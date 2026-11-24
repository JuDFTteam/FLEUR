"""End-to-end test of the embedding workflow on a perfect Cu crystal.

1. A tetragonal 4-atom Cu supercell (ABAB stacking along (001), so the lattice
   period is two atomic planes) is converged with fleur.
2. gfprep cuts one two-plane period as a single layer; gfleur in CBS mode
   (`&CBSMODE cbs=t`) writes the lead embedding potentials embpot11/embpotN2.
   The propagating states of the CBS at Gamma are compared with the bulk band
   structure along (0,0,kz) from a fleur band run on the same density.
3. gfleur computes the transmission (Landauer) of
   a) the same single layer between these leads, and
   b) a stack of two layers of two atomic planes each,
   between the leads. For a perfect crystal the transmission at every k is the
   number of open channels, i.e. the number of right-moving Bloch states of the CBS.

Everything runs for two muffin-tin radii:
- rmt2.35: the standard radius. The embedding planes lie midway between the
  atomic planes (1.705 a.u. from the atoms) and cut the MT spheres.
- rmt1.6: no MT sphere touches the embedding planes (Kmax*RMT=7.2).

Binaries are located via GFPREP_BINARY or <repo>/build*/gfprep; fleur and
gfleur are taken from the same directory. Missing binaries skip the tests.
"""
import math
import os
import re
import shutil
import subprocess
from collections import namedtuple
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
BULK_INPUTS = HERE / "inputfiles" / "cu_supercell"
ENERGIES = (0.1, 0.2)
GF_EN = "l              2\n   0.10000   0.00000   1.00000   0.00000\n   0.20000   0.00000   1.00000   0.00000\n"
# replacements in the bulk inp.xml per MT setup
MT_SETUPS = {
    "rmt2.35": {},
    "rmt1.6": {'radius="2.35"': 'radius="1.6"',
               'Kmax="3.4" Gmax="10.2" GmaxXC="10.2"': 'Kmax="4.5" Gmax="13.5" GmaxXC="11.25"'},
}
# small cutoffs: the tests compare the CBS with the bulk bands of the same density, they
# need not be converged; 2D k-points Gamma, X, M
NAPW = 8
KGRID = "kgrid=2,2"
# k-points of the bulk band run: (0,0,kz) with kz in [0,pi/L]
NKZ = 101
KZ_TOL = 0.01

Lead = namedtuple("Lead", "mt gfprep gfleur bulk base cbs bulk_kz")


def _find_gfprep():
    env = os.environ.get("GFPREP_BINARY")
    if env and Path(env).is_file():
        return Path(env)
    for c in sorted(REPO.glob("build*/gfprep")):
        if c.is_file() and os.access(c, os.X_OK):
            return c
    return None


def _sibling(binary, names):
    for name in names:
        c = binary.parent / name
        if c.is_file() and os.access(c, os.X_OK):
            return c
    return None


def _run(cmd, cwd, timeout=1800):
    proc = subprocess.run([str(cmd)], cwd=cwd, capture_output=True, text=True,
                          timeout=timeout, env=dict(os.environ))
    return proc, proc.stdout + "\n" + proc.stderr


def _replace(path, replacements):
    text = path.read_text()
    for old, new in replacements.items():
        assert text.count(old) == 1, f"{old} not found once in {path}"
        text = text.replace(old, new)
    path.write_text(text)


def _gfprep(gfprep, work, bulk, options):
    work.mkdir()
    # relative path: the xml file name handling of FLEUR is limited to 100 characters
    (work / "gfprep_inp").write_text(f"&GFPREP bulkdir='../{bulk.name}/' {options} {KGRID} /\n")
    proc, log = _run(gfprep, work)
    assert proc.returncode == 0, f"gfprep failed\n{log[-3000:]}"


def _layer_files(work, dest):
    for f in work.iterdir():
        if f.name.startswith("layer") or f.suffix == ".xsd" or f.name in ("gf_setup.hdf", "gf_pot.hdf", "gf_cdn.hdf"):
            shutil.copy(f, dest / f.name)


def _gf_inp(work, mode_line, extra=""):
    """gfprep's &LAYERS block followed by the given mode line."""
    text = (work / "gf_inp").read_text()
    layers = text[:text.index("\n/\n") + 3]
    return layers + mode_line + "\n" + f"&BASIS napw={NAPW} /\n" + extra + "&SC iter=1 /\n"


def _transmission(path):
    """{(kx, ky, E): j} of the first transcurrent column"""
    res = {}
    for line in path.read_text().splitlines():
        p = line.split()
        res[(round(float(p[1]), 4), round(float(p[2]), 4), round(float(p[3]), 4))] = float(p[4])
    return res


def _bulk_period(bulk):
    row3 = re.search(r"<row-3>(.*?)</row-3>", (bulk / "inp.xml").read_text()).group(1)
    return float(row3.split()[2])


def _bulk_band_kz(fleur, bulk, work):
    """{E: sorted kz} of the right-moving bulk states at Gamma for E in ENERGIES:
    crossings of the sorted bands with E for kz in (0,pi/L)."""
    shutil.copytree(bulk, work)
    (work / "out.xml").unlink()
    pts = "".join(f'               <kPoint weight="1.0">0.0 0.0 {0.5 * i / (NKZ - 1):.10f}</kPoint>\n'
                  for i in range(NKZ))
    _replace(work / "kpts.xml", {"</kPointLists>":
                                 f'   <kPointList name="kz" count="{NKZ}" type="path">\n{pts}'
                                 "            </kPointList>\n         </kPointLists>"})
    inp = (work / "inp.xml").read_text()
    inp = re.sub(r'kPointListSelection listName="[^"]*"', 'kPointListSelection listName="kz"', inp)
    (work / "inp.xml").write_text(inp)
    _replace(work / "inp.xml", {'<output dos="F" band="F"': '<output dos="F" band="T"'})
    proc, log = _run(fleur, work)
    assert proc.returncode == 0, f"bulk band run failed\n{log[-3000:]}"
    blocks = re.findall(r"<eigenvaluesAt[^>]*>(.*?)</eigenvaluesAt>", (work / "out.xml").read_text(), re.S)
    bands = [sorted(float(x) for x in b.split()) for b in blocks]
    assert len(bands) == NKZ
    nb = min(len(b) for b in bands)
    dkz = 0.5 * 2 * math.pi / _bulk_period(bulk) / (NKZ - 1)
    kz = {}
    for e in ENERGIES:
        kz[e] = sorted((i + (e - b0[n]) / (b1[n] - b0[n])) * dkz
                       for i, (b0, b1) in enumerate(zip(bands[:-1], bands[1:]))
                       for n in range(nb) if (b0[n] - e) * (b1[n] - e) < 0)
    return kz


@pytest.fixture(scope="module", params=list(MT_SETUPS))
def lead(request, tmp_path_factory):
    """Bulk run, bulk bands at Gamma, CBS of one two-plane period -> lead embedding potentials."""
    gfprep = _find_gfprep()
    if gfprep is None:
        pytest.skip("no gfprep binary found (set GFPREP_BINARY or build the gfprep target)")
    fleur = _sibling(gfprep, ["fleur", "fleur_MPI"])
    gfleur = _sibling(gfprep, ["gfleur", "gfleur_MPI"])
    if fleur is None or gfleur is None:
        pytest.skip(f"fleur/gfleur binaries missing next to {gfprep}")
    base = tmp_path_factory.mktemp(f"transport_{request.param}")
    bulk = base / "bulk"
    shutil.copytree(BULK_INPUTS, bulk)
    _replace(bulk / "inp.xml", MT_SETUPS[request.param])
    proc, log = _run(fleur, bulk)
    assert proc.returncode == 0 and (bulk / "cdn.hdf").is_file(), f"bulk fleur run failed\n{log[-3000:]}"
    bulk_kz = _bulk_band_kz(fleur, bulk, base / "bulk_band")

    one = base / "lead_prep"
    _gfprep(gfprep, one, bulk, "planes_per_layer=2 first_plane=1 last_plane=2")
    cbs = base / "cbs"
    cbs.mkdir()
    _layer_files(one, cbs)
    (cbs / "gf_inp").write_text(_gf_inp(one, "&MODE tmat=t, gmat=t /", "&CBSMODE cbs=t /\n"))
    (cbs / "gf_en").write_text(GF_EN)
    (cbs / "gf_io").write_text("default\n")
    proc, log = _run(gfleur, cbs)
    assert proc.returncode == 0, f"gfleur CBS run failed\n{log[-3000:]}"
    for f in ("embpot11.hdf", "embpotN2.hdf", "CBS_bloch.1"):
        assert (cbs / f).is_file(), f"{f} missing"
    # the defaults of gf_io must not have been written into the user's file
    assert (cbs / "gf_io").read_text().strip() == "default"
    return Lead(request.param, gfprep, gfleur, bulk, base, cbs, bulk_kz)


def _transport(lead, name, prep_options, mode, nlayers):
    prep = lead.base / f"{name}_prep"
    _gfprep(lead.gfprep, prep, lead.bulk, prep_options)
    run = lead.base / name
    run.mkdir()
    _layer_files(prep, run)
    for f in ("embpot11.hdf", "embpotN2.hdf"):
        shutil.copy(lead.cbs / f, run / f)
    (run / "gf_en").write_text(GF_EN)
    (run / "gf_inp").write_text(_gf_inp(prep, f"&MODE tmat=f, gmat=t, {mode} /", "&EMB addemb=t /\n"))
    (run / "gf_io").write_text(
        "2\n&IO mode='E', region=1, side=1, rw='r', id=1 /\n"
        f"&IO mode='E', region={nlayers}, side=2, rw='r', id=2 /\n"
        "&FILE name='embpot11' /\n&FILE name='embpotN2' /\n")
    proc, log = _run(lead.gfleur, run)
    assert proc.returncode == 0, f"gfleur transport run failed\n{log[-3000:]}"
    return _transmission(run / "transcurrent")


def _bloch_states(cbs):
    """{(E, ik): [(kz, current)]} of the Bloch states in CBS_bloch.1"""
    states = {}
    for line in (cbs / "CBS_bloch.1").read_text().splitlines():
        p = line.split()
        states.setdefault((round(float(p[0]), 4), int(p[1])), []).append((float(p[2]), float(p[3])))
    return states


def _open_channels(cbs, ik=1):
    """Number of right-moving Bloch states per energy at k-point ik (CBS_bloch.1 lists +-k pairs)."""
    states = _bloch_states(cbs)
    return {e: len(states.get((e, ik), [])) // 2 for e in ENERGIES}


def _cut_spheres_xfail(lead, energy, request):
    # E=0.1 lies close to the extremum of a flat band (bulk kz=0.15), where the
    # error of the cut spheres turns the propagating state into an evanescent one
    if lead.mt == "rmt2.35" and energy == 0.1:
        request.applymarker(pytest.mark.xfail(
            reason="the embedding planes cut the MT spheres, the CBS is not converged in napw"))


@pytest.mark.parametrize("energy", ENERGIES)
def test_cbs_channels_at_gamma(lead, energy, request):
    """The CBS has as many propagating channels at Gamma as the bulk band structure."""
    _cut_spheres_xfail(lead, energy, request)
    assert _open_channels(lead.cbs)[energy] == len(lead.bulk_kz[energy])


@pytest.mark.parametrize("energy", ENERGIES)
def test_cbs_kz_at_gamma(lead, energy, request):
    """kz of the right-moving CBS states at Gamma agree with the bulk band structure
    (folded into the Brillouin zone of the 4-atom bulk cell)."""
    _cut_spheres_xfail(lead, energy, request)
    period = 2 * math.pi / _bulk_period(lead.bulk)
    kz = sorted(min(k % period, period - k % period)
                for k, current in _bloch_states(lead.cbs).get((energy, 1), []) if current > 0)
    assert kz == pytest.approx(lead.bulk_kz[energy], abs=KZ_TOL)


def test_cbs_spectral_vs_linear_solve(lead):
    """spectral=f (linear solve) and the default spectral=t give the same CBS."""
    run = lead.base / "cbs_nospectral"
    run.mkdir()
    _layer_files(lead.base / "lead_prep", run)
    (run / "gf_inp").write_text(_gf_inp(lead.base / "lead_prep", "&MODE tmat=t, gmat=t, spectral=f /",
                                        "&CBSMODE cbs=t /\n"))
    (run / "gf_en").write_text(GF_EN)
    (run / "gf_io").write_text("default\n")
    proc, log = _run(lead.gfleur, run)
    assert proc.returncode == 0, f"gfleur CBS run with spectral=f failed\n{log[-3000:]}"
    ref, res = _bloch_states(lead.cbs), _bloch_states(run)
    assert res.keys() == ref.keys()
    for key in ref:
        assert sorted(k for k, _ in res[key]) == pytest.approx(sorted(k for k, _ in ref[key]), abs=1e-5), key


def test_transmission_single_layer(lead):
    """One two-plane layer between its own leads transmits every open channel."""
    t = _transport(lead, "stack1", "planes_per_layer=2 first_plane=1 last_plane=2", "curr=-1", 1)
    channels = _open_channels(lead.cbs)
    for e in ENERGIES:
        assert t[(0.0, 0.0, e)] == pytest.approx(channels[e], abs=1e-3)
    assert all(abs(v - round(v)) < 1e-3 for v in t.values() if abs(v) > 1e-3)


@pytest.mark.parametrize("mode", ["curr=3", "curr=16"])
def test_transmission_two_layer_stack(lead, mode):
    """A stack of two two-plane layers (two lattice periods) between the leads: propagation of
    the embedding potentials. (Layers of a single atomic plane are thinner than the muffin-tin
    diameter and are less accurate, they give e.g. 2.65 instead of 3 channels.)"""
    t = _transport(lead, f"stack2_{mode[5:]}", "planes_per_layer=2 first_plane=1 last_plane=4", mode, 2)
    channels = _open_channels(lead.cbs)
    for e in ENERGIES:
        assert t[(0.0, 0.0, e)] == pytest.approx(channels[e], abs=1e-3)
