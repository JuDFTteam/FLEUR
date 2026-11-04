"""End-to-end test of the embedding workflow on a perfect Cu crystal.

1. A tetragonal 4-atom Cu supercell (ABAB stacking along (001), so the lattice
   period is two atomic planes) is converged with fleur.
2. gfprep cuts one two-plane period as a single layer; gfleur in CBS mode
   (`&CBSMODE cbs=t`) writes the lead embedding potentials embpot11/embpotN2.
3. gfleur computes the transmission (Landauer) of
   a) the same single layer between these leads, and
   b) a stack of two layers of two atomic planes each,
   between the leads. For a perfect crystal the transmission at every k is the
   number of open channels, i.e. the number of right-moving Bloch states of the CBS.

Binaries are located via GFPREP_BINARY or <repo>/build*/gfprep; fleur and
gfleur are taken from the same directory. Missing binaries skip the tests.
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
BULK_INPUTS = HERE / "inputfiles" / "cu_supercell"
ENERGIES = (0.1, 0.2)
GF_EN = "l              2\n   0.10000   0.00000   1.00000   0.00000\n   0.20000   0.00000   1.00000   0.00000\n"


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


def _gfprep(gfprep, work, bulk, options):
    work.mkdir()
    # relative path: the xml file name handling of FLEUR is limited to 100 characters
    (work / "gfprep_inp").write_text(f"&GFPREP bulkdir='../{bulk.name}/' {options} /\n")
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
    return layers + mode_line + "\n" + extra + "&SC iter=1 /\n"


def _transmission(path):
    """{(kx, ky, E): j} of the first transcurrent column"""
    res = {}
    for line in path.read_text().splitlines():
        p = line.split()
        res[(round(float(p[1]), 4), round(float(p[2]), 4), round(float(p[3]), 4))] = float(p[4])
    return res


@pytest.fixture(scope="module")
def lead(tmp_path_factory):
    """Bulk run, CBS of one two-plane period -> lead embedding potentials."""
    gfprep = _find_gfprep()
    if gfprep is None:
        pytest.skip("no gfprep binary found (set GFPREP_BINARY or build the gfprep target)")
    fleur = _sibling(gfprep, ["fleur", "fleur_MPI"])
    gfleur = _sibling(gfprep, ["gfleur", "gfleur_MPI"])
    if fleur is None or gfleur is None:
        pytest.skip(f"fleur/gfleur binaries missing next to {gfprep}")
    base = tmp_path_factory.mktemp("transport")
    bulk = base / "bulk"
    shutil.copytree(BULK_INPUTS, bulk)
    proc, log = _run(fleur, bulk)
    assert proc.returncode == 0 and (bulk / "cdn.hdf").is_file(), f"bulk fleur run failed\n{log[-3000:]}"

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
    return gfprep, gfleur, bulk, base, cbs


def _transport(lead, name, prep_options, mode, nlayers):
    gfprep, gfleur, bulk, base, cbs = lead
    prep = base / f"{name}_prep"
    _gfprep(gfprep, prep, bulk, prep_options)
    run = base / name
    run.mkdir()
    _layer_files(prep, run)
    for f in ("embpot11.hdf", "embpotN2.hdf"):
        shutil.copy(cbs / f, run / f)
    (run / "gf_en").write_text(GF_EN)
    (run / "gf_inp").write_text(_gf_inp(prep, f"&MODE tmat=f, gmat=t, {mode} /", "&EMB addemb=t /\n"))
    (run / "gf_io").write_text(
        "2\n&IO mode='E', region=1, side=1, rw='r', id=1 /\n"
        f"&IO mode='E', region={nlayers}, side=2, rw='r', id=2 /\n"
        "&FILE name='embpot11' /\n&FILE name='embpotN2' /\n")
    proc, log = _run(gfleur, run)
    assert proc.returncode == 0, f"gfleur transport run failed\n{log[-3000:]}"
    return _transmission(run / "transcurrent")


def _open_channels(cbs, ik=1):
    """Number of right-moving Bloch states per energy at k-point ik (CBS_bloch.1 lists +-k pairs)."""
    counts = {}
    for line in (cbs / "CBS_bloch.1").read_text().splitlines():
        p = line.split()
        if int(p[1]) == ik:
            counts[round(float(p[0]), 4)] = counts.get(round(float(p[0]), 4), 0) + 1
    return {e: n // 2 for e, n in counts.items()}


def test_cbs_channels_at_gamma(lead):
    """The CBS finds the expected Cu channels (3 and 4 at E=0.1 and 0.2 Htr at Gamma)."""
    assert _open_channels(lead[4]) == {0.1: 3, 0.2: 4}


def _bloch_states(cbs):
    """{(E, ik): sorted kz} of the Bloch states in CBS_bloch.1"""
    states = {}
    for line in (cbs / "CBS_bloch.1").read_text().splitlines():
        p = line.split()
        states.setdefault((round(float(p[0]), 4), int(p[1])), []).append(float(p[2]))
    return {key: sorted(kz) for key, kz in states.items()}


def test_cbs_spectral_vs_linear_solve(lead):
    """spectral=f (linear solve) and the default spectral=t give the same CBS."""
    gfprep, gfleur, bulk, base, cbs = lead
    run = base / "cbs_nospectral"
    run.mkdir()
    _layer_files(base / "lead_prep", run)
    (run / "gf_inp").write_text(_gf_inp(base / "lead_prep", "&MODE tmat=t, gmat=t, spectral=f /",
                                        "&CBSMODE cbs=t /\n"))
    (run / "gf_en").write_text(GF_EN)
    (run / "gf_io").write_text("default\n")
    proc, log = _run(gfleur, run)
    assert proc.returncode == 0, f"gfleur CBS run with spectral=f failed\n{log[-3000:]}"
    ref, res = _bloch_states(cbs), _bloch_states(run)
    assert res.keys() == ref.keys()
    for key in ref:
        assert res[key] == pytest.approx(ref[key], abs=1e-5), key


def test_transmission_single_layer(lead):
    """One two-plane layer between its own leads transmits every open channel."""
    t = _transport(lead, "stack1", "planes_per_layer=2 first_plane=1 last_plane=2", "curr=-1", 1)
    channels = _open_channels(lead[4])
    for e in ENERGIES:
        assert t[(0.0, 0.0, e)] == pytest.approx(channels[e], abs=1e-3)
    assert all(abs(v - round(v)) < 1e-3 for v in t.values() if abs(v) > 1e-3)


@pytest.mark.parametrize("mode", ["curr=3", "curr=16"])
def test_transmission_two_layer_stack(lead, mode):
    """A stack of two two-plane layers (two lattice periods) between the leads: propagation of
    the embedding potentials. (Layers of a single atomic plane are thinner than the muffin-tin
    diameter and are less accurate, they give e.g. 2.65 instead of 3 channels.)"""
    t = _transport(lead, f"stack2_{mode[5:]}", "planes_per_layer=2 first_plane=1 last_plane=4", mode, 2)
    channels = _open_channels(lead[4])
    for e in ENERGIES:
        assert t[(0.0, 0.0, e)] == pytest.approx(channels[e], abs=1e-3)
