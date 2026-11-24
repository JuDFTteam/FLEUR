"""Tests of gfprep for films and of the vacuum embedding of gfleur.

A Cu(001) film of three atomic planes (RMT 1.6, no MT sphere touches an
embedding plane) is converged with fleur. gfprep converts it into
- one layer covering the whole film,
- three one-plane layers,
- one layer with the vacuum embedding planes 1 a.u. beyond the FLEUR vacuum
  boundary z1 (zvac=1).
gfleur embeds each setup in the vacuum on both sides (&EMB surface=T,
vacleft=T) and computes the DOS at Gamma on a line of energies with a small
imaginary part. Its peaks reproduce the eigenvalues of the film from a fleur
run at Gamma on the same density, independent of the split into layers and of
the position of the vacuum embedding planes.

Surface embedding: the upper two planes of the film between a bulk lead (CBS
of a Cu supercell with the same MT radius and cutoffs, potential shifted to
the Fermi energy of the film) and the vacuum.

Binaries are located via GFPREP_BINARY or <repo>/build*/gfprep; fleur and
gfleur are taken from the same directory. Missing binaries skip the tests.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
FILM_INPUTS = HERE / "inputfiles" / "cu_film"
BULK_INPUTS = HERE / "inputfiles" / "cu_supercell"
# the bulk supercell with the MT radius and the cutoffs of the film
BULK_SETUP = {'radius="2.35"': 'radius="1.6"',
              'Kmax="3.4" Gmax="10.2" GmaxXC="10.2"': 'Kmax="4.5" Gmax="13.5" GmaxXC="11.25"',
              'itmax="6"': 'itmax="40"'}
# energies of the DOS: below the d band, Im(E) broadens the peaks
E_MIN, E_MAX, E_STEP, E_IMAG = -0.50, -0.30, 0.002, 0.002
PEAK_TOL = 0.0025
# gfprep options, number of layers and napw per setup (the thicker layer
# with zvac=1 needs more basis functions along z)
SETUPS = {
    "one_layer": ("planes_per_layer=3", 1, 8),
    "three_layers": ("planes_per_layer=1", 3, 8),
    "one_layer_zvac": ("planes_per_layer=3 zvac=1.0", 1, 12),
}


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


def _gamma_eigenvalues(fleur, film, work):
    """Eigenvalues of the film at Gamma from one fleur iteration on the converged density"""
    shutil.copytree(film, work)
    (work / "out.xml").unlink()
    _replace(work / "kpts.xml", {"<kPointLists>":
                                 '<kPointLists>\n      <kPointList name="gamma" count="1" type="unspecified">\n'
                                 '         <kPoint weight="1.0">0.0 0.0 0.0</kPoint>\n      </kPointList>'})
    inp = (work / "inp.xml").read_text()
    inp = re.sub(r'kPointListSelection listName="[^"]*"', 'kPointListSelection listName="gamma"', inp)
    inp = re.sub(r'itmax="\d+"', 'itmax="1"', inp)
    (work / "inp.xml").write_text(inp)
    proc, log = _run(fleur, work)
    assert proc.returncode == 0, f"fleur run at Gamma failed\n{log[-3000:]}"
    block = re.search(r"<eigenvaluesAt[^>]*>(.*?)</eigenvaluesAt>", (work / "out.xml").read_text(), re.S)
    return sorted(float(x) for x in block.group(1).split())


@pytest.fixture(scope="module")
def film(tmp_path_factory):
    """Converge the film and get its eigenvalues at Gamma (once per module)."""
    gfprep = _find_gfprep()
    if gfprep is None:
        pytest.skip("no gfprep binary found (set GFPREP_BINARY or build the gfprep target)")
    fleur = _sibling(gfprep, ["fleur", "fleur_MPI"])
    gfleur = _sibling(gfprep, ["gfleur", "gfleur_MPI"])
    if fleur is None or gfleur is None:
        pytest.skip(f"no fleur or gfleur binary next to {gfprep}")
    base = tmp_path_factory.mktemp("film")
    scf = base / "film"
    shutil.copytree(FILM_INPUTS, scf)
    proc, log = _run(fleur, scf)
    assert proc.returncode == 0, f"film fleur run failed\n{log[-3000:]}"
    assert (scf / "cdn.hdf").is_file()
    eig = _gamma_eigenvalues(fleur, scf, base / "gamma")
    return gfprep, gfleur, base, eig


@pytest.fixture(scope="module", params=list(SETUPS))
def setup(request, film):
    """gfprep and the gfleur DOS run of one setup"""
    gfprep, gfleur, base, eig = film
    options, nl, napw = SETUPS[request.param]
    work = base / request.param
    work.mkdir()
    # relative path: the xml file name handling of FLEUR is limited to 100 characters
    (work / "gfprep_inp").write_text(f"&GFPREP bulkdir='../film/' kgrid=1,1 {options} /\n")
    proc, log = _run(gfprep, work)
    assert proc.returncode == 0, f"gfprep failed\n{log[-3000:]}"
    gfprep_gfinp = (work / "gf_inp").read_text()

    _dos_inputs(work, nl, napw)
    proc, log = _run(gfleur, work)
    assert proc.returncode == 0, f"gfleur failed\n{log[-3000:]}"
    return request.param, nl, work, gfprep_gfinp, eig


def _gf_en(work):
    n = round((E_MAX - E_MIN)/E_STEP) + 1
    (work / "gf_en").write_text(f"l {n:14d}\n" + "".join(
        f"{E_MIN + i*E_STEP:10.5f}{E_IMAG:10.5f}{1.0:10.5f}{0.0:10.5f}\n" for i in range(n)))


def _dos_inputs(work, nl, napw):
    """gf_inp in DOS mode, the energies and the DOS projections of all layers"""
    _replace(work / "gf_inp", {"&MODE tmat=t, gmat=t /": f"&MODE tmat=t, gmat=t, dos=t /\n&BASIS napw={napw} /"})
    _gf_en(work)
    (work / "gf_dos").write_text("".join(f"&LDOS layer={i}, atom=0, l=-1 /\n" for i in range(1, nl + 1)))


def _gfprep(gfprep, work, options):
    work.mkdir()
    (work / "gfprep_inp").write_text(f"&GFPREP {options} kgrid=1,1 /\n")
    proc, log = _run(gfprep, work)
    assert proc.returncode == 0, f"gfprep failed\n{log[-3000:]}"


def _layer_files(work, dest):
    dest.mkdir()
    for f in work.iterdir():
        if f.name.startswith("layer") or f.suffix == ".xsd" or f.name in ("gf_setup.hdf", "gf_pot.hdf", "gf_cdn.hdf",
                                                                         "gf_inp"):
            shutil.copy(f, dest / f.name)


def _fermi_energy(out_xml):
    return float(re.findall(r'<FermiEnergy value="([^"]*)"', out_xml.read_text())[-1])


def _checks(out_text):
    """Self-check lines of gfprep.out: (quantity, layer, region) -> (max|f|, deviation)"""
    pat = re.compile(r"gfprep check (\w+) layer (\d+) (interstitial|muffin-tins): "
                     r"max\|f\|, max deviation:\s+(\S+)\s+(\S+)")
    return {(m[0], int(m[1]), m[2]): (float(m[3]), float(m[4])) for m in pat.findall(out_text)}


def _dos_peaks(path):
    """Positions of the maxima of the summed MT DOS (parabolic refinement)"""
    rows = [[float(x) for x in line.split()] for line in path.read_text().splitlines() if line.strip()]
    e = [r[0] for r in rows]
    y = [sum(r[1:]) for r in rows]
    peaks = []
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] >= y[i + 1] and y[i] > 0.01*max(y):
            a, b, c = y[i - 1], y[i], y[i + 1]
            peaks.append(e[i] + 0.5*(e[1] - e[0])*(a - c)/(a - 2*b + c))
    return peaks


def test_gfprep_film_outputs(setup):
    """gfprep marks both outer planes as vacuum embedding planes and writes the vacuum potentials."""
    _, nl, work, gfinp, _ = setup
    assert f"nlayers={nl}" in gfinp
    assert "&EMB addemb=T, surface=T, vacleft=T /" in gfinp
    planes = (work / "gfprep_planes.txt").read_text()
    assert "# vacuum embedding plane below layer 1" in planes
    assert f"# vacuum embedding plane above layer {nl}" in planes
    pat = re.compile(r"vacuum embedding potential side (\d): (\d+) points.*V\(plane\), V\(end\):\s+(\S+)\s+(\S+)")
    vac = {int(m[0]): (float(m[2]), float(m[3])) for m in pat.findall((work / "gfprep.out").read_text())}
    assert sorted(vac) == [1, 2]
    for v_plane, v_end in vac.values():
        # the potential rises from the plane to the vacuum level (close to 0)
        assert v_plane < v_end
        assert abs(v_end) < 1e-2
    # z-reflection symmetric film
    assert vac[1] == pytest.approx(vac[2], abs=1e-6)


def test_gfprep_film_conversion_accuracy(setup):
    """The converted layer functions reproduce the film, including the vacuum part."""
    _, nl, work, _, _ = setup
    checks = _checks((work / "gfprep.out").read_text())
    assert len(checks) == 4*nl
    for (quantity, layer, region), (fmax, dev) in checks.items():
        if region == "muffin-tins":
            assert dev < 1e-10*fmax, f"{quantity} layer {layer} MT deviation {dev}"
        else:
            assert dev < 1e-2*fmax, f"{quantity} layer {layer} interstitial deviation {dev}"


def test_vacuum_embedding_integration(setup):
    """The integration of the vacuum embedding potential is converged."""
    _, _, work, _, _ = setup
    out = (work / "out").read_text()
    diffs = [float(x) for x in re.findall(r"vacuum embedding side \d:.*step-doubling difference:\s+(\S+)", out)]
    assert len(diffs) == 2
    assert max(diffs) < 1e-6


def test_film_dos_peaks(setup):
    """The DOS of the film embedded in vacuum has its peaks at the fleur eigenvalues."""
    name, _, work, _, eig = setup
    peaks = _dos_peaks(work / "gf_DOS.raw")
    ref = [e for e in eig if E_MIN + 0.01 < e < E_MAX - 0.01]
    assert len(ref) >= 3
    for e in ref:
        nearest = min(peaks, key=lambda p: abs(p - e))
        assert abs(nearest - e) < PEAK_TOL, f"{name}: fleur eigenvalue {e:.4f}, nearest DOS peak {nearest:.4f}"


def test_surface_embedding(film):
    """Bulk lead on the left, vacuum on the right: only the left lead file is read."""
    gfprep, gfleur, base, _ = film
    fleur = _sibling(gfprep, ["fleur", "fleur_MPI"])
    bulk = base / "bulk"
    shutil.copytree(BULK_INPUTS, bulk)
    _replace(bulk / "inp.xml", BULK_SETUP)
    proc, log = _run(fleur, bulk)
    assert proc.returncode == 0, f"bulk fleur run failed\n{log[-3000:]}"
    vshift = _fermi_energy(base / "film" / "out.xml") - _fermi_energy(bulk / "out.xml")

    # lead: the two-plane bulk period that continues the film below its plane 2
    _gfprep(gfprep, base / "lead_prep", f"bulkdir='../bulk/' planes_per_layer=2 first_plane=2 last_plane=3 "
                                        f"vshift={vshift:.6f}")
    cbs = base / "lead_cbs"
    _layer_files(base / "lead_prep", cbs)
    text = (cbs / "gf_inp").read_text()
    (cbs / "gf_inp").write_text(text[:text.index("\n/\n") + 3] +
                                "&MODE tmat=t, gmat=t /\n&BASIS napw=8 /\n&CBSMODE cbs=t /\n&SC iter=1 /\n")
    _gf_en(cbs)
    (cbs / "gf_io").write_text("default\n")
    proc, log = _run(gfleur, cbs)
    assert proc.returncode == 0, f"gfleur CBS run failed\n{log[-3000:]}"

    # surface: planes 2 and 3 of the film, vacuum above
    _gfprep(gfprep, base / "surf_prep", "bulkdir='../film/' planes_per_layer=1 first_plane=2")
    assert "&EMB addemb=T, surface=T, vacleft=F /" in (base / "surf_prep" / "gf_inp").read_text()
    run = base / "surface"
    _layer_files(base / "surf_prep", run)
    shutil.copy(cbs / "embpot11.hdf", run / "embpot11.hdf")
    (run / "gf_io").write_text("default\n")
    _dos_inputs(run, 2, 8)
    proc, log = _run(gfleur, run)
    assert proc.returncode == 0, f"gfleur surface run failed\n{log[-3000:]}"
    out = (run / "out").read_text()
    assert "vacuum embedding side 2" in out
    assert "vacuum embedding side 1" not in out
    # the semi-infinite bulk gives a continuous DOS
    rows = [[float(x) for x in line.split()] for line in (run / "gf_DOS.raw").read_text().splitlines() if line.strip()]
    assert len(rows) == round((E_MAX - E_MIN)/E_STEP) + 1
    assert min(sum(r[1:]) for r in rows) > 0.1
