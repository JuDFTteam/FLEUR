"""Non-collinear (noco) tests of gfprep and gfleur.

Without spin-orbit coupling a global rotation of all moments leaves the
spectrum unchanged, and a noco calculation with all moments along z equals the
collinear one. The tests use this on
- sc Fe with one atom (testing/inputfiles/noco/1atz, 1atx: moment along z or x,
  local orbitals): gfprep cuts one period as one layer; the CBS at Gamma of the
  x and z setups are equal, the z layer run collinearly gives the same states
  (spin up and down together), and the layer between its own CBS leads
  transmits every open channel;
- sc Fe with two atoms of perpendicular moments (noco/FeAPWlo_noco) and the MT
  off-diagonal potential (l_mtNocoPot, no symmetry), rotated as a whole by
  90 degrees: the complex band structure is unchanged;
- the Mn monolayer of noco/MnFilmX (moment along x) and the same film with the
  moment along z, embedded in the vacuum on both sides (2x2 vacuum embedding
  potential): equal DOS peaks, and the z film run collinearly gives the same DOS.

With spin-orbit coupling (first-variation SOC of the noco Hamiltonian) the same Mn film
(moment along z, symmetry reduced for SOC) is compared with FLEUR: the SOC splitting of a
d-level pair at Gamma, and sc Fe with SOC transmits every open channel of its own leads.

Binaries are located via GFPREP_BINARY or <repo>/build*/gfprep; fleur and
gfleur are taken from the same directory. Missing binaries skip the tests.
"""
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent
NOCO_INPUTS = REPO / "testing" / "inputfiles" / "noco"
SCHEMAS = [REPO / "src" / "fleur" / "io" / "xml" / f for f in ("FleurInputSchema.xsd", "FleurOutputSchema.xsd")]
GF_EN = "l              2\n   0.10000   0.00000   1.00000   0.00000\n   0.20000   0.00000   1.00000   0.00000\n"
NAPW = 8
# DOS of the film: d states at Gamma, Im(E) broadens the peaks
E_MIN, E_MAX, E_STEP, E_IMAG = -0.34, -0.07, 0.002, 0.002


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


def _fleur_run(fleur, src, work, replacements, files=("inp.xml",)):
    work.mkdir()
    for f in files:
        shutil.copy(src / f, work / f)
    for f in SCHEMAS:
        shutil.copy(f, work / f.name)
    _replace(work / "inp.xml", replacements)
    proc, log = _run(fleur, work)
    assert proc.returncode == 0 and (work / "cdn.hdf").is_file(), f"fleur run in {work} failed\n{log[-3000:]}"


def _gfprep(gfprep, work, bulk, options):
    work.mkdir()
    for f in SCHEMAS:
        shutil.copy(f, work / f.name)
    # relative path: the xml file name handling of FLEUR is limited to 100 characters
    (work / "gfprep_inp").write_text(f"&GFPREP bulkdir='../{bulk.name}/' kgrid=1,1 {options} /\n")
    proc, log = _run(gfprep, work)
    assert proc.returncode == 0, f"gfprep failed\n{log[-3000:]}"


def _layer_files(work, dest):
    dest.mkdir()
    for f in work.iterdir():
        if f.name.startswith("layer") or f.suffix == ".xsd" or f.name in ("gf_setup.hdf", "gf_pot.hdf", "gf_cdn.hdf"):
            shutil.copy(f, dest / f.name)


def _layers_block(prep):
    text = (prep / "gf_inp").read_text()
    return text[:text.index("\n/\n") + 3]


def _cbs(gfleur, prep, run, gf_en=GF_EN, collinear=False):
    """gfleur in CBS mode on the layers of prep"""
    _layer_files(prep, run)
    if collinear:
        for f in run.glob("layer*_inp.xml"):
            _replace(f, {'l_noco="T"': 'l_noco="F"'})
    (run / "gf_inp").write_text(_layers_block(prep) + f"&MODE tmat=t, gmat=t /\n&BASIS napw={NAPW} /\n"
                                "&CBSMODE cbs=t /\n&SC iter=1 /\n")
    (run / "gf_en").write_text(gf_en)
    (run / "gf_io").write_text("default\n")
    proc, log = _run(gfleur, run)
    assert proc.returncode == 0, f"gfleur CBS run in {run} failed\n{log[-3000:]}"


def _bloch_states(path):
    """{(E, ik): sorted [(kz, current)]} of the Bloch states in a CBS_bloch file"""
    states = {}
    for line in path.read_text().splitlines():
        p = line.split()
        states.setdefault((round(float(p[0]), 4), int(p[1])), []).append((float(p[2]), float(p[3])))
    return {k: sorted(v) for k, v in states.items()}


def _cbs_kz(path):
    """{E: [complex kz]} of all states in a CBS_eigenvalues file"""
    res = {}
    for line in path.read_text().splitlines():
        p = line.split()
        res.setdefault(round(float(p[1]), 4), []).append(complex(float(p[4]), float(p[5])))
    return res


def _checks(out_text):
    """Self-check lines of gfprep.out: (quantity, layer, region) -> (max|f|, deviation)"""
    pat = re.compile(r"gfprep check (\w+) layer (\d+) (interstitial|muffin-tins): "
                     r"max\|f\|, max deviation:\s+(\S+)\s+(\S+)")
    return {(m[0], int(m[1]), m[2]): (float(m[3]), float(m[4])) for m in pat.findall(out_text)}


def _binaries():
    gfprep = _find_gfprep()
    if gfprep is None:
        pytest.skip("no gfprep binary found (set GFPREP_BINARY or build the gfprep target)")
    fleur = _sibling(gfprep, ["fleur", "fleur_MPI"])
    gfleur = _sibling(gfprep, ["gfleur", "gfleur_MPI"])
    if fleur is None or gfleur is None:
        pytest.skip(f"no fleur or gfleur binary next to {gfprep}")
    return gfprep, fleur, gfleur


@pytest.fixture(scope="module")
def fe1(tmp_path_factory):
    """sc Fe with the moment along z and x: bulk, gfprep (one period), CBS at Gamma"""
    gfprep, fleur, gfleur = _binaries()
    base = tmp_path_factory.mktemp("noco_fe1")
    for d in ("z", "x"):
        _fleur_run(fleur, NOCO_INPUTS / f"1at{d}", base / f"bulk_{d}",
                   {'itmax="2"': 'itmax="25"', 'alpha=".05000000"': 'alpha=".10000000"'})
        _gfprep(gfprep, base / f"prep_{d}", base / f"bulk_{d}", "planes_per_layer=1")
        _cbs(gfleur, base / f"prep_{d}", base / f"cbs_{d}")
    return gfprep, gfleur, base


def test_gfprep_noco_layer(fe1):
    """The layer inp.xml carries the noco setup of the bulk, the files the off-diagonal parts."""
    _, _, base = fe1
    for d, beta in (("z", 0.0), ("x", math.pi/2)):
        xml = (base / f"prep_{d}" / "layer1_inp.xml").read_text()
        assert 'l_noco="T"' in xml and 'l_mperp="T"' in xml
        m = re.search(r'<nocoParams\s+alpha="([^"]*)"\s+beta="([^"]*)"', xml)
        assert float(m[1]) == pytest.approx(0.0, abs=1e-10)
        assert float(m[2]) == pytest.approx(beta, abs=1e-10)
        assert xml.count("<lo ") == 2
    h5py = pytest.importorskip("h5py")
    for f in ("gf_pot.hdf", "gf_cdn.hdf"):
        with h5py.File(base / "prep_x" / f, "r") as h:
            assert "pwod" in h["layer-1/pot"] and "vrod" in h["layer-1/pot"]


def test_gfprep_noco_conversion_accuracy(fe1):
    """The converted layer functions, including the off-diagonal (21) parts, reproduce the bulk."""
    _, _, base = fe1
    checks = _checks((base / "prep_x" / "gfprep.out").read_text())
    assert {q for q, _, _ in checks} == {"potential", "density", "potential21", "density21"}
    # the moment along x gives a finite off-diagonal interstitial potential
    assert checks[("potential21", 1, "interstitial")][0] > 1e-3
    for (quantity, layer, region), (fmax, dev) in checks.items():
        if region == "muffin-tins":
            assert dev <= 1e-10*fmax, f"{quantity} layer {layer} MT deviation {dev}"
        else:
            assert dev < 1e-2*fmax, f"{quantity} layer {layer} interstitial deviation {dev}"


def test_cbs_rotation_invariance(fe1):
    """The CBS at Gamma does not depend on the direction of the moment."""
    _, _, base = fe1
    z = _bloch_states(base / "cbs_z" / "CBS_bloch.1")
    x = _bloch_states(base / "cbs_x" / "CBS_bloch.1")
    assert z.keys() == x.keys()
    assert sum(len(v) for v in z.values()) > 0
    for key in z:
        assert len(x[key]) == len(z[key]), key
        for (kz, cz), (kx, cx) in zip(z[key], x[key]):
            assert kx == pytest.approx(kz, abs=1e-5) and cx == pytest.approx(cz, abs=1e-5), key


def test_cbs_collinear_equivalence(fe1):
    """The z layer run collinearly has the noco states: spin up and down together."""
    _, gfleur, base = fe1
    run = base / "cbs_z_collinear"
    _cbs(gfleur, base / "prep_z", run, collinear=True)
    noco = _bloch_states(base / "cbs_z" / "CBS_bloch.1")
    up, down = _bloch_states(run / "CBS_bloch.1"), _bloch_states(run / "CBS_bloch.2")
    for key in noco:
        both = sorted(up.get(key, []) + down.get(key, []))
        assert [k for k, _ in both] == pytest.approx([k for k, _ in noco[key]], abs=1e-6), key
        # spin up and spin down states are both present
        assert up.get(key) and down.get(key), key


def test_transmission_noco(fe1):
    """The x layer between its own CBS leads transmits every open channel."""
    _, gfleur, base = fe1
    run = base / "transport_x"
    _layer_files(base / "prep_x", run)
    for f in ("embpot11.hdf", "embpotN2.hdf"):
        shutil.copy(base / "cbs_x" / f, run / f)
    (run / "gf_en").write_text(GF_EN)
    (run / "gf_inp").write_text(_layers_block(base / "prep_x") + "&MODE tmat=f, gmat=t, curr=-1 /\n"
                                f"&BASIS napw={NAPW} /\n&EMB addemb=t /\n&SC iter=1 /\n")
    (run / "gf_io").write_text("2\n&IO mode='E', region=1, side=1, rw='r', id=1 /\n"
                               "&IO mode='E', region=1, side=2, rw='r', id=2 /\n"
                               "&FILE name='embpot11' /\n&FILE name='embpotN2' /\n")
    proc, log = _run(gfleur, run)
    assert proc.returncode == 0, f"gfleur transport run failed\n{log[-3000:]}"
    states = _bloch_states(base / "cbs_x" / "CBS_bloch.1")
    for line in (run / "transcurrent").read_text().splitlines():
        p = line.split()
        e = round(float(p[3]), 4)
        channels = len(states.get((e, 1), [])) // 2
        assert channels > 0
        assert float(p[4]) == pytest.approx(channels, abs=1e-3), f"E={e}"


def test_cbs_rotation_two_atoms_mtnocopot(tmp_path):
    """Two perpendicular moments with the MT off-diagonal potential (l_mtNocoPot), rotated
    by 90 degrees about y: the complex band structure at Gamma is unchanged. The density is
    not converged (3 iterations); both setups follow the same rotated path."""
    gfprep, fleur, gfleur = _binaries()
    ident = ("<symmetryOperations>\n         <symOp>\n            <row-1>1 0 0 .0000000000</row-1>\n"
             "            <row-2>0 1 0 .0000000000</row-2>\n            <row-3>0 0 1 .0000000000</row-3>\n"
             "         </symOp>\n      </symmetryOperations>")
    n = 4
    pts = "".join(f'            <kPoint weight="1.0">{i/n:.6f} {j/n:.6f} {k/n:.6f}</kPoint>\n'
                  for i in range(n) for j in range(n) for k in range(n))
    kpts = (f'<kPointLists>\n         <kPointList name="default-1" count="{n**3}" type="unspecified">\n'
            f'{pts}         </kPointList>\n      </kPointLists>')
    angles = {"a": ('beta="0.0"', 'beta="1.57079633"'), "b": ('beta="1.570796326795"', 'beta="3.14159265359"')}
    for v, (b1, b2) in angles.items():
        src, bulk = tmp_path / f"src_{v}", tmp_path / f"bulk_{v}"
        src.mkdir()
        text = (NOCO_INPUTS / "FeAPWlo_noco" / "inp.xml").read_text()
        # FLEUR allows l_mtNocoPot only without symmetry operations: full k-point mesh
        text = re.sub(r"<symmetryOperations>.*?</symmetryOperations>", ident, text, flags=re.S)
        text = re.sub(r"<kPointLists>.*?</kPointLists>", kpts, text, flags=re.S)
        (src / "inp.xml").write_text(text)
        _fleur_run(fleur, src, bulk,
                   {'itmax="2"': 'itmax="3"', 'l_mperp="F" l_mtNocoPot="F"': 'l_mperp="T" l_mtNocoPot="T"',
                    '<nocoParams alpha="0.0" beta="0.0"/>': f'<nocoParams alpha="0.0" {b1}/>',
                    '<nocoParams alpha="0.0" beta="1.57079633"/>': f'<nocoParams alpha="0.0" {b2}/>'})
        _gfprep(gfprep, tmp_path / f"prep_{v}", bulk, "planes_per_layer=2")
        checks = _checks((tmp_path / f"prep_{v}" / "gfprep.out").read_text())
        # the moments turn inside the spheres: finite MT off-diagonal potential
        assert checks[("potential21", 1, "muffin-tins")][0] > 1e-6
        _cbs(gfleur, tmp_path / f"prep_{v}", tmp_path / f"cbs_{v}",
             gf_en="l              1\n   0.30000   0.00000   1.00000   0.00000\n")
    a, b = _cbs_kz(tmp_path / "cbs_a" / "CBS_eigenvalues.1"), _cbs_kz(tmp_path / "cbs_b" / "CBS_eigenvalues.1")
    period = 2*math.pi/5.42

    def dist(k1, k2):
        d = k1 - k2
        return abs(complex((d.real + period/2) % period - period/2, d.imag))
    assert a.keys() == b.keys()
    for e in a:
        assert len(a[e]) == len(b[e]) > 0
        assert max(min(dist(k, kb) for kb in b[e]) for k in a[e]) < 1e-5, f"E={e}"


def _gf_en_dos(e_min=E_MIN, e_max=E_MAX, e_step=E_STEP, e_imag=E_IMAG):
    n = round((e_max - e_min)/e_step) + 1
    return f"l {n:14d}\n" + "".join(f"{e_min + i*e_step:10.5f}{e_imag:10.5f}{1.0:10.5f}{0.0:10.5f}\n"
                                    for i in range(n))


def _film_dos(gfleur, prep, run, collinear=False, window=()):
    """DOS at Gamma of the film embedded in the vacuum on both sides; rows of gf_DOS.raw"""
    _layer_files(prep, run)
    if collinear:
        _replace(run / "layer1_inp.xml", {'l_noco="T"': 'l_noco="F"'})
    gfinp = (prep / "gf_inp").read_text()
    assert "&EMB addemb=T, surface=T, vacleft=T /" in gfinp
    (run / "gf_inp").write_text(gfinp.replace("&MODE tmat=t, gmat=t /",
                                              f"&MODE tmat=t, gmat=t, dos=t /\n&BASIS napw={NAPW} /"))
    (run / "gf_en").write_text(_gf_en_dos(*window))
    (run / "gf_dos").write_text("&LDOS layer=1, atom=0, l=-1 /\n")
    proc, log = _run(gfleur, run)
    assert proc.returncode == 0, f"gfleur DOS run in {run} failed\n{log[-3000:]}"
    return [[float(x) for x in line.split()] for line in (run / "gf_DOS.raw").read_text().splitlines()
            if line.strip()]


def _dos_peaks(rows, ncol):
    """Maxima of the total DOS (spin up minus the negative spin down columns)"""
    e = [r[0] for r in rows]
    y = [sum(r[1:1 + ncol]) - sum(r[1 + ncol:1 + 2*ncol]) for r in rows]
    peaks = []
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] >= y[i + 1] and y[i] > 0.01*max(y):
            a, b, c = y[i - 1], y[i], y[i + 1]
            peaks.append(e[i] + 0.5*(e[1] - e[0])*(a - c)/(a - 2*b + c))
    return peaks


def test_film_dos_noco(tmp_path):
    """Mn monolayer with the moment along x and along z in the vacuum: equal DOS peaks; the
    z film run collinearly gives the noco DOS."""
    gfprep, fleur, gfleur = _binaries()
    dos = {}
    for d, beta in (("x", 'beta="0.5*Pi"'), ("z", 'beta="0.0"')):
        _fleur_run(fleur, NOCO_INPUTS / "MnFilmX", tmp_path / f"film_{d}",
                   {'itmax="5"': 'itmax="40"', 'alpha="0.05"': 'alpha="0.15"', 'beta="0.5*Pi"': beta},
                   files=("inp.xml", "kpts.xml", "sym.xml"))
        _gfprep(gfprep, tmp_path / f"prep_{d}", tmp_path / f"film_{d}", "planes_per_layer=1")
        dos[d] = _film_dos(gfleur, tmp_path / f"prep_{d}", tmp_path / f"dos_{d}")
    out = (tmp_path / "dos_x" / "out").read_text()
    diffs = [float(x) for x in re.findall(r"vacuum embedding side \d:.*step-doubling difference:\s+(\S+)", out)]
    assert len(diffs) == 2 and max(diffs) < 1e-6
    # 4 columns (s,p,d,f) per spin, the third spin is the off-diagonal part
    assert all(len(r) == 13 for r in dos["x"])
    px, pz = _dos_peaks(dos["x"], 4), _dos_peaks(dos["z"], 4)
    assert len(px) == len(pz) >= 8
    assert px == pytest.approx(pz, abs=2e-4)
    # moments along the local axes: no off-diagonal MT DOS in the local frame
    assert max(abs(v) for r in dos["x"] for v in r[9:]) < 1e-6*max(abs(v) for r in dos["x"] for v in r[1:9])

    col = _film_dos(gfleur, tmp_path / "prep_z", tmp_path / "dos_z_collinear", collinear=True)
    assert len(col) == len(dos["z"])
    scale = max(abs(v) for r in col for v in r[1:])
    for rc, rn in zip(col, dos["z"]):
        assert rc[1:9] == pytest.approx(rn[1:9], abs=1e-6*scale)


def test_gfprep_rejects_spin_spirals(tmp_path):
    """gfprep stops with a clear message for a spin-spiral bulk."""
    gfprep, _, _ = _binaries()
    bulk = tmp_path / "bulk"
    bulk.mkdir()
    shutil.copy(NOCO_INPUTS / "1atz" / "inp.xml", bulk / "inp.xml")
    _replace(bulk / "inp.xml", {'l_noco="t"': 'l_noco="t" l_ss="t"'})
    work = tmp_path / "prep"
    work.mkdir()
    (work / "gfprep_inp").write_text("&GFPREP bulkdir='../bulk/' kgrid=1,1 /\n")
    proc, log = _run(gfprep, work)
    assert proc.returncode != 0
    assert "Spin spirals are not supported by gfprep" in log
