"""Tests of the substrate embedding potential for an ad-layer of another 2D lattice
(incomm_inp, gf_embtransform/gf_sigtab).

The substrate file is synthetic: analytic, rotation- and time-reversal-invariant
amplitudes sigma_u(p) on a hexagonal 6x6 mesh, stored at the irreducible k-points
only (so the symmetry unfolding is always exercised). gfleur stops after writing the
transformed potential (stop_after=t), so no converged potentials are needed.
The results are compared with the NumPy reference in tools/incomm_ref.py:

* identity: ad-layer cell = substrate cell, the output equals the input;
* 2x2 supercell with a lateral shift and sqrt3 x sqrt3 R30 (theta=30): exact lookups;
* incommensurate (theta=7): interpolated Sigma0, diagonal;
* vacuum: free-electron amplitudes reproduce -1/2 sqrt(|k+G|^2+2(V-E));
* explicit substrate k-list (write_ksub) with exact lookups;
* Umklapp correction (order 2, 3, full) with a synthetic Sigma_up;
* symmetry check: an ad-layer operation that is no substrate operation is rejected.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("h5py")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "tools"))
import incomm_ref as ref  # noqa: E402
from test_gfleur_smoke import SCHEMA_DIR, _find_gfleur  # noqa: E402

TEMPLATE = HERE / "inputfiles" / "2layer_film" / "layer1_inp.xml"
ASUB = 5.0                     # substrate hexagonal lattice constant
NMESH = (6, 6)                 # substrate k-mesh
RK_SUB = 7.0                   # substrate Kmax
KMAX_AD = 3.0                  # ad-layer Kmax
ENERGIES = [complex(-0.2, 0.05), complex(0.1, 0.08)]
SQUARE_C4 = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
             [[-1, 0, 0], [0, -1, 0], [0, 0, 1]], [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]]


@pytest.fixture(scope="module")
def gfleur_binary():
    binary = _find_gfleur()
    if binary is None:
        pytest.skip("no gfleur binary found (set GFLEUR_BINARY or build the gfleur target)")
    return binary


def write_layer(path, amat, kpts, kmax=KMAX_AD, jspins=2, ops=None):
    """layer1_inp.xml of a one-atom film with the in-plane lattice amat (columns)."""
    xml = TEMPLATE.read_text()
    xml = re.sub(r'Kmax="[^"]*"', f'Kmax="{kmax:.8f}"', xml)
    xml = re.sub(r'jspins="\d"', f'jspins="{jspins}"', xml)
    xml = re.sub(r'radius="[^"]*"', 'radius="2.20000000"', xml)
    kp = "".join(f'            <kPoint weight="{1.0 / len(kpts):.10f}">{k[0]:.12f} {k[1]:.12f} 0.0</kPoint>\n'
                 for k in kpts)
    xml = re.sub(r'<kPointList name="default" count="\d+">.*?</kPointList>',
                 f'<kPointList name="default" count="{len(kpts)}">\n{kp}         </kPointList>', xml, flags=re.S)
    ops = ops or [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]
    so = "".join("         <symOp>\n" + "".join(f"            <row-{i + 1}>{r[0]} {r[1]} {r[2]} .0000000000</row-{i + 1}>\n"
                                                for i, r in enumerate(op)) + "         </symOp>\n" for op in ops)
    xml = re.sub(r"<symmetryOperations>.*?</symmetryOperations>",
                 f"<symmetryOperations>\n{so}      </symmetryOperations>", xml, flags=re.S)
    rows = "".join(f"            <row-{i + 1}>{amat[0, i]:.10f} {amat[1, i]:.10f} .0000000000</row-{i + 1}>\n"
                   for i in range(2))
    xml = re.sub(r"<row-1>[^<]*</row-1>\s*<row-2>[^<]*</row-2>(\s*<row-3>.0000000000 .0000000000)",
                 rows.strip() + r"\1", xml)
    path.write_text(xml)


def stage(work, amat, kpts, incomm, ops=None, jspins=2):
    work.mkdir(parents=True, exist_ok=True)
    write_layer(work / "layer1_inp.xml", amat, kpts, ops=ops, jspins=jspins)
    (work / "gf_inp").write_text("&LAYERS nlayers=1 /\n&MODE charge=f,dos=f,tmat=f,gmat=f /\n")
    lines = [f"l{len(ENERGIES):>17d}"] + [f"{e.real:10.5f}{e.imag:10.5f}{1.0:10.5f}{0.0:10.5f}" for e in ENERGIES]
    (work / "gf_en").write_text("\n".join(lines) + "\n")
    (work / "incomm_inp").write_text("&INCOMM " + ", ".join(f"{k}={v}" for k, v in incomm.items()) + " /\n")
    for xsd in ("FleurInputSchema.xsd", "FleurOutputSchema.xsd"):
        shutil.copy(SCHEMA_DIR / xsd, work / xsd)


def run(binary, work, expect_ok=True):
    proc = subprocess.run([str(binary)], cwd=work, capture_output=True, text=True, timeout=600)
    out = (work / "out").read_text() if (work / "out").is_file() else ""
    text = out + "\n" + proc.stdout + "\n" + proc.stderr
    if expect_ok:
        assert proc.returncode == 0, f"gfleur exited with {proc.returncode}\n{text[-4000:]}"
    return proc.returncode, text


def substrate(work, name="sub", model=ref.sigma_model, spins=2, kpts=None, symmetric=True):
    """synthetic substrate file: irreducible k of the 6x6 mesh, C6v + time reversal"""
    asub = ref.hex_amat(ASUB)
    bsub = ref.recip(asub)
    cart = ref.c6v_cartesian()
    if kpts is None:
        kpts = ref.irreducible(ref.mesh(*NMESH), bsub, cart) if symmetric else ref.mesh(*NMESH)
    ops = ref.point_group_ops(asub, cart) if symmetric else None
    gv = ref.glist(bsub, RK_SUB)
    ref.write_embfile(str(work / name), asub, kpts, gv, ENERGIES, spins, model, RK_SUB,
                      ops=ops, timerev=symmetric)
    return asub, gv


def compare_sigma0(work, asub, gv, theta=0.0, shift=(0.0, 0.0), model=ref.sigma_model, rtol=1e-9,
                   outfile="embpot11"):
    out = ref.read_embfile(str(work / outfile))
    ao = out["attrs"]["amat2d"].T
    geo = ref.Geometry(ao, asub, theta, out["attrs"]["rkmax"][0])
    glist = out["gvecs"]
    nspin = out["mat"].shape[2]
    for ie, e in enumerate(ENERGIES):
        for s in range(nspin):
            tab = ref.SigRef(asub, NMESH, gv, RK_SUB, model, e, s + 1)
            for ik, kap in enumerate(out["kpts"]):
                want = ref.sigma0(geo, tab, kap, glist, shift)
                got = out["mat"][ik, ie, s]
                scale = np.max(np.abs(want))
                assert np.max(np.abs(got - want)) <= rtol * scale, \
                    f"E{ie} spin{s} k{kap}: max dev {np.max(np.abs(got - want)):.3e} (scale {scale:.3e})"
    return out, geo


def test_identity(gfleur_binary, tmp_path):
    """Same cell: the output reproduces the (unfolded) substrate matrices."""
    asub, gv = substrate(tmp_path)
    kpts = ref.mesh(*NMESH)[[0, 1, 8, 15, 21, 35]]
    stage(tmp_path, ref.hex_amat(ASUB), kpts, {"subfile": "'sub'", "stop_after": "t"})
    _, text = run(gfleur_binary, tmp_path)
    assert "Lambda_c" in text
    out, geo = compare_sigma0(tmp_path, asub, gv, rtol=1e-11)
    # direct comparison with the analytic matrix of the substrate
    bsub = ref.recip(asub)
    f = ref.sigma_model(ENERGIES[0], 1)
    for ik, k in enumerate(out["kpts"]):
        hl = geo.home(k, out["gvecs"])
        g = out["gvecs"][hl]
        p = k + g[None, :, :]
        u = g[:, None, :] - g[None, :, :]
        pc, puc, uc = p @ bsub, (p + u) @ bsub, u @ bsub
        want = f(np.sum(pc ** 2, -1), np.sum(puc ** 2, -1), np.sum(uc ** 2, -1), np.all(u == 0, -1))
        got = out["mat"][ik, 0, 0][np.ix_(hl, hl)]
        assert np.allclose(got, want, rtol=0, atol=1e-12 * np.max(np.abs(want)))


def test_supercell_shift(gfleur_binary, tmp_path):
    """2x2 supercell with a lateral shift: Lambda_c terms with phases exp(-i u.d)."""
    asub, gv = substrate(tmp_path)
    kpts = ref.mesh(3, 3)
    shift = (0.7, -0.3)
    stage(tmp_path, ref.hex_amat(2 * ASUB), kpts,
          {"subfile": "'sub'", "stop_after": "t", "shift": f"{shift[0]}, {shift[1]}"})
    run(gfleur_binary, tmp_path)
    out, _ = compare_sigma0(tmp_path, asub, gv, shift=shift, rtol=1e-11)
    assert np.count_nonzero(np.abs(out["mat"][0, 0, 0]) > 0) > np.count_nonzero(
        np.abs(np.diag(out["mat"][0, 0, 0])) > 0), "supercell: no off-diagonal Lambda_c couplings"


def test_root3_rotated(gfleur_binary, tmp_path):
    """sqrt3 x sqrt3 R30: the substrate rotated by theta=30 is commensurate."""
    asub, gv = substrate(tmp_path)
    kpts = ref.mesh(2, 2)
    stage(tmp_path, ref.hex_amat(np.sqrt(3.0) * ASUB), kpts, {"subfile": "'sub'", "stop_after": "t", "theta": 30.0})
    run(gfleur_binary, tmp_path)
    compare_sigma0(tmp_path, asub, gv, theta=30.0, rtol=1e-10)


def test_incommensurate(gfleur_binary, tmp_path):
    """Rotated, stretched ad-layer: diagonal Sigma0 from the interpolated amplitudes."""
    asub, gv = substrate(tmp_path)
    kpts = ref.mesh(3, 3)
    stage(tmp_path, ref.hex_amat(5.3), kpts, {"subfile": "'sub'", "stop_after": "t", "theta": 7.0})
    _, text = run(gfleur_binary, tmp_path)
    assert "incommensurate" in text
    out, _ = compare_sigma0(tmp_path, asub, gv, theta=7.0, rtol=1e-9)
    m = out["mat"][0, 0, 0]
    assert np.allclose(m - np.diag(np.diag(m)), 0.0)


def test_vacuum(gfleur_binary, tmp_path):
    """Free-electron substrate: Sigma0 = -1/2 sqrt(|kappa+G|^2 + 2(V-E))."""
    asub, gv = substrate(tmp_path, model=ref.sigma_vacuum, spins=1)
    kpts = ref.mesh(3, 3)
    stage(tmp_path, ref.hex_amat(4.6, 13.0), kpts, {"subfile": "'sub'", "stop_after": "t", "theta": 3.0})
    run(gfleur_binary, tmp_path)
    out = ref.read_embfile(str(tmp_path / "embpot11"))
    bo = ref.recip(out["attrs"]["amat2d"].T)
    geo = ref.Geometry(out["attrs"]["amat2d"].T, asub, 3.0, out["attrs"]["rkmax"][0])
    for ie, e in enumerate(ENERGIES):
        for ik, kap in enumerate(out["kpts"]):
            hl = geo.home(kap, out["gvecs"])
            q = (kap + out["gvecs"][hl]) @ bo
            want = -0.5 * np.sqrt(np.sum(q ** 2, 1) + 2 * (0.2 - e))
            got = np.diag(out["mat"][ik, ie, 0])[hl]
            assert np.max(np.abs(got - want) / np.abs(want)) < 2e-3
    compare_sigma0(tmp_path, asub, gv, theta=3.0, model=ref.sigma_vacuum, rtol=1e-9)


def test_exact_klist(gfleur_binary, tmp_path):
    """write_ksub writes the substrate k-list; a substrate file on that list gives
    Sigma0 by exact lookups (list mode, no interpolation)."""
    kpts = ref.mesh(2, 2)
    stage(tmp_path, ref.hex_amat(5.3), kpts, {"subfile": "'sub'", "theta": 7.0, "write_ksub": "t"})
    asub, gv = substrate(tmp_path)  # only needed for the geometry of the first run
    run(gfleur_binary, tmp_path)
    ksub = ref.read_klist(tmp_path / "kpts_sub")
    assert len(ksub) > len(kpts)
    asub, gv = substrate(tmp_path, kpts=ksub, symmetric=False)
    (tmp_path / "incomm_inp").write_text("&INCOMM subfile='sub', theta=7.0, interp='exact', stop_after=t /\n")
    _, text = run(gfleur_binary, tmp_path)
    assert "list mode" in text
    out = ref.read_embfile(str(tmp_path / "embpot11"))
    geo = ref.Geometry(out["attrs"]["amat2d"].T, asub, 7.0, out["attrs"]["rkmax"][0])
    bsub = ref.recip(asub)
    for ie, e in enumerate(ENERGIES):
        f = ref.sigma_model(e, 1)
        for ik, kap in enumerate(out["kpts"]):
            hl = geo.home(kap, out["gvecs"])
            pc = ((kap + out["gvecs"][hl]) @ geo.m) @ bsub
            p2 = np.sum(pc ** 2, 1)
            want = f(p2, p2, 0.0, True)
            got = np.diag(out["mat"][ik, ie, 0])[hl]
            assert np.allclose(got, want, rtol=1e-10, atol=0)


@pytest.mark.parametrize("order", ["2", "3", "full"])
def test_umklapp(gfleur_binary, tmp_path, order):
    """Umklapp correction with a synthetic Sigma_up at the auxiliary k-points."""
    asub, gv = substrate(tmp_path)
    kpts = ref.mesh(2, 2)
    ao = ref.hex_amat(5.3)
    base = {"subfile": "'sub'", "theta": 7.0, "nshell": 1}
    stage(tmp_path, ao, kpts, dict(base, stop_after="t"))
    run(gfleur_binary, tmp_path)
    glist = ref.read_embfile(str(tmp_path / "embpot11"))["gvecs"]
    (tmp_path / "incomm_inp").write_text("&INCOMM subfile='sub', theta=7.0, nshell=1, write_kaux=t /\n")
    run(gfleur_binary, tmp_path)
    kaux = ref.read_klist(tmp_path / "kpts_aux")
    assert len(kaux) > 0 and (tmp_path / "kaux_map").is_file()
    rk = KMAX_AD
    ref.write_embfile(str(tmp_path / "up"), ao, kaux, glist, ENERGIES, 2, ref.up_model, rk, model_kind="up")
    (tmp_path / "incomm_inp").write_text(
        f"&INCOMM subfile='sub', theta=7.0, nshell=1, order='{order}', upfile='up', stop_after=t /\n")
    _, text = run(gfleur_binary, tmp_path)
    assert "Umklapp correction" in text and "||Delta_u0||/||Sigma0||" in text
    out = ref.read_embfile(str(tmp_path / "embpot11"))
    up = ref.read_embfile(str(tmp_path / "up"))
    geo = ref.Geometry(ao, asub, 7.0, out["attrs"]["rkmax"][0])
    q = ref.shells(geo, 1)
    assert len(q) == 6
    for ie, e in enumerate(ENERGIES):
        for s in range(2):
            tab = ref.SigRef(asub, NMESH, gv, RK_SUB, ref.sigma_model, e, s + 1)
            for ik, kap in enumerate(out["kpts"]):
                hl, want = ref.sigma_eff(geo, tab, kap, glist, q, up, ie, s, {"2": 2, "3": 3, "full": 99}[order])
                got = out["mat"][ik, ie, s][np.ix_(hl, hl)]
                s0 = ref.sigma0(geo, tab, kap, glist)[np.ix_(hl, hl)]
                assert np.max(np.abs(want - s0)) > 1e-6, "the correction vanishes"
                assert np.max(np.abs(got - want)) <= 1e-9 * np.max(np.abs(want))


def test_symmetry_check(gfleur_binary, tmp_path):
    """A square ad-layer with C4 is rejected on the hexagonal substrate."""
    substrate(tmp_path)
    a = 5.2
    stage(tmp_path, np.array([[a, 0.0], [0.0, a]]), ref.mesh(2, 2), {"subfile": "'sub'", "stop_after": "t"},
          ops=SQUARE_C4)
    code, text = run(gfleur_binary, tmp_path, expect_ok=False)
    assert code != 0
    assert "no symmetry of the substrate" in text


# --------------------------------------------------------------------------- writeup (real film)
from test_gfprep_film import _run as _run_bin, film  # noqa: E402,F401  (module fixture: converged Cu film)

FILM_EIGS_WINDOW = (-0.5, -0.3)
FILM_PEAK_TOL = 0.0025


def _writeup_run(gfprep, gfleur, base, planes_per_layer):
    work = base / f"writeup_{planes_per_layer}"
    work.mkdir()
    (work / "gfprep_inp").write_text(f"&GFPREP bulkdir='../film/' kgrid=1,1 planes_per_layer={planes_per_layer} /\n")
    proc, log = _run_bin(gfprep, work)
    assert proc.returncode == 0, f"gfprep failed\n{log[-3000:]}"
    gfinp = (work / "gf_inp").read_text()
    nl = int(re.search(r"nlayers=(\d+)", gfinp).group(1))
    gfinp = gfinp.replace("&MODE tmat=t, gmat=t /", "&MODE tmat=t, gmat=t /\n&BASIS napw=8 /")
    gfinp = re.sub(r"(&EMB[^/]*)/", r"\1, writeup=T /", gfinp)
    (work / "gf_inp").write_text(gfinp)
    n = 101
    (work / "gf_en").write_text(f"l {n:14d}\n" + "".join(
        f"{FILM_EIGS_WINDOW[0] + i * 0.002:10.5f}{0.002:10.5f}{1.0:10.5f}{0.0:10.5f}\n" for i in range(n)))
    (work / "gf_io").write_text(
        "3\n&IO mode='E', region=0, side=2, rw='w', id=1 /\n&IO mode='E', region=1, side=1, rw='w', id=2 /\n"
        f"&IO mode='E', region={nl + 1}, side=1, rw='w', id=3 /\n"
        "&FILE name='sigup' /\n&FILE name='sigleft' /\n&FILE name='sigtop' /\n")
    proc, log = _run_bin(gfleur, work)
    assert proc.returncode == 0, f"gfleur failed\n{log[-3000:]}"
    return [ref.read_embfile(str(work / name)) for name in ("sigup", "sigleft", "sigtop")]


@pytest.mark.parametrize("planes_per_layer", [3, 1])
def test_writeup_film(film, planes_per_layer):
    """Cu(001) film, one layer (3 planes) or three layers: Sigma_up at the plane below
    layer 1 (writeup=t, region 0) gives G_S = (Sigma_left + Sigma_up)^-1 with peaks at the
    eigenvalues of the film (checks G_S^-1 = Sigma_L + Sigma_R); for one layer it equals the
    potential at the top plane seen from below (mirror-symmetric film)."""
    gfprep, gfleur, base, eig = film
    up, left, top = _writeup_run(gfprep, gfleur, base, planes_per_layer)
    energies = up["energies"]
    sup, sl, stop = up["mat"][0, :, 0], left["mat"][0, :, 0], top["mat"][0, :, 0]
    basis = np.where(np.abs(np.diag(sup[0])) > 0)[0]
    dos = np.array([-np.trace(np.linalg.inv((sl[i] + sup[i])[np.ix_(basis, basis)])).imag / np.pi
                    for i in range(len(energies))])
    assert dos.min() > 0.0, "G_S = (Sigma_L + Sigma_up)^-1 is not retarded"
    peaks = [energies[i].real for i in range(1, len(dos) - 1) if dos[i] > dos[i - 1] and dos[i] > dos[i + 1]]
    want = [e for e in eig if FILM_EIGS_WINDOW[0] + 0.01 < e < FILM_EIGS_WINDOW[1] - 0.01]
    assert want, "no film eigenvalue in the energy window"
    for e in want:
        assert min(abs(p - e) for p in peaks) < FILM_PEAK_TOL, f"no peak at the eigenvalue {e}: {peaks}"
    # reciprocity of the non-magnetic Sigma at Gamma, mirror symmetry of the film
    scale = np.abs(sup).max()
    assert np.abs(sup - np.swapaxes(sup, 1, 2)).max() < 1e-6 * scale
    if planes_per_layer == 3:
        # both from the T-matrix of the single layer: exact mirror symmetry
        # (with several layers the composition limits the agreement to ~1e-2)
        assert np.abs(stop - sup).max() < 1e-6 * scale
    # metadata: the planes lie 2 a.u. from the outer atom planes
    assert abs(up["attrs"]["zplane"][0] - top["attrs"]["zplane"][0]) < 1e-8
    assert up["attrs"]["zplane"][0] > 0.0


def _film_sigup(gfprep, gfleur, base, name, kopts):
    """Sigma_up below the film (one layer) on a 4x4 mesh"""
    work = base / name
    work.mkdir()
    (work / "gfprep_inp").write_text(f"&GFPREP bulkdir='../film/' kgrid=4,4 planes_per_layer=3 {kopts} /\n")
    proc, log = _run_bin(gfprep, work)
    assert proc.returncode == 0, f"gfprep failed\n{log[-3000:]}"
    gfinp = (work / "gf_inp").read_text().replace("&MODE tmat=t, gmat=t /", "&MODE tmat=t, gmat=t /\n&BASIS napw=8 /")
    (work / "gf_inp").write_text(re.sub(r"(&EMB[^/]*)/", r"\1, writeup=T /", gfinp))
    energies = [(-0.45, 0.01), (-0.3, 0.02), (-0.1, 0.03)]
    (work / "gf_en").write_text(f"l {len(energies):14d}\n" + "".join(
        f"{a:10.5f}{b:10.5f}{1.0:10.5f}{0.0:10.5f}\n" for a, b in energies))
    (work / "gf_io").write_text("1\n&IO mode='E', region=0, side=2, rw='w', id=1 /\n&FILE name='sigup' /\n")
    proc, log = _run_bin(gfleur, work)
    assert proc.returncode == 0, f"gfleur failed\n{log[-3000:]}"
    return work


def test_unfold_film(film):
    """Sigma_up of the Cu(001) film from the symmetry-reduced 4x4 mesh (C4v and time
    reversal), unfolded by incomm_inp onto the full mesh of the same cell, equals a
    direct run on the full mesh (gfprep kfull=t)."""
    gfprep, gfleur, base, _ = film
    red = _film_sigup(gfprep, gfleur, base, "unfold_red", "")
    full = _film_sigup(gfprep, gfleur, base, "unfold_full", "kfull=t")
    assert len(ref.read_embfile(str(red / "sigup"))["kpts"]) < 16
    shutil.copy(red / "sigup.hdf", full / "subred.hdf")
    (full / "gf_io").unlink()
    (full / "incomm_inp").write_text("&INCOMM subfile='subred', outfile='unfolded', stop_after=t /\n")
    proc, log = _run_bin(gfleur, full)
    assert proc.returncode == 0, f"incomm transform failed\n{log[-3000:]}"
    out = (full / "out").read_text()
    dev = float(re.search(r"max deviation of symmetry images\s+(\S+)", out).group(1))
    assert dev < 1e-8
    a = ref.read_embfile(str(full / "unfolded"))
    b = ref.read_embfile(str(full / "sigup"))
    assert len(b["kpts"]) == 16 and np.allclose(a["kpts"], b["kpts"])
    assert np.abs(a["mat"] - b["mat"]).max() < 1e-6 * np.abs(b["mat"]).max()
