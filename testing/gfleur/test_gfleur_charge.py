"""Self-consistency (charge mode) of gfleur on a perfect Cu crystal.

The 4-atom tetragonal Cu supercell of test_gfleur_transport (ABAB stacking, RMT 1.6) is
converged with fleur. gfprep cuts one two-plane period as the lead; a gfleur CBS run on the
energy contour of gfprep (gf_en.charge: semicircle from E_F-1 Htr to E_F) writes the lead
embedding potentials and the Coulomb boundary values. A stack of two layers of two atomic planes
between these leads is then iterated with `&MODE charge=t`. The leads and the stack are the same
crystal, so the self-consistent result must reproduce the bulk.

Convergence of the first iteration, tested with the fleur bulk density as input (2x2 mesh: the
total charge is 1.2% off, 4x4 and larger: 0.15%, converged in the contour (15, 31, 63 points
identical) and in the cylinder basis (napw 8, 12, 16 identical to 1e-4 e)):

    kgrid            2x2     4x4     6x6     8x8
    MT d charge     8.392   8.246   8.257   8.261   (fleur bulk: 8.180)
    neutrality     -1.370  -0.163  -0.128  -0.232   (4 atoms)

The self-consistent fixed point of the stack approaches the fleur values. The Hartree response
of the stack has a charge-sloshing mode (eigenvalue of the order of 50): plain mixing needs
alpha below 0.04, Broyden mixing is used with alpha=0.02.

Surface (frozen vacuum): the upper four planes of a 7-plane Cu(001) film (testing/gfleur/inputfiles/
cu_film7, planes 4-7, one layer per plane) between a bulk lead below (planes 2-3 of the supercell, its
potential shifted to the Fermi energy of the film) and the vacuum above. The vacuum beyond the last
embedding plane is not updated: gfprep stores its planar Coulomb field at the plane, which is the Neumann
condition of the Coulomb problem for the planar component (G_par=0); the other components decay,
V' = -|G_par| V. The vacuum holds 0.205 e of the electrons of the film (the last layer has 28.795 of 29),
so the stack converges to a charge of +0.205 e. The layers must not be mirror symmetric in z (`zmirror=f`
in gfprep), the Coulomb solver cannot handle planar stars with kz>=0 only. Plain linear mixing is stable
(alpha=0.02); Broyden (imix=7) is not, imix=5 converges in ~25 iterations (distance 7.6 -> 0.06).

Noco (sc Fe, one atom per period, moment along z or x, the semicore 3s3p states need a contour
that starts 4 Htr below E_F): the noco charge mode reproduces the collinear first iteration for the
moment along z exactly (same stack charge, same MT charge, same distance). For the moment along x
the MT d charge differs by 0.37 e and the stack charge by 0.34 e: the first iteration is not yet
rotation invariant (known problem, xfail).

MPI: with GFLEUR_BINARY=<path of gfleur_MPI> and GFLEUR_MPI="mpirun -np 2" the gfleur runs of these
tests use several ranks (the k-points are distributed, the density is summed over the ranks of a
layer). The results must equal the serial ones.

Binaries are located via GFPREP_BINARY or <repo>/build*/gfprep; fleur and gfleur are taken from
the same directory. Missing binaries skip the tests.
"""
import os
import re
import shutil
import subprocess
from collections import namedtuple
from pathlib import Path

import pytest

from test_gfleur_transport import (BULK_INPUTS, MT_SETUPS, _find_gfprep, _layer_files, _replace, _run,
                                   _sibling)

FILM7 = Path(__file__).resolve().parent / "inputfiles" / "cu_film7"
KGRID = "kgrid=4,4"
NAPW = 8
SCF_ITER = 14
SCF_MIX = "imix=7, maxiter=14, alpha=0.02"

Env = namedtuple("Env", "gfprep gfleur bulk base prep lead cbs")


def _run_gf(gfleur, cwd, timeout=3600, mpi=True):
    """gfleur, optionally under a launcher (GFLEUR_MPI="mpirun -np 2"); mpi=False for the serial
    preparation runs (CBS of the leads)"""
    cmd = (os.environ.get("GFLEUR_MPI", "").split() if mpi else []) + [str(gfleur)]
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=dict(os.environ))
    return proc, proc.stdout + "\n" + proc.stderr


def _gfprep(gfprep, work, bulk, options):
    work.mkdir()
    # relative path: the xml file name handling of FLEUR is limited to 100 characters
    (work / "gfprep_inp").write_text(f"&GFPREP bulkdir='../{bulk.name}/' {options} {KGRID} /\n")
    proc, log = _run(gfprep, work)
    assert proc.returncode == 0, f"gfprep failed\n{log[-3000:]}"
    return log


def _layers_block(prep):
    text = (prep / "gf_inp").read_text()
    return text[:text.index("\n/\n") + 3]


def _gfleur_charge(env, name, sc, copy_cdn=True, mpi=True):
    """gfleur charge run of the stack with the given &SC line, returns (run dir, log)"""
    run = env.base / name
    run.mkdir()
    _layer_files(env.prep, run)
    for f in ("embpot11.hdf", "embpotN2.hdf"):
        shutil.copy(env.cbs / f, run / f)
    shutil.copy(env.lead / "gf_en.charge", run / "gf_en")
    (run / "gf_inp").write_text(_layers_block(env.prep) + "&MODE tmat=t, gmat=t, charge=t /\n"
                                f"&BASIS napw={NAPW} /\n&EMB addemb=t /\n{sc}\n")
    nl = len(list(env.prep.glob("layer*_inp.xml")))
    (run / "gf_io").write_text("2\n&IO mode='E', region=1, side=1, rw='r', id=1 /\n"
                               f"&IO mode='E', region={nl}, side=2, rw='r', id=2 /\n"
                               "&FILE name='embpot11' /\n&FILE name='embpotN2' /\n")
    if not copy_cdn:
        (run / "gf_cdn.hdf").unlink()
    proc, log = _run_gf(env.gfleur, run, mpi=mpi)
    return run, proc, log


def _neutrality(out):
    """total electron charge minus nuclear charge per iteration"""
    return [float(a) + float(b) for a, b in re.findall(r"Total charge neutrality:\s*(\S+)\+\s*(\S+)=", out)]


def _distances(out):
    return [float(x) for x in re.findall(r"Mix: total distance =\s*(\S+)", out)]


def _mt_charges(out):
    """(s,p,d,f,total) of the first atom of layer 1 for each iteration"""
    return [tuple(float(x) for x in m) for m in
            re.findall(r"\n\s+1 --> 1\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", out)]


def _bulk_mt_charge(bulk):
    """total MT valence charge of the first atom type of the fleur run (last iteration)"""
    m = re.findall(r'<mtCharge atomType="1" total="(\S+)" s="(\S+)" p="(\S+)" d="(\S+)" f="(\S+)"/>',
                   (bulk / "out.xml").read_text())[-1]
    return tuple(float(x) for x in m[1:]) + (float(m[0]),)


@pytest.fixture(scope="module")
def env(tmp_path_factory):
    """Bulk run, lead CBS on the contour of gfprep, stack for the charge runs."""
    gfprep = _find_gfprep()
    if gfprep is None:
        pytest.skip("no gfprep binary found (set GFPREP_BINARY or build the gfprep target)")
    fleur = _sibling(gfprep, ["fleur", "fleur_MPI"])
    gfleur = _sibling(gfprep, ["gfleur", "gfleur_MPI"])
    if os.environ.get("GFLEUR_BINARY"):
        gfleur = Path(os.environ["GFLEUR_BINARY"])
    if fleur is None or gfleur is None:
        pytest.skip(f"fleur/gfleur binaries missing next to {gfprep}")
    base = tmp_path_factory.mktemp("charge")
    bulk = base / "bulk"
    shutil.copytree(BULK_INPUTS, bulk)
    _replace(bulk / "inp.xml", MT_SETUPS["rmt1.6"])
    proc, log = _run(fleur, bulk)
    assert proc.returncode == 0 and (bulk / "cdn.hdf").is_file(), f"bulk fleur run failed\n{log[-3000:]}"

    lead = base / "lead_prep"
    _gfprep(gfprep, lead, bulk, "planes_per_layer=2 first_plane=1 last_plane=2")
    cbs = base / "cbs"
    cbs.mkdir()
    _layer_files(lead, cbs)
    (cbs / "gf_inp").write_text(_layers_block(lead) + f"&MODE tmat=t, gmat=t /\n&BASIS napw={NAPW} /\n"
                                "&CBSMODE cbs=t /\n&SC iter=1 /\n")
    shutil.copy(lead / "gf_en.charge", cbs / "gf_en")
    (cbs / "gf_io").write_text("default\n")
    proc, log = _run_gf(gfleur, cbs, mpi=False)
    assert proc.returncode == 0, f"gfleur CBS run failed\n{log[-3000:]}"
    for f in ("embpot11.hdf", "embpotN2.hdf"):
        assert (cbs / f).is_file(), f"{f} missing"

    prep = base / "stack_prep"
    _gfprep(gfprep, prep, bulk, "planes_per_layer=2")
    return Env(gfprep, gfleur, bulk, base, prep, lead, cbs)


def test_gfprep_contour_template(env):
    """gf_en.charge: semicircle of 15 points ending at the Fermi energy of the bulk"""
    lines = (env.lead / "gf_en.charge").read_text().splitlines()
    assert lines[0].split() == ["c", "15"]
    ef = float(re.findall(r'<FermiEnergy value="(\S+)"', (env.bulk / "out.xml").read_text())[-1])
    assert float(lines[2].split()[0]) == pytest.approx(ef, abs=1e-4)
    assert float(lines[1].split()[0]) == pytest.approx(ef - 1.0, abs=1e-4)


def test_charge_first_iteration(env):
    """One iteration from the fleur density: neutral to 0.3%, MT charges within 0.1 e of the bulk,
    the potential generated from the density reproduces the one of the bulk."""
    run, proc, log = _gfleur_charge(env, "first", "&SC iter=1, imix=0, alpha=0.0 /")
    assert proc.returncode == 0, f"gfleur charge run failed\n{log[-3000:]}"
    out = (run / "out").read_text()
    assert abs(_neutrality(out)[0]) < 0.3
    s, p, d, f, tot = _mt_charges(out)[0]
    sb, pb, db, fb, totb = _bulk_mt_charge(env.bulk)
    assert tot == pytest.approx(totb, abs=0.1)
    assert d == pytest.approx(db, abs=0.1)
    assert _distances(out)[0] < 10.0
    # gf_vgen replaces gf_pot.hdf; the original is kept as gf_pot1.hdf
    # observed: MT (pot/vr) up to 3e-3, interstitial (pot/pw) up to 0.018 in the planar components,
    # which follow the 0.15% error of the density
    for ds, delta in (("pot/vr", "1e-2"), ("pot/pw", "3e-2")):
        proc = subprocess.run(["h5diff", "-d", delta, str(run / "gf_pot1.hdf"), str(run / "gf_pot.hdf"),
                               f"/layer-1/{ds}"], capture_output=True, text=True)
        assert proc.returncode == 0, f"potential of layer 1 ({ds}) differs from the bulk:\n{proc.stdout[-500:]}"


def test_charge_zero_mixing_is_stationary(env):
    """alpha=0: every iteration reproduces the first one"""
    run, proc, log = _gfleur_charge(env, "stationary", "&SC iter=3, imix=0, alpha=0.0 /")
    assert proc.returncode == 0, f"gfleur charge run failed\n{log[-3000:]}"
    out = (run / "out").read_text()
    q, dist = _neutrality(out), _distances(out)
    assert len(q) == 3
    assert q == pytest.approx([q[0]] * 3, abs=1e-8)
    assert dist == pytest.approx([dist[0]] * 3, rel=1e-8)


def test_charge_converges(env):
    """Broyden mixing converges the stack towards the bulk: the distance falls, the stack is neutral
    and the MT charges approach those of the fleur bulk."""
    run, proc, log = _gfleur_charge(env, "scf", f"&SC iter={SCF_ITER}, {SCF_MIX} /")
    assert proc.returncode == 0, f"gfleur charge run failed\n{log[-3000:]}"
    out = (run / "out").read_text()
    dist, q = _distances(out), _neutrality(out)
    assert len(dist) == SCF_ITER
    assert dist[-1] < dist[0]/10
    assert abs(q[-1]) < 0.02
    d = _mt_charges(out)[-1][2]
    assert d == pytest.approx(_bulk_mt_charge(env.bulk)[2], abs=0.03)


def test_charge_needs_cdn(env):
    """Self-consistency needs the density file."""
    run, proc, log = _gfleur_charge(env, "nocdn", "&SC iter=2, imix=0, alpha=0.02 /", copy_cdn=False)
    assert proc.returncode != 0
    assert "No self-consistency without gf_cdn.hdf" in log


SURF_ITER = 26
SURF_MIX = "imix=5, maxiter=26, alpha=0.02"
# layer 1-4 of the stack are the film planes 4-7: atom types 2, 1, 3, 4 of cu_film7
FILM7_TYPES = (2, 1, 3, 4)


def _film7_mt_charges(film):
    """(s,p,d,f,total) of the atom types of the 7-plane film (last iteration)"""
    blocks = re.findall(r'<mtCharge atomType="(\d)" total="(\S+)" s="(\S+)" p="(\S+)" d="(\S+)" f="(\S+)"/>',
                        (film / "out.xml").read_text())[-4:]
    return {int(b[0]): tuple(float(x) for x in b[2:]) + (float(b[1]),) for b in blocks}


def _fermi(out_xml):
    return float(re.findall(r'<FermiEnergy value="(\S+)"', out_xml.read_text())[-1])


@pytest.fixture(scope="module")
def surf(env):
    """7-plane film, lead from the bulk (shifted to the film), stack of the planes 4-7 on the lead."""
    fleur = _sibling(env.gfprep, ["fleur", "fleur_MPI"])
    film = env.base / "film7"
    shutil.copytree(FILM7, film)
    proc, log = _run(fleur, film)
    assert proc.returncode == 0 and (film / "cdn.hdf").is_file(), f"film fleur run failed\n{log[-3000:]}"
    vshift = _fermi(film / "out.xml") - _fermi(env.bulk / "out.xml")

    # the stack first: its contour is also used for the lead
    prep = env.base / "surf_prep"
    _gfprep(env.gfprep, prep, film, "planes_per_layer=1 first_plane=4 zmirror=f")
    lead = env.base / "surf_lead_prep"
    _gfprep(env.gfprep, lead, env.bulk, f"planes_per_layer=2 first_plane=2 last_plane=3 vshift={vshift:.6f}")
    cbs = env.base / "surf_cbs"
    cbs.mkdir()
    _layer_files(lead, cbs)
    (cbs / "gf_inp").write_text(_layers_block(lead) + f"&MODE tmat=t, gmat=t /\n&BASIS napw={NAPW} /\n"
                                "&CBSMODE cbs=t /\n&SC iter=1 /\n")
    shutil.copy(prep / "gf_en.charge", cbs / "gf_en")
    (cbs / "gf_io").write_text("default\n")
    proc, log = _run_gf(env.gfleur, cbs, mpi=False)
    assert proc.returncode == 0, f"gfleur CBS run of the lead failed\n{log[-3000:]}"
    return film, prep, cbs


def _surface_run(env, surf, name, sc):
    film, prep, cbs = surf
    run = env.base / name
    run.mkdir()
    _layer_files(prep, run)
    shutil.copy(cbs / "embpot11.hdf", run / "embpot11.hdf")
    shutil.copy(prep / "gf_en.charge", run / "gf_en")
    (run / "gf_io").write_text("default\n")
    gfinp = (prep / "gf_inp").read_text()
    assert "&EMB addemb=T, surface=T, vacleft=F /" in gfinp
    (run / "gf_inp").write_text(_layers_block(prep) + "&MODE tmat=t, gmat=t, charge=t /\n"
                                f"&BASIS napw={NAPW} /\n&EMB addemb=T, surface=T, vacleft=F /\n{sc}\n")
    proc, log = _run_gf(env.gfleur, run)
    return run, proc, log


def test_surface_gfprep_outputs(env, surf):
    """The layers of the stack have no z mirror (planar stars of both signs); the Coulomb field of the
    vacuum is stored (positive: the potential rises into the vacuum)."""
    film, prep, cbs = surf
    out = (prep / "gfprep.out").read_text()
    assert "vacuum embedding potential side 2" in out
    proc = subprocess.run(["h5dump", "-a", "/vacemb_side2/vcfield", str(prep / "gf_pot.hdf")],
                          capture_output=True, text=True)
    field = float(re.findall(r"\(0\):\s*(\S+)", proc.stdout)[0])
    assert 0.05 < field < 0.2
    for i in range(1, 5):
        proc = subprocess.run(["h5dump", "-d", f"/layer-{i}/stars/kv3", "-y", "-O", str(prep / "gf_pot.hdf")],
                              capture_output=True, text=True)
        kv = [int(float(x)) for x in re.findall(r"-?\d+", proc.stdout)]
        kz = [kv[j + 2] for j in range(0, len(kv) - 2, 3) if kv[j] == 0 and kv[j + 1] == 0]
        assert min(kz) < 0, f"layer {i} is mirror symmetric"


def test_surface_first_iteration(env, surf):
    """One iteration from the film density: the Coulomb potential between the lead and the vacuum is
    smooth (no spurious net charge), MT charges within 0.2 e of the film, stack charge about +0.2 e."""
    film, prep, cbs = surf
    run, proc, log = _surface_run(env, surf, "surface1", "&SC iter=1, imix=0, alpha=0.0 /")
    assert proc.returncode == 0, f"gfleur surface run failed\n{log[-3000:]}"
    out = (run / "out").read_text()
    assert 0.0 < _neutrality(out)[0] < 0.8
    assert _distances(out)[0] < 12.0
    ref = _film7_mt_charges(film)
    for (s, p, d, f, tot), t in zip(re.findall(r"\n\s+\d --> 1\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", out)[-4:], FILM7_TYPES):
        assert float(tot) == pytest.approx(ref[t][4], abs=0.2)


def test_surface_converges(env, surf):
    """SCF of the surface: the distance falls by two orders of magnitude, the stack converges to the
    charge with the vacuum electrons missing (+0.205 e), the MT charges of all layers equal those of
    the film planes."""
    film, prep, cbs = surf
    run, proc, log = _surface_run(env, surf, "surface_scf", f"&SC iter={SURF_ITER}, {SURF_MIX} /")
    assert proc.returncode == 0, f"gfleur surface run failed\n{log[-3000:]}"
    out = (run / "out").read_text()
    dist, q = _distances(out), _neutrality(out)
    assert len(dist) == SURF_ITER
    assert dist[-1] < dist[0]/50
    assert q[-1] == pytest.approx(0.205, abs=0.04)
    ref = _film7_mt_charges(film)
    layers = re.findall(r"\n\s+\d --> 1\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", out)[-4:]
    for (s, p, d, f, tot), t in zip(layers, FILM7_TYPES):
        assert float(d) == pytest.approx(ref[t][2], abs=0.03)
        assert float(tot) == pytest.approx(ref[t][4], abs=0.04)


def test_charge_rejects_vacuum_on_the_left(env, surf):
    """A vacuum on the left (free film) is not supported in the charge mode."""
    film, prep, cbs = surf
    run = env.base / "vacleft"
    run.mkdir()
    _layer_files(prep, run)
    shutil.copy(prep / "gf_en.charge", run / "gf_en")
    (run / "gf_inp").write_text(_layers_block(prep) + "&MODE tmat=t, gmat=t, charge=t /\n"
                                f"&BASIS napw={NAPW} /\n&EMB addemb=T, surface=T, vacleft=T /\n&SC iter=1 /\n")
    proc, log = _run_gf(env.gfleur, run)
    assert proc.returncode != 0
    assert "charge mode with vacuum on the left is not supported" in log


NOCO_INPUTS = Path(__file__).resolve().parent.parent / "inputfiles" / "noco"
NOCO_PREP = "planes_per_layer=1 ewidth=4.0 nen=31"


def _noco_stack(env, d):
    """sc Fe with the moment along d: bulk, gfprep (one layer), lead CBS on the contour; returns a function
    running the stack (one layer between its own leads) in charge mode"""
    base = env.base / f"noco_{d}"
    base.mkdir()
    fleur = _sibling(env.gfprep, ["fleur", "fleur_MPI"])
    bulk = base / "bulk"
    shutil.copytree(NOCO_INPUTS / f"1at{d}", bulk)
    _replace(bulk / "inp.xml", {'itmax="2"': 'itmax="25"', 'alpha=".05000000"': 'alpha=".10000000"'})
    proc, log = _run(fleur, bulk)
    assert proc.returncode == 0 and (bulk / "cdn.hdf").is_file(), f"bulk fleur run failed\n{log[-3000:]}"
    prep = base / "prep"
    _gfprep(env.gfprep, prep, bulk, NOCO_PREP)
    cbs = base / "cbs"
    cbs.mkdir()
    _layer_files(prep, cbs)
    (cbs / "gf_inp").write_text(_layers_block(prep) + f"&MODE tmat=t, gmat=t /\n&BASIS napw={NAPW} /\n"
                                "&CBSMODE cbs=t /\n&SC iter=1 /\n")
    shutil.copy(prep / "gf_en.charge", cbs / "gf_en")
    (cbs / "gf_io").write_text("default\n")
    proc, log = _run_gf(env.gfleur, cbs, mpi=False)
    assert proc.returncode == 0, f"gfleur CBS run failed\n{log[-3000:]}"

    def run(name, sc, collinear=False):
        work = base / name
        work.mkdir()
        _layer_files(prep, work)
        for f in ("embpot11.hdf", "embpotN2.hdf"):
            shutil.copy(cbs / f, work / f)
        shutil.copy(prep / "gf_en.charge", work / "gf_en")
        if collinear:
            for f in work.glob("layer*_inp.xml"):
                _replace(f, {'l_noco="T"': 'l_noco="F"'})
        (work / "gf_inp").write_text(_layers_block(prep) + "&MODE tmat=t, gmat=t, charge=t /\n"
                                     f"&BASIS napw={NAPW} /\n&EMB addemb=t /\n{sc}\n")
        (work / "gf_io").write_text("2\n&IO mode='E', region=1, side=1, rw='r', id=1 /\n"
                                    "&IO mode='E', region=1, side=2, rw='r', id=2 /\n"
                                    "&FILE name='embpot11' /\n&FILE name='embpotN2' /\n")
        proc, log = _run_gf(env.gfleur, work)
        # the collinear control reads noco-prepared files: it crashes when it writes the new density,
        # after the quantities compared by the tests are in the output
        assert collinear or proc.returncode == 0, f"gfleur charge run failed\n{log[-3000:]}"
        return (work / "out").read_text()
    return run


def _mt_total(out):
    """valence charge in the MT spheres, summed over the spin passes of a collinear run"""
    return sum(float(x) for x in re.findall(r"total valence charge in MT-spheres:\s*(\S+)", out))


@pytest.fixture(scope="module")
def noco_z(env):
    return _noco_stack(env, "z")


@pytest.fixture(scope="module")
def noco_x(env):
    return _noco_stack(env, "x")


def test_noco_charge_equals_collinear(noco_z):
    """Moment along z: the noco charge mode gives the collinear first iteration."""
    sc = "&SC iter=1, imix=0, alpha=0.0 /"
    out_n = noco_z("noco", sc)
    out_c = noco_z("collinear", sc, collinear=True)
    assert _neutrality(out_n)[0] == pytest.approx(_neutrality(out_c)[0], abs=1e-6)
    assert _mt_total(out_n) == pytest.approx(_mt_total(out_c), abs=1e-6)
    assert _distances(out_n)[0] == pytest.approx(_distances(out_c)[0], rel=1e-6)
    # the stack holds most of the 26 electrons of the Fe atom (semicore included)
    assert abs(_neutrality(out_n)[0]) < 3.0


@pytest.mark.xfail(reason="the MT density of the local spin frames is not yet rotation invariant "
                          "(MT d charge 5.64 for x against 5.27 for z)", strict=True)
def test_noco_charge_rotation_invariance(noco_z, noco_x):
    """Without SOC the charge does not depend on the direction of the moment."""
    sc = "&SC iter=1, imix=0, alpha=0.0 /"
    out_z = noco_z("noco_inv", sc)
    out_x = noco_x("noco_inv", sc)
    assert _mt_total(out_x) == pytest.approx(_mt_total(out_z), abs=0.02)
    assert _neutrality(out_x)[0] == pytest.approx(_neutrality(out_z)[0], abs=0.02)


@pytest.mark.skipif(not os.environ.get("GFLEUR_MPI"), reason="set GFLEUR_BINARY=<gfleur_MPI> and GFLEUR_MPI='mpirun -np 2'")
def test_charge_mpi_equals_serial(env):
    """The charge mode on several ranks (k-points distributed, density summed per layer, mixing on
    rank 0 from the files) equals the serial run: stack charge, MT charge and the distances of
    three mixing iterations."""
    sc = "&SC iter=3, imix=0, alpha=0.02 /"
    out = {}
    for name, mpi in (("ser", False), ("par", True)):
        run, proc, log = _gfleur_charge(env, f"mpi_{name}", sc, mpi=mpi)
        assert proc.returncode == 0, f"gfleur charge run ({name}) failed\n{log[-3000:]}"
        out[name] = (run / "out").read_text()
    assert _neutrality(out["par"]) == pytest.approx(_neutrality(out["ser"]), abs=1e-6)
    assert _distances(out["par"]) == pytest.approx(_distances(out["ser"]), rel=1e-6)
    assert _mt_total(out["par"]) == pytest.approx(_mt_total(out["ser"]), abs=1e-6)
