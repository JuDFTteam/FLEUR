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

Binaries are located via GFPREP_BINARY or <repo>/build*/gfprep; fleur and gfleur are taken from
the same directory. Missing binaries skip the tests.
"""
import re
import shutil
import subprocess
from collections import namedtuple
from pathlib import Path

import pytest

from test_gfleur_transport import (BULK_INPUTS, MT_SETUPS, _find_gfprep, _layer_files, _replace, _run,
                                   _sibling)

KGRID = "kgrid=4,4"
NAPW = 8
SCF_ITER = 14
SCF_MIX = "imix=7, maxiter=14, alpha=0.02"

Env = namedtuple("Env", "gfprep gfleur bulk base prep lead cbs")


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


def _gfleur_charge(env, name, sc, copy_cdn=True):
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
    proc, log = _run(env.gfleur, run)
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
    proc, log = _run(gfleur, cbs)
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
