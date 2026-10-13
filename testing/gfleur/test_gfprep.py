"""Tests for gfprep, the conversion of a bulk supercell into the layered
gfleur setup.

A tetragonal 4-atom Cu supercell along (001) is converged briefly with fleur,
gfprep cuts it into four one-plane layers (reading the bulk from a sibling
directory), and gfleur is started on the result.

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
BULK_INPUTS = HERE / "inputfiles" / "cu_supercell"
GF_EN = HERE / "inputfiles" / "2layer_film" / "gf_en"
NLAYERS = 4


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


def _run(cmd, cwd, timeout=900):
    env = dict(os.environ)
    proc = subprocess.run([str(cmd)], cwd=cwd, capture_output=True, text=True,
                          timeout=timeout, env=env)
    return proc, proc.stdout + "\n" + proc.stderr


@pytest.fixture(scope="module")
def gfprep_run(tmp_path_factory):
    """Converge the bulk supercell and run gfprep on it (once per module)."""
    gfprep = _find_gfprep()
    if gfprep is None:
        pytest.skip("no gfprep binary found (set GFPREP_BINARY or build the gfprep target)")
    fleur = _sibling(gfprep, ["fleur", "fleur_MPI"])
    if fleur is None:
        pytest.skip(f"no fleur binary next to {gfprep}")

    base = tmp_path_factory.mktemp("gfprep")
    bulk = base / "bulk"
    shutil.copytree(BULK_INPUTS, bulk)
    proc, log = _run(fleur, bulk)
    assert proc.returncode == 0, f"bulk fleur run failed\n{log[-4000:]}"
    assert (bulk / "cdn.hdf").is_file()

    work = base / "gf"
    work.mkdir()
    (work / "gfprep_inp").write_text("&GFPREP bulkdir='../bulk/' /\n")
    proc, log = _run(gfprep, work)
    assert proc.returncode == 0, f"gfprep failed\n{log[-4000:]}"
    return gfprep, work


def _checks(out_text):
    """Self-check lines of gfprep.out: (quantity, layer, region) -> (max|f|, deviation)"""
    pat = re.compile(r"gfprep check (\w+) layer (\d+) (interstitial|muffin-tins): "
                     r"max\|f\|, max deviation:\s+(\S+)\s+(\S+)")
    return {(m[0], int(m[1]), m[2]): (float(m[3]), float(m[4])) for m in pat.findall(out_text)}


def test_gfprep_outputs(gfprep_run):
    """gfprep writes the layer inputs, the gf_inp template and the HDF files."""
    _, work = gfprep_run
    for i in range(1, NLAYERS + 1):
        xml = (work / f"layer{i}_inp.xml").read_text()
        assert "<filmLattice" in xml
        assert xml.count("<filmPos") == 1
    gfinp = (work / "gf_inp").read_text()
    assert f"nlayers={NLAYERS}" in gfinp
    for f in ("gf_setup.hdf", "gf_pot.hdf", "gf_cdn.hdf", "gfprep_planes.txt"):
        assert (work / f).is_file(), f"{f} missing"


def test_gfprep_hdf_layout(gfprep_run):
    """gf_setup.hdf and gf_pot.hdf contain the groups gfleur reads."""
    _, work = gfprep_run
    h5py = pytest.importorskip("h5py")
    with h5py.File(work / "gf_setup.hdf", "r") as f:
        for i in range(1, NLAYERS + 1):
            assert f"layer-{i}/stepinfo" in f
            assert f"layer-{i}/embdesc/left" in f
            assert abs(f[f"layer-{i}/stepinfo"].attrs["c"][0] - 3.40965) < 1e-5
    with h5py.File(work / "gf_pot.hdf", "r") as f:
        for i in range(1, NLAYERS + 1):
            for d in ("pot/vr", "pot/pw", "pot/starinfo", "coulpot/coul"):
                assert f"layer-{i}/{d}" in f


def test_gfprep_conversion_accuracy(gfprep_run):
    """The converted layer functions reproduce the bulk ones."""
    _, work = gfprep_run
    checks = _checks((work / "gfprep.out").read_text())
    assert len(checks) == 4 * NLAYERS
    for (quantity, layer, region), (fmax, dev) in checks.items():
        if region == "muffin-tins":
            assert dev < 1e-10 * fmax, f"{quantity} layer {layer} MT deviation {dev}"
        elif quantity == "potential":
            assert dev < 1e-2 * fmax, f"potential layer {layer} interstitial deviation {dev}"
        else:
            assert dev < 1e-3 * fmax, f"density layer {layer} interstitial deviation {dev}"


def _run_gfleur(gfprep_run, tmp_path, mode):
    gfprep, work = gfprep_run
    gfleur = _sibling(gfprep, ["gfleur", "gfleur_MPI"])
    if gfleur is None:
        pytest.skip(f"no gfleur binary next to {gfprep}")
    for f in work.iterdir():
        if f.name.startswith("layer") or f.name.startswith("gf_") or f.suffix == ".xsd":
            shutil.copy(f, tmp_path / f.name)
    gfinp = (tmp_path / "gf_inp").read_text().replace("tmat=t, gmat=t", mode)
    (tmp_path / "gf_inp").write_text(gfinp)
    shutil.copy(GF_EN, tmp_path / "gf_en")
    proc, log = _run(gfleur, tmp_path)
    out_file = tmp_path / "out"
    out = out_file.read_text() if out_file.is_file() else ""
    return proc, out, log


def test_gfleur_reads_gfprep_setup(gfprep_run, tmp_path):
    """gfleur initialises all layers of the gfprep output and completes."""
    proc, out, log = _run_gfleur(gfprep_run, tmp_path, "tmat=t, gmat=f")
    assert proc.returncode == 0, f"gfleur exited with {proc.returncode}\n{(out + log)[-4000:]}"
    assert "GF-SETUP done" in out
    assert f"layers : {NLAYERS}" in out


@pytest.mark.xfail(reason="gfleur port: T-matrices are never written (gf_write_tmat has no caller)",
                   strict=False)
def test_gfleur_gmat_with_gfprep_potential(gfprep_run, tmp_path):
    """gfleur in gmat mode, which reads the converted potentials and gf_setup.hdf."""
    proc, out, log = _run_gfleur(gfprep_run, tmp_path, "tmat=t, gmat=t")
    assert proc.returncode == 0, f"gfleur exited with {proc.returncode}\n{(out + log)[-4000:]}"
