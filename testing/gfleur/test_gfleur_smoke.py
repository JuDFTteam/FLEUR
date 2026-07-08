"""Smoke tests for the Green's-function embedding code (gfleur/gfleur_MPI).

These tests are intentionally light-weight and self-contained: they do not
hook into the fleur pytest fixtures (which are built around the fleur/inpgen
binaries). They only check that the gfleur executable runs the full
initialisation + driver chain of a minimal two-layer film setup and shuts
down cleanly - i.e. that per-layer fleur_init, the k x layer x energy setup,
the step-function/2D-matrix stores and gfleur_execute all survive end to end.

The physics kernels (T-matrices, transmission, SCF charge) need converged
per-layer potentials and reference data that are not part of this fixture, so
they are out of scope here; the two-layer input runs in tmat mode with
iter = 0, so the self-consistency loop body is skipped and the run exercises
only the ported infrastructure.

The gfleur binary is located via (in order):
  * the GFLEUR_BINARY environment variable, or
  * <repo>/build*/gfleur  (serial) discovered next to the source tree.
If no binary and no schema files are found the tests skip, so the file is
safe to collect in a plain (non-gfleur) build.
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent  # <repo>/testing/gfleur -> <repo>
INPUTS = HERE / "inputfiles" / "2layer_film"
SCHEMA_DIR = REPO / "src" / "fleur" / "io" / "xml"


def _find_gfleur():
    env = os.environ.get("GFLEUR_BINARY")
    if env and Path(env).is_file():
        return Path(env)
    # search build directories next to the source tree
    candidates = sorted(REPO.glob("build*/gfleur"))
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return c
    return None


def _stage(work: Path):
    work.mkdir(parents=True, exist_ok=True)
    for f in INPUTS.iterdir():
        shutil.copy(f, work / f.name)
    for xsd in ("FleurInputSchema.xsd", "FleurOutputSchema.xsd"):
        src = SCHEMA_DIR / xsd
        if not src.is_file():
            pytest.skip(f"schema {xsd} not found at {src}")
        shutil.copy(src, work / xsd)


@pytest.fixture(scope="module")
def gfleur_binary():
    binary = _find_gfleur()
    if binary is None:
        pytest.skip("no gfleur binary found (set GFLEUR_BINARY or build the gfleur target)")
    return binary


def test_gfleur_2layer_driver_completes(gfleur_binary, tmp_path):
    """gfleur runs gf_init + gfleur_execute for a 2-layer film and finishes."""
    _stage(tmp_path)
    proc = subprocess.run([str(gfleur_binary)], cwd=tmp_path,
                          capture_output=True, text=True, timeout=600)
    out_file = tmp_path / "out"
    out = out_file.read_text() if out_file.is_file() else ""
    combined = out + "\n" + proc.stdout + "\n" + proc.stderr

    assert proc.returncode == 0, f"gfleur exited with {proc.returncode}\n{combined[-4000:]}"
    # the per-layer initialisation finished ...
    assert "GF-SETUP done" in out, "gf_init did not complete\n" + combined[-4000:]
    # ... and the driver printed its setup summary for both layers
    assert "GFLEUR setup summary" in out, "gfleur_execute was not reached\n" + combined[-4000:]
    assert "layers : 2" in out, "unexpected layer count in the driver summary"
    # the step-function cache was produced during setup
    assert (tmp_path / "gf_steps.hdf").is_file(), "gf_steps.hdf was not generated"
