
import shutil
import subprocess
from pathlib import Path

import pytest

"""
Invariants of the occupation-matrix rotation in global/rotMMPmat.f90.

The spin rotation is reached only from flipcdn (lflip with LDA+U) and the
Hubbard-I setup, and none of the test inputs combine lflip with an LDA+U
species. A small standalone Fortran program therefore checks the routine
against an independent reference. It is an EXCLUDE_FROM_ALL CMake target;
this test builds it on demand and runs it.
"""


@pytest.mark.fleur
def test_rotmmpmat_spin_rotation(build_dir):
    """rotMMPmat(...,spin_rotation=.TRUE.) must equal (su x 1) P (su x 1)^dagger.

    rotMMPmat_dwgn updated the three-block matrix in place and read the
    spin off-diagonal element d(1,2) = CONJG(n(m',m,3)) after (m',m) had already
    been rotated. All blocks were wrong for m > m' (up to 0.9 for a flip by pi),
    the spin-diagonal blocks lost Hermiticity and R^-1 R /= 1. The four-block
    branch wrote d(1,2) to (m',m,4) instead of (m,m',4).
    """
    if shutil.which("make") is None:
        pytest.skip("make not available")

    build = subprocess.run(["make", "rotmmpmat_spin_test"], cwd=build_dir,
                           capture_output=True, text=True)
    assert build.returncode == 0, f"could not build rotmmpmat_spin_test:\n{build.stdout}\n{build.stderr}"

    exe = Path(build_dir) / "tools" / "rotmmpmat_spin_test"
    assert exe.is_file(), f"{exe} was not produced"

    run = subprocess.run([str(exe)], capture_output=True, text=True)
    assert run.returncode == 0, f"rotmmpmat_spin_test failed:\n{run.stdout}\n{run.stderr}"
    assert "failures = 0" in run.stdout, run.stdout
