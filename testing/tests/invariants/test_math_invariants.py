
import shutil
import subprocess
from pathlib import Path

import pytest

"""
Invariants of the shared math routines.

These drive a small standalone Fortran program rather than a FLEUR calculation,
because the affected routines are reachable from very few test inputs. The
program is an EXCLUDE_FROM_ALL CMake target, so it costs nothing in a normal
build; this test builds it on demand and runs it.
"""


@pytest.mark.fleur
def test_clebsch_selection_rule(build_dir):
    """clebsch must vanish unless m1+m2 = M (#804).

    The old guard truncated the real difference am+bm-cm, so a mismatch of -1
    landed in the same bucket as 0 and fell through to the factorial evaluation.
    Two of the cases in the driver returned 1/sqrt(3) and sqrt(2) instead of
    zero -- the latter is impossible for a Clebsch-Gordan coefficient.

    Reachable in production from the Wannier SOC projections in
    wannierlib_rad_twd, which sweep m over a shell at fixed jm.
    """
    if shutil.which("make") is None:
        pytest.skip("make not available")

    build = subprocess.run(["make", "clebsch_test"], cwd=build_dir,
                           capture_output=True, text=True)
    assert build.returncode == 0, f"could not build clebsch_test:\n{build.stdout}\n{build.stderr}"

    exe = Path(build_dir) / "tools" / "clebsch_test"
    assert exe.is_file(), f"{exe} was not produced"

    run = subprocess.run([str(exe)], capture_output=True, text=True)
    assert run.returncode == 0, f"clebsch_test failed:\n{run.stdout}\n{run.stderr}"
    assert "CLEBSCH SELECTION RULE TEST: PASS" in run.stdout, run.stdout
