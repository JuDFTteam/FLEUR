"""
Regression test for <wannierization numIter="0">, which the schema documents as
"skips the spread minimisation and keeps the projection gauge".

It did not. The adapter guarded the option with `IF (this%num_iter > 0)`, using zero
as the sentinel for "the input said nothing", so a requested zero never reached
Wannier90 and Wannier90 applied its own default of a hundred iterations. An input
asking for no minimisation got a hundred steps of it and said nothing about it: the
only trace is the iteration table in the run's own stdout, which nothing reads. The
schema and the code disagreed and the code won, silently.

What makes it testable is that zero iterations leaves a signature nothing else
produces: the Wannier functions at the end are the ones at the start. Wannier90
prints that pair itself, as the Initial State and Final State blocks,

  Initial State
   Sum of centres and spreads (  0.089688, -0.269063, -0.179375 )     8.59701030
  Final State
   Sum of centres and spreads (  0.089688, -0.269063, -0.179375 )     8.59701030

so the check is an identity between two blocks of the same run and needs no stored
number. It compares the centres too, not just the spread: a spread can come back to
its starting value while the functions have moved, three numbers in step cannot.

The test is written against the real text and not a guess at it. Both blocks print
the same "Sum of centres and spreads" line, which is why they are paired by their
own label rather than by counting lines; and the label is separated from that line
by one row per Wannier function, which is why the match has to span them without
being allowed to run past the next label.

Needs neither numpy nor h5py, so unlike its neighbour in this directory it runs in
CI like any other test.
"""

import re

import pytest

BASE = "./inputfiles/wannier/WannFeBccNoIter"
#> Silicon rather than iron for the disentangled half: eight k-points and eight bands
#> against four Wannier functions, which is the cheapest case in the house with a frozen
#> window genuinely narrower than the outer one. A case whose two windows coincide runs
#> zero disentanglement iterations whatever is asked of it, and would pass this test
#> without testing anything.
BASE_DIS = "./inputfiles/wannier/WannSiDiamondNoIter"

#> The line is printed to six and eight decimals, so it can be compared at the
#> precision of the printing and no finer. Anything that moves the gauge is O(1e-2).
TOL = 1e-7

#> " Initial State \n  Sum of centres and spreads (  x,  y,  z )   8.59701030 ".
#> Keeping the label in the match is what lets the blocks be paired; both blocks
#> print the same line, so counting them separately gets 2n against n.
STATE = re.compile(
    r"(Initial|Final) State\b(?:(?!State\b).)*?Sum of centres and spreads\s*"
    r"\(([^)]*)\)\s+([-\d.]+)", re.S)
#> An iteration row: "      0     0.860E+01   0.0000000000    8.5970102960   0.00  <-- CONV".
#> The border and header lines carry the same <-- CONV tag, hence the leading integer.
ITER_ROW = re.compile(r"^\s*(\d+)\s+[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+.*<-- CONV", re.M)
#> The same, for the disentanglement: its table is tagged <-- DIS and is counted
#> separately because a disentangled run prints a Sum of centres and spreads line on
#> every one of its iterations, not only in the two State blocks. Counting lines rather
#> than tags is how that output gets misread.
DIS_ROW = re.compile(r"^\s*(\d+)\s+[-\d.E+]+\s+[-\d.E+]+\s+[-\d.E+]+.*<-- DIS", re.M)


@pytest.mark.wannierlib
@pytest.mark.bulk
def test_num_iter_zero_keeps_the_projection_gauge(execute_fleur):
    """numIter="0" must reach Wannier90 as zero, not fall through to its default."""
    res = execute_fleur(BASE, rm_files=["."], mpi_procs=2)
    assert "out" in res, "the run wrote no out file; there is nothing to read"
    out = open(res["out"]).read()

    #> jspins=2 without spinors, so the whole wannierisation runs once per channel and
    #> every block below appears twice. Both have to obey this, not just the first.
    estados = STATE.findall(out)
    assert len(estados) >= 2 and len(estados) % 2 == 0, (
        f"found {len(estados)} State blocks; the run did not get as far as "
        f"wannierising both spin channels")

    #> The whole claim, stated directly: nothing past iteration zero ever ran.
    iters = [int(n) for n in ITER_ROW.findall(out)]
    assert iters, "no iteration rows at all; this test cannot tell what happened"
    assert max(iters) == 0, (
        f'numIter="0" was asked for and Wannier90 ran up to iteration {max(iters)}. '
        f"That is the old defect: the adapter guarded the option with > 0, so a "
        f"requested zero never reached Wannier90 and its own default was applied")

    #> And the consequence, which is what the schema actually promises: the functions
    #> at the end are the functions at the start, centres and spread alike.
    for ch, (ini, fin) in enumerate(zip(estados[::2], estados[1::2]), start=1):
        assert (ini[0], fin[0]) == ("Initial", "Final"), (
            f"spin channel {ch}: the blocks came out as {ini[0]}/{fin[0]}, so this "
            f"test is not reading the pair it thinks it is")

        c_ini = [float(v) for v in ini[1].split(",")]
        c_fin = [float(v) for v in fin[1].split(",")]
        assert len(c_ini) == len(c_fin) == 3, f"spin channel {ch}: centres unreadable"
        for axis, (a, b) in enumerate(zip(c_ini, c_fin)):
            assert abs(b - a) < TOL * max(1.0, abs(a)), (
                f"spin channel {ch}: the sum of centres moved along axis {axis} from "
                f"{a} to {b} although no minimisation was asked for")

        s_ini, s_fin = float(ini[2]), float(fin[2])
        assert abs(s_fin - s_ini) < TOL * max(1.0, abs(s_ini)), (
            f"spin channel {ch}: the spread moved from {s_ini:.8f} to {s_fin:.8f} "
            f"although no minimisation was asked for, so the projection gauge was "
            f"not kept")


@pytest.mark.wannierlib
@pytest.mark.bulk
def test_dis_num_iter_zero_keeps_the_projection_subspace(execute_fleur):
    """<disentanglement numIter="0"> must be expressible, and must reach Wannier90.

    The schema used to type this attribute as positiveInteger, so zero could not even be
    written down. That was stricter than the library it drives: Wannier90 rejects
    dis_num_iter < 0, not zero, and its loop runs from 1, so zero simply does not
    execute while dis_project and dis_proj_froz -- the steps that actually build the
    subspace out of the projections and the frozen window -- run before the loop and
    always. Zero therefore names a real object: the projection subspace, unoptimised.

    This case asks for zero on both counts, so neither Omega_I nor the spread is
    minimised and the run is the pure projection answer from end to end.
    """
    res = execute_fleur(BASE_DIS, rm_files=["."], mpi_procs=2)
    assert "out" in res, "the run wrote no out file; there is nothing to read"
    out = open(res["out"]).read()

    #> If the two windows coincided there would be nothing to disentangle and the count
    #> below would be zero for a reason that has nothing to do with numIter.
    assert "Frozen states" in out or "<-- DIS" in out, (
        "this run shows no sign of a disentanglement stage at all, so a count of zero "
        "iterations says nothing")

    dis = [int(n) for n in DIS_ROW.findall(out)]
    assert not dis, (
        f'<disentanglement numIter="0"> was asked for and Wannier90 ran '
        f"{len(dis)} iterations, up to {max(dis)}. The zero did not reach it")

    #> The wannierise count is zero here too, so the same signature as the other test
    #> has to hold: what comes out is what went in.
    iters = [int(n) for n in ITER_ROW.findall(out)]
    assert iters and max(iters) == 0, (
        f"the wannierise stage ran up to iteration {max(iters) if iters else None} "
        f'although numIter="0" was asked for')

    estados = STATE.findall(out)
    assert len(estados) >= 2 and len(estados) % 2 == 0, (
        f"found {len(estados)} State blocks; the run did not get as far as wannierising")
    for ch, (ini, fin) in enumerate(zip(estados[::2], estados[1::2]), start=1):
        assert (ini[0], fin[0]) == ("Initial", "Final"), (
            f"block {ch}: came out as {ini[0]}/{fin[0]}, not the pair this expects")
        s_ini, s_fin = float(ini[2]), float(fin[2])
        assert abs(s_fin - s_ini) < TOL * max(1.0, abs(s_ini)), (
            f"block {ch}: the spread moved from {s_ini:.8f} to {s_fin:.8f} although "
            f"nothing was asked to be minimised")
