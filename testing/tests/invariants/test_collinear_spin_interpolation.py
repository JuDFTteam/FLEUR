"""The spin operator interpolated on a collinear, two-channel run.

Two spin channels are two eigenproblems wannierised one after the other, and most
operators cannot cross them: the cross-spin block needs both gauges at once, so it
exists only as the combined 2N real-space file. sigma_z is the exception, because it
is diagonal in the channel index -- each channel carries its own -- and until this was
written the code refused it anyway, aborting on <operator name="spin"/> inside
<interpolation> whenever jspins=2 and the run was not spinorial.

WHAT CARRIES THE INFORMATION HERE. Within one channel sigma_z is that channel's own
overlap, so the coarse block is the unit matrix and the interpolated expectation value
is +-1 at every k', exactly, whatever the mesh and whatever the gauge. That is the
assertion: not a stored number but a value fixed by orthonormality, which is why the
case may be re-converged, re-wannierised or re-meshed without touching the tolerance.

It is not a tautology on the FLEUR side. The block is an INTEGRAL -- muffin-tin spheres
plus interstitial plus, on a film, the vacuum -- so a region dropped from the sum, a
normalisation lost, or a provider handing back a transposed block all move the diagonal
off one. The run also carries a Fortran-side gate at every k that aborts before any of
this is written; what reaches here has already passed it, and the printed line below is
that gate reporting.

WHAT IT DOES NOT CATCH, stated so nobody reads more into a green than is there:

  * the R-set weights. O(R) is delta_{R,0} times the identity, so sum_R e^{ikR} O(R)
    is the identity whatever the degeneracies are; a wrong ndegen is invisible to this
    test and needs the Hamiltonian, which has structure in R.
  * the transverse components as a physical result. They are filled with zeros by
    construction, because a spin-pure eigenstate has none. Asserting them is still
    worth doing, but for a different reason: together with the +-1 in z they pin the
    COMPONENT ORDER of the dataset, so a permuted axis turns red instead of silently
    relabelling x as z.
  * the split of the two allocations behind it. Whether the cross-spin block is also
    reserved for a run that never exports it is not observable from outside; what this
    pins is that the interpolation-only request runs at all.

The path domain is deliberately kept. Interpolating back onto the coarse mesh is an
exact round trip of the discrete transform and would hold even if the transform were
wrong in both directions; the path asks for k' the gauge was never built at.
"""
import pytest

BASE = "./inputfiles/wannier/WannFeBccInterpSpin"

#> The Fortran gate uses 1e-3 and the line it prints has six decimals, so the tightest
#> thing that can be read back from it is 1e-6. Anything that actually breaks this is
#> O(1), not O(1e-5).
TOL_PRINTED = 1.0e-4
#> On the datasets there is no printing in the way: these are the stored reals.
TOL_DATA = 1.0e-8


def _channel_deviations(out_text):
    """{channel: (max||diag|-1|, largest off-diagonal)} as the Fortran gate printed them."""
    found = {}
    for line in out_text.splitlines():
        if "max | |diag|-1 |" not in line:
            continue
        head, rest = line.split(":", 1)
        parts = rest.split("=")
        found[int(head.split()[1])] = (float(parts[1].split()[0]), float(parts[2].split()[0]))
    return found


@pytest.mark.fleur
@pytest.mark.wannierlib
@pytest.mark.hdf
def test_collinear_spin_interpolation_runs_and_blocks_are_unit(execute_fleur):
    """The run completes, and its own per-k gate says both channel blocks are the identity.

    Completing IS the differential half: against the code before this change the same
    input aborts in wannierlib_request with "the spin operator cannot be interpolated",
    and execute_fleur fails a run that does not report success.
    """
    res = execute_fleur(BASE, rm_files=['.'])
    out = open(res["out"], encoding="utf-8", errors="replace").read()

    assert "cannot be interpolated" not in out, "the refusal is back"

    dev = _channel_deviations(out)
    assert set(dev) == {1, 2}, f"expected both channels to report, got {sorted(dev)}"
    for ch, (diag, offdiag) in dev.items():
        assert diag < TOL_PRINTED, f"channel {ch}: |diag| is {diag} away from 1"
        assert offdiag < TOL_PRINTED, f"channel {ch}: off-diagonal reaches {offdiag}"

    for name in ("WF1_interp.hdf", "WF2_interp.hdf"):
        assert name in res, f"{name} was not written"


@pytest.mark.fleur
@pytest.mark.wannierlib
@pytest.mark.hdf
def test_collinear_interpolated_spin_is_plus_minus_one(execute_fleur):
    """<sigma_z> = +1 over channel 1 and -1 over channel 2, at every interpolated k.

    Skipped where h5py is absent, which includes the CI venv; the test above runs there
    and already covers the refusal and the gate.
    """
    h5py = pytest.importorskip("h5py")
    numpy = pytest.importorskip("numpy")

    res = execute_fleur(BASE, rm_files=['.'])

    for name, sign in (("WF1_interp.hdf", +1.0), ("WF2_interp.hdf", -1.0)):
        checked = 0
        with h5py.File(res[name], "r") as fh:
            for group in fh:
                if "spin" not in fh[group]:
                    continue
                #> Fortran writes (3, numWann, nk) and h5py sees it reversed, so the
                #> Cartesian component is the LAST axis here, not the first.
                s = numpy.asarray(fh[group]["spin"])
                assert s.ndim == 3 and s.shape[-1] == 3, f"{name}/{group}: shape {s.shape}"
                worst_z = numpy.max(numpy.abs(s[..., 2] - sign))
                worst_t = numpy.max(numpy.abs(s[..., 0:2]))
                assert worst_z < TOL_DATA, \
                    f"{name}/{group}: <sigma_z> is {worst_z} away from {sign:+.0f}"
                assert worst_t < TOL_DATA, \
                    f"{name}/{group}: a transverse component reaches {worst_t}"
                checked += s[..., 2].size
        assert checked > 0, f"{name} carries no spin dataset"
