"""The quantum conductance that <transport> asks Wannier90 for.

A perfect crystal has no scatterer, so its Landauer transmission is not a smooth curve: it
counts the bands that propagate at that energy. T(E) is an integer staircase in units of the
conductance quantum, and it is an integer however fine the energy mesh is -- which is what
makes it a check and not a stored number.
"""
import pytest

#> The <transport> energies are given in the inp.xml convention -- Hartree, absolute -- and
#> the adapter converts them, because Wannier90 speaks eV throughout. Asserting the ends of
#> the energy axis is what pins that conversion.
HTR_TO_EV = 27.211386245988
E_MIN_HTR = 0.025
E_MAX_HTR = 0.085

#> One file per spin channel, because the case is jspins="2" and Wannier90 names its
#> output after the seedname, which the adapter suffixes per channel. Before that suffix
#> existed the second channel overwrote the first, so a single bare file was all that was
#> left and this test read channel 2 believing it was the whole system.
QC_FILES = ["fleur_wlib_internal.1_qc.dat", "fleur_wlib_internal.2_qc.dat"]
BASE = "./inputfiles/wannier/WannAlNTransport"


def read_qc(path):
    """Read the two-column (energy, transmission) file Wannier90 writes."""
    points = []
    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) == 2:
                try:
                    points.append((float(fields[0]), float(fields[1])))
                except ValueError:
                    pass
    return points


@pytest.mark.wannierlib
@pytest.mark.bulk
def test_transport_transmission_is_quantised(execute_fleur):
    """Run AlN with <transport> and check the transmission counts whole channels."""
    res = execute_fleur(BASE, rm_files=["."], mpi_procs=2)
    for qc_file in QC_FILES:
        assert qc_file in res, f"transport was requested but {qc_file} was not written"
        check_channel(res[qc_file], qc_file)


def check_channel(path, name):
    """Every assertion below holds for each spin channel on its own."""
    points = read_qc(path)
    assert len(points) > 100, f"only {len(points)} energies in {name}"

    assert points[0][0] == pytest.approx(E_MIN_HTR * HTR_TO_EV, abs=5e-3)
    assert points[-1][0] == pytest.approx(E_MAX_HTR * HTR_TO_EV, abs=5e-3)

    trans = [point[1] for point in points]
    assert min(trans) > -1e-8, f"negative transmission in {name}: {min(trans)}"

    #> Wannier90 derives its translation centre from the atomic positions, which the library
    #> API has no way to receive -- atoms_frac is a block, and w90_set_option takes no array
    #> of strings. The adapter names the centre instead. Without that, the centre is 0/0, the
    #> Wannier centres are all NaN and the transmission is identically zero: a run that
    #> reports success and conducts nothing. This is the assertion that sees it.
    conducting = [value for value in trans if value > 0.5]
    assert len(conducting) > 50, f"nothing conducts anywhere in {name}; are the Wannier centres NaN?"

    on_step = [value for value in trans if abs(value - round(value)) < 1e-3]
    assert len(on_step) >= 0.9 * len(trans), \
        f"only {len(on_step)} of {len(trans)} energies sit on an integer step in {name}"
