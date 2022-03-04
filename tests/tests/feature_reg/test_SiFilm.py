
import pytest

@pytest.mark.serial
@pytest.mark.film
@pytest.mark.plot
@pytest.mark.xml
def test_SiFilmPlotXML(execute_fleur, grep_number, grep_exists):
    """Fleur Si film slice plot

    Simple test of Fleur with two steps:
    1.Generate a starting density and perform a single iteration.
    2.Generate a slice plot.
    """
    test_file_folder = './inputfiles/SiFilmPlotXML/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml', 'sym.out']) #kpts.xml 
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert ('cdn.hdf' in res_file_names) or ('cdn1' in res_file_names)

    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp-2.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['out', 'denIn.xsf']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    # unit cell
    assert grep_exists(res_files['denIn.xsf'], ".0000000 11.36143")
    assert grep_exists(res_files['denIn.xsf'], "1.93419")
    # atom positions
    assert grep_exists(res_files['denIn.xsf'], "3.35012")
    assert grep_exists(res_files['denIn.xsf'], ".39481")
    assert grep_exists(res_files['denIn.xsf'], "2.23341")
    # density values
    assert grep_exists(res_files['denIn.xsf'], "1.9231") # line 24
    assert grep_exists(res_files['denIn.xsf'], "7.9875") # line 521
    assert grep_exists(res_files['denIn.xsf'], "10.2889") # line 523

@pytest.mark.serial
@pytest.mark.film
@pytest.mark.plot
@pytest.mark.xml
def test_SiFilmSlicePlotXML(execute_fleur, grep_number, grep_exists):
    """Fleur Si film slice plot

    Simple test of Fleur with three steps:
    1.Generate a starting density and perform a single iteration.
    2.Generate a charge density slice.
    3.Generate a plot from the slice.
    """
    test_file_folder = './inputfiles/SiFilmSlicePlotXML/'

    # Stage 1
    res_files = execute_fleur(test_file_folder, only_copy=['inp.xml'])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert ('cdn.hdf' in res_file_names) or ('cdn1' in res_file_names)

    # Stage 2
    res_files = execute_fleur(test_file_folder, only_copy=[['inp-2.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    # now test output cdn_slice
    assert ('cdn_slice' in res_file_names) or ('cdn_slice.hdf' in res_file_names)

    # Stage 3
    res_files = execute_fleur(test_file_folder, only_copy=[['inp-3.xml', 'inp.xml']])
    res_file_names = list(res_files.keys())
    should_files = ['out', 'slice.xsf']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    # unit cell
    assert grep_exists(res_files['slice.xsf'], ".0000000 11.36143")
    assert grep_exists(res_files['slice.xsf'], "1.93419")
    # atom positions
    assert grep_exists(res_files['slice.xsf'], "3.35012")
    assert grep_exists(res_files['slice.xsf'], ".39481")
    assert grep_exists(res_files['slice.xsf'], "2.23341")
    # density values
    assert grep_exists(res_files['slice.xsf'], "1.43968")  # line 24
    assert grep_exists(res_files['slice.xsf'], "2.54654")  # line 289
    assert grep_exists(res_files['slice.xsf'], "2.86363")  # line 521
    assert grep_exists(res_files['slice.xsf'], "0.20168")  # line 523
    assert grep_exists(res_files['slice.xsf'], "3.10848")  # line 3018
