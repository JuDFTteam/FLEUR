import pytest
import os
from pytest_plugins.pytest_dependency import depends

@pytest.mark.fleur_parser
@pytest.mark.masci_tools
@pytest.mark.dos
@pytest.mark.hdf
def test_banddos_parser_default_dos(request, fleur_test_name, test_file, parser_testdir):
    pytest.importorskip('masci_tools',minversion='0.4.0')
    from masci_tools.io.parsers.hdf5 import HDF5Reader
    from masci_tools.io.parsers.hdf5.recipes import FleurDOS
    depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')

    banddosfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'banddos.hdf'))

    assert os.path.isfile(banddosfilepath)

    with HDF5Reader(banddosfilepath) as h5reader:
        data, attributes = h5reader.read(recipe=FleurDOS)

    assert data is not None
    assert isinstance(data, dict)
    assert data != {}

    assert attributes is not None
    assert isinstance(attributes, dict)
    assert attributes != {}

@pytest.mark.fleur_parser
@pytest.mark.masci_tools
@pytest.mark.jdos
@pytest.mark.hdf
def test_banddos_parser_jdos(request, fleur_test_name, test_file, parser_testdir):
    pytest.importorskip('masci_tools',minversion='0.4.0')
    from masci_tools.io.parsers.hdf5 import HDF5Reader
    from masci_tools.io.parsers.hdf5.recipes import FleurJDOS
    depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')

    banddosfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'banddos.hdf'))

    assert os.path.isfile(banddosfilepath)

    with HDF5Reader(banddosfilepath) as h5reader:
        data, attributes = h5reader.read(recipe=FleurJDOS)

    assert data is not None
    assert isinstance(data, dict)
    assert data != {}

    assert attributes is not None
    assert isinstance(attributes, dict)
    assert attributes != {}

@pytest.mark.fleur_parser
@pytest.mark.masci_tools
@pytest.mark.orbcomp
@pytest.mark.hdf
def test_banddos_parser_orbcomp(request, fleur_test_name, test_file, parser_testdir):
    pytest.importorskip('masci_tools',minversion='0.4.0')
    from masci_tools.io.parsers.hdf5 import HDF5Reader
    from masci_tools.io.parsers.hdf5.recipes import FleurORBCOMP
    depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')

    banddosfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'banddos.hdf'))

    assert os.path.isfile(banddosfilepath)

    with HDF5Reader(banddosfilepath) as h5reader:
        data, attributes = h5reader.read(recipe=FleurORBCOMP)

    assert data is not None
    assert isinstance(data, dict)
    assert data != {}

    assert attributes is not None
    assert isinstance(attributes, dict)
    assert attributes != {}




@pytest.mark.fleur_parser
@pytest.mark.masci_tools
@pytest.mark.mcd
@pytest.mark.hdf
def test_banddos_parser_mcd(request, fleur_test_name, test_file, parser_testdir):
    pytest.importorskip('masci_tools',minversion='0.4.0')
    from masci_tools.io.parsers.hdf5 import HDF5Reader
    from masci_tools.io.parsers.hdf5.recipes import FleurMCD
    depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')

    banddosfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'banddos.hdf'))

    assert os.path.isfile(banddosfilepath)

    with HDF5Reader(banddosfilepath) as h5reader:
        data, attributes = h5reader.read(recipe=FleurMCD)

    assert data is not None
    assert isinstance(data, dict)
    assert data != {}

    assert attributes is not None
    assert isinstance(attributes, dict)
    assert attributes != {}

@pytest.mark.fleur_parser
@pytest.mark.masci_tools
@pytest.mark.band
@pytest.mark.hdf
def test_banddos_parser_default_band(request, fleur_test_name, test_file, parser_testdir):
    pytest.importorskip('masci_tools',minversion='0.4.0')
    from masci_tools.io.parsers.hdf5 import HDF5Reader
    from masci_tools.io.parsers.hdf5.recipes import FleurBands
    depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')

    banddosfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'banddos.hdf'))

    assert os.path.isfile(banddosfilepath)

    with HDF5Reader(banddosfilepath) as h5reader:
        data, attributes = h5reader.read(recipe=FleurBands)

    assert data is not None
    assert isinstance(data, dict)
    assert data != {}

    assert attributes is not None
    assert isinstance(attributes, dict)
    assert attributes != {}
