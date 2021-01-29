import pytest
import os
from pytest_dependency import depends

@pytest.mark.fleur_parser
@pytest.mark.masci_tools
def test_fleur_mt_inpxml_parser(request, fleur_test_name, test_file, parser_testdir):
    """
    For each folder (parametrization happens in conftest.py) in parser test,
    try if the fleur parser in masci-tools can handle the parsing without
    crashing, successful and an empty parser log
    """
    pytest.importorskip('masci_tools',minversion='0.4.0-dev3')
    from masci_tools.io.parsers.fleur import inpxml_parser

    depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')

    inpxmlfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'inp.xml'))
    assert os.path.isfile(inpxmlfilepath)

    parser_info = {'parser_warnings': []}
    inp_dict = inpxml_parser(inpxmlfilepath, parser_info_out=parser_info)

    assert inp_dict is not None
    assert isinstance(inp_dict, dict)
    assert inp_dict != {}
    assert parser_info['parser_warnings'] == []


@pytest.mark.fleur_parser
@pytest.mark.masci_tools
def test_fleur_mt_outxml_parser(request, fleur_test_name, test_file, parser_testdir):
    """
    For each folder (parametrization happens in conftest.py) in parser test,
    try if the fleur parser in masci-tools can handle the parsing without
    crashing, successful and an empty parser log
    """

    #These warnings are expected to appear at the moment
    KNOWN_WARNINGS = {'No text found for tag targetStructureClass'} #Will be removed in the future

    pytest.importorskip('masci_tools',minversion='0.4.0-dev3')
    from masci_tools.io.parsers.fleur import outxml_parser
    depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')

    outxmlfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'out.xml'))

    assert os.path.isfile(outxmlfilepath)

    parser_info = {'parser_warnings': []}
    out_dict = outxml_parser(outxmlfilepath, parser_info_out=parser_info)

    assert out_dict is not None
    assert isinstance(out_dict, dict)
    assert out_dict != {}
    assert set(parser_info['parser_warnings']).difference(KNOWN_WARNINGS) == set() #At the moment there is always at least one warning