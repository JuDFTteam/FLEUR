import pytest
import os
import warnings
from pytest_plugins.pytest_dependency import depends

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

    parser_info = {}
    inp_dict = inpxml_parser(inpxmlfilepath, parser_info_out=parser_info)

    if any("Schema available for version" in warning for warning in parser_info['parser_warnings']):
        for warning in parser_info['parser_warnings'].copy():
            if 'Input file does not validate against the schema' in warning:
                parser_info['parser_warnings'].remove(warning)
            if 'Schema available for version' in warning:
                parser_info['parser_warnings'].remove(warning)

    assert inp_dict is not None
    assert isinstance(inp_dict, dict)
    assert inp_dict != {}
    assert parser_info['parser_errors'] == []
    assert parser_info['parser_critical'] == []
    assert parser_info['parser_warnings'] == []


@pytest.mark.fleur_parser
@pytest.mark.masci_tools
def test_fleur_mt_outxml_parser(request, fleur_test_name, test_file, expected_failure, parser_testdir):
    """
    For each folder (parametrization happens in conftest.py) in parser test,
    try if the fleur parser in masci-tools can handle the parsing without
    crashing, successful and an empty parser log
    """
    #These warnings are expected to appear at the moment
    KNOWN_WARNINGS = {'No values found for attribute l_f'}

    pytest.importorskip('masci_tools',minversion='0.4.0-dev3')
    from masci_tools.io.parsers.fleur import outxml_parser
    depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')

    outxmlfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'out.xml'))

    assert os.path.isfile(outxmlfilepath)

    parser_info = {}
    out_dict = outxml_parser(outxmlfilepath, parser_info_out=parser_info)

    if any("Schema available for version" in warning for warning in parser_info['parser_warnings']):
        for warning in parser_info['parser_warnings'].copy():
            if 'Output file does not validate against the schema' in warning:
                parser_info['parser_warnings'].remove(warning)
            if 'Schema available for version' in warning:
                parser_info['parser_warnings'].remove(warning)

    assert out_dict is not None
    assert isinstance(out_dict, dict)
    assert out_dict != {}
    if expected_failure:
        assert parser_info['parser_errors'] != []
    else:
        assert parser_info['parser_errors'] == []
    assert parser_info['parser_critical'] == [] #Critical errors should not happen in any case
    warning_set = set(parser_info['parser_warnings']).difference(KNOWN_WARNINGS)
    if warning_set != set():
        print(warning_set)#Maybe better to write to a file
        warnings.warn(f'The outxml_parser encountered {len(warning_set)} unexpected warning(s) for test {fleur_test_name}')
