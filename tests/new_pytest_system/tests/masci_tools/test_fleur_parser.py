import pytest
from helpers.utils import RUN_PARSER_TESTS

@pytest.mark.skip('test not implemented')
#@pytest.mark.skipif(not RUN_PARSER_TESTS)
#@pytest.mark.parameterize()
@pytest.mark.fleur_parser
def test_fleur_mt_parser():
    """
    For each folder in parser test, try if the fleur parser in masci-tools can handle the parsing without
    crashing, successful and an empty parser log
    """
    assert False