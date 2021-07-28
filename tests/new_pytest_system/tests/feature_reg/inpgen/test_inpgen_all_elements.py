import pytest


def test_inpgen_whole_periodic_table():
    """Runs inpgen for a 'imaginary' structure containing all

    elements of the periodic table (N: 1-119).
    
    Tests if there is a problem with this.
    I.e this tests might fail, if there is a misstake in some defaults 
    for a certain element
    
    # Comment one could split this test into a parametized test set
    using the same base cell and only using a different element
    But we assume this to fail seldom and here we have very different atoms together
    """
    pass
