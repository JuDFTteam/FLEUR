

def test_inpgen_scaling_performance():
    """Run inpgen for a very large structure 80 000 atoms, with namelists and time it. 
    
    Even on a very slow machine this should not take longer then a minute.

    Notice, the time constrains is not so important if you add new feature and it is a bit longer, 
    adjust the time.

    This test exists, because once inpgen scaled very bad (exponential) and 
    running it on even 1000 atoms would have taken weeks if not longer.
    Also Fleur is a MaX flagship and getting ready for Exascale means also being able to simulate very
    large systems, possible millions. 
    """
    pass