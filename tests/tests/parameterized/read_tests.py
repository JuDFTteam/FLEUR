
import pytest

def read_tests(testset):
    """
    Read test.md and return list of tests in the given testset
    """
    mark_list = {
        "band":pytest.mark.band,
        "fast":pytest.mark.fast,
        "bulk":pytest.mark.bulk,
        "film":pytest.mark.film,
        "soc":pytest.mark.soc,
        "orbcomp":pytest.mark.orbcomp,
        "mcd":pytest.mark.mcd,
        "ldau":pytest.mark.ldau,
        "libxc":pytest.mark.libxc,
        "forcetheorem":pytest.mark.forcetheorem,
        "dos":pytest.mark.dos
    }
    import os
    import re
    test_list=[]
    with open(os.path.dirname(os.path.abspath(__file__))+"/tests.md") as file:
        for s in file:
            m=s.split("|")
            if len(m)==3:
                dir=m[1]
                desc=m[0]
                testsetname,id=dir.split("/")
                marks=[]
                for mark in m[2].split(","):
                    mark=mark.strip()
                if mark:
                    marks.append(mark_list[mark])
                    if testsetname==testset:
                        test_list.append(pytest.param(dir,desc,marks=marks,id=id))    
    return test_list                        

