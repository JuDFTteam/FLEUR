
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
        "plot":pytest.mark.plot,
        "wannier":pytest.mark.wannier,
        "eels":pytest.mark.eels,
        "dos":pytest.mark.dos,
        "edsolver":pytest.mark.edsolver,
        "spinspiral":pytest.mark.spinspiral,
        "hdf":pytest.mark.hdf
    }
    import os
    import re
    test_list=[]
    with open(os.path.dirname(os.path.abspath(__file__))+"/tests.md") as file:
        for s in file:
            #all active tests should be in lines starting with |+|
            if s.startswith("|+|"):
                m=s.split("|")
                #Split the table: 
                #m[0]: empty, m[1]=+ due to start of line
                #m[2]: description, m[3]:directory name, m[4]: marks for pytest
                dir=m[3]
                desc=m[2]
                #use name of directory for testsetname & test-id 
                testsetname,id=dir.split("/")
                #extract marks
                marks=[]
                for mark in m[4].split(","):
                    mark=mark.strip()
                    if len(mark)>1:
                        marks.append(mark_list[mark])
                #only use this test if it matches the testset-name given
                if testsetname==testset:
                    test_list.append(pytest.param(dir,desc,marks=marks,id=id))    
    return test_list                        

