
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
        "wannierlib":pytest.mark.wannierlib,
        "eels":pytest.mark.eels,
        "dos":pytest.mark.dos,
        "edsolver":pytest.mark.edsolver,
        "spinspiral":pytest.mark.spinspiral,
        "hdf":pytest.mark.hdf,
        "dfpt":pytest.mark.dfpt,
        "mpi":pytest.mark.mpi,
        "elpa":pytest.mark.elpa,
        "scalapack":pytest.mark.scalapack,
        "chase":pytest.mark.chase
        }
    import os
    import re
    test_list=[]
    with open(os.path.dirname(os.path.abspath(__file__))+"/tests.md") as file:
        for s in file:
            m=s.split("|")
            testMarker = " "
            if (len(m) > 1):
                testMarker = m[1].strip()
            #all active tests should be in lines starting with |+|
            if testMarker == "+":
                #Split the table: 
                #m[0]: empty, m[1]=+ due to start of line
                #m[2]: description, m[3]:directory name, m[4]: marks for pytest, m[5]: Remarks, m[6]: cmdline, m[7]: mpi_procs
                dir=m[3].strip()
                desc=m[2].strip()
                if (len(m)>6):
                    cmdline=m[6].split(",")
                    for i in range(len(cmdline)):
                       cmdline[i] = cmdline[i].strip()
                    if (len(cmdline)==1 and cmdline[0].strip()==""):
                        cmdline=None
                else:
                    cmdline=None
                if(len(m)>7):
                    if (m[7].strip().isdigit()):
                        mpi_procs=int(m[7].strip())
                    else:
                        mpi_procs=None
                else:
                    mpi_procs=None
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
                    test_list.append(pytest.param(dir,desc,cmdline,mpi_procs,marks=marks,id=id))    
    return test_list                        

