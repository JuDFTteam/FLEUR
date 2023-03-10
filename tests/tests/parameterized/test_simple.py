
import pytest

#Here noco tests can be simply added by:
# creating a new entry in the list
# first string is the directory of the input
# second string is the short descriptions
# in addition marks can be set.
# it is also a good idea to add an explicit id to avoid complex directory names
basic_tests = []
noco_tests=[]
extra_tests=[]
film_tests=[]
hybrid_tests=[]
force_tests=[]

testsets={
    "basic":basic_tests,
    "noco":noco_tests,
    "extra":extra_tests,
    "film":film_tests,
    "hybrid":hybrid_tests,
    "forces":force_tests
}
mark_list = {
    "band":pytest.mark.band,
    "fast":pytest.mark.fast,
    "bulk":pytest.mark.bulk
}
import os
import re
with open(os.path.dirname(os.path.abspath(__file__))+"/tests.md") as file:
    for s in file:
        m=s.split("|")
        if len(m)==3:
            dir=m[1]
            desc=m[0]
            print("DIR:",dir)
            testsetname,id=dir.split("/")
            marks=[]
            for mark in m[2].split(","):
                try:
                    marks.append(mark_list[mark])
                except:
                    print("Unkown mark",mark)
            testsets[testsetname].append(pytest.param(dir,desc,marks=marks,id=id))            

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.parametrize(("dir","desc"), basic_tests)
def test_basic(dir,desc,default_fleur_test):
    """
    """
    
    assert default_fleur_test(dir)

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.parametrize(("dir","desc"), noco_tests)
def test_noco(dir,desc,default_fleur_test):
    """
    """
    
    assert default_fleur_test(dir)

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.extra
@pytest.mark.parametrize(("dir","desc"), extra_tests)
def test_extra(dir,desc,default_fleur_test):
    """
    """
    
    assert default_fleur_test(dir)

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.film
@pytest.mark.parametrize(("dir","desc"), film_tests)
def test_film(dir,desc,default_fleur_test):
    """
    """
    
    assert default_fleur_test(dir)

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.hybrid
@pytest.mark.parametrize(("dir","desc"), hybrid_tests)
def test_hybrid(dir,desc,default_fleur_test):
    """
    """
    
    assert default_fleur_test(dir)

@pytest.mark.serial
@pytest.mark.fleur
@pytest.mark.force
@pytest.mark.parametrize(("dir","desc"), force_tests)
def test_force(dir,desc,default_fleur_test):
    """
    """
    
    assert default_fleur_test(dir)
