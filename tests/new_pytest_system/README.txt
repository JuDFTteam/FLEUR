###################
# How to run tests:
# execute in the 'new_pytest_system' folder:
# To run all tests:
$pytest

#or instead of pytest, since it ships with python
$python3 -m pytest

# execute all tests in a given subtest folder or file
$pytest tests/bar
$pytest tests/bar/test_foo.py

# stop at first failure
$pytest -

#checkout
pytest --help

Pytest provides two ways to run the subset of the test suite.

    Select tests to run based on substring matching of test names.
    Select tests groups to run based on the markers applied.

pytest -k <substring> -v
pytest -m <markername>

# example markers are:
bulk,film,xml,collinear,soc,lo,ldau,non-collinear,spinspiral,inversion,forces,...
 

# –duration=N
##################
# Folder layout:

- tests: here and under the subfolders go all further test folders, tests are organized through the functionality they test
- work: here we run the tests
- failed_run: folder where all failed tests are preserved of the last test sesstion
- inputfiles: here all the input files for the tests go into a folder

- conftest.py
This is a central pytest file, here all fixtures i.e helpers go (what former has partly been taken care of by the 'libtest.py' file or the 'scripts/*')


####################
# pytest basics

https://docs.pytest.org/en/stable/

pytest automaticly discovers all tests, as long as they are in a python file starting with 'test' and have a function or class defined inside starting with 'test'

fixtures are helper functions, which can be put in the .py file of the tests to be available there or in the conftest.py file to be available globally.
We use them for prepartion or teardown code, also basicly to make anything availble we need throughout the different tests
fixtures with scope 'session' are executed one a session (for example the binaries, because we use the same fleur exe for all tests)
fixtures with scope 'function' are auto available or auto executed for every test

usefull pytest decorators
@pytest.mark.skip() # skips the test (shown in report)
@pytest.mark.skipif() # skips the test under a certain condition
@pytest.mark.parameterize() # creates a test case for each 'parameter'
@pytest.mark.<markername> # mark test if a given label to allow subset execution

# pytest plugins, there are a bunch of them, we use two

pytest-datadir
pytest-regression


#################
# how to create a new test

1. Create a new folder under inputfiles where you put all the input files needed for your test
2. Optional add a new test_<some_name>.py file somewhere under the tests directory where your test functions goes, or add you test function to an already existing file


example test function:

@pytest.mark.bulk                       # markers to group tests
@pytest.mark.collinear                  # markers to group tests
def test_foo_bar():                     # what fixures you want to include
    """Describe what you test
    """
    


# TODO: How to do coverage of fortan, enable what we had before
# TODO?: if wanted we can move also the inputfiles of the tests in the tests sub folders, I made it like this, because pytest will create there data subfolder and to be very clear what is autogerated and what has the developer put in




