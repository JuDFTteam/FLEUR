# -*- coding: utf-8 -*-
"""
In this file go all implemented of fixtures for fleur tests which are useful by more then one test

Some code from taken and adapted from redis libtest.py
We stay close to the old test workflow, which is serial, everything is executed in the work dir
Maybe one could also execute tests in parallel, since fleur tests take quite some time
"""
import os
import re
import pytest
import logging
import shutil
pytest_plugins = []

######### Helpers ############
# C: By using os.path instead of pathlib, this will prob fail on Windows

# TODO allow User to specify some of these over the cmd line when executing pytest
# because build dir can have different names
# see https://stackoverflow.com/questions/36141024/how-to-pass-environment-variables-to-pytest#39162893
# or https://adamj.eu/tech/2020/10/13/how-to-mock-environment-variables-with-pytest/

# TODO test log, and better reporting
# TODO time out tests,
# TODO MPI
# TODO set environment i.e open mp
# TODO make tests work on CI
# TODO smoke tests, i.e early end test session if no executable, and so on...
# TODO instead of making all functions fixtures, most of them could be put into a libtest.py and imported
# Through I am not sure how to do this without installing anything. Prob works as long as it is a sub path
# TODO generate docs for devs on webside of test system. Should be done from README.txt. Ideal if possible generate also docs for all tests
# and fixtures from docstring


def test_dir():
    """Get path to the parent test directory other paths are relative to this"""
    test_dir_path = os.path.dirname(os.path.abspath(__file__))
    return test_dir_path

def build_dir():
    """Return directory path where to look for executables"""
    path = "../../build/"
    build_path = os.path.abspath(os.path.join(test_dir(), path))
    return build_path

def work_dir():
    """Return directory path where execute tests in."""
    path = "./work/"
    work_dir_path = os.path.abspath(os.path.join(test_dir(), path))
    return work_dir_path

def failed_dir():
    """Return directory path where execute tests in."""
    path = "./failed_test_results/"
    failed_dir_path = os.path.abspath(os.path.join(test_dir(), path))
    return failed_dir_path


##### Add some markers to pytest to group tests

def pytest_configure(config):
    """
    Here you can add things by a pytest config, could be also part of a seperate file
    So far we add some markers here to be able to execute a certain group of tests
    """
    config.addinivalue_line("markers", "bulk: test with bulk")
    config.addinivalue_line("markers", "film: test with film")
    config.addinivalue_line("markers", "inpgen: test running inpgen")
    config.addinivalue_line("markers", "serial: test running fleur serial")
    config.addinivalue_line("markers", "mpi: test running fleur in parallel")
    config.addinivalue_line("markers", "fast: tests which take < 1 sec to execute")
    config.addinivalue_line("markers", "slow: tests which take < 1 min to execute")
    config.addinivalue_line("markers", "very_slow: tests which take > 1 min to execute")
    config.addinivalue_line("markers", "masci_tools: tests which use function from masci-tools repo")
    config.addinivalue_line("markers", "soc: tests with soc")
    config.addinivalue_line("markers", "lda: tests with ldau")
    config.addinivalue_line("markers", "lo: tests with lo")
    config.addinivalue_line("markers", "forces: tests with forces")
    config.addinivalue_line("markers", "xml: test with xml")
    config.addinivalue_line("markers", "noxml: test with no xml")
    config.addinivalue_line("markers", "collinear: test with collinear")
    config.addinivalue_line("markers", "non-collinear: test with non-collinear")
    config.addinivalue_line("markers", "spinspiral: test with spinspiral")

####################################
########### fixtures
# useful ones from pytest:
# capsys, capsysbinary, capfd
# plugins: data-dir, data-regression, file-regression


##### fixtures for whole test session ####

@pytest.fixture(scope='session')
def fleur_test_session():
    """
    Setup a fleur test session
    - initiate a global logger
    - look for executables to use within that session for all tests
    """
    # put session preparation code here

    yield # now all the tests run

    # put session tear down code here

@pytest.fixture
def set_environment():
    """
    Fixture to set given environment variables
    """
    pass

##### fixtures active EACH test case ####

@pytest.hookimpl(tryfirst=True, hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """This is needed for the base_test_case fixture to work
    """
    # execute all other hooks to obtain the report object
    outcome = yield
    rep = outcome.get_result()

    # set a report attribute for each phase of a call, which can
    # be "setup", "call", "teardown"

    setattr(item, "rep_" + rep.when, rep)


@pytest.fixture(scope='function', autouse=True)
def base_test_case(request, clean_workdir):
    """
    Base fixture for every test case to execute cleanup code after a test
    Write testlog.
    """
    
    workdir = work_dir()
    faildir = failed_dir()

    yield # test is running
    
    # clean up code goes here:
   
    method_name = request.node.name
    if request.node.rep_call.failed:
        #log('Test {} failed :('.format(method_name))
        # if failed move test result to failed dir, will replace dir if existent
        destination = os.path.abspath(os.path.join(faildir, method_name))
        if os.path.isdir(destination):
            shutil.rmtree(destination)
        shutil.move(workdir, destination)
        os.mkdir(workdir)
    else:
        clean_workdir()


##### other fixtures ####

@pytest.fixture
def execute_inpgen(set_environment, inpgen_binary):
    """
    Fixture which returns an execute_inpgen function
    """
    def _execute_inpgen(cmdline_param=None, test_file_folder=None, exclude=[]):
        """
        Function which copies the input files
        executes inpgen with the given cmdline_param
        and returns a path to the folder with the results
        
        :param cmdline_param:
        :param test_file_folder: string relative path to a folder containing all files, to copy for the test
        :param exclude: list of strings file names to exclude form copy
        
        :return: a dictionary of the form 'filename :filepath'
        """
        import subprocess
        
        workdir = work_dir()
        testdir = test_dir()
        if cmdline_param is None:
            cmdline_param = []

        if test_file_folder is not None:
            abspath = os.path.abspath(os.path.join(testdir, test_file_folder))
            source = os.listdir(abspath)
            for files in source:
                if files not in exclude:
                    shutil.copy(os.path.abspath(os.path.join(abspath, files)), workdir)
        
        #print(inpgen_binary)
        arg_list = [inpgen_binary] + cmdline_param
        #print(arg_list)
        os.chdir(workdir)
        with open(f"{workdir}/stdout", "w") as f_stdout:
            with open(f"{workdir}/stderr", "w") as f_stderr:
                subprocess.run(arg_list + ["-no_send"], stdout=f_stdout, stderr=f_stderr, check=True)

        result_files = {}
        source = os.listdir(workdir) # Notice this is simple and not recursive,
        # TODO if we have output directories use os.walk or so instead
        for files in source:
            result_files[files] = os.path.abspath(os.path.join(workdir, files))      
        os.chdir(testdir)

        return result_files
        
    return _execute_inpgen


@pytest.fixture
def execute_fleur(set_environment, fleur_binary):
    """
    Fixture which returns an execute_fleur function
    """
    def _execute_fleur(cmdline_param=None, test_file_folder=None, exclude=[]):
        """
        Function which copies the input files
        executes fleur with the given cmdline_param
        and returns a path to the folder with the results
        
        :param cmdline_param:
        :return:
        """
        
        import subprocess
        
        workdir = work_dir()
        testdir = test_dir()
        if cmdline_param is None:
            cmdline_param = []

        if test_file_folder is not None:
            abspath = os.path.abspath(os.path.join(testdir, test_file_folder))
            source = os.listdir(abspath)
            for files in source:
                if files not in exclude:
                    shutil.copy(os.path.abspath(os.path.join(abspath, files)), workdir)
        
        fleur, parallel = fleur_binary
        mpiruncmd = []
        arg_list = mpiruncmd + [fleur] + cmdline_param
        #print(arg_list)
        os.chdir(workdir)
        with open(f"{workdir}/stdout", "w") as f_stdout:
            with open(f"{workdir}/stderr", "w") as f_stderr:
                subprocess.run(arg_list, stdout=f_stdout, stderr=f_stderr, check=True)

        result_files = {}
        source = os.listdir(workdir) # Notice this is simple and not recursive,
        # TODO if we have output directories use os.walk or so instead
        for files in source:
            result_files[files] = os.path.abspath(os.path.join(workdir, files))      
        os.chdir(testdir)

        return result_files
        

    return _execute_fleur

@pytest.fixture
def grep_exists():
    """returns the grep_exits function
    """
    def _grep_exists(filepath, expression):
        """for an expression in a file
        Args:
            filepath ([string, path]): path to the file to search for
            expression (string): 'string' to look for
    
        :return: Bool, if exists
        """
        exists = False
        with open(filepath, "r") as file1:
            for line in file1:
                if re.search(expression, line):
                    exists = True
                    break
        return exists
    return _grep_exists


@pytest.fixture
def grep_number():
    """returns the grep number function
    """
    def _grep_number(filepath, expression, split, index=1, first_only=True):
        
        numbers = []
        with open(filepath, "r") as file1:
            for line in file1:
                if re.search(expression, line):
                    res = line.split(expression)[index]
                    try:
                        number = float(res)
                    except ValueError as exc: # There is still something after the number
                        number = float(re.findall("\d+\.\d+", res)[0])
                    numbers.append(number)
        if first_only:
            return numbers[0]
        else:
            return numbers
        
    return _grep_number
    
def check_value_outfile(filepath, before_str, after_str, expected, delta):
    """
    Check a value in the out file, coped form libtest.py
    """
    exp_idx = 0
    found = False
    passed = False
    with open(filepath, "r") as f:
        for line in f.readlines():
            if(before_str in line):
                value_string = line.split(before_str)[-1]
                value_string = value_string.split(after_str)
                # remove empty strings
                value_string = [i for i in value_string if i is not ""][0]
                value = float(value_string)
             
                if(expected[exp_idx] is not None):
                    if(abs(value - expected[exp_idx]) < delta):
                        #log_info(f"PASSED [{exp_idx}]: {before_str} found: {value} expected: {expected[exp_idx]}")
                        passed = True
                    else:
                        #log_info(f"FAILED [{exp_idx}]: {before_str} found: {value} expected: {expected[exp_idx]}")
                        errors += 1
                exp_idx += 1
                found = True
    if(len(expected) != exp_idx):
        #log_error("Number of expected values disagree with found values.")
        #log_error(f"before_str: '{before_str}''")
        passed = False
    if(errors != 0):
        passed = False
    return found, passed

def parse_inp_xml():
    pass

def parse_out_xml():
    pass

def parse_stderr():
    pass

def parse_stdout():
    pass

def parse_file():
    pass

@pytest.fixture
def clean_workdir():
    """
    Fixture which returns a function to delete files in the test work dir
    the default is to delete all
    """
    def _clean_workdir(filelist=['all']):
        """
        Delete files in the test work dir the default is to delete all
        :param filelist (list, optional): filesnames to delete. Defaults to ['all'].
        """
        workdir = work_dir()
        if 'all' in filelist:
            shutil.rmtree(workdir) # this remove the work dir too.
            os.mkdir(workdir)
        else:
            for files in filelist:
                os.remove(os.path.abspath(os.path.join(workdir, files)))
        
    return _clean_workdir





@pytest.fixture(scope='session')
def inpgen_binary():
    """
    Fixture returning the path to a inpgen executable
    """
    fleur_dir = build_dir()
    if(fleur_dir[-1] == "/"):
        fleur_dir = fleur_dir[:-1]
    
    if(os.path.isfile(f"{fleur_dir}/inpgen")):
        binary = f"{fleur_dir}/inpgen"
        #logging.info(f"Use {self.binary} as executable")
    else:
        binary = None
    #else:
    #   logging.warning("Can not find any executables")

    return binary

@pytest.fixture(scope='session')
def fleur_binary():
    """
    Fixture returning the path to a fleur executable
    """
    fleur_dir = build_dir()
    parallel = False
    if(fleur_dir[-1] == "/"):
        fleur_dir = fleur_dir[:-1]

    if(os.path.isfile(f"{fleur_dir}/fleur_MPI")):
        binary = f"{fleur_dir}/fleur_MPI"
        #logging.info(f"Use {self.binary} as executable")
        parallel = True

    elif(os.path.isfile(f"{fleur_dir}/fleur")):
        binary = f"{fleur_dir}/fleur"
        #logging.info(f"Use {self.binary} as executable")
    else:
        binary = None
        

    #else:
    #   logging.warning("Can not find any executables")

    return binary, parallel


@pytest.fixture
def inpxml_etree():
    """Returns the etree generator"""
    def _get_etree(filepath):
        from lxml import etree
        with open(filepath, 'r') as inpxmlfile:
            tree = etree.parse(inpxmlfile)
        return tree
    return _get_etree



