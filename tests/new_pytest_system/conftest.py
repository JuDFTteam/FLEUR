# -*- coding: utf-8 -*-
"""
In this file contains all implemented fixtures for fleur tests which are useful by more then one test

Some code was taken and adapted from redis libtest.py
We stay close to the old test workflow, which is serial, everything is executed in the work dir.
Also so far we did not rename or restructure the tests, in a new directory tree
Maybe one could also execute tests in parallel, since fleur tests take quite some time
"""
import os
import re
import sys
import time
import pytest
import shlex
import logging
import shutil
sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
sys.path.append(os.path.join(os.path.dirname(__file__), 'pytest_plugins'))
# Now we can import everything that is in helpers, but be careful about name clashing
from helpers.utils import RUN_PARSER_TESTS, MASCI_TOOLS_VERSION_STR

pytest_plugins = ("pytest_plugins.pytest_dependency",
                  "pytest_plugins.pytest_modify_terminal_report")

LOGGER = logging.getLogger(__name__)

# pytest logfile: --log-file=path
# --log-file-level
#TODO other optional aiida tests
# test log, and better reporting https://docs.pytest.org/en/stable/logging.html
# python test_out.py | tee myoutput.log, # maybe with capsys, capsysbinary, capfd
# use | tee pytest_stdout to pipe output
# TODO time out tests,
# TODO smoke tests, i.e early end test session if no executable, and so on...
# TODO generate docs for devs on webside of test system. Should be done from README.txt. Ideal if possible generate also docs for all tests
# and fixtures from docstring
# TODO: install pytest-xdist to run pytest -n=2 to execute tests in parallel
# pytest --durations=0 for showing slowest tests
# we could also create a timeing fixture around each test for example see
# https://stackoverflow.com/questions/51490166/how-to-time-tests-with-pytest
# TODO: All hdf tests do not test for values in hdf files, this is bad
# Write grep so that if hdf file, reads it and greps in there
# Speed up slowest test by staging the last cdn and running onlyone iteration

# C: maybe instead of running everything in the work dir, which is a problem, for parallel test execution
# use the native tmp dir instead and implement something to allow to copy temp dir results of a test to a next test.

# C: One can also run non python tests, so we could if we want to run the old tests as they are...
# https://docs.pytest.org/en/stable/example/nonpython.html



######### Helpers ############
# C: By using os.path instead of pathlib, this will prob fail on Windows

# The current fleur test workflow is as follows, create Testing dir in build dir
# there is the workdir in which the test runs, after the test is finished, they copy the work dir
# they create logs there.


def test_dir():
    """Get path to the parent test directory defined by the position of this conftest.py file
    other paths are relative to this"""
    test_dir_path = os.path.dirname(os.path.abspath(__file__))
    return test_dir_path

def get_build_dir(pytestconfig):
    """Return directory path where to look for executables, some other paths are relative to this"""
    path = pytestconfig.getoption("build_dir")
    if os.path.isabs(path):
        build_path = path
    else:
        build_path = os.path.abspath(os.path.join(test_dir(), path))
    return build_path

def get_work_dir(build_dir):
    """Return directory path where execute tests in"""
    path = "./Testing/work/"
    work_dir_path = os.path.abspath(os.path.join(build_dir, path))
    return work_dir_path

def get_stage_dir(build_dir):
    """Return directory path where tests results can be stage in"""
    path = "./Testing/stage_dir/"
    work_dir_path = os.path.abspath(os.path.join(build_dir, path))
    return work_dir_path

def get_failed_dir(build_dir):
    """Return directory path where execute tests in"""
    path = "./Testing/failed_test_results/"
    failed_dir_path = os.path.abspath(os.path.join(build_dir, path))
    return failed_dir_path

def get_parser_testdir(build_dir):
    """Return directory path where execute tests in"""
    path = "./Testing/parser_testdir/"
    parser_testdir_path = os.path.abspath(os.path.join(build_dir, path))
    return parser_testdir_path

@pytest.fixture(scope='session')
def build_dir(pytestconfig):
    """Return directory path where to look for executables, some other paths are relative to this"""
    return get_build_dir(pytestconfig)

@pytest.fixture(scope='session')
def cleanup(pytestconfig):
    """Flag if only failed test results are stored or all tests"""
    return not pytestconfig.getoption("--no-cleanup")


@pytest.fixture(scope='session')
def work_dir(build_dir):
    """Return directory path where execute tests in"""
    return get_work_dir(build_dir)

@pytest.fixture(scope='session')
def stage_dir(build_dir):
    """Return directory path where tests results can be stage in"""
    return get_stage_dir(build_dir)

@pytest.fixture(scope='session')
def failed_dir(build_dir):
    """Return directory path where execute tests in"""
    return get_failed_dir(build_dir)

@pytest.fixture(scope='session')
def parser_testdir(build_dir):
    """Return directory path where execute tests in"""
    return get_parser_testdir(build_dir)

##### change pytest configuration and execution:

def pytest_addoption(parser):
    """We add an option to pytest to parse the build dir
    """
    parser.addoption("--build_dir", action="store", default="../../build/",
                    help='Path to the build dir with fleur exe')
    parser.addoption("--no-cleanup", action="store_true",
                     help='Move every test result to failed dir and not delete any data')
    parser.addoption("--runevery", action="store", default=None, help='Run every x test')
    parser.addoption("--testoffset", action="store", default=None, help='Do not run first x tests')
    parser.addoption("--skipmarkers", action="store",
                     default="", help="skip tests with these markers")
    #parser.addoption("--testing_dir", action="store", default="")


# Does not work with the fixtures
@pytest.mark.trylast
def pytest_report_header(config):#libs):
    """Add further information to header"""
    # ggf add version info

    build_dir = get_build_dir(config)
    work_dir = get_work_dir(build_dir)
    failed_dir = get_failed_dir(build_dir)
    path_fleur, para = get_fleur_binary(build_dir)
    path_inp = get_inpgen_binary(build_dir)
    mpiruncmd = get_mpi_command(dict(os.environ), '{mpi_procs}', para)
    mpiruncmd = ' '.join(mpiruncmd)
    libs = []

    reporter = config.pluginmanager.getplugin("terminalreporter")
    reporter.write_sep('-',title='Fleur test session')
    add_header_strings = [f"Fleur exe: {path_fleur}",
                          f"Inpgen exe: {path_inp}",
                          f"NOT linked libraries: {libs}",
                          f"Running tests in: {work_dir}",
                          f"Failed tests will be copied to: {failed_dir}",
                          f"Parser tests {'will' if RUN_PARSER_TESTS else 'will not'} be run (masci-tools version {MASCI_TOOLS_VERSION_STR})",
                          f"Default MPI command: {mpiruncmd if para else 'None (serial build)'}",
                           "Documentation for the test system can be found here:",
                           "    https://iffgit.fz-juelich.de/fleur/fleur/-/wikis/Testing/Pytest-test-system\n",
                           "Now cleaning work, failed and parser_test directories...\n",]

    return add_header_strings

def read_cmake_config(configfilepath):
    """
    reads build path
    and which markers to ignore by cmake
    """
    marker_list = []

    if os.path.exists(configfilepath):
        with open(configfilepath, 'r') as conf:
            content = conf.readlines()
        for line in content:
            if 'excl_flags=' in line:
                marker_string = str(line.split('excl_flags=')[1])
                marker_string = marker_string.replace('"', '')
                marker_string = marker_string.strip()
                marker_list = marker_string.split()

    marker_list.append('noci') # these tests are always ignored
    return marker_list


_parser_tests_collected = False
fleur_tests = set()
inpgen_tests = set()

def pytest_generate_tests(metafunc):
    """Generate tests during collection

        This is modified to be able to use parametrized tests for the parsers and automatically performing the tests
        for all available fleur/inpgen tests

        Fleur and inpgen tests are detected based on their usage of the `execute_fleur` or `execute_ingpen` fixture
        They are added to (module level) sets as tuples (function name, function module file, all given pytest markers)
        and when a parser test is encountered (detection based on the fixture `fleur_test_name`) the corresponding fixture
        is parametrized with all collected tests

        Selection of tests for parsers based on markers is also possible but not yet used

    """
    #This is needed to be able to issue warnings if fleur tests are missed for the parsers
    global _parser_tests_collected

    if 'execute_fleur' in metafunc.fixturenames:
        markers = tuple({mark.name for mark in metafunc.function.pytestmark})
        fleur_tests.add((metafunc.function.__name__,
                         metafunc.module.__file__.replace(os.path.dirname(os.path.abspath(__file__))+'/','')) + markers)
        if _parser_tests_collected:
            metafunc.config.issue_config_time_warning(UserWarning('Fleur test collected after parser test. Missed: '\
                                                                  f"{metafunc.function.__name__}"),2)
    if 'execute_inpgen' in metafunc.fixturenames:
        markers = tuple({mark.name for mark in metafunc.function.pytestmark})
        inpgen_tests.add((metafunc.function.__name__,
                          metafunc.module.__file__.replace(os.path.dirname(os.path.abspath(__file__))+'/','')) + markers)
        if _parser_tests_collected:
            metafunc.config.issue_config_time_warning(UserWarning('Inpgen test collected after parser test. Missed: '\
                                                                  f"{metafunc.function.__name__}"),2)

    if 'fleur_test_name' in metafunc.fixturenames:
        _parser_tests_collected = True
        if 'inpxml' in metafunc.function.__name__:
            test_info = fleur_tests.union(inpgen_tests)
        else:
            test_info = fleur_tests

        default_markers = set(['fleur_parser', 'masci_tools', 'hdf']) #These markers should be ignored for selection

        markers = {mark.name for mark in metafunc.function.pytestmark}

        required_markers = markers - default_markers
        #Here we could select tests based on the markers of the Test (at the moment we just discard the marker info here)
        #This is useful for the eventual tests of banddos parsers, nmmpmat parser, ...
        test_info = {(info[0], info[1]) for info in test_info if all(marker in info[2:] for marker in required_markers)}

        metafunc.parametrize('fleur_test_name, test_file', test_info, ids=[info[0] for info in test_info])


# To modify the collected tests AFTER collections
def pytest_collection_modifyitems(session, config, items):
    """After test collection modify collection.

    Depending on how fleur was compiled we mark some tests
    with certain markers to be skiped.

    """
    from _pytest.mark import deselect_by_keyword, deselect_by_mark
    # I have not found any other way, but to import from a protective method...
    deselect_by_keyword(items, config)
    deselect_by_mark(items, config)

    markers_cmd_to_skip = config.getoption("--skipmarkers").split(',')
    filename = 'pytest_incl.py'
    path = config.getoption("build_dir")
    test_dir_path = os.path.dirname(os.path.abspath(__file__))

    confile = os.path.abspath(os.path.join(test_dir_path, path))
    confile = os.path.join(confile, filename)
    marker_list = read_cmake_config(confile)
    marker_list = marker_list + markers_cmd_to_skip
    marker_list = sorted(list(set(marker_list)))
    print("\nExcluding tests with the following markers in 'pytest_incl.py': ", marker_list)
    run_every = config.getoption("runevery")
    if run_every is None:
        run_every = 1

    testoffset = config.getoption("testoffset")
    if testoffset is None:
        testoffset = 0
    print(f'Running every {run_every}st test with offset {testoffset}, others will be skiped.')

    if not RUN_PARSER_TESTS:
        marker_list = marker_list + ['masci_tools']

    skip_unselected = pytest.mark.skip(reason='This test was unselected by commandline arguments given.')
    #Add to all tests marked with masci_tools
    #pytest.importorskip("masci_tools")

    # only left with selected items.
    deselection_items = []
    for i, item in enumerate(items):
        # be careful with this because it also applies for test classes,
        # and tests which might depend on each other
        mod = (i+int(testoffset))%int(run_every)
        if mod !=0:
            item.add_marker(skip_unselected)
        for marker in marker_list:
            if marker in item.keywords:
                deselection_items.append(item)

        '''
        for marker in marker_to_skip:
            if marker in item.keywords:
                item.add_marker(skip_markers[marker])
        '''
    items[:] = [item for item in items if item not in deselection_items]
    config.hook.pytest_deselected(items=deselection_items)



##### Add some markers to pytest to group tests
# control skipping test on command line options, for test collection
# https://docs.pytest.org/en/stable/example/simple.html?highlight=pytest_configure

def pytest_configure(config):
    """
    Here you can add things by a pytest config, could be also part of a separate file
    So far we add some markers here to be able to execute a certain group of tests
    We make them all lowercaps as convention
    """
    # run mode markers
    config.addinivalue_line("markers", "inpgen: test running inpgen")
    config.addinivalue_line("markers", "fleur: test running fleur")
    config.addinivalue_line("markers", "serial: test running fleur serial")
    config.addinivalue_line("markers", "mpi: test running fleur in parallel")
    config.addinivalue_line("markers", "mpionly: test which need a MPI capable fleur to run.")
    config.addinivalue_line("markers", "fast: tests which take < 1 sec to execute")
    config.addinivalue_line("markers", "slow: tests which take < 1 min to execute")
    config.addinivalue_line("markers", "very_slow: tests which take > 1 min to execute")
    config.addinivalue_line("markers", "xml: test with xml")
    config.addinivalue_line("markers", "noxml: test with no xml")
    config.addinivalue_line("markers", "gpu: this test will a GPU capbale fleur version")

    config.addinivalue_line("markers", "noci: this test will not be run on CI ")
    # the reason for this is that it is not run in the old set.

    # feature makers
    config.addinivalue_line("markers", "bulk: test with bulk")
    config.addinivalue_line("markers", "film: test with film")
    config.addinivalue_line("markers", "band: testing bandstructure")
    config.addinivalue_line("markers", "dos: testing DOS")
    config.addinivalue_line("markers", "jdos: testing jDOS")
    config.addinivalue_line("markers", "orbcomp: testing orbcomp")
    config.addinivalue_line("markers", "mcd: testing mcd")
    config.addinivalue_line("markers", "hybrid: testing hybrid functionals")
    config.addinivalue_line("markers", "eigpara: test testing eig para")
    config.addinivalue_line("markers", "soc: tests with soc")
    config.addinivalue_line("markers", "ldau: tests with LDA+U")
    config.addinivalue_line("markers", "lo: tests with LOs")
    config.addinivalue_line("markers", "forces: tests with forces")
    config.addinivalue_line("markers", "relaxation: tests involving relaxation")
    config.addinivalue_line("markers", "collinear: test with collinear")
    config.addinivalue_line("markers", "non_collinear: test with non-collinear")
    config.addinivalue_line("markers", "spinspiral: test with spinspiral calculations")
    config.addinivalue_line("markers", "greensfunction: test with greensfunction")
    config.addinivalue_line("markers", "edsolver: test for fleur using the edsolver library")
    config.addinivalue_line("markers", "magnetism: test with magnetism")
    config.addinivalue_line("markers", "plot: tests testing a plot feature")
    config.addinivalue_line("markers", "eels: test with eels")
    config.addinivalue_line("markers", "gw: test for gw interface")
    config.addinivalue_line("markers", "interface: tests testing some interface")

    # main libs
    config.addinivalue_line("markers", "hdf: tests needing hdf")
    config.addinivalue_line("markers", "libxc: test for fleur using libxc")
    config.addinivalue_line("markers", "wannier: test for fleur using wannier") # TODO account for differnet wannier versions?
    config.addinivalue_line("markers", "wannier4: test for fleur using wannier 4D calculations")
    config.addinivalue_line("markers", "wannier5: test for fleur using wannier 5D calculations")
    config.addinivalue_line("markers", "masci_tools: tests which use functions from masci-tools repo")
    config.addinivalue_line("markers", "fleur_parser: tests testing fleur parsers or generate files for them")

    # solvers, ffts and other libs
    config.addinivalue_line("markers", "edsolver: test needing the edsolver")
    config.addinivalue_line("markers", "cusolver: test needing the cusolver")
    config.addinivalue_line("markers", "fftmkl: test needing the fftmkl")
    config.addinivalue_line("markers", "fftw: test needing the fftw")
    config.addinivalue_line("markers", "spfft: test needing the spfft")
    config.addinivalue_line("markers", "magma: test needing the magma")
    config.addinivalue_line("markers", "progthread: test needing the progthread")
    config.addinivalue_line("markers", "elpa: test needing the elpa")
    config.addinivalue_line("markers", "elpaonenode: test needing the elpaonenode")
    config.addinivalue_line("markers", "chase: test needing the chase")
    config.addinivalue_line("markers", "scalapack: test needing the scalapack")


# Replace for now with an even shorter version in pytest_plugins/pytest_modify_terminal_report.py
# def pytest_runtest_logreport(report):
#     """
#     Short path names in report
#     """
#     # remove the tests/ in all paths printed
#     report.nodeid = report.nodeid[5:]

####################################
########### fixtures
# useful ones from pytest:
# capsys, capsysbinary, capfd
# plugins: data-dir, data-regression, file-regression


##### fixtures for whole test session ####

@pytest.fixture(scope='session', autouse=True)
def fleur_test_session(parser_testdir, failed_dir, work_dir, inpgen_binary, fleur_binary, request):
    """
    Setup a fleur test session
    - cleanup directories
    - print some info
    - look for executables to use within that session for all tests
    """
    # This fixture gets executed just before the first test runs

    # put session preparation code here

    parser_testdir_path = parser_testdir
    failed_dir_path = failed_dir
    work_dir_path = work_dir
    libs = []
    capmanager = request.config.pluginmanager.getplugin("capturemanager")

    # First we print out some info, which will show at the end of session:
    path_fleur, para = fleur_binary
    path_inp = inpgen_binary
    add_header_string = "\n#########################################\n"
    add_header_string += "Fleur test session information:\n\n"
    add_header_string += "Fleur exe: {} \n".format(path_fleur)
    add_header_string += "Inpgen exe: {} \n".format(path_inp)
    add_header_string += "NOT linked libraries: {} \n".format(str(libs))
    add_header_string += "Running tests in: {} \n".format(work_dir_path)
    add_header_string += "Failed tests will be copied to: {} \n\n".format(failed_dir_path)
    add_header_string += "Cleaning now work, failed and parser_test directories...\n"
    add_header_string += "#########################################\n"
    #LOGGER.info(add_header_string)
    #with capmanager.global_and_fixture_disabled():
        #print(add_header_string)

    # We clean before and not after test session, because that way
    # people can investigate the files

    # create work dir if not existent
    if not os.path.isdir(work_dir_path):
        os.makedirs(work_dir_path)  # Will create also Testing dir if not existent
    else:
        shutil.rmtree(work_dir_path)
        os.mkdir(work_dir_path)

    # clean parser_testdir
    if os.path.isdir(parser_testdir_path):
        shutil.rmtree(parser_testdir_path)
    os.mkdir(parser_testdir_path)

    # clean failed_dir
    if os.path.isdir(failed_dir_path):
        shutil.rmtree(failed_dir_path)
    os.mkdir(failed_dir_path)

    time.sleep(0.5) # otherwise somehow the first test fails

    yield # now all the tests run

    # put session tear down code here
    #print("Fleur test session ended.")

@pytest.fixture
def set_environment_session():
    """
    Fixture to set given environment variables
    """
    # So far this is done on the outside, and for the whole test set
    # or on an individual test basis
    yield

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
def base_test_case(request, work_dir, failed_dir, clean_workdir, cleanup):
    """
    Base fixture for every test case to execute cleanup code after a test
    Write testlog.
    """

    if 'fleur_parser' not in request.keywords:

        workdir = work_dir
        faildir = failed_dir

        yield # test is running

        # clean up code goes here:

        method_name = request.node.name
        if request.node.rep_call.failed or not cleanup:
            #log('Test {} failed :('.format(method_name))
            # if failed move test result to failed dir, will replace dir if existent
            destination = os.path.abspath(os.path.join(faildir, method_name))
            if os.path.isdir(destination):
                shutil.rmtree(destination)
            shutil.copytree(workdir, destination)
        clean_workdir()
    else:
        yield


##### other fixtures ####

@pytest.fixture
def execute_inpgen(inpgen_binary, work_dir):
    """
    Fixture which returns an execute_inpgen function
    """
    def _execute_inpgen(test_file_folder=None, cmdline_param=None, exclude=[], only_copy=[], rm_files=[]):
        """
        Function which copies the input files
        executes inpgen with the given cmdline_param
        and returns a path to the folder with the results

        :param test_file_folder: string relative path to a folder containing all files, to copy for the test
        :param cmdline_param: list of strings, containing cmdline args, for example ['-inp' , 'simple_inp']
        :param exclude: list of strings file names to exclude from copy
        :param only_copy: list of string file names, or length 2 to change file name. example ['inp', ['kpts2', 'kpts']]
        in which the file 'kpts2' in the source dir will be renamed to kpts in the destination dir.
        :param rm_files: list of strings files in the workdir to be removed, will be executed before copy

        :return: a dictionary of the form 'filename :filepath'
        """
        import subprocess

        workdir = str(work_dir)
        testdir = test_dir()
        if cmdline_param is None:
            cmdline_param = []


        # Prepare only copy list, since we allow for name changes.

        new_only_copy_list = {}
        for entry in only_copy:
            if isinstance(entry, list):
                new_only_copy_list[entry[0]] = entry[1]
            else:
                new_only_copy_list[entry] = entry



        if test_file_folder is not None:
            abspath = os.path.abspath(os.path.join(testdir, test_file_folder))
            source = os.listdir(abspath)
            for files in source:
                if files not in exclude:
                    if new_only_copy_list != {}:
                        if files in list(new_only_copy_list.keys()):
                            #os.remove(os.path.abspath(os.path.join(workdir, new_only_copy_list[files])))
                            shutil.copy(os.path.abspath(os.path.join(abspath, new_only_copy_list[files])), workdir)
                    else:
                        shutil.copy(os.path.abspath(os.path.join(abspath, files)), workdir)
        #print(inpgen_binary)
        if inpgen_binary is None:
            print('No Inpgen binary found')
            return {}
        arg_list = [inpgen_binary] + cmdline_param
        #print(arg_list)

        os.chdir(workdir)
        with open(f"{workdir}/stdout", "w") as f_stdout:
            with open(f"{workdir}/stderr", "w") as f_stderr:
                p1 = subprocess.run(arg_list + ["-no_send"], stdout=f_stdout, stderr=f_stderr)
        # Check per hand if successful:
        with open(f"{workdir}/stderr", "r") as f_stderr:
                error_content= f_stderr.read()

        if 'Run finished successfully' not in error_content:
            # failure
            print('Inpgen execution failed.')
            print('======================================')
            print('The following was printed to stdout:')
            with open(f"{workdir}/stdout", "r") as f_stdout:
                print(f_stdout.read())
            print('======================================')
            print('The following was printed to stderr:')
            print(error_content)
            raise RuntimeError('Inpgen Execution failed')

        result_files = {}
        for root, dirs, files in os.walk(workdir):
            for file in files:
                rel_path = os.path.relpath(os.path.join(root, file), workdir)
                rel_path = rel_path.lstrip('./')
                abs_path = os.path.abspath(os.path.join(root, file))
                result_files[rel_path] = abs_path

        os.chdir(testdir)

        return result_files

    return _execute_inpgen

def get_mpi_command(env, mpi_procs, parallel):
    """
    Function returning the mpirun command for fleur with the given number of processes.
    The environment variables ``juDFT_MPI`` and ``juDFT_NPROCS`` can be used to customize it
    to the environment

    ``juDFT_MPI`` should not define the number of mpi processes. Either the task should be
    added at the end without inserting the number or the text ``{mpi_procs}`` should be added instead
    of the number

    :param env: dictionary with the defined environment variables
    :param mpi_procs: number of mpi processes to run
    :param parallel: boolean, if True the fleur executable is assumed to be fleur_MPI

    :returns: list of string with the mpi command to prepend the fleur call
    """
    import string
    import warnings

    mpiruncmd = env.get('juDFT_MPI', None)

    if env.get('juDFT_NPROCS', ''):
        mpi_procs = env.get('juDFT_NPROCS', '')

    if mpiruncmd is None and parallel:
        mpiruncmd = 'mpirun -n {mpi_procs} '

    if mpiruncmd is not None:
        if mpiruncmd.strip() != 'time' and len(mpiruncmd.strip()) > 0:
            if len([val[0] for val in string.Formatter().parse(mpiruncmd)]) != 0:
                try:
                    mpiruncmd = mpiruncmd.format(mpi_procs=mpi_procs)
                except (ValueError, KeyError) as exc:
                    raise KeyError("mpirun command could not be constructed: {}".format(exc))

            else:
                warnings.warn('The number of mpi processes will be appended to the mpicommand:\n {}'
                              'If this should not happen enter the text {{mpi_procs}} into the'
                              'command at the right place. For overriding the number of mpi processes use "juDFT_NPROCS"'.format(mpiruncmd))

                mpiruncmd += ' {} '.format(mpi_procs)

        mpiruncmd = mpiruncmd.split()
    else:
        mpiruncmd = []

    return mpiruncmd


@pytest.fixture
def mpi_command():
    """
    Fixture which returns the mpi commadn to run fleur
    """
    def _mpi_command(env, mpi_procs, parallel):
        """
        Function returning the mpirun command for fleur with the given number of processes
        The environment variables ``juDFT_MPI`` and ``juDFT_NPROCS`` can be used to customize it
        to the environment

        ``juDFT_MPI`` should not define the number of mpi processes. Either the task should be
        added at the end without inserting the number or the text ``{mpi_procs}`` should be added instead
        of the number

        :param env: dictionary with the defined environment variables
        :param mpi_procs: number of mpi processes to run
        :param parallel: boolean, if True the fleur executable is assumed to be fleur_MPI

        :returns: list of string with the mpi command to prepend the fleur call
        """
        return get_mpi_command(env, mpi_procs, parallel)

    return _mpi_command

@pytest.fixture
def execute_fleur(fleur_binary, work_dir, mpi_command):
    """
    Fixture which returns an execute_fleur function
    """
    def _execute_fleur(test_file_folder=None, cmdline_param=None, exclude=[], only_copy=[], rm_files=[], env={}, stderr='stderr', stdout='stdout', sub_dir=None, mpi_procs=2):
        """
        Function which copies the input files
        executes fleur with the given cmdline_param
        and returns a path to the folder with the results

        :param test_file_folder: string relative path to a folder containing all files, to copy for the test
        :param cmdline_param: list of strings, containing cmdline args, for example ['-inp' , 'simple_inp']
        :param exclude: list of strings file names to exclude from copy
        :param only_copy: list of string file names, or length 2 to change file name. example ['inp.xml', ['kpts2', 'kpts']]
        in which the file 'kpts2' in the source dir will be renamed to kpts in the destination dir.
        :param rm_files: list of strings files in the workdir to be removed, will be executed before copy

        :return: a dictionary of the form 'filename :filepath'
        """
        import subprocess

        if sub_dir is not None:
            workdir = os.path.abspath(os.path.join(work_dir, sub_dir))
        else:
            workdir = work_dir
        testdir = test_dir()
        if cmdline_param is None:
            cmdline_param = []

        osenv = dict(os.environ)
        run_env = osenv # This creates a complete env for the fleur, but we keep all session env
        # Not sure if this is good
        run_env['OMP_NUM_THREADS'] = osenv.get('OMP_NUM_THREADS', '1')
        run_env.update(env) # apply custom user changes

        # Prepare only copy list, since we allow for name changes.
        # but each file name is allowed only once
        new_only_copy_list = {}
        for entry in only_copy:
            if isinstance(entry, list):
                new_only_copy_list[entry[0]] = entry[1]
            else:
                new_only_copy_list[entry] = entry

        files_work_dir = os.listdir(workdir)
        for entry in rm_files:
            path = os.path.abspath(os.path.join(workdir, entry))
            if os.path.isfile(path):
                os.remove(path)
            else: # Either the file does not exits, or an expression was given
                for filename in files_work_dir:
                    if re.search(entry, filename):
                        path = os.path.abspath(os.path.join(workdir, filename))
                        if os.path.isfile(path):
                             os.remove(path)

        if test_file_folder is not None: # Does it even make sense to not give a folder?
            abspath = os.path.abspath(os.path.join(testdir, test_file_folder))
            source = os.listdir(abspath)
            for files in source:
                if files not in exclude:
                    if new_only_copy_list != {}:
                        if files in list(new_only_copy_list.keys()):
                            srcpath = os.path.abspath(os.path.join(abspath, files))
                            destpath = os.path.abspath(os.path.join(workdir, new_only_copy_list[files]))
                            shutil.copy(srcpath, destpath)
                    else:
                        shutil.copy(os.path.abspath(os.path.join(abspath, files)), workdir)

        fleur, parallel = fleur_binary
        if fleur is None:
            print('No Fleur binary found')
            return {}
        mpiruncmd = mpi_command(run_env, mpi_procs, parallel)

        arg_list = mpiruncmd + [fleur] + cmdline_param
        #print(arg_string)
        os.chdir(workdir)
        #t0 = time.perf_counter()
        with open(f"{workdir}/{stdout}", "bw") as f_stdout:
            with open(f"{workdir}/{stderr}", "bw") as f_stderr:
                # we parse the whole string and execute in shell,
                # otherwise popen thinks 'mpirun -np 2 /path/fleur' is the path to the executable...
                p1 = subprocess.run(arg_list, env=run_env, stdout=f_stdout, stderr=f_stderr)#check=True
        
        # Check per hand if successful:
        with open(f"{workdir}/stderr", "r") as f_stderr:
                error_content= f_stderr.read()

        if 'Run finished successfully' not in error_content:
            # failure
            print('Fleur execution failed.')
            print('======================================')
            print('The following was printed to stdout:')
            with open(f"{workdir}/stdout", "r") as f_stdout:
                print(f_stdout.read())
            print('======================================')
            print('The following was printed to stderr:')
            print(error_content)
            raise RuntimeError('Fleur Execution failed')

        result_files = {}
        for root, dirs, files in os.walk(workdir):
            for file in files:
                rel_path = os.path.relpath(os.path.join(root, file), workdir)
                rel_path = rel_path.lstrip('./')
                abs_path = os.path.abspath(os.path.join(root, file))
                result_files[rel_path] = abs_path
        os.chdir(testdir)

        return result_files

    return _execute_fleur

# Consider maybe running this after each fleur execution?
@pytest.fixture(scope='function')#, autouse=True) # make this available in every test
def validate_out_xml_file(execute_fleur):
    """
    return function validate_out_xml_file_f
    """
    def _validate_out_xml_file(file_path, schema_path=None):
        """
        Validates and outxml file via a fleur execution
        Maybe we also want to validate the out.xml file outside of fleur with python and lxml instead?
        Which is probably faster.
        So far we stay to the validatation test. Still one should test the fleur validatation feature at least once
        """
        import subprocess
        import shutil

        if 'out.xml' not in file_path:
            raise ValueError('No out.xml file given for validation.')
        root = file_path.split('out.xml')[0]
        if schema_path is None:
            schema_path = os.path.join(root,'FleurOutputSchema.xsd')
        print(f"Test validating outputfile: {file_path}")
        if not os.path.isfile(schema_path):
            msg = "No OutputSchema found"
            # the original test just continued
            print(msg)
            return True
            #raise ValueError(msg)
        # this fails as validation fails
        xmllint = shutil.which('xmllint')
        if xmllint is None:
            msg = "No xmllint executable found"
            # the original test just continued
            print(msg)
            return True

        with open(f"{root}/xmllintOut", "bw") as f_stdout:
            with open(f"{root}/xmllintErrors", "bw") as f_stderr:
                try:
                    subprocess.run([xmllint, '--schema', f'{schema_path}', f'{file_path}'],
                                   stderr=f_stderr, stdout=f_stdout, check=True)
                except Exception as e:
                    print(f"Failed validating outputfile: {e}")
                    return False
        return True

    return _validate_out_xml_file

# Comment: Instead of implementing grep in python one could also just execute grep via subprocess
# This might be rather unsave someone not nice could put 'grep x y; rm -rf /' in a test...
# The same goes for fleur and inpgen execute
# Pro: The rexpressions stay close that what people are used to.
@pytest.fixture
def grep_exists():
    """returns the grep_exits function
    """
    def _grep_exists(filepath, expression):
        """for an expression in a file
        Args:
            filepath ([str, path]): path to the file to search for
            expression (str): 'string' to look for

        :return: Bool, if exists
        """
        #t0 = time.perf_counter()
        exists = False
        #regex_pattern = "(" + expression + ")"
        #pattern = re.compile(regex_pattern) # This might be unsave,
        # but we do it to allow for same expressions as grep does
        with open(filepath, "r") as file1:
            for line in file1:
                if re.search(expression, line):#pattern.search(line):
                    # re.search allows also for patters, complains about some
                    #print(line)
                    exists = True
                    break
        #t1 = time.perf_counter()
        #print(f'Executing grep exits took {t1 - t0:0.4f} seconds')

        return exists

    return _grep_exists


@pytest.fixture
def grep_number():
    """returns the grep number function
    """
    def _grep_number(filepath, expression, split=None, line_index=1, res_index=-1):
        """Implements grep for a float number in a file

        :param filepath (str): path to the dile
        :param expression (str): a python expression to look for in lines (everything that is can be used by re.search())
        :param split (str, optional): if the expression is in line, split for this expression. if
None then the given expression will be used.
        :param line_index (int, optional): After the split where to look for the number, Default 1, therefore the first number ofter the split string
        :param res_index (int, optional): If expression matches several lines, which one to use. Defaults to -1, last number only
        if res_index=None a list of all numbers is returned

        :return: float, list of floats
        """
        #t0 = time.perf_counter()
        numbers = []
        with open(filepath, "r") as file1:
            for line in file1:
                if re.search(expression, line):
                    if split is not None:
                        res = line.split(split)[line_index]
                    else:
                        res = line.split(expression)[line_index]
                    try:
                        number = float(res)
                    except ValueError as exc: # There is still something after the number
                        number = float(re.findall(r"[+-]?\d+\.\d+", res)[0])
                    numbers.append(number)
        #t1 = time.perf_counter()
        #print(f'Executing grep number took {t1 - t0:0.4f} seconds')

        if len(numbers) == 0:
            raise ValueError(f'Number for "{expression}" was not found in {filepath}')
        elif len(numbers) == 1:
            return numbers[0]

        if res_index is not None:
            return numbers[res_index]
        else:
            return numbers
    return _grep_number

@pytest.fixture
def check_value_outfile():
    """Fixture which returns a check_value_outfile function
    """
    def _check_value_outfile(filepath, before_str, after_str, expected, delta):
        """
        Check a value in the out file, copied from libtest.py
        """
        exp_idx = 0
        found = False
        passed = False
        errors = 0
        with open(filepath, "r") as f:
            for line in f.readlines():
                if(before_str in line):
                    value_string = line.split(before_str)[-1]
                    value_string = value_string.split(after_str)
                    # remove empty strings
                    value_string = [i for i in value_string if i != ""][0]
                    value = float(value_string)
                    if (expected[exp_idx] is not None):
                        if(abs(value - expected[exp_idx]) < delta):
                            #log_info(f"PASSED [{exp_idx}]: {before_str} found: {value} expected: {expected[exp_idx]}")
                            passed = True
                        else:
                            #log_info(f"FAILED [{exp_idx}]: {before_str} found: {value} expected: {expected[exp_idx]}")
                            errors += 1
                    exp_idx = exp_idx + 1
                    found = True
                    if exp_idx > len(expected):
                        break # failed, so we stop
        if(len(expected) != exp_idx):
            #log_error("Number of expected values disagree with found values.")
            #log_error(f"before_str: '{before_str}''")
            passed = False
        if(errors != 0):
            passed = False
        return found, passed
    return _check_value_outfile

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
def clean_workdir(work_dir):
    """
    Fixture which returns a function to delete files in the test work dir
    the default is to delete all
    """
    def _clean_workdir(filelist=['all']):
        """
        Delete files in the test work dir the default is to delete all
        :param filelist (list, optional): filesnames to delete. Defaults to ['all'].
        """
        workdir = work_dir
        if 'all' in filelist:
            for entry in os.listdir(workdir):
                path = os.path.abspath(os.path.join(workdir, entry))
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)
        else:
            for file in filelist:
                os.remove(os.path.abspath(os.path.join(workdir, file)))

    return _clean_workdir

@pytest.fixture(scope='function')
def stage_workdir(work_dir, stage_dir, request):
    """
    A fixture when used will cause a test to copy the workdir content to the stage folder.
    Needed for tests which depend on other tests. Or can be used during development to shorten tests
    """
    def _stage_workdir(foldername=None, index=0):
        """
        A method causing a test to copy the workdir content to the stage folder.
        Needed for tests which depend on other tests. Or can be used during development to shorten tests        """
        if foldername is None:
            method_name = request.node.name + f'_{index}'
            foldername = os.path.abspath(os.path.join(stage_dir, method_name))
        else:
            foldername = os.path.abspath(os.path.join(work_dir, foldername))
        if os.path.isdir(foldername):
            shutil.rmtree(foldername)
        shutil.copytree(work_dir, foldername)

    return _stage_workdir

@pytest.fixture(scope='function')
def load_stage(workdir, clean_workdir, stage_dir, request):
    """
    A fixture when used will cause a test to copy all content of a certain dir in the stage folder to the
    workdir. This will overwrite all files in the current workdir. Needed for tests which depend on other tests.
    """
    def _use_stage(foldername=None, index=0):
        """
        A method to copy all content of a certain dir in the stage folder to the
        workdir. This will overwrite all files in the current workdir. Needed for tests which depend on other tests.
        """
        if foldername is None:
            method_name = request.node.name + f'_{index}'
            foldername = os.path.abspath(os.path.join(stage_dir, method_name))
        #if not os.path.isdir(foldername):
        #    #throw error
        shutil.copytree(foldername, workdir)

    return _use_stage

@pytest.fixture(scope='function', autouse=True)
def stage_for_parser_test(request, work_dir, parser_testdir):
    """
    Fixture to copy the result files of a test to the parser test folder.
    To later autogenerate tests for the fleur parsers
    """

    if 'fleur_parser' not in request.keywords:

        parsertestdir = parser_testdir
        workdir = work_dir

        yield # test is running

        # clean up code goes here:

        method_name = request.node.name
        if RUN_PARSER_TESTS and not request.node.rep_call.failed:
            # if not failed move test result to parsertestdir, will replace dir if existent
            destination = os.path.abspath(os.path.join(parsertestdir, method_name))
            if os.path.isdir(destination):
                shutil.rmtree(destination)
            #os.mkdir(destination)
            shutil.copytree(workdir, destination)
    else:
        yield

def get_inpgen_binary(fleur_dir):
    """
    Fixture returning the path to a inpgen executable
    """
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
def inpgen_binary(build_dir):
    """
    Fixture returning the path to a inpgen executable
    """
    return get_inpgen_binary(build_dir)


def get_fleur_binary(fleur_dir):
    """
    Fixture returning the path to a fleur executable
    """
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

@pytest.fixture(scope='session')
def fleur_binary(build_dir):
    """
    Fixture returning the path to a fleur executable
    """
    return get_fleur_binary(build_dir)


@pytest.fixture
def inpxml_etree():
    """Returns the etree generator"""
    def _get_etree(filepath):
        from lxml import etree
        with open(filepath, 'r') as inpxmlfile:
            tree = etree.parse(inpxmlfile)
        return tree
    return _get_etree



@pytest.fixture()
def collect_all_judft_messages_f():
    return collect_all_judft_messages

def collect_all_judft_messages():
    """Helper function, to create a list of all judft messages within the fleur source code,
    by grepping for them.
    """
    testdir = test_dir()
    rel_fleur_source = '../../'
    fleur_source_dir = os.path.abspath(os.path.join(testdir, rel_fleur_source))
    # source code is top dir, to much other stuff in there, thats why hardcode source dirs for speed
    # and to avoid problems with binaries and so on.
    src_folders = ['cdn', 'cdn_mt', 'core', 'diagonalization', 'dos', 'eels', 'eigen',
    'eigen_secvar', 'eigen_soc', 'external', 'fermi', 'fft', 'fleurinput', 'force',
    'forcetheorem', 'global', 'greensf', 'hybrid', 'include', 'init', 'inpgen2',
    'io', 'juDFT', 'kpoints', 'ldahia', 'ldau', 'main', 'math', 'mix', 'mpi', 'optional', 'orbdep',
    'rdmft', 'tetra', 'types', 'vgen', 'wannier', 'xc-pot'
     ]

    grep_results = []
    grep_string = '(judft_error|error_output)'
    # fortran is not case sensitive, sometimes output is programmed line before.
    # maybe use real grep instead of this python implementation...
    for folder in src_folders:
        folder_path = os.path.join(fleur_source_dir, folder)
        for root, dirs, files in os.walk(folder_path):
            for filename in files:
                if not (filename.endswith('f90') or filename.endswith('F90')):
                    continue
                file_path = os.path.join(root, filename)
                try:
                    with open(file_path, encoding="utf-8") as file1:
                        for line in file1:
                            try:
                                if re.search(
                                    grep_string, line, re.IGNORECASE):
                                    grep_results.append(line)
                            except (TypeError):
                                pass
                except (IOError, OSError, UnicodeDecodeError):
                    pass
    # Construct list
    all_messages = []
    # there are all combinations of strings all over the place in the fleur src
    for judft_string in grep_results:
        ju_str = judft_string.split('"')
        if len(ju_str) == 1:
            ju_str = judft_string.split("'")
            if len(ju_str) == 1:
                # ignore this one
                continue
        string = ''
        for s in judft_string:
            string = string + s
        # if calledby is used the split is always different, therefore we split calledby.
        if re.search('calledby', string):
            ju_str1 = string.split('calledby')[0]
            ju_str = ju_str1.split('"')
            if len(ju_str) == 1:
                ju_str = ju_str1.split("'")
                if len(ju_str) == 1:
                    # ignore this one
                    continue
            try:
                message = ju_str[-2]
            except IndexError: # In the source code both strings sings are used
                # We are missing some, but currently we do not care
                print(ju_str)
                message = None
                continue

        else:
            message = ju_str[-2]
        all_messages.append(message)
    return list(set(all_messages))
