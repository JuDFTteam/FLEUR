**Notice:** This only works with **python version >3.5** and **NOT with python2**. If run with python2 you get an error message stating that you pytest version is prior 5.0.

# The (python) test system explained
The python test system relies on pytest (https://docs.pytest.org/en/stable/), which is a standard for python and per default installed on most our machines.
Also looking at:
```
$pytest --help
```
may provide you with some basic understanding. 
If the `pytest` command is not available try out:
```
$python3 -m pytest masci-tools
```
instead, maybe pytest is already installed within python and you can use it as a module. If this is the case you have to replace `pytest` in all further commands with this.
If this command fails, you might need to install pytest using:
```
pip install -U pytest
```   

# How to run tests:
To run the tests pytest has to know two things: where to find the `conftest.py` file, which is under `<fleur_source_folder>/tests/new_pytest_system` and which build dir to test for.

## Execution from tests source dir
One way is to go to the `tests/new_pytest_system` folder and to run all tests execute:
```
$pytest
```
(default assumed build dir is `<fleur_source_folder>/build`)
pytest automatically discovers all tests in the any sub directory structure of the head `conftest.py` file, by looking for certain identifiers per default, like files and functions starting with `test_` or ending with `_test` (https://docs.pytest.org/en/stable/goodpractices.html#test-discovery) . During execution pytest will print a standard report for the test run.

## Execution from the build dir
If you execute pytest from any other directory, you have to provide the path to the folder with the head `conftest.py` and a path to the build dir with the fleur executables, either relative to the `contest.py` file or absolute.
For example if you build dir is `<fleur_source_folder>/build.123`
```
$pytest ../tests/new_pytest_system --build_dir=../../build.123
```
or
```
$pytest ../tests/new_pytest_system --build_dir=`<fleur_source_folder>/build.123`
```

The build process also generates a shell script in the build folder to perform the right pytest command for the current build directory. The arguments to this script are forwarded to pytest.
```
./run_tests.sh
```
If you need a non-default python executable for the tests, you can specify this for the `run_tests.sh` script by setting the `juDFT_PYTHON` environment variable with the right executable


## Executing a subset off tests

Pytest ways to run a subset of the test suite:

* Select tests to run based on substring matching of test names.
* Select tests groups to run based on the markers applied.
```
pytest -k <substring> -v
pytest -m <markername>
```
example markers are:
```
bulk,film,xml,collinear,soc,lo,ldau,non-collinear,spinspiral,forces,hybrid,slow,fast... 
```
You can register markers in the `conftest.py` file. Also pytest will print a warning if a it does not understand a certain marker.
 
To not execute all tests but instead exit at first or xth failure, execute:
```
$pytest -x           # exit after first failure
$pytest --maxfail=2  # exit after two failures
```

Per default pytest times all tests but only displays the total runtime in the report.
adding the `–duration=N` option will print the times of the `N` slowest tests during the run. 

To execute all tests in a given (sub)folder (`tests/bar`) or file `tests/bar/test_foo.py` run (will only work from the `tests/new_pytest_system` ):
```
$pytest tests/bar
$pytest tests/bar/test_foo.py
```
# What is added to the build dir:
Under `Testing/` a test session run will create several folders:
- `work`: here we run the tests, i.e execute inpgen and fleur. This folder is cleaned after every test.
- `failed_test_results`: a folder where for all failed tests the content of `work` is preserved of the last test session. This folder is cleaned at the beginning of a test session.
- `parser_testdir`: a folder where all scheduled parser tests go, this is cleaned before each test session. Soon this is cleaned during a session, and failed parser tests are also moved to `failed_test_results`
- (only on CI) `pytest_session.stdout`, `pytest_summmary.out`: On the CI we also write the pytest report and the short summary to files (using tee, check the ci.yml file)

Also cmake writes out a file for pytests with information on the fleur build, to automaticly exclude tests for specific fleur build, i.e specifc libraries like magma, libxc, ...

# How to read the test log:

Pytest writes a log for each test session, which is usually outputed into the terminal. We also pipe (via `| tee <filepath>`) it to the `pytest_session.stdout` file on the CI.

A example (shorted) log from a test session run with line numbers may look like this (usually colored):
```
 0 $ ./run_tests.sh | tee $CI_PROJECT_DIR/build/Testing/pytest_session.stdout
 1 ============================= test session starts ==============================
 2 platform linux -- Python 3.8.5, pytest-6.2.4, py-1.10.0, pluggy-0.13.1
 3 ------------------------------ Fleur test session ------------------------------
 4 Fleur exe: /builds/fleur/build/fleur_MPI
 5 Inpgen exe: /builds/fleur/build/inpgen
 6 NOT linked libraries: []
 7 Running tests in: /builds/fleur/build/Testing/work
 8 Failed tests will be copied to: /builds/fleur/build/Testing/failed_test_results
 9 Parser tests will be run (masci-tools version 0.4.8)
10 Default MPI command: mpirun -n {mpi_procs} --allow-run-as-root --mca btl vader,self
11 Now cleaning work, failed and parser_test directories...
12 rootdir: /builds/fleur/tests/new_pytest_system, configfile: setup.cfg
13 Excluding tests with the following markers in 'pytest_incl.py':  ['chase', 'cusolver', 'edsolver', 'elpa', 'elpaonenode', 'fftmkl', 'gpu', 'magma', 'noci', 'progthread', 'spfft', 'wannier4', 'wannier5']
14 Running every 1st test with offset 0, others will be skiped.
15 collected 210 items / 4 deselected / 206 selected
16 ../tests/feature_reg/test_AlLibxcPbe.py .         [  0%]
17 ../tests/feature_reg/test_Co.py ..                [  1%]
18 ../tests/feature_reg/test_CrystalFieldOutput.py . [  1%]
19 ../tests/feature_reg/test_CuBulk.py .....         [  4%]
20-80.....
81 ../tests/libxc/test_libx.py s                     [ 31%]
82../tests/masci_tools/test_banddos_parser.py ....   [ 33%]
83 ../tests/masci_tools/test_fleur_parser.py ..s... [ 38%]
84 s...s..........s.......s............s...s...........s...ss.....s..s....s [ 73%]
85 ........s.......s............s...s.........F.s...ss...                   [ 99%]
86 ../tests/new_pytest_system/tests/masci_tools/test_judft_errors.py .      [100%]
87 =========================== short test summary info ============================
88 FAILED ../test_fleur_mt_outxml_parser[test_CwannXML]
89 ===== 1 failed, 178 passed, 27 skipped, 4 deselected in 1408.04s (0:23:28) =====
```
## Log: Test session header
The log starts with a session header, going in this case from line 1 to 15.
In this header contains default pytest output of versions (line 1), and the output of how many tests pytest has discovered, selected and deselected (line 15).
Further the session header contains information we put there. I.e which fleur executable used (line 4), inpgen executable used (line 5), libraries not linked (line 6), in which folder the tests will run (line 7), where files from failed tests will be copied to (line 8), if parser tests will be run, and for which masci-tools version (line 9), the default mpi command to execute fleur (line 10), at the start of the sessions these folders are cleared and/or created (mentioned on line 11), the rootdir (all other paths below are relative to this) and the pytest configfile (line 12).  Line 13 list all markers cmake has written into `pytest_incl.py` within the build dir for compilation related test deselection. Line 14 states if only every x test is run and if the test session has an offset.

## Log: Running tests
Then in the lines 16-86 information on the running tests is outputed. Each line contains information on tests from which file are run, followed by a '.' (dot) for each passed test. Skipped tests are marked with an 's', failed tests with and 'F' and tests which had unexpected errors with an 'E'.
The line ends with a procentage, for how far we are into the test session.

## Log: Stacktraces, understanding what failed
After the running tests (not shown in this example) follow long stack traces containing information on were the test failed, and what was recorded in stderr.
This starts after a line containing ```=========== FAILURES ========```.
Further some examples of stack traces explained (also from different test session):
### Example 1
A full long stackstrace of a parser test failure:
```
______________________________________ test_fleur_mt_outxml_parser[test_CwannXML] ______________________________________

request = <FixtureRequest for <Function test_fleur_mt_outxml_parser[test_CwannXML]>>, fleur_test_name = 'test_CwannXML'
test_file = 'tests/feature_reg/test_Cwann.py', parser_testdir = '/builds/fleur/build/Testing/parser_testdir'

    @pytest.mark.fleur_parser
    @pytest.mark.masci_tools
    def test_fleur_mt_outxml_parser(request, fleur_test_name, test_file, parser_testdir):
        """
        For each folder (parametrization happens in conftest.py) in parser test,
        try if the fleur parser in masci-tools can handle the parsing without
        crashing, successful and an empty parser log
        """
        #Note:
        #   You might notice that a lot of the output parser tests fail either because the out file does not validate
        #   Or there are warnings
        #   1. The validation errors occur since the schemas in fleur (develop) and masci-tools are slightly out of sync
        #      at the moment (fixed in the raise_fleur_file_version branch)
        #   2. There are a couple of output differences not yet accounted for in the output parser (Some could be solved in fleur)
        #           - bandgap output only in hist mode
    
        #These warnings are expected to appear at the moment
        KNOWN_WARNINGS = {'No values found for attribute l_f'}
    
        pytest.importorskip('masci_tools',minversion='0.4.0-dev3')
        from masci_tools.io.parsers.fleur import outxml_parser
        depends(request, [f'{test_file}::{fleur_test_name}'], scope='session')
    
        outxmlfilepath = os.path.abspath(os.path.join(parser_testdir, f'./{fleur_test_name}/', 'out.xml'))
    
        assert os.path.isfile(outxmlfilepath)
    
        parser_info = {}
        out_dict = outxml_parser(outxmlfilepath, parser_info_out=parser_info)
    
        if any("Schema available for version" in warning for warning in parser_info['parser_warnings']):
            for warning in parser_info['parser_warnings'].copy():
                if 'Output file does not validate against the schema' in warning:
                    parser_info['parser_warnings'].remove(warning)
                if 'Schema available for version' in warning:
                    parser_info['parser_warnings'].remove(warning)
    
        assert out_dict is not None
        assert isinstance(out_dict, dict)
        assert out_dict != {}
        assert parser_info['parser_errors'] == []
        assert parser_info['parser_critical'] == []
>       assert set(parser_info['parser_warnings']).difference(KNOWN_WARNINGS) == set() #At the moment there is always at least one warning
E       assert {"[Iteration ...tribute", ...} == set()
E         Extra items in the left set:
E         '[Iteration 1] No values found for attribute distance'
E         '[Iteration 1] No values found for attribute interstitial at tag spinDependentCharge'
E         '[Iteration 1] No values found for attribute value at tag sumOfEigenvalues'
E         "[Iteration 1] Failed to evaluate singleValue from tag chargeDenXCDenIntegral: Has no 'value' attribute"
E         '[Iteration 1] No values found for attribute value at tag valenceElectrons'
E         '[Iteration 1] No values found for attribute mtSpheres at tag spinDependentCharge'
E         "[Iteration 1] Failed to evaluate singleValue from tag totalCharge: Has no 'value' attribute"
E         '[Iteration 1] No values found for attribute total at tag spinDependentCharge'
E         '[Iteration 1] convert_total_energy cannot convert None to eV'
E         "[Iteration 1] Failed to evaluate singleValue from tag totalEnergy: Has no 'units' attribute"
E         '[Iteration 1] No values found for attribute value at tag totalCharge'
E         "[Iteration 1] Failed to evaluate singleValue from tag totalEnergy: Has no 'value' attribute"
E         '[Iteration 1] No values found for attribute value at tag chargeDenXCDenIntegral'
E         Use -v to get the full diff

/builds/fleur/tests/new_pytest_system/tests/masci_tools/test_fleur_parser.py:81: AssertionError
-------------------------------------------------- Captured log call ---------------------------------------------------
INFO     masci_tools.io.parsers.fleur.fleur_outxml_parser:fleur_outxml_parser.py:97 Masci-Tools Fleur out.xml Parser v0.5.0
WARNING  masci_tools.io.parsers.fleur.fleur_outxml_parser:schema_dict.py:395 No Output Schema available for version '0.35'; falling back to '0.34'
INFO     masci_tools.io.parsers.fleur.fleur_outxml_parser:fleur_outxml_parser.py:201 Found fleur out file with the versions out: 0.34; inp: 0.34
WARNING  masci_tools.io.parsers.fleur.fleur_outxml_parser:fleur_outxml_parser.py:212 Output file does not validate against the schema: 
Line 2: Element 'fleurOutput', attribute 'fleurOutputVersion': [facet 'enumeration'] The value '0.35' is not an element of the set {'0.34'}.
...
```
The first lines tells you from which test this is, below is information for which previous run test case this parser test is from (in this case `test_CwannXML`). This information is followed by a print of the complete python code of the test until the line, where the first exception was thrown (indicated by the `>` on the start of the line). This way one sees, what is run, what parsed already.
In the lines following the line starting with the `>` comes the probably most important information often telling us what exactly went wrong. In this case we tested that `assert set(parser_info['parser_warnings']).difference(KNOWN_WARNINGS) == set()`, i.e we throw an `AssertionError` if there are any unsupected parser warnings. The next line shows us what the values on each side of the `==` were. Then at the end of the report any caputered output or captured logging is also printed.
### Example 2
Most important lines of a stacktrace of a tests which failed due to a value not beeing as tested for:
```
 tenergy = grep_number(res_files['out'], "total energy=", "=")
>       assert abs(tenergy - -1270.4886) <= 0.0001
E       assert 7.6435286517000804 <= 0.0001
E        +  where 7.6435286517000804 = abs((-1278.1321286517 - -1270.4886))
/builds/fleur/tests/new_pytest_system/tests/feature_reg/test_Noncollinear_downward_comp.py:28: AssertionError
```
This means we tested for that the total energy is -1270.4886 but instead it turned out to be -1278.1321286517.

### Example 3
Often we grep in files and expect a certain expression to be there, if it is not this will look in a stack trace like that:
```
>       assert grep_exists(res_files['out'], "it=  1  is completed")
E       AssertionError: assert False
E        +  where False = <function grep_exists.<locals>._grep_exists at 0x7f3b7aba21f0>('/home/build/Testing/work/out', 'it=  1  is completed')

/tests/new_pytest_system/tests/feature_reg/test_CuBulk.py:26: AssertionError
```
here "it=  1  is completed" is not in the "out" file.

### Example 4
A failure of a fleur execution.

## Log: Test session Summary
At the end of the test session, a short summary of the test session is outputed (line 87-89).
It contains a single line (line 88 in this case) for each failed tests with the error message (truncated to terminal width).
And a line summing up how many tests failed, passed, where skipped, deselected, errored and how long the full session was.

# Tests (source) folder layout:

It makes sense to organize and group tests already in a directory layout, which will also make it easier to execute a certain test subset. Besides folders with tests there are other folders and files in the top directory which are important:

- `tests`: here and under the sub-folders go all further test folders, tests are organized through test type and the functionality they test
- `inputfiles`: To separate all tests input files and data files from the test code, we put them in this folder. The folder names do not matter, through it would be nice if the name would tell to which test(s) the folder belongs to. In the future we might have test which have (autogenerated) files which they test against, and we want to keep the input files separated from them.
- `helpers`: Folder where one can place python code which will be in the python path for a pytest run and can therefore be imported for test code.
- `conftest.py`: This is a (are) central pytest file(s), here all fixtures i.e helpers go (what former has partly been taken care of by the 'libtest.py' file or the 'scripts/*') which can be used within the tests, or automatically do something before or after a test run. For more on fixtures read the pytest basics section below or take a look at (https://docs.pytest.org/en/stable/fixture.html#fixture)


# Implementing new tests:

## Python basics:

### Decorators:

the stuff with the `@` above the functions are so called decorators, which are executed before the actual function they decorate is executed. So they wrap around the functions. Implementing it this way makes the code easier to read and understand. For example see:
https://pythonbasics.org/decorators/

## Pytest basics

https://docs.pytest.org/en/stable/

Fixtures are helper functions, which can be put in the `.py` file of the tests to be available there or in the conftest.py file(s) to be available globally.
We use them for preparation or teardown code, also basicly to make anything available we need throughout the different tests.
Fixtures with scope `session` are executed one a test session level (for example the binaries, because we use the same fleur exe for all tests)
Fixtures with scope `function` can be used within a test, and could with `autouse=True` automatically be executed for every test.
(for more on this look at the `conftest.py`file)

usefull pytest decorators
```
@pytest.mark.skip() # skips the test (shown in report)
@pytest.mark.skipif() # skips the test under a certain condition
@pytest.mark.parameterize() # creates a test case for each 'parameter'
@pytest.mark.<markername> # mark test if a given label to allow subset execution
```
### Pytest plugins we use:

there are a bunch of them, we use two so far
```
pytest-dependency
```
If you want to add a plugin, please add it its source code, that one does not have to install it to run the tests everywhere.
i.e add it under `tests/new_pytest_system/pytest_plugins`

# How to create a new test

As an example and good starting point look at the `tests/feature_reg/test_CuBulk.py` file which contains several tests on a Cu bulk system.

1. Check if the thing you want to test, is not covered by any other test. (Since fleur regression tests are slow). Maybe you can extend a further test, or add a stage, or depend on an other test on which to continue. 
1. Create a new folder under the `inputfiles` folder where you put all the input files needed for your test
2. Optional add a new `test_<some_name>.py` or `<some_name>_test.py` file somewhere under the `tests` directory where your test functions/code goes, or add you test function to an already existing file. Usually, one executes all tests within a file, so if you want to execute just your test alone, it needs its own file, or a special marker.
3. If the test function executes without throwing an `exception` pytest will mark it as `PASSED` during execution, Otherwise, the test will fail on the first error (python exception) thrown. Therefore, if you want to test for a simple comparison of something you can use the `assert` statement to do so, which will throw an `AssertionError` if what comes after `assert` is `False`.
examples:
```
assert ('out' in res_file_names), 'The out file is missing'
assert fermi == 0.4233
assert abs(fermi - 0.4233) <= 0.001
assert grep_exists(res_files['out'], "it=  1  is completed")
```
You can specify some error message behind the statement like in the first example above. Through this is not needed and sometimes less transparent since pytest per default prints a stack trace for failed tests, where you see which lines where executes and for example what the value of `fermi` from above turned out to be.

Any other exception will due of course, for example this is also fine
```
validate('inp.xml')  # this will throw some libxml errors if it does not validate
```

example test function explained:
```
@pytest.mark.bulk                  # markers to group tests
@pytest.mark.xml                   # markers to group tests
@pytest.mark.fast                  # markers to group tests
def test_CuBulkXML(execute_fleur, grep_exists, grep_number, stage_for_parser_test): # arguments are what fixures you want to explicitly include, in this case execute_fleur, grep_number and grep_exists
    """ # Describe what your test does here
    Simple test of Fleur with XML input with one step:
    1.Generate a starting density and run 1 iteration and compare fermi-energy & total energy
    """
    # here the python code for the test follows.
    
    # execute fleur needs to know the path to your test files.
    # a relative path to the folder with the top conftest.py file is fine
    test_file_folder = './inputfiles/CuBulkXML/' 

    res_files = execute_fleur(test_file_folder)   # execute fleur in work dir
    should_files = ['out.xml', 'out', 'usage.json']
    res_file_names = list(res_files.keys())
    
    # Test if all files are there
    for file1 in should_files:
        assert file1 in res_file_names
    
    # some more tests
    assert grep_exists(res_files['out'], "it=  1  is completed")

    # get some numbers from the out file
    fermi = grep_number(res_files['out'], "new fermi energy", ":")
    tenergy = grep_number(res_files['out'], "total energy=", "=")
    dist = grep_number(res_files['out'], "distance of charge densitie", "1:")
    
    # test if these number are what you expect
    assert abs(fermi - 0.4233) <= 0.001
    assert abs(tenergy - -3305.016) <= 0.001
    assert abs(dist - 45.8259) <= 0.001
```


There are certain things, which are still missing in the pytest system. Feel free to implement them or to open feature request issues that we can improve.
