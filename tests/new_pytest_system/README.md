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
