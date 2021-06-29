
import pytest
import os
from _pytest.terminal import TerminalReporter
from _pytest.config import create_terminal_writer
from contextlib import contextmanager

@contextmanager
def replace_attribute(obj, attr_name, attr_value):

    old_value = getattr(obj, attr_name)

    setattr(obj, attr_name, attr_value)
    try:
        yield
    finally:
        setattr(obj, attr_name, old_value)


def pytest_runtest_logreport(report):
    report.nodeid = report.nodeid.split('::')[-1]


def pytest_addoption(parser):
    parser.addoption("--test-summary-file",
                     help="Redirect the long test summary to this file",
                     default=None, action="store")
    parser.addoption("--overwrite-terminal-width",
                     help="Overwrite the automatic terminal width for the TerminalWriter",
                     default=None, action="store")

class FleurTestsTerminalReporter(TerminalReporter):

    def __init__(self, config, file=None):

        super().__init__(config, file=file)

        try:
            self.startpath = self.config.rootpath
        except AttributeError:
            self.startpath = self.config.rootdir

        force_width = config.getoption("--overwrite-terminal-width")
        if force_width is not None:
            self._tw.fullwidth = int(force_width)

        long_test_summary_file = config.getoption("--test-summary-file")
        if long_test_summary_file is None:
            self._summary_tw = self._tw
        else:
            filepath = os.path.abspath(long_test_summary_file)
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            long_test_summary_file = open(filepath, 'w+')
            self._summary_tw = create_terminal_writer(config, long_test_summary_file)
            self._summary_tw.hasmarkup = False
            self._summary_tw.code_highlight = False
            self._summary_tw.fullwidth = 120


    @pytest.hookimpl(hookwrapper=True)
    def pytest_terminal_summary(self):
        with replace_attribute(self, '_tw', self._summary_tw):
            self.summary_errors()
            self.summary_failures()
            self.summary_warnings()
            self.summary_passes()
            self.flush()
        yield
        self.short_test_summary()
        # Display any extra warnings from teardown here (if any).
        self.summary_warnings()


@pytest.mark.trylast
def pytest_configure(config):
    vanilla_reporter = config.pluginmanager.getplugin("terminalreporter")
    my_reporter = FleurTestsTerminalReporter(config)
    config.pluginmanager.unregister(vanilla_reporter)
    config.pluginmanager.register(my_reporter, "terminalreporter")