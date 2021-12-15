# -*- coding: utf-8 -*-
MASCI_TOOLS_VERSION_STR = 'Not installed'
try:
    import masci_tools
    MASCI_TOOLS_VERSION_STR = masci_tools.__version__
    version = tuple(int(val) for val in MASCI_TOOLS_VERSION_STR.replace('-','.').split('.')[:3])
except ImportError:
    RUN_PARSER_TESTS = False
except AttributeError:
    MASCI_TOOLS_VERSION_STR = 'too old'
    RUN_PARSER_TESTS = False
else:
    if version >= (0,4,0):
        RUN_PARSER_TESTS = True
