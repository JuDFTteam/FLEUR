RUN_PARSER_TESTS = True
try:
    import masci_tools
    from masci_tools.io.parsers.fleur import inpxml_parser, outxml_parser
    # from masci_tools.parser import parse_fleur_banddos, parse_fleur_nnm_mat
except ImportError:
    print('Required masci-tools version not installed (>=0.4.0). Skipping all parser tests and other test requiring it.')
    RUN_PARSER_TESTS = False
