RUN_PARSER_TESTS = True
try:
    import masci_tools
    # from masci_tools.parser import parse_fleur_outxml, parse_fleur_inpxml, parse_fleur_banddos, parse_fleur_nnm_mat
except ImportError:
    print('Masci-tools not installed. Skipping all parser tests and other test requireing it.')
    RUN_PARSER_TESTS = False
