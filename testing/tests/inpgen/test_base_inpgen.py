"""
Collection of tests running only inpgen
"""

import pytest


@pytest.mark.inpgen
@pytest.mark.bulk
@pytest.mark.fast
def test_inpgen_Si_simple_interface_files(execute_inpgen, clean_workdir):
    """Simple test of Inpgen:
    1.Generate input files where no namelists are given & compare with expected files
    2.Generate input files with explicit option & compare with expected files
    """
    test_file_folder = './inputfiles/inpgen/Si_no_para/'

    # Test if inpgen runs and required files are there
    cmd_params = ['-f', 'inp_simple']
    res_files = execute_inpgen(test_file_folder, cmd_params)
    should_files = ['inp.xml', 'struct.xsf', 'default.econfig', 'out', 'usage.json', 'sym.xml', 'kpts.xml']
    res_file_names = list(res_files.keys())
    
    for file1 in should_files:
        assert file1 in res_file_names

    clean_workdir()

    # Test if inpgen can run with --explicit and different includes are right in inp.xml
    cmd_params = ['-explicit', '-inc', 'all', '-f', 'inp_simple']
    res_files = execute_inpgen(test_file_folder, cmd_params)
    should_files = ['inp.xml', 'struct.xsf', 'default.econfig', 'out', 'usage.json']
    not_there = ['sym.xml', 'kpts.xml']
    res_file_names = list(res_files.keys())
    for file1 in should_files:
        assert file1 in res_file_names
    for file1 in not_there:
        assert file1 not in res_file_names


@pytest.mark.inpgen
@pytest.mark.bulk
@pytest.mark.fast
def test_inpgen_Si_inputfile_content(execute_inpgen, clean_workdir):
    """Simple test of Inpgen:
    1.Generate input files where no namelists are given & compare with expected files
    Be very careful when changing the inpgen input file, it is parsed by various other tools:

    Masci-tools (used by aiida-fleur): https://github.com/JuDFTteam/masci-tools/tree/develop/masci_tools/io/io_fleur_inpgen.py
    Phonopy: https://github.com/phonopy/phonopy/blob/develop/phonopy/interface/fleur.py
    # ASE
    """
    test_file_folder = './inputfiles/inpgen/Si_full_para/'

    # Test if inpgen runs and required files are there
    cmd_params = ['-f', 'inp_simple']
    res_files = execute_inpgen(test_file_folder, cmd_params)
    should_files = ['inp.xml', 'struct.xsf', 'default.econfig', 'out', 'usage.json', 'sym.xml', 'kpts.xml']
    res_file_names = list(res_files.keys())
    
    # if it would complain about a list the inp.xml would not be there
    for file1 in should_files:
        assert file1 in res_file_names

    # Now we test if the namelists are really set in the inp.xml


@pytest.mark.inpgen
@pytest.mark.bulk
@pytest.mark.fast
def test_inpgen_wann_projection_counts(execute_inpgen, clean_workdir):
    """-wann ELEMENT:channels writes the trial orbitals and sizes the window.

    The two silicon atoms are one species with neq=2, so the count has to pick up
    the atoms of the species and not the species: s and p give (1+3)*2 = 8, and d
    alone gives 5*2 = 10. The window is twice that, so the disentanglement has
    states to choose from, and #nbands replaces it without touching the
    projections. minBand, maxBand and the numbands cutoff all follow the window.
    """
    folder = './inputfiles/inpgen/Si_no_para/'

    res = execute_inpgen(folder, ['-f', 'inp_simple', '-wann', 'Si:sp'])
    with open(res['inp.xml']) as f:
        content = f.read()
    assert content.count('<wannierproj l=') == 2
    assert '<wannierproj l="0" m="0" spin=""/>' in content
    assert '<wannierproj l="1" m="0" spin=""/>' in content
    assert 'numBands="16"' in content
    assert 'minBand="1" maxBand="16"' in content
    assert 'numbands="16"' in content          # the <cutoffs> follow maxBand
    clean_workdir()

    res = execute_inpgen(folder, ['-f', 'inp_simple', '-wann', 'Si:d'])
    with open(res['inp.xml']) as f:
        content = f.read()
    assert content.count('<wannierproj l=') == 1
    assert '<wannierproj l="2" m="0" spin=""/>' in content
    assert 'numBands="20"' in content
    clean_workdir()

    res = execute_inpgen(folder, ['-f', 'inp_simple', '-wann', 'Si:sp#nbands=24'])
    with open(res['inp.xml']) as f:
        content = f.read()
    assert content.count('<wannierproj l=') == 2
    assert 'numBands="24"' in content
    clean_workdir()

    # The bare number keeps its old meaning: a window and no projections.
    res = execute_inpgen(folder, ['-f', 'inp_simple', '-wann', '24'])
    with open(res['inp.xml']) as f:
        content = f.read()
    assert '<wannierproj l=' not in content
    assert 'numBands="24"' in content
    clean_workdir()

    # And without the flag at all. The species loop asks whether a projection
    # list was given, which has to hold for a run that never looked at one.
    res = execute_inpgen(folder, ['-f', 'inp_simple'])
    with open(res['inp.xml']) as f:
        content = f.read()
    assert '<wannierproj l=' not in content
    assert '<wannierlib' not in content

    clean_workdir()


@pytest.mark.inpgen
@pytest.mark.bulk
@pytest.mark.fast
def test_inpgen_wann_projections_reach_every_species(execute_inpgen, clean_workdir):
    """One key reaches every species of its element.

    In CoF2 the two cobalt sites are separate species, written by non-integer ids,
    while the four fluorines are one species with neq=4. A single Co key therefore
    has to land in both cobalt species, and the F key has to count four atoms:
    3 + 3 + 2 = 8 entries and 9 + 9 + 4*4 = 34 Wannier functions. That is what
    makes the count independent of how the species came out, which a single-species
    structure cannot show.
    """
    folder = './inputfiles/inpgen/CoF2_two_species/'

    res = execute_inpgen(folder, ['-f', 'inp_simple', '-wann', 'Co:spd,F:sp'])
    with open(res['inp.xml']) as f:
        content = f.read()

    species = content.split('<species name=')
    assert len(species) == 4          # the split leaves the preamble in front
    assert species[1].count('<wannierproj l=') == 3
    assert species[2].count('<wannierproj l=') == 3
    assert species[3].count('<wannierproj l=') == 2
    assert content.count('<wannierproj l=') == 8
    assert 'numBands="68"' in content          # twice the 34 Wannier functions

    clean_workdir()


@pytest.mark.inpgen
@pytest.mark.bulk
@pytest.mark.fast
def test_inpgen_wann_projections_refused(execute_inpgen, clean_workdir):
    """A key that cannot be honoured stops inpgen instead of writing a short answer.

    Without these, a typo leaves a species with no projections and the run only
    fails much later, when wannierlib_amn finds nothing to project onto. The
    fixture reports a non-zero inpgen with pytest.fail, whose exception derives
    from BaseException rather than Exception.
    """
    folder = './inputfiles/inpgen/Si_no_para/'

    for bad in ('Fw:sp',            # not an element symbol
                'Cu:sp',            # an element, but not in this structure
                'Si:sp,Si:d',       # the same element twice
                'Si:spx',           # not an s, p, d or f
                'Si:',              # no channels at all
                'Si:sp#nbands=4',   # fewer bands than Wannier functions
                'Si:sp#kmax=4'):    # not a modifier
        with pytest.raises(BaseException):
            execute_inpgen(folder, ['-f', 'inp_simple', '-wann', bad])
        clean_workdir()
