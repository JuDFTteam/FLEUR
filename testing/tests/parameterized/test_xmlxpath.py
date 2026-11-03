import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree

import pytest

"""
Tests of the -xmlXPath command line option (issue #822).
"""

INPUT = "./inputfiles/basic/CuBulk"


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.fast
def test_xmlxpath_sets_attribute(execute_fleur):
    """An attribute present in inp.xml is changed for this run"""
    res_files = execute_fleur(INPUT, rm_files=['.'], mpi_procs=1,
                              cmdline_param=["-xmlXPath", "/fleurInput/calculationSetup/scfLoop/@itmax=2"])
    root = ElementTree.parse(res_files["out.xml"]).getroot()
    assert len(root.findall(".//scfLoop/iteration")) == 2


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.xml
@pytest.mark.fast
def test_xmlxpath_missing_attribute(fleur_binary, tmp_path):
    """An XPath without a match in inp.xml stops FLEUR with an error instead of a segfault"""
    fleur, _ = fleur_binary
    shutil.copy(Path(INPUT) / "inp.xml", tmp_path)
    run = subprocess.Popen([fleur, "-xmlXPath", "/fleurInput/output/bandDOS/@orbcomp=T"],
                           cwd=tmp_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        stdout, stderr = run.communicate(timeout=120)
    except subprocess.TimeoutExpired:
        # no wait(): a crashed process may hang while exiting
        run.kill()
        pytest.fail("FLEUR did not stop within 120 s")
    output = stdout + stderr
    assert "Segmentation fault" not in output
    assert "Attribute value could not be set" in output
    assert "Run finished successfully" not in output
