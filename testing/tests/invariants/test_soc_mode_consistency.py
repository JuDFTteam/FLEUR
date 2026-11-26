import math
import re
import shutil
from pathlib import Path
from xml.etree import ElementTree

import pytest

"""
Consistency of the three SOC implementations (issue #846).

FLEUR evaluates spin-orbit coupling in three independent ways:

  FV   first variation, noco + SOC (hsmt_soc_offdiag), moments along n
  SV   second variation, collinear SOC (secvar_soc / t_matelements_soc), SOC axis n
  DMI  first-order force theorem (ssomat), q = 0, moments along z, SOC axis n

For the same unperturbed state all three must give the same first-order SOC
energy at a single k-point,

  dE/dlambda = sum_{b occ} <b|H_so|b> = n . V,

where lambda scales the SOC strength (special/@socscale). FV and SV are
differentiated numerically (central difference in lambda), the force theorem
yields the first-order term directly. The film has no inversion symmetry and
the k-point is general, so V is a generic vector and an error in the
orientation of the SOC axis (e.g. phi -> -phi) changes the result.

All runs start from the same density: the noco starting density for moments
along z. FV and DMI generate it themselves (rotated to n for FV); the collinear
SV runs read it, because the collinear starting density differs slightly.
"""

BASE = Path(__file__).parents[2] / "inputfiles" / "noco" / "FePt_film_SOCModes" / "inp.xml"
LAMBDA = 0.005
#: SOC axes (theta, phi): along x, along y and a generic direction
AXES = {"x": (math.pi / 2, 0.0), "y": (math.pi / 2, math.pi / 2), "g": (1.0, 2.0)}
#: absolute tolerance on dE/dlambda in Htr; the values are of order 1e-3
TOL = 2.0e-6


def _make_input(mode, theta=0.0, phi=0.0, lam=0.0, ef_shift=0.0):
    s = BASE.read_text()
    soc = lam != 0.0
    if mode == "FV":
        mag = '<magnetism jspins="2" l_noco="T" l_ss="F" swsp="F" lflip="F">'
        socl = f'<soc theta="0.0" phi="0.0" l_soc="{"T" if soc else "F"}" spav="F"/>'
        alpha, beta = phi, theta
    elif mode == "SV":
        mag = '<magnetism jspins="2" l_noco="F" l_ss="F" swsp="F" lflip="F">'
        socl = f'<soc theta="{theta!r}" phi="{phi!r}" l_soc="{"T" if soc else "F"}" spav="F"/>'
        alpha, beta = 0.0, 0.0
    else:  # DMI force theorem with q = 0
        mag = '<magnetism jspins="2" l_noco="T" l_ss="T" swsp="F" lflip="F">'
        socl = '<soc theta="0.0" phi="0.0" l_soc="F" spav="F"/>'
        alpha, beta = 0.0, 0.0
    s = re.sub(r'<magnetism [^>]*>', mag, s)
    s = re.sub(r'<soc [^>]*/>', socl, s)
    s = re.sub(r'<nocoParams [^>]*/>', f'<nocoParams alpha="{alpha!r}" beta="{beta!r}"/>', s)
    s = re.sub(r'<special socscale="[^"]*"/>', f'<special socscale="{lam if soc else 1.0!r}"/>', s)
    if mode == "DMI":
        thetas = " ".join(repr(t) for t, _ in AXES.values())
        phis = " ".join(repr(p) for _, p in AXES.values())
        s = s.replace("</output>", "</output>\n   <forceTheorem>\n"
                      f'      <DMI theta="{thetas}" phi="{phis}" ef_shift="{ef_shift!r}">\n'
                      "         <qVectors>\n            <q> 0.0 0.0 0.0 </q>\n         </qVectors>\n"
                      "      </DMI>\n   </forceTheorem>")
    return s


def _occupied_sum(outxml, nocc):
    """Sum of the lowest nocc eigenvalues; SV stores its eigenvalues for both spins."""
    root = ElementTree.parse(outxml).getroot()
    blocks = [[float(v) for v in e.text.split()] for e in root.findall(".//eigenvaluesAt")]
    if len(blocks) == 2 and blocks[0] == blocks[1]:
        blocks = blocks[:1]
    eig = sorted(v for b in blocks for v in b)
    return sum(eig[:nocc]), eig


def _valence_electrons(outxml):
    root = ElementTree.parse(outxml).getroot()
    return int(round(float(root.find(".//bzIntegration").attrib["valenceElectrons"])))


@pytest.mark.fleur
@pytest.mark.noco
@pytest.mark.soc
@pytest.mark.film
@pytest.mark.spinspiral
@pytest.mark.forcetheorem
@pytest.mark.hdf
@pytest.mark.slow
def test_soc_modes_first_order_consistent(execute_fleur, work_dir):
    """First-order SOC energy agrees between first variation, second variation and DMI.

    Before the fix of #846 second variation used the axis (theta, -phi): its
    value for the y axis has the opposite sign of the other two modes.
    """
    work = Path(work_dir)
    keep = work / "socmodes"  # inputs and out.xml of all runs, kept for inspection
    keep.mkdir(exist_ok=True)

    def run(label, inp, cdn=None):
        for f in work.iterdir():
            if f.is_file():
                f.unlink()
            elif f != keep:
                shutil.rmtree(f)
        (work / "inp.xml").write_text(inp)
        if cdn is not None:
            shutil.copy(cdn, work / "cdn.hdf")
        res = execute_fleur(mpi_procs=1)
        shutil.copy(work / "inp.xml", keep / f"{label}_inp.xml")
        shutil.copy(res["out.xml"], keep / f"{label}_out.xml")
        return keep / f"{label}_out.xml"

    # unperturbed noco state, moments along z; its starting density is kept for SV
    start_xml = run("start", _make_input("FV"))
    start_cdn = keep / "start_cdn.hdf"
    shutil.copy(work / "cdn_old.hdf", start_cdn)
    nocc = _valence_electrons(start_xml)
    e0, eig0 = _occupied_sum(start_xml, nocc)
    gap = eig0[nocc] - eig0[nocc - 1]
    assert gap > 2.0e-3, f"no gap at the Fermi level of the test k-point: {gap}"

    # the same unperturbed state in all modes
    sv0, _ = _occupied_sum(run("SV_nosoc", _make_input("SV"), cdn=start_cdn), nocc)
    assert abs(sv0 - e0) < 1.0e-5, f"SV starts from a different state: {sv0} vs {e0}"

    first_order = {}
    for mode in ("FV", "SV"):
        for name, (theta, phi) in AXES.items():
            cdn = start_cdn if mode == "SV" else None
            ep, _ = _occupied_sum(run(f"{mode}_{name}_p", _make_input(mode, theta, phi, LAMBDA), cdn), nocc)
            em, _ = _occupied_sum(run(f"{mode}_{name}_m", _make_input(mode, theta, phi, -LAMBDA), cdn), nocc)
            first_order[(mode, name)] = (ep - em) / (2 * LAMBDA)

    # force theorem: Fermi level moved to the middle of the gap, so that exactly
    # nocc states are occupied (hist mode puts it on the highest occupied state)
    dmi_xml = run("DMI", _make_input("DMI", ef_shift=round(gap / 2, 6)))
    root = ElementTree.parse(dmi_xml).getroot()
    hso = [float(e.attrib["H_so"]) for e in root.findall(".//Forcetheorem_DMI/allAtoms")]
    assert len(hso) == len(AXES)
    for name, value in zip(AXES, hso):
        first_order[("DMI", name)] = value

    table = "\n".join(f"{name}: FV {first_order[('FV', name)]: .8f}  SV {first_order[('SV', name)]: .8f}"
                      f"  DMI {first_order[('DMI', name)]: .8f}" for name in AXES)
    print("first-order SOC energy dE/dlambda (Htr)\n" + table)

    for name in AXES:
        fv, sv, dmi = (first_order[(m, name)] for m in ("FV", "SV", "DMI"))
        assert abs(fv) > 1.0e-3, f"axis {name}: first-order term too small to test\n{table}"
        assert abs(dmi - fv) < TOL, f"axis {name}: DMI force theorem differs from first variation\n{table}"
        assert abs(sv - fv) < TOL, f"axis {name}: second variation differs from first variation\n{table}"
