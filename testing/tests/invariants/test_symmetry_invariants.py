
import re
from xml.etree import ElementTree

import pytest

"""
Differential regression tests for the symmetry treatment.

Each test runs the same system once with its symmetry operations and once
without (identity only, one atom group per atom, full k-point mesh) and
asserts that both give the same result. They fail against the pre-fix code
rather than against stored reference numbers. See issues #838 to #842.
"""


def _iterations(outxml):
    return ElementTree.parse(outxml).getroot().findall("./scfLoop/iteration")


def _forces(iteration):
    return {int(f.attrib["atomType"]): tuple(float(f.attrib[c]) for c in ("F_x", "F_y", "F_z"))
            for f in iteration.findall("totalForcesOnRepresentativeAtoms/forceTotal")}


def _scf_data(outxml):
    data = []
    for it in _iterations(outxml):
        mt = {int(m.attrib["atomType"]): float(m.attrib["total"])
              for m in it.findall("valenceDensity/mtCharges/mtCharge")}
        vac = [float(v) for k, v in it.find("valenceDensity/spinDependentCharge").attrib.items()
               if k.startswith("vacuum")]
        data.append({"etot": float(it.find("totalEnergy").attrib["value"]),
                     "dist": float(it.find("densityConvergence/chargeDensity").attrib["distance"]),
                     "mt": mt, "vac": vac, "forces": _forces(it)})
    return data


def _xsf_values(path, name):
    """Data of the 2D grid block called name; one xsf file holds all plots of a run."""
    text = open(path).read()
    match = re.search(r"BEGIN_BLOCK_DATAGRID_2D\s*\n\s*" + name + r"\s*\n\s*BEGIN_DATAGRID_2D\S*\s*\n(.*?)END_DATAGRID",
                      text, re.S)
    assert match, f"plot block '{name}' missing in {path}"
    return [float(v) for v in match.group(1).split()[11:]]  # skip dims, origin and spanning vectors


@pytest.mark.fleur
@pytest.mark.film
@pytest.mark.forces
@pytest.mark.plot
@pytest.mark.serial
def test_film_nonsymmorphic_matches_nosym(execute_fleur):
    """A nonsymmorphic film must give the same SCF path and plots without symmetry.

    Layer group with an in-plane mirror {m_y|(0,1/3)}, a z-glide {m_z|(1/2,0)}
    and a 2_1 screw axis: the vacua are related by an operation that is neither
    inversion nor z-reflection, and the 2D star phases are complex. Checks the
    vacuum density (complex star phases in vacden), the vacuum mixing metric,
    the film Kerker preconditioner (VYukawaFilm) and the plot of the lower
    vacuum (outcdn), with two plots in one run.
    """
    input_dir = "./inputfiles/film/FilmNonsymmorphicInvariant"

    def scf(name):
        res = execute_fleur(input_dir, rm_files=['.'],
                            only_copy=[[f"inp{name}.xml", "inp.xml"]], mpi_procs=1)
        return _scf_data(res["out.xml"])

    def plot_lower(name):
        scf(name)
        # no rm_files: the plot run reads the density of the SCF run
        res = execute_fleur(input_dir, only_copy=[[f"inp_plot{name}.xml", "inp.xml"]], mpi_procs=1)
        return _xsf_values(res["denIn.xsf"], "lower")  # the second plot, below the film

    sym, nosym = scf(""), scf("_nosym")
    assert len(sym) == len(nosym) == 3
    atom_type = [1, 1, 2, 2]  # nosym atom -> sym atom type
    for s, n in zip(sym, nosym):
        assert abs(s["etot"] - n["etot"]) < 1.0e-7, f"total energy: {s['etot']} vs {n['etot']}"
        assert abs(s["dist"] - n["dist"]) < 1.0e-5, f"distance: {s['dist']} vs {n['dist']}"
        for atom, t in enumerate(atom_type, start=1):
            assert abs(s["mt"][t] - n["mt"][atom]) < 5.0e-6, f"MT charge of atom {atom}"
        assert len(s["vac"]) == 1 and len(n["vac"]) == 2
        assert all(abs(s["vac"][0] - v) < 1.0e-6 for v in n["vac"]), f"vacuum charge: {s['vac']} vs {n['vac']}"
        for t, atom in ((1, 1), (2, 3)):  # representative atoms
            assert max(abs(a - b) for a, b in zip(s["forces"][t], n["forces"][atom])) < 1.0e-5, \
                f"force on atom {atom}: {s['forces'][t]} vs {n['forces'][atom]}"

    lower = {name: plot_lower(name) for name in ("", "_nosym")}
    assert len(lower[""]) == len(lower["_nosym"]) > 0
    scale = max(abs(v) for v in lower["_nosym"])
    assert max(abs(a - b) for a, b in zip(lower[""], lower["_nosym"])) < 1.0e-4 * scale, \
        "density plotted in the lower vacuum differs"


@pytest.mark.fleur
@pytest.mark.bulk
@pytest.mark.forces
@pytest.mark.serial
def test_coretail_force_level1_matches_nosym(execute_fleur):
    """Core-tail forces of force level 1 must not depend on the symmetry.

    Trigonal Se (P3_121) off equilibrium: three equivalent atoms related by a
    3_1 screw axis, so the star phases are complex. The muffin-tin part of the
    core-tail force has to be summed over all equivalent atoms; before the fix
    only the representative atom was used, weighted with the star phase.
    """
    input_dir = "./inputfiles/forces/SeCoreTailForceInvariant"
    forces = {}
    for name in ("", "_nosym"):
        res = execute_fleur(input_dir, rm_files=['.'],
                            only_copy=[[f"inp{name}.xml", "inp.xml"]], mpi_procs=1)
        forces[name] = _forces(_iterations(res["out.xml"])[-1])

    sym, nosym = forces[""][1], forces["_nosym"][1]
    assert abs(sym[0]) > 0.1, f"test needs a sizeable force: {sym}"
    assert max(abs(a - b) for a, b in zip(sym, nosym)) < 1.0e-5, f"{sym} vs {nosym}"
    # the 2-fold axis through the atom forces F_y = F_z = 0
    assert abs(sym[1]) < 1.0e-6 and abs(sym[2]) < 1.0e-6, f"symmetry-forbidden force components: {sym}"
