"""
FLEURiste - A Textual-based GUI editor for FLEUR inp.xml files.

This package provides a terminal-based graphical editor for FLEUR input files,
using the FleurInputSchema.xsd to provide structure-aware editing capabilities.
Also includes a K-Point manager for managing k-point lists in inp.xml files,
an Inpgen mode for generating inp.xml from input files or structure files,
and a Job Generator mode for creating SLURM job scripts.
"""

__version__ = "0.1.0"
__author__ = "FLEUR Team"

from .schema_parser import XSDParser, SchemaElement
from .xml_editor import XMLDocument, XMLNode
from .app import FLEURisteApp
from .cli import cli, main as cli_main
from .kpoint_manager import InpXMLManager, KPoint, KPointList, KPointMode, KPointModifiers
from .kpoint_gui import KPointPlot, CreateMeshDialog, CreateDensityDialog, CreatePathDialog
from .inpgen_gui import (
    FilePickerDialog, SupercellDialog, MagneticMomentDialog, SurfaceDialog,
    INPGEN_PROFILES, STRUCTURE_FORMATS, SURFACE_PRESETS, ELEMENT_LATTICE_CONSTANTS,
    ase_to_fleur_input, create_supercell, parse_namelist_to_atoms,
    get_elements_in_input, add_magnetic_moments, extract_magnetic_moments,
    generate_surface
)

# Lazy-import Jupyter interface (requires ipywidgets)
def jupyter(*args, **kwargs):
    """Create and return a FLEURisteJupyter GUI instance.

    Shorthand for::

        from fleuriste.jupyter_gui import FLEURisteJupyter
        gui = FLEURisteJupyter(...)
        gui.display()
    """
    from .jupyter_gui import FLEURisteJupyter
    gui = FLEURisteJupyter(*args, **kwargs)
    gui.display()
    return gui

__all__ = [
    "cli",
    "cli_main",
    "XSDParser",
    "SchemaElement", 
    "XMLDocument",
    "XMLNode",
    "FLEURisteApp",
    "InpXMLManager",
    "KPoint",
    "KPointList",
    "KPointMode",
    "KPointModifiers",
    "KPointPlot",
    "CreateMeshDialog",
    "CreateDensityDialog",
    "CreatePathDialog",
    "FilePickerDialog",
    "SupercellDialog",
    "MagneticMomentDialog",
    "INPGEN_PROFILES",
    "STRUCTURE_FORMATS",
    "ase_to_fleur_input",
    "create_supercell",
    "parse_namelist_to_atoms",
    "get_elements_in_input",
    "add_magnetic_moments",
    "jupyter",
]
