"""
Inpgen GUI components for FLEURiste.
Provides widgets for generating FLEUR inp.xml files from input files or structure files.
"""

from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical, VerticalScroll
from textual.widgets import (
    Button, Label, Input, Static, Select, TextArea, Rule, DirectoryTree
)
from textual.binding import Binding
from textual.screen import ModalScreen
from textual import on
from pathlib import Path
from typing import Optional, Tuple, Dict, List
import os
import numpy as np

# Optional imports
try:
    from FleurInpgen import InpgenInterface
    HAS_INPGEN = True
except ImportError:
    HAS_INPGEN = False
    InpgenInterface = None

try:
    from ase.io import read as ase_read
    from ase.build import make_supercell
    ASE_AVAILABLE = True
except ImportError:
    ASE_AVAILABLE = False


def ase_to_fleur_input(atoms, **kwargs) -> str:
    """
    Convert ASE Atoms object to FLEUR namelist input format.
    
    Args:
        atoms: ASE Atoms object
        **kwargs: Additional parameters:
            - film: bool, whether this is a film calculation
            - magnetic_moments: dict mapping atomic number to moment specification
    
    Returns:
        Namelist input string for FLEUR
    
    Note:
        For film calculations, the z-coordinate is in absolute units (Bohr)
        and atoms are centered around z=0.
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE library is required for this functionality")
    
    lines = []
    
    # Input namelist
    film = kwargs.get('film', False)
    film_str = "t" if film else "f"
    lines.append(f"&input film={film_str} /")
    
    # Bravais matrix - convert from Angstrom to Bohr
    angstrom_to_bohr = 1.88972612456506
    cell = atoms.get_cell() * angstrom_to_bohr
    
    # Output the three lattice vectors (in Bohr)
    for i in range(3):
        vec = cell[i]
        lines.append(f"{vec[0]:.8f} {vec[1]:.8f} {vec[2]:.8f}")
    
    # Scale factor
    lines.append("1.0")
    
    # Number of atoms
    lines.append(str(len(atoms)))
    
    # Get magnetic moments if provided
    magnetic_moments = kwargs.get('magnetic_moments', {})
    
    if film:
        # For films: x,y in scaled coordinates, z in absolute Bohr centered around 0
        scaled_positions = atoms.get_scaled_positions()
        positions_cart = atoms.get_positions() * angstrom_to_bohr  # in Bohr
        
        # Center z around 0
        z_coords = positions_cart[:, 2]
        z_center = (z_coords.max() + z_coords.min()) / 2.0
        z_centered = z_coords - z_center
        
        # Atom positions
        for i, atom in enumerate(atoms):
            elem_z = atom.number
            # x, y in scaled; z in absolute Bohr
            line = f"{elem_z}  {scaled_positions[i][0]:.8f} {scaled_positions[i][1]:.8f} {z_centered[i]:.8f}"
            
            # Add magnetic moment if specified for this element
            if elem_z in magnetic_moments:
                line += f" : {magnetic_moments[elem_z]}"
            
            lines.append(line)
    else:
        # For bulk: all coordinates in scaled units
        scaled_positions = atoms.get_scaled_positions()
        
        # Atom positions
        for i, atom in enumerate(atoms):
            elem_z = atom.number
            scaled_pos = scaled_positions[i]
            line = f"{elem_z}  {scaled_pos[0]:.8f} {scaled_pos[1]:.8f} {scaled_pos[2]:.8f}"
            
            # Add magnetic moment if specified for this element
            if elem_z in magnetic_moments:
                line += f" : {magnetic_moments[elem_z]}"
            
            lines.append(line)
     
    return "\n".join(lines) + "\n"


def create_supercell(atoms, nx: int, ny: int, nz: int):
    """
    Create a supercell from an ASE Atoms object.
    
    Uses ASE's make_supercell function with a diagonal transformation matrix.
    
    Args:
        atoms: ASE Atoms object
        nx: Number of repetitions along x (a1)
        ny: Number of repetitions along y (a2)
        nz: Number of repetitions along z (a3)
    
    Returns:
        ASE Atoms object of the supercell
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE library is required for supercell creation")
    
    # Create diagonal transformation matrix
    P = np.array([[nx, 0, 0],
                  [0, ny, 0],
                  [0, 0, nz]])
    
    # Use ASE's make_supercell
    supercell = make_supercell(atoms, P)
    
    return supercell


def parse_namelist_to_atoms(input_text: str):
    """
    Parse FLEUR namelist input to an ASE Atoms object.
    
    Args:
        input_text: FLEUR namelist input string
    
    Returns:
        ASE Atoms object or None if parsing fails
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE library is required")
    
    from ase import Atoms
    
    lines = [l.strip() for l in input_text.strip().split('\n') if l.strip()]
    
    # Find the lattice vectors (skip &input line)
    idx = 0
    # Skip namelist lines
    while idx < len(lines) and (lines[idx].startswith('&') or lines[idx].startswith('/')):
        idx += 1
    
    if idx + 5 > len(lines):
        return None
    
    try:
        # Read lattice vectors (in Bohr)
        bohr_to_angstrom = 0.529177210903
        cell = []
        for i in range(3):
            vec = [float(x) * bohr_to_angstrom for x in lines[idx + i].split()]
            cell.append(vec)
        
        # Scale factor
        scale = float(lines[idx + 3])
        cell = np.array(cell) * scale
        
        # Number of atoms
        natoms = int(lines[idx + 4])
        
        # Read atom positions
        symbols = []
        positions = []
        atomic_numbers = {
            1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
            9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
            16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
            23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
            30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
            37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
            44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
            51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
            72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
            79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi',
        }
        
        for i in range(natoms):
            parts = lines[idx + 5 + i].split()
            z = int(parts[0])
            symbol = atomic_numbers.get(z, 'X')
            symbols.append(symbol)
            # Scaled positions
            scaled_pos = [float(parts[1]), float(parts[2]), float(parts[3])]
            positions.append(scaled_pos)
        
        # Create Atoms object with scaled positions
        atoms = Atoms(symbols=symbols, scaled_positions=positions, cell=cell, pbc=True)
        return atoms
        
    except (ValueError, IndexError) as e:
        return None


def extract_magnetic_moments(input_text: str) -> Dict[int, str]:
    """
    Extract magnetic moment specifications from FLEUR namelist input.
    
    Args:
        input_text: FLEUR namelist input string
    
    Returns:
        Dict mapping atomic number to moment specification string
    """
    magnetic_moments = {}
    lines = [l.strip() for l in input_text.strip().split('\n') if l.strip()]
    
    # Find atom lines
    in_atoms = False
    natoms = 0
    atom_count = 0
    
    for line in lines:
        if line.startswith('&') or line.startswith('/'):
            continue
        
        parts = line.split()
        
        # Check for number of atoms line
        if not in_atoms and len(parts) == 1 and parts[0].isdigit():
            natoms = int(parts[0])
            in_atoms = True
            continue
        
        if in_atoms and atom_count < natoms:
            # Check for magnetic moment (after colon)
            if ':' in line:
                try:
                    z = int(parts[0])
                    moment_part = line.split(':')[1].strip()
                    if moment_part:
                        magnetic_moments[z] = moment_part
                except (ValueError, IndexError):
                    pass
            atom_count += 1
    
    return magnetic_moments


class SupercellDialog(ModalScreen):
    """Modal dialog for supercell parameters."""
    
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    
    CSS = """
    #supercell-dialog {
        width: 50;
        height: auto;
        background: $panel;
        border: thick $primary;
        padding: 1 2;
    }
    
    #supercell-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 1;
        text-align: center;
    }
    
    .supercell-row {
        height: auto;
        margin-bottom: 1;
    }
    
    .supercell-row Label {
        width: 15;
    }
    
    .supercell-row Input {
        width: 10;
    }
    
    #supercell-info {
        margin: 1 0;
        color: $text-muted;
    }
    
    #supercell-buttons {
        margin-top: 1;
        height: auto;
    }
    """
    
    def __init__(self, current_atoms: int = 0):
        super().__init__()
        self.current_atoms = current_atoms
    
    def compose(self) -> ComposeResult:
        with Container(id="supercell-dialog"):
            yield Label("Create Supercell", id="supercell-title")
            yield Rule()
            
            with Horizontal(classes="supercell-row"):
                yield Label("Repeat along a₁:")
                yield Input(value="2", id="input-nx", type="integer")
            
            with Horizontal(classes="supercell-row"):
                yield Label("Repeat along a₂:")
                yield Input(value="2", id="input-ny", type="integer")
            
            with Horizontal(classes="supercell-row"):
                yield Label("Repeat along a₃:")
                yield Input(value="2", id="input-nz", type="integer")
            
            info = f"Current atoms: {self.current_atoms}" if self.current_atoms > 0 else "Load a structure first"
            yield Static(info, id="supercell-info")
            
            with Horizontal(id="supercell-buttons"):
                yield Button("Create", variant="primary", id="btn-create-supercell")
                yield Button("Cancel", variant="default", id="btn-cancel-supercell")
    
    @on(Input.Changed)
    def on_input_changed(self, event: Input.Changed):
        """Update info when input changes."""
        if self.current_atoms > 0:
            try:
                nx = int(self.query_one("#input-nx", Input).value or "1")
                ny = int(self.query_one("#input-ny", Input).value or "1")
                nz = int(self.query_one("#input-nz", Input).value or "1")
                total = self.current_atoms * nx * ny * nz
                self.query_one("#supercell-info", Static).update(
                    f"Current: {self.current_atoms} atoms → Supercell: {total} atoms"
                )
            except ValueError:
                pass
    
    @on(Button.Pressed, "#btn-create-supercell")
    def on_create(self):
        try:
            nx = int(self.query_one("#input-nx", Input).value or "1")
            ny = int(self.query_one("#input-ny", Input).value or "1")
            nz = int(self.query_one("#input-nz", Input).value or "1")
            
            if nx < 1 or ny < 1 or nz < 1:
                self.notify("All values must be >= 1", severity="error")
                return
            
            self.dismiss((nx, ny, nz))
        except ValueError:
            self.notify("Please enter valid integers", severity="error")
    
    @on(Button.Pressed, "#btn-cancel-supercell")
    def on_cancel(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)


class FilePickerDialog(ModalScreen):
    """Modal dialog for picking a file."""
    
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    
    CSS = """
    #file-picker-dialog {
        width: 80%;
        height: 80%;
        background: $panel;
        border: thick $primary;
        padding: 1;
    }
    
    #file-picker-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 1;
        text-align: center;
    }
    
    #file-tree {
        height: 1fr;
        border: solid $primary;
    }
    
    #file-picker-buttons {
        margin-top: 1;
        height: auto;
    }
    
    #selected-file-label {
        margin-top: 1;
        color: $text-muted;
    }
    """
    
    def __init__(self, start_path: str = ".", title: str = "Select File"):
        super().__init__()
        self.start_path = start_path
        self.dialog_title = title
        self.selected_path: Optional[Path] = None
    
    def compose(self) -> ComposeResult:
        with Container(id="file-picker-dialog"):
            yield Label(self.dialog_title, id="file-picker-title")
            yield DirectoryTree(self.start_path, id="file-tree")
            yield Label("No file selected", id="selected-file-label")
            with Horizontal(id="file-picker-buttons"):
                yield Button("Select", variant="primary", id="btn-select-file")
                yield Button("Cancel", variant="default", id="btn-cancel-file")
    
    @on(DirectoryTree.FileSelected)
    def on_file_selected(self, event: DirectoryTree.FileSelected):
        """Handle file selection in tree."""
        self.selected_path = event.path
        label = self.query_one("#selected-file-label", Label)
        label.update(f"Selected: {event.path.name}")
    
    @on(Button.Pressed, "#btn-select-file")
    def on_select(self):
        if self.selected_path:
            self.dismiss(self.selected_path)
        else:
            self.notify("Please select a file first", severity="warning")
    
    @on(Button.Pressed, "#btn-cancel-file")
    def on_cancel(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)


# List of profiles that can be used with inpgen
INPGEN_PROFILES = [
    ("default_profile", "Default Profile"),
    ("careful_profile", "Careful Profile"),
    ("fast_profile", "Fast Profile"),
]

# Common structure file formats for ASE
STRUCTURE_FORMATS = [
    ("auto", "Auto-detect"),
    ("vasp", "VASP POSCAR"),
    ("cif", "CIF"),
    ("xyz", "XYZ"),
    ("pdb", "PDB"),
    ("extxyz", "Extended XYZ"),
]


# Atomic number to symbol mapping
ATOMIC_SYMBOLS = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
    9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
    16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
    23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu',
    30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr',
    37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc',
    44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn',
    51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
    58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd',
    65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu',
    72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt',
    79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At',
    86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
}

# Symbol to atomic number mapping
SYMBOL_TO_Z = {v: k for k, v in ATOMIC_SYMBOLS.items()}


def get_elements_in_input(input_text: str) -> List[Tuple[int, str]]:
    """
    Get list of unique elements (Z, symbol) in the input.
    
    Returns:
        List of (atomic_number, symbol) tuples
    """
    elements = set()
    lines = [l.strip() for l in input_text.strip().split('\n') if l.strip()]
    
    # Find atom lines (after header)
    in_atoms = False
    natoms = 0
    atom_count = 0
    
    for line in lines:
        if line.startswith('&') or line.startswith('/'):
            continue
        
        parts = line.split()
        if len(parts) == 1 and parts[0].isdigit():
            # This is the number of atoms line
            natoms = int(parts[0])
            in_atoms = True
            continue
        
        if in_atoms and atom_count < natoms:
            try:
                z = int(parts[0])
                symbol = ATOMIC_SYMBOLS.get(z, f'Z{z}')
                elements.add((z, symbol))
                atom_count += 1
            except (ValueError, IndexError):
                pass
    
    return sorted(list(elements))


def add_magnetic_moments(input_text: str, element_z: int, moment_spec: str) -> str:
    """
    Add magnetic moment specification to atoms of a specific element.
    
    Args:
        input_text: FLEUR namelist input string
        element_z: Atomic number of the element to magnetize
        moment_spec: Moment specification (e.g., "2.4", "up", "down", "0.0 2.4 0.0")
    
    Returns:
        Modified input string with magnetic moments
    """
    lines = input_text.strip().split('\n')
    result = []
    
    in_atoms = False
    natoms = 0
    atom_count = 0
    
    for line in lines:
        stripped = line.strip()
        
        if stripped.startswith('&') or stripped.startswith('/'):
            result.append(line)
            continue
        
        parts = stripped.split()
        
        # Check for number of atoms line
        if not in_atoms and len(parts) == 1 and parts[0].isdigit():
            natoms = int(parts[0])
            in_atoms = True
            result.append(line)
            continue
        
        if in_atoms and atom_count < natoms:
            try:
                z = int(parts[0])
                # Check if this is the target element
                if z == element_z:
                    # Remove any existing moment specification
                    if ':' in stripped:
                        base_line = stripped.split(':')[0].strip()
                    else:
                        # Get position part only (first 4 values: Z x y z)
                        base_parts = parts[:4]
                        base_line = ' '.join(base_parts)
                    # Add new moment
                    new_line = f"{base_line} : {moment_spec}"
                    result.append(new_line)
                else:
                    result.append(line)
                atom_count += 1
            except (ValueError, IndexError):
                result.append(line)
        else:
            result.append(line)
    
    return '\n'.join(result) + '\n'


class MagneticMomentDialog(ModalScreen):
    """Modal dialog for setting magnetic moments."""
    
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    
    CSS = """
    #magnetic-dialog {
        width: 60;
        height: auto;
        background: $panel;
        border: thick $primary;
        padding: 1 2;
    }
    
    #magnetic-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 1;
        text-align: center;
    }
    
    .magnetic-row {
        height: auto;
        margin-bottom: 1;
    }
    
    .magnetic-row Label {
        width: 18;
    }
    
    .magnetic-row Input {
        width: 1fr;
    }
    
    .magnetic-row Select {
        width: 1fr;
    }
    
    #magnetic-info {
        margin: 1 0;
        color: $text-muted;
    }
    
    #noncollinear-inputs {
        display: none;
    }
    
    #noncollinear-inputs.visible {
        display: block;
    }
    
    #collinear-inputs {
        display: block;
    }
    
    #collinear-inputs.hidden {
        display: none;
    }
    
    #magnetic-buttons {
        margin-top: 1;
        height: auto;
    }
    """
    
    def __init__(self, elements: List[Tuple[int, str]]):
        """
        Initialize with list of elements in the structure.
        
        Args:
            elements: List of (atomic_number, symbol) tuples
        """
        super().__init__()
        self.elements = elements
        self._is_noncollinear = False
    
    def compose(self) -> ComposeResult:
        with Container(id="magnetic-dialog"):
            yield Label("Set Magnetic Moments", id="magnetic-title")
            yield Rule()
            
            with Horizontal(classes="magnetic-row"):
                yield Label("Element:")
                if self.elements:
                    options = [(f"{sym} (Z={z})", z) for z, sym in self.elements]
                    yield Select(options, value=self.elements[0][0], id="select-element")
                else:
                    yield Select([("No elements found", 0)], value=0, id="select-element")
            
            with Horizontal(classes="magnetic-row"):
                yield Label("Type:")
                yield Select(
                    [("Collinear", "collinear"), ("Non-collinear", "noncollinear")],
                    value="collinear",
                    id="select-magtype"
                )
            
            # Collinear inputs
            with Vertical(id="collinear-inputs"):
                with Horizontal(classes="magnetic-row"):
                    yield Label("Moment:")
                    yield Select(
                        [("Custom value", "custom"), ("up (default)", "up"), ("down (default)", "down")],
                        value="custom",
                        id="select-moment-type"
                    )
                with Horizontal(classes="magnetic-row"):
                    yield Label("Value (μB):")
                    yield Input(value="2.0", id="input-moment", type="number")
            
            # Non-collinear inputs
            with Vertical(id="noncollinear-inputs"):
                with Horizontal(classes="magnetic-row"):
                    yield Label("mx (or keyword):")
                    yield Input(value="0.0", id="input-mx")
                with Horizontal(classes="magnetic-row"):
                    yield Label("my:")
                    yield Input(value="0.0", id="input-my")
                with Horizontal(classes="magnetic-row"):
                    yield Label("mz:")
                    yield Input(value="2.0", id="input-mz")
                yield Static("Use 'up'/'down' keywords or numeric values", id="magnetic-info")
            
            with Horizontal(id="magnetic-buttons"):
                yield Button("Apply", variant="primary", id="btn-apply-magnetic")
                yield Button("Cancel", variant="default", id="btn-cancel-magnetic")
    
    @on(Select.Changed, "#select-magtype")
    def on_magtype_changed(self, event: Select.Changed):
        """Toggle between collinear and non-collinear inputs."""
        is_nc = event.value == "noncollinear"
        self._is_noncollinear = is_nc
        
        collinear = self.query_one("#collinear-inputs", Vertical)
        noncollinear = self.query_one("#noncollinear-inputs", Vertical)
        
        if is_nc:
            collinear.add_class("hidden")
            noncollinear.add_class("visible")
        else:
            collinear.remove_class("hidden")
            noncollinear.remove_class("visible")
    
    @on(Button.Pressed, "#btn-apply-magnetic")
    def on_apply(self):
        element_z = self.query_one("#select-element", Select).value
        
        if element_z == 0:
            self.notify("No valid element selected", severity="error")
            return
        
        if self._is_noncollinear:
            # Non-collinear: mx my mz
            mx = self.query_one("#input-mx", Input).value.strip()
            my = self.query_one("#input-my", Input).value.strip()
            mz = self.query_one("#input-mz", Input).value.strip()
            moment_spec = f"{mx} {my} {mz}"
        else:
            # Collinear
            moment_type = self.query_one("#select-moment-type", Select).value
            if moment_type == "up":
                moment_spec = "up"
            elif moment_type == "down":
                moment_spec = "down"
            else:
                moment_spec = self.query_one("#input-moment", Input).value.strip()
        
        self.dismiss((element_z, moment_spec))
    
    @on(Button.Pressed, "#btn-cancel-magnetic")
    def on_cancel(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)


# Surface generation presets
SURFACE_PRESETS = [
    ("fcc100", "FCC (100)"),
    ("fcc110", "FCC (110)"),
    ("fcc111", "FCC (111)"),
    ("bcc100", "BCC (100)"),
    ("bcc110", "BCC (110)"),
    ("bcc111", "BCC (111)"),
    ("hcp0001", "HCP (0001)"),
    ("hcp10m10", "HCP (10-10)"),
    ("diamond100", "Diamond (100)"),
    ("diamond111", "Diamond (111)"),
    ("custom", "Custom Miller indices"),
]

# Common elements with default lattice constants (in Angstrom)
ELEMENT_LATTICE_CONSTANTS = {
    "Fe": (2.87, "bcc"),
    "Ni": (3.52, "fcc"),
    "Cu": (3.61, "fcc"),
    "Ag": (4.09, "fcc"),
    "Au": (4.08, "fcc"),
    "Pt": (3.92, "fcc"),
    "Pd": (3.89, "fcc"),
    "Al": (4.05, "fcc"),
    "Co": (2.51, "hcp"),  # a parameter
    "Ti": (2.95, "hcp"),  # a parameter
    "W": (3.16, "bcc"),
    "Mo": (3.15, "bcc"),
    "Cr": (2.88, "bcc"),
    "V": (3.02, "bcc"),
    "Si": (5.43, "diamond"),
    "Ge": (5.66, "diamond"),
    "C": (3.57, "diamond"),
}


def generate_surface(element: str, surface_type: str, layers: int, vacuum: float,
                     lattice_const: float = None, size: Tuple[int, int, int] = (1, 1, 1),
                     orthogonal: bool = True, miller: Tuple[int, int, int] = None):
    """
    Generate a surface slab using ASE's surface building functions.
    
    Args:
        element: Chemical symbol (e.g., 'Fe', 'Cu')
        surface_type: Surface preset (e.g., 'fcc111', 'bcc100') or 'custom'
        layers: Number of atomic layers
        vacuum: Vacuum thickness in Angstrom
        lattice_const: Lattice constant in Angstrom (optional, uses default if not given)
        size: Supercell size (nx, ny, 1) - z is typically 1 for slabs
        orthogonal: Whether to use orthogonal cell
        miller: Miller indices for custom surface (h, k, l)
    
    Returns:
        ASE Atoms object with the surface slab
    """
    if not ASE_AVAILABLE:
        raise ImportError("ASE library is required for surface generation")
    
    from ase.build import (
        fcc100, fcc110, fcc111,
        bcc100, bcc110, bcc111,
        hcp0001, hcp10m10,
        diamond100, diamond111,
        surface
    )
    from ase.data import atomic_numbers, covalent_radii
    
    # Get lattice constant
    if lattice_const is None:
        if element in ELEMENT_LATTICE_CONSTANTS:
            lattice_const = ELEMENT_LATTICE_CONSTANTS[element][0]
        else:
            # Estimate from covalent radius
            z = atomic_numbers.get(element, 1)
            lattice_const = 2.5 * covalent_radii[z]
    
    # Map surface type to builder function
    surface_builders = {
        'fcc100': lambda: fcc100(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'fcc110': lambda: fcc110(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'fcc111': lambda: fcc111(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'bcc100': lambda: bcc100(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'bcc110': lambda: bcc110(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'bcc111': lambda: bcc111(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'hcp0001': lambda: hcp0001(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'hcp10m10': lambda: hcp10m10(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'diamond100': lambda: diamond100(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        'diamond111': lambda: diamond111(element, size=size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
    }
    
    if surface_type == 'custom' and miller is not None:
        # Use general surface function
        from ase.build import bulk
        # Determine crystal structure from element
        if element in ELEMENT_LATTICE_CONSTANTS:
            struct = ELEMENT_LATTICE_CONSTANTS[element][1]
        else:
            struct = 'fcc'  # default
        
        # Build bulk first
        if struct == 'hcp':
            bulk_atoms = bulk(element, struct, a=lattice_const, c=lattice_const * 1.633)
        else:
            bulk_atoms = bulk(element, struct, a=lattice_const)
        
        # Create surface
        slab = surface(bulk_atoms, miller, layers, vacuum=vacuum)
        
        # Apply size
        if size[0] > 1 or size[1] > 1:
            slab = slab.repeat((size[0], size[1], 1))
        
        return slab
    elif surface_type in surface_builders:
        # Update size to include layers
        build_size = (size[0], size[1], layers)
        surface_builders_with_layers = {
            'fcc100': lambda: fcc100(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'fcc110': lambda: fcc110(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'fcc111': lambda: fcc111(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'bcc100': lambda: bcc100(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'bcc110': lambda: bcc110(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'bcc111': lambda: bcc111(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'hcp0001': lambda: hcp0001(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'hcp10m10': lambda: hcp10m10(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'diamond100': lambda: diamond100(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
            'diamond111': lambda: diamond111(element, size=build_size, a=lattice_const, vacuum=vacuum, orthogonal=orthogonal),
        }
        return surface_builders_with_layers[surface_type]()
    else:
        raise ValueError(f"Unknown surface type: {surface_type}")


class SurfaceDialog(ModalScreen):
    """Modal dialog for generating surface/film structures."""
    
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    
    CSS = """
    #surface-dialog {
        width: 65;
        height: auto;
        background: $panel;
        border: thick $primary;
        padding: 1 2;
    }
    
    #surface-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 1;
        text-align: center;
    }
    
    .surface-row {
        height: auto;
        margin-bottom: 1;
    }
    
    .surface-row Label {
        width: 20;
    }
    
    .surface-row Input {
        width: 1fr;
    }
    
    .surface-row Select {
        width: 1fr;
    }
    
    #miller-inputs {
        display: none;
    }
    
    #miller-inputs.visible {
        display: block;
    }
    
    #surface-info {
        margin: 1 0;
        color: $text-muted;
    }
    
    #surface-buttons {
        margin-top: 1;
        height: auto;
    }
    """
    
    def __init__(self):
        super().__init__()
        self._show_miller = False
    
    def compose(self) -> ComposeResult:
        # Create element options from known elements
        element_options = [(f"{sym} ({struct}, a={a:.2f}Å)", sym) 
                          for sym, (a, struct) in ELEMENT_LATTICE_CONSTANTS.items()]
        
        with Container(id="surface-dialog"):
            yield Label("Generate Surface/Film", id="surface-title")
            yield Rule()
            
            with Horizontal(classes="surface-row"):
                yield Label("Element:")
                yield Select(element_options, value="Fe", id="select-element")
            
            with Horizontal(classes="surface-row"):
                yield Label("Surface type:")
                yield Select(
                    [(name, value) for value, name in SURFACE_PRESETS],
                    value="bcc100",
                    id="select-surface-type"
                )
            
            # Miller indices for custom surface
            with Vertical(id="miller-inputs"):
                yield Label("Miller indices (h k l):", classes="surface-section")
                with Horizontal(classes="surface-row"):
                    yield Label("h:")
                    yield Input(value="1", id="input-miller-h", type="integer")
                with Horizontal(classes="surface-row"):
                    yield Label("k:")
                    yield Input(value="0", id="input-miller-k", type="integer")
                with Horizontal(classes="surface-row"):
                    yield Label("l:")
                    yield Input(value="0", id="input-miller-l", type="integer")
            
            yield Rule()
            
            with Horizontal(classes="surface-row"):
                yield Label("Lattice constant (Å):")
                yield Input(value="2.87", id="input-lattice", type="number")
            
            with Horizontal(classes="surface-row"):
                yield Label("Number of layers:")
                yield Input(value="5", id="input-layers", type="integer")
            
            with Horizontal(classes="surface-row"):
                yield Label("Vacuum (Å):")
                yield Input(value="10.0", id="input-vacuum", type="number")
            
            yield Rule()
            yield Label("Supercell size (in-plane):", classes="surface-section")
            
            with Horizontal(classes="surface-row"):
                yield Label("Size x:")
                yield Input(value="1", id="input-size-x", type="integer")
            
            with Horizontal(classes="surface-row"):
                yield Label("Size y:")
                yield Input(value="1", id="input-size-y", type="integer")
            
            yield Static("", id="surface-info")
            
            with Horizontal(id="surface-buttons"):
                yield Button("Generate", variant="primary", id="btn-generate-surface")
                yield Button("Cancel", variant="default", id="btn-cancel-surface")
    
    @on(Select.Changed, "#select-element")
    def on_element_changed(self, event: Select.Changed):
        """Update lattice constant when element changes."""
        element = event.value
        if element in ELEMENT_LATTICE_CONSTANTS:
            lattice, struct = ELEMENT_LATTICE_CONSTANTS[element]
            self.query_one("#input-lattice", Input).value = str(lattice)
            
            # Update surface type to match structure
            surface_select = self.query_one("#select-surface-type", Select)
            if struct == 'bcc':
                surface_select.value = 'bcc100'
            elif struct == 'fcc':
                surface_select.value = 'fcc100'
            elif struct == 'hcp':
                surface_select.value = 'hcp0001'
            elif struct == 'diamond':
                surface_select.value = 'diamond100'
    
    @on(Select.Changed, "#select-surface-type")
    def on_surface_type_changed(self, event: Select.Changed):
        """Show/hide Miller indices for custom surface."""
        is_custom = event.value == "custom"
        self._show_miller = is_custom
        
        miller_inputs = self.query_one("#miller-inputs", Vertical)
        if is_custom:
            miller_inputs.add_class("visible")
        else:
            miller_inputs.remove_class("visible")
    
    @on(Button.Pressed, "#btn-generate-surface")
    def on_generate(self):
        try:
            element = self.query_one("#select-element", Select).value
            surface_type = self.query_one("#select-surface-type", Select).value
            lattice = float(self.query_one("#input-lattice", Input).value or "2.87")
            layers = int(self.query_one("#input-layers", Input).value or "5")
            vacuum = float(self.query_one("#input-vacuum", Input).value or "10.0")
            size_x = int(self.query_one("#input-size-x", Input).value or "1")
            size_y = int(self.query_one("#input-size-y", Input).value or "1")
            
            miller = None
            if self._show_miller:
                h = int(self.query_one("#input-miller-h", Input).value or "1")
                k = int(self.query_one("#input-miller-k", Input).value or "0")
                l = int(self.query_one("#input-miller-l", Input).value or "0")
                miller = (h, k, l)
            
            self.dismiss({
                'element': element,
                'surface_type': surface_type,
                'lattice_const': lattice,
                'layers': layers,
                'vacuum': vacuum,
                'size': (size_x, size_y, 1),
                'miller': miller,
            })
        except ValueError as e:
            self.notify(f"Invalid input: {e}", severity="error")
    
    @on(Button.Pressed, "#btn-cancel-surface")
    def on_cancel(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)
