"""
K-Point Manager for FLEUR inp.xml files
Provides functionality to read, write, and modify k-point lists in inp.xml
"""

import xml.etree.ElementTree as ET
from typing import List, Dict, Optional, Tuple, Union
from pathlib import Path
from dataclasses import dataclass, field
from enum import Enum

from .inpgen_loader import create_inpgen_interface


class KPointMode(Enum):
    """K-Point generation modes."""
    GRID = "grid"
    DENSITY = "den"
    NUMBER = "nk"
    BAND = "band"


@dataclass
class KPointModifiers:
    """Modifiers for k-point generation."""
    gamma: bool = False  # Include gamma point
    gauss: bool = False  # Use Gaussian integration
    tria: bool = False   # Use triangular method
    tetra: bool = False  # Use tetrahedral method
    soc: bool = False    # Compatible with spin-orbit coupling
    
    def to_string(self) -> str:
        """Convert modifiers to string for inpgen."""
        parts = []
        if self.gamma:
            parts.append("gamma")
        if self.gauss:
            parts.append("gauss")
        if self.tria:
            parts.append("tria")
        if self.tetra:
            parts.append("tetra")
        if self.soc:
            parts.append("soc")
        return "@".join(parts) if parts else ""


@dataclass
class KPoint:
    """Represents a single k-point with coordinates and weight."""
    x: Union[float, str]
    y: Union[float, str]
    z: Union[float, str]
    weight: Union[float, str]
    label: str = ""
    
    def to_xml(self, parent: ET.Element) -> ET.Element:
        """Create XML element for this k-point."""
        kpoint = ET.SubElement(parent, "kPoint")
        kpoint.text = f"{self.x} {self.y} {self.z}"
        kpoint.set("weight", str(self.weight))
        if self.label:
            kpoint.set("label", self.label)
        return kpoint
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'KPoint':
        """Create KPoint from XML element."""
        coord_strings = element.text.strip().split()
        if len(coord_strings) != 3:
            raise ValueError(f"Expected 3 coordinates, got {len(coord_strings)}")
        
        # Try to parse as float, otherwise keep as string
        coords = []
        for coord_str in coord_strings:
            try:
                coords.append(float(coord_str))
            except ValueError:
                coords.append(coord_str)
        
        weight_str = element.get("weight", "1.0")
        try:
            weight = float(weight_str)
        except ValueError:
            weight = weight_str
        
        label = element.get("label", "")
        
        return cls(x=coords[0], y=coords[1], z=coords[2], weight=weight, label=label)


@dataclass
class KPointList:
    """Represents a list of k-points with metadata."""
    name: str
    kpoints: List[KPoint] = field(default_factory=list)
    type: str = "unspecified"
    count: Optional[int] = None
    nx: Optional[int] = None
    ny: Optional[int] = None
    nz: Optional[int] = None
    nkq_pairs: Optional[int] = None
    
    def to_xml(self, parent: ET.Element) -> ET.Element:
        """Create XML element for this k-point list."""
        kplist = ET.SubElement(parent, "kPointList")
        kplist.set("name", self.name)
        kplist.set("type", self.type)
        
        if self.count is not None:
            kplist.set("count", str(self.count))
        if self.nx is not None:
            kplist.set("nx", str(self.nx))
        if self.ny is not None:
            kplist.set("ny", str(self.ny))
        if self.nz is not None:
            kplist.set("nz", str(self.nz))
        if self.nkq_pairs is not None:
            kplist.set("nkq_pairs", str(self.nkq_pairs))
        
        for kpoint in self.kpoints:
            kpoint.to_xml(kplist)
        
        return kplist
    
    @classmethod
    def from_xml(cls, element: ET.Element) -> 'KPointList':
        """Create KPointList from XML element."""
        name = element.get("name", "default")
        type_ = element.get("type", "unspecified")
        count = int(element.get("count")) if element.get("count") else None
        nx = int(element.get("nx")) if element.get("nx") else None
        ny = int(element.get("ny")) if element.get("ny") else None
        nz = int(element.get("nz")) if element.get("nz") else None
        nkq_pairs = int(element.get("nkq_pairs")) if element.get("nkq_pairs") else None
        
        kpoints = []
        for kp_elem in element.findall("kPoint"):
            kpoints.append(KPoint.from_xml(kp_elem))
        
        return cls(
            name=name,
            kpoints=kpoints,
            type=type_,
            count=count,
            nx=nx,
            ny=ny,
            nz=nz,
            nkq_pairs=nkq_pairs
        )
    
    def add_kpoint(self, kpoint: KPoint):
        """Add a k-point to this list."""
        self.kpoints.append(kpoint)
        if self.count is not None:
            self.count += 1
    
    def remove_kpoint(self, index: int) -> Optional[KPoint]:
        """Remove k-point at given index."""
        if 0 <= index < len(self.kpoints):
            kpoint = self.kpoints.pop(index)
            if self.count is not None:
                self.count -= 1
            return kpoint
        return None


def _process_xincludes(elem: ET.Element, base_dir: Path) -> None:
    """Process XInclude directives in-place.

    Each ``<xi:include href="..."/>`` element whose referenced file exists is
    replaced with the parsed root of that file.  If the file cannot be found
    (or the href is a remote URL / text include) the xi:include element is left
    untouched so that round-trip serialisation preserves it.
    """
    XI_INCLUDE = "{http://www.w3.org/2001/XInclude}include"

    i = 0
    while i < len(elem):
        child = elem[i]
        if child.tag == XI_INCLUDE:
            href = child.get("href", "")
            parse = child.get("parse", "xml")
            replaced = False
            # Only handle local XML includes (skip URLs and text-mode includes)
            if href and parse == "xml" and not href.startswith(("http://", "https://", "ftp://")):
                include_path = base_dir / href
                if include_path.exists():
                    try:
                        included_root = ET.parse(str(include_path)).getroot()
                        # Preserve the whitespace tail that was on the xi:include tag
                        included_root.tail = child.tail
                        elem.remove(child)
                        elem.insert(i, included_root)
                        # Recurse so nested xi:includes are also resolved
                        _process_xincludes(included_root, include_path.parent)
                        replaced = True
                    except ET.ParseError:
                        pass  # keep the xi:include element on parse failure
            if not replaced:
                i += 1
        else:
            _process_xincludes(child, base_dir)
            i += 1


class InpXMLManager:
    """Manager for FLEUR inp.xml files focusing on k-point operations."""
    
    def __init__(self, xml_path: str = "inp.xml", quiet: bool = False):
        """Initialize manager with path to inp.xml file.
        
        Args:
            xml_path: Path to the inp.xml file
            quiet: If True, suppress console output (messages stored in last_messages)
        """
        self.xml_path = Path(xml_path)
        self.tree: Optional[ET.ElementTree] = None
        self.root: Optional[ET.Element] = None
        self.kpoint_lists: Dict[str, KPointList] = {}
        self.selected_list: Optional[str] = None
        self.quiet = quiet
        self.last_messages: List[str] = []  # Messages from last operation
        
        if self.xml_path.exists():
            self.load()
    
    def _log(self, message: str):
        """Log a message. Stores in last_messages and optionally prints."""
        self.last_messages.append(message)
        if not self.quiet:
            print(message)
    
    def get_last_messages(self) -> str:
        """Get messages from last operation as a single string."""
        return "\n".join(self.last_messages)
    
    def clear_messages(self):
        """Clear collected messages."""
        self.last_messages.clear()
    
    def load(self) -> None:
        """Load and parse inp.xml file."""
        if not self.xml_path.exists():
            raise FileNotFoundError(f"File not found: {self.xml_path}")
        
        self.tree = ET.parse(self.xml_path)
        # Resolve xi:include directives; missing files are left as xi:include nodes
        _process_xincludes(self.tree.getroot(), self.xml_path.parent)
        self.root = self.tree.getroot()
        self._parse_kpoints()
    
    def _parse_kpoints(self) -> None:
        """Parse k-points from loaded XML."""
        self.kpoint_lists.clear()
        
        # Find kPointLists element
        cell = self.root.find(".//cell")
        if cell is None:
            return
        
        bz_integration = cell.find("bzIntegration")
        if bz_integration is None:
            return
        
        # Get selected k-point list
        kpoint_selection = bz_integration.find("kPointListSelection")
        if kpoint_selection is not None:
            self.selected_list = kpoint_selection.get("listName")
        
        # Parse all k-point lists
        kpoint_lists_elem = bz_integration.find("kPointLists")
        if kpoint_lists_elem is not None:
            for kplist_elem in kpoint_lists_elem.findall("kPointList"):
                kplist = KPointList.from_xml(kplist_elem)
                self.kpoint_lists[kplist.name] = kplist
    
    def save(self, output_path: Optional[str] = None) -> None:
        """Save XML back to file."""
        if self.tree is None or self.root is None:
            raise RuntimeError("No XML loaded. Load or create XML first.")
        
        # Update k-points in XML
        self._update_kpoints_in_xml()
        
        # Write to file
        save_path = Path(output_path) if output_path else self.xml_path
        self._pretty_print(self.root)
        self.tree.write(save_path, encoding='utf-8', xml_declaration=True)
    
    def _update_kpoints_in_xml(self) -> None:
        """Update k-point data in XML tree."""
        cell = self.root.find(".//cell")
        if cell is None:
            raise ValueError("No cell element found in XML")
        
        bz_integration = cell.find("bzIntegration")
        if bz_integration is None:
            bz_integration = ET.SubElement(cell, "bzIntegration")
        
        # Update selection
        kpoint_selection = bz_integration.find("kPointListSelection")
        if kpoint_selection is None and self.selected_list:
            kpoint_selection = ET.SubElement(bz_integration, "kPointListSelection")
        if kpoint_selection is not None and self.selected_list:
            kpoint_selection.set("listName", self.selected_list)
        
        # Update lists
        kpoint_lists_elem = bz_integration.find("kPointLists")
        if kpoint_lists_elem is None:
            kpoint_lists_elem = ET.SubElement(bz_integration, "kPointLists")
        
        # Clear and rebuild k-point lists
        for elem in list(kpoint_lists_elem):
            kpoint_lists_elem.remove(elem)
        
        for kplist in self.kpoint_lists.values():
            kplist.to_xml(kpoint_lists_elem)
    
    def _pretty_print(self, elem: ET.Element, level: int = 0) -> None:
        """Add indentation for pretty printing."""
        indent = "\n" + "  " * level
        if len(elem):
            if not elem.text or not elem.text.strip():
                elem.text = indent + "  "
            if not elem.tail or not elem.tail.strip():
                elem.tail = indent
            for child in elem:
                self._pretty_print(child, level + 1)
            if not child.tail or not child.tail.strip():
                child.tail = indent
        else:
            if level and (not elem.tail or not elem.tail.strip()):
                elem.tail = indent
    
    def list_kpoint_lists(self) -> List[Tuple[str, int, str]]:
        """Get list of k-point list names with metadata.
        
        Returns:
            List of tuples (name, count, type)
        """
        result = []
        for name, kplist in self.kpoint_lists.items():
            count = len(kplist.kpoints)
            selected = " [SELECTED]" if name == self.selected_list else ""
            result.append((name, count, kplist.type, selected))
        return result
    
    def get_kpoint_list(self, name: str) -> Optional[KPointList]:
        """Get k-point list by name."""
        return self.kpoint_lists.get(name)
    
    def add_kpoint_list(self, kplist: KPointList) -> None:
        """Add a new k-point list."""
        self.kpoint_lists[kplist.name] = kplist
    
    def remove_kpoint_list(self, name: str) -> bool:
        """Remove k-point list by name."""
        if name in self.kpoint_lists:
            del self.kpoint_lists[name]
            if self.selected_list == name:
                self.selected_list = None
            return True
        return False
    
    def select_kpoint_list(self, name: str) -> bool:
        """Set the active k-point list."""
        if name in self.kpoint_lists:
            self.selected_list = name
            return True
        return False
    
    def get_selected_list(self) -> Optional[str]:
        """Get the name of the currently selected (active) k-point list."""
        return self.selected_list
    
    def add_kpoint_to_list(self, list_name: str, kpoint: KPoint) -> bool:
        """Add k-point to specific list."""
        kplist = self.kpoint_lists.get(list_name)
        if kplist:
            kplist.add_kpoint(kpoint)
            return True
        return False
    
    def create_kpoints(self, name: str, mode: KPointMode, 
                      modifiers: Optional[KPointModifiers] = None,
                      grid: Optional[Tuple[int, int, int]] = None,
                      density: Optional[float] = None,
                      num_kpoints: Optional[int] = None,
                      use_symmetry: bool = True,
                      lib_path: Optional[str] = None) -> KPointList:
        """Create k-points using FLEUR inpgen library with flexible modes.
        
        Args:
            name: Name for the k-point list (must be unique)
            mode: K-point generation mode (GRID, DENSITY, NUMBER, or BAND)
            modifiers: Optional KPointModifiers for gamma, gauss, tria, etc.
            grid: For GRID mode: (nx, ny, nz) divisions
            density: For DENSITY mode: k-point density
            num_kpoints: For NUMBER or BAND mode: number of k-points
            use_symmetry: If True, k-points will be reduced by symmetry
            lib_path: Optional path to inpgen library
        
        Returns:
            KPointList with generated k-points
        
        Raises:
            ValueError: If parameters are invalid or name already exists
        
        Examples:
            # Grid mode with gamma-centered mesh
            create_kpoints("mesh", KPointMode.GRID, 
                         modifiers=KPointModifiers(gamma=True),
                         grid=(8, 8, 8))
            
            # Density mode with triangular integration
            create_kpoints("dense", KPointMode.DENSITY,
                         modifiers=KPointModifiers(tria=True),
                         density=0.1)
            
            # Number mode requesting ~100 k-points
            create_kpoints("nk100", KPointMode.NUMBER,
                         num_kpoints=100)
        """
        # Check that the name is unique
        if name in self.kpoint_lists:
            raise ValueError(f"K-point list with name '{name}' already exists. "
                           f"Existing lists: {list(self.kpoint_lists.keys())}")
        
        # Set default modifiers
        if modifiers is None:
            modifiers = KPointModifiers()
        
        # Build the k-point specification string: name#modifier@mode=details
        parts = [name]
        
        # Add modifiers
        modifier_str = modifiers.to_string()
        if modifier_str:
            parts.append(modifier_str)
        
        # Add mode and details
        if mode == KPointMode.GRID:
            if grid is None:
                raise ValueError("grid parameter required for GRID mode")
            nx, ny, nz = grid
            mode_str = f"grid={nx},{ny},{nz}"
        elif mode == KPointMode.DENSITY:
            if density is None:
                raise ValueError("density parameter required for DENSITY mode")
            mode_str = f"den={density}"
        elif mode == KPointMode.NUMBER:
            if num_kpoints is None:
                raise ValueError("num_kpoints parameter required for NUMBER mode")
            mode_str = f"nk={num_kpoints}"
        elif mode == KPointMode.BAND:
            if num_kpoints is None:
                raise ValueError("num_kpoints parameter required for BAND mode")
            mode_str = f"band={num_kpoints}"
        else:
            raise ValueError(f"Unknown mode: {mode}")
        
        parts.append(mode_str)
        
        # Join with appropriate separators: name#modifier@mode
        if len(parts) == 3:  # name, modifier, mode
            kpts_str = f"{parts[0]}#{parts[1]}@{parts[2]}"
        else:  # name, mode (no modifiers)
            kpts_str = f"{parts[0]}#{parts[1]}"
        
        # In FLEUR inpgen: nosym=True means no symmetry reduction
        nosym = not use_symmetry
        
        # Clear messages for this operation
        self.clear_messages()
        
        try:
            inpgen = create_inpgen_interface(lib_path=lib_path, quiet=True)
            inpgen.add_kpoints(kpts_str, kpts_path="", nosym=nosym)
            # Collect messages from inpgen
            for msg in inpgen.messages:
                self._log(msg)
            symmetry_msg = "without symmetry" if nosym else "with symmetry"
            self._log(f"✓ Generated k-points ({mode.value}) {symmetry_msg} using inpgen library")
        except Exception as e:
            self._log(f"✗ Error calling inpgen library: {e}")
            raise RuntimeError(f"Cannot generate k-points:\n{e}") from e
        
        # Reload the XML to get the updated k-points
        self.load()
        
        # Find the newly created k-point list
        return self._find_or_rename_kplist(name)
    
    def create_mesh_kpoints(self, name: str, nx: int, ny: int, nz: int, 
                           gamma_centered: bool = True, use_symmetry: bool = True,
                           lib_path: Optional[str] = None) -> KPointList:
        """Create a k-point mesh using FLEUR inpgen library.
        
        This is a convenience wrapper around create_kpoints() for backward compatibility.
        
        Args:
            name: Name for the k-point list (must be unique)
            nx, ny, nz: Number of divisions in each direction
            gamma_centered: If True, mesh is centered at Gamma point (0,0,0)
            use_symmetry: If True, k-points will be reduced by symmetry
            lib_path: Optional path to inpgen library
        
        Returns:
            KPointList with generated k-points
        """
        modifiers = KPointModifiers(gamma=gamma_centered)
        return self.create_kpoints(
            name=name,
            mode=KPointMode.GRID,
            modifiers=modifiers,
            grid=(nx, ny, nz),
            use_symmetry=use_symmetry,
            lib_path=lib_path
        )
    
    def _find_or_rename_kplist(self, name: str) -> KPointList:
        """Find newly created k-point list or rename from 'default'."""
        # The inpgen library typically creates a list named "default" or similar
        # We need to rename it if a specific name was requested
        if name != "default" and "default" in self.kpoint_lists:
            kplist = self.kpoint_lists["default"]
            del self.kpoint_lists["default"]
            kplist.name = name
            self.kpoint_lists[name] = kplist
            
            # Update the selected list if it was "default"
            if self.selected_list == "default":
                self.selected_list = name
            
            return kplist
        elif name in self.kpoint_lists:
            return self.kpoint_lists[name]
        else:
            # If we can't find the list, return the first one
            if self.kpoint_lists:
                return list(self.kpoint_lists.values())[0]
            else:
                raise RuntimeError("No k-point list created by inpgen library")
    
    def create_kpoint_path(self, name: str, special_points: List[Tuple[str, Tuple[float, float, float]]], 
                          num_points: int = 100, path_string: Optional[str] = None,
                          lib_path: Optional[str] = None) -> KPointList:
        """Create a k-point path for band structure calculations using FLEUR inpgen library.
        
        Args:
            name: Name for the k-point list (must be unique)
            special_points: List of (label, (x, y, z)) tuples defining the path.
                          If path_string is provided, this is used to define custom points.
            num_points: Total number of k-points to generate along the path
            path_string: Optional path specification string (e.g., "G-X-W-L-G").
                       If None, will be constructed from special_points labels.
                       Can also use format "gamma=0,0,0;x=0.5,0.5,0.5" for custom points.
            lib_path: Optional path to inpgen library
        
        Returns:
            KPointList with path k-points
        
        Raises:
            ValueError: If a k-point list with the given name already exists or
                       if fewer than 2 special points are provided
        
        Note: This uses the FLEUR inpgen library to generate k-point paths.
              The inp.xml file is modified in place.
              Band mode does not support modifiers.
        
        Examples:
            # Using pre-defined special points
            create_kpoint_path("bands", [
                ("Γ", (0.0, 0.0, 0.0)),
                ("X", (0.5, 0.0, 0.5)),
                ("W", (0.5, 0.25, 0.75)),
                ("L", (0.5, 0.5, 0.5)),
                ("Γ", (0.0, 0.0, 0.0))
            ], num_points=200)
            
            # Using custom path string with defined points
            create_kpoint_path("bands", [], num_points=200,
                path_string="gamma=0,0,0;x=0.5,0.5,0.5;l=0.5,0.5,0.5;gamma=0,0,0")
        """
        # Check that the name is unique
        if name in self.kpoint_lists:
            raise ValueError(f"K-point list with name '{name}' already exists. "
                           f"Existing lists: {list(self.kpoint_lists.keys())}")
        
        # Build path string
        if path_string is None:
            if len(special_points) < 2:
                raise ValueError("Need at least 2 special points to define a path")
            # Build path string from special points labels
            path_labels = "-".join(label for label, _ in special_points)
        else:
            path_labels = path_string
        
        # Build k-point specification string
        kpts_str = f"{name}#band={num_points}"
        
        # Clear messages for this operation
        self.clear_messages()
        
        # Call inpgen library to add k-points (band mode always uses nosym=True)
        try:
            inpgen = create_inpgen_interface(lib_path=lib_path, quiet=True)
            inpgen.add_kpoints(kpts_str, kpts_path=path_labels, nosym=True)
            # Collect messages from inpgen
            for msg in inpgen.messages:
                self._log(msg)
            self._log(f"✓ Generated k-point path '{path_labels}' with {num_points} points using inpgen library")
        except Exception as e:
            self._log(f"✗ Error calling inpgen library: {e}")
            if special_points:
                self._log(f"  Falling back to manual path generation...")
                # Fall back to manual generation if inpgen fails and we have special_points
                return self._create_kpoint_path_manual(name, special_points, num_points)
            else:
                raise RuntimeError(f"Cannot create path: inpgen failed and no special_points provided\n{e}") from e
        
        # Reload the XML to get the updated k-points
        self.load()
        
        # Find the newly created k-point list
        return self._find_or_rename_kplist(name)
    
    def _create_kpoint_path_manual(self, name: str, special_points: List[Tuple[str, Tuple[float, float, float]]], 
                                   num_points: int = 100) -> KPointList:
        """Fallback method to manually create a k-point path.
        
        This is used when the inpgen library is not available or fails.
        """
        import math
        
        kplist = KPointList(name=name, type="path")
        
        # Calculate total path length
        total_length = 0.0
        segment_lengths = []
        
        for i in range(len(special_points) - 1):
            p1 = special_points[i][1]
            p2 = special_points[i+1][1]
            length = math.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2)
            segment_lengths.append(length)
            total_length += length
        
        # Generate k-points along the path
        points_generated = 0
        
        for seg_idx, (start_pt, end_pt) in enumerate(zip(special_points[:-1], special_points[1:])):
            start_label, start_coords = start_pt
            end_label, end_coords = end_pt
            
            # Number of points for this segment (proportional to length)
            if seg_idx == len(special_points) - 2:
                # Last segment gets remaining points
                seg_points = num_points - points_generated
            else:
                seg_length = segment_lengths[seg_idx]
                seg_points = max(2, int(num_points * seg_length / total_length))
            
            # Generate points along this segment
            for i in range(seg_points):
                t = i / (seg_points - 1) if seg_points > 1 else 0.0
                x = start_coords[0] + t * (end_coords[0] - start_coords[0])
                y = start_coords[1] + t * (end_coords[1] - start_coords[1])
                z = start_coords[2] + t * (end_coords[2] - start_coords[2])
                
                # Label only at special points
                label = ""
                if i == 0:
                    label = start_label
                elif i == seg_points - 1 and seg_idx == len(special_points) - 2:
                    label = end_label
                
                kplist.add_kpoint(KPoint(x, y, z, weight=1.0, label=label))
                points_generated += 1
        
        kplist.count = len(kplist.kpoints)
        self.add_kpoint_list(kplist)
        return kplist
    
    def get_summary(self) -> Dict:
        """Get summary of k-point configuration."""
        return {
            "file": str(self.xml_path),
            "num_lists": len(self.kpoint_lists),
            "selected_list": self.selected_list,
            "total_kpoints": sum(len(kpl.kpoints) for kpl in self.kpoint_lists.values()),
            "lists": {name: len(kpl.kpoints) for name, kpl in self.kpoint_lists.items()}
        }
