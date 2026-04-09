#!/usr/bin/env python3
"""
FLEUR Input File Analyzer

Parses FLEUR inp.xml files to extract key computational parameters
and estimate the number of basis functions for memory/resource planning.
"""

import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple


def _cross(a: List[float], b: List[float]) -> List[float]:
    """Cross product of two 3D vectors."""
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ]


def _dot(a: List[float], b: List[float]) -> float:
    """Dot product of two vectors."""
    return sum(x * y for x, y in zip(a, b))


@dataclass
class FleurAnalysisResult:
    """Results from analyzing a FLEUR inp.xml file."""
    
    # Raw extracted values
    kmax: float
    jspins: int
    has_soc: bool
    has_noco: bool
    has_inversion: bool
    num_symmetry_ops: int
    
    # Cell information
    cell_vectors: List[List[float]]  # 3x3 matrix of lattice vectors
    cell_volume: float  # in Bohr^3
    lattice_scale: float
    
    # Atom information
    num_atoms: int
    num_atom_types: int
    atom_species: List[str]
    
    # K-point information
    num_kpoints: int
    
    # Calculated estimates
    n_basis_per_kpoint: int  # Pure basis function count
    matrix_dimension: int  # Matrix dimension (n_basis × 2 for noco)
    is_complex: bool  # True if complex storage (no inversion or noco)
    
    # File info
    file_path: str
    
    def __str__(self) -> str:
        """Human-readable summary."""
        lines = [
            "=" * 50,
            "FLEUR Input Analysis",
            "=" * 50,
            f"File: {self.file_path}",
            "",
            "--- Computational Parameters ---",
            f"  Kmax: {self.kmax:.2f} a.u.⁻¹",
            f"  Cell volume: {self.cell_volume:.2f} Bohr³",
            f"  Spins (jspins): {self.jspins} ({'spin-polarized' if self.jspins == 2 else 'non-magnetic'})",
            f"  SOC enabled: {'Yes' if self.has_soc else 'No'}",
            f"  Non-collinear: {'Yes' if self.has_noco else 'No'}",
            f"  Has inversion: {'Yes' if self.has_inversion else 'No'}",
            f"  Symmetry ops: {self.num_symmetry_ops}",
            f"  K-points: {self.num_kpoints}",
            "",
            "--- Atoms ---",
            f"  Total atoms: {self.num_atoms}",
            f"  Atom types: {self.num_atom_types}",
            f"  Species: {', '.join(self.atom_species)}",
            "",
            "--- Basis Function Estimate ---",
            f"  Base formula: Kmax³ × Volume / (6π²)",
            f"  N_basis: {self.n_basis_per_kpoint:,}",
            f"  Noco factor: ×2" if self.has_noco else "",
            f"  >>> Matrix dimension: {self.matrix_dimension:,}",
            "",
            "--- Memory Estimates (per k-point) ---",
            f"  Storage type: {'complex' if self.is_complex else 'real'} {'(no inversion or noco)' if self.is_complex else '(has inversion, no noco)'}",
            f"  Hamiltonian matrix: {self.matrix_dimension:,} × {self.matrix_dimension:,}",
            f"  Memory: {self._format_memory(self.matrix_dimension ** 2 * (16 if self.is_complex else 8))}",
            "=" * 50,
        ]
        return "\n".join(lines)
    
    @staticmethod
    def _format_memory(bytes_count: int) -> str:
        """Format bytes to human-readable size."""
        for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
            if bytes_count < 1024:
                return f"{bytes_count:.2f} {unit}"
            bytes_count /= 1024
        return f"{bytes_count:.2f} PB"
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON export."""
        return {
            "kmax": self.kmax,
            "jspins": self.jspins,
            "has_soc": self.has_soc,
            "has_noco": self.has_noco,
            "has_inversion": self.has_inversion,
            "num_symmetry_ops": self.num_symmetry_ops,
            "cell_volume": self.cell_volume,
            "lattice_scale": self.lattice_scale,
            "num_atoms": self.num_atoms,
            "num_atom_types": self.num_atom_types,
            "atom_species": self.atom_species,
            "num_kpoints": self.num_kpoints,
            "n_basis_per_kpoint": self.n_basis_per_kpoint,
            "matrix_dimension": self.matrix_dimension,
            "is_complex": self.is_complex,
            "file_path": self.file_path,
        }


class FleurInputAnalyzer:
    """Analyzer for FLEUR inp.xml files."""
    
    def __init__(self, file_path: str):
        """Initialize with path to inp.xml file."""
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        self.tree = ET.parse(self.file_path)
        self.root = self.tree.getroot()
    
    def analyze(self) -> FleurAnalysisResult:
        """Perform full analysis of the inp.xml file."""
        # Extract all parameters
        kmax = self._get_kmax()
        jspins = self._get_jspins()
        has_soc = self._get_has_soc()
        has_noco = self._get_has_noco()
        
        cell_vectors, scale = self._get_cell_vectors()
        cell_volume = self._calculate_volume(cell_vectors)
        
        has_inversion = self._check_inversion_symmetry()
        num_sym_ops = self._count_symmetry_ops()
        
        num_atoms, num_types, species = self._get_atom_info()
        num_kpoints = self._get_num_kpoints()
        
        # Calculate basis function estimate
        # Formula: N = (4/3) * π * kmax^3 * V / (2π)^3 = kmax^3 * V / (6π²)
        n_basis = int((kmax ** 3) * cell_volume / (6.0 * math.pi ** 2))
        
        # Matrix dimension: noco doubles it (2-component spinors)
        matrix_dimension = n_basis * 2 if has_noco else n_basis
        
        # Storage type: complex if no inversion OR noco, otherwise real
        is_complex = (not has_inversion) or has_noco
        
        return FleurAnalysisResult(
            kmax=kmax,
            jspins=jspins,
            has_soc=has_soc,
            has_noco=has_noco,
            has_inversion=has_inversion,
            num_symmetry_ops=num_sym_ops,
            cell_vectors=cell_vectors,
            cell_volume=cell_volume,
            lattice_scale=scale,
            num_atoms=num_atoms,
            num_atom_types=num_types,
            atom_species=species,
            num_kpoints=num_kpoints,
            n_basis_per_kpoint=n_basis,
            matrix_dimension=matrix_dimension,
            is_complex=is_complex,
            file_path=str(self.file_path),
        )
    
    def _get_kmax(self) -> float:
        """Extract Kmax from cutoffs element."""
        cutoffs = self.root.find(".//cutoffs")
        if cutoffs is None:
            raise ValueError("No <cutoffs> element found in inp.xml")
        
        kmax_str = cutoffs.get("Kmax")
        if kmax_str is None:
            raise ValueError("No Kmax attribute in <cutoffs>")
        
        return float(kmax_str)
    
    def _get_jspins(self) -> int:
        """Extract jspins from magnetism element."""
        magnetism = self.root.find(".//magnetism")
        if magnetism is None:
            return 1  # Default to non-magnetic
        
        jspins_str = magnetism.get("jspins", "1")
        return int(jspins_str)
    
    def _get_has_soc(self) -> bool:
        """Check if spin-orbit coupling is enabled."""
        soc = self.root.find(".//soc")
        if soc is None:
            return False
        
        l_soc = soc.get("l_soc", "F")
        return l_soc.upper() in ("T", "TRUE", ".TRUE.")
    
    def _get_has_noco(self) -> bool:
        """Check if non-collinear mode is enabled."""
        # l_noco is an attribute on the <magnetism> element
        magnetism = self.root.find(".//magnetism")
        if magnetism is not None:
            l_noco = magnetism.get("l_noco", "F")
            if l_noco.upper() in ("T", "TRUE", ".TRUE."):
                return True
        
        # Also check if SOC implies noco
        return self._get_has_soc()
    
    def _get_cell_vectors(self) -> Tuple[List[List[float]], float]:
        """Extract cell vectors from bulkLattice or filmLattice."""
        # Try bulkLattice first
        bulk_lattice = self.root.find(".//bulkLattice")
        film_lattice = self.root.find(".//filmLattice")
        
        if bulk_lattice is not None:
            return self._get_bulk_vectors(bulk_lattice)
        elif film_lattice is not None:
            return self._get_film_vectors(film_lattice)
        else:
            raise ValueError("No bulkLattice or filmLattice found")
    
    def _get_bulk_vectors(self, lattice) -> Tuple[List[List[float]], float]:
        """Extract cell vectors from bulkLattice."""
        scale = float(lattice.get("scale", "1.0"))
        
        bravais = lattice.find("bravaisMatrix")
        if bravais is None:
            raise ValueError("No bravaisMatrix found in bulkLattice")
        
        vectors: List[List[float]] = []
        for row_name in ["row-1", "row-2", "row-3"]:
            row_elem = bravais.find(row_name)
            if row_elem is None:
                raise ValueError(f"Missing {row_name} in bravaisMatrix")
            
            values = row_elem.text.strip().split()
            row = [float(v) * scale for v in values[:3]]
            vectors.append(row)
        
        return vectors, scale
    
    def _get_film_vectors(self, lattice) -> Tuple[List[List[float]], float]:
        """
        Extract cell vectors from filmLattice.
        
        Film lattices have a 2D bravaisMatrixFilm and use dTilda for the 
        third direction (perpendicular to the film).
        """
        scale = float(lattice.get("scale", "1.0"))
        dTilda = float(lattice.get("dTilda", "1.0"))
        
        bravais = lattice.find("bravaisMatrixFilm")
        if bravais is None:
            raise ValueError("No bravaisMatrixFilm found in filmLattice")
        
        vectors: List[List[float]] = []
        
        # Read the 2D lattice vectors (row-1 and row-2)
        for row_name in ["row-1", "row-2"]:
            row_elem = bravais.find(row_name)
            if row_elem is None:
                raise ValueError(f"Missing {row_name} in bravaisMatrixFilm")
            
            values = row_elem.text.strip().split()
            # 2D vectors: extend to 3D with z=0
            row = [float(values[0]) * scale, float(values[1]) * scale, 0.0]
            vectors.append(row)
        
        # Third vector is perpendicular to the film: (0, 0, dTilda)
        vectors.append([0.0, 0.0, dTilda * scale])
        
        return vectors, scale
    
    def _calculate_volume(self, vectors: List[List[float]]) -> float:
        """Calculate cell volume from lattice vectors."""
        # Volume = |a · (b × c)|
        a, b, c = vectors[0], vectors[1], vectors[2]
        volume = abs(_dot(a, _cross(b, c)))
        return volume
    
    def _check_inversion_symmetry(self) -> bool:
        """Check if inversion symmetry (-I) is present."""
        sym_ops = self.root.find(".//symmetryOperations")
        if sym_ops is None:
            return False
        
        for sym_op in sym_ops.findall("symOp"):
            # Extract the 3x3 rotation matrix
            matrix: List[List[float]] = []
            try:
                for row_name in ["row-1", "row-2", "row-3"]:
                    row_elem = sym_op.find(row_name)
                    if row_elem is None:
                        continue
                    values = row_elem.text.strip().split()
                    matrix.append([float(v) for v in values[:3]])
                
                if len(matrix) != 3:
                    continue
                
                # Check if this is -I (inversion: diagonal = -1, off-diagonal = 0)
                is_inversion = True
                for i in range(3):
                    for j in range(3):
                        expected = -1.0 if i == j else 0.0
                        if abs(matrix[i][j] - expected) > 1e-6:
                            is_inversion = False
                            break
                    if not is_inversion:
                        break
                
                if is_inversion:
                    return True
            except (ValueError, TypeError):
                continue
        
        return False
    
    def _count_symmetry_ops(self) -> int:
        """Count the number of symmetry operations."""
        sym_ops = self.root.find(".//symmetryOperations")
        if sym_ops is None:
            return 1  # Identity only
        
        return len(sym_ops.findall("symOp"))
    
    def _get_atom_info(self) -> Tuple[int, int, List[str]]:
        """Get atom count, type count, and species names."""
        atom_groups = self.root.find(".//atomGroups")
        if atom_groups is None:
            return 0, 0, []
        
        num_atoms = 0
        species_list = []
        
        for group in atom_groups.findall("atomGroup"):
            species = group.get("species", "unknown")
            if species not in species_list:
                species_list.append(species)
            
            # Count atoms in this group (relPos or absPos elements)
            for pos_type in ["relPos", "absPos"]:
                num_atoms += len(group.findall(pos_type))
        
        return num_atoms, len(species_list), species_list
    
    def _get_num_kpoints(self) -> int:
        """Get the number of k-points from the active k-point list."""
        # Find the selected k-point list name
        selection = self.root.find(".//kPointListSelection")
        if selection is not None:
            list_name = selection.get("listName", "")
        else:
            list_name = ""
        
        # Find k-point lists
        kpoint_lists = self.root.find(".//kPointLists")
        if kpoint_lists is None:
            return 1
        
        # Find the matching list or use the first one
        for kpl in kpoint_lists.findall("kPointList"):
            name = kpl.get("name", "")
            if name == list_name or not list_name:
                count = kpl.get("count")
                if count:
                    return int(count)
                # If no count attribute, count the kPoint elements
                return len(kpl.findall("kPoint"))
        
        return 1


def analyze_fleur_input(file_path: str) -> FleurAnalysisResult:
    """Convenience function to analyze a FLEUR inp.xml file."""
    analyzer = FleurInputAnalyzer(file_path)
    return analyzer.analyze()


# CLI interface
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python fleur_analyzer.py <inp.xml>")
        sys.exit(1)
    
    result = analyze_fleur_input(sys.argv[1])
    print(result)
