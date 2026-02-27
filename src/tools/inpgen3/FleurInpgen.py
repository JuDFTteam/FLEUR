"""
Python interface for FLEUR inpgen library (inpgen_lib.F90)
Provides functions to generate inp.xml and add k-points from Python.
"""

import ctypes
import os
from pathlib import Path
from typing import Optional, Union, List


class InpgenInterface:
    """Interface to FLEUR inpgen Fortran library."""
    
    def __init__(self, lib_path: Optional[Union[str, Path]] = None, quiet: bool = False):
        """
        Initialize the inpgen interface.
        
        Args:
            lib_path: Path to the shared library (.so or .dylib).
                     If None, searches in common locations.
            quiet: If True, suppress console output (messages still collected in self.messages)
        """
        self.quiet = quiet
        self.messages: List[str] = []  # Collected log messages
        
        if lib_path is None:
            lib_path = self._find_library()
        
        self.lib = ctypes.CDLL(str(lib_path))
        self._setup_functions()
        self._log(f"✓ Loaded library: {lib_path}")
    
    def _log(self, message: str):
        """Log a message. Stores in self.messages and optionally prints."""
        self.messages.append(message)
        if not self.quiet:
            print(message)
    
    def get_messages(self) -> str:
        """Get all collected messages as a single string."""
        return "\n".join(self.messages)
    
    def clear_messages(self):
        """Clear the collected messages."""
        self.messages.clear()
    
    def _find_library(self) -> Path:
        """Search for the inpgen library in common locations."""
        possible_names = [
            'libfleurinpgen.so',
            'libfleurinpgen.dylib',
            'libinpgen3.so',
            'libinpgen3.dylib',
            'fleurinpgen.so',
            'fleurinpgen.dylib',
            'inpgen3.so',
            'inpgen3.dylib'
        ]
        
        search_paths = [
            Path.cwd(),
            Path.cwd() / 'build' / 'lib',
            Path.cwd() / 'build' / 'src' / 'tools' / 'inpgen3',
            Path.cwd() / 'lib',
            Path(__file__).parent / 'lib',
            Path(__file__).parent.parent.parent.parent / 'build' / 'lib',
        ]
        
        # Add build directory from environment (set by fleuriste wrapper)
        builddir = os.environ.get('FLEUR_BUILDDIR')
        if builddir:
            search_paths.insert(0, Path(builddir) / 'src' / 'tools' / 'inpgen3')
            search_paths.insert(1, Path(builddir))
        
        for path in search_paths:
            if not path.exists():
                continue
            for name in possible_names:
                lib_file = path / name
                if lib_file.exists():
                    return lib_file
        
        raise FileNotFoundError(
            f"Could not find inpgen library.\n"
            f"Searched in: {[str(p) for p in search_paths]}\n"
            f"Looking for: {possible_names}\n\n"
            "Please build the library first:\n"
            "  cd build && cmake .. && make inpgen3\n"
            "Or specify lib_path explicitly."
        )
    
    def _setup_functions(self):
        """Setup ctypes function signatures for Fortran routines."""
        
        # say_hello function - test function
        # SUBROUTINE say_hello() BIND(C, name="say_hello")
        self.lib.say_hello.argtypes = [ctypes.c_int]  # number of times to say hello
        self.lib.say_hello.restype = None
        
        # make_inp function
        # SUBROUTINE make_inp_py(len_simple_input, simple_input, len_profileName, 
        #                        profileName, nosym) BIND(C, name="make_inp")
        self.lib.make_inp.argtypes = [
            ctypes.c_char_p,      # simple_input (CHARACTER(KIND=c_char,len=*))
            ctypes.c_int,      # len_simple_input (INTEGER)
            ctypes.c_char_p,      # profileName (CHARACTER(KIND=c_char,len=*))
            ctypes.c_int,      # len_profileName (INTEGER)
            ctypes.c_bool         # nosym (LOGICAL)
        ]
        self.lib.make_inp.restype = None
        
        # make_kpt function
        # SUBROUTINE make_kpt_py(len_kpts_str, kpts_str, len_kpts_path, 
        #                        kpts_path, nosym) BIND(C, name="make_kpt")
        self.lib.make_kpt.argtypes = [
            ctypes.c_int,      # len_kpts_str (INTEGER(C_SIZE_T))
            ctypes.c_char_p,      # kpts_str (CHARACTER(KIND=c_char,len=*))
            ctypes.c_int,      # len_kpts_path (INTEGER(C_SIZE_T))
            ctypes.c_char_p,      # kpts_path (CHARACTER(KIND=c_char,len=*))
            ctypes.c_bool         # nosym (LOGICAL)
        ]
        self.lib.make_kpt.restype = None
    
    def say_hello(self,Num: int =2 ) -> None:
        """Test function to verify library connection."""
        self.lib.say_hello(Num)
    
    def make_inp(self, 
                 simple_input: str, 
                 profile_name: str = "default_profile",
                 nosym: bool = False) -> None:
        """
        Generate inp.xml from simple input string.
        
        Args:
            simple_input: FLEUR simple input format string
            profile_name: Profile name for inpgen defaults
            nosym: If True, disable symmetry detection
        
        Examples:
            Simple format:
            >>> inpgen = InpgenInterface()
            >>> simple_input = "Cu fcc\\n3.61\\n"
            >>> inpgen.make_inp(simple_input)
            
            Namelist format:
            >>> simple_input = '''&input film=f /
            ... &lattice latsys='bcc' a=2.87 /
            ... 1
            ... 26  0.0 0.0 0.0
            ... &atom element="Fe" /
            ... &comp /
            ... &kpt div1=8 div2=8 div3=8 /
            ... '''
            >>> inpgen.make_inp(simple_input)
        """
        # Encode strings to bytes
        c_simple_input = ctypes.c_char_p(simple_input.encode('utf-8'))
        c_profile_name = ctypes.c_char_p(profile_name.encode('utf-8'))

        # Get lengths (excluding null terminator)
        len_simple = ctypes.c_int(len(simple_input))
        len_profile = ctypes.c_int(len(profile_name))
        
        self._log(f"Calling make_inp with {len_simple.value} chars input, profile='{profile_name}'")
        
        try:
            # Call Fortran function with lengths
            self.lib.make_inp(
                c_simple_input,
                len_simple,
                c_profile_name,
                len_profile,
                ctypes.c_bool(nosym)
            )
            self._log("✓ inp.xml generated successfully")
        except Exception as e:
            self._log(f"✗ Error generating inp.xml: {e}")
            raise
    
    def make_inp_from_file(self,
                          input_file: Union[str, Path],
                          profile_name: str = "default",
                          nosym: bool = False) -> None:
        """
        Generate inp.xml from simple input file.
        
        Args:
            input_file: Path to simple input file
            profile_name: Profile name for inpgen defaults
            nosym: If True, disable symmetry detection
        
        Example:
            >>> inpgen = InpgenInterface()
            >>> inpgen.make_inp_from_file("input.txt")
        """
        input_path = Path(input_file)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        with open(input_path, 'r') as f:
            simple_input = f.read()
        
        self._log(f"Reading input from: {input_path}")
        self.make_inp(simple_input, profile_name, nosym)
    
    def add_kpoints(self,
                    kpts_str: str,
                    kpts_path: str = "",
                    nosym: bool = False) -> None:
        """
        Add k-points to existing inp.xml.
        
        Args:
            kpts_str: K-point specification string (e.g., "8 8 8" or "band=100")
            kpts_path: Path specification for band structure (e.g., "GXWLG", "default")
            nosym: If True, disable symmetry in k-points
        
        Examples:
            >>> inpgen = InpgenInterface()
            >>> inpgen.add_kpoints("8 8 8")
            >>> inpgen.add_kpoints("band=100", "GXWLG")
            >>> inpgen.add_kpoints("band=240", "default")
        """
        # Encode strings to bytes
        c_kpts_str = kpts_str.encode('utf-8')
        c_kpts_path = kpts_path.encode('utf-8')
        
        # Get lengths (excluding null terminator)
        len_kpts_str = len(c_kpts_str)
        len_kpts_path = len(c_kpts_path)
        
        self._log(f"Adding k-points: '{kpts_str}' path='{kpts_path}'")
        
        try:
            # Call Fortran function with lengths
            self.lib.make_kpt(
                ctypes.c_int(len_kpts_str),
                ctypes.c_char_p(c_kpts_str),
                ctypes.c_int(len_kpts_path),
                ctypes.c_char_p(c_kpts_path),
                ctypes.c_bool(nosym)
            )
            self._log(f"✓ K-points added successfully")
        except Exception as e:
            self._log(f"✗ Error adding k-points: {e}")
            raise


# Convenience functions
def make_inp(simple_input: str,
             profile_name: str = "default",
             nosym: bool = False,
             lib_path: Optional[str] = None) -> None:
    """
    Convenience function to generate inp.xml.
    
    Args:
        simple_input: FLEUR simple input format string
        profile_name: Profile name for inpgen defaults
        nosym: If True, disable symmetry detection
        lib_path: Optional path to shared library
    
    Example:
        >>> from FleurInpgen import make_inp
        >>> make_inp("Cu fcc\\n3.61\\n")
    """
    inpgen = InpgenInterface(lib_path)
    inpgen.make_inp(simple_input, profile_name, nosym)


def make_inp_from_file(input_file: Union[str, Path],
                       profile_name: str = "default_profile",
                       nosym: bool = False,
                       lib_path: Optional[str] = None) -> None:
    """
    Convenience function to generate inp.xml from file.
    
    Args:
        input_file: Path to simple input file
        profile_name: Profile name for inpgen defaults
        nosym: If True, disable symmetry detection
        lib_path: Optional path to shared library
    
    Example:
        >>> from FleurInpgen import make_inp_from_file
        >>> make_inp_from_file("input.txt")
    """
    inpgen = InpgenInterface(lib_path)
    inpgen.make_inp_from_file(input_file, profile_name, nosym)


def add_kpoints(kpts_str: str,
                kpts_path: str = "",
                nosym: bool = False,
                lib_path: Optional[str] = None) -> None:
    """
    Convenience function to add k-points to inp.xml.
    
    Args:
        kpts_str: K-point specification string
        kpts_path: Path specification for band structure
        nosym: If True, disable symmetry in k-points
        lib_path: Optional path to shared library
    
    Example:
        >>> from FleurInpgen import add_kpoints
        >>> add_kpoints("8 8 8")
        >>> add_kpoints("band=100", "GXWLG")
    """
    inpgen = InpgenInterface(lib_path)
    inpgen.add_kpoints(kpts_str, kpts_path, nosym)


# Example usage
if __name__ == "__main__":
    import sys
    
    print("FLEUR Inpgen Python Interface")
    print("=" * 60)
    
    try:
        # Initialize interface
        inpgen = InpgenInterface()
        
        # Test library connection
        print("\n1. Testing library connection...")
        inpgen.say_hello()
        print("✓ Library connection successful\n")
        
        
        # Example 4: Namelist format (Fe bcc) 
        print("2. Example: Fe bcc (namelist format)")
        namelist_input = """&input film=f /
&lattice latsys='fcc' a=2.87 /
1
26  0.0 0.0 0.0
&atom element="Fe" jri=565 lmax=8 lnonsph=6 /
&comp kmax=3.4 gmaxxc=10.2 gmax=8.5 /
&kpt div1=8 div2=8 div3=8 tkb=0.001 /
"""
        inpgen.make_inp(namelist_input)
        print()
        
        print("=" * 60)
        print("✓ All examples completed successfully!")
        
    except FileNotFoundError as e:
        print(f"\n✗ Error: {e}", file=sys.stderr)
        print("\nPlease build the inpgen library first:", file=sys.stderr)
        print("  cd build", file=sys.stderr)
        print("  cmake ..", file=sys.stderr)
        print("  make inpgen3", file=sys.stderr)
        sys.exit(1)
        
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)