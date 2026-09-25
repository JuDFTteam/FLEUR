# to be executed in the fleur/ directory
# Find Fortran source files that are not referenced in any CMakeLists.txt
# (i.e., orphaned files not compiled into the build).
# Usage: python3 cmake/scripts/find_dead_files.py

from glob import glob
import re
import os

# Source files to check (exclude vendored external/ and bundled site-packages).
files = []
for ext in ["f90", "F90", "f", "F"]:
    for f in glob(f"src/**/*.{ext}", recursive=True):
        # Skip vendored third-party and venv/site-packages trees.
        if "site-packages" not in f and "venv" not in f and "masci_tools" not in f:
            files.append(f)

cmake_files = glob("src/**/CMakeLists.txt", recursive=True)
cmake_files.extend(glob("cmake/*.txt", recursive=True))


def rm_comments(line):
    """Remove CMake comment (# to end of line)."""
    comm_pos = line.find("#")
    if comm_pos == -1:
        return line
    else:
        return line[:comm_pos]


def basename_in_cmake(cmake_files, file):
    """Check if the basename of 'file' is referenced in any CMakeLists.txt.

    Uses word-boundary regex to avoid false positives like 'kpoints.f90'
    matching inside 'make_kpoints.f90'.
    """
    base = os.path.basename(file)
    # Regex: word boundary + basename + word boundary (case-insensitive).
    pattern = r"\b" + re.escape(base) + r"\b"

    for cmf in cmake_files:
        try:
            with open(cmf, encoding="utf-8", errors="ignore") as f:
                content = f.read()
                # Remove comments from the entire file content.
                lines = content.split("\n")
                stripped = "\n".join(rm_comments(line) for line in lines)

                if re.search(pattern, stripped, re.IGNORECASE):
                    return True
        except (IOError, OSError):
            continue

    return False


def main():
    dead_files = []
    for f in files:
        if not basename_in_cmake(cmake_files, f):
            dead_files.append(f)

    if dead_files:
        print(f"Found {len(dead_files)} potentially orphaned files:")
        for f in sorted(dead_files):
            print(f)
    else:
        print("No orphaned files found.")


if __name__ == "__main__":
    main()