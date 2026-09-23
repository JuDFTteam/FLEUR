"""Layering invariants of the matrix-element / wannierisation split.

These do not run FLEUR: they read the source and check which module uses which.

The one that matters is that matrixelements/ does not depend on wannierlib/. That is what
makes it a layer other code can use -- secvar_soc already does -- rather than a private part
of the wannierisation. It holds today by accident of how the code was written; this makes it
hold on purpose, so that a USE added in the wrong direction fails here instead of quietly
turning the two directories into one.

The second is that postproc/ does not depend on the driver above it. It takes a gauge and
writes files, and it says so by importing nothing from wannierlib/ itself -- which keeps the
stack three layers deep rather than two directories that happen to sit apart.
"""
import os
import re

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "src", "fleur"))

ZONES = {
    "matrixelements": ["matrixelements"],
    "postproc": [os.path.join("wannierlib", "postproc")],
    "wannierlib": ["wannierlib"],
}

MODULE_RE = re.compile(r"^\s*MODULE\s+([A-Za-z_]\w*)\s*$", re.IGNORECASE)
USE_RE = re.compile(r"^\s*USE\s+([A-Za-z_]\w*)", re.IGNORECASE)


def _sources(zone):
    """(path, filename) of every Fortran source in a zone, without recursing into others."""
    for rel in ZONES[zone]:
        d = os.path.join(SRC, rel)
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn.endswith((".F90", ".f90")):
                yield os.path.join(d, fn), fn


def _module_owner():
    """module name (lower case) -> (zone, filename) for every module in the two zones."""
    owner = {}
    for zone in ZONES:
        for path, fn in _sources(zone):
            with open(path, errors="ignore") as fh:
                for line in fh:
                    m = MODULE_RE.match(line)
                    if m:
                        owner[m.group(1).lower()] = (zone, fn)
    return owner


def _imports_from(zone, forbidden):
    """Every USE in `zone` of a module owned by one of `forbidden`, as printable lines."""
    owner = _module_owner()
    assert owner, f"no Fortran modules found under {SRC} -- the paths in this test are stale"

    offenders = []
    for path, fn in _sources(zone):
        with open(path, errors="ignore") as fh:
            for lineno, line in enumerate(fh, 1):
                m = USE_RE.match(line)
                if not m:
                    continue
                used = m.group(1).lower()
                where = owner.get(used, ("", ""))
                if where[0] in forbidden:
                    offenders.append(
                        f"  {zone}/{fn}:{lineno} USE {m.group(1)}"
                        f"  (lives in {where[0]}/{where[1]})")
    return offenders


def test_matrixelements_does_not_use_wannierlib():
    """No module under matrixelements/ may USE one from wannierlib/ or its postproc/."""
    offenders = _imports_from("matrixelements", {"wannierlib", "postproc"})

    assert not offenders, (
        "matrixelements/ must not depend on wannierlib/: it is a layer other code uses "
        "(secvar_soc does), not a private part of the wannierisation. Offending imports:\n"
        + "\n".join(offenders)
        + "\n\nMove what is shared to a place both can see, or keep the wannier-specific "
          "part on the wannierlib side.")


def test_postproc_does_not_use_the_driver():
    """No module under wannierlib/postproc/ may USE one from wannierlib/ itself."""
    offenders = _imports_from("postproc", {"wannierlib"})

    assert not offenders, (
        "wannierlib/postproc/ must not depend on the driver above it: it takes the gauge "
        "and writes the files, and everything it needs is passed in. Offending imports:\n"
        + "\n".join(offenders)
        + "\n\nPass what is missing as an argument, or move the shared part down into "
          "matrixelements/.")
