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

The third is that only the driver reaches the consumers in export/. What wannierlib/ is for is
producing U: the overlaps, the projections and the run of Wannier90. Writing the basis, the
gauge, the Bloch coefficients, the plots and the C/F tensors is what somebody does with U
afterwards, and a routine on the way to U that starts importing one of those has turned the
directory back into a bag of everything.
"""
import os
import re

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "src", "fleur"))

ZONES = {
    "matrixelements": ["matrixelements"],
    "postproc": [os.path.join("wannierlib", "postproc")],
    "export": [os.path.join("wannierlib", "export")],
    "wannierlib": ["wannierlib"],
}

# The one file in wannierlib/ that is allowed to name a consumer: it is the driver, and
# calling them in order is what it does.
DRIVER = "wannierlib_main.F90"

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


def test_matrixelements_does_not_name_the_exposure_tables():
    """The exposure tables say how the wannierisation spells things, so matrixelements/
    must not read them. They live in fleurinput/, out of reach of the zone check above,
    which is why this is asserted by name."""
    offenders = []
    for path, fn in _sources("matrixelements"):
        with open(path, errors="ignore") as fh:
            for lineno, line in enumerate(fh, 1):
                for table in ("WANNIERLIB_INTERP", "WANNIERLIB_OPR"):
                    if table in line:
                        offenders.append(f"  matrixelements/{fn}:{lineno} names {table}")

    assert not offenders, (
        "matrixelements/ must not read the wannierisation's exposure tables: what reaches "
        "it is the catalogue entry a name needs built, resolved by whoever owns the "
        "vocabulary. Offending lines:\n" + "\n".join(offenders))


def test_postproc_does_not_use_the_driver():
    """No module under wannierlib/postproc/ may USE one from wannierlib/ or its export/."""
    offenders = _imports_from("postproc", {"wannierlib", "export"})

    assert not offenders, (
        "wannierlib/postproc/ must not depend on the driver above it: it takes the gauge "
        "and writes the files, and everything it needs is passed in. Offending imports:\n"
        + "\n".join(offenders)
        + "\n\nPass what is missing as an argument, or move the shared part down into "
          "matrixelements/.")


def test_only_the_driver_uses_the_consumers():
    """Of everything in wannierlib/, only wannierlib_main.F90 may USE a module from export/.

    The rest of the directory exists to produce U, and it has to be readable as that."""
    offenders = [o for o in _imports_from("wannierlib", {"export"}) if DRIVER + ":" not in o]

    assert not offenders, (
        "only the driver may import from wannierlib/export/: the rest of the directory is "
        "the path to U, and a consumer imported halfway along it is a dependency that runs "
        "backwards. Offending imports:\n" + "\n".join(offenders)
        + "\n\nCall it from " + DRIVER + " with what it needs, or move the shared part down "
          "into postproc/.")
