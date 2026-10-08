#!/usr/bin/env python3
"""Extract magnetic moments from FLEUR out.xml and compare with the SCAN
reference values in reference_values.csv.

Called by run_mgga_tests.sh; usable on its own:

    ./compare_moments.py --workdir work fm/Fe af/CrSb2

Quantities, following reference_values.csv:

* FM  cases: total spin moment per formula unit, i.e. the cell moment
  (interstitial + all muffin tins) divided by the number of formula units.
  The number of formula units is derived from <atomsInCell nat=".."> in
  out.xml and the stoichiometry implied by the directory name.

* AFM cases: the muffin-tin spin moment on the transition-metal atom. Note
  the paper quotes a *Bader* moment for these; the two integrate over
  different volumes and will not agree exactly. See README.md caveat 4.
"""
import argparse
import csv
import os
import re
import sys

MU = "μ"


def formula_atoms(name):
    """Atoms per formula unit from a chemical formula, e.g. Ni3Al -> 4."""
    parts = re.findall(r"([A-Z][a-z]?)(\d*)", name)
    return sum(int(n) if n else 1 for _, n in parts if _)


def last_iteration(text):
    i = text.rfind("<iteration")
    return text[i:] if i >= 0 else text


def parse_outxml(path):
    s = open(path, errors="replace").read()
    li = last_iteration(s)
    out = {}

    m = re.search(r'<atomsInCell nat="(\d+)" ntype="(\d+)"', s)
    if m:
        out["nat"], out["ntype"] = int(m.group(1)), int(m.group(2))

    # element per atom type, in declaration order
    out["elements"] = re.findall(r'<species name="[^"]*" element="(\w+)"', s)

    # cell-integrated charge per spin -> total moment
    ch = {int(a): float(b) for a, b in
          re.findall(r'<spinDependentCharge spin="(\d)" total="([-\d.Ee+]+)"', li)}
    out["cell_moment"] = ch[1] - ch[2] if 1 in ch and 2 in ch else None

    # muffin-tin spin moments; orbitalMomentsInMTSpheres reuses localMagMoment,
    # so restrict to the magnetic block
    blk = re.search(r"<magneticMomentsInMTSpheres.*?</magneticMomentsInMTSpheres>",
                    li, re.S)
    mts = {}
    if blk:
        for a, _x, _y, z in re.findall(
                r'<localMagMoment atomType="(\d+)" vec="\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)"',
                blk.group(0)):
            mts[int(a)] = float(z)
    out["mt_moments"] = mts

    d = re.findall(r'<chargeDensity spin="\d" distance="([\d.Ee+-]+)"', li)
    out["distance"] = max(float(x) for x in d) if d else None
    out["iterations"] = len(re.findall(r"<iteration", s))
    return out


def read_reference(suite):
    ref = {}
    with open(os.path.join(suite, "reference_values.csv")) as fh:
        for row in csv.DictReader(fh):
            ref[row["material"]] = row
    return ref


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--suite", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("cases", nargs="*")
    args = ap.parse_args()

    ref = read_reference(args.suite)
    cases = args.cases or sorted(
        os.path.join(c, m)
        for c in ("fm", "af")
        for m in sorted(os.listdir(os.path.join(args.suite, c)))
        if os.path.isdir(os.path.join(args.suite, c, m)))

    hdr = ("%-12s %-4s %-9s %9s %9s %9s  %s"
           % ("case", "type", "converged", "FLEUR", "SCAN ref", "diff", "note"))
    print(hdr)
    print("-" * len(hdr))

    missing = 0
    for rel in cases:
        material = os.path.basename(rel)
        r = ref.get(material)
        out_xml = os.path.join(args.workdir, rel, "out.xml")

        if not os.path.exists(out_xml):
            print("%-12s %-4s %-9s %9s %9s %9s  %s"
                  % (rel, r["category"] if r else "?", "-", "-",
                     r["scan_reference_muB"] if r else "-", "-", "no out.xml"))
            missing += 1
            continue

        d = parse_outxml(out_xml)
        target = float(r["scan_reference_muB"]) if r else None
        note = ""

        if r and r["category"] == "FM":
            apf = formula_atoms(material)
            nfu = (d.get("nat") or apf) / apf
            value = d["cell_moment"] / nfu if d["cell_moment"] is not None else None
            note = "%d atom(s)/cell = %g f.u." % (d.get("nat", 0), nfu)
        else:
            # AFM: muffin-tin moment on the transition-metal atom
            tm = r["tm_atom"] if r else None
            idx = [i + 1 for i, e in enumerate(d["elements"]) if e == tm]
            vals = [d["mt_moments"][i] for i in idx if i in d["mt_moments"]]
            value = max((abs(v) for v in vals), default=None)
            note = "MT moment on %s (paper quotes Bader)" % tm

        if d["distance"] is None:
            conv = "?"
        elif d["distance"] < 1e-5:
            conv = "yes"
        else:
            conv = "NO"
        diff = (value - target) if (value is not None and target is not None) else None

        print("%-12s %-4s %-9s %9s %9s %9s  %s"
              % (rel,
                 r["category"] if r else "?",
                 "%s" % conv,
                 "%.3f" % value if value is not None else "-",
                 "%.2f" % target if target is not None else "-",
                 "%+.3f" % diff if diff is not None else "-",
                 note))
        if conv == "NO":
            print("%-12s   ^ density distance %.2e after %d iterations - not converged, "
                  "moment not meaningful" % ("", d["distance"], d["iterations"]))
        elif conv == "?":
            print("%-12s   ^ no density distance in out.xml - run did not complete an "
                  "iteration" % "")

    print()
    print("Reference: Tran et al., Phys. Rev. B 102, 024407 (2020), SCAN column.")
    print("Differences are physics, not regressions - read README.md before")
    print("treating a mismatch as a FLEUR bug (k-convergence in particular).")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
