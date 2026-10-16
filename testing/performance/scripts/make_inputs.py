"""Generate the MnGe performance test inputs with inpgen.

Usage: python make_inputs.py <path to inpgen>

Each test is written to inputfiles/<name>/ as inp.xml together with the
inpgen input it was created from. All tests run a single SCF iteration.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

# SOC quantization axis; with it only the identity remains as symmetry
THETA = 0.1
PHI = 0.2

tests = {
    #name                    atoms  k-mesh      magnetism
    "MnGe-128-FFN-SOC-1k":   (128, (1, 1, 1), "ffn"),
    "MnGe-128-FFN-SOC-16k":  (128, (2, 2, 4), "ffn"),
    "MnGe-016-FFN-SOC-1k":   (16,  (1, 1, 1), "ffn"),
    "MnGe-016-FFN-SOC-128k": (16,  (4, 8, 4), "ffn"),
    "MnGe-128-NM-SOC-1k":    (128, (1, 1, 1), "nm"),
}


def edit_inpxml(xml, mag):
    xml = re.sub(r'itmax="\d+"', 'itmax="1"', xml)
    # LDA as noco is not recommended with GGA
    xml = re.sub(r'<xcFunctional name="[^"]*"', '<xcFunctional name="vwn"', xml)
    # drop the band structure path
    xml = re.sub(r'\s*<kPointList name="path-2".*?</kPointList>', '', xml, flags=re.S)
    if mag == "ffn":
        xml = xml.replace('l_noco="F"', 'l_noco="T"', 1)
        xml = xml.replace('l_mperp="F"', 'l_mperp="T"', 1)
        xml = xml.replace('l_mtNocoPot="F"', 'l_mtNocoPot="T"', 1)
        # in noco the SOC axis is not used, the moments carry the direction
        xml = re.sub(r'theta="[^"]*" phi="[^"]*"', 'theta="0.0" phi="0.0"', xml, count=1)
        xml = re.sub(r'<nocoParams\s+alpha="[^"]*" beta="[^"]*"/>',
                     f'<nocoParams alpha="{PHI}" beta="{THETA}"/>', xml)
    else:
        xml = xml.replace('jspins="2"', 'jspins="1"', 1)
    return xml


def make_test(inpgen, name, natoms, kmesh, mag):
    basedir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "inputfiles")
    src = os.path.join(basedir, "MnGe", f"MnGe-{natoms:04d}", f"inpMnGe-{natoms:04d}.txt")
    with open(src) as f:
        inp = f.read().rstrip() + "\n\n"
    inp += f"&soc {THETA} {PHI} /\n"
    inp += f"&kpt div1={kmesh[0]} div2={kmesh[1]} div3={kmesh[2]} /\n"

    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "inp.txt"), "w") as f:
            f.write(inp)
        args = [inpgen, "-f", "inp.txt", "-inc", "+all", "-warn_only"]
        if mag == "ffn":
            args.append("-noco")
        res = subprocess.run(args, cwd=tmp, capture_output=True, text=True)
        if res.returncode != 0 or not os.path.isfile(os.path.join(tmp, "inp.xml")):
            print(res.stdout, res.stderr)
            sys.exit(f"inpgen failed for {name}")
        with open(os.path.join(tmp, "inp.xml")) as f:
            xml = f.read()

    outdir = os.path.join(basedir, name)
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "inp.txt"), "w") as f:
        f.write(inp)
    with open(os.path.join(outdir, "inp.xml"), "w") as f:
        f.write(edit_inpxml(xml, mag))
    print("Created:", outdir)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit("Usage: python make_inputs.py <path to inpgen>")
    inpgen = os.path.abspath(sys.argv[1])
    for name, (natoms, kmesh, mag) in tests.items():
        make_test(inpgen, name, natoms, kmesh, mag)
