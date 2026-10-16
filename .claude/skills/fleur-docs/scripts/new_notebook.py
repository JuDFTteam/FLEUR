#!/usr/bin/env python3
"""Create a tutorial notebook skeleton in the fleur_tutorial_container layout.

  new_notebook.py tutorial/3.AdvancedTopics/A9/DensityMatrix.ipynb --title "Band-resolved density matrices" \
      --system Cu --input inp_Cu.txt --objective "..." --objective "..." \
      --docs Calculations/densityMatrixDOS --docs Reference/densityMatrixDOSSetup --register

Writes a bash-kernel notebook with the standard sections (Objectives, Introduction,
hands-on steps, Learn more), empty outputs and TODO markers that check_notebook.py
reports until they are replaced. --system creates the input subfolder next to the
notebook; --register adds the '***A9***: [Title](A9/Name.ipynb)' line to the chapter's
Outline.ipynb, which tutorial_to_webpage.pl in www.flapw.de parses.
"""
import argparse
import json
import os
import re
import sys
import uuid

KERNELS = {
    "bash": {
        "kernelspec": {"display_name": "Bash", "language": "bash", "name": "bash"},
        "language_info": {"codemirror_mode": "shell", "file_extension": ".sh", "mimetype": "text/x-sh", "name": "bash"},
    },
    "sos": {
        "kernelspec": {"display_name": "SoS", "language": "sos", "name": "sos"},
        "language_info": {"codemirror_mode": "sos", "file_extension": ".sos", "mimetype": "text/x-sos",
                          "name": "sos", "nbconvert_exporter": "sos_notebook.converter.SoS_Exporter",
                          "pygments_lexer": "sos"},
        "sos": {"kernels": [["Bash", "bash", "bash", "", "shell"], ["Python3", "python3", "python3", "", ""]],
                "default_kernel": "Bash"},
    },
}
LINK_LINE = re.compile(r"^\*\*\*[^*]+\*\*\*: *\[[^\]]*\]\([^)]*\)")


def lines(text):
    parts = text.split("\n")
    return [p + "\n" for p in parts[:-1]] + ([parts[-1]] if parts[-1] else [])


def cell(kind, text, sos_kernel=None):
    c = {"cell_type": kind, "id": uuid.uuid4().hex[:8], "metadata": {}, "source": lines(text)}
    if kind == "code":
        c["execution_count"] = None
        c["outputs"] = []
        if sos_kernel:
            c["metadata"]["kernel"] = sos_kernel
    return c


def skeleton(args):
    sos = "Bash" if args.kernel == "sos" else None
    objectives = "\n".join("- " + o for o in args.objective) or "- TODO: what the reader learns"
    learn = "\n".join("- [TODO: page title](https://www.flapw.de/rel/documentation/%s/)" % d for d in args.docs) \
        or "- TODO: links to the howto and reference pages on www.flapw.de"
    setup = []
    if args.system:
        setup.append("cd %s" % args.system)
    if args.input:
        setup.append("cat %s" % args.input)
    cells = [
        cell("markdown", "# %s\n\n## Objectives\n\n%s\n\n## Introduction\n\n"
             "TODO: what the feature computes and why it is useful; link the theory page if there is one." % (args.title, objectives)),
        cell("code", "\n".join(setup) or "ls", sos),
        cell("markdown", "## Setting up the calculation\n\n"
             "TODO: explain the inpgen input above, then generate inp.xml."),
        cell("code", ("inpgen -f %s" % args.input) if args.input else "inpgen -f TODO", sos),
        cell("markdown", "```{attention} Modify inp.xml\n"
             "TODO: list the changes as XPaths, e.g. set `/fleurInput/output/@dos` to `T`.\n```\n\n"
             "## Running FLEUR\n\nTODO: what this run does and how long it takes."),
        cell("code", "fleur_MPI", sos),
        cell("markdown", "## Inspecting the results\n\nTODO: which files and out.xml entries hold the result, and what to look for."),
        cell("code", "", sos),
        cell("markdown", "## Learn more\n\n" + learn),
        cell("code", "", sos),
    ]
    meta = json.loads(json.dumps(KERNELS[args.kernel]))
    return {"cells": cells, "metadata": meta, "nbformat": 4, "nbformat_minor": 5}


def write(path, nb):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
        f.write("\n")


def register(nb_path, title):
    nb_dir = os.path.dirname(os.path.abspath(nb_path))
    code, name = os.path.basename(nb_dir), os.path.basename(nb_path)
    outline = os.path.join(os.path.dirname(nb_dir), "Outline.ipynb")
    if not os.path.exists(outline):
        sys.exit("no Outline.ipynb in %s - register the tutorial by hand" % os.path.dirname(nb_dir))
    with open(outline, encoding="utf-8") as f:
        nb = json.load(f)
    entry = "***%s***: [%s](%s/%s)" % (code, title, code, name)
    for c in nb["cells"]:
        if c["cell_type"] != "markdown":
            continue
        src = "".join(c["source"]).split("\n")
        idx = [i for i, line in enumerate(src) if LINK_LINE.match(line)]
        if not idx:
            continue
        if any("](%s/%s)" % (code, name) in line for line in src):
            print("already registered in %s" % outline)
            return
        prefix = re.match(r"[A-Za-z]*", code).group(0)
        same = [i for i in idx if re.match(r"^\*\*\*%s\d" % re.escape(prefix), src[i])] if prefix else []
        after = (same or idx)[-1]
        src[after + 1:after + 1] = ["", entry]
        c["source"] = lines("\n".join(src))
        write(outline, nb)
        print("registered in %s:\n  %s" % (outline, entry))
        return
    sys.exit("no '***Xn***: [..](..)' line found in %s - register the tutorial by hand" % outline)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("notebook")
    p.add_argument("--title", required=True)
    p.add_argument("--kernel", choices=sorted(KERNELS), default="bash",
                   help="bash (default); sos when Python plotting cells are needed")
    p.add_argument("--system", help="input subfolder next to the notebook, e.g. Cu")
    p.add_argument("--input", help="inpgen input file inside --system, e.g. inp_Cu.txt")
    p.add_argument("--objective", action="append", default=[])
    p.add_argument("--docs", action="append", default=[], help="www.flapw.de page as <Section>/<page>, e.g. Calculations/dos (path under docs/documentation without .md)")
    p.add_argument("--register", action="store_true", help="add the entry to the chapter Outline.ipynb")
    args = p.parse_args()

    if os.path.exists(args.notebook):
        sys.exit("%s exists - not overwriting" % args.notebook)
    os.makedirs(os.path.dirname(os.path.abspath(args.notebook)), exist_ok=True)
    if args.system:
        os.makedirs(os.path.join(os.path.dirname(os.path.abspath(args.notebook)), args.system), exist_ok=True)
    write(args.notebook, skeleton(args))
    print("wrote %s" % args.notebook)
    if args.register:
        register(args.notebook, args.title)


if __name__ == "__main__":
    main()
