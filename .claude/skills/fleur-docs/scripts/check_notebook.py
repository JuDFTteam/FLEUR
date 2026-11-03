#!/usr/bin/env python3
"""Check a fleur_tutorial_container notebook, and optionally run it against a local FLEUR build.

  check_notebook.py NOTEBOOK [--www WWW_ROOT]                       static checks
  check_notebook.py NOTEBOOK --execute --build-dir ../fleur/build   also execute a copy
  check_notebook.py NOTEBOOK --clear                                 clear outputs in place

Static checks: valid nbformat 4 JSON, bash/sos kernel, '# Title' first, sections
Objectives/Introduction/Learn more, cleared outputs, no TODO left, balanced ``` in markdown,
local links/images and 'cd' targets exist, registered in the chapter Outline.ipynb, and
(with --www) links to www.flapw.de/rel/documentation/<Section>/<page>/ have a page.

--execute copies the notebook's folder to a scratch directory and runs
'jupyter nbconvert --execute' there, or 'sos convert --execute' (sos-papermill) for a SoS
notebook, which runs each cell in its own sub-kernel (Bash, Python3). A PATH shim maps
fleur_MPI, fleur and inpgen to the build directory, the stack limit is raised, and
OMP_NUM_THREADS=2 as in the container. --python PY uses that interpreter (a venv with
nbconvert + bash_kernel, plus the sos packages for SoS notebooks) and only its kernels
instead of the 'jupyter' on PATH. Exit code 2 means execution was skipped (no Jupyter,
no kernel or no build). The original notebook is never modified, except by --clear.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)\)")
DOCLINK = re.compile(r"flapw\.de/[\w.-]+/documentation/([\w./-]+?)/?(?=[#)\s>\]]|$)")
CD = re.compile(r"(?:^|[;&|]\s*)cd\s+([^\s;&|]+)")
MKDIR = re.compile(r"(?:^|[;&|]\s*)mkdir\s+(?:-p\s+)?([^\s;&|]+)")
INPGEN_F = re.compile(r"inpgen\s+(?:-\S+\s+)*-f\s+([^\s;&|]+)")
BARE_FLEUR =re.compile(r"(?:^|[\s;&|])fleur(?:\s|$)")
FENCE = re.compile(r"^\s*```(?!.*```)")  # not a one-line ```span```
RUNTIME_ERRORS = re.compile(r"juDFT-Error|command not found|No such file or directory|Segmentation fault|Traceback|Abort")


def src(c):
    return "".join(c.get("source", []))


def static_checks(path, nb, www, report):
    err = lambda m: report.append(("E", m))  # noqa: E731
    warn = lambda m: report.append(("W", m))  # noqa: E731
    nb_dir = os.path.dirname(os.path.abspath(path))
    cells = nb.get("cells", [])
    if nb.get("nbformat") != 4:
        err("nbformat is %s, expected 4" % nb.get("nbformat"))
    kernel = nb.get("metadata", {}).get("kernelspec", {}).get("name")
    if kernel not in ("bash", "sos"):
        warn("kernel '%s' - FLEUR tutorials use bash (or sos for Python plots)" % kernel)
    if not cells or cells[0]["cell_type"] != "markdown" or not src(cells[0]).startswith("# "):
        err("first cell must be markdown starting with '# <Title>'")
    md = "\n".join(src(c) for c in cells if c["cell_type"] == "markdown")
    for section in ("## Objectives", "## Introduction", "## Learn more"):
        if not re.search("^" + re.escape(section) + r"\s*$", md, re.M):
            warn("missing section '%s'" % section)

    cwd = nb_dir
    made = set()
    for i, c in enumerate(cells):
        text = src(c)
        if "TODO" in text:
            warn("cell %d: TODO left" % i)
        if c["cell_type"] == "markdown":
            if sum(1 for line in text.splitlines() if FENCE.match(line)) % 2:
                err("cell %d: unbalanced ``` fence" % i)
            for target in LINK.findall(text):
                if re.match(r"^[a-z]+:", target) or target.startswith("#"):
                    continue
                if not os.path.exists(os.path.join(nb_dir, target.split("#")[0])):
                    err("cell %d: broken link %s" % (i, target))
            for page in DOCLINK.findall(text):
                if "/" not in page:
                    warn("cell %d: link .../documentation/%s/ uses the old flat layout - use "
                         ".../documentation/<Calculations|Reference|Theory>/%s/" % (i, page, page))
                elif www and not os.path.exists(os.path.join(www, "docs", "documentation", page + ".md")):
                    err("cell %d: no page docs/documentation/%s.md in %s" % (i, page, www))
            continue
        if c.get("outputs") or c.get("execution_count") is not None:
            err("cell %d: outputs not cleared (run with --clear)" % i)
        for line in text.splitlines():
            if line.lstrip().startswith("#"):
                continue
            for d in MKDIR.findall(line):
                made.add(os.path.normpath(os.path.join(cwd, d)))
            if cwd is None:
                break  # lost track of the working directory ($VAR, cd -, ~)
            for d in CD.findall(line):
                if re.search(r"[$~`]", d) or d == "-":
                    cwd = None
                    break
                target = os.path.normpath(os.path.join(cwd, d))
                if not (os.path.isdir(target) or target in made):
                    warn("cell %d: 'cd %s' - no such directory next to the notebook" % (i, d))
                cwd = target
            if cwd is None:
                break
            for f in INPGEN_F.findall(line):
                if not os.path.exists(os.path.join(cwd, f)) and cwd not in made:
                    warn("cell %d: 'inpgen -f %s' - file not found in %s" % (i, f, os.path.relpath(cwd, nb_dir)))
            if BARE_FLEUR.search(line):
                warn("cell %d: calls 'fleur' - the container provides 'fleur_MPI'" % i)

    outline = os.path.join(os.path.dirname(nb_dir), "Outline.ipynb")
    entry = "](%s/%s)" % (os.path.basename(nb_dir), os.path.basename(path))
    if not os.path.exists(outline):
        warn("no chapter Outline.ipynb at %s" % outline)
    else:
        with open(outline, encoding="utf-8") as f:
            ol = "\n".join(src(c) for c in json.load(f)["cells"])
        if not re.search(r"^\*\*\*[^*]+\*\*\*: *\[[^\]]*\]" + re.escape(entry[1:]), ol, re.M):
            err("not registered in %s (line '***%s***: [Title]%s')" % (outline, os.path.basename(nb_dir), entry[1:]))


def clear(path, nb):
    for c in nb.get("cells", []):
        if c["cell_type"] == "code":
            c["outputs"], c["execution_count"] = [], None
    with open(path, "w", encoding="utf-8") as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)
        f.write("\n")
    print("cleared outputs in %s" % path)


def jupyter_env(python):
    """With a venv interpreter, its own kernelspecs take precedence over user-level ones."""
    env = dict(os.environ)
    if python:
        bindir = os.path.dirname(os.path.abspath(python))
        prefix = os.path.dirname(bindir)
        env["JUPYTER_PATH"] = os.pathsep.join(filter(None, [os.path.join(prefix, "share", "jupyter"), env.get("JUPYTER_PATH")]))
        # kernels started by name (SoS sub-kernels) must find the venv's python first
        env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
    return env


def sos_kernels(nb):
    """Jupyter kernels of the SoS sub-kernels used by the code cells."""
    names = {k[0]: k[1] for k in nb.get("metadata", {}).get("sos", {}).get("kernels", []) if len(k) > 1}
    used = {c.get("metadata", {}).get("kernel") for c in nb["cells"] if c["cell_type"] == "code"}
    return sorted(names.get(u, u) for u in used if u and u != "SoS")


def sos_cmd(python):
    """Command for 'sos convert', or None with a message if sos-papermill is missing."""
    if python:
        probe = subprocess.run([python, "-c", "import sos, sos_notebook, sos_papermill"], capture_output=True, text=True)
        if probe.returncode != 0:
            print("skip execution: SoS notebook needs sos, sos-notebook, sos-papermill, sos-bash and sos-python in the venv")
            return None
        return [python, "-m", "sos"]
    if not shutil.which("sos"):
        print("skip execution: SoS notebook needs 'sos' on PATH (pip install sos sos-notebook sos-papermill sos-bash sos-python)")
        return None
    return ["sos"]


def papermill_cell(c):
    """Error banners that sos-papermill inserts into the executed notebook."""
    return c["cell_type"] == "markdown" and "papermill" in src(c) and src(c).lstrip().startswith("<span")


def kernels(python):
    """Installed kernel names, or None with a message if Jupyter does not work."""
    if python:
        probe = [python, "-c", "import json, jupyter_client.kernelspec as k; print(json.dumps(sorted(k.find_kernel_specs())))"]
    else:
        probe = ["jupyter", "kernelspec", "list", "--json"]
    try:
        out = subprocess.run(probe, capture_output=True, text=True, env=jupyter_env(python))
    except OSError as e:
        print("skip execution: %s" % e)
        return None
    if out.returncode != 0:
        print("skip execution: listing Jupyter kernels failed:\n%s" % out.stderr.strip()[-400:])
        return None
    data = json.loads(out.stdout)
    return data if python else list(data.get("kernelspecs", {}))


def nbconvert_cmd(python):
    return [python, "-m", "nbconvert"] if python else ["jupyter", "nbconvert"]


def execute(path, nb, args, report):
    """Returns False if execution had to be skipped."""
    kernel = nb.get("metadata", {}).get("kernelspec", {}).get("name", "bash")
    installed = kernels(args.python)
    if installed is None:
        return False
    if kernel not in installed:
        print("skip execution: Jupyter kernel '%s' is not installed (python -m bash_kernel.install --sys-prefix)" % kernel)
        return False
    if kernel == "sos":
        missing = [k for k in sos_kernels(nb) if k not in installed]
        if missing:
            print("skip execution: SoS sub-kernel(s) %s not installed" % ", ".join(missing))
            return False
        runner = sos_cmd(args.python)
        if runner is None:
            return False
    build =os.path.abspath(args.build_dir or os.environ.get("FLEUR_BUILD_DIR", ""))
    fleur = next((os.path.join(build, b) for b in ("fleur_MPI", "fleur") if os.path.isfile(os.path.join(build, b))), None)
    inpgen = os.path.join(build, "inpgen")
    if not fleur or not os.path.isfile(inpgen):
        print("skip execution: no fleur/fleur_MPI and inpgen in build dir '%s'" % build)
        return False
    code = "\n".join(src(c) for c in nb["cells"] if c["cell_type"] == "code")
    for tool in ("mpirun", "aiida-fleur", "verdi", "masci-tools"):
        if re.search(r"(?:^|[\s;&|])%s\b" % re.escape(tool), code) and not shutil.which(tool):
            report.append(("W", "notebook uses '%s', which is not on PATH here - those cells will fail locally" % tool))

    work = tempfile.mkdtemp(prefix="nbrun-", dir=args.workdir)
    shim = os.path.join(work, "bin")
    os.makedirs(shim)
    for name, target in (("fleur_MPI", fleur), ("fleur", fleur), ("inpgen", inpgen)):
        os.symlink(target, os.path.join(shim, name))
    nb_dir = os.path.dirname(os.path.abspath(path))
    copy = os.path.join(work, os.path.basename(nb_dir))
    shutil.copytree(nb_dir, copy)
    env = jupyter_env(args.python)
    os.makedirs(os.path.join(work, "jupyter-config"))  # ignore personal nbconvert config
    env.update(PATH=shim + os.pathsep + env.get("PATH", ""), OMP_NUM_THREADS="2",
               JUPYTER_CONFIG_DIR=os.path.join(work, "jupyter-config"))
    if kernel == "sos":
        # sos-papermill continues after failing cells and exits 1; it has no per-cell timeout
        tool = runner + ["convert", os.path.basename(path), "executed.ipynb", "--execute"]
    else:
        tool = nbconvert_cmd(args.python) + [
            "--to", "notebook", "--execute", "--allow-errors",
            "--ExecutePreprocessor.timeout=%d" % args.timeout,
            "--ExecutePreprocessor.kernel_name=%s" % kernel,
            "--output", "executed.ipynb", os.path.basename(path)]
    cmd = ["bash", "-c", 'ulimit -s unlimited 2>/dev/null || ulimit -s 65520 2>/dev/null; exec "$@"', "_"] + tool
    print("executing in %s (fleur: %s) ..." % (copy, fleur))
    run = subprocess.run(cmd, cwd=copy, env=env, capture_output=True, text=True)
    executed = os.path.join(copy, "executed.ipynb")
    if not os.path.exists(executed) or (run.returncode != 0 and kernel != "sos"):
        report.append(("E", "%s failed (exit %d):\n%s" % ("sos convert" if kernel == "sos" else "nbconvert", run.returncode,
                                                          (run.stderr or run.stdout).strip()[-1500:])))
        return True
    with open(executed, encoding="utf-8") as f:
        result = json.load(f)
    cells = [c for c in result["cells"] if not papermill_cell(c)]
    for i, c in enumerate(cells):
        outs = c.get("outputs", [])
        text = "".join("".join(o.get("text", [])) for o in outs if o.get("output_type") == "stream")
        for line in text.splitlines():
            if RUNTIME_ERRORS.search(line):
                report.append(("E", "cell %d output: %s" % (i, line.strip()[:200])))
        for o in outs:
            if o.get("output_type") == "error":
                lines = [line for line in text.strip().splitlines() if not line.startswith("[IPKernelApp]")]
                tail = " | ".join(line.strip() for line in lines[-3:])
                what = (o.get("ename") + ": ") if o.get("ename") else "exit status "
                report.append(("E", "cell %d failed (%s%s): %s" % (i, what, o.get("evalue", "")[:100], tail[:300])))
    print("executed notebook: %s" % executed)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("notebook")
    p.add_argument("--www", help="www.flapw.de checkout/worktree to resolve Learn-more links")
    p.add_argument("--clear", action="store_true", help="clear all outputs in place, then check")
    p.add_argument("--execute", action="store_true")
    p.add_argument("--build-dir", help="FLEUR build directory (default: $FLEUR_BUILD_DIR)")
    p.add_argument("--python", help="interpreter of a venv with nbconvert and bash_kernel (SoS notebooks: also the sos packages)")
    p.add_argument("--timeout", type=int, default=1800, help="seconds per cell (default 1800; not applied to SoS notebooks)")
    p.add_argument("--workdir", help="parent of the scratch run directory (default: system temp)")
    args = p.parse_args()

    with open(args.notebook, encoding="utf-8") as f:
        try:
            nb = json.load(f)
        except json.JSONDecodeError as e:
            sys.exit("E  not valid JSON: %s" % e)
    if args.clear:
        clear(args.notebook, nb)
    report = []
    static_checks(args.notebook, nb, args.www and os.path.abspath(args.www), report)
    ran = execute(args.notebook, nb, args, report) if args.execute else True
    for level, msg in report:
        print("%s  %s" % (level, msg))
    n_err = sum(1 for r in report if r[0] == "E")
    print("%d error(s), %d warning(s)%s" % (n_err, len(report) - n_err, "" if ran else ", execution skipped"))
    sys.exit(1 if n_err else (0 if ran else 2))


if __name__ == "__main__":
    main()
