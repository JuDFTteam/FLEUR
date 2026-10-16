#!/usr/bin/env python3
"""Check new or changed www.flapw.de pages before committing them.

  check_pages.py WWW_ROOT [PAGE ...]

Without PAGE arguments, all Markdown files under docs/ that differ from origin/develop
(including untracked ones) are checked. Errors (exit code 1):
  - page missing from the mkdocs.yml nav
  - broken relative link or image
  - image name used more than once under docs/documentation/ (the PDF build copies all
    images into one directory)
  - unbalanced ``` fences, $$ or \\begin{equation}
  - admonition (!!!) without indented body
Warnings:
  - new page under docs/documentation/ outside Theory/, Calculations/ or Reference/
  - page under docs/documentation/ missing from latex/include.list (not in the PDF manual)
  - inline HTML, unknown #anchor, no heading at the top
  - xml snippet that is not well formed, or uses elements/attributes unknown to the
    current FLEUR schemas; inline XPath such as `output/bandDOS/@sigma` that the schemas lack
"""
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import schema_diff  # noqa: E402  (also checks that this Python can parse XML)

ET = schema_diff.ET

LINK = re.compile(r"!?\[[^\]]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
FENCE = re.compile(r"^\s*```(?!.*```)\s*(\S*)")  # opening/closing fence, not a one-line ```span```
HTML = re.compile(r"<(div|span|img|br|p|table|tr|td|font|center|b|i|u|sup|sub|details)\b", re.I)
SECTIONS = ("documentation/Theory/", "documentation/Calculations/", "documentation/Reference/")
XPATH = re.compile(r"`{1,3}(/?[A-Za-z_][\w]*(?:/[A-Za-z_][\w]*)*(?:/@[A-Za-z_]\w*)?)`{1,3}")


def changed_pages(www):
    def git(*a):
        out = subprocess.run(["git", "-C", www] + list(a), capture_output=True, text=True)
        return out.stdout.split() if out.returncode == 0 else []

    files = set(git("diff", "--name-only", "origin/develop", "--", "docs"))
    files |= set(git("ls-files", "--others", "--exclude-standard", "--", "docs"))
    return sorted(f for f in files if f.endswith(".md") and os.path.exists(os.path.join(www, f)))


def slug(heading):
    s = re.sub(r"[^\w\s-]", "", heading.strip().lower())
    return re.sub(r"[\s]+", "-", s).strip("-")


def anchors(path):
    with open(path, encoding="utf-8") as f:
        return {slug(m.group(1)) for m in re.finditer(r"^#+\s+(.*?)\s*#*$", f.read(), re.M)}


def schema_names():
    """element name -> set of attribute names, and all XPaths, from the current schemas."""
    attrs, paths = {}, []
    for kind in ("input", "output"):
        for path, e in schema_diff.load(kind).items():
            name = path.rsplit("/", 1)[1]
            attrs.setdefault(name, set()).update(e["attrs"])
            paths.append((path, set(e["attrs"])))
    return attrs, paths


def strip_code(text):
    """Text with fenced blocks and inline code removed, plus the list of fenced blocks."""
    blocks, out, cur, lang = [], [], None, None
    for line in text.splitlines():
        m = FENCE.match(line)
        if m and cur is None:
            cur, lang = [], m.group(1).lower()
            out.append("")
            continue
        if m and cur is not None:
            blocks.append((lang, "\n".join(cur)))
            cur = None
            out.append("")
            continue
        if cur is not None:
            cur.append(line)
            out.append("")
        else:
            out.append(re.sub(r"`[^`]*`", "", line))
    return "\n".join(out), blocks


def image_names(www):
    """basename -> list of PNG/JPG/SVG paths under docs/documentation/."""
    found = {}
    for root, _, files in os.walk(os.path.join(www, "docs", "documentation")):
        for f in files:
            if f.lower().endswith((".png", ".jpg", ".jpeg", ".svg", ".gif")):
                found.setdefault(f, []).append(os.path.relpath(os.path.join(root, f), www))
    return found


def check(www, rel, nav, names, paths, report, layout):
    path = os.path.join(www, rel)
    pdir = os.path.dirname(path)
    with open(path, encoding="utf-8") as f:
        text = f.read()
    lines = text.splitlines()
    err = lambda msg: report.append(("E", rel, msg))  # noqa: E731
    warn = lambda msg: report.append(("W", rel, msg))  # noqa: E731

    site = rel[len("docs/"):] if rel.startswith("docs/") else rel
    if rel.startswith("docs/") and site not in nav:
        err("not in the mkdocs.yml nav (add '- <Title>: %s')" % site)
    if site.startswith("documentation/"):
        if not site.startswith(SECTIONS) and site not in layout["existing"]:
            warn("feature pages belong in documentation/Theory/, Calculations/ or Reference/")
        if layout["include"] is not None and site not in layout["include"]:
            warn("not in latex/include.list - the page will be missing from the PDF manual")
    if sum(1 for line in lines if FENCE.match(line)) % 2:
        err("unbalanced ``` fences")

    prose, blocks = strip_code(text)
    if prose.count("$$") % 2:
        err("odd number of $$")
    if prose.count("\\begin{equation}") != prose.count("\\end{equation}"):
        err("\\begin{equation} and \\end{equation} do not match")
    first = next((line for line in lines if line.strip()), "")
    if not first.startswith("#"):
        warn("page does not start with a heading")

    plines = prose.splitlines()
    html_lines = []
    for i, line in enumerate(plines):
        if line.startswith("!!!"):
            body = next((x for x in plines[i + 1:] if x.strip()), "")
            if not body.startswith("    "):
                err("line %d: admonition without a 4-space indented body" % (i + 1))
        if HTML.search(line):
            html_lines.append(str(i + 1))
        for target in LINK.findall(line):
            if re.match(r"^[a-z]+:", target) or target.startswith("#"):
                if target.startswith("#") and target[1:] not in anchors(path):
                    warn("line %d: no heading for anchor %s" % (i + 1, target))
                continue
            file, _, anchor = target.partition("#")
            tpath = os.path.normpath(os.path.join(pdir, file))
            if not os.path.exists(tpath):
                err("line %d: broken link %s" % (i + 1, target))
            elif len(layout["images"].get(os.path.basename(file), [])) > 1:
                err("line %d: image name %s is not unique under docs/documentation/ (%s)"
                    % (i + 1, os.path.basename(file), ", ".join(layout["images"][os.path.basename(file)])))
            elif anchor and tpath.endswith(".md") and anchor not in anchors(tpath):
                warn("line %d: no heading for anchor %s in %s" % (i + 1, anchor, file))
    if html_lines:
        more = ", ..." if len(html_lines) > 8 else ""
        warn("inline HTML on %d line(s) (%s%s) - use Markdown" % (len(html_lines), ", ".join(html_lines[:8]), more))

    if names is None:
        return
    for lang, code in blocks:
        if lang != "xml":
            continue
        try:
            tree = ET.fromstring("<snippet>%s</snippet>" % re.sub(r"<\?xml[^>]*\?>", "", code))
        except ET.ParseError as e:
            warn("xml snippet not well formed (%s): %s" % (e, code.strip().splitlines()[0][:60]))
            continue
        for el in tree.iter():
            if el is tree:
                continue
            if el.tag not in names:
                warn("xml snippet: element <%s> is not in the FLEUR schemas" % el.tag)
                continue
            for a in el.attrib:
                if a not in names[el.tag]:
                    warn("xml snippet: <%s %s=...> is not in the FLEUR schemas" % (el.tag, a))
    for m in XPATH.finditer(text):
        xp = m.group(1)
        if "/" not in xp:
            continue
        elem, _, attr = xp.partition("/@")
        suffix = "/" + elem.strip("/")
        hits = [a for p, a in paths if p.endswith(suffix)]
        if not hits:
            warn("XPath %s: no such element in the FLEUR schemas" % xp)
        elif attr and not any(attr in a for a in hits):
            warn("XPath %s: element exists but has no attribute %s" % (xp, attr))


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    www = os.path.abspath(sys.argv[1])
    mk = os.path.join(www, "mkdocs.yml")
    if not os.path.exists(mk):
        sys.exit("%s has no mkdocs.yml" % www)
    with open(mk, encoding="utf-8") as f:
        nav = set(re.findall(r":\s*(\S+\.md)\s*$", f.read(), re.M))
    pages = [p if os.path.exists(os.path.join(www, p)) else os.path.relpath(os.path.abspath(p), www)
             for p in sys.argv[2:]] or changed_pages(www)
    if not pages:
        sys.exit("no changed pages found (compare against origin/develop) - pass them explicitly")
    try:
        names, paths = schema_names()
    except SystemExit:
        names, paths = None, None
        print("W  FLEUR schemas not found - xml/XPath checks skipped")
    inc = os.path.join(www, "latex", "include.list")
    include = None
    if os.path.exists(inc):
        with open(inc, encoding="utf-8") as f:
            include = {line.strip() for line in f if line.strip()}
    existing = set(subprocess.run(["git", "-C", www, "ls-tree", "-r", "--name-only", "origin/develop", "docs/"],
                                  capture_output=True, text=True).stdout.split())
    layout = {"include": include, "images": image_names(www),
              "existing": {p[len("docs/"):] for p in existing}}
    report = []
    for rel in pages:
        check(www, rel, nav, names, paths, report, layout)
    for level, rel, msg in report:
        print("%s  %s: %s" % (level, rel, msg))
    n_err = sum(1 for r in report if r[0] == "E")
    print("%d page(s), %d error(s), %d warning(s)" % (len(pages), n_err, len(report) - n_err))
    sys.exit(1 if n_err else 0)


if __name__ == "__main__":
    main()
