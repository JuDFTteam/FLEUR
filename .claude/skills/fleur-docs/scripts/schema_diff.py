#!/usr/bin/env python3
"""Elements and attributes of the FLEUR XML schemas (inp.xml, out.xml) by XPath.

  schema_diff.py diff BASE [HEAD]          what changed between two git refs
                                           (HEAD omitted: the working tree)
  schema_diff.py show XPATH [--ref REF]    one element as a Markdown reference table
  schema_diff.py find NAME [--ref REF]     XPaths of all elements/attributes called NAME

--schema input|output|both selects FleurInputSchema.xsd, FleurOutputSchema.xsd or both
(default: both). The out.xml copy of the input (/fleurOutput/fleurInput) is skipped.
Standard library only; also imported by check_pages.py.
"""
import argparse
import os
import re
import subprocess
import sys

try:
    import xml.etree.ElementTree as ET
    ET.fromstring("<a/>")
except ImportError:
    sys.exit("This Python cannot parse XML (broken pyexpat). Run the script with /usr/bin/python3 or python3.13.")

XS = "{http://www.w3.org/2001/XMLSchema}"
SCHEMA_DIR = "src/fleur/io/xml"
SCHEMAS = {"input": "FleurInputSchema.xsd", "output": "FleurOutputSchema.xsd"}
SKIP_CHILDREN = {"output": {"fleurInput"}, "input": set()}

TYPE_NAMES = {
    "FleurBool": "logical (T/F)",
    "xsd:boolean": "logical",
    "FleurDouble": "real",
    "xsd:double": "real",
    "xsd:string": "string",
    "xsd:integer": "integer",
    "xsd:int": "integer",
    "xsd:nonNegativeInteger": "integer >= 0",
    "xsd:positiveInteger": "integer > 0",
}


def fleur_root():
    here = os.path.dirname(os.path.abspath(__file__))
    out = subprocess.run(["git", "-C", here, "rev-parse", "--show-toplevel"], capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit("cannot locate the fleur repository from " + here)
    return out.stdout.strip()


def read_file(root, ref, relpath):
    """File content at a git ref, or in the working tree if ref is None. None if absent."""
    if ref is None:
        path = os.path.join(root, relpath)
        if not os.path.exists(path):
            return None
        with open(path, encoding="utf-8") as f:
            return f.read()
    out = subprocess.run(["git", "-C", root, "show", "%s:%s" % (ref, relpath)], capture_output=True, text=True)
    return out.stdout if out.returncode == 0 else None


def local(tag):
    return tag[len(XS):] if tag.startswith(XS) else tag


def doc(node):
    if node is None:
        return ""
    d = node.find(XS + "annotation/" + XS + "documentation")
    if d is None or d.text is None:
        return ""
    return " ".join(d.text.split())


class Schema:
    def __init__(self, kind, ref=None, root=None):
        self.kind = kind
        self.root = root or fleur_root()
        self.ctypes, self.stypes, self.groups, self.roots = {}, {}, {}, []
        self.found = self._load(ref, SCHEMAS[kind], True, set())

    def _load(self, ref, filename, main, seen):
        if filename in seen:
            return True
        seen.add(filename)
        text = read_file(self.root, ref, SCHEMA_DIR + "/" + filename)
        if text is None:
            return False
        for child in ET.fromstring(text):
            tag, name = local(child.tag), child.get("name")
            if tag == "include":
                self._load(ref, child.get("schemaLocation"), False, seen)
            elif tag == "complexType":
                self.ctypes[name] = child
            elif tag == "simpleType":
                self.stypes[name] = child
            elif tag == "group":
                self.groups[name] = child
            elif tag == "element" and main:
                self.roots.append(child)
        return True

    def enum(self, typename=None, node=None, depth=0):
        """Allowed values of an enumerated simple type, else None."""
        if node is None:
            node = self.stypes.get(typename)
        if node is None or depth > 5:
            return None
        values = [e.get("value") for e in node.iter(XS + "enumeration")]
        for u in node.iter(XS + "union"):
            for member in (u.get("memberTypes") or "").split():
                values += self.enum(member, depth=depth + 1) or []
        return values or None

    def attributes(self, tnode, depth=0):
        if tnode is None or depth > 5:
            return []
        attrs = []
        for node in tnode:
            tag = local(node.tag)
            if tag == "attribute" and node.get("name"):
                attrs.append(node)
            elif tag in ("simpleContent", "complexContent", "extension", "restriction"):
                if node.get("base") in self.ctypes:
                    attrs += self.attributes(self.ctypes[node.get("base")], depth + 1)
                attrs += self.attributes(node, depth + 1)
        return attrs

    def children(self, tnode):
        """Direct child elements of a complex type, through sequence/all/choice/group."""
        result = []

        def scan(node):
            for c in node:
                tag = local(c.tag)
                if tag == "element":
                    result.append(c)
                elif tag in ("sequence", "all", "choice", "complexContent", "extension"):
                    if tag == "extension" and c.get("base") in self.ctypes:
                        scan(self.ctypes[c.get("base")])
                    scan(c)
                elif tag == "group" and c.get("ref") in self.groups:
                    scan(self.groups[c.get("ref")])

        if tnode is not None:
            scan(tnode)
        return result

    def attr_info(self, a):
        typename = a.get("type")
        inline = a.find(XS + "simpleType")
        values = self.enum(node=inline) if inline is not None else self.enum(typename)
        if typename is None and inline is not None and not values:
            r = inline.find(XS + "restriction")
            typename = (r.get("base") + " (restricted)") if r is not None and r.get("base") else "restricted"
        if typename in TYPE_NAMES:
            values = None
        base = None
        if typename in self.stypes:
            r = self.stypes[typename].find(XS + "restriction")
            base = r.get("base") if r is not None else None
        return {
            "type": typename or "enumeration",
            "base": base,
            "use": a.get("use", "optional"),
            "default": a.get("default", a.get("fixed")),
            "enum": values,
            "doc": doc(a),
        }

    def walk(self):
        """Map XPath -> element info, attributes included."""
        out = {}
        skip = SKIP_CHILDREN[self.kind]

        def visit(path, el, stack):
            typename = el.get("type")
            tnode = el.find(XS + "complexType")
            if tnode is None:
                tnode = self.ctypes.get(typename)
            out[path] = {
                "type": typename or ("(inline)" if tnode is not None else ""),
                "min": el.get("minOccurs", "1"),
                "max": el.get("maxOccurs", "1"),
                "doc": doc(el) or doc(tnode),
                "attrs": {a.get("name"): self.attr_info(a) for a in self.attributes(tnode)},
                "children": [],
            }
            key = typename or id(tnode)
            if tnode is None or key in stack:
                return
            for c in self.children(tnode):
                name = c.get("name") or c.get("ref")
                if name is None or name in skip:
                    continue
                out[path]["children"].append(name)
                visit(path + "/" + name, c, stack | {key})

        for r in self.roots:
            visit("/" + r.get("name"), r, frozenset())
        return out


def load(kind, ref=None, root=None):
    s = Schema(kind, ref, root)
    return s.walk() if s.found else {}


def kinds(arg):
    return ["input", "output"] if arg == "both" else [arg]


def fmt_attr(info):
    parts = [info["type"], info["use"]]
    if info["default"] is not None:
        parts.append('default="%s"' % info["default"])
    if info["enum"]:
        parts.append("one of " + ",".join(info["enum"]))
    return "  ".join(parts)


def cmd_diff(args):
    root = fleur_root()
    for kind in kinds(args.schema):
        base, head = load(kind, args.base, root), load(kind, args.head, root)
        lines = []
        for path in sorted(set(head) - set(base)):
            e = head[path]
            lines.append("+ element   %s  (%s, minOccurs=%s)" % (path, e["type"], e["min"]))
            if e["doc"]:
                lines.append("      " + e["doc"][:200])
            for name, info in sorted(e["attrs"].items()):
                lines.append("      @%-16s %s" % (name, fmt_attr(info)))
        for path in sorted(set(base) - set(head)):
            lines.append("- element   %s" % path)
        for path in sorted(set(base) & set(head)):
            b, h = base[path]["attrs"], head[path]["attrs"]
            for name in sorted(set(h) - set(b)):
                lines.append("+ attribute %s/@%s  %s" % (path, name, fmt_attr(h[name])))
                if h[name]["doc"]:
                    lines.append("      " + h[name]["doc"][:200])
            for name in sorted(set(b) - set(h)):
                lines.append("- attribute %s/@%s" % (path, name))
            for name in sorted(set(b) & set(h)):
                for field in ("type", "use", "default", "enum") + (("doc",) if args.docs else ()):
                    if b[name][field] != h[name][field]:
                        lines.append("~ attribute %s/@%s  %s: %s -> %s" % (path, name, field, b[name][field], h[name][field]))
        label = args.head or "working tree"
        print("== %s schema (%s): %s -> %s" % (kind, SCHEMAS[kind], args.base, label))
        print("\n".join(lines) if lines else "   no changes")


def match(paths, xpath):
    xpath = xpath.rstrip("/")
    if xpath.startswith("/"):
        return [p for p in paths if p == xpath]
    return [p for p in paths if p == "/" + xpath or p.endswith("/" + xpath)]


def type_label(info):
    if info["enum"]:
        return "one of " + ", ".join("`%s`" % v for v in info["enum"])
    if info["type"] in TYPE_NAMES:
        return TYPE_NAMES[info["type"]]
    if info.get("base"):
        return TYPE_NAMES.get(info["base"], info["base"])
    return info["type"].replace(" (restricted)", "")


def cmd_show(args):
    root = fleur_root()
    hits = 0
    for kind in kinds(args.schema):
        tree = load(kind, args.ref, root)
        for path in match(tree, args.xpath):
            hits += 1
            e = tree[path]
            name = path.rsplit("/", 1)[1]
            print("<!-- %s  (%s schema, type %s, minOccurs=%s maxOccurs=%s) -->" % (path, kind, e["type"], e["min"], e["max"]))
            if e["doc"]:
                print(e["doc"] + "\n")
            attrs = sorted(e["attrs"].items())
            example = " ".join('%s="%s"' % (n, i["default"] if i["default"] is not None else "...") for n, i in attrs)
            print("```xml\n<%s%s%s/>\n```\n" % (name, " " if example else "", example))
            if attrs:
                print("|Attribute|Type|Default|Description|")
                print("|---|---|---|---|")
                for n, i in attrs:
                    default = "required" if i["use"] == "required" else (i["default"] if i["default"] is not None else "-")
                    print("|%s|%s|%s|%s|" % (n, type_label(i), default, i["doc"].replace("|", "\\|")))
                print()
            if e["children"]:
                print("Child elements: " + ", ".join("`%s`" % c for c in e["children"]) + "\n")
    if not hits:
        sys.exit("no element matches %s" % args.xpath)


def cmd_find(args):
    root = fleur_root()
    for kind in kinds(args.schema):
        tree = load(kind, args.ref, root)
        for path, e in sorted(tree.items()):
            if path.rsplit("/", 1)[1] == args.name:
                print("%-6s element   %s" % (kind, path))
            if args.name in e["attrs"]:
                print("%-6s attribute %s/@%s  %s" % (kind, path, args.name, fmt_attr(e["attrs"][args.name])))


def main():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--schema", choices=["input", "output", "both"], default="both")
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("diff", parents=[common])
    d.add_argument("base")
    d.add_argument("head", nargs="?")
    d.add_argument("--docs", action="store_true", help="also report changed documentation strings")
    s = sub.add_parser("show", parents=[common])
    s.add_argument("xpath")
    s.add_argument("--ref")
    f = sub.add_parser("find", parents=[common])
    f.add_argument("name")
    f.add_argument("--ref")
    args = p.parse_args()
    {"diff": cmd_diff, "show": cmd_show, "find": cmd_find}[args.cmd](args)


if __name__ == "__main__":
    main()
