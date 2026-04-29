#!/usr/bin/env python3
"""
FLEURiste MCP Server – Model Context Protocol server for FLEUR inp.xml editing.

Exposes the FLEURiste CLI/library capabilities as MCP tools so that LLM-based
agents can search the FLEUR XML schema, read/edit inp.xml files, manage k-point
lists, validate input, and analyse FLEUR calculations.

Transport: Streamable-HTTP (default) or SSE, suitable for running inside a
container and connecting from any MCP client.
"""

import json
import logging
import os
import sys
import textwrap
from pathlib import Path
from typing import Any

from mcp.server.fastmcp import FastMCP

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("fleuriste-mcp")

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
SCHEMA_PATH = os.environ.get(
    "FLEURISTE_SCHEMA",
    "/app/FleurInputSchema.xsd",
)
WORKDIR = Path(os.environ.get("FLEURISTE_WORKDIR", "/work"))

# ---------------------------------------------------------------------------
# Lazy-loaded singletons
# ---------------------------------------------------------------------------
_schema_parser = None
_xml_document = None


def _get_schema():
    """Return (and cache) the XSD schema parser."""
    global _schema_parser
    if _schema_parser is None:
        from fleuriste.schema_parser import XSDParser
        if not Path(SCHEMA_PATH).exists():
            raise FileNotFoundError(
                f"Schema file not found at {SCHEMA_PATH}. "
                "Set FLEURISTE_SCHEMA or mount it into the container."
            )
        _schema_parser = XSDParser(SCHEMA_PATH)
        log.info("Loaded schema from %s", SCHEMA_PATH)
    return _schema_parser


def _get_document(inp_xml: str | None = None) -> "XMLDocument":
    """Return (and cache) the XML document, optionally (re-)loading *inp_xml*."""
    global _xml_document
    from fleuriste.xml_editor import XMLDocument

    if _xml_document is None:
        _xml_document = XMLDocument(schema=_get_schema())

    if inp_xml is not None:
        path = _resolve_path(inp_xml)
        _xml_document.load(path)
        log.info("Loaded inp.xml from %s", path)
    elif _xml_document.root is None:
        # Try default location
        default = WORKDIR / "inp.xml"
        if default.exists():
            _xml_document.load(default)
            log.info("Auto-loaded %s", default)
    return _xml_document


def _resolve_path(rel: str) -> Path:
    """Resolve a path relative to WORKDIR."""
    p = Path(rel)
    if not p.is_absolute():
        p = WORKDIR / p
    return p


# ---------------------------------------------------------------------------
# MCP Application
# ---------------------------------------------------------------------------
_host = os.environ.get("MCP_HOST", "0.0.0.0")
_port = int(os.environ.get("MCP_PORT", "8000"))

mcp = FastMCP(
    "FLEURiste",
    instructions=textwrap.dedent("""\
        FLEURiste MCP server – tools for working with FLEUR DFT inp.xml files.

        Typical workflow:
        1. Use `schema_search` to discover available XML elements & attributes.
        2. Use `read_inp_xml` to inspect the current inp.xml.
        3. Use `set_attribute`, `add_element`, `remove_element` to edit.
        4. Use `validate` to check for schema errors.
        5. Use `save_inp_xml` to persist changes.
        6. Use `kpoints_*` tools to manage k-point lists.
        7. Use `analyze_fleur_input` for parallelisation insights.

        The server operates on files inside the /work directory by default.
    """),
    host=_host,
    port=_port,
)


# ── Schema tools ─────────────────────────────────────────────────────────────

@mcp.tool()
def schema_search(query: str, max_results: int = 20) -> str:
    """Search the FLEUR XML schema for elements, attributes, or documentation.

    Returns matching paths, element names, and documentation snippets.
    Use this to discover which XML elements and attributes are available
    and what they mean.
    """
    schema = _get_schema()
    results = schema.search(query, max_results=max_results)
    if not results:
        return f"No schema matches for '{query}'."
    lines = []
    for r in results:
        lines.append(
            f"- **{r.path}** ({r.match_type}): {r.match_text}"
            + (f"  — {r.context}" if r.context else "")
        )
    return "\n".join(lines)


@mcp.tool()
def schema_element_info(path: str) -> str:
    """Get detailed schema information for an element path.

    *path* is a ``/``-separated path like ``calculationSetup/cutoffs``.
    Returns attributes (name, type, required, docs) and allowed child elements.
    """
    schema = _get_schema()
    parts = [p for p in path.strip("/").split("/") if p]
    elem = schema.get_element_for_path(parts)
    if elem is None:
        return f"Element not found for path: {path}"

    lines = [f"## {elem.name}"]
    if elem.documentation:
        lines.append(elem.documentation)

    if elem.attributes:
        lines.append("\n### Attributes")
        for name, attr in elem.attributes.items():
            req = "required" if attr.is_required else "optional"
            typ = schema.get_type_description(attr.type_name)
            default = f", default={attr.default}" if attr.default else ""
            doc = f"  {attr.documentation}" if attr.documentation else ""
            lines.append(f"- **{name}** ({typ}, {req}{default}){doc}")

    if elem.children:
        lines.append("\n### Child elements")
        for ch in elem.children:
            occ = "required" if ch.is_required else "optional"
            if ch.is_unbounded:
                occ += ", repeatable"
            lines.append(f"- **{ch.name}** ({occ})")

    return "\n".join(lines)


# ── inp.xml reading / writing ────────────────────────────────────────────────

@mcp.tool()
def read_inp_xml(file: str = "inp.xml") -> str:
    """Read and return the full contents of an inp.xml file.

    The file path is relative to /work unless an absolute path is given.
    """
    doc = _get_document(file)
    return doc.to_string()


@mcp.tool()
def get_xml_tree(file: str | None = None, max_depth: int = 3) -> str:
    """Return a compact tree view of the inp.xml structure.

    Shows element names, key attributes and nesting up to *max_depth* levels.
    Useful for getting an overview before targeted edits.
    """
    doc = _get_document(file)
    if doc.root is None:
        return "No document loaded."

    lines: list[str] = []

    def _walk(node, depth):
        if depth > max_depth:
            return
        indent = "  " * depth
        attrs = " ".join(f'{k}="{v}"' for k, v in list(node.attributes.items())[:4])
        text_hint = ""
        if node.text and node.text.strip():
            t = node.text.strip()
            text_hint = f'  "{t[:60]}{"…" if len(t) > 60 else ""}"'
        lines.append(f"{indent}<{node.display_tag} {attrs}>{text_hint}")
        for ch in node.children:
            _walk(ch, depth + 1)

    _walk(doc.root, 0)
    return "\n".join(lines)


@mcp.tool()
def get_element(xpath: str, file: str | None = None) -> str:
    """Return the XML subtree for elements matching *xpath* (simple tag path).

    *xpath* is a ``/``-separated tag path, e.g.
    ``fleurInput/calculationSetup/cutoffs``.
    """
    doc = _get_document(file)
    if doc.root is None:
        return "No document loaded."

    parts = [p for p in xpath.strip("/").split("/") if p]
    # Walk the tree
    nodes = [doc.root] if doc.root.tag == parts[0] else []
    for part in parts[1:]:
        next_nodes = []
        for n in nodes:
            next_nodes.extend(n.get_children_by_tag(part))
        nodes = next_nodes

    if not nodes:
        return f"No element found for path: {xpath}"

    from fleuriste.xml_editor import XMLNode
    import xml.etree.ElementTree as ET
    from xml.dom import minidom

    results = []
    for node in nodes:
        elem = node.to_element()
        raw = ET.tostring(elem, encoding="unicode")
        pretty = minidom.parseString(raw).toprettyxml(indent="  ")
        # Drop XML declaration line
        results.append("\n".join(pretty.split("\n")[1:]))
    return "\n---\n".join(results)


@mcp.tool()
def set_attribute(xpath: str, attribute: str, value: str, file: str | None = None) -> str:
    """Set an attribute on the element(s) at *xpath*.

    Example: ``set_attribute("fleurInput/calculationSetup/cutoffs", "Kmax", "4.5")``
    """
    doc = _get_document(file)
    if doc.root is None:
        return "No document loaded."

    parts = [p for p in xpath.strip("/").split("/") if p]
    nodes = [doc.root] if doc.root.tag == parts[0] else []
    for part in parts[1:]:
        next_nodes = []
        for n in nodes:
            next_nodes.extend(n.get_children_by_tag(part))
        nodes = next_nodes

    if not nodes:
        return f"No element found for path: {xpath}"

    for n in nodes:
        old = n.attributes.get(attribute, "(unset)")
        n.attributes[attribute] = value
    doc.mark_modified()

    return f"Set {attribute}={value} on {len(nodes)} element(s) (was {old})."


@mcp.tool()
def set_text(xpath: str, text: str, file: str | None = None) -> str:
    """Set the text content of the element(s) at *xpath*."""
    doc = _get_document(file)
    if doc.root is None:
        return "No document loaded."

    parts = [p for p in xpath.strip("/").split("/") if p]
    nodes = [doc.root] if doc.root.tag == parts[0] else []
    for part in parts[1:]:
        next_nodes = []
        for n in nodes:
            next_nodes.extend(n.get_children_by_tag(part))
        nodes = next_nodes

    if not nodes:
        return f"No element found for path: {xpath}"

    for n in nodes:
        n.text = text
    doc.mark_modified()
    return f"Set text on {len(nodes)} element(s)."


@mcp.tool()
def add_element(parent_xpath: str, tag: str, attributes: str = "{}", text: str = "", file: str | None = None) -> str:
    """Add a new child element under *parent_xpath*.

    *attributes* is a JSON object string, e.g. ``{"name": "default", "count": "10"}``.
    If the schema is loaded the element is created with required defaults.
    """
    doc = _get_document(file)
    if doc.root is None:
        return "No document loaded."

    parts = [p for p in parent_xpath.strip("/").split("/") if p]
    nodes = [doc.root] if doc.root.tag == parts[0] else []
    for part in parts[1:]:
        next_nodes = []
        for n in nodes:
            next_nodes.extend(n.get_children_by_tag(part))
        nodes = next_nodes

    if not nodes:
        return f"No parent element found for path: {parent_xpath}"

    attrs = json.loads(attributes) if attributes else {}
    from fleuriste.xml_editor import XMLNode

    count = 0
    for parent in nodes:
        # Try schema-aware creation first
        new_node = None
        if parent.schema_element:
            for ch_schema in parent.schema_element.children:
                if ch_schema.name == tag:
                    new_node = doc.create_element_from_schema(ch_schema)
                    break
        if new_node is None:
            new_node = XMLNode(tag=tag)
        # Override with user-supplied attributes
        new_node.attributes.update(attrs)
        if text:
            new_node.text = text
        parent.add_child(new_node)
        count += 1

    doc.mark_modified()
    return f"Added <{tag}> to {count} parent(s)."


@mcp.tool()
def remove_element(xpath: str, file: str | None = None) -> str:
    """Remove element(s) at *xpath* from the document.

    Removes only the *last* matching element if multiple exist.
    """
    doc = _get_document(file)
    if doc.root is None:
        return "No document loaded."

    parts = [p for p in xpath.strip("/").split("/") if p]
    nodes = [doc.root] if doc.root.tag == parts[0] else []
    for part in parts[1:]:
        next_nodes = []
        for n in nodes:
            next_nodes.extend(n.get_children_by_tag(part))
        nodes = next_nodes

    if not nodes:
        return f"No element found for path: {xpath}"

    target = nodes[-1]
    if target.parent is None:
        return "Cannot remove the root element."

    target.parent.remove_child(target)
    doc.mark_modified()
    return f"Removed <{target.tag}> from <{target.parent.tag if target.parent else '?'}>."


@mcp.tool()
def validate(file: str | None = None) -> str:
    """Validate the current inp.xml against the FLEUR schema.

    Returns a list of errors or a success message.
    """
    doc = _get_document(file)
    if doc.root is None:
        return "No document loaded."

    all_errors: list[str] = []

    def _check(node):
        errs = doc.validate_node(node)
        for e in errs:
            all_errors.append(f"{node.get_path_string()}: {e}")
        for ch in node.children:
            _check(ch)

    _check(doc.root)

    if not all_errors:
        return "✓ No validation errors."
    return "Validation errors:\n" + "\n".join(f"- {e}" for e in all_errors)


@mcp.tool()
def save_inp_xml(file: str = "inp.xml") -> str:
    """Save the current (possibly modified) inp.xml back to disk."""
    doc = _get_document()
    if doc.root is None:
        return "No document loaded."
    path = _resolve_path(file)
    doc.save(path)
    return f"Saved to {path}."


# ── K-point tools ────────────────────────────────────────────────────────────

@mcp.tool()
def kpoints_list(file: str = "inp.xml") -> str:
    """List all k-point sets defined in an inp.xml file."""
    from fleuriste.kpoint_manager import InpXMLManager

    path = _resolve_path(file)
    mgr = InpXMLManager(str(path), quiet=True)
    summary = mgr.get_summary()

    lines = [
        f"File: {summary['file']}",
        f"K-point lists: {summary['num_lists']}",
        f"Active list: {summary['selected_list'] or '(none)'}",
        "",
    ]
    for name, kplist in mgr.kpoint_lists.items():
        marker = " ← active" if name == summary["selected_list"] else ""
        lines.append(f"  {name} ({len(kplist.kpoints)} pts){marker}")
    return "\n".join(lines)


@mcp.tool()
def kpoints_show(name: str, file: str = "inp.xml") -> str:
    """Show the k-points in a specific k-point list."""
    from fleuriste.kpoint_manager import InpXMLManager

    path = _resolve_path(file)
    mgr = InpXMLManager(str(path), quiet=True)
    kplist = mgr.get_kpoint_list(name)
    if kplist is None:
        avail = ", ".join(mgr.kpoint_lists.keys())
        return f"'{name}' not found. Available: {avail}"

    lines = [f"K-point list: {name}", f"Points: {len(kplist.kpoints)}", ""]
    header = f"{'#':>5}  {'kx':>12}  {'ky':>12}  {'kz':>12}  {'weight':>10}  label"
    lines.append(header)
    lines.append("─" * 70)
    for i, kp in enumerate(kplist.kpoints):
        kx = float(kp.x) if not isinstance(kp.x, str) else kp.x
        ky = float(kp.y) if not isinstance(kp.y, str) else kp.y
        kz = float(kp.z) if not isinstance(kp.z, str) else kp.z
        w = float(kp.weight) if not isinstance(kp.weight, str) else kp.weight
        lines.append(f"{i:>5}  {kx:>12}  {ky:>12}  {kz:>12}  {w:>10}  {kp.label or ''}")
    return "\n".join(lines)


@mcp.tool()
def kpoints_select(name: str, file: str = "inp.xml") -> str:
    """Set a k-point list as the active list in inp.xml."""
    from fleuriste.kpoint_manager import InpXMLManager

    path = _resolve_path(file)
    mgr = InpXMLManager(str(path), quiet=True)
    if not mgr.select_kpoint_list(name):
        avail = ", ".join(mgr.kpoint_lists.keys())
        return f"Cannot select '{name}'. Available: {avail}"
    mgr.save()
    return f"'{name}' is now the active k-point list."


@mcp.tool()
def kpoints_delete(name: str, file: str = "inp.xml") -> str:
    """Delete a k-point list from inp.xml."""
    from fleuriste.kpoint_manager import InpXMLManager

    path = _resolve_path(file)
    mgr = InpXMLManager(str(path), quiet=True)
    if name not in mgr.kpoint_lists:
        avail = ", ".join(mgr.kpoint_lists.keys())
        return f"'{name}' not found. Available: {avail}"
    if name == mgr.get_selected_list():
        return f"Cannot delete '{name}' – it is the active list."
    mgr.remove_kpoint_list(name)
    mgr.save()
    return f"Deleted '{name}'."


# ── Analysis tools ───────────────────────────────────────────────────────────

@mcp.tool()
def analyze_fleur_input(file: str = "inp.xml") -> str:
    """Analyse a FLEUR inp.xml and return a summary of the calculation setup.

    Includes atom count, k-points, matrix size, memory estimate, SOC/noco flags, etc.
    """
    from fleuriste.pyjob.fleur_analyzer import FleurInputAnalyzer

    path = _resolve_path(file)
    analyzer = FleurInputAnalyzer(str(path))
    a = analyzer.analyze()

    elem_size = 16 if a.is_complex else 8
    mem_gb = a.matrix_dimension ** 2 * elem_size / (1024 ** 3)
    storage = "complex" if a.is_complex else "real"

    return "\n".join([
        f"File      : {file}",
        f"Kmax      : {a.kmax:.2f} a.u.⁻¹",
        f"Volume    : {a.cell_volume:.1f} Bohr³",
        f"Spins     : {a.jspins}",
        f"SOC       : {'yes' if a.has_soc else 'no'}",
        f"Noco      : {'yes' if a.has_noco else 'no'}",
        f"Inversion : {'yes' if a.has_inversion else 'no'}",
        f"K-points  : {a.num_kpoints}",
        f"Atoms     : {a.num_atoms} ({a.num_atom_types} types: {', '.join(a.atom_species)})",
        f"N_basis   : {a.n_basis_per_kpoint:,}",
        f"Matrix    : {a.matrix_dimension:,}{'  (×2 noco)' if a.has_noco else ''}",
        f"Storage   : {storage}",
        f"Memory/k  : {mem_gb:.3f} GB",
    ])


# ── Resource: list files in workdir ──────────────────────────────────────────

@mcp.tool()
def list_workdir(subdir: str = ".") -> str:
    """List files in the /work directory (or a subdirectory)."""
    target = WORKDIR / subdir
    if not target.is_dir():
        return f"{target} is not a directory."
    entries = sorted(target.iterdir())
    lines = [f"Contents of {target}:", ""]
    for e in entries:
        kind = "dir" if e.is_dir() else "file"
        lines.append(f"  [{kind}] {e.name}")
    return "\n".join(lines) if entries else "Directory is empty."


@mcp.tool()
def read_file_content(path: str) -> str:
    """Read the text content of any file in /work."""
    p = _resolve_path(path)
    if not p.exists():
        return f"File not found: {p}"
    if not p.is_file():
        return f"Not a regular file: {p}"
    return p.read_text(errors="replace")


# ── Entry point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    transport = os.environ.get("MCP_TRANSPORT", "streamable-http")

    log.info("Starting FLEURiste MCP server  transport=%s  %s:%d", transport, _host, _port)

    if transport == "sse":
        mcp.run(transport="sse")
    else:
        mcp.run(transport="streamable-http")
