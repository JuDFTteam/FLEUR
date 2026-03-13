"""
FLEURiste JupyterLab Interface

Provides an interactive ipywidgets-based GUI for FLEUR inp.xml editing
in JupyterLab / Jupyter Notebook environments.

Mirrors the functionality of the Textual TUI:
  - XML Editor: schema-aware tree editing of inp.xml
  - K-Point Manager: create / inspect / visualise k-point lists
  - Inpgen: generate inp.xml from namelist or structure files
  - Job Generator: build SLURM job scripts with auto-parallelisation

Usage (inside a notebook):

    from fleuriste.jupyter_gui import FLEURisteJupyter
    gui = FLEURisteJupyter(schema="FleurInputSchema.xsd", input_file="inp.xml")
    gui.display()
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import ipywidgets as widgets
from IPython.display import display, HTML, clear_output

# ── Internal imports ─────────────────────────────────────────────────────────
from .schema_parser import XSDParser, SchemaElement, AttributeInfo
from .xml_editor import XMLDocument, XMLNode
from .kpoint_manager import (
    InpXMLManager, KPoint, KPointList, KPointMode, KPointModifiers,
)

# Optional heavy imports ---------------------------------------------------
try:
    from .inpgen_gui import (
        ase_to_fleur_input, create_supercell, parse_namelist_to_atoms,
        get_elements_in_input, add_magnetic_moments, extract_magnetic_moments,
        generate_surface,
        INPGEN_PROFILES, STRUCTURE_FORMATS,
    )
    _HAS_INPGEN_GUI = True
except ImportError:
    _HAS_INPGEN_GUI = False

try:
    from FleurInpgen import InpgenInterface
    _HAS_INPGEN = True
except ImportError:
    _HAS_INPGEN = False

try:
    from ase.io import read as ase_read
    _HAS_ASE = True
except ImportError:
    _HAS_ASE = False

try:
    import matplotlib
    import matplotlib.pyplot as plt
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False

try:
    from .pyjob.slurm_generator import (
        MachineConfig, SlurmJobConfig, SlurmJobGenerator,
        load_machine_configs, get_machines_directory,
    )
    from .pyjob.fleur_analyzer import FleurInputAnalyzer, FleurAnalysisResult
    from .pyjob.parallelization import suggest_parallelization
    _HAS_PYJOB = True
except ImportError:
    _HAS_PYJOB = False


# ── Helpers ──────────────────────────────────────────────────────────────────

def _make_btn(description: str, *, icon: str = "", tooltip: str = "",
              style: str = "", layout: dict | None = None, **kw) -> widgets.Button:
    """Shorthand button factory."""
    btn = widgets.Button(
        description=description,
        icon=icon,
        tooltip=tooltip or description,
        layout=widgets.Layout(**(layout or {})),
        **kw,
    )
    if style:
        btn.button_style = style
    return btn


def _html_tag(tag: str, text: str, **attrs: str) -> str:
    a = " ".join(f'{k}="{v}"' for k, v in attrs.items())
    return f"<{tag} {a}>{text}</{tag}>" if a else f"<{tag}>{text}</{tag}>"


# ═══════════════════════════════════════════════════════════════════════════
#  XML Editor Panel
# ═══════════════════════════════════════════════════════════════════════════

class XMLEditorPanel:
    """Schema-aware XML tree editor using ipywidgets."""

    def __init__(self, doc: XMLDocument):
        self.doc = doc
        self._selected_node: Optional[XMLNode] = None

        # ── Widgets ──
        self._tree_select = widgets.Select(
            options=[], rows=25,
            layout=widgets.Layout(width="100%"),
        )
        self._tree_select.observe(self._on_tree_select, names="value")

        self._editor_area = widgets.VBox(layout=widgets.Layout(width="100%"))
        self._status = widgets.HTML("")

        self._btn_save = _make_btn("Save", icon="save", style="success")
        self._btn_undo = _make_btn("Undo", icon="undo", style="warning")
        self._btn_redo = _make_btn("Redo", icon="repeat", style="warning")
        self._btn_add  = _make_btn("Add child", icon="plus", style="info")
        self._btn_del  = _make_btn("Delete", icon="trash", style="danger")
        self._btn_dup  = _make_btn("Duplicate", icon="clone", style="info")

        self._btn_save.on_click(self._on_save)
        self._btn_undo.on_click(self._on_undo)
        self._btn_redo.on_click(self._on_redo)
        self._btn_add.on_click(self._on_add)
        self._btn_del.on_click(self._on_delete)
        self._btn_dup.on_click(self._on_duplicate)

        self._search = widgets.Text(placeholder="Search schema…",
                                    layout=widgets.Layout(width="100%"))
        self._search.observe(self._on_search, names="value")
        self._search_results = widgets.Select(options=[], rows=5,
                                              layout=widgets.Layout(width="100%", display="none"))
        self._search_results.observe(self._on_search_select, names="value")

        self._refresh_tree()

    # ── Public ────────────────────────────────────────────────────────
    @property
    def widget(self) -> widgets.Widget:
        toolbar = widgets.HBox([
            self._btn_save, self._btn_undo, self._btn_redo,
            self._btn_add, self._btn_del, self._btn_dup,
        ])
        left = widgets.VBox([
            self._search,
            self._search_results,
            self._tree_select,
        ], layout=widgets.Layout(width="40%", overflow_y="auto"))
        right = widgets.VBox([
            self._editor_area,
        ], layout=widgets.Layout(width="60%", padding="0 0 0 12px"))
        body = widgets.HBox([left, right], layout=widgets.Layout(width="100%"))
        return widgets.VBox([toolbar, body, self._status])

    # ── Tree helpers ──────────────────────────────────────────────────
    def _flatten_tree(self, node: XMLNode, depth: int = 0) -> List[Tuple[str, XMLNode]]:
        indent = "  " * depth
        label = f"{indent}{node.get_display_name()}"
        items: List[Tuple[str, XMLNode]] = [(label, node)]
        for child in node.children:
            items.extend(self._flatten_tree(child, depth + 1))
        return items

    def _refresh_tree(self):
        if not self.doc.root:
            self._tree_select.options = [("(empty document)", None)]
            return
        items = self._flatten_tree(self.doc.root)
        # Store tuples (label, node) – Select widget uses value
        self._tree_select.options = [(lbl, node) for lbl, node in items]

    # ── Event handlers ────────────────────────────────────────────────
    def _on_tree_select(self, change):
        node = change["new"]
        if node is None:
            return
        self._selected_node = node
        self._render_editor(node)

    def _render_editor(self, node: XMLNode):
        children: List[widgets.Widget] = []

        # Title & path
        children.append(widgets.HTML(
            f"<h3 style='color:#4488cc'>{node.display_tag}</h3>"
            f"<code style='color:gray'>{node.get_path_string()}</code>"
        ))

        # Schema doc
        if node.schema_element and node.schema_element.documentation:
            children.append(widgets.HTML(
                f"<blockquote style='color:#888;border-left:3px solid #4488cc;padding-left:8px'>"
                f"{node.schema_element.documentation}</blockquote>"
            ))

        # Attributes
        if node.attributes or (node.schema_element and node.schema_element.attributes):
            children.append(widgets.HTML("<b>Attributes</b>"))
            all_attrs = dict(node.attributes)
            # Merge schema-defined but unset attrs
            if node.schema_element:
                for aname, ainfo in node.schema_element.attributes.items():
                    if aname not in all_attrs:
                        all_attrs[aname] = ainfo.default or ""
            for attr_name in sorted(all_attrs):
                value = all_attrs[attr_name]
                attr_widget = self._make_attr_widget(node, attr_name, value)
                children.append(attr_widget)

        # Text content
        if node.schema_element and node.schema_element.is_simple_content:
            children.append(widgets.HTML("<b>Content</b>"))
            txt = widgets.Textarea(value=node.text or "",
                                   layout=widgets.Layout(width="100%", height="60px"))
            txt.observe(lambda c, n=node: self._on_text_change(n, c), names="value")
            children.append(txt)
        elif node.text:
            children.append(widgets.HTML("<b>Text</b>"))
            txt = widgets.Textarea(value=node.text or "",
                                   layout=widgets.Layout(width="100%", height="60px"))
            txt.observe(lambda c, n=node: self._on_text_change(n, c), names="value")
            children.append(txt)

        # Validation
        errors = self.doc.validate_node(node)
        if errors:
            err_html = "<br>".join(f"⚠️ {e}" for e in errors)
            children.append(widgets.HTML(f"<div style='color:red;margin-top:8px'>{err_html}</div>"))

        self._editor_area.children = children

    def _make_attr_widget(self, node: XMLNode, attr_name: str, value: str) -> widgets.Widget:
        schema_attr: Optional[AttributeInfo] = None
        if node.schema_element:
            schema_attr = node.schema_element.attributes.get(attr_name)

        label_text = attr_name
        if schema_attr and schema_attr.is_required:
            label_text += " *"

        label = widgets.Label(label_text, layout=widgets.Layout(width="40%"))
        tooltip_parts: List[str] = []
        if schema_attr:
            tooltip_parts.append(f"Type: {schema_attr.type_name}")
            if schema_attr.documentation:
                tooltip_parts.append(schema_attr.documentation)

        # Choose widget type based on schema
        enum_values: List[str] = []
        is_bool = False
        if schema_attr and self.doc.schema:
            enum_values = self.doc.schema.get_enum_values(schema_attr.type_name)
            is_bool = self.doc.schema.is_boolean_type(schema_attr.type_name)

        if is_bool:
            w = widgets.Dropdown(
                options=["T", "F"],
                value=value if value in ("T", "F") else "F",
                layout=widgets.Layout(width="60%"),
            )
            w.observe(lambda c, n=node, a=attr_name: self._on_attr_change(n, a, c), names="value")
        elif enum_values:
            opts = enum_values if value in enum_values else [value] + enum_values
            w = widgets.Dropdown(
                options=opts, value=value,
                layout=widgets.Layout(width="60%"),
            )
            w.observe(lambda c, n=node, a=attr_name: self._on_attr_change(n, a, c), names="value")
        else:
            w = widgets.Text(
                value=value,
                layout=widgets.Layout(width="60%"),
            )
            w.observe(lambda c, n=node, a=attr_name: self._on_attr_change(n, a, c), names="value")

        row = widgets.HBox([label, w])
        return row

    # ── Data mutation callbacks ───────────────────────────────────────
    def _on_attr_change(self, node: XMLNode, attr: str, change):
        self.doc.save_undo_state()
        node.attributes[attr] = str(change["new"])
        self.doc.mark_modified()
        self._status.value = f"<span style='color:orange'>Modified: {node.display_tag}/@{attr}</span>"

    def _on_text_change(self, node: XMLNode, change):
        self.doc.save_undo_state()
        node.text = change["new"]
        self.doc.mark_modified()
        self._status.value = f"<span style='color:orange'>Modified text of {node.display_tag}</span>"

    def _on_save(self, _):
        try:
            self.doc.save()
            self._status.value = f"<span style='color:green'>✓ Saved to {self.doc.file_path}</span>"
        except Exception as e:
            self._status.value = f"<span style='color:red'>✗ Save failed: {e}</span>"

    def _on_undo(self, _):
        if self.doc.undo():
            self._refresh_tree()
            self._status.value = "<span style='color:blue'>↩ Undo</span>"

    def _on_redo(self, _):
        if self.doc.redo():
            self._refresh_tree()
            self._status.value = "<span style='color:blue'>↪ Redo</span>"

    def _on_add(self, _):
        if not self._selected_node:
            return
        addable = self.doc.get_addable_children(self._selected_node)
        if not addable:
            self._status.value = "<span style='color:gray'>No addable children for this element</span>"
            return
        # Show dialog
        dd = widgets.Dropdown(options=[(s.name, s) for s in addable],
                              description="Element:")
        btn = _make_btn("Add", style="success")
        out = widgets.Output()

        def do_add(_):
            schema_elem = dd.value
            self.doc.save_undo_state()
            new_node = self.doc.create_element_from_schema(schema_elem)
            self._selected_node.add_child(new_node)
            self.doc.mark_modified()
            self._refresh_tree()
            with out:
                clear_output()
                print(f"✓ Added {schema_elem.name}")

        btn.on_click(do_add)
        self._editor_area.children = [dd, btn, out]

    def _on_delete(self, _):
        node = self._selected_node
        if not node or not node.parent:
            return
        self.doc.save_undo_state()
        node.parent.remove_child(node)
        self.doc.mark_modified()
        self._selected_node = None
        self._refresh_tree()
        self._status.value = f"<span style='color:red'>Deleted {node.display_tag}</span>"

    def _on_duplicate(self, _):
        node = self._selected_node
        if not node or not node.parent:
            return
        self.doc.save_undo_state()
        clone = node.clone()
        idx = node.parent.children.index(node)
        node.parent.insert_child(idx + 1, clone)
        self.doc.mark_modified()
        self._refresh_tree()
        self._status.value = f"<span style='color:green'>Duplicated {node.display_tag}</span>"

    # ── Search ────────────────────────────────────────────────────────
    def _on_search(self, change):
        query = change["new"].strip()
        if not query or not self.doc.schema:
            self._search_results.layout.display = "none"
            return
        results = self.doc.schema.search(query)
        if results:
            self._search_results.options = [
                (f"{r.path}  ({r.match_type}: {r.match_text})", r) for r in results[:20]
            ]
            self._search_results.layout.display = ""
        else:
            self._search_results.options = [("(no results)", None)]
            self._search_results.layout.display = ""

    def _on_search_select(self, change):
        result = change["new"]
        if result is None:
            return
        # Try to navigate to the node matching the path
        path_parts = [p for p in result.path.split("/") if p]
        node = self.doc.root
        if node and path_parts:
            for part in path_parts[1:]:  # skip root
                child = node.get_child_by_tag(part)
                if child:
                    node = child
                else:
                    break
            self._selected_node = node
            self._render_editor(node)


# ═══════════════════════════════════════════════════════════════════════════
#  K-Point Manager Panel
# ═══════════════════════════════════════════════════════════════════════════

class KPointPanel:
    """Interactive k-point manager for JupyterLab."""

    def __init__(self, inp_xml: str = "inp.xml"):
        self._inp_xml = inp_xml
        self._mgr: Optional[InpXMLManager] = None
        self._output = widgets.Output()
        self._plot_output = widgets.Output()

        # ── Sidebar widgets ──
        self._list_select = widgets.Select(options=[], rows=8,
                                           layout=widgets.Layout(width="100%"))
        self._list_select.observe(self._on_list_select, names="value")

        self._btn_refresh = _make_btn("Refresh", icon="refresh", style="info")
        self._btn_select = _make_btn("Set Active", icon="check", style="success")
        self._btn_delete = _make_btn("Delete", icon="trash", style="danger")
        self._btn_save = _make_btn("Save", icon="save", style="success")

        self._btn_refresh.on_click(self._on_refresh)
        self._btn_select.on_click(self._on_set_active)
        self._btn_delete.on_click(self._on_delete_list)
        self._btn_save.on_click(self._on_save)

        # Create controls
        self._create_name = widgets.Text(value="new_kpts", description="Name:",
                                         layout=widgets.Layout(width="100%"))
        self._create_mode = widgets.Dropdown(
            options=[("Grid", "grid"), ("Density", "density"),
                     ("Number", "number"), ("Path", "path")],
            value="grid", description="Mode:",
            layout=widgets.Layout(width="100%"),
        )
        self._create_mode.observe(self._on_mode_change, names="value")

        self._grid_nx = widgets.IntText(value=8, description="nx:", layout=widgets.Layout(width="100%"))
        self._grid_ny = widgets.IntText(value=8, description="ny:", layout=widgets.Layout(width="100%"))
        self._grid_nz = widgets.IntText(value=8, description="nz:", layout=widgets.Layout(width="100%"))
        self._grid_box = widgets.VBox([self._grid_nx, self._grid_ny, self._grid_nz])

        self._density_val = widgets.FloatText(value=0.1, description="Density:",
                                              layout=widgets.Layout(width="100%"))
        self._density_box = widgets.VBox([self._density_val])
        self._density_box.layout.display = "none"

        self._number_val = widgets.IntText(value=100, description="Target N:",
                                           layout=widgets.Layout(width="100%"))
        self._number_box = widgets.VBox([self._number_val])
        self._number_box.layout.display = "none"

        self._path_str = widgets.Text(
            value="gamma=0,0,0;x=0.5,0,0.5;gamma=0,0,0",
            description="Path:", layout=widgets.Layout(width="100%"),
        )
        self._path_npts = widgets.IntText(value=100, description="Points:",
                                          layout=widgets.Layout(width="100%"))
        self._path_box = widgets.VBox([self._path_str, self._path_npts])
        self._path_box.layout.display = "none"

        self._chk_gamma = widgets.Checkbox(value=False, description="Gamma-centred")
        self._chk_sym = widgets.Checkbox(value=True, description="Use symmetry")

        self._btn_create = _make_btn("Create", icon="plus", style="primary")
        self._btn_create.on_click(self._on_create)

        # Projection selector for plot
        self._proj_select = widgets.Dropdown(
            options=["xy", "xz", "yz"], value="xy",
            description="Projection:",
            layout=widgets.Layout(width="160px"),
        )
        self._proj_select.observe(lambda _: self._plot_selected(), names="value")

        self._status = widgets.HTML("")

        # Initial load
        self._load_manager()

    # ── Public ────────────────────────────────────────────────────────
    @property
    def widget(self) -> widgets.Widget:
        create_section = widgets.VBox([
            widgets.HTML("<b>Create k-points</b>"),
            self._create_name,
            self._create_mode,
            self._grid_box,
            self._density_box,
            self._number_box,
            self._path_box,
            self._chk_gamma,
            self._chk_sym,
            self._btn_create,
        ])
        sidebar = widgets.VBox([
            widgets.HTML("<b>K-Point Lists</b>"),
            self._list_select,
            widgets.HBox([self._btn_refresh, self._btn_select]),
            widgets.HBox([self._btn_delete, self._btn_save]),
            widgets.HTML("<hr>"),
            create_section,
        ], layout=widgets.Layout(width="280px", padding="8px",
                                  border_right="2px solid #4488cc"))
        content = widgets.VBox([
            widgets.HBox([self._proj_select]),
            self._plot_output,
            self._output,
            self._status,
        ], layout=widgets.Layout(width="calc(100% - 280px)", padding="8px"))
        return widgets.HBox([sidebar, content], layout=widgets.Layout(width="100%"))

    # ── Internal ──────────────────────────────────────────────────────
    def _load_manager(self):
        path = Path(self._inp_xml)
        if path.exists():
            self._mgr = InpXMLManager(str(path), quiet=True)
            self._refresh_list_widget()
        else:
            self._status.value = f"<span style='color:orange'>⚠ {self._inp_xml} not found</span>"

    def _refresh_list_widget(self):
        if not self._mgr:
            return
        opts = []
        for name, kplist in self._mgr.kpoint_lists.items():
            marker = " ← active" if name == self._mgr.selected_list else ""
            opts.append((f"{name} ({len(kplist.kpoints)} pts){marker}", name))
        self._list_select.options = opts

    def _on_list_select(self, change):
        self._show_list_details(change["new"])
        self._plot_selected()

    def _show_list_details(self, name: str):
        if not self._mgr or not name:
            return
        kplist = self._mgr.get_kpoint_list(name)
        if not kplist:
            return
        with self._output:
            clear_output(wait=True)
            # Build HTML table
            rows = []
            from .kpoint_gui import evaluate_coordinate
            for i, kp in enumerate(kplist.kpoints):
                rows.append(
                    f"<tr><td>{i}</td>"
                    f"<td>{evaluate_coordinate(kp.x):.8f}</td>"
                    f"<td>{evaluate_coordinate(kp.y):.8f}</td>"
                    f"<td>{evaluate_coordinate(kp.z):.8f}</td>"
                    f"<td>{evaluate_coordinate(kp.weight):.6f}</td>"
                    f"<td>{kp.label or ''}</td></tr>"
                )
            html = (
                f"<h4>{name} — {len(kplist.kpoints)} k-points (type: {kplist.type})</h4>"
                "<table border='1' style='border-collapse:collapse;font-size:12px'>"
                "<tr><th>#</th><th>kx</th><th>ky</th><th>kz</th>"
                "<th>weight</th><th>label</th></tr>"
                + "".join(rows)
                + "</table>"
            )
            display(HTML(html))

    def _plot_selected(self):
        name = self._list_select.value
        if not self._mgr or not name:
            return
        kplist = self._mgr.get_kpoint_list(name)
        if not kplist or not kplist.kpoints:
            return
        proj = self._proj_select.value
        from .kpoint_gui import evaluate_coordinate

        if proj == "xy":
            xs = [evaluate_coordinate(kp.x) for kp in kplist.kpoints]
            ys = [evaluate_coordinate(kp.y) for kp in kplist.kpoints]
            xl, yl = "kx", "ky"
        elif proj == "xz":
            xs = [evaluate_coordinate(kp.x) for kp in kplist.kpoints]
            ys = [evaluate_coordinate(kp.z) for kp in kplist.kpoints]
            xl, yl = "kx", "kz"
        else:
            xs = [evaluate_coordinate(kp.y) for kp in kplist.kpoints]
            ys = [evaluate_coordinate(kp.z) for kp in kplist.kpoints]
            xl, yl = "ky", "kz"

        with self._plot_output:
            clear_output(wait=True)
            if _HAS_MPL:
                fig, ax = plt.subplots(figsize=(5, 4))
                ax.scatter(xs, ys, s=12, marker="x")
                ax.set_xlabel(xl)
                ax.set_ylabel(yl)
                ax.set_title(f"{name} ({len(kplist.kpoints)} pts, {proj.upper()})")
                ax.set_aspect("equal", adjustable="datalim")
                fig.tight_layout()
                plt.show()
            else:
                display(HTML(f"<em>Install matplotlib for k-point plots</em>"))

    # ── Actions ───────────────────────────────────────────────────────
    def _on_refresh(self, _):
        self._load_manager()
        self._status.value = "<span style='color:green'>Refreshed</span>"

    def _on_set_active(self, _):
        name = self._list_select.value
        if self._mgr and name:
            self._mgr.select_kpoint_list(name)
            self._mgr.save()
            self._refresh_list_widget()
            self._status.value = f"<span style='color:green'>✓ '{name}' is now active</span>"

    def _on_delete_list(self, _):
        name = self._list_select.value
        if not self._mgr or not name:
            return
        if name == self._mgr.selected_list:
            self._status.value = "<span style='color:red'>Cannot delete the active list</span>"
            return
        self._mgr.remove_kpoint_list(name)
        self._mgr.save()
        self._refresh_list_widget()
        self._status.value = f"<span style='color:green'>✓ Deleted '{name}'</span>"

    def _on_save(self, _):
        if self._mgr:
            self._mgr.save()
            self._status.value = "<span style='color:green'>✓ Saved</span>"

    def _on_mode_change(self, change):
        mode = change["new"]
        self._grid_box.layout.display = "" if mode == "grid" else "none"
        self._density_box.layout.display = "" if mode == "density" else "none"
        self._number_box.layout.display = "" if mode == "number" else "none"
        self._path_box.layout.display = "" if mode == "path" else "none"

    def _on_create(self, _):
        if not self._mgr:
            self._status.value = "<span style='color:red'>No inp.xml loaded</span>"
            return
        name = self._create_name.value.strip()
        if not name:
            self._status.value = "<span style='color:red'>Enter a list name</span>"
            return
        mode_str = self._create_mode.value
        try:
            if mode_str == "grid":
                self._mgr.create_kpoints(
                    name=name, mode=KPointMode.GRID,
                    modifiers=KPointModifiers(gamma=self._chk_gamma.value),
                    grid=(self._grid_nx.value, self._grid_ny.value, self._grid_nz.value),
                    use_symmetry=self._chk_sym.value,
                )
            elif mode_str == "density":
                self._mgr.create_kpoints(
                    name=name, mode=KPointMode.DENSITY,
                    modifiers=KPointModifiers(gamma=self._chk_gamma.value),
                    density=self._density_val.value,
                    use_symmetry=self._chk_sym.value,
                )
            elif mode_str == "number":
                self._mgr.create_kpoints(
                    name=name, mode=KPointMode.NUMBER,
                    modifiers=KPointModifiers(gamma=self._chk_gamma.value),
                    num_kpoints=self._number_val.value,
                    use_symmetry=self._chk_sym.value,
                )
            elif mode_str == "path":
                self._mgr.create_kpoint_path(
                    name=name, special_points=[],
                    num_points=self._path_npts.value,
                    path_string=self._path_str.value,
                )
            self._mgr.save()
            self._refresh_list_widget()
            nk = len(self._mgr.get_kpoint_list(name).kpoints) if self._mgr.get_kpoint_list(name) else "?"
            self._status.value = f"<span style='color:green'>✓ Created '{name}' → {nk} k-points</span>"
        except Exception as e:
            self._status.value = f"<span style='color:red'>✗ {e}</span>"


# ═══════════════════════════════════════════════════════════════════════════
#  Inpgen Panel
# ═══════════════════════════════════════════════════════════════════════════

# Lazy-import surface helpers (they live in inpgen_gui and need ASE)
def _surface_presets():
    try:
        from .inpgen_gui import SURFACE_PRESETS, ELEMENT_LATTICE_CONSTANTS
        return SURFACE_PRESETS, ELEMENT_LATTICE_CONSTANTS
    except Exception:
        return [], {}


class InpgenPanel:
    """Generate inp.xml from namelist text or structure files.

    Mirrors the TUI Inpgen mode: load files from disk or upload, convert
    structures via ASE, generate surfaces, create supercells, set magnetic
    moments, and run inpgen to produce ``inp.xml``.
    """

    def __init__(self):
        self._output = widgets.Output()
        # Keeps a reference to the last ASE Atoms loaded (for supercell etc.)
        self._current_atoms = None

        # ── Generation options ────────────────────────────────────────
        profiles = list(INPGEN_PROFILES) if _HAS_INPGEN_GUI else [("default_profile", "Default Profile")]
        self._profile = widgets.Dropdown(
            options=[(label, val) for val, label in profiles],
            value=profiles[0][0],
            description="Profile:", layout=widgets.Layout(width="100%"),
        )

        fmt_options = [("FLEUR namelist (none)", "")]
        if _HAS_INPGEN_GUI:
            fmt_options += [(label, val) for val, label in STRUCTURE_FORMATS]
        self._format = widgets.Dropdown(
            options=fmt_options, value="",
            description="ASE format:", layout=widgets.Layout(width="100%"),
        )

        self._chk_nosym = widgets.Checkbox(value=False, description="No symmetry")
        self._chk_film  = widgets.Checkbox(value=False, description="Film mode")

        # ── Input sources ─────────────────────────────────────────────
        # 1) Local file path (works on HPC / when notebook runs next to files)
        self._file_path = widgets.Text(
            value="", placeholder="path/to/input_file or structure",
            description="File path:", layout=widgets.Layout(width="100%"),
        )
        self._btn_load_file = _make_btn("Load from path", icon="folder-open", style="primary")
        self._btn_load_file.on_click(self._on_load_file)

        # 2) Browser upload (works when notebook runs remotely)
        self._file_upload = widgets.FileUpload(
            accept=".txt,.inp,.cif,.vasp,.xyz,.pdb,.xsf,.json,POSCAR",
            multiple=False, description="Upload file",
            layout=widgets.Layout(width="100%"),
        )
        self._file_upload.observe(self._on_file_upload, names="value")

        self._btn_process_upload = _make_btn(
            "⬆ Process uploaded file", icon="upload", style="warning")
        self._btn_process_upload.on_click(self._on_process_upload)

        # ── Text area for namelist input ──────────────────────────────
        self._textarea = widgets.Textarea(
            placeholder=(
                "Paste FLEUR namelist input here, or use one of the load\n"
                "buttons above.\n\n"
                "If you load a structure file (CIF, POSCAR, …), set the\n"
                "ASE format dropdown and it will be auto-converted to a\n"
                "FLEUR namelist."
            ),
            layout=widgets.Layout(width="100%", height="300px"),
        )

        # ── Output dir ────────────────────────────────────────────────
        self._outdir = widgets.Text(value=".", description="Output dir:",
                                    layout=widgets.Layout(width="100%"))

        # ── Action buttons ────────────────────────────────────────────
        self._btn_generate = _make_btn("▶ Generate inp.xml", icon="cog", style="success")
        self._btn_generate.on_click(self._on_generate)

        self._btn_preview = _make_btn("Preview namelist", icon="eye", style="info")
        self._btn_preview.on_click(self._on_preview)

        # ── Structure tools ───────────────────────────────────────────
        self._btn_surface   = _make_btn("🎬 Generate Surface", icon="globe", style="")
        self._btn_supercell = _make_btn("🔲 Supercell", icon="th", style="")
        self._btn_magnetic  = _make_btn("🧲 Magnetic Moments", icon="magnet", style="")

        self._btn_surface.on_click(self._on_surface)
        self._btn_supercell.on_click(self._on_supercell)
        self._btn_magnetic.on_click(self._on_magnetic)

        # Surface sub-panel (hidden by default)
        self._surface_panel = self._build_surface_panel()
        self._surface_panel.layout.display = "none"

        # Supercell sub-panel (hidden by default)
        self._supercell_panel = self._build_supercell_panel()
        self._supercell_panel.layout.display = "none"

        # Magnetic moments sub-panel (hidden by default)
        self._magnetic_panel = self._build_magnetic_panel()
        self._magnetic_panel.layout.display = "none"

        self._status = widgets.HTML("")

    # ── Layout ────────────────────────────────────────────────────────
    @property
    def widget(self) -> widgets.Widget:
        sidebar = widgets.VBox([
            widgets.HTML("<b style='color:#4488cc'>Generate inp.xml</b>"),
            self._btn_generate,
            widgets.HTML("<hr>"),

            widgets.HTML("<b style='color:#4488cc'>Generation Options</b>"),
            self._profile,
            self._chk_nosym,
            widgets.HTML("<hr>"),

            widgets.HTML("<b style='color:#4488cc'>Input Source</b>"),
            self._file_path,
            self._btn_load_file,
            self._file_upload,
            self._btn_process_upload,
            self._format,
            self._chk_film,
            widgets.HTML("<hr>"),

            widgets.HTML("<b style='color:#4488cc'>Structure Tools</b>"),
            self._btn_surface,
            self._btn_supercell,
            self._btn_magnetic,
            widgets.HTML("<hr>"),

            self._outdir,
            self._btn_preview,
            self._status,
        ], layout=widgets.Layout(width="300px", padding="8px",
                                  border_right="2px solid #4488cc",
                                  overflow_y="auto"))

        content = widgets.VBox([
            widgets.HTML("<b>Input Content</b> (edit namelist input below):"),
            self._textarea,
            self._surface_panel,
            self._supercell_panel,
            self._magnetic_panel,
            widgets.HTML("<b>Output Log:</b>"),
            self._output,
        ], layout=widgets.Layout(width="calc(100% - 300px)", padding="8px"))

        return widgets.HBox([sidebar, content], layout=widgets.Layout(width="100%"))

    # ══════════════════════════════════════════════════════════════════
    #  Loading files
    # ══════════════════════════════════════════════════════════════════

    def _load_content_from_path(self, filepath: Path):
        """Read a file from disk, converting via ASE if a format is selected."""
        fmt = self._format.value  # "" means plain namelist

        if fmt:
            # Structure file → convert to namelist via ASE
            if not _HAS_ASE:
                self._status.value = (
                    "<span style='color:red'>ASE is required for structure files. "
                    "Install with: pip install ase</span>")
                return
            try:
                ase_fmt = fmt if fmt != "auto" else None
                atoms = ase_read(str(filepath), format=ase_fmt)
                self._current_atoms = atoms
                namelist = ase_to_fleur_input(atoms, film=self._chk_film.value)
                self._textarea.value = namelist
                self._status.value = (
                    f"<span style='color:green'>Loaded structure: {filepath.name} — "
                    f"{atoms.get_chemical_formula()}, {len(atoms)} atoms</span>")
            except Exception as e:
                self._status.value = f"<span style='color:red'>ASE error: {e}</span>"
        else:
            # Plain text (namelist / input file)
            try:
                text = filepath.read_text()
                self._textarea.value = text
                self._current_atoms = None
                # Try to parse as namelist for later tools
                if _HAS_INPGEN_GUI:
                    atoms = parse_namelist_to_atoms(text)
                    if atoms is not None:
                        self._current_atoms = atoms
                self._status.value = f"<span style='color:green'>Loaded {filepath.name}</span>"
            except Exception as e:
                self._status.value = f"<span style='color:red'>Read error: {e}</span>"

    def _on_load_file(self, _):
        """Load a file from the local filesystem using the path text field."""
        raw = self._file_path.value.strip()
        if not raw:
            self._status.value = "<span style='color:red'>Enter a file path first</span>"
            return
        filepath = Path(raw).expanduser()
        if not filepath.exists():
            self._status.value = f"<span style='color:red'>Not found: {filepath}</span>"
            return
        self._load_content_from_path(filepath)

    def _process_uploaded_file(self, fname: str, raw_content, fmt: str):
        """Process a single uploaded file (called from _on_file_upload or button)."""
        import tempfile

        # Ensure we have bytes (ipywidgets 8.x gives memoryview)
        if isinstance(raw_content, memoryview):
            content_bytes = bytes(raw_content)
        elif isinstance(raw_content, bytes):
            content_bytes = raw_content
        else:
            content_bytes = bytes(raw_content)

        if fmt:
            # Structure file — write to temp file so ASE can read it
            suffix = Path(fname).suffix or f".{fmt}"
            with tempfile.NamedTemporaryFile(
                mode="wb", suffix=suffix, delete=False
            ) as tmp:
                tmp.write(content_bytes)
                tmp_path = tmp.name
            try:
                self._load_content_from_path(Path(tmp_path))
                # Override status to show original filename
                if "color:green" in self._status.value:
                    self._status.value = self._status.value.replace(
                        Path(tmp_path).name, fname)
            finally:
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass
        else:
            # Plain text — decode directly
            try:
                text = content_bytes.decode("utf-8", errors="replace")
            except AttributeError:
                text = str(content_bytes)
            self._textarea.value = text
            self._current_atoms = None
            if _HAS_INPGEN_GUI:
                atoms = parse_namelist_to_atoms(text)
                if atoms is not None:
                    self._current_atoms = atoms
            self._status.value = f"<span style='color:green'>Loaded {fname}</span>"

    def _on_file_upload(self, change):
        """Handle browser file upload widget (fires automatically on file selection)."""
        uploaded = change.get("new", None)
        if not uploaded:
            return

        fmt = self._format.value

        # ipywidgets 8.x: value is a tuple of dicts with 'name', 'content' (memoryview)
        # ipywidgets 7.x: value is a dict mapping filename → bytes
        if isinstance(uploaded, (list, tuple)):
            for item in uploaded:
                if isinstance(item, dict):
                    fname = item.get("name", "uploaded_file")
                    raw_content = item.get("content", b"")
                    self._process_uploaded_file(fname, raw_content, fmt)
                    break  # only first file
        elif isinstance(uploaded, dict):
            # ipywidgets 7.x style: {filename: content}
            for fname, raw_content in uploaded.items():
                self._process_uploaded_file(fname, raw_content, fmt)
                break  # only first file
        else:
            self._status.value = (
                "<span style='color:red'>Unexpected upload format — "
                "check ipywidgets version</span>")

    def _on_process_upload(self, _):
        """Explicit button to process uploaded file (in case observe doesn't trigger)."""
        uploaded = self._file_upload.value
        if not uploaded:
            self._status.value = (
                "<span style='color:orange'>No file uploaded yet — "
                "click the Upload button first to select a file</span>")
            return

        fmt = self._format.value

        if isinstance(uploaded, (list, tuple)):
            for item in uploaded:
                if isinstance(item, dict):
                    fname = item.get("name", "uploaded_file")
                    raw_content = item.get("content", b"")
                    self._process_uploaded_file(fname, raw_content, fmt)
                    break
        elif isinstance(uploaded, dict):
            for fname, raw_content in uploaded.items():
                self._process_uploaded_file(fname, raw_content, fmt)
                break

    # ══════════════════════════════════════════════════════════════════
    #  Preview & Generate
    # ══════════════════════════════════════════════════════════════════

    def _get_namelist_content(self) -> Optional[str]:
        """Return the namelist content ready for inpgen.

        If the textarea already holds namelist text, return it directly.
        """
        text = self._textarea.value.strip()
        if not text:
            self._status.value = "<span style='color:red'>No input content</span>"
            return None
        return text

    def _on_preview(self, _):
        content = self._get_namelist_content()
        if content is None:
            return
        with self._output:
            clear_output(wait=True)
            print(content)

    def _on_generate(self, _):
        if not _HAS_INPGEN:
            self._status.value = (
                "<span style='color:red'>FleurInpgen library not available. "
                "Build FLEUR with Python bindings.</span>")
            return

        content = self._get_namelist_content()
        if content is None:
            return

        output_dir = Path(self._outdir.value).resolve()
        output_dir.mkdir(parents=True, exist_ok=True)
        original_dir = os.getcwd()

        try:
            os.chdir(str(output_dir))
            inpgen = InpgenInterface(quiet=True)
            inpgen.make_inp(content, self._profile.value, self._chk_nosym.value)
            messages = inpgen.get_messages()
            target = output_dir / "inp.xml"

            with self._output:
                clear_output(wait=True)
                if messages.strip():
                    print(messages)
                if target.exists():
                    print(f"\n✓ Generated {target}")
                else:
                    print("\n✗ inp.xml was not created — check input")

            if target.exists():
                self._status.value = f"<span style='color:green'>✓ Generated {target}</span>"
            else:
                self._status.value = "<span style='color:red'>✗ inp.xml was not created</span>"
        except Exception as e:
            self._status.value = f"<span style='color:red'>✗ {e}</span>"
            with self._output:
                clear_output(wait=True)
                print(f"Error: {e}")
        finally:
            os.chdir(original_dir)

    # ══════════════════════════════════════════════════════════════════
    #  Surface Generator
    # ══════════════════════════════════════════════════════════════════

    def _build_surface_panel(self) -> widgets.Widget:
        presets, elem_consts = _surface_presets()

        elem_opts = [(f"{sym} ({st}, a={a:.2f} Å)", sym)
                     for sym, (a, st) in elem_consts.items()] if elem_consts else [("Fe", "Fe")]

        self._surf_element = widgets.Dropdown(
            options=elem_opts, description="Element:",
            layout=widgets.Layout(width="100%"))
        self._surf_type = widgets.Dropdown(
            options=[(label, val) for val, label in presets] if presets else [("BCC (100)", "bcc100")],
            description="Surface:", layout=widgets.Layout(width="100%"))
        self._surf_lattice = widgets.FloatText(value=2.87, description="a (Å):",
                                                layout=widgets.Layout(width="100%"))
        self._surf_layers = widgets.IntText(value=5, description="Layers:",
                                             layout=widgets.Layout(width="100%"))
        self._surf_vacuum = widgets.FloatText(value=10.0, description="Vacuum (Å):",
                                               layout=widgets.Layout(width="100%"))
        self._surf_size_x = widgets.IntText(value=1, description="Size x:",
                                             layout=widgets.Layout(width="50%"))
        self._surf_size_y = widgets.IntText(value=1, description="Size y:",
                                             layout=widgets.Layout(width="50%"))
        # Miller indices (for custom)
        self._surf_miller_h = widgets.IntText(value=1, description="h:",
                                               layout=widgets.Layout(width="33%"))
        self._surf_miller_k = widgets.IntText(value=0, description="k:",
                                               layout=widgets.Layout(width="33%"))
        self._surf_miller_l = widgets.IntText(value=0, description="l:",
                                               layout=widgets.Layout(width="33%"))
        self._surf_miller_box = widgets.HBox(
            [self._surf_miller_h, self._surf_miller_k, self._surf_miller_l])
        self._surf_miller_box.layout.display = "none"

        self._surf_type.observe(self._on_surf_type_change, names="value")
        self._surf_element.observe(self._on_surf_element_change, names="value")

        btn_go = _make_btn("Generate surface", style="primary")
        btn_go.on_click(self._on_surf_generate)
        btn_close = _make_btn("Close", style="")
        btn_close.on_click(lambda _: setattr(self._surface_panel.layout, "display", "none"))

        return widgets.VBox([
            widgets.HTML("<b>🎬 Generate Surface / Film</b>"),
            self._surf_element, self._surf_type,
            self._surf_miller_box,
            self._surf_lattice, self._surf_layers, self._surf_vacuum,
            widgets.HBox([self._surf_size_x, self._surf_size_y]),
            widgets.HBox([btn_go, btn_close]),
        ], layout=widgets.Layout(
            border="1px solid #4488cc", padding="8px", margin="8px 0"))

    def _on_surf_type_change(self, change):
        self._surf_miller_box.layout.display = "" if change["new"] == "custom" else "none"

    def _on_surf_element_change(self, change):
        _, elem_consts = _surface_presets()
        sym = change["new"]
        if sym in elem_consts:
            a, struct = elem_consts[sym]
            self._surf_lattice.value = a
            mapping = {"bcc": "bcc100", "fcc": "fcc100", "hcp": "hcp0001", "diamond": "diamond100"}
            if struct in mapping:
                self._surf_type.value = mapping[struct]

    def _on_surface(self, _):
        self._surface_panel.layout.display = ""
        self._supercell_panel.layout.display = "none"
        self._magnetic_panel.layout.display = "none"

    def _on_surf_generate(self, _):
        if not _HAS_ASE or not _HAS_INPGEN_GUI:
            self._status.value = "<span style='color:red'>ASE required for surface generation</span>"
            return
        try:
            miller = None
            if self._surf_type.value == "custom":
                miller = (self._surf_miller_h.value,
                          self._surf_miller_k.value,
                          self._surf_miller_l.value)

            slab = generate_surface(
                element=self._surf_element.value,
                surface_type=self._surf_type.value,
                layers=self._surf_layers.value,
                vacuum=self._surf_vacuum.value,
                lattice_const=self._surf_lattice.value,
                size=(self._surf_size_x.value, self._surf_size_y.value, 1),
                orthogonal=True,
                miller=miller,
            )
            self._current_atoms = slab
            namelist = ase_to_fleur_input(slab, film=True)
            self._textarea.value = namelist
            self._chk_film.value = True

            miller_str = (f"({miller[0]}{miller[1]}{miller[2]})"
                          if miller else self._surf_type.value)
            self._status.value = (
                f"<span style='color:green'>✓ {self._surf_element.value} "
                f"{miller_str} surface — {len(slab)} atoms</span>")
        except Exception as e:
            self._status.value = f"<span style='color:red'>✗ {e}</span>"

    # ══════════════════════════════════════════════════════════════════
    #  Supercell
    # ══════════════════════════════════════════════════════════════════

    def _build_supercell_panel(self) -> widgets.Widget:
        self._sc_nx = widgets.IntText(value=2, description="a₁ ×:",
                                       layout=widgets.Layout(width="33%"))
        self._sc_ny = widgets.IntText(value=2, description="a₂ ×:",
                                       layout=widgets.Layout(width="33%"))
        self._sc_nz = widgets.IntText(value=2, description="a₃ ×:",
                                       layout=widgets.Layout(width="33%"))
        self._sc_info = widgets.HTML("")

        btn_go = _make_btn("Create supercell", style="primary")
        btn_go.on_click(self._on_sc_generate)
        btn_close = _make_btn("Close", style="")
        btn_close.on_click(lambda _: setattr(self._supercell_panel.layout, "display", "none"))

        return widgets.VBox([
            widgets.HTML("<b>🔲 Create Supercell</b>"),
            widgets.HBox([self._sc_nx, self._sc_ny, self._sc_nz]),
            self._sc_info,
            widgets.HBox([btn_go, btn_close]),
        ], layout=widgets.Layout(
            border="1px solid #4488cc", padding="8px", margin="8px 0"))

    def _on_supercell(self, _):
        self._supercell_panel.layout.display = ""
        self._surface_panel.layout.display = "none"
        self._magnetic_panel.layout.display = "none"
        # Try to parse current text to get atom count
        if self._current_atoms is not None:
            self._sc_info.value = f"Current structure: {len(self._current_atoms)} atoms"
        else:
            text = self._textarea.value.strip()
            if text and _HAS_INPGEN_GUI:
                atoms = parse_namelist_to_atoms(text)
                if atoms is not None:
                    self._current_atoms = atoms
                    self._sc_info.value = f"Current structure: {len(atoms)} atoms"
                else:
                    self._sc_info.value = "<em>Could not parse structure from textarea</em>"
            else:
                self._sc_info.value = "<em>Load a structure first</em>"

    def _on_sc_generate(self, _):
        if not _HAS_ASE or not _HAS_INPGEN_GUI:
            self._status.value = "<span style='color:red'>ASE required for supercell</span>"
            return

        # Parse current input if we don't have atoms yet
        if self._current_atoms is None:
            text = self._textarea.value.strip()
            if text:
                self._current_atoms = parse_namelist_to_atoms(text)
        if self._current_atoms is None:
            self._status.value = "<span style='color:red'>Load a structure first</span>"
            return

        nx, ny, nz = self._sc_nx.value, self._sc_ny.value, self._sc_nz.value
        if nx < 1 or ny < 1 or nz < 1:
            self._status.value = "<span style='color:red'>All values must be ≥ 1</span>"
            return

        try:
            # Preserve magnetic moments from current text
            mag_moments = {}
            text = self._textarea.value.strip()
            if text:
                mag_moments = extract_magnetic_moments(text)

            sc = create_supercell(self._current_atoms, nx, ny, nz)
            self._current_atoms = sc
            namelist = ase_to_fleur_input(sc, film=self._chk_film.value,
                                          magnetic_moments=mag_moments)
            self._textarea.value = namelist
            self._status.value = (
                f"<span style='color:green'>✓ {nx}×{ny}×{nz} supercell — "
                f"{len(sc)} atoms</span>")
        except Exception as e:
            self._status.value = f"<span style='color:red'>✗ {e}</span>"

    # ══════════════════════════════════════════════════════════════════
    #  Magnetic Moments
    # ══════════════════════════════════════════════════════════════════

    def _build_magnetic_panel(self) -> widgets.Widget:
        self._mag_element = widgets.Dropdown(
            options=[("(load a structure first)", 0)],
            description="Element:", layout=widgets.Layout(width="100%"),
        )
        self._mag_type = widgets.Dropdown(
            options=[("Collinear", "col"), ("Non-collinear", "nc")],
            value="col", description="Type:",
            layout=widgets.Layout(width="100%"),
        )
        self._mag_type.observe(self._on_mag_type_change, names="value")

        self._mag_moment = widgets.Text(value="2.0", description="Moment (μB):",
                                         layout=widgets.Layout(width="100%"))
        self._mag_preset = widgets.Dropdown(
            options=[("Custom value", "custom"), ("up", "up"), ("down", "down")],
            value="custom", description="Preset:",
            layout=widgets.Layout(width="100%"),
        )
        self._mag_col_box = widgets.VBox([self._mag_preset, self._mag_moment])

        self._mag_mx = widgets.Text(value="0.0", description="mx:",
                                     layout=widgets.Layout(width="33%"))
        self._mag_my = widgets.Text(value="0.0", description="my:",
                                     layout=widgets.Layout(width="33%"))
        self._mag_mz = widgets.Text(value="2.0", description="mz:",
                                     layout=widgets.Layout(width="33%"))
        self._mag_nc_box = widgets.HBox([self._mag_mx, self._mag_my, self._mag_mz])
        self._mag_nc_box.layout.display = "none"

        btn_apply = _make_btn("Apply", style="primary")
        btn_apply.on_click(self._on_mag_apply)
        btn_close = _make_btn("Close", style="")
        btn_close.on_click(lambda _: setattr(self._magnetic_panel.layout, "display", "none"))

        return widgets.VBox([
            widgets.HTML("<b>🧲 Set Magnetic Moments</b>"),
            self._mag_element, self._mag_type,
            self._mag_col_box, self._mag_nc_box,
            widgets.HBox([btn_apply, btn_close]),
        ], layout=widgets.Layout(
            border="1px solid #4488cc", padding="8px", margin="8px 0"))

    def _on_mag_type_change(self, change):
        is_nc = change["new"] == "nc"
        self._mag_col_box.layout.display = "none" if is_nc else ""
        self._mag_nc_box.layout.display = "" if is_nc else "none"

    def _on_magnetic(self, _):
        self._magnetic_panel.layout.display = ""
        self._surface_panel.layout.display = "none"
        self._supercell_panel.layout.display = "none"

        text = self._textarea.value.strip()
        if not text or not _HAS_INPGEN_GUI:
            self._status.value = "<span style='color:orange'>Load a structure first</span>"
            return
        elements = get_elements_in_input(text)
        if elements:
            self._mag_element.options = [
                (f"{sym} (Z={z})", z) for z, sym in elements]
        else:
            self._mag_element.options = [("(no atoms found)", 0)]

    def _on_mag_apply(self, _):
        element_z = self._mag_element.value
        if not element_z:
            return
        text = self._textarea.value.strip()
        if not text:
            return

        if self._mag_type.value == "nc":
            moment_spec = f"{self._mag_mx.value} {self._mag_my.value} {self._mag_mz.value}"
        else:
            preset = self._mag_preset.value
            if preset in ("up", "down"):
                moment_spec = preset
            else:
                moment_spec = self._mag_moment.value.strip()

        try:
            new_text = add_magnetic_moments(text, element_z, moment_spec)
            self._textarea.value = new_text
            from .inpgen_gui import ATOMIC_SYMBOLS
            sym = ATOMIC_SYMBOLS.get(element_z, f"Z{element_z}")
            self._status.value = (
                f"<span style='color:green'>✓ Set moment for {sym}: {moment_spec}</span>")
        except Exception as e:
            self._status.value = f"<span style='color:red'>✗ {e}</span>"


# ═══════════════════════════════════════════════════════════════════════════
#  Job Generator Panel
# ═══════════════════════════════════════════════════════════════════════════

class JobGeneratorPanel:
    """SLURM job script generator for FLEUR calculations."""

    _LABEL_W = "150px"
    _FIELD_W = "calc(100% - 158px)"

    def __init__(self, input_file: str = "inp.xml"):
        self._input_file = input_file
        self._output = widgets.Output()
        self._preview_area = widgets.Textarea(
            value="", placeholder="Job script preview…",
            layout=widgets.Layout(width="100%", height="400px"),
        )

        # ── Load machines ──────────────────────────────────────────────
        self._machines = {}
        if _HAS_PYJOB:
            try:
                self._machines = load_machine_configs()
            except Exception:
                pass

        # ── Helper to build a labeled row ─────────────────────────────
        def _row(label, widget):
            lbl = widgets.HTML(
                f"<div style='line-height:32px;font-weight:500'>{label}</div>",
                layout=widgets.Layout(width=self._LABEL_W, flex_shrink="0"),
            )
            return widgets.HBox([lbl, widget],
                                layout=widgets.Layout(width="100%", margin="2px 0"))

        # ── Section header helper stored as instance method ──────────
        def _section(title):
            return widgets.HTML(
                f"<div style='background:#2255884d;padding:4px 8px;border-radius:4px;"
                f"font-weight:bold;margin:8px 0 4px'>{title}</div>")
        self._section = _section

        # ── Machine & Partition ───────────────────────────────────────
        machine_opts = [("(custom)", "")] + [(n, n) for n in sorted(self._machines.keys())]
        self._machine_select = widgets.Dropdown(
            options=machine_opts, value="",
            layout=widgets.Layout(width=self._FIELD_W),
        )
        self._partition_select = widgets.Dropdown(
            options=[("(select machine first)", "")], value="",
            layout=widgets.Layout(width=self._FIELD_W),
        )
        self._machine_info = widgets.HTML("")
        self._machine_select.observe(self._on_machine_change, names="value")
        self._partition_select.observe(self._on_partition_change, names="value")

        # ── Job Identity ──────────────────────────────────────────────
        self._job_name = widgets.Text(value="fleur_job",
                                      layout=widgets.Layout(width=self._FIELD_W))
        self._account = widgets.Text(placeholder="(optional)",
                                     layout=widgets.Layout(width=self._FIELD_W))
        self._time = widgets.Text(value="01:00:00", placeholder="HH:MM:SS",
                                  layout=widgets.Layout(width=self._FIELD_W))

        # ── Resources ────────────────────────────────────────────────
        self._nodes = widgets.IntText(value=1,
                                      layout=widgets.Layout(width=self._FIELD_W))
        self._ntasks = widgets.IntText(value=1,
                                       layout=widgets.Layout(width=self._FIELD_W))
        self._cpus_per_task = widgets.IntText(value=1,
                                              layout=widgets.Layout(width=self._FIELD_W))
        self._memory = widgets.Text(value="4G", placeholder="e.g. 4G, 8000M",
                                    layout=widgets.Layout(width=self._FIELD_W))
        self._gpus = widgets.Text(placeholder="e.g. 1  or  a100:2  (optional)",
                                  layout=widgets.Layout(width=self._FIELD_W))

        # ── Notifications ─────────────────────────────────────────────
        self._mail_user = widgets.Text(placeholder="user@example.com (optional)",
                                       layout=widgets.Layout(width=self._FIELD_W))
        self._mail_type = widgets.Dropdown(
            options=["END,FAIL", "BEGIN,END,FAIL", "ALL", "NONE"],
            value="END,FAIL",
            layout=widgets.Layout(width=self._FIELD_W),
        )

        # ── Modules & Commands ────────────────────────────────────────
        self._modules_area = widgets.Textarea(
            placeholder="One module per line",
            layout=widgets.Layout(width="100%", height="80px"),
        )
        self._commands_area = widgets.Textarea(
            placeholder="Shell commands to execute",
            layout=widgets.Layout(width="100%", height="80px"),
        )

        # ── Buttons & status ──────────────────────────────────────────
        self._btn_analyze  = _make_btn("Analyze inp.xml", icon="search",    style="info")
        self._btn_generate = _make_btn("Generate",        icon="file-text", style="success")
        self._btn_save     = _make_btn("Save .slurm",     icon="save",      style="primary")
        self._btn_analyze.on_click(self._on_analyze)
        self._btn_generate.on_click(self._on_generate)
        self._btn_save.on_click(self._on_save)

        self._status       = widgets.HTML("")
        self._analysis_html = widgets.HTML("")
        self._par_html     = widgets.HTML("")

        self._fleur_result: Optional[Any] = None
        self._machine_config: Optional[Any] = None

        # ── Build form layout ─────────────────────────────────────────
        self._form = widgets.VBox([
            self._section("🖥️ Machine Configuration"),
            _row("Machine:", self._machine_select),
            _row("Partition:", self._partition_select),
            self._machine_info,

            self._section("📋 Job Identity"),
            _row("Job name:", self._job_name),
            _row("Account:", self._account),
            _row("Time limit:", self._time),

            self._section("⚡ Parallelization & Resources"),
            _row("Nodes:", self._nodes),
            _row("Tasks:", self._ntasks),
            _row("CPUs/task:", self._cpus_per_task),
            _row("Memory:", self._memory),
            _row("GPUs:", self._gpus),

            self._section("📧 Notifications"),
            _row("Email:", self._mail_user),
            _row("Notify on:", self._mail_type),

            self._section("📦 Modules"),
            self._modules_area,

            self._section("🔧 Commands"),
            self._commands_area,

            widgets.HTML("<hr style='margin:8px 0'>"),
            widgets.HBox([self._btn_analyze, self._btn_generate, self._btn_save]),
            self._status,
        ], layout=widgets.Layout(width="360px", padding="8px",
                                  border_right="2px solid #4488cc",
                                  overflow_y="auto"))

    @property
    def widget(self) -> widgets.Widget:
        right = widgets.VBox([
            self._section("🔍 FLEUR Input Analysis"),
            self._analysis_html,
            self._par_html,
            widgets.HTML("<b>Script Preview</b>"),
            self._preview_area,
            self._output,
        ], layout=widgets.Layout(flex="1", padding="8px", min_width="0"))
        return widgets.HBox([self._form, right],
                            layout=widgets.Layout(width="100%", align_items="flex-start"))

    # ── Callbacks ─────────────────────────────────────────────────────
    def _on_machine_change(self, change):
        name = change["new"]
        if name and name in self._machines:
            self._machine_config = self._machines[name]
            mc = self._machine_config
            # Populate partition dropdown
            parts = mc.list_partitions() if hasattr(mc, "list_partitions") else []
            if parts:
                self._partition_select.options = [("(default)", "")] + [(p, p) for p in parts]
                self._partition_select.value = ""
            else:
                self._partition_select.options = [("(no partitions)", "")]
                self._partition_select.value = ""
            info = (f"<div style='border:1px solid #4488cc;padding:6px;border-radius:4px'>"
                    f"<b>{mc.name}</b><br>"
                    f"Cores/node: {mc.cores_per_node} | "
                    f"Mem: {mc.memory_per_node_gb} GB/node")
            if mc.gpus_per_node:
                info += f" | GPUs: {mc.gpus_per_node}/node"
            info += "</div>"
            self._machine_info.value = info
            # Pre-fill modules
            if hasattr(mc, "modules_needed") and mc.modules_needed:
                self._modules_area.value = "\n".join(mc.modules_needed)
        else:
            self._machine_config = None
            self._machine_info.value = ""
            self._partition_select.options = [("(select machine first)", "")]
            self._partition_select.value = ""

    def _on_partition_change(self, change):
        """Update machine info / re-auto-suggest when partition changes."""
        if self._machine_config and self._fleur_result and _HAS_PYJOB:
            self._auto_suggest()

    def _auto_suggest(self):
        """Fill resource fields from parallelization suggestion."""
        if not (self._machine_config and self._fleur_result):
            return
        try:
            par = suggest_parallelization(
                fleur_result=self._fleur_result,
                machine=self._machine_config,
                partition=self._partition_select.value or None,
                initial_nodes=self._nodes.value,
            )
            self._nodes.value        = par.num_nodes
            self._ntasks.value       = par.num_tasks
            self._cpus_per_task.value = par.cpus_per_task
            self._memory.value       = f"{int(par.memory_per_node_gb * 1.2) + 1}G"
            self._par_html.value = (
                f"<div style='border:1px solid #00aa44;padding:6px;border-radius:4px'>"
                f"<b>💡 Auto-suggestion</b><br>"
                f"Nodes: {par.num_nodes} | Tasks: {par.num_tasks}<br>"
                f"CPUs/task: {par.cpus_per_task} | "
                f"Mem/node: {par.memory_per_node_gb:.1f} GB"
                f"</div>"
            )
        except Exception as e:
            self._par_html.value = f"<span style='color:orange'>Suggestion failed: {e}</span>"

    def _on_analyze(self, _):
        if not _HAS_PYJOB:
            self._status.value = "<span style='color:red'>pyjob module not available</span>"
            return
        path = Path(self._input_file)
        if not path.exists():
            self._status.value = f"<span style='color:red'>{self._input_file} not found</span>"
            return
        try:
            analyzer = FleurInputAnalyzer(str(path))
            result = analyzer.analyze()
            self._fleur_result = result
            elem_size = 16 if result.is_complex else 8
            mem_gb = result.matrix_dimension ** 2 * elem_size / (1024 ** 3)
            self._analysis_html.value = (
                f"<div style='border:1px solid #cc8800;padding:6px;border-radius:4px'>"
                f"<b>FLEUR Analysis</b><br>"
                f"Kmax: {result.kmax:.2f} | Atoms: {result.num_atoms}<br>"
                f"K-pts: {result.num_kpoints} | Spins: {result.jspins}<br>"
                f"Matrix dim: {result.matrix_dimension:,}<br>"
                f"Est. mem/k: {mem_gb:.3f} GB ({'complex' if result.is_complex else 'real'})"
                f"</div>"
            )
            if self._machine_config:
                self._auto_suggest()
            self._status.value = "<span style='color:green'>✓ Analysis complete</span>"
        except Exception as e:
            self._status.value = f"<span style='color:red'>✗ {e}</span>"

    def _on_generate(self, _):
        if not _HAS_PYJOB:
            self._status.value = "<span style='color:red'>pyjob not available</span>"
            return
        try:
            mc  = self._machine_config
            partition = self._partition_select.value or None

            gpus_val = self._gpus.value.strip() or None
            gpu_syntax = "gpus"
            use_mem_option = True
            modules = []
            commands = []
            if mc:
                gpu_syntax = mc.get_effective_value("gpu_syntax", partition) or "gpus"
                use_mem = mc.get_effective_value("use_mem_option", partition)
                use_mem_option = use_mem if use_mem is not None else True
                cmd = mc.get_effective_value("command", partition)
                if cmd:
                    commands = [cmd]
                modules_cfg = mc.get_effective_value("modules_needed", partition) or []
                if modules_cfg:
                    modules = list(modules_cfg)
            # Manual module/command overrides from text areas
            if self._modules_area.value.strip():
                modules = [l for l in self._modules_area.value.splitlines() if l.strip()]
            if self._commands_area.value.strip():
                commands = [l for l in self._commands_area.value.splitlines() if l.strip()]

            mail_user = self._mail_user.value.strip() or None
            mail_type = self._mail_type.value if mail_user else None

            config = SlurmJobConfig(
                job_name      = self._job_name.value,
                nodes         = self._nodes.value,
                ntasks        = self._ntasks.value,
                cpus_per_task = self._cpus_per_task.value,
                memory        = self._memory.value,
                time          = self._time.value,
                partition     = partition,
                account       = self._account.value or None,
                gpus          = gpus_val,
                gpu_syntax    = gpu_syntax,
                use_mem_option = use_mem_option,
                modules       = modules,
                commands      = commands,
                mail_user     = mail_user,
                mail_type     = mail_type,
            )
            generator = SlurmJobGenerator(config, mc)
            script = generator.generate()
            self._preview_area.value = script
            self._status.value = "<span style='color:green'>✓ Script generated</span>"
        except Exception as e:
            self._status.value = f"<span style='color:red'>✗ {e}</span>"

    def _on_save(self, _):
        script = self._preview_area.value.strip()
        if not script:
            self._status.value = "<span style='color:red'>Generate a script first</span>"
            return
        fname = f"{self._job_name.value}.slurm"
        Path(fname).write_text(script)
        self._status.value = f"<span style='color:green'>✓ Saved {fname}</span>"



# ═══════════════════════════════════════════════════════════════════════════
#  Main Application Widget
# ═══════════════════════════════════════════════════════════════════════════

class FLEURisteJupyter:
    """
    Top-level FLEURiste interface for JupyterLab.

    Usage::

        from fleuriste.jupyter_gui import FLEURisteJupyter
        gui = FLEURisteJupyter(schema="FleurInputSchema.xsd", input_file="inp.xml")
        gui.display()
    """

    def __init__(
        self,
        schema: Optional[str] = None,
        input_file: Optional[str] = None,
    ):
        self._schema_path = schema
        self._input_file = input_file or "inp.xml"

        # Auto-detect schema
        if self._schema_path is None:
            from .cli import find_schema
            found = find_schema()
            if found:
                self._schema_path = str(found)

        # Load documents
        self._schema: Optional[XSDParser] = None
        if self._schema_path and Path(self._schema_path).exists():
            self._schema = XSDParser(self._schema_path)

        self._doc = XMLDocument(schema=self._schema)
        if Path(self._input_file).exists():
            self._doc.load(self._input_file)

        # Build panels
        self._xml_panel = XMLEditorPanel(self._doc)
        self._kpoint_panel = KPointPanel(self._input_file)
        self._inpgen_panel = InpgenPanel()
        self._job_panel = JobGeneratorPanel(self._input_file)

        # Tab container
        self._tabs = widgets.Tab()
        self._tabs.children = [
            self._inpgen_panel.widget,
            self._xml_panel.widget,
            self._kpoint_panel.widget,
            self._job_panel.widget,
        ]
        self._tabs.set_title(0, "⚙️ Inpgen")
        self._tabs.set_title(1, "📝 XML Editor")
        self._tabs.set_title(2, "🔑 K-Points")
        self._tabs.set_title(3, "🖥️ Job Generator")

        # Header
        self._header = widgets.HTML(
            "<div style='background:linear-gradient(90deg,#225588,#4488cc);"
            "color:white;padding:8px 16px;border-radius:6px 6px 0 0;"
            "font-size:16px;font-weight:bold'>"
            "🌱 FLEURiste — grow your FLEUR input"
            "<span style='float:right;font-size:12px;font-weight:normal'>"
            "JupyterLab Edition v0.1.0</span></div>"
        )

    # ── Public API ────────────────────────────────────────────────────

    def display(self):
        """Show the full GUI in the notebook output."""
        display(widgets.VBox([
            self._header,
            self._tabs,
        ], layout=widgets.Layout(border="1px solid #4488cc", border_radius="6px")))

    # Convenience accessors ------------------------------------------------

    @property
    def document(self) -> XMLDocument:
        """Access the underlying XMLDocument."""
        return self._doc

    @property
    def kpoint_manager(self) -> Optional[InpXMLManager]:
        """Access the k-point manager."""
        return self._kpoint_panel._mgr

    def reload(self):
        """Reload inp.xml and refresh all panels."""
        if Path(self._input_file).exists():
            self._doc.load(self._input_file)
            self._xml_panel._refresh_tree()
            self._kpoint_panel._load_manager()
