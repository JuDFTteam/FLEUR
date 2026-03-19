"""
FLEURiste - Textual GUI Application for FLEUR inp.xml editing.

A terminal-based graphical editor using the Textual framework.
Includes XML editor, K-Point manager, Inpgen, and Job Generator modes.
"""

from pathlib import Path
import os
from typing import Optional, List, Dict

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Container, Horizontal, Vertical, VerticalScroll
from textual.widgets import (
    Header, Footer, Static, Tree, Input, Label, 
    Select, Rule, Button, DataTable, TabbedContent, TabPane, TextArea
)
from textual import on

from .schema_parser import XSDParser, SchemaElement, AttributeInfo, SearchResult
from .xml_editor import XMLDocument, XMLNode
from .kpoint_gui import (
    KPointPlot, CreateDensityDialog, CreatePathDialog, 
    CreateMeshDialog, STANDARD_KPOINTS
)
from .kpoint_manager import InpXMLManager, KPointMode, KPointModifiers
from .inpgen_gui import (
    FilePickerDialog, SupercellDialog, MagneticMomentDialog, SurfaceDialog,
    INPGEN_PROFILES, STRUCTURE_FORMATS,
    ase_to_fleur_input, create_supercell, parse_namelist_to_atoms,
    get_elements_in_input, add_magnetic_moments, extract_magnetic_moments,
    generate_surface, ASE_AVAILABLE
)
from .inpgen_loader import create_inpgen_interface

# Import pyjob components for Job Generator mode
try:
    from .pyjob.slurm_generator import (
        MachineConfig, Partition, SlurmJobConfig, SlurmJobGenerator,
        load_machine_configs, get_machines_directory
    )
    from .pyjob.fleur_analyzer import FleurInputAnalyzer, FleurAnalysisResult
    from .pyjob.parallelization import ParallelizationStrategy, ParallelizationResult, suggest_parallelization
    HAS_PYJOB = True
except ImportError:
    HAS_PYJOB = False
    MachineConfig = None


class FLEURisteApp(App):
    """FLEURiste - Main application for editing FLEUR inp.xml files."""
    
    TITLE = "FLEURiste - grow your FLEUR input"
    
    CSS = """
    Screen {
        layout: vertical;
    }
    
    /* Menu bar at the top */
    #menu-bar {
        height: 3;
        background: $primary-darken-1;
        padding: 0 1;
        layout: horizontal;
    }
    
    #menu-bar Button {
        margin: 0 1;
        min-width: 16;
    }
    
    #menu-bar .active-mode {
        background: $primary;
        text-style: bold;
    }
    
    /* Main content area */
    #main-content {
        height: 1fr;
        layout: horizontal;
    }
    
    /* XML Editor styles */
    #xml-editor-container {
        width: 100%;
        height: 100%;
        layout: horizontal;
    }
    
    #tree-panel {
        width: 40%;
        height: 100%;
        border-right: solid $primary;
    }
    
    #search-container {
        height: auto;
        padding: 1;
        background: $surface;
    }
    
    #search-input {
        width: 100%;
    }
    
    #search-results {
        height: auto;
        max-height: 15;
        display: none;
        background: $surface-darken-1;
        border: solid $primary;
        padding: 1;
        overflow-y: auto;
    }
    
    #search-results.visible {
        display: block;
    }
    
    .search-result {
        padding: 0 1;
        margin-bottom: 1;
    }
    
    .search-result:hover {
        background: $primary-darken-2;
    }
    
    .search-result-path {
        color: $accent;
    }
    
    .search-result-type {
        color: $text-muted;
        text-style: italic;
    }
    
    .search-result-context {
        color: $text-muted;
    }
    
    #xml-tree {
        height: 1fr;
    }
    
    #editor-panel {
        width: 60%;
        height: 100%;
        padding: 1 2;
    }
    
    .node-title {
        text-style: bold;
        margin-bottom: 1;
        color: $accent;
        text-align: center;
    }
    
    .node-path {
        color: $text-muted;
        margin-bottom: 1;
    }
    
    .section-title {
        margin-top: 1;
        margin-bottom: 1;
        text-style: bold;
        color: $secondary;
    }
    
    .attr-row {
        height: auto;
        margin-bottom: 1;
        padding: 0 1;
    }
    
    .attr-label {
        width: 40%;
        padding-top: 1;
    }
    
    .attr-input {
        width: 60%;
    }
    
    .schema-info {
        color: $text-muted;
        margin: 1;
    }
    
    .doc-text {
        color: $text;
        margin: 1;
        padding: 1;
        background: $surface-darken-1;
        border: solid $primary-darken-2;
    }
    
    .attr-doc {
        color: $text-muted;
        margin-left: 4;
        margin-bottom: 1;
        text-style: italic;
    }
    
    #no-selection {
        text-align: center;
        margin-top: 10;
        color: $text-muted;
    }
    
    #tree-buttons {
        height: auto;
        dock: bottom;
        padding: 1;
        background: $surface;
        layout: grid;
        grid-size: 4;
        grid-gutter: 1;
    }
    
    #tree-buttons Button {
        min-width: 10;
    }
    
    #no-add-elements {
        color: $text-muted;
        padding: 1;
        column-span: 4;
    }
    
    #add-elements-label {
        color: $text;
        padding: 0 1;
        text-style: bold;
        column-span: 4;
    }
    
    /* K-Point Manager styles */
    #kpoint-manager-container {
        width: 100%;
        height: 100%;
        layout: horizontal;
        display: none;
    }
    
    #kpoint-manager-container.visible {
        display: block;
        layout: horizontal;
    }
    
    #kpoint-sidebar {
        width: 30;
        background: $panel;
        padding: 1;
        border-right: solid $primary;
    }
    
    #kpoint-content {
        width: 1fr;
        padding: 1;
    }
    
    #kpoint-info-panel {
        height: auto;
        margin-bottom: 1;
        padding: 1;
        background: $panel;
        border: solid $primary;
    }
    
    #kpoint-table {
        height: 1fr;
    }
    
    #kpoint-visualization-panel {
        height: 28;
        margin-bottom: 1;
        padding: 1;
        background: $panel;
        border: solid $primary;
    }
    
    #kpoint-plot-controls {
        height: 3;
        margin-bottom: 1;
    }
    
    #kpoint-plot {
        height: 1fr;
        overflow: auto;
    }
    
    #kpoint-sidebar Button {
        margin-bottom: 1;
        width: 100%;
    }
    
    /* Inpgen mode styles */
    #inpgen-container {
        width: 100%;
        height: 100%;
        layout: horizontal;
        display: none;
    }
    
    #inpgen-container.visible {
        display: block;
        layout: horizontal;
    }
    
    #inpgen-sidebar {
        width: 35;
        height: 100%;
        background: $panel;
        padding: 1;
        border-right: solid $primary;
    }
    
    #inpgen-content {
        width: 1fr;
        padding: 1;
    }
    
    #inpgen-input-area {
        height: 1fr;
        border: solid $primary;
    }
    
    #inpgen-output-area {
        height: 12;
        margin-top: 1;
        padding: 1;
        background: $surface-darken-1;
        border: solid $primary;
        overflow-y: auto;
    }
    
    #inpgen-sidebar Button {
        margin-bottom: 1;
        width: 100%;
    }
    
    #inpgen-sidebar Select {
        margin-bottom: 1;
        width: 100%;
    }
    
    .inpgen-section-title {
        text-style: bold;
        color: $secondary;
        margin-top: 1;
        margin-bottom: 1;
    }
    
    /* Job Generator mode styles */
    #job-generator-container {
        width: 100%;
        height: 100%;
        layout: horizontal;
        display: none;
    }
    
    #job-generator-container.visible {
        display: block;
        layout: horizontal;
    }
    
    #job-sidebar {
        width: 45;
        height: 100%;
        background: $panel;
        padding: 1;
        border-right: solid $primary;
    }
    
    #job-content {
        width: 1fr;
        height: 100%;
        padding: 1;
    }
    
    #job-form-container {
        padding: 1 2;
    }
    
    .job-form-section-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 0;
    }
    
    .job-form-row {
        height: 3;
        layout: horizontal;
    }
    
    .job-form-row Label {
        width: 18;
        height: 1;
        content-align: left middle;
    }
    
    .job-form-row Input {
        width: 1fr;
        height: 3;
    }
    
    .job-form-row Select {
        width: 1fr;
        height: 3;
    }
    
    #job-machine-info {
        height: 8;
        border: solid $primary-darken-2;
        padding: 0 1;
        margin: 1 0;
        overflow-y: auto;
    }
    
    #job-preview-area {
        height: 1fr;
        border: solid $success;
        margin: 1;
    }
    
    #job-command-list {
        height: 8;
        border: solid $primary-darken-2;
    }
    
    #job-module-list {
        height: 5;
        border: solid $primary-darken-2;
    }
    
    .job-button-bar {
        height: 4;
        align: center middle;
        padding: 0 1;
    }
    
    .job-button-bar Button {
        margin: 0 1;
    }
    
    #job-fleur-analysis {
        height: auto;
        max-height: 16;
        border: solid $warning;
        padding: 0 1;
        margin: 1 0;
        overflow-y: auto;
    }
    
    #job-parallelization-info {
        height: auto;
        max-height: 12;
        border: solid $success;
        padding: 0 1;
        margin: 1 0;
        overflow-y: auto;
    }
    
    /* Status bar */
    #status-bar {
        dock: bottom;
        height: 1;
        background: $primary-darken-2;
        padding: 0 1;
    }
    """
    
    BINDINGS = [
        Binding("ctrl+s", "save", "Save"),
        Binding("ctrl+z", "undo", "Undo"),
        Binding("ctrl+y", "redo", "Redo"),
        Binding("ctrl+q", "quit", "Quit"),
        Binding("f5", "refresh", "Refresh"),
        Binding("ctrl+f", "focus_search", "Search"),
        Binding("delete", "delete_element", "Delete"),
        Binding("escape", "clear_search", "Clear Search", show=False),
        Binding("f1", "switch_to_xml_editor", "XML Editor", show=False),
        Binding("f2", "switch_to_kpoint_manager", "K-Points", show=False),
        Binding("f3", "switch_to_inpgen", "Inpgen", show=False),
        Binding("f4", "switch_to_job_generator", "Job Generator", show=False),
    ]
    
    def __init__(
        self, 
        schema_path: Optional[str | Path] = None,
        input_file: Optional[str | Path] = None,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.schema: Optional[XSDParser] = None
        self.document = XMLDocument()
        self.schema_path = schema_path
        self.input_file = input_file
        self._node_map: Dict[str, XMLNode] = {}
        self._current_node: Optional[XMLNode] = None
        self._editor_version = 0  # Counter for unique widget IDs
        self._search_version = 0  # Counter for unique search result IDs
        
        # Auto-detect inp.xml in current directory if no file specified
        if not input_file:
            default_inp = Path.cwd() / "inp.xml"
            if default_inp.exists():
                self.input_file = default_inp
        
        # Current view: "xml_editor", "kpoint_manager", or "inpgen"
        # Default to inpgen if no input file, otherwise xml_editor
        self._active_view = "inpgen" if not self.input_file else "xml_editor"
        
        # K-Point manager state
        self._kpoint_manager: Optional[InpXMLManager] = None
        self._selected_kpoint_list: Optional[str] = None
        self._current_projection: str = "xy"
        
        # Inpgen state
        self._inpgen_profile = "default_profile"
        self._inpgen_nosym = False
        self._inpgen_film = False
        self._inpgen_format = "auto"
        self._overwrite_confirmed = False
        
        # Job generator state
        self._job_machine: Optional[MachineConfig] = None
        self._job_available_machines = load_machine_configs() if HAS_PYJOB else {}
        self._job_fleur_analysis = None
        self._job_parallelization = None
        self._job_updating_from_parallelization = False
        self._job_last_machine_command = ""
        self._job_last_machine_modules = ""
        
        if schema_path:
            self.schema = XSDParser(schema_path)
            self.document.schema = self.schema
    
    def compose(self) -> ComposeResult:
        yield Header()
        
        # Menu bar for switching modes
        with Horizontal(id="menu-bar"):
            yield Button("🔧 Inpgen", id="menu-inpgen")
            yield Button("📄 XML Editor", id="menu-xml-editor")
            yield Button("📊 K-Points", id="menu-kpoint-manager")
            yield Button("🖥️ Job", id="menu-job-generator")
            yield Static("", id="menu-spacer")
        
        # Main content area
        with Container(id="main-content"):
            # Inpgen container (shown by default when no input file)
            with Horizontal(id="inpgen-container"):
                with VerticalScroll(id="inpgen-sidebar"):
                    yield Label("[bold cyan]Generate inp.xml[/bold cyan]")
                    yield Button("▶ Generate inp.xml", id="inpgen-btn-generate", variant="success")
                    yield Button("📄 Open Generated File", id="inpgen-btn-open", variant="primary")
                    
                    yield Rule()
                    yield Label("[bold cyan]Generation Options[/bold cyan]")
                    yield Label("Profile:", classes="inpgen-section-title")
                    yield Select(
                        [(name, value) for value, name in INPGEN_PROFILES],
                        value="default_profile",
                        id="inpgen-profile-select"
                    )
                    yield Button("☐ No Symmetry", id="inpgen-btn-nosym", variant="default")
                    
                    yield Rule()
                    yield Label("[bold cyan]Input Source[/bold cyan]")
                    yield Button("📂 Load Input File", id="inpgen-btn-load-input", variant="primary")
                    yield Button("📂 Load Structure (ASE)", id="inpgen-btn-load-structure", variant="default")
                    yield Label("ASE Format:", classes="inpgen-section-title")
                    yield Select(
                        [(name, value) for value, name in STRUCTURE_FORMATS],
                        value="auto",
                        id="inpgen-format-select"
                    )
                    
                    yield Rule()
                    yield Label("[bold cyan]Structure Tools[/bold cyan]")
                    yield Button("🎬 Generate Surface/Film", id="inpgen-btn-surface", variant="default")
                    yield Button("🔲 Create Supercell", id="inpgen-btn-supercell", variant="default")
                    yield Button("🧲 Set Magnetic Moments", id="inpgen-btn-magnetic", variant="default")
                
                with Vertical(id="inpgen-content"):
                    yield Label("[bold]Input Content[/bold] (edit namelist input below):")
                    yield TextArea(id="inpgen-input-area", language="text")
                    yield Label("[bold]Output Log:[/bold]")
                    yield Static("Ready to generate inp.xml", id="inpgen-output-area")
            
            # XML Editor container
            with Container(id="xml-editor-container"):
                with Container(id="tree-panel"):
                    with Container(id="search-container"):
                        yield Input(placeholder="Search schema...", id="search-input")
                    with VerticalScroll(id="search-results"):
                        pass  # Results will be added dynamically
                    yield Tree("Document", id="xml-tree")
                    with Horizontal(id="tree-buttons"):
                        yield Static("Select an element", id="no-add-elements")
                
                with VerticalScroll(id="editor-panel"):
                    yield Static("Select an element from the tree", id="no-selection")
            
            # K-Point Manager container (hidden by default)
            with Horizontal(id="kpoint-manager-container"):
                with Vertical(id="kpoint-sidebar"):
                    yield Label("[bold cyan]K-Point Lists[/bold cyan]")
                    yield Tree("K-Point Lists", id="kpoint-tree")
                    
                    yield Label("\n[bold cyan]Actions:[/bold cyan]")
                    yield Button("Refresh", id="kp-btn-refresh", variant="primary")
                    yield Button("Select as Active", id="kp-btn-select-list", variant="success")
                    yield Button("Create Mesh", id="kp-btn-create-mesh")
                    yield Button("Create Density", id="kp-btn-create-density")
                    yield Button("Create Path", id="kp-btn-create-path")
                    yield Button("Delete List", id="kp-btn-delete-list", variant="error")
                
                with Vertical(id="kpoint-content"):
                    yield Static(id="kpoint-info-panel")
                    
                    with Container(id="kpoint-visualization-panel"):
                        with Horizontal(id="kpoint-plot-controls"):
                            yield Label("Projection (in internal coordinates):")
                            yield Button("XY", id="kp-btn-proj-xy", variant="primary")
                            yield Button("XZ", id="kp-btn-proj-xz", variant="default")
                            yield Button("YZ", id="kp-btn-proj-yz", variant="default")
                        yield KPointPlot(id="kpoint-plot")
                    
                    yield DataTable(id="kpoint-table")
            
            # Job Generator container (hidden by default)
            with Horizontal(id="job-generator-container"):
                with VerticalScroll(id="job-sidebar"):
                    yield Static("[bold]Machine Configuration[/bold]", classes="job-form-section-title")
                    with Horizontal(classes="job-form-row"):
                        yield Label("Machine:")
                        yield Select(
                            self._get_job_machine_choices(),
                            value="(None)",
                            id="job-machine-select"
                        )
                    with Horizontal(classes="job-form-row"):
                        yield Label("Partition:")
                        yield Select(
                            [("(select machine first)", "")],
                            value="",
                            id="job-partition"
                        )
                    yield Static("No machine selected", id="job-machine-info")
                    
                    yield Rule()
                    
                    yield Static("[bold]Parallelization & Resources[/bold]", classes="job-form-section-title")
                    yield Static("Load machine config and FLEUR input to see suggestions", id="job-parallelization-info")
                    
                    with Horizontal(classes="job-form-row"):
                        yield Label("Job Name:")
                        yield Input(value="fleur_job", id="job-name")
                    with Horizontal(classes="job-form-row"):
                        yield Label("Account:")
                        yield Input(placeholder="(optional)", id="job-account")
                    with Horizontal(classes="job-form-row"):
                        yield Label("Nodes:")
                        yield Input(value="1", id="job-nodes", type="integer")
                    with Horizontal(classes="job-form-row"):
                        yield Label("Tasks:")
                        yield Input(value="1", id="job-ntasks", type="integer")
                    with Horizontal(classes="job-form-row"):
                        yield Label("CPUs/Task:")
                        yield Input(value="1", id="job-cpus-per-task", type="integer")
                    with Horizontal(classes="job-form-row"):
                        yield Label("Memory:")
                        yield Input(value="4G", id="job-memory", placeholder="e.g., 4G, 8000M")
                    with Horizontal(classes="job-form-row"):
                        yield Label("Time Limit:")
                        yield Input(value="01:00:00", id="job-time-limit", placeholder="HH:MM:SS")
                    with Horizontal(classes="job-form-row"):
                        yield Label("GPUs:")
                        yield Input(placeholder="e.g., 1, a100:2", id="job-gpus")
                    
                    with Horizontal(classes="job-button-bar"):
                        yield Button("Generate", variant="primary", id="job-btn-generate")
                        yield Button("Save", variant="success", id="job-btn-save")
                    
                    yield Rule()
                    
                    yield Static("[bold]FLEUR Input Analysis[/bold]", classes="job-form-section-title")
                    yield Static("No FLEUR input loaded", id="job-fleur-analysis")
                    
                    yield Rule()
                    
                    yield Static("[bold]Notifications[/bold]", classes="job-form-section-title")
                    with Horizontal(classes="job-form-row"):
                        yield Label("Email:")
                        yield Input(placeholder="user@example.com", id="job-mail-user")
                    with Horizontal(classes="job-form-row"):
                        yield Label("Notify on:")
                        yield Select(
                            [("END,FAIL", "END,FAIL"), ("BEGIN,END,FAIL", "BEGIN,END,FAIL"), 
                             ("ALL", "ALL"), ("NONE", "NONE")],
                            value="END,FAIL",
                            id="job-mail-type"
                        )
                    
                    yield Rule()
                    
                    yield Static("[bold]Modules[/bold]", classes="job-form-section-title")
                    yield TextArea(id="job-module-list", language=None)
                    yield Label("[dim]One module per line[/dim]")
                    
                    yield Rule()
                    
                    yield Static("[bold]Commands[/bold]", classes="job-form-section-title")
                    yield TextArea(id="job-command-list", language="bash")
                    yield Label("[dim]Shell commands to execute[/dim]")
                
                with Vertical(id="job-content"):
                    yield Static("[bold]SLURM Script Preview[/bold]", id="job-preview-title")
                    yield TextArea(id="job-preview-area", language="bash", read_only=True)
        
        yield Static("Ready", id="status-bar")
        yield Footer()
    
    def on_mount(self) -> None:
        """Called when app is mounted."""
        # Set initial view visibility
        self._apply_view_visibility()
        
        if self.input_file:
            self._load_file(Path(self.input_file))
            # Also pre-load for Job Generator mode
            if HAS_PYJOB:
                self._load_job_fleur_input(Path(self.input_file))
    
    def _apply_view_visibility(self):
        """Apply visibility based on current active view."""
        # Hide all containers first
        self.query_one("#xml-editor-container", Container).styles.display = "none"
        self.query_one("#kpoint-manager-container", Horizontal).styles.display = "none"
        self.query_one("#inpgen-container", Horizontal).styles.display = "none"
        self.query_one("#job-generator-container", Horizontal).styles.display = "none"
        
        # Remove active-mode from all menu buttons
        self.query_one("#menu-xml-editor", Button).remove_class("active-mode")
        self.query_one("#menu-kpoint-manager", Button).remove_class("active-mode")
        self.query_one("#menu-inpgen", Button).remove_class("active-mode")
        self.query_one("#menu-job-generator", Button).remove_class("active-mode")
        
        # Show the active container and highlight menu button
        if self._active_view == "xml_editor":
            self.query_one("#xml-editor-container", Container).styles.display = "block"
            self.query_one("#menu-xml-editor", Button).add_class("active-mode")
        elif self._active_view == "kpoint_manager":
            self.query_one("#kpoint-manager-container", Horizontal).styles.display = "block"
            self.query_one("#menu-kpoint-manager", Button).add_class("active-mode")
        elif self._active_view == "job_generator":
            self.query_one("#job-generator-container", Horizontal).styles.display = "block"
            self.query_one("#menu-job-generator", Button).add_class("active-mode")
        else:  # inpgen
            self.query_one("#inpgen-container", Horizontal).styles.display = "block"
            self.query_one("#menu-inpgen", Button).add_class("active-mode")
    
    # ==================== Mode Switching ====================
    
    @on(Button.Pressed, "#menu-xml-editor")
    def on_menu_xml_editor(self):
        """Switch to XML Editor mode."""
        self.action_switch_to_xml_editor()
    
    @on(Button.Pressed, "#menu-kpoint-manager")
    def on_menu_kpoint_manager(self):
        """Switch to K-Point Manager mode."""
        self.action_switch_to_kpoint_manager()
    
    @on(Button.Pressed, "#menu-inpgen")
    def on_menu_inpgen(self):
        """Switch to Inpgen mode."""
        self.action_switch_to_inpgen()
    
    @on(Button.Pressed, "#menu-job-generator")
    def on_menu_job_generator(self):
        """Switch to Job Generator mode."""
        self.action_switch_to_job_generator()
    
    def action_switch_to_xml_editor(self):
        """Switch to XML Editor mode."""
        if self._active_view == "xml_editor":
            return
        
        self._active_view = "xml_editor"
        self._apply_view_visibility()
        
        # Check if we have an inp.xml file
        if self.input_file and Path(self.input_file).exists():
            # Reload the file if tree is empty
            tree = self.query_one("#xml-tree", Tree)
            if not tree.root.children:
                self._load_file(Path(self.input_file))
            self._update_status("XML Editor mode")
        else:
            # No inp.xml file - show hint to use Inpgen mode
            self.query_one("#no-selection", Static).update(
                "[yellow]No inp.xml file loaded.[/yellow]\n\n"
                "Use [bold cyan]Inpgen mode (F3)[/bold cyan] to generate an inp.xml file first,\n"
                "or load an existing file."
            )
            self._update_status("XML Editor mode - No file loaded")
    
    def action_switch_to_kpoint_manager(self):
        """Switch to K-Point Manager mode."""
        if self._active_view == "kpoint_manager":
            return
        
        self._active_view = "kpoint_manager"
        self._apply_view_visibility()
        
        # Check if we have an inp.xml file
        if self.input_file and Path(self.input_file).exists():
            # Initialize k-point manager if not already done
            if self._kpoint_manager is None:
                self._load_kpoint_manager()
            self._refresh_kpoint_display()
            self._update_status("K-Point Manager mode")
        else:
            # No inp.xml file - show hint to use Inpgen mode
            info_panel = self.query_one("#kpoint-info-panel", Static)
            info_panel.update(
                "[yellow]No inp.xml file loaded.[/yellow]\n\n"
                "Use [bold cyan]Inpgen mode (F3)[/bold cyan] to generate an inp.xml file first,\n"
                "or load an existing file."
            )
            self._update_status("K-Point Manager mode - No file loaded")
    
    def action_switch_to_inpgen(self):
        """Switch to Inpgen mode."""
        if self._active_view == "inpgen":
            return
        
        self._active_view = "inpgen"
        self._apply_view_visibility()
        self._update_status("Inpgen mode")
    
    def action_switch_to_job_generator(self):
        """Switch to Job Generator mode."""
        if self._active_view == "job_generator":
            return
        
        self._active_view = "job_generator"
        self._apply_view_visibility()
        
        # Check if we have an inp.xml file (same logic as K-Point Manager)
        if self.input_file and Path(self.input_file).exists():
            # Load FLEUR analysis if not already done
            if self._job_fleur_analysis is None:
                self._load_job_fleur_input(Path(self.input_file))
            self._update_job_preview()
            self._update_status("Job Generator mode")
        else:
            # No inp.xml file - show hint to use Inpgen mode
            fleur_display = self.query_one("#job-fleur-analysis", Static)
            fleur_display.update(
                "[yellow]No inp.xml file loaded.[/yellow]\n\n"
                "Use [bold cyan]Inpgen mode (F3)[/bold cyan] to generate an inp.xml file first,\n"
                "or click 'Load inp.xml' to load an existing file."
            )
            self._update_job_preview()
            self._update_status("Job Generator mode - No file loaded")
    
    # ==================== Inpgen Methods ====================
    
    @on(Button.Pressed, "#inpgen-btn-nosym")
    def on_inpgen_toggle_nosym(self, event: Button.Pressed):
        """Toggle no symmetry option."""
        self._inpgen_nosym = not self._inpgen_nosym
        btn = event.button
        if self._inpgen_nosym:
            btn.label = "☑ No Symmetry"
            btn.variant = "success"
        else:
            btn.label = "☐ No Symmetry"
            btn.variant = "default"
    
    @on(Select.Changed, "#inpgen-profile-select")
    def on_inpgen_profile_changed(self, event: Select.Changed):
        """Handle profile selection change."""
        if event.value:
            self._inpgen_profile = str(event.value)
    
    @on(Select.Changed, "#inpgen-format-select")
    def on_inpgen_format_changed(self, event: Select.Changed):
        """Handle format selection change."""
        if event.value:
            self._inpgen_format = str(event.value)
    
    @on(Button.Pressed, "#inpgen-btn-load-input")
    def on_inpgen_load_input(self):
        """Load an input file."""
        def handle_result(path):
            if path:
                try:
                    content = path.read_text()
                    self.query_one("#inpgen-input-area", TextArea).text = content
                    self._update_inpgen_output(f"Loaded input file: {path.name}")
                except Exception as e:
                    self._update_inpgen_output(f"Error loading file: {e}")
        
        self.push_screen(FilePickerDialog(start_path=".", title="Select Input File"), handle_result)
    
    @on(Button.Pressed, "#inpgen-btn-load-structure")
    def on_inpgen_load_structure(self):
        """Load a structure file using ASE."""
        if not ASE_AVAILABLE:
            self.notify("ASE library not available. Install with: pip install ase", severity="error")
            return
        
        def handle_result(path):
            if path:
                try:
                    # Read structure using ASE
                    from ase.io import read as ase_read
                    fmt = None if self._inpgen_format == "auto" else self._inpgen_format
                    atoms = ase_read(str(path), format=fmt)
                    
                    # Convert to FLEUR input
                    fleur_input = ase_to_fleur_input(atoms, film=self._inpgen_film)
                    self.query_one("#inpgen-input-area", TextArea).text = fleur_input
                    
                    self._update_inpgen_output(
                        f"Loaded structure: {path.name}\n"
                        f"Formula: {atoms.get_chemical_formula()}\n"
                        f"Atoms: {len(atoms)}"
                    )
                except Exception as e:
                    self._update_inpgen_output(f"Error loading structure: {e}")
        
        self.push_screen(FilePickerDialog(start_path=".", title="Select Structure File"), handle_result)
    
    @on(Button.Pressed, "#inpgen-btn-generate")
    def on_inpgen_generate(self):
        """Generate inp.xml from the input."""
        input_text = self.query_one("#inpgen-input-area", TextArea).text
        if not input_text.strip():
            self.notify("No input content. Load a file or enter namelist input.", severity="warning")
            return
        
        # Check if inp.xml already exists
        existing_file = Path.cwd() / "inp.xml"
        if existing_file.exists():
            self._pending_generate_input = input_text
            self.notify(
                "inp.xml already exists! Press Generate again to overwrite.",
                severity="warning",
                timeout=5
            )
            self._update_inpgen_output(
                "[yellow]⚠ Warning: inp.xml already exists![/yellow]\n"
                "Press 'Generate inp.xml' again to overwrite the existing file."
            )
            # Set a flag to allow overwrite on next click
            if hasattr(self, '_overwrite_confirmed') and self._overwrite_confirmed:
                self._overwrite_confirmed = False
                self._do_generate(input_text)
            else:
                self._overwrite_confirmed = True
            return
        
        self._do_generate(input_text)
    
    def _do_generate(self, input_text: str):
        """Actually perform the inp.xml generation."""
        try:
            self._update_inpgen_output("Generating inp.xml...")
            
            # Use quiet=True to prevent console output, get messages via get_messages()
            inpgen = create_inpgen_interface(quiet=True)
            inpgen.make_inp(input_text, self._inpgen_profile, self._inpgen_nosym)
            inpgen_output = inpgen.get_messages()
            
            # Set the generated file as the active file for other modes
            generated_path = Path.cwd() / "inp.xml"
            if generated_path.exists():
                self.input_file = generated_path
                # Pre-load for XML Editor mode
                self._load_file(generated_path)
                # Reset k-point manager so it reloads on next switch
                self._kpoint_manager = None
                # Reset job generator analysis so it reloads on next switch
                self._job_fleur_analysis = None
            
            self._update_inpgen_output(
                f"✓ Successfully generated inp.xml\n"
                f"Profile: {self._inpgen_profile}\n"
                f"No symmetry: {self._inpgen_nosym}\n\n"
                f"[bold]Library output:[/bold]\n{inpgen_output}\n"
                f"File available in XML Editor, K-Point Manager, and Job Generator modes."
            )
            self.notify("inp.xml generated successfully!", severity="success")
            
        except Exception as e:
            self._update_inpgen_output(f"✗ Error generating inp.xml:\n{e}")
            self.notify(f"Error: {e}", severity="error")
    
    @on(Button.Pressed, "#inpgen-btn-open")
    def on_inpgen_open_generated(self):
        """Open the generated inp.xml file."""
        inp_path = Path("inp.xml")
        if inp_path.exists():
            self.input_file = inp_path
            self._load_file(inp_path)
            self.action_switch_to_xml_editor()
            self.notify("Opened inp.xml", severity="success")
        else:
            self.notify("inp.xml not found. Generate it first.", severity="warning")
    
    @on(Button.Pressed, "#inpgen-btn-surface")
    def on_inpgen_surface(self):
        """Generate a surface/film structure using ASE."""
        if not ASE_AVAILABLE:
            self.notify("ASE library required for surface generation", severity="error")
            return
        
        def handle_surface(result):
            if result is None:
                return
            
            try:
                element = result['element']
                surface_type = result['surface_type']
                layers = result['layers']
                vacuum = result['vacuum']
                lattice_const = result['lattice_const']
                size = result['size']
                miller = result.get('miller')
                
                # Generate surface
                slab = generate_surface(
                    element=element,
                    surface_type=surface_type,
                    layers=layers,
                    vacuum=vacuum,
                    lattice_const=lattice_const,
                    size=size,
                    orthogonal=True,
                    miller=miller
                )
                
                # Convert to FLEUR input (film mode for surfaces)
                fleur_input = ase_to_fleur_input(slab, film=True)
                self.query_one("#inpgen-input-area", TextArea).text = fleur_input
                
                miller_str = f"({miller[0]}{miller[1]}{miller[2]})" if miller else surface_type
                self._update_inpgen_output(
                    f"✓ Generated {element} {miller_str} surface\n"
                    f"Layers: {layers}\n"
                    f"Vacuum: {vacuum} Å\n"
                    f"Size: {size[0]}×{size[1]}\n"
                    f"Total atoms: {len(slab)}"
                )
                self.notify(f"Generated {element} surface with {len(slab)} atoms", severity="success")
                
            except Exception as e:
                self._update_inpgen_output(f"✗ Error generating surface:\n{e}")
                self.notify(f"Error: {e}", severity="error")
        
        self.push_screen(SurfaceDialog(), handle_surface)
    
    @on(Button.Pressed, "#inpgen-btn-supercell")
    def on_inpgen_supercell(self):
        """Create a supercell from the current input."""
        if not ASE_AVAILABLE:
            self.notify("ASE library required for supercell creation", severity="error")
            return
        
        input_text = self.query_one("#inpgen-input-area", TextArea).text
        if not input_text.strip():
            self.notify("Load a structure first", severity="warning")
            return
        
        # Try to parse the current input to get atom count
        atoms = parse_namelist_to_atoms(input_text)
        natoms = len(atoms) if atoms else 0
        
        # Extract magnetic moments to preserve them
        magnetic_moments = extract_magnetic_moments(input_text)
        
        def handle_supercell(result):
            if result is None:
                return
            
            nx, ny, nz = result
            
            if atoms is None:
                self.notify("Could not parse input. Load a structure file.", severity="error")
                return
            
            try:
                supercell_atoms = create_supercell(atoms, nx, ny, nz)
                # Convert back to FLEUR input format, preserving magnetic moments
                film = self._inpgen_film
                new_input = ase_to_fleur_input(
                    supercell_atoms, 
                    film=film, 
                    magnetic_moments=magnetic_moments
                )
                self.query_one("#inpgen-input-area", TextArea).text = new_input
                
                moment_info = ""
                if magnetic_moments:
                    moment_info = f"\nMagnetic moments preserved: {list(magnetic_moments.keys())}"
                
                self._update_inpgen_output(
                    f"✓ Created {nx}×{ny}×{nz} supercell\n"
                    f"Original: {len(atoms)} atoms\n"
                    f"Supercell: {len(supercell_atoms)} atoms{moment_info}"
                )
                self.notify(f"Created {nx}×{ny}×{nz} supercell with {len(supercell_atoms)} atoms", severity="success")
            except Exception as e:
                self._update_inpgen_output(f"✗ Error creating supercell:\n{e}")
                self.notify(f"Error: {e}", severity="error")
        
        self.push_screen(SupercellDialog(natoms), handle_supercell)
    
    @on(Button.Pressed, "#inpgen-btn-magnetic")
    def on_inpgen_magnetic(self):
        """Set magnetic moments for atoms of a specific element."""
        input_text = self.query_one("#inpgen-input-area", TextArea).text
        if not input_text.strip():
            self.notify("Load a structure first", severity="warning")
            return
        
        # Get elements in the structure
        elements = get_elements_in_input(input_text)
        if not elements:
            self.notify("No atoms found in input", severity="warning")
            return
        
        def handle_magnetic(result):
            if result is None:
                return
            
            element_z, moment_spec = result
            
            try:
                from .inpgen_gui import ATOMIC_SYMBOLS
                symbol = ATOMIC_SYMBOLS.get(element_z, f'Z{element_z}')
                
                new_input = add_magnetic_moments(input_text, element_z, moment_spec)
                self.query_one("#inpgen-input-area", TextArea).text = new_input
                
                self._update_inpgen_output(
                    f"✓ Set magnetic moments for {symbol} (Z={element_z})\n"
                    f"Moment specification: {moment_spec}"
                )
                self.notify(f"Set magnetic moment for {symbol}: {moment_spec}", severity="success")
            except Exception as e:
                self._update_inpgen_output(f"✗ Error setting magnetic moments:\n{e}")
                self.notify(f"Error: {e}", severity="error")
        
        self.push_screen(MagneticMomentDialog(elements), handle_magnetic)
    
    def _update_inpgen_output(self, text: str):
        """Update the inpgen output area."""
        self.query_one("#inpgen-output-area", Static).update(text)
    
    # ==================== K-Point Manager Methods ====================
    
    def _load_kpoint_manager(self):
        """Load the k-point manager for the current file."""
        if not self.input_file:
            return
        try:
            if Path(self.input_file).exists():
                self._kpoint_manager = InpXMLManager(str(self.input_file), quiet=True)
                self.notify(f"K-Point Manager: Loaded {self.input_file}", severity="information")
            else:
                self.notify(f"File not found: {self.input_file}", severity="warning")
        except Exception as e:
            self.notify(f"Error loading k-points: {e}", severity="error")

    def _save_kpoints_and_sync(self) -> None:
        """Save k-point manager state to disk and reload XMLDocument to keep views in sync."""
        if self._kpoint_manager is None:
            return
        self._kpoint_manager.save()
        # Reload the XML Editor model so it reflects the saved kpoint changes
        if self.document.file_path and self.document.file_path.exists():
            try:
                self.document.load(self.document.file_path)
            except Exception:
                pass
    
    def _refresh_kpoint_display(self):
        """Refresh the k-point display."""
        self._update_kpoint_info_panel()
        self._update_kpoint_tree()
        self._update_kpoint_table()
    
    def _update_kpoint_info_panel(self):
        """Update the k-point information panel."""
        info = self.query_one("#kpoint-info-panel", Static)
        
        if self._kpoint_manager is None:
            info.update("[error]No inp.xml loaded[/error]")
            return
        
        summary = self._kpoint_manager.get_summary()
        
        text = f"""[bold]File:[/bold] {summary['file']}
[bold]K-Point Lists:[/bold] {summary['num_lists']}
[bold]Selected List:[/bold] {summary['selected_list'] or 'None'}"""
        
        info.update(text)
    
    def _show_library_output(self, output: str):
        """Display library output in the k-point info panel or via notification.
        
        Args:
            output: The captured stdout from FleurInpgen library calls
        """
        if not output.strip():
            return
        
        # Update the info panel to include library output
        try:
            info = self.query_one("#kpoint-info-panel", Static)
            current_text = str(info.renderable) if info.renderable else ""
            
            # Add library output below current content
            new_text = f"{current_text}\n\n[bold]Library output:[/bold]\n{output.strip()}"
            info.update(new_text)
        except Exception:
            # Fallback: show as notification
            # Truncate if too long
            if len(output) > 200:
                output = output[:200] + "..."
            self.notify(f"Library: {output.strip()}", severity="information", timeout=5)
    
    def _update_kpoint_tree(self):
        """Update the k-point lists tree."""
        tree = self.query_one("#kpoint-tree", Tree)
        tree.clear()
        
        if self._kpoint_manager is None:
            return
        
        tree.root.expand()
        
        for name, kplist in self._kpoint_manager.kpoint_lists.items():
            selected = " ✓" if name == self._kpoint_manager.selected_list else ""
            node_label = f"{name} ({len(kplist.kpoints)}){selected}"
            node = tree.root.add(node_label, data=name)
            node.allow_expand = True
    
    @on(Tree.NodeSelected, "#kpoint-tree")
    def on_kpoint_tree_node_selected(self, event: Tree.NodeSelected):
        """Handle k-point tree node selection."""
        if event.node.data:
            self._selected_kpoint_list = event.node.data
            self._update_kpoint_table()
    
    def _update_kpoint_table(self):
        """Update the k-points table and visualization."""
        table = self.query_one("#kpoint-table", DataTable)
        table.clear(columns=True)
        
        if self._kpoint_manager is None or self._selected_kpoint_list is None:
            # Clear visualization
            plot = self.query_one("#kpoint-plot", KPointPlot)
            plot.set_kpoint_list(None)
            return
        
        kplist = self._kpoint_manager.get_kpoint_list(self._selected_kpoint_list)
        
        if kplist is None:
            return
        
        # Update visualization
        plot = self.query_one("#kpoint-plot", KPointPlot)
        plot.set_kpoint_list(kplist, self._current_projection)
        
        # Set up columns
        table.add_columns("Index", "X", "Y", "Z", "Weight", "Label")
        
        def format_value(val):
            if isinstance(val, float):
                return f"{val:.8f}"
            return str(val)
        
        # Add rows
        for i, kp in enumerate(kplist.kpoints):
            table.add_row(
                str(i),
                format_value(kp.x),
                format_value(kp.y),
                format_value(kp.z),
                format_value(kp.weight),
                kp.label
            )
    
    # K-Point projection buttons
    @on(Button.Pressed, "#kp-btn-proj-xy")
    def on_kpoint_projection_xy(self):
        self._current_projection = "xy"
        self._update_projection_buttons()
        self._update_kpoint_table()
    
    @on(Button.Pressed, "#kp-btn-proj-xz")
    def on_kpoint_projection_xz(self):
        self._current_projection = "xz"
        self._update_projection_buttons()
        self._update_kpoint_table()
    
    @on(Button.Pressed, "#kp-btn-proj-yz")
    def on_kpoint_projection_yz(self):
        self._current_projection = "yz"
        self._update_projection_buttons()
        self._update_kpoint_table()
    
    def _update_projection_buttons(self):
        """Update projection button variants."""
        self.query_one("#kp-btn-proj-xy", Button).variant = "primary" if self._current_projection == "xy" else "default"
        self.query_one("#kp-btn-proj-xz", Button).variant = "primary" if self._current_projection == "xz" else "default"
        self.query_one("#kp-btn-proj-yz", Button).variant = "primary" if self._current_projection == "yz" else "default"
    
    # K-Point action buttons
    @on(Button.Pressed, "#kp-btn-refresh")
    def on_kpoint_refresh(self):
        """Refresh k-point data from file."""
        self._load_kpoint_manager()
        self._refresh_kpoint_display()
    
    @on(Button.Pressed, "#kp-btn-select-list")
    def on_kpoint_select_list(self):
        """Select the current list as active."""
        if self._kpoint_manager is None:
            self.notify("No inp.xml loaded", severity="error")
            return
        
        if self._selected_kpoint_list is None:
            self.notify("Please select a k-point list first", severity="warning")
            return
        
        try:
            if self._kpoint_manager.select_kpoint_list(self._selected_kpoint_list):
                self.notify(f"Set '{self._selected_kpoint_list}' as active k-point list", severity="success")
                self._save_kpoints_and_sync()
                self._refresh_kpoint_display()
            else:
                self.notify(f"Could not select '{self._selected_kpoint_list}'", severity="error")
        except Exception as e:
            self.notify(f"Error selecting list: {e}", severity="error")
    
    @on(Button.Pressed, "#kp-btn-create-mesh")
    def on_kpoint_create_mesh(self):
        """Open create mesh dialog."""
        if self._kpoint_manager is None:
            self.notify("No inp.xml loaded", severity="error")
            return
        
        def handle_result(result):
            if result:
                name, nx, ny, nz, modifiers_dict, use_symmetry = result
                try:
                    modifiers = KPointModifiers(**modifiers_dict)
                    self._kpoint_manager.create_kpoints(
                        name=name,
                        mode=KPointMode.GRID,
                        modifiers=modifiers,
                        grid=(nx, ny, nz),
                        use_symmetry=use_symmetry
                    )
                    lib_output = self._kpoint_manager.get_last_messages()
                    mod_str = modifiers.to_string().replace('@', '+') if modifiers.to_string() else "standard"
                    self.notify(f"Created mesh '{name}' ({nx}x{ny}x{nz}) [{mod_str}]", severity="success")
                    self._show_library_output(lib_output)
                    self._save_kpoints_and_sync()
                    self._refresh_kpoint_display()
                except Exception as e:
                    self.notify(f"Error creating mesh: {e}", severity="error")
        
        self.push_screen(CreateMeshDialog(), handle_result)
    
    @on(Button.Pressed, "#kp-btn-create-density")
    def on_kpoint_create_density(self):
        """Open create density/number dialog."""
        if self._kpoint_manager is None:
            self.notify("No inp.xml loaded", severity="error")
            return
        
        def handle_result(result):
            if result:
                name, mode, value, modifiers_dict, use_symmetry = result
                try:
                    modifiers = KPointModifiers(**modifiers_dict)
                    
                    if mode == "density":
                        self._kpoint_manager.create_kpoints(
                            name=name,
                            mode=KPointMode.DENSITY,
                            modifiers=modifiers,
                            density=value,
                            use_symmetry=use_symmetry
                        )
                        mod_str = modifiers.to_string().replace('@', '+') if modifiers.to_string() else "standard"
                        self.notify(f"Created density-based '{name}' (density={value}) [{mod_str}]", 
                                   severity="success")
                    else:
                        self._kpoint_manager.create_kpoints(
                            name=name,
                            mode=KPointMode.NUMBER,
                            modifiers=modifiers,
                            num_kpoints=value,
                            use_symmetry=use_symmetry
                        )
                        mod_str = modifiers.to_string().replace('@', '+') if modifiers.to_string() else "standard"
                        self.notify(f"Created number-based '{name}' (nk={value}) [{mod_str}]", 
                                   severity="success")
                    
                    lib_output = self._kpoint_manager.get_last_messages()
                    self._show_library_output(lib_output)
                    self._save_kpoints_and_sync()
                    self._refresh_kpoint_display()
                except Exception as e:
                    self.notify(f"Error creating k-points: {e}", severity="error")
        
        self.push_screen(CreateDensityDialog(), handle_result)
    
    @on(Button.Pressed, "#kp-btn-create-path")
    def on_kpoint_create_path(self):
        """Open create path dialog."""
        if self._kpoint_manager is None:
            self.notify("No inp.xml loaded", severity="error")
            return
        
        def handle_result(result):
            if result:
                name, lattice, path_str, npoints = result
                try:
                    points_dict = STANDARD_KPOINTS.get(lattice, STANDARD_KPOINTS['fcc'])
                    
                    # Parse path
                    if '-' in path_str:
                        path_labels = path_str.split('-')
                    else:
                        path_labels = list(path_str)
                    
                    # Build special points list
                    special_points = []
                    for label in path_labels:
                        label = label.strip()
                        if label not in points_dict:
                            self.notify(f"Unknown k-point: {label}", severity="error")
                            return
                        special_points.append((label, points_dict[label]))
                    
                    # Create path
                    self._kpoint_manager.create_kpoint_path(name, special_points, npoints)
                    lib_output = self._kpoint_manager.get_last_messages()
                    
                    path_display = '-'.join(path_labels)
                    self.notify(f"Created path '{name}': {path_display} ({npoints} points)", 
                               severity="success")
                    self._show_library_output(lib_output)
                    self._save_kpoints_and_sync()
                    self._refresh_kpoint_display()
                    
                except Exception as e:
                    self.notify(f"Error creating path: {e}", severity="error")
        
        self.push_screen(CreatePathDialog(), handle_result)
    
    @on(Button.Pressed, "#kp-btn-delete-list")
    def on_kpoint_delete_list(self):
        """Delete selected k-point list."""
        if self._kpoint_manager is None or self._selected_kpoint_list is None:
            self.notify("No list selected", severity="warning")
            return
        
        # Prevent deletion of the active k-point list
        if self._selected_kpoint_list == self._kpoint_manager.get_selected_list():
            self.notify(f"Cannot delete '{self._selected_kpoint_list}': it is the active k-point list", severity="error")
            return
        
        try:
            self._kpoint_manager.remove_kpoint_list(self._selected_kpoint_list)
            self.notify(f"Deleted list '{self._selected_kpoint_list}'", severity="success")
            self._selected_kpoint_list = None
            self._save_kpoints_and_sync()
            self._refresh_kpoint_display()
        except Exception as e:
            self.notify(f"Error deleting list: {e}", severity="error")
    
    # ==================== XML Editor Methods ====================

    def _set_schema(self, schema_path: Path) -> bool:
        """Load and set schema for XML editor/document linking."""
        try:
            self.schema = XSDParser(schema_path)
            self.document.schema = self.schema
            self.schema_path = schema_path
            return True
        except Exception as exc:
            self.notify(f"Failed to load schema '{schema_path}': {exc}", severity="warning")
            return False

    def _schema_candidates_for_input(self, xml_path: Path) -> List[Path]:
        """Return potential schema paths in priority order."""
        candidates: List[Path] = []

        if self.schema_path:
            candidates.append(Path(self.schema_path))

        candidates.extend([
            xml_path.parent / "FleurInputSchema.xsd",
            Path.cwd() / "FleurInputSchema.xsd",
            Path.cwd() / "build" / "FleurInputSchema.xsd",
        ])

        unique: List[Path] = []
        seen: set[str] = set()
        for candidate in candidates:
            key = str(candidate.resolve(strict=False))
            if key not in seen:
                seen.add(key)
                unique.append(candidate)
        return unique

    def _ensure_schema_for_xml_editor(self, xml_path: Path):
        """Ensure schema is loaded; generate it for inp.xml if missing."""
        if self.schema is not None:
            return

        for candidate in self._schema_candidates_for_input(xml_path):
            if candidate.exists() and self._set_schema(candidate):
                return

        if xml_path.name != "inp.xml":
            return

        try:
            inpgen = create_inpgen_interface(quiet=True)
            original_cwd = Path.cwd()
            os.chdir(xml_path.parent)
            try:
                inpgen.dropxmlschema()
            finally:
                os.chdir(original_cwd)
        except Exception as exc:
            self.notify(
                f"No FleurInputSchema.xsd found and schema generation failed: {exc}",
                severity="warning",
            )
            return

        generated_schema = xml_path.parent / "FleurInputSchema.xsd"
        if generated_schema.exists() and self._set_schema(generated_schema):
            self.notify(f"Generated schema: {generated_schema}", severity="information")
        else:
            self.notify(
                "Schema generation completed, but FleurInputSchema.xsd was not found.",
                severity="warning",
            )
    
    def _load_file(self, path: Path):
        """Load an XML file."""
        try:
            self._ensure_schema_for_xml_editor(path)
            self.document.load(path)
            self._populate_tree()
            self._update_status(f"Loaded: {path.name}")
            self.title = f"FLEURiste - {path.name}"
        except Exception as e:
            self._update_status(f"Error: {e}")
    
    def _populate_tree(self):
        """Populate the tree from the document."""
        tree = self.query_one("#xml-tree", Tree)
        tree.clear()
        
        # Remember current node path if any
        current_path = None
        if self._current_node:
            current_path = self._current_node.get_path_string()
        
        self._node_map.clear()
        self._current_node = None
        
        if self.document.root:
            self._add_node_to_tree(self.document.root, tree.root)
            tree.root.expand()
            
            # Try to restore selection by path
            if current_path:
                for node_id, node in self._node_map.items():
                    if node.get_path_string() == current_path:
                        self._current_node = node
                        break
    
    def _add_node_to_tree(self, xml_node: XMLNode, tree_node):
        """Recursively add XML nodes to the tree."""
        display_name = xml_node.get_display_name()
        node_id = str(xml_node.uid)
        
        child_tree_node = tree_node.add(display_name, data=node_id)
        self._node_map[node_id] = xml_node
        
        for child_xml in xml_node.children:
            self._add_node_to_tree(child_xml, child_tree_node)
    
    def _update_status(self, message: str):
        """Update the status bar."""
        self.query_one("#status-bar", Static).update(message)
    
    def on_input_submitted(self, event: Input.Submitted) -> None:
        """Handle search input submission (Enter key)."""
        if event.input.id == "search-input":
            self._perform_search(event.value)
    
    def _perform_search(self, query: str):
        """Perform search and display results."""
        results_container = self.query_one("#search-results", VerticalScroll)
        
        # Increment search version for unique IDs
        self._search_version += 1
        
        # Clear previous results
        for child in list(results_container.children):
            child.remove()
        
        if not query or len(query) < 2:
            results_container.remove_class("visible")
            return
        
        if not self.schema:
            results_container.mount(Static("No schema loaded", classes="search-result"))
            results_container.add_class("visible")
            return
        
        # Perform search
        results = self.schema.search(query)
        
        if not results:
            results_container.mount(Static(f"No results for '{query}'", classes="search-result"))
            results_container.add_class("visible")
            return
        
        # Display results
        for i, result in enumerate(results):
            result_id = f"search-result-v{self._search_version}-{i}"
            type_label = {
                "element": "[elem]",
                "attribute": "[attr]",
                "documentation": "[doc]",
                "attribute_doc": "[attr doc]"
            }.get(result.match_type, "")
            
            # Check if path exists in current document
            search_path = result.path
            if "/@" in search_path:
                search_path = search_path.split("/@")[0]
            
            path_exists = any(
                node.get_path_string() == search_path 
                for node in self._node_map.values()
            )
            
            # Style based on whether path exists
            if path_exists:
                label = f"{type_label} {result.path}"
                variant = "default"
            else:
                label = f"[dim]{type_label} {result.path}[/dim]"
                variant = "default"
            
            result_widget = Button(
                label,
                id=result_id,
                variant=variant,
                classes="search-result"
            )
            # Store the path and context as data
            result_widget._search_path = result.path
            result_widget._search_context = result.context
            result_widget._path_exists = path_exists
            results_container.mount(result_widget)
        
        results_container.add_class("visible")
        self._update_status(f"Found {len(results)} result(s)")
    
    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses for adding elements and search results."""
        button_id = event.button.id
        
        # Handle search result clicks
        if button_id and button_id.startswith("search-result-"):
            if hasattr(event.button, "_search_path"):
                # Check if path exists in document
                if hasattr(event.button, "_path_exists") and not event.button._path_exists:
                    # Path doesn't exist - show documentation instead
                    self._show_schema_documentation(
                        event.button._search_path,
                        getattr(event.button, "_search_context", "")
                    )
                else:
                    self._navigate_to_path(event.button._search_path)
            return
        
        # Handle add element buttons
        if button_id and button_id.startswith("add-elem-"):
            element_name = button_id.replace("add-elem-", "")
            self._add_element_by_name(element_name)
    
    def _show_schema_documentation(self, path: str, context: str):
        """Show documentation for a schema element that doesn't exist in the document."""
        self._current_node = None
        self._editor_version += 1
        
        panel = self.query_one("#editor-panel", VerticalScroll)
        
        # Remove all existing children
        for child in list(panel.children):
            child.remove()
        
        # Extract element name from path
        element_name = path.split("/")[-1]
        if element_name.startswith("@"):
            element_name = element_name[1:]  # Remove @ for attributes
        
        # Title
        panel.mount(Static(f"[bold]{element_name}[/bold]", classes="node-title"))
        panel.mount(Static("[italic](not in current document)[/italic]", classes="node-path"))
        panel.mount(Rule())
        
        # Documentation
        if context:
            panel.mount(Static("[bold]Documentation[/bold]", classes="section-title"))
            panel.mount(Static(context, classes="doc-text"))
            panel.mount(Rule())
        
        # Path info
        panel.mount(Static("[bold]Schema Path:[/bold]", classes="section-title"))
        panel.mount(Static(path, classes="node-path"))
        
        self._update_status(f"Showing schema info for {path}")

    def _navigate_to_path(self, path: str):
        """Navigate to a node by its XML path."""
        # Remove attribute part if present (e.g., /@attrName)
        if "/@" in path:
            path = path.split("/@")[0]
        
        # Find the node with this path
        for node_id, node in self._node_map.items():
            if node.get_path_string() == path:
                # Select this node in the tree
                tree = self.query_one("#xml-tree", Tree)
                # Find the tree node and select it
                self._select_tree_node_by_id(tree.root, node_id)
                
                # Update editor
                self._current_node = node
                self._editor_version += 1
                self._update_add_button()
                self.call_after_refresh(self._rebuild_editor)
                
                # Hide search results
                self.query_one("#search-results", VerticalScroll).remove_class("visible")
                self._update_status(f"Navigated to {path}")
                return
        
        self._update_status(f"Path not found in document: {path}")
    
    def _select_tree_node_by_id(self, tree_node, target_id: str) -> bool:
        """Recursively find and select a tree node by its data ID."""
        if tree_node.data == target_id:
            tree_node.expand()
            tree = self.query_one("#xml-tree", Tree)
            tree.select_node(tree_node)
            return True
        
        for child in tree_node.children:
            if self._select_tree_node_by_id(child, target_id):
                tree_node.expand()
                return True
        
        return False

    def on_tree_node_highlighted(self, event: Tree.NodeHighlighted) -> None:
        """Handle tree cursor movement."""
        # Only handle events from our XML tree, not from DirectoryTree in dialogs
        if event.control.id != "xml-tree":
            return
        self._schedule_editor_update(event.node.data)
    
    def on_tree_node_selected(self, event: Tree.NodeSelected) -> None:
        """Handle tree node selection (Enter key)."""
        # Only handle events from our XML tree, not from DirectoryTree in dialogs
        if event.control.id != "xml-tree":
            return
        # Don't update on select - highlighted already handles it
        pass
    
    def _schedule_editor_update(self, node_id: Optional[str]):
        """Schedule an editor update, avoiding duplicates."""
        if not node_id or node_id not in self._node_map:
            return
        
        # Skip if same node
        node = self._node_map[node_id]
        if self._current_node and self._current_node.uid == node.uid:
            return
        
        self._current_node = node
        self._editor_version += 1  # Increment for unique IDs
        self._update_add_button()
        self.call_after_refresh(self._rebuild_editor)
    
    def _update_add_button(self):
        """Update the Add buttons to show available child elements."""
        try:
            buttons_container = self.query_one("#tree-buttons", Horizontal)
            
            # Remove all existing children
            for child in list(buttons_container.children):
                child.remove()
            
            if not self._current_node:
                buttons_container.mount(Static("Select an element", id="no-add-elements"))
                return
            
            # Check if schema is linked
            if not self._current_node.schema_element:
                buttons_container.mount(Static("No schema info available", id="no-add-elements"))
                return
            
            # Get addable children (respecting maxOccurs constraints)
            addable = self.document.get_addable_children(self._current_node)
            
            if not addable:
                buttons_container.mount(Static("No elements can be added", id="no-add-elements"))
                return
            
            # Add header label
            buttons_container.mount(Static("Add XML elements:", id="add-elements-label"))
            
            # Create a button for each addable element
            for schema_elem in addable:
                btn_id = f"add-elem-{schema_elem.name}"
                label = schema_elem.name
                btn = Button(label, id=btn_id, variant="primary")
                buttons_container.mount(btn)
                
        except Exception as e:
            self._update_status(f"Button error: {e}")
            pass
    
    def _rebuild_editor(self):
        """Rebuild the editor panel for the current node."""
        if not self._current_node:
            return
            
        node = self._current_node
        ver = self._editor_version  # Use version for unique IDs
        
        # Get the editor panel and rebuild its contents
        panel = self.query_one("#editor-panel", VerticalScroll)
        
        # Remove all existing children
        for child in list(panel.children):
            child.remove()
        
        # Title - use display_tag to resolve namespace prefixes
        panel.mount(Static(f"[bold]{node.display_tag}[/bold]", classes="node-title"))
        
        # Get schema info for later use
        schema_elem = node.schema_element
        
        # Show documentation first if available
        if schema_elem and schema_elem.documentation:
            panel.mount(Static("[bold]Documentation[/bold]", classes="section-title"))
            panel.mount(Static(schema_elem.documentation, classes="doc-text"))
            panel.mount(Rule())
        
        # Attributes
        has_attrs = bool(node.attributes) or (schema_elem and schema_elem.attributes)
        
        if has_attrs:
            panel.mount(Static("[bold]Attributes[/bold]", classes="section-title"))
            
            shown_attrs = set()
            
            # Schema-defined attributes first
            if schema_elem and schema_elem.attributes:
                for attr_name, attr_info in schema_elem.attributes.items():
                    current_val = node.attributes.get(attr_name, attr_info.default or "")
                    self._mount_attribute_row(panel, attr_name, current_val, attr_info, ver)
                    shown_attrs.add(attr_name)
            
            # Extra attributes not in schema
            for attr_name, attr_val in node.attributes.items():
                if attr_name not in shown_attrs:
                    self._mount_attribute_row(panel, attr_name, attr_val, None, ver)
        
        # Text content - show if has text OR schema says it's simple content
        has_text = node.text and node.text.strip()
        is_simple = schema_elem and schema_elem.is_simple_content
        
        if has_text or is_simple:
            panel.mount(Rule())
            panel.mount(Static("[bold]Text Content[/bold]", classes="section-title"))
            placeholder = ""
            if is_simple and schema_elem.simple_content_type:
                placeholder = f"Type: {schema_elem.simple_content_type}"
            panel.mount(Input(
                value=node.text or "", 
                id=f"text-content-{ver}",
                placeholder=placeholder
            ))
        
        # Schema info
        if schema_elem:
            panel.mount(Rule())
            panel.mount(Static("[bold]Schema Info[/bold]", classes="section-title"))
            
            info = []
            info.append(f"• {'Optional' if schema_elem.is_optional else 'Required'} element")
            if schema_elem.is_unbounded:
                info.append("• Multiple occurrences allowed")
            if schema_elem.children:
                names = [c.name for c in schema_elem.children[:5]]
                info.append(f"• Children: {', '.join(names)}{'...' if len(schema_elem.children) > 5 else ''}")
            
            panel.mount(Static("\n".join(info), classes="schema-info"))
        else:
            # No schema for this element (e.g., XInclude)
            panel.mount(Rule())
            panel.mount(Static("[bold]Element Info[/bold]", classes="section-title"))
            info = ["• Element not defined in schema (e.g., XInclude directive)"]
            if node.children:
                info.append(f"• Has {len(node.children)} child element(s)")
            panel.mount(Static("\n".join(info), classes="schema-info"))
        
        # XML Path at the end
        panel.mount(Rule())
        panel.mount(Static(f"[bold]Path:[/bold] {node.get_path_string()}", classes="node-path"))
    
    def _mount_attribute_row(self, container, attr_name: str, value: str, attr_info: Optional[AttributeInfo], ver: int):
        """Mount an attribute editor row."""
        required = "*" if attr_info and attr_info.is_required else ""
        
        row = Horizontal(classes="attr-row")
        label = Label(f"{attr_name}{required}:", classes="attr-label")
        
        # Use version in ID to ensure uniqueness
        widget_id = f"attr-{ver}-{attr_name}"
        
        # Get placeholder text - use documentation or type description
        placeholder = ""
        if attr_info:
            if attr_info.documentation:
                placeholder = attr_info.documentation
            elif self.schema:
                placeholder = self.schema.get_type_description(attr_info.type_name)
        
        # Determine input widget type
        if attr_info and self.schema:
            type_name = attr_info.type_name
            
            if self.schema.is_boolean_type(type_name):
                input_widget = Select(
                    [("T", "T"), ("F", "F")],
                    value=value.upper() if value else "F",
                    id=widget_id,
                    classes="attr-input"
                )
            elif self.schema.get_enum_values(type_name):
                values = self.schema.get_enum_values(type_name)
                options = [(v, v) for v in values]
                input_widget = Select(
                    options,
                    value=value if value in values else (values[0] if values else ""),
                    id=widget_id,
                    classes="attr-input"
                )
            else:
                input_widget = Input(
                    value=value,
                    id=widget_id,
                    placeholder=placeholder,
                    classes="attr-input"
                )
        else:
            input_widget = Input(value=value, id=widget_id, classes="attr-input")
        
        # Mount the row with its children
        container.mount(row)
        row.mount(label)
        row.mount(input_widget)
        
        # Show attribute documentation below the input if available
        if attr_info and attr_info.documentation:
            container.mount(Static(f"↳ {attr_info.documentation}", classes="attr-doc"))
    
    def _parse_attr_id(self, widget_id: str) -> Optional[str]:
        """Parse attribute name from widget ID like 'attr-123-attrName'."""
        if not widget_id or not widget_id.startswith("attr-"):
            return None
        # Format: attr-{version}-{attr_name}
        parts = widget_id.split("-", 2)  # Split into at most 3 parts
        if len(parts) >= 3:
            return parts[2]  # Return the attribute name
        return None
    
    def on_input_changed(self, event: Input.Changed) -> None:
        """Handle input field changes."""
        if not self._current_node:
            return
        
        input_id = event.input.id
        if input_id and input_id.startswith("text-content-"):
            self._current_node.text = event.value
            self.document.mark_modified()
            self._update_status("Modified")
        else:
            attr_name = self._parse_attr_id(input_id)
            if attr_name:
                self._current_node.attributes[attr_name] = event.value
                self.document.mark_modified()
                self._update_status("Modified")
    
    def on_select_changed(self, event: Select.Changed) -> None:
        """Handle select field changes."""
        if not self._current_node or event.value is None:
            return
        
        attr_name = self._parse_attr_id(event.select.id)
        if attr_name:
            self._current_node.attributes[attr_name] = str(event.value)
            self.document.mark_modified()
            self._update_status("Modified")
    
    def action_save(self) -> None:
        """Save the document."""
        if self._active_view == "kpoint_manager":
            if self._kpoint_manager is None:
                self._update_status("No k-points loaded")
                return
            try:
                self._save_kpoints_and_sync()
                self._update_status(f"Saved: {self._kpoint_manager.xml_path.name}")
            except Exception as e:
                self._update_status(f"Error: {e}")
            return
        if not self.document.root:
            self._update_status("No document to save")
            return
        try:
            self.document.save()
            self._update_status(f"Saved: {self.document.file_path.name}")
        except Exception as e:
            self._update_status(f"Error: {e}")
    
    def action_undo(self) -> None:
        if self.document.undo():
            self._populate_tree()
            self._update_status("Undo")
    
    def action_redo(self) -> None:
        if self.document.redo():
            self._populate_tree()
            self._update_status("Redo")
    
    def action_refresh(self) -> None:
        self._populate_tree()
        self._update_status("Refreshed")
    
    def action_focus_search(self) -> None:
        """Focus the search input."""
        search_input = self.query_one("#search-input", Input)
        search_input.focus()
    
    def action_clear_search(self) -> None:
        """Clear search results and input."""
        search_input = self.query_one("#search-input", Input)
        search_input.value = ""
        results_container = self.query_one("#search-results", VerticalScroll)
        for child in list(results_container.children):
            child.remove()
        results_container.remove_class("visible")
    
    def action_delete_element(self) -> None:
        """Delete the currently selected element."""
        if not self._current_node:
            self._update_status("No element selected")
            return
        
        # Cannot delete root element
        if self._current_node.parent is None:
            self._update_status("Cannot delete root element")
            return
        
        parent = self._current_node.parent
        deleted_tag = self._current_node._get_readable_tag()
        
        # Save state for undo
        self.document.save_undo_state()
        
        # Remove the node from its parent
        parent.remove_child(self._current_node)
        self.document.mark_modified()
        
        # Store parent ID for selection after tree refresh
        parent_id = str(parent.uid)
        
        # Refresh the tree
        self._populate_tree()
        
        # Select the parent node
        tree = self.query_one("#xml-tree", Tree)
        self._select_tree_node_by_id(tree.root, parent_id)
        
        # Update the current node and editor
        if parent_id in self._node_map:
            self._current_node = self._node_map[parent_id]
            self._rebuild_editor()
            self._update_add_button()
        
        self._update_status(f"Deleted: {deleted_tag}")

    def _add_element_by_name(self, element_name: str) -> None:
        """Add a child element by name."""
        if not self._current_node:
            self._update_status("Select an element first")
            return
        
        # Get addable children and find the matching schema element
        addable = self.document.get_addable_children(self._current_node)
        for schema_elem in addable:
            if schema_elem.name == element_name:
                self.document.save_undo_state()
                new_node = self.document.create_element_from_schema(schema_elem)
                self._current_node.add_child(new_node)
                self.document.mark_modified()
                
                # Store the new node's ID for selection after tree refresh
                new_node_id = str(new_node.uid)
                
                # Refresh the tree
                self._populate_tree()
                
                # Select the newly added node in the tree
                tree = self.query_one("#xml-tree", Tree)
                self._select_tree_node_by_id(tree.root, new_node_id)
                
                # Update the current node and editor
                if new_node_id in self._node_map:
                    self._current_node = self._node_map[new_node_id]
                    self._rebuild_editor()
                    self._update_add_button()
                
                self._update_status(f"Added: {element_name}")
                return
        
        self._update_status(f"Cannot add: {element_name}")

    # ==================== Job Generator Methods ====================
    
    def _get_job_machine_choices(self):
        """Build machine selection dropdown choices."""
        choices = [("(None)", "(None)")]
        for name in sorted(self._job_available_machines.keys()):
            choices.append((name, name))
        choices.append(("(Load from file...)", "(Load from file...)"))
        return choices
    
    def _get_job_config(self) -> "SlurmJobConfig":
        """Build SlurmJobConfig from form inputs."""
        if not HAS_PYJOB:
            return None
        
        # Parse modules (one per line)
        modules_text = self.query_one("#job-module-list", TextArea).text
        modules = [m.strip() for m in modules_text.split("\n") if m.strip()]
        
        # Parse commands (one per line)
        commands_text = self.query_one("#job-command-list", TextArea).text
        commands = [c for c in commands_text.split("\n") if c.strip()]
        
        # Get form values
        job_name = self.query_one("#job-name", Input).value or "fleur_job"
        account = self.query_one("#job-account", Input).value or None
        partition_select = self.query_one("#job-partition", Select)
        partition = partition_select.value if partition_select.value else None
        
        nodes = int(self.query_one("#job-nodes", Input).value or "1")
        ntasks = int(self.query_one("#job-ntasks", Input).value or "1")
        cpus = int(self.query_one("#job-cpus-per-task", Input).value or "1")
        memory = self.query_one("#job-memory", Input).value or "4G"
        time_limit = self.query_one("#job-time-limit", Input).value or "01:00:00"
        gpus = self.query_one("#job-gpus", Input).value or None
        
        mail_user = self.query_one("#job-mail-user", Input).value or None
        mail_type_select = self.query_one("#job-mail-type", Select)
        mail_type = mail_type_select.value if mail_user else None
        
        # Get gpu_syntax and use_mem_option from machine config
        gpu_syntax = "gpus"
        use_mem_option = True
        if self._job_machine:
            gpu_syntax = self._job_machine.get_effective_value("gpu_syntax", partition) or "gpus"
            use_mem = self._job_machine.get_effective_value("use_mem_option", partition)
            use_mem_option = use_mem if use_mem is not None else True
        
        return SlurmJobConfig(
            job_name=job_name,
            nodes=nodes,
            ntasks=ntasks,
            cpus_per_task=cpus,
            memory=memory,
            time=time_limit,
            partition=partition,
            account=account,
            gpus=gpus,
            gpu_syntax=gpu_syntax,
            use_mem_option=use_mem_option,
            mail_user=mail_user,
            mail_type=mail_type,
            modules=modules,
            commands=commands,
        )
    
    def _update_job_preview(self):
        """Generate and update the job script preview."""
        if not HAS_PYJOB:
            return
        try:
            config = self._get_job_config()
            generator = SlurmJobGenerator(config, self._job_machine)
            script = generator.generate()
            
            preview = self.query_one("#job-preview-area", TextArea)
            preview.load_text(script)
            
            self._update_status(f"Preview updated | Job: {config.job_name}")
        except Exception as e:
            self._update_status(f"Error: {e}")
    
    def _update_job_machine_info(self):
        """Update machine info display."""
        if not HAS_PYJOB:
            return
        info_panel = self.query_one("#job-machine-info", Static)
        
        if self._job_machine:
            # Update partition dropdown
            partition_select = self.query_one("#job-partition", Select)
            partition_choices = []
            for p in self._job_machine.partitions:
                partition_choices.append((p.name, p.name))
            if partition_choices:
                partition_select.set_options(partition_choices)
                partition_select.value = partition_choices[0][1]
            else:
                partition_select.set_options([("(no partitions)", "")])
            
            # Suggest modules from selected/default partition
            partition_value = partition_select.value if partition_select.value else None
            self._update_job_module_suggestion(partition_value)
            
            # Update command suggestion
            self._update_job_command_suggestion()
            
            # Update partition info
            self._update_job_partition_info()
        else:
            machines_dir = get_machines_directory() if HAS_PYJOB else Path(".")
            if not self._job_available_machines:
                info_panel.update(f"No machine configs found in:\n{machines_dir}")
            else:
                info_panel.update("No machine selected")
            partition_select = self.query_one("#job-partition", Select)
            partition_select.set_options([("(select machine first)", "")])
    
    def _update_job_partition_info(self):
        """Update partition info display."""
        if not HAS_PYJOB or self._job_machine is None:
            return
        
        info_panel = self.query_one("#job-machine-info", Static)
        partition_select = self.query_one("#job-partition", Select)
        partition = partition_select.value if partition_select.value else None
        
        cores = self._job_machine.get_effective_value("cores_per_node", partition)
        gpus = self._job_machine.get_effective_value("gpus_per_node", partition) or 0
        memory = self._job_machine.get_effective_value("memory_per_node_gb", partition)
        max_runtime = self._job_machine.get_effective_value("max_runtime", partition)
        
        partition_name = partition if partition else "(none)"
        info = f"""[bold]{self._job_machine.name}[/bold] - [cyan]{partition_name}[/cyan]
[dim]Cores:[/dim] [green]{cores}[/green]/node  [dim]GPUs:[/dim] [green]{gpus}[/green]/node  [dim]Memory:[/dim] [green]{memory}[/green]GB"""
        
        if max_runtime:
            info += f"\n[dim]Max runtime:[/dim] {max_runtime}"
        
        info_panel.update(info)

        self._update_job_module_suggestion(partition)

        self._update_job_command_suggestion()

    def _update_job_module_suggestion(self, partition):
        """Update module textarea with machine/partition modules without overwriting user edits."""
        if not HAS_PYJOB or self._job_machine is None:
            return

        suggested_modules = self._job_machine.get_effective_value("modules_needed", partition) or []
        suggested_text = "\n".join(suggested_modules[:3]) if suggested_modules else ""

        modules_area = self.query_one("#job-module-list", TextArea)
        current_text = modules_area.text.strip()
        if (not current_text) or (current_text == self._job_last_machine_modules):
            modules_area.load_text(suggested_text)
            self._job_last_machine_modules = suggested_text
    
    def _update_job_command_suggestion(self):
        """Update command textarea with machine/partition command."""
        if not HAS_PYJOB or self._job_machine is None:
            return
        
        partition_select = self.query_one("#job-partition", Select)
        partition = partition_select.value if partition_select.value else None
        command = self._job_machine.get_effective_value("command", partition)
        
        if command:
            command_area = self.query_one("#job-command-list", TextArea)
            current_text = command_area.text.strip()
            if not current_text or current_text == self._job_last_machine_command:
                command_area.load_text(command)
                self._job_last_machine_command = command
    
    @on(Select.Changed, "#job-machine-select")
    def on_job_machine_changed(self, event: Select.Changed):
        """Handle machine selection change."""
        if not HAS_PYJOB:
            return
        value = event.value
        
        if value == "(None)":
            self._job_machine = None
        elif value == "(Load from file...)":
            self._action_load_job_machine()
            return
        elif value in self._job_available_machines:
            self._job_machine = self._job_available_machines[value]
            self._update_status(f"Loaded machine: {value}")
        
        self._update_job_machine_info()
        self._update_job_parallelization()
        self._update_job_preview()
    
    @on(Select.Changed, "#job-partition")
    def on_job_partition_changed(self, event: Select.Changed):
        """Handle partition selection change."""
        self._update_job_partition_info()
        self._update_job_parallelization()
        self._update_job_preview()
    
    @on(Input.Changed, "#job-nodes")
    def on_job_nodes_changed(self, event: Input.Changed):
        """Handle nodes input change."""
        if self._job_updating_from_parallelization:
            return
        self._update_job_parallelization()
        self._update_job_preview()
    
    @on(Input.Changed, "#job-ntasks")
    def on_job_ntasks_changed(self, event: Input.Changed):
        """Handle tasks input change."""
        if self._job_updating_from_parallelization:
            return
        self._update_job_cpus_per_task()
        self._update_job_preview()
    
    def _update_job_cpus_per_task(self):
        """Update cpus_per_task based on tasks and cores."""
        if not HAS_PYJOB or self._job_machine is None:
            return
        
        try:
            partition_select = self.query_one("#job-partition", Select)
            partition = partition_select.value if partition_select.value else None
            cores_per_node = self._job_machine.get_effective_value("cores_per_node", partition)
            
            nodes = int(self.query_one("#job-nodes", Input).value or "1")
            ntasks = int(self.query_one("#job-ntasks", Input).value or "1")
            
            tasks_per_node = max(1, ntasks // nodes) if nodes > 0 else ntasks
            cpus_per_task = max(1, cores_per_node // tasks_per_node) if tasks_per_node > 0 else 1
            
            self.query_one("#job-cpus-per-task", Input).value = str(cpus_per_task)
        except (ValueError, ZeroDivisionError):
            pass
    
    @on(Button.Pressed, "#job-btn-generate")
    def on_job_generate(self):
        """Generate/refresh the preview."""
        self._update_job_preview()
    
    @on(Button.Pressed, "#job-btn-save")
    def on_job_save(self):
        """Save the script."""
        if not HAS_PYJOB:
            self.notify("pyjob module not available", severity="error")
            return
        try:
            config = self._get_job_config()
            default_name = f"{config.job_name}.slurm"
            
            # Simple save - write to current directory
            generator = SlurmJobGenerator(config, self._job_machine)
            generator.save(default_name)
            self._update_status(f"Saved to: {default_name}")
            self.notify(f"Script saved to {default_name}", severity="success")
        except Exception as e:
            self.notify(f"Error: {e}", severity="error")
    
    def _action_load_job_machine(self):
        """Load machine config from file."""
        if not HAS_PYJOB:
            return
        
        def handle_file(path):
            if path and path.suffix == ".json":
                try:
                    self._job_machine = MachineConfig.load(str(path))
                    self._update_job_machine_info()
                    self._update_job_preview()
                    self._update_status(f"Loaded machine config: {path.name}")
                except Exception as e:
                    self.notify(f"Error loading config: {e}", severity="error")
        
        self.push_screen(FilePickerDialog(title="Select Machine Config (JSON)"), handle_file)
    
    def _load_job_fleur_input(self, path: Path):
        """Load and analyze FLEUR inp.xml file."""
        if not HAS_PYJOB:
            return
        try:
            analyzer = FleurInputAnalyzer(str(path))
            self._job_fleur_analysis = analyzer.analyze()
            self._update_job_fleur_display()
            self._update_status(f"Analyzed FLEUR input: {path.name}")
            self.notify(f"Loaded {path.name}", severity="information")
        except Exception as e:
            self.notify(f"Error analyzing FLEUR input: {e}", severity="error")
    
    def _update_job_fleur_display(self):
        """Update the FLEUR analysis display."""
        display = self.query_one("#job-fleur-analysis", Static)
        
        if self._job_fleur_analysis is None:
            display.update("No FLEUR input loaded")
            return
        
        a = self._job_fleur_analysis
        element_size = 16 if a.is_complex else 8
        matrix_elements = a.matrix_dimension ** 2
        matrix_memory_gb = matrix_elements * element_size / (1024**3)
        storage_type = "complex" if a.is_complex else "real"
        
        info = f"""[bold yellow]Analysis Results[/bold yellow]
[dim]Kmax:[/dim] [green]{a.kmax:.2f}[/green] a.u.⁻¹  [dim]Volume:[/dim] [green]{a.cell_volume:.1f}[/green] Bohr³
[dim]Spins:[/dim] {a.jspins}  [dim]SOC:[/dim] {'[green]Yes[/green]' if a.has_soc else 'No'}  [dim]Noco:[/dim] {'[green]Yes[/green]' if a.has_noco else 'No'}
[dim]Inversion:[/dim] {'[green]Yes[/green]' if a.has_inversion else '[red]No[/red]'}  [dim]K-points:[/dim] {a.num_kpoints}
[dim]Atoms:[/dim] {a.num_atoms}  [dim]Types:[/dim] {a.num_atom_types} ({', '.join(a.atom_species)})

[bold cyan]N_basis:[/bold cyan] [green]{a.n_basis_per_kpoint:,}[/green]  [bold cyan]Matrix:[/bold cyan] [green]{a.matrix_dimension:,}[/green]{'  (×2 noco)' if a.has_noco else ''}
[dim]Storage:[/dim] {storage_type}  [dim]Memory/k-point:[/dim] [yellow]{matrix_memory_gb:.3f} GB[/yellow]"""
        
        display.update(info)
        self._update_job_parallelization()
    
    def _update_job_parallelization(self):
        """Update the parallelization suggestion display."""
        if not HAS_PYJOB:
            return
        display = self.query_one("#job-parallelization-info", Static)
        
        if self._job_machine is None or self._job_fleur_analysis is None:
            display.update("Load machine config and FLEUR input to see suggestions")
            self._job_parallelization = None
            return
        
        try:
            partition_select = self.query_one("#job-partition", Select)
            partition = partition_select.value if partition_select.value else None
            
            nodes_input = self.query_one("#job-nodes", Input)
            initial_nodes = int(nodes_input.value) if nodes_input.value else 1
            
            effective_gpus = self._job_machine.get_effective_value("gpus_per_node", partition) or 0
            use_gpu = effective_gpus > 0
            
            self._job_parallelization = suggest_parallelization(
                fleur_result=self._job_fleur_analysis,
                machine=self._job_machine,
                partition=partition,
                initial_nodes=initial_nodes,
                use_gpu=use_gpu,
            )
            
            p = self._job_parallelization
            
            info = f"""[bold green]Recommended Settings[/bold green]
[dim]Nodes:[/dim] [green]{p.num_nodes}[/green]  [dim]Tasks:[/dim] [green]{p.num_tasks}[/green]  [dim]Tasks/node:[/dim] [green]{p.tasks_per_node}[/green]
[dim]CPUs/task:[/dim] [green]{p.cpus_per_task}[/green]  [dim]K-pts/task:[/dim] [green]{p.kpoints_per_task}[/green]
[dim]Memory/node:[/dim] [yellow]{p.memory_per_node_gb:.2f} GB[/yellow]"""
            
            if p.uses_gpu and p.tasks_per_gpu:
                info += f"\n[dim]Tasks/GPU:[/dim] [green]{p.tasks_per_gpu}[/green]  [dim]Mem/GPU:[/dim] [yellow]{p.memory_per_gpu_gb:.2f} GB[/yellow]"
            
            if p.notes:
                info += "\n[dim italic]" + "; ".join(p.notes[:2]) + "[/dim italic]"
            
            display.update(info)
            self._apply_job_parallelization_to_form()
            
        except Exception as e:
            display.update(f"[red]Error: {e}[/red]")
            self._job_parallelization = None
    
    def _apply_job_parallelization_to_form(self):
        """Apply parallelization suggestions to form fields."""
        if not HAS_PYJOB or self._job_parallelization is None:
            return
        
        if self._job_updating_from_parallelization:
            return
        
        self._job_updating_from_parallelization = True
        try:
            p = self._job_parallelization
            
            nodes_input = self.query_one("#job-nodes", Input)
            if int(nodes_input.value or "1") != p.num_nodes:
                nodes_input.value = str(p.num_nodes)
            
            self.query_one("#job-ntasks", Input).value = str(p.num_tasks)
            self.query_one("#job-cpus-per-task", Input).value = str(p.cpus_per_task)
            
            if self._job_machine:
                partition_select = self.query_one("#job-partition", Select)
                partition = partition_select.value if partition_select.value else None
                
                gpus_input = self.query_one("#job-gpus", Input)
                gpus_per_node = self._job_machine.get_effective_value("gpus_per_node", partition) or 0
                if gpus_per_node > 0:
                    gpu_type = self._job_machine.get_effective_value("gpu_type", partition) or ""
                    if gpu_type:
                        gpus_input.value = f"{gpu_type}:{gpus_per_node}"
                    else:
                        gpus_input.value = str(gpus_per_node)
                else:
                    gpus_input.value = ""
            
            memory_input = self.query_one("#job-memory", Input)
            suggested_mem_gb = int(p.memory_per_node_gb * 1.2) + 1
            memory_input.value = f"{suggested_mem_gb}G"
        finally:
            self._job_updating_from_parallelization = False


def run_editor(
    schema_path: Optional[str | Path] = None,
    input_file: Optional[str | Path] = None
):
    """Run the FLEUR input editor application."""
    app = FLEURisteApp(schema_path=schema_path, input_file=input_file)
    app.run()
