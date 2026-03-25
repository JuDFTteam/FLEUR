#!/usr/bin/env python3
"""
Textual GUI for SLURM Job Submission File Generator

A terminal-based graphical interface for creating SLURM batch job
submission scripts with configurable options for HPC cluster scheduling.
"""

from pathlib import Path
from typing import Optional, List

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Container, Horizontal, Vertical, VerticalScroll, Grid
from textual.widgets import (
    Header, Footer, Static, Input, Label, Select, Rule, 
    Button, TextArea, Checkbox, TabbedContent, TabPane, DirectoryTree
)
from textual.screen import ModalScreen
from textual import on

# Support both direct execution and package import
try:
    from .slurm_generator import (
        MachineConfig, Partition, SlurmJobConfig, SlurmJobGenerator,
        load_machine_configs, get_machines_directory
    )
    from .fleur_analyzer import FleurInputAnalyzer, FleurAnalysisResult
    from .parallelization import ParallelizationStrategy, ParallelizationResult, suggest_parallelization
except ImportError:
    from slurm_generator import (
        MachineConfig, Partition, SlurmJobConfig, SlurmJobGenerator,
        load_machine_configs, get_machines_directory
    )
    from fleur_analyzer import FleurInputAnalyzer, FleurAnalysisResult
    from parallelization import ParallelizationStrategy, ParallelizationResult, suggest_parallelization


# ==================== Dialogs ====================

class FilePickerDialog(ModalScreen[Optional[Path]]):
    """Dialog for picking a file."""
    
    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]
    
    CSS = """
    FilePickerDialog {
        align: center middle;
    }
    
    FilePickerDialog > Vertical {
        width: 70;
        height: 30;
        border: thick $primary;
        background: $surface;
        padding: 1 2;
    }
    
    FilePickerDialog DirectoryTree {
        height: 20;
        border: solid $primary-darken-2;
    }
    
    FilePickerDialog .dialog-buttons {
        height: 3;
        align: center middle;
    }
    """
    
    def __init__(self, start_path: str = ".", title: str = "Select File", 
                 file_filter: Optional[str] = None):
        super().__init__()
        self.start_path = start_path
        self.title_text = title
        self.file_filter = file_filter
        self.selected_path: Optional[Path] = None
    
    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label(f"[bold]{self.title_text}[/bold]")
            yield Rule()
            yield DirectoryTree(self.start_path, id="file-tree")
            yield Label("", id="selected-label")
            with Horizontal(classes="dialog-buttons"):
                yield Button("Select", variant="primary", id="btn-select")
                yield Button("Cancel", variant="default", id="btn-cancel")
    
    @on(DirectoryTree.FileSelected)
    def on_file_selected(self, event: DirectoryTree.FileSelected):
        self.selected_path = event.path
        label = self.query_one("#selected-label", Label)
        label.update(f"Selected: {event.path.name}")
    
    @on(Button.Pressed, "#btn-select")
    def on_select(self):
        if self.selected_path:
            self.dismiss(self.selected_path)
        
    @on(Button.Pressed, "#btn-cancel")
    def on_cancel(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)


class SaveDialog(ModalScreen[Optional[str]]):
    """Dialog for saving a file."""
    
    BINDINGS = [
        Binding("escape", "cancel", "Cancel"),
    ]
    
    CSS = """
    SaveDialog {
        align: center middle;
    }
    
    SaveDialog > Vertical {
        width: 60;
        height: 12;
        border: thick $primary;
        background: $surface;
        padding: 1 2;
    }
    
    SaveDialog .dialog-buttons {
        height: 3;
        align: center middle;
        margin-top: 1;
    }
    """
    
    def __init__(self, default_name: str = "job.slurm"):
        super().__init__()
        self.default_name = default_name
    
    def compose(self) -> ComposeResult:
        with Vertical():
            yield Label("[bold]Save SLURM Script[/bold]")
            yield Rule()
            yield Label("Filename:")
            yield Input(value=self.default_name, id="filename-input")
            with Horizontal(classes="dialog-buttons"):
                yield Button("Save", variant="primary", id="btn-save")
                yield Button("Cancel", variant="default", id="btn-cancel")
    
    @on(Button.Pressed, "#btn-save")
    def on_save(self):
        filename = self.query_one("#filename-input", Input).value
        if filename:
            self.dismiss(filename)
        
    @on(Button.Pressed, "#btn-cancel")
    def on_cancel(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)


# ==================== Main Application ====================

class SlurmGeneratorApp(App):
    """Textual GUI for SLURM Job Submission File Generator."""
    
    TITLE = "SLURM Job Generator"
    SUB_TITLE = "Generate HPC Job Submission Scripts"
    
    CSS = """
    Screen {
        layout: horizontal;
    }
    
    #left-panel {
        width: 45;
        height: 100%;
        border-right: solid $primary;
        background: $surface-darken-1;
    }
    
    #right-panel {
        width: 1fr;
        height: 100%;
    }
    
    #form-container {
        padding: 1 2;
    }
    
    .form-section {
        margin-bottom: 1;
    }
    
    .form-section-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 0;
    }
    
    .form-row {
        height: 3;
        layout: horizontal;
    }
    
    .form-row Label {
        width: 18;
        height: 1;
        content-align: left middle;
    }
    
    .form-row Input {
        width: 1fr;
        height: 3;
    }
    
    .form-row Select {
        width: 1fr;
        height: 3;
    }
    
    #machine-info {
        height: 8;
        border: solid $primary-darken-2;
        padding: 0 1;
        margin: 1 0;
        overflow-y: auto;
    }
    
    #preview-area {
        height: 1fr;
        border: solid $success;
        margin: 1;
    }
    
    #command-list {
        height: 8;
        border: solid $primary-darken-2;
    }
    
    #module-list {
        height: 5;
        border: solid $primary-darken-2;
    }
    
    .button-bar {
        height: 4;
        align: center middle;
        padding: 0 1;
    }
    
    .button-bar Button {
        margin: 0 1;
    }
    
    #status-bar {
        height: 1;
        dock: bottom;
        background: $primary-darken-2;
        padding: 0 1;
    }
    
    .checkbox-row {
        height: 3;
        layout: horizontal;
    }
    
    .checkbox-row Checkbox {
        width: auto;
        margin-right: 2;
    }
    
    #fleur-analysis {
        height: auto;
        max-height: 16;
        border: solid $warning;
        padding: 0 1;
        margin: 1 0;
        overflow-y: auto;
    }
    
    #fleur-analysis .analysis-title {
        text-style: bold;
        color: $warning;
    }
    
    #fleur-analysis .analysis-value {
        color: $success;
    }
    
    #parallelization-info {
        height: auto;
        max-height: 12;
        border: solid $success;
        padding: 0 1;
        margin: 1 0;
        overflow-y: auto;
    }
    """
    
    BINDINGS = [
        Binding("ctrl+s", "save_script", "Save Script"),
        Binding("ctrl+g", "generate_preview", "Generate Preview"),
        Binding("ctrl+l", "load_machine", "Load Machine"),
        Binding("ctrl+f", "load_fleur_input", "Load FLEUR inp.xml"),
        Binding("ctrl+q", "quit", "Quit"),
        Binding("f5", "generate_preview", "Refresh"),
    ]
    
    def __init__(self):
        super().__init__()
        self.machine: Optional[MachineConfig] = None
        self.available_machines = load_machine_configs()
        self.fleur_analysis: Optional[FleurAnalysisResult] = None
        self.parallelization: Optional[ParallelizationResult] = None
        self._last_machine_modules = ""
    
    def _get_machine_choices(self):
        """Build machine selection dropdown choices."""
        choices = [("(None)", "(None)")]
        for name in sorted(self.available_machines.keys()):
            choices.append((name, name))
        choices.append(("(Load from file...)", "(Load from file...)"))
        return choices
    
    def compose(self) -> ComposeResult:
        yield Header(show_clock=True)
        
        with Container(id="left-panel"):
            with VerticalScroll(id="form-container"):
                # 1. Machine Selection + Partition
                yield Static("[bold]Machine Configuration[/bold]", classes="form-section-title")
                with Horizontal(classes="form-row"):
                    yield Label("Machine:")
                    yield Select(
                        self._get_machine_choices(),
                        value="(None)",
                        id="machine-select"
                    )
                
                with Horizontal(classes="form-row"):
                    yield Label("Partition:")
                    yield Select(
                        [("(select machine first)", "")],
                        value="",
                        id="partition"
                    )
                yield Static("No machine selected", id="machine-info")
                
                yield Rule()
                
                # 2. Recommended Settings + Resources
                yield Static("[bold]Parallelization & Resources[/bold]", classes="form-section-title")
                yield Static("Load machine config and FLEUR input to see suggestions", id="parallelization-info")
                
                with Horizontal(classes="form-row"):
                    yield Label("Job Name:")
                    yield Input(value="my_job", id="job-name")
                
                with Horizontal(classes="form-row"):
                    yield Label("Account:")
                    yield Input(placeholder="(optional)", id="account")
                
                with Horizontal(classes="form-row"):
                    yield Label("Nodes:")
                    yield Input(value="1", id="nodes", type="integer")
                
                with Horizontal(classes="form-row"):
                    yield Label("Tasks:")
                    yield Input(value="1", id="ntasks", type="integer")
                
                with Horizontal(classes="form-row"):
                    yield Label("CPUs/Task:")
                    yield Input(value="1", id="cpus-per-task", type="integer")
                
                with Horizontal(classes="form-row"):
                    yield Label("Memory:")
                    yield Input(value="4G", id="memory", placeholder="e.g., 4G, 8000M")
                
                with Horizontal(classes="form-row"):
                    yield Label("Time Limit:")
                    yield Input(value="01:00:00", id="time-limit", placeholder="HH:MM:SS")
                
                with Horizontal(classes="form-row"):
                    yield Label("GPUs:")
                    yield Input(placeholder="e.g., 1, a100:2", id="gpus")
                
                # 3. Save Button
                with Horizontal(classes="button-bar"):
                    yield Button("Generate", variant="primary", id="btn-generate")
                    yield Button("Save", variant="success", id="btn-save")
                
                yield Rule()
                
                # 4. Analysis Results
                yield Static("[bold]FLEUR Input Analysis[/bold]", classes="form-section-title")
                with Horizontal(classes="form-row"):
                    yield Button("Load inp.xml", variant="warning", id="btn-load-fleur")
                yield Static("No FLEUR input loaded", id="fleur-analysis")
                
                yield Rule()
                
                # 5. Notifications
                yield Static("[bold]Notifications[/bold]", classes="form-section-title")
                with Horizontal(classes="form-row"):
                    yield Label("Email:")
                    yield Input(placeholder="user@example.com", id="mail-user")
                
                with Horizontal(classes="form-row"):
                    yield Label("Notify on:")
                    yield Select(
                        [("END,FAIL", "END,FAIL"), ("BEGIN,END,FAIL", "BEGIN,END,FAIL"), 
                         ("ALL", "ALL"), ("NONE", "NONE")],
                        value="END,FAIL",
                        id="mail-type"
                    )
                
                yield Rule()
                
                # 6. Modules
                yield Static("[bold]Modules[/bold]", classes="form-section-title")
                yield TextArea(id="module-list", language=None)
                yield Label("[dim]One module per line, e.g., python/3.10[/dim]")
                
                yield Rule()
                
                # 7. Commands
                yield Static("[bold]Commands[/bold]", classes="form-section-title")
                yield TextArea(id="command-list", language="bash")
                yield Label("[dim]Enter shell commands to execute[/dim]")
        
        with Container(id="right-panel"):
            yield Static("[bold]SLURM Script Preview[/bold]", id="preview-title")
            yield TextArea(id="preview-area", language="bash", read_only=True)
        
        yield Static("Ready", id="status-bar")
        yield Footer()
    
    def on_mount(self):
        """Initialize after mount."""
        self._update_preview()
    
    def _update_status(self, message: str):
        """Update status bar."""
        self.query_one("#status-bar", Static).update(message)
    
    def _get_job_config(self) -> SlurmJobConfig:
        """Build SlurmJobConfig from form inputs."""
        # Parse modules (one per line)
        modules_text = self.query_one("#module-list", TextArea).text
        modules = [m.strip() for m in modules_text.split("\n") if m.strip()]
        
        # Parse commands (one per line)
        commands_text = self.query_one("#command-list", TextArea).text
        commands = [c for c in commands_text.split("\n") if c.strip()]
        
        # Get form values
        job_name = self.query_one("#job-name", Input).value or "my_job"
        account = self.query_one("#account", Input).value or None
        partition_select = self.query_one("#partition", Select)
        partition = partition_select.value if partition_select.value else None
        
        nodes = int(self.query_one("#nodes", Input).value or "1")
        ntasks = int(self.query_one("#ntasks", Input).value or "1")
        cpus = int(self.query_one("#cpus-per-task", Input).value or "1")
        memory = self.query_one("#memory", Input).value or "4G"
        time_limit = self.query_one("#time-limit", Input).value or "01:00:00"
        gpus = self.query_one("#gpus", Input).value or None
        
        mail_user = self.query_one("#mail-user", Input).value or None
        mail_type_select = self.query_one("#mail-type", Select)
        mail_type = mail_type_select.value if mail_user else None
        
        # Get gpu_syntax and use_mem_option from machine config (partition or machine default)
        gpu_syntax = "gpus"  # Default to modern syntax
        use_mem_option = True  # Default to using --mem
        if self.machine:
            gpu_syntax = self.machine.get_effective_value("gpu_syntax", partition) or "gpus"
            use_mem = self.machine.get_effective_value("use_mem_option", partition)
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
    
    def _update_preview(self):
        """Generate and update the preview."""
        try:
            config = self._get_job_config()
            generator = SlurmJobGenerator(config, self.machine)
            script = generator.generate()
            
            preview = self.query_one("#preview-area", TextArea)
            preview.load_text(script)
            
            self._update_status(f"Preview updated | Job: {config.job_name}")
        except Exception as e:
            self._update_status(f"Error: {e}")
    
    def _update_machine_info(self):
        """Update machine info display."""
        info_panel = self.query_one("#machine-info", Static)
        
        if self.machine:
            # Update partition dropdown choices - first partition is default
            partition_select = self.query_one("#partition", Select)
            partition_choices = []
            for p in self.machine.partitions:
                partition_choices.append((p.name, p.name))
            if partition_choices:
                partition_select.set_options(partition_choices)
                partition_select.value = partition_choices[0][1]  # Select first partition
            else:
                partition_select.set_options([("(no partitions)", "")])
            
            # Suggest modules from selected/default partition
            partition_value = partition_select.value if partition_select.value else None
            self._update_module_suggestion(partition_value)
            
            # Suggest command from machine config (will be updated per partition)
            self._update_command_suggestion()
            
            # Update the partition info display
            self._update_partition_info()
        else:
            machines_dir = get_machines_directory()
            if not self.available_machines:
                info_panel.update(f"No machine configs found in:\\n{machines_dir}\\nAdd .json files or use 'Load from file...'")
            else:
                info_panel.update("No machine selected")
            # Reset partition dropdown
            partition_select = self.query_one("#partition", Select)
            partition_select.set_options([("(select machine first)", "")])
    
    def _update_partition_info(self):
        """Update the partition info display with effective values."""
        info_panel = self.query_one("#machine-info", Static)
        
        if self.machine is None:
            info_panel.update("No machine selected")
            return
        
        # Get current partition
        partition_select = self.query_one("#partition", Select)
        partition = partition_select.value if partition_select.value else None
        
        # Get effective values for this partition
        cores = self.machine.get_effective_value("cores_per_node", partition)
        gpus = self.machine.get_effective_value("gpus_per_node", partition) or 0
        memory = self.machine.get_effective_value("memory_per_node_gb", partition)
        max_runtime = self.machine.get_effective_value("max_runtime", partition)
        
        partition_name = partition if partition else "(none)"
        info = f"""[bold]{self.machine.name}[/bold] - [cyan]{partition_name}[/cyan]
[dim]Cores:[/dim] [green]{cores}[/green]/node  [dim]GPUs:[/dim] [green]{gpus}[/green]/node  [dim]Memory:[/dim] [green]{memory}[/green]GB"""
        
        if max_runtime:
            info += f"\n[dim]Max runtime:[/dim] {max_runtime}"
        
        info_panel.update(info)
        
        self._update_module_suggestion(partition)

        # Update command suggestion for this partition
        self._update_command_suggestion()

    def _update_module_suggestion(self, partition):
        """Update module textarea with machine/partition modules without overwriting user edits."""
        if self.machine is None:
            return

        suggested_modules = self.machine.get_effective_value("modules_needed", partition) or []
        suggested_text = "\n".join(suggested_modules[:3]) if suggested_modules else ""

        modules_area = self.query_one("#module-list", TextArea)
        current_text = modules_area.text.strip()
        if (not current_text) or (current_text == self._last_machine_modules):
            modules_area.load_text(suggested_text)
            self._last_machine_modules = suggested_text
    
    def _update_command_suggestion(self):
        """Update the command textarea with machine/partition command if user hasn't modified it."""
        if self.machine is None:
            return
        
        # Get current partition
        partition_select = self.query_one("#partition", Select)
        partition = partition_select.value if partition_select.value else None
        
        # Get effective command for this partition
        command = self.machine.get_effective_value("command", partition)
        
        if command:
            command_area = self.query_one("#command-list", TextArea)
            # Only update if the area is empty or contains the previous machine command
            current_text = command_area.text.strip()
            if not current_text or current_text == getattr(self, '_last_machine_command', ''):
                command_area.load_text(command)
                self._last_machine_command = command
    
    @on(Select.Changed, "#machine-select")
    def on_machine_changed(self, event: Select.Changed):
        """Handle machine selection change."""
        value = event.value
        
        if value == "(None)":
            self.machine = None
        elif value == "(Load from file...)":
            self.action_load_machine()
            return
        elif value in self.available_machines:
            self.machine = self.available_machines[value]
            self._update_status(f"Loaded machine: {value}")
        
        self._update_machine_info()
        self._update_parallelization()
        self._update_preview()
    
    @on(Select.Changed, "#partition")
    def on_partition_changed(self, event: Select.Changed):
        """Handle partition selection change - recalculate parallelization."""
        self._update_partition_info()
        self._update_parallelization()
        self._update_preview()
    
    @on(Input.Changed, "#nodes")
    def on_nodes_changed(self, event: Input.Changed):
        """Handle nodes input change - recalculate parallelization."""
        # Skip if we're auto-updating from parallelization
        if hasattr(self, '_updating_from_parallelization') and self._updating_from_parallelization:
            return
        self._update_parallelization()
        self._update_preview()
    
    @on(Input.Changed, "#ntasks")
    def on_ntasks_changed(self, event: Input.Changed):
        """Handle tasks input change - recalculate cpus_per_task."""
        # Skip if we're auto-updating from parallelization
        if hasattr(self, '_updating_from_parallelization') and self._updating_from_parallelization:
            return
        self._update_cpus_per_task()
        self._update_preview()
    
    def _update_cpus_per_task(self):
        """Update cpus_per_task so that cpus_per_task * tasks_per_node = cores_per_node."""
        if self.machine is None:
            return
        
        try:
            partition_select = self.query_one("#partition", Select)
            partition = partition_select.value if partition_select.value else None
            
            cores_per_node = self.machine.get_effective_value("cores_per_node", partition)
            
            nodes_input = self.query_one("#nodes", Input)
            nodes = int(nodes_input.value) if nodes_input.value else 1
            
            ntasks_input = self.query_one("#ntasks", Input)
            ntasks = int(ntasks_input.value) if ntasks_input.value else 1
            
            # Calculate tasks per node
            tasks_per_node = max(1, ntasks // nodes) if nodes > 0 else ntasks
            
            # Calculate cpus_per_task = cores_per_node / tasks_per_node
            cpus_per_task = max(1, cores_per_node // tasks_per_node) if tasks_per_node > 0 else 1
            
            cpus_input = self.query_one("#cpus-per-task", Input)
            cpus_input.value = str(cpus_per_task)
        except (ValueError, ZeroDivisionError):
            pass

    @on(Button.Pressed, "#btn-generate")
    def on_generate(self):
        """Generate/refresh the preview."""
        self._update_preview()
    
    @on(Button.Pressed, "#btn-save")
    def on_save(self):
        """Save the script."""
        self.action_save_script()
    
    @on(Input.Changed)
    @on(TextArea.Changed)
    @on(Checkbox.Changed)
    def on_form_changed(self, event):
        """Auto-update preview on form changes."""
        # Debounce would be nice here, but for simplicity update immediately
        self._update_preview()
    
    def action_generate_preview(self):
        """Generate preview action."""
        self._update_preview()
    
    def action_save_script(self):
        """Save script action."""
        try:
            config = self._get_job_config()
            default_name = f"{config.job_name}.slurm"
            
            def handle_save(filename: Optional[str]):
                if filename:
                    generator = SlurmJobGenerator(config, self.machine)
                    generator.save(filename)
                    self._update_status(f"Saved to: {filename}")
                    self.notify(f"Script saved to {filename}", severity="success")
            
            self.push_screen(SaveDialog(default_name), handle_save)
        except Exception as e:
            self.notify(f"Error: {e}", severity="error")
    
    def action_load_machine(self):
        """Load machine config from file."""
        def handle_file(path: Optional[Path]):
            if path and path.suffix == ".json":
                try:
                    self.machine = MachineConfig.load(str(path))
                    self._update_machine_info()
                    self._update_preview()
                    self._update_status(f"Loaded machine config: {path.name}")
                    
                    # Update select to show loaded machine
                    # (ideally would add it to the list, but that's complex)
                except Exception as e:
                    self.notify(f"Error loading config: {e}", severity="error")
        
        self.push_screen(FilePickerDialog(title="Select Machine Config (JSON)"), handle_file)
    
    @on(Button.Pressed, "#btn-load-fleur")
    def on_load_fleur_button(self):
        """Handle FLEUR input load button."""
        self.action_load_fleur_input()
    
    def action_load_fleur_input(self):
        """Load and analyze FLEUR inp.xml file."""
        def handle_file(path: Optional[Path]):
            if path and path.suffix == ".xml":
                try:
                    analyzer = FleurInputAnalyzer(str(path))
                    self.fleur_analysis = analyzer.analyze()
                    self._update_fleur_display()
                    self._update_status(f"Analyzed FLEUR input: {path.name}")
                    self.notify(f"Loaded {path.name}", severity="information")
                except Exception as e:
                    self.notify(f"Error analyzing FLEUR input: {e}", severity="error")
        
        self.push_screen(FilePickerDialog(title="Select FLEUR inp.xml"), handle_file)
    
    def _update_fleur_display(self):
        """Update the FLEUR analysis display."""
        display = self.query_one("#fleur-analysis", Static)
        
        if self.fleur_analysis is None:
            display.update("No FLEUR input loaded")
            return
        
        a = self.fleur_analysis
        
        # Calculate memory estimate (complex if no inversion or noco)
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
        
        # Also update parallelization when FLEUR analysis changes
        self._update_parallelization()
        
        # Also update parallelization when FLEUR analysis changes
        self._update_parallelization()
    
    def _update_parallelization(self):
        """Update the parallelization suggestion display."""
        display = self.query_one("#parallelization-info", Static)
        
        # Need both machine and FLEUR analysis
        if self.machine is None or self.fleur_analysis is None:
            display.update("Load machine config and FLEUR input to see suggestions")
            self.parallelization = None
            return
        
        try:
            # Get current values from UI
            partition_select = self.query_one("#partition", Select)
            partition = partition_select.value if partition_select.value else None
            
            nodes_input = self.query_one("#nodes", Input)
            initial_nodes = int(nodes_input.value) if nodes_input.value else 1
            
            # Auto-use GPUs if machine has them
            effective_gpus = self.machine.get_effective_value("gpus_per_node", partition) or 0
            use_gpu = effective_gpus > 0
            
            # Calculate parallelization
            self.parallelization = suggest_parallelization(
                fleur_result=self.fleur_analysis,
                machine=self.machine,
                partition=partition,
                initial_nodes=initial_nodes,
                use_gpu=use_gpu,
            )
            
            p = self.parallelization
            
            # Build display
            info = f"""[bold green]Recommended Settings[/bold green]
[dim]Nodes:[/dim] [green]{p.num_nodes}[/green]  [dim]Tasks:[/dim] [green]{p.num_tasks}[/green]  [dim]Tasks/node:[/dim] [green]{p.tasks_per_node}[/green]
[dim]CPUs/task:[/dim] [green]{p.cpus_per_task}[/green]  [dim]K-pts/task:[/dim] [green]{p.kpoints_per_task}[/green]
[dim]Memory/node:[/dim] [yellow]{p.memory_per_node_gb:.2f} GB[/yellow]"""
            
            if p.uses_gpu and p.tasks_per_gpu:
                info += f"\n[dim]Tasks/GPU:[/dim] [green]{p.tasks_per_gpu}[/green]  [dim]Mem/GPU:[/dim] [yellow]{p.memory_per_gpu_gb:.2f} GB[/yellow]"
            
            if p.notes:
                info += "\n[dim italic]" + "; ".join(p.notes[:2]) + "[/dim italic]"
            
            display.update(info)
            
            # Auto-update the form fields with suggested values
            self._apply_parallelization_to_form()
            
        except Exception as e:
            display.update(f"[red]Error: {e}[/red]")
            self.parallelization = None
    
    def _apply_parallelization_to_form(self):
        """Apply parallelization suggestions to form fields (without triggering updates)."""
        if self.parallelization is None:
            return
        
        p = self.parallelization
        
        # Update form fields - use _updating flag to prevent infinite loop
        if not hasattr(self, '_updating_from_parallelization'):
            self._updating_from_parallelization = False
        
        if self._updating_from_parallelization:
            return
        
        self._updating_from_parallelization = True
        try:
            # Only update if values differ significantly
            nodes_input = self.query_one("#nodes", Input)
            if int(nodes_input.value or "1") != p.num_nodes:
                nodes_input.value = str(p.num_nodes)
            
            ntasks_input = self.query_one("#ntasks", Input)
            ntasks_input.value = str(p.num_tasks)
            
            cpus_input = self.query_one("#cpus-per-task", Input)
            cpus_input.value = str(p.cpus_per_task)
            
            # Always set GPUs based on partition's gpus_per_node (if > 0)
            if self.machine:
                # Get partition to determine effective GPU settings
                partition_select = self.query_one("#partition", Select)
                partition = partition_select.value if partition_select.value else None
                
                gpus_input = self.query_one("#gpus", Input)
                gpus_per_node = self.machine.get_effective_value("gpus_per_node", partition) or 0
                if gpus_per_node > 0:
                    gpu_type = self.machine.get_effective_value("gpu_type", partition) or ""
                    if gpu_type:
                        gpus_input.value = f"{gpu_type}:{gpus_per_node}"
                    else:
                        gpus_input.value = str(gpus_per_node)
                else:
                    # Clear GPU field if partition has no GPUs
                    gpus_input.value = ""
            
            # Update memory suggestion
            memory_input = self.query_one("#memory", Input)
            suggested_mem_gb = int(p.memory_per_node_gb * 1.2) + 1  # 20% margin
            memory_input.value = f"{suggested_mem_gb}G"
        finally:
            self._updating_from_parallelization = False


def main():
    """Main entry point."""
    app = SlurmGeneratorApp()
    app.run()


if __name__ == "__main__":
    main()
