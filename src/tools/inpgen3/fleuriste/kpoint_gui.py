"""
K-Point Manager GUI components for PyInpEdit.
Provides widgets and screens for managing k-points in FLEUR inp.xml files.
"""

from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical, ScrollableContainer
from textual.widgets import (
    Button, Label, Input, Static, 
    DataTable, Select, Tree, Rule
)
from textual.binding import Binding
from textual.screen import ModalScreen
from textual import on
from pathlib import Path
from typing import Optional, Union
from fractions import Fraction

try:
    from textual_plotext import PlotextPlot
    HAS_PLOTEXT = True
except ImportError:
    HAS_PLOTEXT = False
    PlotextPlot = None

from .kpoint_manager import InpXMLManager, KPoint, KPointList, KPointMode, KPointModifiers


def evaluate_coordinate(coord: Union[float, str]) -> float:
    """Evaluate a k-point coordinate that may be a float or string expression.
    
    Args:
        coord: Coordinate value as float or string expression (e.g., "1/3", "0.5")
    
    Returns:
        Float value of the coordinate
    """
    if isinstance(coord, (int, float)):
        return float(coord)
    
    coord_str = str(coord).strip()
    
    try:
        return float(coord_str)
    except ValueError:
        pass
    
    try:
        return float(Fraction(coord_str))
    except (ValueError, ZeroDivisionError):
        pass
    
    try:
        allowed_chars = set('0123456789.+-*/()')
        if all(c in allowed_chars or c.isspace() for c in coord_str):
            return float(eval(coord_str, {"__builtins__": {}}, {}))
    except (ValueError, SyntaxError, ZeroDivisionError):
        pass
    
    return 0.0


if HAS_PLOTEXT:
    class KPointPlot(PlotextPlot):
        """Widget for visualizing k-points using textual-plotext."""
        
        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            self.kpoint_list: Optional[KPointList] = None
            self.projection: str = "xy"
        
        def on_mount(self) -> None:
            """Initialize the plot when mounted."""
            self._update_plot()
        
        def set_kpoint_list(self, kplist: Optional[KPointList], projection: str = "xy"):
            """Set the k-point list to visualize."""
            self.kpoint_list = kplist
            self.projection = projection
            self._update_plot()
        
        def _update_plot(self):
            """Update the plot with current k-point data."""
            plt = self.plt
            plt.clear_figure()
            
            if self.kpoint_list is None or len(self.kpoint_list.kpoints) == 0:
                plt.title("No k-points to display")
                self.refresh()
                return
            
            try:
                kpoints = self.kpoint_list.kpoints
                
                if self.projection == "xy":
                    x_coords = [evaluate_coordinate(kp.x) for kp in kpoints]
                    y_coords = [evaluate_coordinate(kp.y) for kp in kpoints]
                    x_label, y_label = "kx", "ky"
                elif self.projection == "xz":
                    x_coords = [evaluate_coordinate(kp.x) for kp in kpoints]
                    y_coords = [evaluate_coordinate(kp.z) for kp in kpoints]
                    x_label, y_label = "kx", "kz"
                else:
                    x_coords = [evaluate_coordinate(kp.y) for kp in kpoints]
                    y_coords = [evaluate_coordinate(kp.z) for kp in kpoints]
                    x_label, y_label = "ky", "kz"
                
                plt.scatter(x_coords, y_coords, marker="x")
                plt.xlabel(x_label)
                plt.ylabel(y_label)
                plt.title(f"{self.kpoint_list.name} ({len(kpoints)} k-points, {self.projection.upper()} plane)")
                self.refresh()
                
            except Exception as e:
                plt.clear_figure()
                plt.title(f"Error creating plot: {e}")
                self.refresh()
else:
    # Fallback when plotext is not available
    class KPointPlot(Static):
        """Fallback widget when textual-plotext is not available."""
        
        def __init__(self, id: str = None, **kwargs):
            super().__init__("Install textual-plotext for k-point visualization", id=id, **kwargs)
            self.kpoint_list = None
            self.projection = "xy"
        
        def set_kpoint_list(self, kplist, projection="xy"):
            self.kpoint_list = kplist
            self.projection = projection
            if kplist:
                self.update(f"K-points: {len(kplist.kpoints)} ({projection.upper()} view)")
            else:
                self.update("No k-points")


class CreateDensityDialog(ModalScreen):
    """Modal dialog for creating k-points by density or number."""
    
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    
    CSS = """
    #kpoint-dialog {
        width: 60;
        height: auto;
        background: $panel;
        border: thick $primary;
        padding: 2;
    }
    
    #kpoint-dialog-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 1;
    }
    
    #kpoint-dialog-buttons {
        margin-top: 1;
        height: auto;
    }
    """
    
    def compose(self) -> ComposeResult:
        with Container(id="kpoint-dialog"):
            yield Label("Create K-Points (Density/Number)", id="kpoint-dialog-title")
            yield Label("List name:")
            yield Input(placeholder="kpoints", id="input-name")
            yield Label("Mode:")
            yield Select(
                [("Density", "density"), ("Number", "number")],
                value="density",
                id="select-mode"
            )
            yield Label("Value:")
            yield Input(placeholder="0.1 (density) or 100 (number)", id="input-value")
            yield Label("Options:")
            with Horizontal():
                yield Button("Gamma", id="btn-gamma", variant="default")
                yield Button("Symmetry ✓", id="btn-symmetry", variant="success")
            with Horizontal():
                yield Button("Gauss", id="btn-gauss", variant="default")
                yield Button("Tria", id="btn-tria", variant="default")
            with Horizontal():
                yield Button("Tetra", id="btn-tetra", variant="default")
                yield Button("SOC", id="btn-soc", variant="default")
            
            with Horizontal(id="kpoint-dialog-buttons"):
                yield Button("Create", variant="primary", id="btn-create")
                yield Button("Cancel", variant="default", id="btn-cancel")
    
    def __init__(self):
        super().__init__()
        self.use_gamma = False
        self.use_symmetry = True
        self.use_gauss = False
        self.use_tria = False
        self.use_tetra = False
        self.use_soc = False
    
    def _toggle_button(self, event, attr_name, label_base):
        """Helper to toggle a button state."""
        current = getattr(self, attr_name)
        setattr(self, attr_name, not current)
        btn = event.button
        if getattr(self, attr_name):
            btn.label = f"{label_base} ✓"
            btn.variant = "success"
        else:
            btn.label = label_base
            btn.variant = "default"
    
    @on(Button.Pressed, "#btn-gamma")
    def toggle_gamma(self, event):
        self._toggle_button(event, "use_gamma", "Gamma")
    
    @on(Button.Pressed, "#btn-symmetry")
    def toggle_symmetry(self, event):
        self._toggle_button(event, "use_symmetry", "Symmetry")
    
    @on(Button.Pressed, "#btn-gauss")
    def toggle_gauss(self, event):
        self._toggle_button(event, "use_gauss", "Gauss")
    
    @on(Button.Pressed, "#btn-tria")
    def toggle_tria(self, event):
        self._toggle_button(event, "use_tria", "Tria")
    
    @on(Button.Pressed, "#btn-tetra")
    def toggle_tetra(self, event):
        self._toggle_button(event, "use_tetra", "Tetra")
    
    @on(Button.Pressed, "#btn-soc")
    def toggle_soc(self, event):
        self._toggle_button(event, "use_soc", "SOC")
    
    @on(Button.Pressed, "#btn-create")
    def on_create_button(self):
        try:
            name = self.query_one("#input-name", Input).value or "kpoints"
            mode = self.query_one("#select-mode", Select).value
            value_str = self.query_one("#input-value", Input).value
            
            if not value_str:
                raise ValueError("Value is required")
            
            if mode == "density":
                value = float(value_str)
                if value <= 0:
                    raise ValueError("Density must be positive")
            else:
                value = int(value_str)
                if value <= 0:
                    raise ValueError("Number of k-points must be positive")
            
            modifiers = {
                'gamma': self.use_gamma,
                'gauss': self.use_gauss,
                'tria': self.use_tria,
                'tetra': self.use_tetra,
                'soc': self.use_soc
            }
            
            self.dismiss((name, mode, value, modifiers, self.use_symmetry))
        
        except ValueError as e:
            self.notify(f"Invalid input: {e}", severity="error")
    
    @on(Button.Pressed, "#btn-cancel")
    def on_cancel_button(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)


class CreatePathDialog(ModalScreen):
    """Modal dialog for creating a k-point path."""
    
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    
    CSS = """
    #kpoint-dialog {
        width: 60;
        height: auto;
        background: $panel;
        border: thick $primary;
        padding: 2;
    }
    
    #kpoint-dialog-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 1;
    }
    
    #kpoint-dialog-buttons {
        margin-top: 1;
        height: auto;
    }
    """
    
    def compose(self) -> ComposeResult:
        with Container(id="kpoint-dialog"):
            yield Label("Create K-Point Path", id="kpoint-dialog-title")
            yield Label("Path name:")
            yield Input(placeholder="bands", id="input-name")
            yield Label("Lattice type:")
            yield Select(
                [(label, label) for label in ["fcc", "bcc", "sc"]],
                value="fcc",
                id="select-lattice"
            )
            yield Label("K-point path (e.g., G-X-W-L-G):")
            yield Input(placeholder="G-X-W-L-G", id="input-path")
            yield Label("Number of k-points:")
            yield Input(placeholder="100", id="input-npoints")
            
            with Horizontal(id="kpoint-dialog-buttons"):
                yield Button("Create", variant="primary", id="btn-create")
                yield Button("Cancel", variant="default", id="btn-cancel")
    
    @on(Button.Pressed, "#btn-create")
    def on_create_button(self):
        try:
            name = self.query_one("#input-name", Input).value or "bands"
            lattice = self.query_one("#select-lattice", Select).value
            path = self.query_one("#input-path", Input).value or "G-X-W-L-G"
            npoints = int(self.query_one("#input-npoints", Input).value or "100")
            
            if npoints <= 0:
                raise ValueError("Number of points must be positive")
            
            self.dismiss((name, lattice, path, npoints))
        
        except ValueError as e:
            self.notify(f"Invalid input: {e}", severity="error")
    
    @on(Button.Pressed, "#btn-cancel")
    def on_cancel_button(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)


class CreateMeshDialog(ModalScreen):
    """Modal dialog for creating a k-point mesh."""
    
    BINDINGS = [Binding("escape", "cancel", "Cancel")]
    
    CSS = """
    #kpoint-dialog {
        width: 60;
        height: auto;
        background: $panel;
        border: thick $primary;
        padding: 2;
    }
    
    #kpoint-dialog-title {
        text-style: bold;
        color: $primary;
        margin-bottom: 1;
    }
    
    #kpoint-dialog-buttons {
        margin-top: 1;
        height: auto;
    }
    """
    
    def compose(self) -> ComposeResult:
        with Container(id="kpoint-dialog"):
            yield Label("Create K-Point Mesh", id="kpoint-dialog-title")
            yield Label("Mesh name:")
            yield Input(placeholder="mesh", id="input-name")
            yield Label("Number of divisions in X:")
            yield Input(placeholder="8", id="input-nx")
            yield Label("Number of divisions in Y:")
            yield Input(placeholder="8", id="input-ny")
            yield Label("Number of divisions in Z:")
            yield Input(placeholder="8", id="input-nz")
            yield Label("Options:")
            with Horizontal():
                yield Button("Gamma ✓", id="btn-gamma", variant="success")
                yield Button("Symmetry ✓", id="btn-symmetry", variant="success")
            with Horizontal():
                yield Button("Gauss", id="btn-gauss", variant="default")
                yield Button("Tria", id="btn-tria", variant="default")
            with Horizontal():
                yield Button("Tetra", id="btn-tetra", variant="default")
                yield Button("SOC", id="btn-soc", variant="default")
            
            with Horizontal(id="kpoint-dialog-buttons"):
                yield Button("Create", variant="primary", id="btn-create")
                yield Button("Cancel", variant="default", id="btn-cancel")
    
    def __init__(self):
        super().__init__()
        self.gamma_centered = True
        self.use_symmetry = True
        self.use_gauss = False
        self.use_tria = False
        self.use_tetra = False
        self.use_soc = False
    
    def _toggle_button(self, event, attr_name, label_base):
        """Helper to toggle a button state."""
        current = getattr(self, attr_name)
        setattr(self, attr_name, not current)
        btn = event.button
        if getattr(self, attr_name):
            btn.label = f"{label_base} ✓"
            btn.variant = "success"
        else:
            btn.label = label_base
            btn.variant = "default"
    
    @on(Button.Pressed, "#btn-gamma")
    def toggle_gamma(self, event):
        self._toggle_button(event, "gamma_centered", "Gamma")
    
    @on(Button.Pressed, "#btn-symmetry")
    def toggle_symmetry(self, event):
        self._toggle_button(event, "use_symmetry", "Symmetry")
    
    @on(Button.Pressed, "#btn-gauss")
    def toggle_gauss(self, event):
        self._toggle_button(event, "use_gauss", "Gauss")
    
    @on(Button.Pressed, "#btn-tria")
    def toggle_tria(self, event):
        self._toggle_button(event, "use_tria", "Tria")
    
    @on(Button.Pressed, "#btn-tetra")
    def toggle_tetra(self, event):
        self._toggle_button(event, "use_tetra", "Tetra")
    
    @on(Button.Pressed, "#btn-soc")
    def toggle_soc(self, event):
        self._toggle_button(event, "use_soc", "SOC")
    
    @on(Button.Pressed, "#btn-create")
    def on_create_button(self):
        try:
            name = self.query_one("#input-name", Input).value or "mesh"
            nx = int(self.query_one("#input-nx", Input).value or "8")
            ny = int(self.query_one("#input-ny", Input).value or "8")
            nz = int(self.query_one("#input-nz", Input).value or "8")
            
            if nx <= 0 or ny <= 0 or nz <= 0:
                raise ValueError("Divisions must be positive")
            
            modifiers = {
                'gamma': self.gamma_centered,
                'gauss': self.use_gauss,
                'tria': self.use_tria,
                'tetra': self.use_tetra,
                'soc': self.use_soc
            }
            self.dismiss((name, nx, ny, nz, modifiers, self.use_symmetry))
        
        except ValueError as e:
            self.notify(f"Invalid input: {e}", severity="error")
    
    @on(Button.Pressed, "#btn-cancel")
    def on_cancel_button(self):
        self.dismiss(None)
    
    def action_cancel(self):
        self.dismiss(None)


# Standard k-point paths for different lattice types
STANDARD_KPOINTS = {
    'fcc': {
        'G': (0.0, 0.0, 0.0), 'Γ': (0.0, 0.0, 0.0),
        'X': (0.5, 0.0, 0.5), 'W': (0.5, 0.25, 0.75),
        'K': (0.375, 0.375, 0.75), 'L': (0.5, 0.5, 0.5),
        'U': (0.625, 0.25, 0.625)
    },
    'bcc': {
        'G': (0.0, 0.0, 0.0), 'Γ': (0.0, 0.0, 0.0),
        'H': (0.5, -0.5, 0.5), 'N': (0.0, 0.0, 0.5),
        'P': (0.25, 0.25, 0.25)
    },
    'sc': {
        'G': (0.0, 0.0, 0.0), 'Γ': (0.0, 0.0, 0.0),
        'X': (0.5, 0.0, 0.0), 'M': (0.5, 0.5, 0.0),
        'R': (0.5, 0.5, 0.5)
    }
}
