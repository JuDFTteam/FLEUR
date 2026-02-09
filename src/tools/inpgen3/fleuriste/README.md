# FLEURiste - FLEUR Input File Editor

A terminal-based graphical editor for FLEUR `inp.xml` files, using the [Textual](https://textual.textualize.io/) framework.

## Features

### Inpgen Mode (Default when no inp.xml)
- **Generate inp.xml**: Create FLEUR input files from scratch
- **Load input files**: Load namelist input files for editing
- **Load structures**: Import structure files (POSCAR, CIF, XYZ, etc.) using ASE
- **Profile selection**: Choose inpgen profiles (default, careful, fast)
- **Options**: Toggle symmetry detection and film mode
- **Direct editing**: Edit namelist input in a text area before generation

### XML Editor Mode
- **Schema-aware editing**: Uses `FleurInputSchema.xsd` to provide intelligent editing features
- **Tree-based navigation**: Navigate the XML structure with an interactive tree view
- **Attribute editing**: Edit element attributes with type-aware input fields
  - Boolean attributes show as T/F dropdowns
  - Enumerated types show as selection dropdowns
  - Numeric and text fields with type hints
- **Validation**: Real-time validation against the XSD schema
- **Undo/Redo**: Full undo/redo support for all edits
- **Element management**: Add, delete, and duplicate elements respecting schema constraints
- **Search**: Search schema documentation and navigate to elements

### K-Point Manager Mode
- **View k-point lists**: Browse all k-point lists defined in inp.xml
- **Create k-points**: Create meshes, paths, or density-based k-point sets
- **Visualize k-points**: 2D projection plots (XY, XZ, YZ planes)
- **Manage lists**: Select active list, delete lists, save changes
- **Integration with inpgen**: Uses FLEUR inpgen library for k-point generation

### Job Generator Mode
- **SLURM script generation**: Create HPC job submission scripts
- **Machine configurations**: Load machine configs from `~/.fleuriste/machines/`
- **FLEUR analysis**: Analyze inp.xml to suggest parallelization settings
- **Auto-parallelization**: Automatic task/node/CPU suggestions based on problem size
- **GPU support**: Configure GPU jobs with proper SLURM directives

## Installation

### From the fleuriste directory:

```bash
cd fleuriste
pip install -e .

# For k-point visualization (optional):
pip install -e ".[plotting]"

# For structure file support (optional):
pip install ase
```

### Or install dependencies directly:

```bash
pip install textual>=0.47.0
```

## Usage

### Command Line

```bash
# Start with empty editor (Inpgen mode)
fleuriste

# Edit an existing inp.xml file
fleuriste inp.xml

# Use a custom schema file
fleuriste -s /path/to/FleurInputSchema.xsd inp.xml

# Start with schema but no input file
fleuriste --schema-only
```

### As a Python module

```bash
python -m fleuriste inp.xml
```

### From Python code

```python
from fleuriste import FLEURisteApp

app = FLEURisteApp(
    schema_path="/path/to/FleurInputSchema.xsd",
    input_file="inp.xml"
)
app.run()
```

## Keyboard Shortcuts

| Key | Action |
|-----|--------|
| `Ctrl+S` | Save file |
| `Ctrl+Z` | Undo |
| `Ctrl+Y` | Redo |
| `Ctrl+F` | Focus search |
| `Ctrl+Q` | Quit |
| `F1` | Switch to XML Editor mode |
| `F2` | Switch to K-Point Manager mode |
| `F3` | Switch to Inpgen mode |
| `F4` | Switch to Job Generator mode |
| `F5` | Refresh |
| `Delete` | Delete selected element |
| `Escape` | Clear search |

## Schema Auto-Detection

FLEURiste automatically searches for `FleurInputSchema.xsd` in these locations:
1. `../src/fleur/io/xml/FleurInputSchema.xsd` (relative to fleuriste)
2. `../build/FleurInputSchema.xsd` (relative to fleuriste)
3. Current working directory
4. `./build/` subdirectory

## Machine Configurations

Machine configurations for the Job Generator are stored in `~/.fleuriste/machines/` as JSON files.

Example machine config (`~/.fleuriste/machines/my_cluster.json`):
```json
{
  "name": "my_cluster",
  "cores_per_node": 48,
  "memory_per_node_gb": 192,
  "modules_needed": ["intel", "mpi"],
  "partitions": [
    {
      "name": "batch",
      "max_runtime": "24:00:00",
      "command": "srun ./fleur"
    },
    {
      "name": "gpu",
      "gpus_per_node": 4,
      "gpu_type": "a100",
      "command": "srun ./fleur_gpu"
    }
  ]
}
```

## Project Structure

```
fleuriste/
├── __init__.py       # Package exports
├── __main__.py       # CLI entry point
├── app.py            # Main Textual application
├── schema_parser.py  # XSD schema parsing
├── xml_editor.py     # XML document handling
├── kpoint_gui.py     # K-point visualization widgets
├── kpoint_manager.py # K-point management logic
├── inpgen_gui.py     # Inpgen mode widgets
└── pyjob/            # Job generator components
    ├── slurm_generator.py
    ├── fleur_analyzer.py
    └── parallelization.py
```

## Requirements

- Python >= 3.10
- textual >= 0.47.0

### Optional Dependencies

- `textual-plotext` - For k-point visualization plots
- `ase` - For loading structure files (POSCAR, CIF, XYZ, etc.)
- `FleurInpgen` - For inp.xml generation (from FLEUR build)
