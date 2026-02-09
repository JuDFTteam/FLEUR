"""
FLEURiste CLI - Command-line interface using Click.

Provides commands for launching the TUI editor, generating inp.xml files,
managing k-points, and generating SLURM job scripts.
"""

import sys
from pathlib import Path
from typing import Optional

import click


def find_schema() -> Optional[Path]:
    """Try to find the FleurInputSchema.xsd file."""
    possible_paths = [
        Path(__file__).parent.parent / "src" / "fleur" / "io" / "xml" / "FleurInputSchema.xsd",
        Path(__file__).parent.parent / "build" / "FleurInputSchema.xsd",
        Path.cwd() / "FleurInputSchema.xsd",
        Path.cwd() / "build" / "FleurInputSchema.xsd",
    ]
    for path in possible_paths:
        if path.exists():
            return path.resolve()
    return None


@click.group(invoke_without_command=True)
@click.option('--schema', '-s', type=click.Path(exists=True), default=None,
              help='Path to the FleurInputSchema.xsd file.')
@click.option('--input', '-i', 'input_file', type=click.Path(), default=None,
              help='Path to inp.xml file (auto-detected in cwd if omitted).')
@click.version_option(version='0.1.0', prog_name='fleuriste')
@click.pass_context
def cli(ctx, schema, input_file):
    """FLEURiste - Terminal UI and CLI tools for FLEUR inp.xml editing."""
    ctx.ensure_object(dict)
    ctx.obj['schema'] = schema
    ctx.obj['input_file'] = input_file

    if ctx.invoked_subcommand is None:
        ctx.invoke(tui, schema=schema, input_file=input_file)


# ── TUI ──────────────────────────────────────────────────────────────────────

@cli.command()
@click.option('--schema', '-s', type=click.Path(exists=True), default=None,
              help='Path to the FleurInputSchema.xsd file.')
@click.option('--input', '-i', 'input_file', type=click.Path(), default=None,
              help='Path to inp.xml file.')
@click.pass_context
def tui(ctx, schema, input_file):
    """Launch the interactive terminal UI editor."""
    schema = schema or ctx.obj.get('schema')
    input_file = input_file or ctx.obj.get('input_file')

    # Auto-detect schema
    if schema is None:
        found = find_schema()
        if found:
            click.echo(f"Using schema: {found}")
            schema = str(found)
        else:
            click.echo("Warning: FleurInputSchema.xsd not found. "
                        "Running without schema validation.")

    from .app import run_editor
    run_editor(schema_path=schema, input_file=input_file)


# ── Generate ─────────────────────────────────────────────────────────────────

@cli.command()
@click.argument('input_source', type=click.Path(exists=True))
@click.option('--format', '-F', 'fmt', default=None,
              help='Structure file format for ASE (e.g., cif, vasp, xyz, pdb, extxyz). '
                   'If omitted, the input is treated as a FLEUR namelist for inpgen.')
@click.option('--profile', '-p', default='default_profile',
              help='Inpgen profile to use.')
@click.option('--nosym', is_flag=True, default=False,
              help='Disable symmetry detection.')
@click.option('--output', '-o', type=click.Path(), default=None,
              help='Output directory (default: current directory).')
@click.option('--force', '-f', is_flag=True, default=False,
              help='Overwrite existing inp.xml without prompting.')
@click.option('--film', is_flag=True, default=False,
              help='Treat structure as a film (only with --format).')
def generate(input_source, fmt, profile, nosym, output, force, film):
    """Generate inp.xml from an input or structure file via inpgen.

    INPUT_SOURCE is the path to the input file. Without --format it is
    treated as a FLEUR namelist.  With --format, ASE is used to read
    the structure and convert it to namelist input first.

    Examples:
      fleuriste generate input.txt
      fleuriste generate POSCAR -F vasp
      fleuriste generate struct.cif -F cif
      fleuriste generate struct.xyz -F xyz --film
    """
    try:
        from FleurInpgen import InpgenInterface
    except ImportError:
        raise click.ClickException(
            "FleurInpgen not available. Build FLEUR with Python bindings."
        )

    import os

    input_path = Path(input_source)
    output_dir = Path(output) if output else Path.cwd()
    output_dir.mkdir(parents=True, exist_ok=True)

    target = output_dir / "inp.xml"
    if target.exists() and not force:
        if not click.confirm(f"'{target}' already exists. Overwrite?"):
            raise SystemExit("Aborted.")

    # Read the input content
    if fmt is not None:
        # Use ASE to read the structure file, then convert to namelist
        try:
            from ase.io import read as ase_read
        except ImportError:
            raise click.ClickException(
                "ASE is required for reading structure files. "
                "Install it with: pip install ase"
            )
        from .inpgen_gui import ase_to_fleur_input

        ase_format = fmt if fmt != "auto" else None
        try:
            atoms = ase_read(str(input_path), format=ase_format)
        except Exception as e:
            raise click.ClickException(f"ASE could not read '{input_path}' (format={fmt}): {e}")

        content = ase_to_fleur_input(atoms, film=film)
        click.echo(f"Read structure via ASE (format={fmt}): "
                    f"{len(atoms)} atoms, {'film' if film else 'bulk'}")
    else:
        content = input_path.read_text()

    click.echo(f"Generating inp.xml from {input_path}")
    click.echo(f"  Profile : {profile}")
    click.echo(f"  No-sym  : {nosym}")

    original_dir = os.getcwd()
    try:
        os.chdir(str(output_dir))
        inpgen = InpgenInterface(quiet=True)
        inpgen.make_inp(content, profile, nosym)
        messages = inpgen.get_messages()
        if messages.strip():
            click.echo(f"\n{messages}")
    finally:
        os.chdir(original_dir)

    if target.exists():
        click.secho(f"✓ Generated {target}", fg='green')
    else:
        raise click.ClickException("inp.xml was not created – check input file.")


# ── K-Points ─────────────────────────────────────────────────────────────────

def _resolve_input(ctx) -> str:
    """Return the inp.xml path from context or cwd, or abort."""
    input_file = ctx.obj.get('input_file')
    if input_file is None:
        default = Path.cwd() / "inp.xml"
        if default.exists():
            return str(default)
        raise click.ClickException("No inp.xml found. Specify with -i/--input.")
    return input_file


@cli.group()
@click.option('--input', '-i', 'input_file', type=click.Path(exists=True),
              default=None, help='Path to inp.xml file.')
@click.pass_context
def kpoints(ctx, input_file):
    """Manage k-point lists in inp.xml."""
    ctx.ensure_object(dict)
    if input_file:
        ctx.obj['input_file'] = input_file
    # Resolve once so sub-commands don't have to
    ctx.obj['input_file'] = _resolve_input(ctx)


@kpoints.command('list')
@click.pass_context
def kpoints_list(ctx):
    """List all k-point sets defined in the inp.xml file."""
    from .kpoint_manager import InpXMLManager

    mgr = InpXMLManager(ctx.obj['input_file'], quiet=True)
    summary = mgr.get_summary()

    click.echo(f"File          : {summary['file']}")
    click.echo(f"K-point lists : {summary['num_lists']}")
    click.echo(f"Active list   : {summary['selected_list'] or '(none)'}")
    click.echo()
    for name, kplist in mgr.kpoint_lists.items():
        marker = " ← active" if name == summary['selected_list'] else ""
        click.echo(f"  {name} ({len(kplist.kpoints)} pts){marker}")


@kpoints.command('show')
@click.argument('name')
@click.pass_context
def kpoints_show(ctx, name):
    """Show details of a specific k-point list."""
    from .kpoint_manager import InpXMLManager

    mgr = InpXMLManager(ctx.obj['input_file'], quiet=True)
    kplist = mgr.get_kpoint_list(name)
    if kplist is None:
        avail = ', '.join(mgr.kpoint_lists.keys())
        raise click.ClickException(f"'{name}' not found. Available: {avail}")

    click.echo(f"K-point list : {name}")
    click.echo(f"Points       : {len(kplist.kpoints)}")
    click.echo()
    from .kpoint_gui import evaluate_coordinate

    click.echo(f"{'#':>5}  {'kx':>12}  {'ky':>12}  {'kz':>12}  {'weight':>10}  label")
    click.echo("─" * 70)
    for i, kp in enumerate(kplist.kpoints):
        click.echo(
            f"{i:>5}  {evaluate_coordinate(kp.x):>12.8f}  {evaluate_coordinate(kp.y):>12.8f}  "
            f"{evaluate_coordinate(kp.z):>12.8f}  {evaluate_coordinate(kp.weight):>10.6f}  {kp.label or ''}"
        )


@kpoints.command('select')
@click.argument('name')
@click.pass_context
def kpoints_select(ctx, name):
    """Set a k-point list as the active list."""
    from .kpoint_manager import InpXMLManager

    mgr = InpXMLManager(ctx.obj['input_file'], quiet=True)
    if not mgr.select_kpoint_list(name):
        avail = ', '.join(mgr.kpoint_lists.keys())
        raise click.ClickException(f"Cannot select '{name}'. Available: {avail}")
    mgr.save()
    click.secho(f"✓ '{name}' is now the active k-point list.", fg='green')


@kpoints.command('create')
@click.argument('name')
@click.option('--grid', '-g', nargs=3, type=int, default=None,
              help='Grid mode: nx ny nz dimensions.')
@click.option('--density', '-d', type=float, default=None,
              help='Density mode: k-point density value.')
@click.option('--number', '-n', type=int, default=None,
              help='Number mode: target number of k-points.')
@click.option('--path', '-p', type=str, default=None,
              help='Path mode: path with special points (e.g., "gamma=0,0,0;x=0.5,0,0.5;gamma=0,0,0").')
@click.option('--custom', '-c', type=str, default=None,
              help='Custom mode: direct inpgen string (e.g., "grid=10,10,6" or "gamma@den=0.03").')
@click.option('--gamma/--no-gamma', default=False,
              help='Use gamma-centred grid (grid/density/number modes).')
@click.option('--symmetry/--no-symmetry', default=True,
              help='Use symmetry to reduce k-points.')
@click.option('--points', type=int, default=100,
              help='Number of points for path mode (default: 100).')
@click.pass_context
def kpoints_create(ctx, name, grid, density, number, path, custom, gamma, symmetry, points):
    """Create k-points in various modes (specify exactly one mode).

    Modes:
      -g/--grid      : Monkhorst-Pack mesh (e.g., -g 12 12 8)
      -d/--density   : K-point density (e.g., -d 0.05)
      -n/--number    : Target number of k-points (e.g., -n 100)
      -p/--path      : Band structure path (e.g., -p "gamma=0,0,0;x=0.5,0,0.5")
      -c/--custom    : Direct inpgen string (e.g., -c "gamma@grid=10,10,6")

    Examples:
      fleuriste kpoints create fine -g 12 12 8 --gamma
      fleuriste kpoints create dense -d 0.05
      fleuriste kpoints create nk100 -n 100
      fleuriste kpoints create bands -p "gamma=0,0,0;x=0.5,0,0.5;l=0.5,0.5,0.5;gamma=0,0,0" --points 200
      fleuriste kpoints create custom -c "tria@den=0.03"
    """
    from .kpoint_manager import InpXMLManager, KPointMode, KPointModifiers
    
    # Count which modes are specified
    modes_specified = sum([
        grid is not None,
        density is not None,
        number is not None,
        path is not None,
        custom is not None,
    ])
    
    if modes_specified == 0:
        raise click.ClickException(
            "No mode specified. Use one of: -g/--grid, -d/--density, "
            "-n/--number, -p/--path, or -c/--custom"
        )
    
    if modes_specified > 1:
        raise click.ClickException(
            "Multiple modes specified. Use exactly one of: -g/--grid, -d/--density, "
            "-n/--number, -p/--path, or -c/--custom"
        )
    
    mgr = InpXMLManager(ctx.obj['input_file'], quiet=True)
    
    # Custom mode - direct inpgen string
    if custom is not None:
        if name in mgr.kpoint_lists:
            raise click.ClickException(
                f"K-point list '{name}' already exists. "
                f"Available: {', '.join(mgr.kpoint_lists.keys())}"
            )
        
        try:
            from FleurInpgen import InpgenInterface
        except ImportError:
            raise click.ClickException("FleurInpgen not available.")
        
        full_string = f"{name}#{custom}"
        nosym = not symmetry
        
        mgr.clear_messages()
        try:
            inpgen = InpgenInterface(quiet=True)
            inpgen.add_kpoints(full_string, kpts_path="", nosym=nosym)
            messages = inpgen.get_messages()
            if messages.strip():
                click.echo(messages)
        except Exception as e:
            raise click.ClickException(f"Inpgen error: {e}")
        
        mgr.load()
        mgr.save()
        
        kplist = mgr.get_kpoint_list(name)
        nk = len(kplist.kpoints) if kplist else "?"
        click.secho(f"✓ Created '{name}' ({custom}) → {nk} k-points", fg='green')
        return
    
    # Path mode - band structure
    if path is not None:
        mgr.create_kpoint_path(
            name=name,
            special_points=[],
            num_points=points,
            path_string=path,
        )
        mgr.save()
        
        kplist = mgr.get_kpoint_list(name)
        nk = len(kplist.kpoints) if kplist else "?"
        click.secho(f"✓ Created path '{name}' ({path}) → {nk} k-points", fg='green')
        return
    
    # Grid mode
    if grid is not None:
        mgr.create_kpoints(
            name=name,
            mode=KPointMode.GRID,
            modifiers=KPointModifiers(gamma=gamma),
            grid=tuple(grid),
            use_symmetry=symmetry,
        )
        mgr.save()
        
        kplist = mgr.get_kpoint_list(name)
        nk = len(kplist.kpoints) if kplist else "?"
        click.secho(
            f"✓ Created mesh '{name}' ({grid[0]}×{grid[1]}×{grid[2]}) → {nk} k-points",
            fg='green',
        )
        return
    
    # Density mode
    if density is not None:
        mgr.create_kpoints(
            name=name,
            mode=KPointMode.DENSITY,
            modifiers=KPointModifiers(gamma=gamma),
            density=density,
            use_symmetry=symmetry,
        )
        mgr.save()
        
        kplist = mgr.get_kpoint_list(name)
        nk = len(kplist.kpoints) if kplist else "?"
        click.secho(f"✓ Created density '{name}' (d={density}) → {nk} k-points", fg='green')
        return
    
    # Number mode
    if number is not None:
        mgr.create_kpoints(
            name=name,
            mode=KPointMode.NUMBER,
            modifiers=KPointModifiers(gamma=gamma),
            num_kpoints=number,
            use_symmetry=symmetry,
        )
        mgr.save()
        
        kplist = mgr.get_kpoint_list(name)
        nk = len(kplist.kpoints) if kplist else "?"
        click.secho(f"✓ Created '{name}' (target ~{number}) → {nk} k-points", fg='green')
        return


@kpoints.command('delete')
@click.argument('name')
@click.option('--force', '-f', is_flag=True, help='Skip confirmation.')
@click.pass_context
def kpoints_delete(ctx, name, force):
    """Delete a k-point list from inp.xml."""
    from .kpoint_manager import InpXMLManager

    mgr = InpXMLManager(ctx.obj['input_file'], quiet=True)

    if name not in mgr.kpoint_lists:
        avail = ', '.join(mgr.kpoint_lists.keys())
        raise click.ClickException(f"'{name}' not found. Available: {avail}")

    if name == mgr.get_selected_list():
        raise click.ClickException(f"Cannot delete '{name}' – it is the active list.")

    if not force:
        nk = len(mgr.get_kpoint_list(name).kpoints)
        if not click.confirm(f"Delete '{name}' ({nk} k-points)?"):
            raise SystemExit("Aborted.")

    mgr.remove_kpoint_list(name)
    mgr.save()
    click.secho(f"✓ Deleted '{name}'.", fg='green')


# ── Job ──────────────────────────────────────────────────────────────────────

@cli.group()
def job():
    """SLURM job script generation and FLEUR analysis tools."""
    pass


@job.command('analyze')
@click.argument('input_file', type=click.Path(exists=True), required=False,
                default=None)
def job_analyze(input_file):
    """Analyze a FLEUR inp.xml and print a summary."""
    from .pyjob.fleur_analyzer import FleurInputAnalyzer

    if input_file is None:
        default = Path.cwd() / "inp.xml"
        if default.exists():
            input_file = str(default)
        else:
            raise click.ClickException("No inp.xml found. Pass a file path.")

    analyzer = FleurInputAnalyzer(input_file)
    a = analyzer.analyze()

    elem_size = 16 if a.is_complex else 8
    mem_gb = a.matrix_dimension ** 2 * elem_size / (1024 ** 3)
    storage = "complex" if a.is_complex else "real"

    click.echo(f"File      : {input_file}")
    click.echo(f"Kmax      : {a.kmax:.2f} a.u.⁻¹")
    click.echo(f"Volume    : {a.cell_volume:.1f} Bohr³")
    click.echo(f"Spins     : {a.jspins}")
    click.echo(f"SOC       : {'yes' if a.has_soc else 'no'}")
    click.echo(f"Noco      : {'yes' if a.has_noco else 'no'}")
    click.echo(f"Inversion : {'yes' if a.has_inversion else 'no'}")
    click.echo(f"K-points  : {a.num_kpoints}")
    click.echo(f"Atoms     : {a.num_atoms} ({a.num_atom_types} types: "
               f"{', '.join(a.atom_species)})")
    click.echo(f"N_basis   : {a.n_basis_per_kpoint:,}")
    click.echo(f"Matrix    : {a.matrix_dimension:,}"
               f"{'  (×2 noco)' if a.has_noco else ''}")
    click.echo(f"Storage   : {storage}")
    click.echo(f"Memory/k  : {mem_gb:.3f} GB")


@job.command('generate')
@click.option('--input', '-i', 'input_file', type=click.Path(exists=True),
              default=None, help='Path to inp.xml for FLEUR analysis.')
@click.option('--machine', '-m', type=str, default=None,
              help='Machine name or path to JSON config file.')
@click.option('--partition', '-p', type=str, default=None,
              help='Partition name on the machine.')
@click.option('--name', '-n', 'job_name', default='fleur_job',
              help='SLURM job name.')
@click.option('--nodes', type=int, default=None,
              help='Number of nodes (auto-suggested if omitted).')
@click.option('--time', '-t', 'time_limit', default='01:00:00',
              help='Wall-time limit (HH:MM:SS).')
@click.option('--account', '-A', default=None, help='SLURM account.')
@click.option('--output', '-o', type=click.Path(), default=None,
              help='Output file (default: <job_name>.slurm).')
@click.option('--preview', is_flag=True, help='Print to stdout instead of saving.')
def job_generate(input_file, machine, partition, job_name, nodes, time_limit,
                 account, output, preview):
    """Generate a SLURM job script for a FLEUR calculation."""
    from .pyjob.slurm_generator import (
        MachineConfig, SlurmJobConfig, SlurmJobGenerator,
        load_machine_configs,
    )
    from .pyjob.fleur_analyzer import FleurInputAnalyzer
    from .pyjob.parallelization import suggest_parallelization

    # Auto-detect inp.xml
    if input_file is None:
        default = Path.cwd() / "inp.xml"
        if default.exists():
            input_file = str(default)

    # ── Machine config ────────────────────────────────────────────────
    machine_config = None
    if machine:
        mp = Path(machine)
        if mp.exists() and mp.suffix == '.json':
            machine_config = MachineConfig.load(str(mp))
        else:
            available = load_machine_configs()
            if machine in available:
                machine_config = available[machine]
            else:
                names = ', '.join(sorted(available.keys())) or '(none)'
                raise click.ClickException(
                    f"Machine '{machine}' not found. Available: {names}")
        click.echo(f"Machine  : {machine_config.name}")

    # ── Analyse FLEUR input ───────────────────────────────────────────
    fleur_result = None
    if input_file:
        try:
            analyzer = FleurInputAnalyzer(input_file)
            fleur_result = analyzer.analyze()
            click.echo(f"Input    : {Path(input_file).name}")
            click.echo(f"  Atoms={fleur_result.num_atoms}  "
                        f"K-pts={fleur_result.num_kpoints}  "
                        f"Matrix={fleur_result.matrix_dimension}")
        except Exception as exc:
            click.echo(f"Warning: cannot analyse input – {exc}", err=True)

    # ── Parallelisation suggestion ────────────────────────────────────
    par = None
    if machine_config and fleur_result:
        gpus_avail = (
            machine_config.get_effective_value("gpus_per_node", partition) or 0
        )
        par = suggest_parallelization(
            fleur_result=fleur_result,
            machine=machine_config,
            partition=partition,
            initial_nodes=nodes or 1,
            use_gpu=gpus_avail > 0,
        )
        click.echo(f"  → {par.num_nodes} nodes, {par.num_tasks} tasks, "
                    f"{par.cpus_per_task} cpus/task")

    # ── Build SLURM config ────────────────────────────────────────────
    num_nodes = nodes or (par.num_nodes if par else 1)
    ntasks    = par.num_tasks if par else num_nodes
    cpus      = par.cpus_per_task if par else 1
    memory    = (f"{int(par.memory_per_node_gb * 1.2) + 1}G"
                 if par else "4G")

    modules = []
    command = None
    gpu_syntax = "gpus"
    use_mem_option = True
    gpus = None

    if machine_config:
        modules = list(machine_config.modules_needed or [])
        command = machine_config.get_effective_value("command", partition)
        gpu_syntax = (
            machine_config.get_effective_value("gpu_syntax", partition)
            or "gpus"
        )
        use_mem = machine_config.get_effective_value("use_mem_option", partition)
        use_mem_option = use_mem if use_mem is not None else True

        if par and par.uses_gpu:
            gpn = (
                machine_config.get_effective_value("gpus_per_node", partition)
                or 0
            )
            gtype = (
                machine_config.get_effective_value("gpu_type", partition) or ""
            )
            if gpn > 0:
                gpus = f"{gtype}:{gpn}" if gtype else str(gpn)

    commands = [command] if command else []

    config = SlurmJobConfig(
        job_name=job_name,
        nodes=num_nodes,
        ntasks=ntasks,
        cpus_per_task=cpus,
        memory=memory,
        time=time_limit,
        partition=partition,
        account=account,
        gpus=gpus,
        gpu_syntax=gpu_syntax,
        use_mem_option=use_mem_option,
        modules=modules,
        commands=commands,
    )

    generator = SlurmJobGenerator(config, machine_config)
    script = generator.generate()

    if preview:
        click.echo()
        click.echo(script)
    else:
        out_file = output or f"{job_name}.slurm"
        generator.save(out_file)
        click.secho(f"✓ Saved job script: {out_file}", fg='green')


# ── Entry point ──────────────────────────────────────────────────────────────

def main():
    """Entry point for the CLI."""
    cli()
