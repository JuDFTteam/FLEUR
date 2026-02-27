#!/usr/bin/env python3
"""
SLURM Job Submission File Generator

A Python application to generate SLURM batch job submission scripts
with configurable options for HPC cluster job scheduling.
"""

import json
import click
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Tuple
from datetime import datetime


@dataclass
class Partition:
    """
    Describes a SLURM partition/queue.
    
    Partition-level settings override machine defaults when specified (not None).
    
    Attributes:
        name: Partition name (required)
        max_nodes: Maximum nodes available in this partition
        max_runtime: Maximum walltime (HH:MM:SS or D-HH:MM:SS)
        cores_per_node: Override machine's cores_per_node (None = use machine default)
        gpus_per_node: Override machine's gpus_per_node (None = use machine default)
        memory_per_node_gb: Override machine's memory (None = use machine default)
        gpu_type: Override machine's gpu_type (None = use machine default)
        command: Command to run (None = use machine default)
        shell_commands: Additional shell commands to run before main execution (None = use machine default)
        default: Whether this is the default partition
        description: Human-readable description
    """
    
    name: str
    max_nodes: int = 1
    max_runtime: str = "24:00:00"  # HH:MM:SS or D-HH:MM:SS
    cores_per_node: Optional[int] = None  # None = use machine default
    gpus_per_node: Optional[int] = None  # None = use machine default
    memory_per_node_gb: Optional[int] = None  # None = use machine default
    gpu_type: Optional[str] = None  # None = use machine default
    gpu_syntax: Optional[str] = None  # "gpus" (--gpus=N) or "gres" (--gres=gpu:N)
    use_mem_option: Optional[bool] = None  # None = use machine default; whether to use --mem in SLURM
    command: Optional[str] = None  # Command to run (None = use machine default)
    shell_commands: Optional[list] = None  # Additional shell commands (None = use machine default)
    default: bool = False
    description: str = ""
    
    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict) -> "Partition":
        """Create from dictionary, handling key name variations."""
        # Handle alternate key names from JSON
        mapping = {
            "max_runime": "max_runtime",  # typo in example JSON
            "max_time": "max_runtime",
            "max_cpus_per_node": "cores_per_node",
            "max_gpus_per_node": "gpus_per_node",
            "max_memory_gb": "memory_per_node_gb",
        }
        normalized = {}
        for key, value in data.items():
            norm_key = mapping.get(key, key)
            if norm_key in cls.__dataclass_fields__:
                normalized[norm_key] = value
        return cls(**normalized)


@dataclass
class MachineConfig:
    """
    Describes the properties of an HPC cluster/machine.
    
    Machine-level settings serve as defaults that can be overridden
    at the partition level. When a partition specifies a value (not None),
    it takes precedence over the machine default.
    
    Example JSON structure:
    {
        "name": "gpu_cluster",
        "description": "Generic GPU compute cluster",
        "account": "default",
        "cores_per_node": 64,
        "gpus_per_node": 4,
        "memory_per_node_gb": 256,
        "gpu_type": "a100",
        "max_runtime": "24:00:00",
        "max_nodes": 16,
        "fleur_path": "/home/fleur_mpi",
        "modules_needed": ["cuda/12.0", "python/3.10"],
        "partitions": [
            {
                "name": "gpu",
                "max_nodes": 16,
                "max_runtime": "24:00:00",
                "cores_per_node": null,       // null = use machine default (64)
                "gpus_per_node": null,        // null = use machine default (4)
                "memory_per_node_gb": null,   // null = use machine default (256)
                "gpu_type": "a100",           // can override machine default
                "fleur_path": "/alt/fleur",   // can override machine default
                "description": "GPU partition"
            }
        ]
    }
    
    Attributes:
        name: Machine/cluster identifier (required)
        description: Human-readable description
        account: Default SLURM account for job submission
        cores_per_node: Default CPU cores per node
        gpus_per_node: Default GPUs per node (0 for CPU-only)
        memory_per_node_gb: Default RAM in GB per node
        gpu_type: Default GPU type (e.g., "a100", "v100")
        gpu_syntax: SLURM syntax for GPU requests: "gpus" (--gpus=N) or "gres" (--gres=gpu:N)
        use_mem_option: Whether to include --mem in SLURM directives (some clusters don't use it)
        max_runtime: Default maximum walltime (HH:MM:SS)
        max_nodes: Default maximum nodes per job
        partitions: List of available partitions (can override defaults)
        modules_needed: List of modules to load for FLEUR jobs
        command: Default command to run in job script
        shell_commands: Additional shell commands to run before main execution
    """
    
    # Machine identification
    name: str
    description: str = ""
    account: Optional[str] = None
    
    # Node hardware specs (defaults, can be overridden per partition)
    cores_per_node: int = 1
    gpus_per_node: int = 0
    memory_per_node_gb: int = 1
    gpu_type: Optional[str] = None
    gpu_syntax: str = "gpus"  # "gpus" (--gpus=N) or "gres" (--gres=gpu:N)
    use_mem_option: bool = True  # Whether to include --mem in SLURM directives
    
    # Default limits
    max_runtime: str = "24:00:00"
    max_nodes: int = 1
    
    # Partitions (can override machine defaults)
    partitions: list = field(default_factory=list)
    
    # Additional properties
    modules_needed: list = field(default_factory=list)
    command: Optional[str] = None  # Default command to run
    shell_commands: list = field(default_factory=list)  # Additional shell commands
    
    def __post_init__(self):
        """Convert partition dicts to Partition objects if needed."""
        self.partitions = [
            Partition.from_dict(p) if isinstance(p, dict) else p
            for p in self.partitions
        ]
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        data = {
            "name": self.name,
            "description": self.description,
            "account": self.account,
            "cores_per_node": self.cores_per_node,
            "gpus_per_node": self.gpus_per_node,
            "memory_per_node_gb": self.memory_per_node_gb,
            "gpu_type": self.gpu_type,
            "gpu_syntax": self.gpu_syntax,
            "use_mem_option": self.use_mem_option,
            "max_runtime": self.max_runtime,
            "max_nodes": self.max_nodes,
            "partitions": [p.to_dict() for p in self.partitions],
            "modules_needed": self.modules_needed,
            "command": self.command,
            "shell_commands": self.shell_commands,
        }
        return data
    
    @classmethod
    def from_dict(cls, data: dict) -> "MachineConfig":
        """Create from dictionary, handling key name variations."""
        # Handle alternate key names from JSON
        mapping = {
            "default_max_runtime": "max_runtime",
            "default_max_nodes": "max_nodes",
            "modules_available": "modules_needed",
        }
        normalized = {}
        for key, value in data.items():
            norm_key = mapping.get(key, key)
            if norm_key in cls.__dataclass_fields__:
                normalized[norm_key] = value
        return cls(**normalized)
    
    def get_effective_value(self, attr: str, partition_name: Optional[str] = None):
        """
        Get effective value for an attribute, considering partition override.
        
        If partition specifies a non-None value, use it; otherwise use machine default.
        """
        machine_value = getattr(self, attr, None)
        
        if partition_name:
            partition = self.get_partition(partition_name)
            if partition:
                part_value = getattr(partition, attr, None)
                if part_value is not None:
                    return part_value
        
        return machine_value
    
    def save(self, filepath: str) -> None:
        """Save machine configuration to a JSON file."""
        path = Path(filepath)
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)
        print(f"Machine configuration saved to: {filepath}")
    
    @classmethod
    def load(cls, filepath: str) -> "MachineConfig":
        """Load machine configuration from a JSON file."""
        path = Path(filepath)
        with open(path, "r") as f:
            data = json.load(f)
        return cls.from_dict(data)
    
    def get_partition(self, name: str) -> Optional[Partition]:
        """Get a partition by name."""
        for p in self.partitions:
            if p.name == name:
                return p
        return None
    
    def get_default_partition(self) -> Optional[Partition]:
        """Get the default partition."""
        for p in self.partitions:
            if p.default:
                return p
        return self.partitions[0] if self.partitions else None
    
    def list_partitions(self) -> list:
        """Get list of partition names."""
        return [p.name for p in self.partitions]
    
    def validate_resources(self, nodes: int = 1, cpus_per_task: int = 1,
                          gpus: int = 0, memory_gb: int = 0,
                          partition: Optional[str] = None) -> list:
        """
        Validate requested resources against machine limits.
        Returns a list of warning/error messages.
        """
        issues = []
        
        # Get effective values (partition overrides machine defaults)
        max_cores = self.get_effective_value("cores_per_node", partition)
        max_gpus = self.get_effective_value("gpus_per_node", partition)
        max_mem = self.get_effective_value("memory_per_node_gb", partition)
        
        # Get partition for node limits
        part = self.get_partition(partition) if partition else None
        max_nodes = part.max_nodes if part else self.max_nodes
        
        # Check cores
        if cpus_per_task > max_cores:
            issues.append(f"Warning: Requested {cpus_per_task} CPUs exceeds {max_cores} cores per node")
        
        # Check GPUs
        if gpus > max_gpus:
            issues.append(f"Warning: Requested {gpus} GPUs exceeds {max_gpus} GPUs per node")
        
        # Check memory
        if memory_gb > max_mem:
            issues.append(f"Warning: Requested {memory_gb}GB exceeds {max_mem}GB per node")
        
        # Check nodes
        if nodes > max_nodes:
            issues.append(f"Warning: Requested {nodes} nodes exceeds partition limit of {max_nodes}")
        
        return issues
    
    def get_info(self) -> str:
        """Get a formatted string with machine information."""
        lines = [
            f"Machine: {self.name}",
            f"Description: {self.description}" if self.description else "",
            "",
            "Hardware per Node (defaults):",
            f"  - Cores: {self.cores_per_node}",
            f"  - Memory: {self.memory_per_node_gb} GB",
            f"  - GPUs: {self.gpus_per_node}" + (f" ({self.gpu_type})" if self.gpu_type else ""),
            "",
            f"Max Runtime: {self.max_runtime}",
            f"Max Nodes: {self.max_nodes}",
        ]
        
        if self.partitions:
            lines.append("")
            lines.append("Partitions (can override defaults):")
            for p in self.partitions:
                default_marker = " (default)" if p.default else ""
                lines.append(f"  - {p.name}{default_marker}:")
                lines.append(f"      Max Nodes: {p.max_nodes}, Max Time: {p.max_runtime}")
                # Show overrides
                overrides = []
                if p.cores_per_node is not None:
                    overrides.append(f"cores={p.cores_per_node}")
                if p.gpus_per_node is not None:
                    overrides.append(f"gpus={p.gpus_per_node}")
                if p.memory_per_node_gb is not None:
                    overrides.append(f"mem={p.memory_per_node_gb}GB")
                if p.gpu_type is not None:
                    overrides.append(f"gpu_type={p.gpu_type}")
                if overrides:
                    lines.append(f"      Overrides: {', '.join(overrides)}")
                if p.description:
                    lines.append(f"      {p.description}")
        
        if self.modules_needed:
            lines.append("")
            lines.append(f"Modules: {', '.join(self.modules_needed)}")
        
        return "\n".join(line for line in lines if line is not None)
        
        return "\n".join(line for line in lines if line is not None)


def get_machines_directory() -> Path:
    """Get the path to the machines configuration directory."""
    return Path.home() / ".fleuriste" / "machines"


def load_machine_configs() -> dict:
    """
    Load machine configurations from ~/.fleuriste/machines/ directory.
    
    Scans the directory for JSON files and loads each as a MachineConfig.
    Returns a dictionary mapping machine names to MachineConfig objects.
    
    Returns:
        dict: Machine name -> MachineConfig mapping
    """
    configs = {}
    machines_dir = get_machines_directory()
    
    if not machines_dir.exists():
        # Create the directory if it doesn't exist
        machines_dir.mkdir(parents=True, exist_ok=True)
        return configs
    
    # Load all JSON files in the directory
    for json_file in machines_dir.glob("*.json"):
        try:
            machine = MachineConfig.load(str(json_file))
            configs[machine.name] = machine
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            print(f"Warning: Could not load machine config {json_file.name}: {e}")
            continue
    
    return configs


@dataclass
class SlurmJobConfig:
    """Configuration for a SLURM job submission."""
    
    # Required
    job_name: str
    
    # Resource allocation
    nodes: int = 1
    ntasks: int = 1
    cpus_per_task: int = 1
    memory: str = "4G"  # e.g., "4G", "8000M"
    time: str = "01:00:00"  # Format: HH:MM:SS or D-HH:MM:SS
    
    # Partition and QOS
    partition: Optional[str] = None
    qos: Optional[str] = None
    account: Optional[str] = None
    
    # Output files
    output_file: Optional[str] = None  # %j = job ID, %x = job name
    error_file: Optional[str] = None
    
    # GPU options
    gpus: Optional[str] = None  # e.g., "1", "a100:2"
    gpu_constraint: Optional[str] = None
    gpu_syntax: str = "gpus"  # "gpus" (--gpus=N) or "gres" (--gres=gpu:N)
    
    # Memory option
    use_mem_option: bool = True  # Whether to include --mem directive
    
    # Email notifications
    mail_user: Optional[str] = None
    mail_type: Optional[str] = None  # BEGIN, END, FAIL, ALL, NONE
    
    # Job dependencies
    dependency: Optional[str] = None  # e.g., "afterok:12345"
    
    # Array jobs
    array: Optional[str] = None  # e.g., "1-10", "1-100%5"
    
    # Working directory
    chdir: Optional[str] = None
    
    # Environment
    export: Optional[str] = None  # e.g., "ALL", "NONE", "VAR1,VAR2"
    
    # Constraints
    constraint: Optional[str] = None
    exclude: Optional[str] = None
    nodelist: Optional[str] = None
    
    # Additional SBATCH options
    extra_options: list = field(default_factory=list)
    
    # Script content
    modules: list = field(default_factory=list)  # Modules to load
    environment_setup: list = field(default_factory=list)  # Setup commands
    commands: list = field(default_factory=list)  # Main commands to run
    
    # Script settings
    shell: str = "/bin/bash"
    use_strict_mode: bool = True


class SlurmJobGenerator:
    """Generator for SLURM job submission scripts."""
    
    def __init__(self, config: SlurmJobConfig, machine: Optional[MachineConfig] = None):
        self.config = config
        self.machine = machine
        self._apply_machine_defaults()
    
    def _apply_machine_defaults(self) -> None:
        """Apply machine defaults and validate resources."""
        if not self.machine:
            return
        
        # If no partition specified, use default partition
        if not self.config.partition:
            default_part = self.machine.get_default_partition()
            if default_part:
                self.config.partition = default_part.name
        
        # Validate and warn about resource issues
        gpu_count = 0
        if self.config.gpus:
            # Parse GPU count (handles "1", "a100:2", etc.)
            gpu_str = self.config.gpus
            if ":" in gpu_str:
                gpu_count = int(gpu_str.split(":")[1])
            else:
                gpu_count = int(gpu_str)
        
        # Parse memory in GB
        mem_gb = 0
        mem_str = self.config.memory.upper()
        if mem_str.endswith("G"):
            mem_gb = int(mem_str[:-1])
        elif mem_str.endswith("M"):
            mem_gb = int(mem_str[:-1]) // 1024
        elif mem_str.endswith("T"):
            mem_gb = int(mem_str[:-1]) * 1024
        
        issues = self.machine.validate_resources(
            nodes=self.config.nodes,
            cpus_per_task=self.config.cpus_per_task,
            gpus=gpu_count,
            memory_gb=mem_gb,
            partition=self.config.partition,
        )
        
        for issue in issues:
            print(issue)
    
    def generate(self) -> str:
        """Generate the complete SLURM submission script."""
        lines = []
        
        # Shebang
        lines.append(f"#!{self.config.shell}")
        lines.append("")
        
        # Header comment
        lines.append(f"# SLURM Job Submission Script")
        lines.append(f"# Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append(f"# Job Name: {self.config.job_name}")
        if self.machine:
            lines.append(f"# Target Machine: {self.machine.name}")
        lines.append("")
        
        # SBATCH directives
        lines.append("#" + "=" * 60)
        lines.append("# SBATCH DIRECTIVES")
        lines.append("#" + "=" * 60)
        lines.append("")
        
        lines.extend(self._generate_sbatch_directives())
        lines.append("")
        
        # Environment setup
        if self.config.use_strict_mode:
            lines.append("#" + "=" * 60)
            lines.append("# SHELL OPTIONS")
            lines.append("#" + "=" * 60)
            lines.append("")
            lines.append("set -euo pipefail  # Exit on error, undefined vars, pipe failures")
            lines.append("")
        
        # Load modules
        if self.config.modules:
            lines.append("#" + "=" * 60)
            lines.append("# LOAD MODULES")
            lines.append("#" + "=" * 60)
            lines.append("")
            for module in self.config.modules:
                lines.append(f"module load {module}")
            lines.append("")
        
        # Environment setup commands
        if self.config.environment_setup:
            lines.append("#" + "=" * 60)
            lines.append("# ENVIRONMENT SETUP")
            lines.append("#" + "=" * 60)
            lines.append("")
            for cmd in self.config.environment_setup:
                lines.append(cmd)
            lines.append("")
        
        # Shell commands from machine config (machine + partition combined)
        if self.machine:
            partition = self.config.partition
            # Get machine-level shell_commands
            shell_commands = list(self.machine.shell_commands) if self.machine.shell_commands else []
            # Add partition-level shell_commands (if any)
            if partition:
                partition_obj = self.machine.get_partition(partition)
                if partition_obj and partition_obj.shell_commands:
                    shell_commands.extend(partition_obj.shell_commands)
            if shell_commands:
                lines.append("#" + "=" * 60)
                lines.append("# ADDITIONAL SHELL COMMANDS")
                lines.append("#" + "=" * 60)
                lines.append("")
                for cmd in shell_commands:
                    lines.append(cmd)
                lines.append("")
        
        # Main commands
        lines.append("#" + "=" * 60)
        lines.append("# MAIN EXECUTION")
        lines.append("#" + "=" * 60)
        lines.append("")
        
        if self.config.commands:
            for cmd in self.config.commands:
                lines.append(cmd)
        else:
            lines.append("# Add your commands here")
            lines.append("echo 'No commands specified'")
        
        lines.append("")
        
        return "\n".join(lines)
    
    def _generate_sbatch_directives(self) -> list:
        """Generate SBATCH directive lines."""
        directives = []
        c = self.config
        
        # Job name (required)
        directives.append(f"#SBATCH --job-name={c.job_name}")
        
        # Resource allocation
        directives.append(f"#SBATCH --nodes={c.nodes}")
        directives.append(f"#SBATCH --ntasks={c.ntasks}")
        directives.append(f"#SBATCH --cpus-per-task={c.cpus_per_task}")
        if c.use_mem_option:
            directives.append(f"#SBATCH --mem={c.memory}")
        directives.append(f"#SBATCH --time={c.time}")
        
        # Partition and account
        if c.partition:
            directives.append(f"#SBATCH --partition={c.partition}")
        if c.qos:
            directives.append(f"#SBATCH --qos={c.qos}")
        if c.account:
            directives.append(f"#SBATCH --account={c.account}")
        
        # Output files
        output = c.output_file or f"{c.job_name}_%j.out"
        directives.append(f"#SBATCH --output={output}")
        if c.error_file:
            directives.append(f"#SBATCH --error={c.error_file}")
        
        # GPU options
        if c.gpus:
            if c.gpu_syntax == "gres":
                # Legacy syntax: --gres=gpu:type:N or --gres=gpu:N
                if ":" in c.gpus:
                    # Format is "type:N" -> "gpu:type:N"
                    directives.append(f"#SBATCH --gres=gpu:{c.gpus}")
                else:
                    # Format is just "N" -> "gpu:N"
                    directives.append(f"#SBATCH --gres=gpu:{c.gpus}")
            else:
                # Modern syntax: --gpus=type:N or --gpus=N
                directives.append(f"#SBATCH --gpus={c.gpus}")
        if c.gpu_constraint:
            directives.append(f"#SBATCH --gpu-bind={c.gpu_constraint}")
        
        # Email notifications
        if c.mail_user:
            directives.append(f"#SBATCH --mail-user={c.mail_user}")
            mail_type = c.mail_type or "END,FAIL"
            directives.append(f"#SBATCH --mail-type={mail_type}")
        
        # Job dependencies
        if c.dependency:
            directives.append(f"#SBATCH --dependency={c.dependency}")
        
        # Array jobs
        if c.array:
            directives.append(f"#SBATCH --array={c.array}")
        
        # Working directory
        if c.chdir:
            directives.append(f"#SBATCH --chdir={c.chdir}")
        
        # Environment export
        if c.export:
            directives.append(f"#SBATCH --export={c.export}")
        
        # Node constraints
        if c.constraint:
            directives.append(f"#SBATCH --constraint={c.constraint}")
        if c.exclude:
            directives.append(f"#SBATCH --exclude={c.exclude}")
        if c.nodelist:
            directives.append(f"#SBATCH --nodelist={c.nodelist}")
        
        # Extra options
        for opt in c.extra_options:
            directives.append(f"#SBATCH {opt}")
        
        return directives
    
    def save(self, filepath: str) -> None:
        """Save the generated script to a file."""
        script = self.generate()
        path = Path(filepath)
        path.write_text(script)
        # Make executable
        path.chmod(path.stat().st_mode | 0o755)
        print(f"SLURM script saved to: {filepath}")


def interactive_mode(machine: Optional[MachineConfig] = None) -> SlurmJobConfig:
    """Interactively gather job configuration from user."""
    print("\n" + "=" * 60)
    print("SLURM Job Submission File Generator - Interactive Mode")
    print("=" * 60 + "\n")
    
    # Show machine info if available
    if machine:
        print(f"Target Machine: {machine.name}")
        print(f"  Cores/node: {machine.cores_per_node}, GPUs/node: {machine.gpus_per_node}, Memory/node: {machine.memory_per_node_gb}GB")
        if machine.partitions:
            print(f"  Partitions: {', '.join(machine.list_partitions())}")
        print()
    
    # Required: Job name
    job_name = input("Job name [my_job]: ").strip() or "my_job"
    
    # Resources with machine-aware defaults
    print("\n--- Resource Allocation ---")
    if machine:
        print(f"(Machine limits: {machine.cores_per_node} cores, {machine.memory_per_node_gb}GB memory per node)")
    
    nodes = int(input("Number of nodes [1]: ").strip() or "1")
    ntasks = int(input("Number of tasks [1]: ").strip() or "1")
    
    default_cpus = min(4, machine.cores_per_node) if machine else 1
    cpus_per_task = int(input(f"CPUs per task [{default_cpus}]: ").strip() or str(default_cpus))
    
    default_mem = "4G"
    memory = input(f"Memory (e.g., 4G, 8000M) [{default_mem}]: ").strip() or default_mem
    
    default_time = machine.default_max_runtime if machine else "01:00:00"
    time = input(f"Time limit (HH:MM:SS) [{default_time}]: ").strip() or default_time
    
    # Partition
    print("\n--- Partition & Account ---")
    if machine and machine.partitions:
        default_part = machine.get_default_partition()
        default_name = default_part.name if default_part else ""
        print(f"Available partitions: {', '.join(machine.list_partitions())}")
        partition = input(f"Partition [{default_name}]: ").strip() or default_name or None
    else:
        partition = input("Partition (leave empty for default): ").strip() or None
    
    account = input("Account (leave empty for default): ").strip() or None
    
    # GPU
    print("\n--- GPU Options ---")
    if machine and machine.gpus_per_node > 0:
        print(f"(Machine has {machine.gpus_per_node} {machine.gpu_type or 'GPU'}(s) per node)")
    gpus = input("GPUs (e.g., 1, a100:2, leave empty for none): ").strip() or None
    
    # Email
    print("\n--- Notifications ---")
    mail_user = input("Email for notifications (leave empty for none): ").strip() or None
    
    # Modules
    print("\n--- Modules ---")
    if machine and machine.modules_available:
        print(f"Common modules: {', '.join(machine.modules_available)}")
    modules_input = input("Modules to load (comma-separated, e.g., python/3.9,cuda/11.0): ").strip()
    modules = [m.strip() for m in modules_input.split(",")] if modules_input else []
    
    # Commands
    print("\n--- Commands ---")
    print("Enter commands to run (one per line, empty line to finish):")
    commands = []
    while True:
        cmd = input("> ")
        if not cmd:
            break
        commands.append(cmd)
    
    return SlurmJobConfig(
        job_name=job_name,
        nodes=nodes,
        ntasks=ntasks,
        cpus_per_task=cpus_per_task,
        memory=memory,
        time=time,
        partition=partition,
        account=account,
        gpus=gpus,
        mail_user=mail_user,
        modules=modules,
        commands=commands,
    )


@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
    """SLURM Job Submission File Generator.
    
    Generate SLURM batch job submission scripts with configurable options
    for HPC cluster job scheduling.
    
    \b
    Examples:
      # Interactive mode
      slurm-generator generate -i
      
      # Quick generation with common options
      slurm-generator generate -n my_job -t 02:00:00 --nodes 2 --mem 8G -c "python train.py"
      
      # GPU job with modules
      slurm-generator generate -n gpu_job --gpus 1 -m python/3.9 -m cuda/11.0 -c "python train.py"
      
      # Using a machine configuration file
      slurm-generator generate -n my_job --machine cluster.json -c "python train.py"
      
      # Create example machine config
      slurm-generator machine create gpu_cluster
      
      # Show machine info
      slurm-generator machine show cluster.json
    """
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())


@cli.group()
def machine():
    """Manage machine configurations."""
    pass


@machine.command("list")
def machine_list():
    """List available machine configurations from ~/.fleuriste/machines/."""
    machines_dir = get_machines_directory()
    configs = load_machine_configs()
    
    if not configs:
        click.echo(f"No machine configs found in: {machines_dir}")
        click.echo("Add .json files to this directory to configure machines.")
        return
    
    click.echo(f"Machine configs in: {machines_dir}")
    click.echo("")
    for name, machine_config in sorted(configs.items()):
        click.echo(f"  {name}: {machine_config.description or '(no description)'}")


@machine.command("show")
@click.argument("name_or_file")
def machine_show(name_or_file: str):
    """Display machine configuration.
    
    NAME_OR_FILE can be a machine name (from ~/.fleuriste/machines/) or a path to a JSON file.
    """
    # First try to load from available machines
    configs = load_machine_configs()
    if name_or_file in configs:
        click.echo(configs[name_or_file].get_info())
        return
    
    # Try as a file path
    try:
        machine_config = MachineConfig.load(name_or_file)
        click.echo(machine_config.get_info())
    except FileNotFoundError:
        raise click.ClickException(f"Machine '{name_or_file}' not found. Use 'machine list' to see available machines.")
    except json.JSONDecodeError:
        raise click.ClickException(f"Invalid JSON in machine config: {name_or_file}")


@cli.command()
@click.option("-i", "--interactive", is_flag=True,
              help="Run in interactive mode")
@click.option("-o", "--output", type=click.Path(), default=None,
              help="Output file path (default: <job_name>.slurm)")
@click.option("--machine", "machine_file", type=click.Path(exists=True), default=None,
              help="Path to machine config JSON file")
@click.option("-n", "--name", "job_name", default="my_job",
              help="Job name (default: my_job)")
@click.option("--account", default=None,
              help="Account/allocation name")
@click.option("--nodes", type=int, default=1,
              help="Number of nodes (default: 1)")
@click.option("--ntasks", type=int, default=1,
              help="Number of tasks (default: 1)")
@click.option("--cpus", type=int, default=1,
              help="CPUs per task (default: 1)")
@click.option("--mem", default="4G",
              help="Memory per node (default: 4G)")
@click.option("-t", "--time", "time_limit", default="01:00:00",
              help="Time limit HH:MM:SS (default: 01:00:00)")
@click.option("-p", "--partition", default=None,
              help="Partition/queue name")
@click.option("--gpus", default=None,
              help="Number/type of GPUs (e.g., 1, a100:2)")
@click.option("--mail", default=None,
              help="Email address for notifications")
@click.option("--mail-type", "mail_type", default="END,FAIL",
              help="When to send email (default: END,FAIL)")
@click.option("--array", default=None,
              help="Array job specification (e.g., 1-10)")
@click.option("--dependency", default=None,
              help="Job dependencies (e.g., afterok:12345)")
@click.option("-m", "--module", "modules", multiple=True,
              help="Module to load (can specify multiple times)")
@click.option("-c", "--command", "commands", multiple=True,
              help="Command to run (can specify multiple times)")
def generate(
    interactive: bool,
    output: Optional[str],
    machine_file: Optional[str],
    job_name: str,
    account: Optional[str],
    nodes: int,
    ntasks: int,
    cpus: int,
    mem: str,
    time_limit: str,
    partition: Optional[str],
    gpus: Optional[str],
    mail: Optional[str],
    mail_type: str,
    array: Optional[str],
    dependency: Optional[str],
    modules: Tuple[str, ...],
    commands: Tuple[str, ...],
):
    """Generate a SLURM job submission script."""
    # Load machine config if specified
    machine = None
    if machine_file:
        try:
            machine = MachineConfig.load(machine_file)
            click.echo(f"Using machine config: {machine.name}")
        except json.JSONDecodeError:
            raise click.ClickException(f"Invalid JSON in machine config: {machine_file}")
    
    if interactive:
        config = interactive_mode(machine)
    else:
        config = SlurmJobConfig(
            job_name=job_name,
            nodes=nodes,
            ntasks=ntasks,
            cpus_per_task=cpus,
            memory=mem,
            time=time_limit,
            partition=partition,
            account=account,
            gpus=gpus,
            mail_user=mail,
            mail_type=mail_type if mail else None,
            array=array,
            dependency=dependency,
            modules=list(modules),
            commands=list(commands),
        )
    
    # Generate script with machine config
    generator = SlurmJobGenerator(config, machine)
    script = generator.generate()
    
    # Output
    output_path = output or f"{config.job_name}.slurm"
    
    if output or interactive:
        generator.save(output_path)
    else:
        # Print to stdout if no output specified and not interactive
        click.echo(script)
        click.echo("\n" + "-" * 60)
        click.echo(f"To save to file, use: -o {output_path}")


def main():
    """Main entry point."""
    cli()


if __name__ == "__main__":
    main()
