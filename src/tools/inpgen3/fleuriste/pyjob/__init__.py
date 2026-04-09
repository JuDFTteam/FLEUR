"""
PyJob - SLURM Job Submission and FLEUR Analysis Tools.

Provides machine configuration, SLURM script generation,
FLEUR input analysis, and parallelization strategy.
"""

from .slurm_generator import (
    MachineConfig,
    Partition,
    SlurmJobConfig,
    SlurmJobGenerator,
    load_machine_configs,
    get_machines_directory,
)

from .fleur_analyzer import (
    FleurInputAnalyzer,
    FleurAnalysisResult,
    analyze_fleur_input,
)

from .parallelization import (
    ParallelizationStrategy,
    ParallelizationResult,
    suggest_parallelization,
)

from .slurm_machine_discovery import (
    SlurmDiscoveryError,
    discover_machine_config,
    write_discovered_machine_config,
)

__all__ = [
    "MachineConfig",
    "Partition",
    "SlurmJobConfig",
    "SlurmJobGenerator",
    "load_machine_configs",
    "get_machines_directory",
    "FleurInputAnalyzer",
    "FleurAnalysisResult",
    "analyze_fleur_input",
    "ParallelizationStrategy",
    "ParallelizationResult",
    "suggest_parallelization",
    "SlurmDiscoveryError",
    "discover_machine_config",
    "write_discovered_machine_config",
]
