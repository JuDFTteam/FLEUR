#!/usr/bin/env python3
"""
FLEUR Parallelization Strategy

Implements an algorithm to suggest optimal parallelization settings
for FLEUR calculations based on the input file analysis and machine
configuration.
"""

import math
from dataclasses import dataclass
from typing import Optional, List, Tuple

# Support both direct execution and package import
try:
    from .fleur_analyzer import FleurAnalysisResult
    from .slurm_generator import MachineConfig
except ImportError:
    from fleur_analyzer import FleurAnalysisResult
    from slurm_generator import MachineConfig


def get_factors(n: int) -> List[int]:
    """Get all factors of n in ascending order."""
    factors = []
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            factors.append(i)
            if i != n // i:
                factors.append(n // i)
    return sorted(factors)


def get_common_factors(a: int, b: int) -> List[int]:
    """Get all common factors of a and b in descending order."""
    factors_a = set(get_factors(a))
    factors_b = set(get_factors(b))
    common = factors_a.intersection(factors_b)
    return sorted(common, reverse=True)


@dataclass
class ParallelizationResult:
    """Result of parallelization strategy calculation."""
    
    # Recommended settings
    num_nodes: int
    num_tasks: int  # Total MPI tasks
    tasks_per_node: int
    tasks_per_gpu: Optional[int]  # None if not using GPUs
    cpus_per_task: int
    
    # Memory estimates
    memory_per_kpoint_gb: float
    memory_per_node_gb: float
    memory_per_gpu_gb: Optional[float]
    
    # Parallelization info
    kpoints_per_task: int
    uses_gpu: bool
    
    # Warnings/notes
    notes: List[str]
    
    def __str__(self) -> str:
        """Human-readable summary."""
        lines = [
            "=" * 50,
            "Parallelization Recommendation",
            "=" * 50,
            "",
            "--- Resource Allocation ---",
            f"  Nodes: {self.num_nodes}",
            f"  Total MPI tasks: {self.num_tasks}",
            f"  Tasks per node: {self.tasks_per_node}",
            f"  CPUs per task: {self.cpus_per_task}",
        ]
        
        if self.uses_gpu:
            lines.append(f"  Tasks per GPU: {self.tasks_per_gpu}")
        
        lines.extend([
            "",
            "--- K-point Distribution ---",
            f"  K-points per task: {self.kpoints_per_task}",
            "",
            "--- Memory Estimates ---",
            f"  Per k-point: {self.memory_per_kpoint_gb:.3f} GB",
            f"  Per node: {self.memory_per_node_gb:.2f} GB",
        ])
        
        if self.uses_gpu and self.memory_per_gpu_gb:
            lines.append(f"  Per GPU: {self.memory_per_gpu_gb:.2f} GB")
        
        if self.notes:
            lines.extend(["", "--- Notes ---"])
            for note in self.notes:
                lines.append(f"  • {note}")
        
        lines.append("=" * 50)
        return "\n".join(lines)


class ParallelizationStrategy:
    """
    Calculates optimal parallelization settings for FLEUR jobs.
    
    Algorithm:
    1. Start with k-point parallelization
    2. Try with initial number of nodes (default 1)
    3. Estimate memory per k-point as ~10 × matrix_size (conservative)
    4. Find common factors between num_kpoints and total_cores
    5. Try largest factor as number of tasks
    6. Calculate memory needed per node/GPU
    7. If memory exceeds limits, try smaller factor
    8. If even single task fails, increase number of nodes
    """
    
    # Memory multiplier: accounts for Hamiltonian, overlap, eigenvectors, work arrays
    MEMORY_MULTIPLIER = 10
    
    def __init__(
        self,
        fleur_result: FleurAnalysisResult,
        machine: MachineConfig,
        partition: Optional[str] = None,
        initial_nodes: int = 1,
        use_gpu: bool = False,
    ):
        """
        Initialize the parallelization strategy.
        
        Args:
            fleur_result: Analysis result from FleurInputAnalyzer
            machine: Machine configuration
            partition: Optional partition name (uses machine defaults if None)
            initial_nodes: Starting number of nodes to try
            use_gpu: Whether to use GPU parallelization
        """
        self.fleur = fleur_result
        self.machine = machine
        self.partition = partition
        self.initial_nodes = initial_nodes
        self.use_gpu = use_gpu and machine.gpus_per_node > 0
        
        # Get effective machine values (partition can override)
        self.cores_per_node = machine.get_effective_value("cores_per_node", partition)
        self.gpus_per_node = machine.get_effective_value("gpus_per_node", partition) or 0
        self.memory_per_node = machine.get_effective_value("memory_per_node_gb", partition)
        
        # Estimate memory per k-point
        # Matrix size = n_basis² × 8 bytes (real double)
        # Total memory ≈ MEMORY_MULTIPLIER × matrix_size
        matrix_bytes = self.fleur.n_basis_per_kpoint ** 2 * 8
        self.memory_per_kpoint_gb = (matrix_bytes * self.MEMORY_MULTIPLIER) / (1024**3)
        
        # For GPU, assume GPU memory is roughly 1/4 of node memory (typical)
        # This should be made configurable in MachineConfig
        self.gpu_memory_gb = self.memory_per_node / 4 if self.use_gpu else 0
    
    def calculate(self) -> ParallelizationResult:
        """
        Calculate the optimal parallelization strategy.
        
        Returns:
            ParallelizationResult with recommended settings
        """
        num_kpoints = self.fleur.num_kpoints
        notes = []
        
        # Start with initial number of nodes
        num_nodes = self.initial_nodes
        max_nodes = self.machine.max_nodes
        
        # Try to find a valid configuration
        result = None
        while num_nodes <= max_nodes:
            result = self._try_configuration(num_nodes, notes)
            if result is not None:
                break
            
            # Increase nodes
            num_nodes += 1
            notes.append(f"Increased to {num_nodes} nodes due to memory constraints")
        
        if result is None:
            # Fallback: use maximum resources
            notes.append("WARNING: Could not find optimal configuration, using fallback")
            result = self._fallback_configuration(max_nodes, notes)
        
        return result
    
    def _try_configuration(
        self, 
        num_nodes: int, 
        notes: List[str]
    ) -> Optional[ParallelizationResult]:
        """
        Try to find a valid configuration for the given number of nodes.
        
        Returns None if no valid configuration found.
        """
        num_kpoints = self.fleur.num_kpoints
        
        if self.use_gpu:
            total_gpus = num_nodes * self.gpus_per_node
            # Find common factors between k-points and total GPUs
            common = get_common_factors(num_kpoints, total_gpus)
        else:
            total_cores = num_nodes * self.cores_per_node
            # Find common factors between k-points and total cores
            common = get_common_factors(num_kpoints, total_cores)
        
        if not common:
            common = [1]
        
        # Try each factor from largest to smallest
        for num_tasks in common:
            if self.use_gpu:
                tasks_per_gpu = num_tasks // (num_nodes * self.gpus_per_node)
                if tasks_per_gpu < 1:
                    tasks_per_gpu = 1
                    # Recalculate actual tasks
                    num_tasks = tasks_per_gpu * num_nodes * self.gpus_per_node
                
                tasks_per_node = tasks_per_gpu * self.gpus_per_node
                
                # Memory check for GPU
                memory_per_gpu = tasks_per_gpu * self.memory_per_kpoint_gb
                if memory_per_gpu > self.gpu_memory_gb:
                    continue  # Try smaller factor
                
                memory_per_node = tasks_per_node * self.memory_per_kpoint_gb
                
            else:
                tasks_per_node = num_tasks // num_nodes
                if tasks_per_node < 1:
                    tasks_per_node = 1
                    num_tasks = tasks_per_node * num_nodes
                
                tasks_per_gpu = None
                memory_per_gpu = None
                
                # Memory check for node
                memory_per_node = tasks_per_node * self.memory_per_kpoint_gb
                if memory_per_node > self.memory_per_node:
                    continue  # Try smaller factor
            
            # Valid configuration found
            kpoints_per_task = num_kpoints // num_tasks
            if kpoints_per_task < 1:
                kpoints_per_task = 1
            
            # Calculate CPUs per task
            if self.use_gpu:
                cpus_per_task = self.cores_per_node // tasks_per_node
            else:
                cpus_per_task = max(1, self.cores_per_node // tasks_per_node)
            
            return ParallelizationResult(
                num_nodes=num_nodes,
                num_tasks=num_tasks,
                tasks_per_node=tasks_per_node,
                tasks_per_gpu=tasks_per_gpu,
                cpus_per_task=cpus_per_task,
                memory_per_kpoint_gb=self.memory_per_kpoint_gb,
                memory_per_node_gb=memory_per_node,
                memory_per_gpu_gb=memory_per_gpu,
                kpoints_per_task=kpoints_per_task,
                uses_gpu=self.use_gpu,
                notes=notes.copy(),
            )
        
        return None
    
    def _fallback_configuration(
        self, 
        num_nodes: int, 
        notes: List[str]
    ) -> ParallelizationResult:
        """Generate a fallback configuration when optimal search fails."""
        num_kpoints = self.fleur.num_kpoints
        
        if self.use_gpu:
            tasks_per_gpu = 1
            tasks_per_node = self.gpus_per_node
            num_tasks = num_nodes * tasks_per_node
            memory_per_gpu = self.memory_per_kpoint_gb
        else:
            # One task per node for memory-constrained cases
            tasks_per_node = 1
            num_tasks = num_nodes
            tasks_per_gpu = None
            memory_per_gpu = None
        
        memory_per_node = tasks_per_node * self.memory_per_kpoint_gb
        kpoints_per_task = max(1, num_kpoints // num_tasks)
        cpus_per_task = self.cores_per_node // tasks_per_node
        
        return ParallelizationResult(
            num_nodes=num_nodes,
            num_tasks=num_tasks,
            tasks_per_node=tasks_per_node,
            tasks_per_gpu=tasks_per_gpu,
            cpus_per_task=cpus_per_task,
            memory_per_kpoint_gb=self.memory_per_kpoint_gb,
            memory_per_node_gb=memory_per_node,
            memory_per_gpu_gb=memory_per_gpu,
            kpoints_per_task=kpoints_per_task,
            uses_gpu=self.use_gpu,
            notes=notes,
        )


def suggest_parallelization(
    fleur_result: FleurAnalysisResult,
    machine: MachineConfig,
    partition: Optional[str] = None,
    initial_nodes: int = 1,
    use_gpu: bool = False,
) -> ParallelizationResult:
    """
    Convenience function to calculate parallelization strategy.
    
    Args:
        fleur_result: Analysis result from FleurInputAnalyzer
        machine: Machine configuration
        partition: Optional partition name
        initial_nodes: Starting number of nodes to try
        use_gpu: Whether to use GPU parallelization
    
    Returns:
        ParallelizationResult with recommended settings
    """
    strategy = ParallelizationStrategy(
        fleur_result=fleur_result,
        machine=machine,
        partition=partition,
        initial_nodes=initial_nodes,
        use_gpu=use_gpu,
    )
    return strategy.calculate()


# CLI interface
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 3:
        print("Usage: python parallelization.py <inp.xml> <machine.json> [--gpu]")
        sys.exit(1)
    
    from fleur_analyzer import analyze_fleur_input
    
    inp_xml = sys.argv[1]
    machine_json = sys.argv[2]
    use_gpu = "--gpu" in sys.argv
    
    fleur_result = analyze_fleur_input(inp_xml)
    print(fleur_result)
    print()
    
    machine = MachineConfig.load(machine_json)
    
    result = suggest_parallelization(
        fleur_result=fleur_result,
        machine=machine,
        use_gpu=use_gpu,
    )
    print(result)
