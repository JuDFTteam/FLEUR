class resources:
    nodes=1
    mpi=1
    omp=1
    gpus=0
    time=0
    mem=0
    efficiency=1.0
    timestring="24:00:00"

    def description(self):
        return (
            f"Nodes       :{self.nodes}\n"
            f"Total MPI   :{self.mpi}\n"
            f"MPI per Node:{self.mpi//self.nodes}\n"
            f"OMP per MPI :{self.omp}\n"
            f"Est. Efficiency :{self.efficiency}\n"
        )
    def __init__(self):
        pass

    def __eq__(self, other): 
        if not isinstance(other, resources):
            return NotImplemented

        return (self.nodes==other.nodes and
                self.mpi==other.mpi     and
                self.mpi==other.mpi     and
                self.omp==other.omp     and
                self.gpus==other.gpus   and
                self.time==other.time   and
                self.mem==other.mem)
    
    def cost(self):
        return self.nodes*self.time
    
    def slurm_output(self):
        taskpernode=self.mpi//self.nodes
        slurm=(
         f"#SBATCH --nodes={self.nodes}\n"
         f"#SBATCH --ntasks={self.mpi}\n"
         f"#SBATCH --ntasks-per-node={taskpernode}\n"
         f"#SBATCH --cpus-per-task={self.omp}\n"
         f"#SBATCH --time={self.timestring}\n"
         f"export OMP_NUM_THREADS={self.omp}\n"
        )
        return slurm
         
    def aiida_resources(self):
        return {
            "num_machines":self.nodes,
            "tot_num_mpiprocs":self.mpi,
            "num_mpiprocs_per_machine":self.mpi//self.nodes,
            "num_cores_per_mpiproc":self.omp
        }