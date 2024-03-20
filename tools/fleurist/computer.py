class computer:
    #basic properties
    name="default"
    maxnodes=1
    gpupernode=0 #not used so far
    corepernode=1
    fleurcommand=""       

    def __init__(self,conf):
        self.maxnodes=conf["params"]["maxnodes"]
        self.corepernode=conf["params"]["corepernode"]
        self.name=conf["name"]
        self.fleurcommand=conf["params"]["fleurcommand"]
     
    def submit(self,resources):
        pass

    def estimate_runtime(self,input, resources):
        """
        This routine should be overwritten
        """
        resources.time=0
        return resources.time
    
    def find_best_config(self,input,efficiency):
        from resources import resources
        import copy
        debug=0
        res=resources()
        #First try to find best performance on a single node
        min_time=1E99
        best_res=None
        for res.mpi in range(1,self.corepernode+1):
            if(debug): print(f"Checking on node with {res.mpi} mpi tasks")
            res.omp=self.corepernode//res.mpi
            runtime=self.estimate_runtime(input,res)
            if(debug): print(f"time:{runtime}")
            if runtime<min_time:
                min_time=runtime
                min_cost=res.cost()
                best_res=copy.deepcopy(res)
        if(debug): print(f"Best time: {min_time}")
        if(debug): print(best_res)        
        #Now try to use more nodes
        for res.nodes in range(2,self.maxnodes+1):
            for mpi in range(1,self.corepernode+1):
                if(debug): print(f"Checking {res.nodes} node with {mpi} mpi tasks")
                res.omp=self.corepernode//mpi      
                res.mpi=mpi*res.nodes
                runtime=self.estimate_runtime(input,res)
                if(debug): print(runtime,res.time)
                if res.time<best_res.time and (min_cost/res.cost()>efficiency):
                    best_res=copy.deepcopy(res)
        best_res.efficiency=min_cost/best_res.cost()            
        return best_res


    

class basic_cluster(computer):
    """
    Simple cluster
    Computing time estimated by a simple amdahl's law
    The serial part for the two parallel levels (mpi,omp) can be set seperately
    An additional constant can be added to the omp scaling to set a preference
    For an eigenvalue parallelization the same MPI scaling is used, but the problem size can be scaled
    """
    omp_fixed_serial=100
    mpi_fixed_serial=5
    mpi_fixed_serial_ev=100
    omp_parallel=0.7
    mpi_parallel=0.9
    mpi_parallel_ev=0.6
    ev_parallel=False

    def __init__(self,conf):
        if conf["type"] != "basic_cluster": 
            raise ValueError('Wrong type of computer')
        super().__init__(conf)
        if "omp_fixed_serial" in conf["params"]: self.omp_fixed_serial=conf["params"]["omp_fixed_serial"]
        if "mpi_fixed_serial" in conf["params"]: self.mpi_fixed_serial=conf["params"]["mpi_fixed_serial"]
        if "mpi_fixed_serial_ev" in conf["params"]: self.mpi_fixed_serials_ev=conf["params"]["mpi_fixed_serial_ev"]
        if "omp_parallel" in conf["params"]: self.omp_parallel=conf["params"]["omp_parallel"]
        if "mpi_parallel" in conf["params"]: self.mpi_parallel=conf["params"]["mpi_parallel"]
        if "mpi_parallel_ev" in conf["params"]: self.mpi_parallel_ev=conf["params"]["mpi_parallel_ev"]
        if "ev_parallel" in conf["params"]: self.ev_parallel=conf["params"]["ev_parallel"]
        
        
    def amdahl_scaling(self,workload,processes,parallel_percentage,fixed_seq):
        """
        create a scaling using amdahl's law. 
        """
        t_p=(parallel_percentage*workload)/processes
        t_s=(1.-parallel_percentage)*workload+fixed_seq
        scale=(fixed_seq+workload)/(t_s+t_p)
        return scale
    


    def estimate_runtime(self,input,resources):
        nkpt=input["nkpt"]
        basis_size=input["basis"]
        
        if nkpt%resources.mpi==0:
            #No EV-Parallelism needed
            s_mpi=self.amdahl_scaling(nkpt*basis_size**3,resources.mpi,self.mpi_parallel,self.mpi_fixed_serial)
            s_omp=self.amdahl_scaling(basis_size**3,resources.omp,self.omp_parallel,self.omp_fixed_serial)
            scale=s_omp*s_mpi
        else:
            if self.ev_parallel:
                import math
                mpi_k=math.gcd(nkpt,resources.mpi)
                mpi_ev=resources.mpi/mpi_k
                s_mpi=self.amdahl_scaling(nkpt*basis_size**3,mpi_k,self.mpi_parallel,self.mpi_fixed_serial)
                s_omp=self.amdahl_scaling(basis_size**3,resources.omp,self.omp_parallel,self.omp_fixed_serial)
                s_mpi_ev=self.amdahl_scaling(basis_size**3,mpi_ev,self.mpi_parallel_ev,self.mpi_fixed_serial_ev)
                scale=s_mpi*s_mpi_ev*s_omp
            else:
                scale=None    

        if scale:
            resources.time=basis_size**3/scale
        else:
            resources.time=1E99
        return resources.time
   

class slurm_cluster(basic_cluster):
    """ 
    Simple cluster with no k-point parallelism
    """
    slurm_begin=""
    slurm_end=""

    def __init__(self,conf):
        if conf["type"] != "basic_slurm_cluster": 
            raise ValueError('Wrong type of computer')
        conf["type"]="basic_cluster"
        super().__init__(conf)
        if "slurm_begin" in conf["params"]: self.slurm_begin=conf["params"]["slurm_begin"]
        if "slurm_end" in conf["params"]: self.slurm_begin=conf["params"]["slurm_end"]
                  

    def submit(self,resources,createonly=False):
        import subprocess

        with open("fleurist.job","w") as file:
            if self.slurm_begin: file.write(self.slurm_begin+"\n")
            file.write(resources.slurm_output())
            if self.slurm_end: file.write(self.slurm_end+"\n")
            
            file.write(f"\n\nsrun {self.fleurcommand}\n")

        if not createonly:
            subprocess.run(["sbatch","fleurist.job"])

class localhost(basic_cluster):
    """ 
    Determine parameters for localhost
    """
    def __init__(self,conf):
        import os
        if conf["type"] != "localhost": 
            raise ValueError('Wrong type of computer')

        self.corepernode=os.cpu_count()
        self.name="localhost"    
        self.name=conf["name"]
        self.fleurcommand=conf["params"]["fleurcommand"]

    def submit(self,resources,createonly=False):
        import os
        command=f"OMP_NUM_THREADS={resources.omp} mpirun -n {resources.mpi} {self.fleurcommand}"
        if createonly:
            print("You should execute:")
            print(command)
        else:
            print(f"Running: {command}")    
            os.run(command)
        

def get_computer_by_config(conf):
    classes=[slurm_cluster,localhost]
    for c in classes:
        try:
            return c(conf)
        except:
            pass
