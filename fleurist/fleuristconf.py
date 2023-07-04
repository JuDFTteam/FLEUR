import json
import click
import os
import sys
conffile=os.getenv("HOME")+"/.fleurist.conf"

default_configs=[
    {"name":"localhost",
     "type":"localhost",
     "params": {"maxnodes":1,"corepernode":1,"fleurcommand":None,"ev_parallel":None}},
    {"name":"iffslurm: th1-2022-64 (no eigenvalue parallelism)",
     "type":"basic_slurm_cluster",
     "params": {"maxnodes":10,"corepernode":64,"fleurcommand":None,"ev_parallel":False,"slurm_begin":"#SBATCH -p th1-2022-64"}},   
    {"name":"iffslurm: th1-2022-32 (no eigenvalue parallelism)",
     "type":"basic_slurm_cluster",
     "params":{"maxnodes":10,"corepernode":32,"fleurcommand":None,"ev_parallel":False,"slurm_begin":"#SBATCH -p th1-2022-32"}},    
    {"name":"jureca-dc-cpu (eigenvalue parallelism)",
     "type":"basic_slurm_cluster",
     "params": {"maxnodes":64,"corepernode":128,"fleurcommand":None,"ev_parallel":True,"slurm_begin":None}},    
    {"name":"Slurm cluster",
     "type":"basic_slurm_cluster",
     "params": {"maxnodes":10,"corepernode":32,"fleurcommand":None,"ev_parallel":None,"slurm_begin":None,"slurm_end":None}}    
]

parameter_descriptions={
     "ev_parallel": ["Can FLEUR be used for a distributed eigenvalue problem?","logical"],
     "maxnodes": ["Maximal number of compute nodes to use","integer"],
     "corepernode": ["Number of computational cores per node","integer"],
     "fleurcommand": ["The command to run FLEUR (could be just fleur_MPI or the full path if needed)","string"],
     "omp_fixed_serial":["omp_f_s","real"],
     "omp_parallel":["omp_p","real"],
     "mpi_fixed_serial":["mpi_f_s","real"],
     "mpi_parallel":["mpi_p","real"],
     "mpi_fixed_serial_ev":["mpi_f_s_ev","real"],
     "mpi_parallel_ev":["mpi_p_ev","real"],
     "slurm_begin":["Extra commands for header of slurm file, e.g. account or queue settings","string"],
     "slurm_end":["Extra commands after header of slurm file before calling FLEUR, e.g. modules to load","string"],
}
def load_conf(allow_fail=False):
    try:
        with open(conffile,"r") as file:
            conf=json.load(file)
        return conf
    except:
        if allow_fail:
             return []
        click.secho("FLEURist is not yet usable. Please use 'FLEURist computer add' to configure at least one computer",fg='red')
        sys.exit()


def get_computer(name):
    import computer

    conf=load_conf()
    machines=[]
    for c in conf:
        if c['name'] in name:
            machines.append(computer.get_computer_by_config(c))
    return machines


def add():
    import click
    import os
    click.secho("The following computers can be defined:",fg='red')
    for i in range(len(default_configs)):
             c=default_configs[i]
             click.secho(f"{i}: {c['name']}")
    click.secho("Please choose number of type or computer",fg='red')
    type_no=int(click.prompt("Type"))
    if not (type_no in range(len(default_configs))):
        import sys
        click.secho("Invalid type choosen")
        sys.exit()
 

    new_config=None
    import copy
    new_conf=copy.deepcopy(default_configs[type_no])
    if new_conf['name']!="localhost": new_conf['name']=click.prompt("Name of the new computer")
    for p in new_conf['params']:
        if new_conf['params'][p] in ["",0,None]:
            if p in parameter_descriptions:
                description=parameter_descriptions[p]
            else:
                click.secho(f"Unkown parameter {p}",fg="red")
                sys.exit()
            s_in=click.prompt(description[0])
            if description[1]=="integer":    
                new_conf['params'][p]=int(s_in)
            elif description[1]=="logical":
                new_conf['params'][p]=s_in.lower() in ['true', '1', 't', 'y', 'yes']
            elif description[1]=="real":
                new_conf['params'][p]=float(s_in)
            else:
                new_conf['params'][p]=s_in
            
        
    conf=load_conf(allow_fail=True)
    
    conf.append(new_conf)

    with open(conffile,"w") as file:
            json.dump(conf,file)

              
def remove(name):
    conf=load_conf()
    
    found=False
    for c in conf:
        if  c['name']==name:
             conf.remove(c)
             found=True
    
    if found:
        with open(conffile,"w") as file:
            json.dump(conf,file)
    else:
         click.secho(f"No computer with name '{name}' found",fg='red')

def list():
    import click
    conf=load_conf()
    
    for c in conf:
        click.secho(c["name"])