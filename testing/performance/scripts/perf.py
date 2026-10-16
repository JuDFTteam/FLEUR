import json
import os
import sys
import shutil

timerlist={
    "uname": "",
    "mpi-tasks": 0,
    "omp-threads" : 0,
    "gpus"        : 0,
    "Total":0.0,
    "Initialization":0.0,
    "Iteration":0.0,
    "generation of potential":0.0,
    "H generation and diagonalizati":0.0,
    "Setup of H&S matrices":0.0,
    "Diagonalization":0.0,
    "generation of new charge densi":0.0,
    "Charge Density Mixing":0.0,
}

def process_subtimers(st):
    for i in range(len(st)):
        name=st[i]["timername"]
        if name in timerlist:
            timerlist[name]=st[i]["totaltime"]
        try:
            process_subtimers(st[i]["subtimers"])
        except:
            pass #No subtimers

def process_judft_times(dir):
    all_times=json.load(open(dir+"/juDFT_times.json","r"))
    #Get the description of the total run
    timerlist["uname"]=all_times["uname"]
    timerlist["mpi-tasks"]=all_times["mpi-tasks"]
    timerlist["omp-threads"]=all_times["omp-threads"]
    timerlist["gpus"]=all_times["gpus"]
    timerlist["Total"]=all_times["totaltime"]    
    #the selected subtimers
    process_subtimers(all_times["subtimers"])
    with open(dir+"/perf.json","w") as f:
        f.write(json.dumps(timerlist,indent=4))
    write_BMF(dir,timerlist)
    return dir+"/perf.json"    


def write_BMF(dir,timerlist):
    benchmarks=["Total","Initialization","Iteration","generation of potential","H generation and diagonalizati","Setup of H&S matrices",
        "Diagonalization","generation of new charge densi","Charge Density Mixing"]
    bmf={}

    for b in benchmarks:
        bmf[b]={"runtime":{"value":timerlist[b]}}
    with open(dir+"/bencher.json","w") as f:
        f.write(json.dumps(bmf,indent=4))



def run_fleur(testdir,env):
    import subprocess
    import calendar
    import time
    #check for executable
    if os.path.isfile("fleur"):
        fleur="fleur"
    elif os.path.isfile("fleur_MPI"):
        fleur="fleur_MPI"
    else:
        sys.exit("No FLEUR executable found")
    #create a directory
    current_GMT = time.gmtime()
    dir="Testing/performance/"+os.path.basename(os.path.normpath(testdir))+"/"+str(calendar.timegm(current_GMT))
    if not os.path.isdir(dir):
        os.makedirs(dir)
    else:
        sys.exit("Test already exists?")
    #copy inp.xml
    shutil.copy(testdir+"/inp.xml",dir)

    #run fleur
    cwd=os.getcwd()
    os.chdir(dir)
    if env:
        env={**os.environ,**env}
    else:
        env=os.environ.copy()
    if os.environ.get("juDFT_MPI"):
            mpi=os.environ.get("juDFT_MPI")
    else:
        mpi=""        
                
    subprocess.run([f"{mpi} {cwd}/{fleur}"],env=env,shell=True)  
    os.chdir(cwd)
     
    #postprocess
    return process_judft_times(dir)

def available_tests():
    """All directories in inputfiles containing an inp.xml"""
    inputdir=os.path.dirname(os.path.abspath(__file__))+"/../inputfiles"
    return sorted(d for d in os.listdir(inputdir) if os.path.isfile(f"{inputdir}/{d}/inp.xml"))

def run_test(name=None,env=None):
    #Get the directory of the script:
    scriptdir=os.path.dirname(__file__)+"/../inputfiles/"
    if (name):
        scriptdir=scriptdir+"/"+name
    else:
        #default test
        scriptdir=scriptdir+"/Noco"
    if not os.path.isfile(scriptdir+"/inp.xml"):
        sys.exit(f"Unknown performance test: {name}. Available: {' '.join(available_tests())}")
    return run_fleur(scriptdir,env)

if __name__=='__main__':

    if len(sys.argv)>1:
        if sys.argv[1]=="bencher":
            process_judft_times(os.getcwd())
            exit()
        if sys.argv[1]=="list":
            print("\n".join(available_tests()))
            exit()
        names=available_tests() if sys.argv[1]=="all" else sys.argv[1:]
        for name in names:
            print(run_test(name))
    else:
        print(run_test())

