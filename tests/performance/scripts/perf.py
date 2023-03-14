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

def run_fleur(testdir):
    import subprocess
    import calendar
    import time
    #check for executable
    if os.path.isfile("fleur"): 
        fleur="fleur"
    elif os.path.isfile("fleur_MPI"): 
        fleur="fleur_MPI"
    else:
        print("No FLEUR executable found")
        sys.exit    
    #create a directory
    current_GMT = time.gmtime()
    dir="Testing/performance/"+str(calendar.timegm(current_GMT))
    if not os.path.isdir(dir):
        os.makedirs(dir)
    else:
        print("Test already exists?")
        sys.exit
    #copy inp.xml
    shutil.copy(testdir+"/inp.xml",dir)

    #run fleur
    cwd=os.getcwd()
    os.chdir(dir)
    subprocess.run([cwd+"/"+fleur])  
    os.chdir(cwd)
     
    #postprocess
    process_judft_times(dir)

def run_test(name=None):
    #Get the directory of the script:
    scriptdir=os.path.dirname(__file__)+"/../inputfiles/"
    if (name):
        scriptdir=scriptdir+"/"+name
    else:    
        #default test
        scriptdir=scriptdir+"/Noco"
    run_fleur(scriptdir)

if __name__=='__main__': run_test()

