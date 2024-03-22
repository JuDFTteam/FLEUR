from perf import run_test
import os
import shutil
import subprocess
import sys
import json

def configure(source_dir,name,env,options):
    if os.path.isdir(f"build.{name}"): return 
    if (env):
        my_env={**os.environ,**env}
    else:
        my_env=os.environ.copy()

    if options:
        options=[f"{source_dir}/configure.sh","-l",f"{name}"]+options
    else:
        options=[f"{source_dir}/configure.sh","-l",f"{name}"]

    print(options)
    return subprocess.run(options,env=my_env)        

def make(name):
    workdir=os.getcwd()+f"/build.{name}"
    cpus=os.cpu_count()
    return subprocess.run(["make",f"-j{cpus}"],cwd=workdir)    

def run_perf(name,iter,env):
    cwd=os.getcwd()
    os.chdir(f"build.{name}")
    json_file=run_test(env=env)
    shutil.copyfile(json_file,f"{cwd}/{name}.{iter}.json")
    os.chdir(cwd)

def test_configs(filename):
    """Tests all configs defined in filename"""

    source_dir=os.path.dirname(os.path.realpath(__file__))+"/../../../"
    if os.path.isfile(filename):
        filename=filename
    elif os.path.isfile(f"{filename}.json"):
        filename=f"{filename}.json"
    elif os.path.isfile(f"{source_dir}tests/performance/{filename}"):
        filename=f"{source_dir}tests/performance/{filename}"
    elif os.path.isfile(f"{source_dir}tests/performance/{filename}.json"):
        filename=f"{source_dir}tests/performance/{filename}.json"
    else:
        print("Config-file not found")
        print(filename)
        print(f"{filename}.json")
        print(f"{source_dir}tests/performance/{filename}")
        print(f"{source_dir}tests/performance/{filename}.json")
        sys.exit()
    

    with open(filename,"r") as f:
        configs=json.load(f)

    for name in configs:
        if "env" in configs[name]: 
            env=configs[name]["env"]
        else:
            env=None
        
        if "opt" in configs[name]:
            opt=configs[name]["opt"]
        else:
            opt=None

        configure(source_dir,name,env,opt)
        make(name)
        i=0
        for timer_opts in configs[name]["timers"]:
            run_perf(f"{name}",i,timer_opts)
            i=i+1

if __name__ == "__main__":
    if len(sys.argv)<2:
        print("Error, please provide filename of the config file")
        sys.exit
    config_file=sys.argv[1]    
    test_configs(config_file)    

