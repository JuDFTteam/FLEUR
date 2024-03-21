import subprocess
import os
#Simple script to run several configurations and tests the performance

configs={
    "basic":["FC=gfortran",""],
    "ffast-math":["FC=gfortran","-flags -ffast-math"]
}

for conf in configs:
    print("Testing:",conf)
    cmd_arg=configs[conf][0]+" ./configure.sh -l "+conf+" "+configs[conf][1]
    #configure and make    
    print("Configure:",cmd_arg)
    subprocess.run(cmd_arg,shell=True)
    os.chdir("build."+conf)
    print("Make")
    subprocess.run("make",shell=True)
    print("Performance Test")
    subprocess.run("./run_tests.sh -perf",shell=True)
    os.chdir("..")
    