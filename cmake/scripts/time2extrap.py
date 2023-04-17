levels=10
def print_timer(tt,path,parameter,scale,max_level=0):
    from json import dumps
    from copy import deepcopy
    name=tt["timername"]
    time=tt["totaltime"]
    line=deepcopy(parameter)
    line["callpath"]=f"{path}:{name}"
    line["value"]=time*scale
    print(dumps(line))
    #recursive call
    if (max_level<levels):
        if ("subtimers" in tt):
            for stt in tt["subtimers"]:
                print_timer(stt,f"{path}:{name}",parameter,scale,max_level+1)


def process_output(path):
    #first read out.xml to find number of MPI and number of OpenMP
    from xml.etree import ElementTree
    outxml=ElementTree.parse(f"{path}/out.xml")
    e2=outxml.findall(".//mpi")[0]
    mpi=int(e2.attrib["mpiProcesses"])
    e2=outxml.findall(".//openMP")[0]
    omp=int(e2.attrib["ompThreads"])
    parameter={"params":{"p":mpi,"o":omp},"metric":"time"}
    #process timer-files
    from json import load
    with open(f"{path}/juDFT_times.json") as f:
        tt=load(f)
    print_timer(tt,"",parameter,mpi*omp)

import sys
if __name__ == "__main__":
    process_output(sys.argv[1])

