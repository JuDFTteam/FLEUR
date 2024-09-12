import fleuristconf
import sys

def suggest(input,computer,nodes,efficiency):
    #First try to get the computer
    if computer:
        computernames=computer.split(",")
    else:
        computernames=["localhost"]
    machines=fleuristconf.get_computer(computernames)
    if not machines:
        print("ERROR: no computer found")
        sys.exit()
    
    #There might be several computers found:
    comp=[]
    res=[]
    for m in machines:
        if nodes:
            if m.maxnodes>nodes:
                m.maxnodes=0  #this machine is not available
            else:
                m.maxnodes=nodes
                m.minnodes=nodes
        #Find default config
        defaultconfig=m.find_best_config(input,efficiency=0.9)
        #Find possible less efficient config
        altconfig=m.find_best_config(input,efficiency=efficiency)
        if defaultconfig==altconfig:
            comp.append(m)
            res.append(defaultconfig)
        else:
            comp.append(m)
            res.append(defaultconfig)
            comp.append(m)
            res.append(altconfig)
    return comp,res

