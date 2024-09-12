import sys 
import re 

def end_sub(line):
    re_end_sub = re.compile(r".*end\s*subroutine.*", re.IGNORECASE)
    return re_end_sub.match(line) is not None

def sub(line):
    if("end" in line.lower()):
        return None
    re_end_sub = re.compile(r".*subroutine(.*)\(", re.IGNORECASE)
    m = re_end_sub.match(line)
    if(m):
        return m.group(1).strip()
    return None

def timestart(line):
    return "call timestart(" in line.lower()

def timestop(line):
    return "call timestop(" in line.lower()

def process_file(filename):

    subroutine = False
    name_sub = ""
    timers = 0
    with open(filename) as f:
        for line in f:
            if(line[0] != "c"):
                if(timestart(line)):
                    if(subroutine):
                        timers += 1 
                        tot_timer += 1
                    else:
                        print("start timer outside of subroutine")
                    
                if(timestop(line)):
                    if(subroutine):
                        timers -= 1 
                    else:
                        print("stop timer outside of subroutine")


                tmp = sub(line)
                if(tmp is not None):
                    subroutine = True 
                    name_sub = tmp
                    tot_timer = 0

                if(end_sub(line)):
                    if(subroutine):
                        subroutine = False 
                    else:
                        print("ended sub without starting it")

                    if(timers != 0):
                        print(f"{name_sub} isn't balanced")

                    if(tot_timer == 0):
                        print(f" tot_timer = {tot_timer} -> {name_sub}")

for filename in sys.argv[1:]:
    process_file(filename)