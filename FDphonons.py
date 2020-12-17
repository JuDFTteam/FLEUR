#!/usr/local/bin/python

# Script to automate finite-displacement calculations with FLEUR of FZ Juelich
# and PHON by Alfe.
# Uses only an input file for the inpgen of FLEUR.
# Call with python 3.5 upwards.
# Klueppelberg, Aug 2015

# Changing around a lot for a new Fleur version, new Cluster and new PHON.
# Neukirchen, Dec 2020

import subprocess
import sys
import os
import math
import time
import matplotlib
matplotlib.use('Agg') #Needed when operating the script from ssh with no XWindows. AN
import matplotlib.cbook as cbook
import matplotlib.pyplot as mpl
import decimal

# Define a function that properly converts decimal numbers to strings.
ctx = decimal.Context()
ctx.prec = 20

def float_to_str(f):
    """
    Convert the given float to a string,
    without resorting to scientific notation
    """
    d1 = ctx.create_decimal(repr(f))
    return format(d1, 'f')

os.system('source ~/.bashrc')

### Please specify parameters (all inputs required by user can be parsed with
### us_inp).

### us_inp: Name of (and path to) input file for inpgen.
inputfilename='./inp.sc'

### us_inp: Size of the supercell (repetitions of primitive unit cell in x-, y-
### and z-direction)

### Keeps the lattice system from the primitive unit cell to the supercell (i.e.
### if the primitive cell is cubic, the supercell has to be cubic as well).
scsizex=2
scsizey=2
scsizez=2

### us_inp: Path to FLEUR input-file generator as string
inpgen='/Users/neukirchen/fleurnov/build.mpintel/inpgen -f inp.sc'
### us_inp: Path to FLEUR executable as string (and call to queuing script)
FLEUR='sbatch ./scf.sh'
FORCERUN='sbatch ./force.sh'
fleursolo='/Users/neukirchen/fleurnov/build.mpintel/fleur_MPI'
### us_inp: Level of force implementation after Klueppelberg
### 0: Yu et al, 1: coretails, 2: kinetic energy surface, 3: disc potden)
level=3
folder='level'+str(level)

### us_inp: Path to PHON code
PHON='/Users/neukirchen/phon.1.47/src/phon'

### us_inp: Parameters for the cluster, i.e. a tag to name the calculation,
### the partition that is to be used, the number of nodes to request and how
### many cores are available per node.
calc_title='PHON'
partition='oscar'
node_count=1
#tasks_per_node=1
#cpus_per_task=12
proc_per_node=12

### us_inp: Convergence thresholds for the charge density/forces.
eps=0.00001
epsforce=0.00001

### Maximal number of available processors
maxproc=node_count*proc_per_node
### us_inp: Maximal number of scf-iteration cycles before script stops
### (failsafe threshold); maximal iteration count will be 15*fthreshold.
fsthreshold=5

### au to angstrom
autoangst=0.529177210903

### Some stuff needed later
##### Check if directory exists, if not, generate it
def ensure_dir(f):
  d = os.path.abspath(f)
  if not os.path.exists(d):
    os.mkdir(d)
    
### Step 1: Read input file for the input-file generator & generate INPHON and POSCAR files for PHON, which suggests displacements and generates supercell configuration
f=open(inputfilename)
data=f.read()
dataline=data.splitlines()

##### Are the coordinates in the input file in cartesian or internal coordinates?
##### The awk call selects the line containing the string 'cartesian' from the input file and selects the substring of length 1
##### starting two characters after 'cartesian' (i.e. the boolean defining if the coordinates are given in cartesian or internal coordinates)
cart=subprocess.Popen("awk '/cartesian/ {print substr($0,index($0,\"cartesian\")+length(\"cartesian\")+1,1)}' "+inputfilename,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
cart.wait()
##### cart contains the output of the awk call with a line break in stdout and possible error messages in stderr. use cartout and carterr as variables
cartout,carterr=cart.communicate()
##### remove the line break from the output
cartout=cartout.decode('ASCII').strip()
print("coordinates in",inputfilename,"are cartesian:",cartout)
if cartout!='t' and cartout!='T' and cartout!='f' and cartout!='F':
  sys.exit('unable to determine if cartesian or internal coordinates are used')

##### What kind of lattice is in place?
for line in range(len(dataline)):
##### Search the line containing the lattice system
  if 'latsys' in dataline[line]:
    entries=dataline[line].split()
##### Search the entry in the line containing the lattice system
    for entry in range(len(entries)):
      if 'latsys' in entries[entry]:
##### Read out the lattice system (format in input file is supposed to be "latsys='type'", no spaces in between)
        latsys=entries[entry][8:-1]
        break
    break
  if line == len(dataline):
    sys.exit('unable to determine lattice type, please use latsys specifier in '+inputfilename)
print("lattice system is",latsys)
##### use dictionary to find bravais matrix to convert cartesian to internal coordinates and vice versa

##### How many atoms of which type are present where? (list of atoms start with number of atoms, then list of according length containing atom coordinates)
atype=[]
apos=[]
xfac=1.0
yfac=1.0
zfac=1.0
for line in range(len(dataline)):
  entries=dataline[line].split()
##### Format of atom list in input file is #atoms in a separate line, followed by #atoms lines containing the charge number and the x, y, and z coordinates of the atom
##### Find line with one argument
  if len(entries) == 1:
    lastline=dataline[line+int(entries[0])].split()
##### Check if last line of potential atom list has the correct format and a potential charge number as first entry (118=ununoctium), 119 is excluded in range()
    if len(lastline) == 4 and int(lastline[0]) in range(119):
##### Read in atom types and coordinates
      for line2 in range(int(entries[0])):
        entries2=dataline[line+line2+1].split()
        atype.append(int(entries2[0]))
        apos.append([float(entries2[1]),float(entries2[2]),float(entries2[3])])
      break
##### modify entries if line &factor is present in input file
  if 'factor' in dataline[line]:
    xfac = float(entries[1])
    yfac = float(entries[2])
    zfac = float(entries[3])

for i in range(len(atype)):
  apos[i][0] = apos[i][0]/xfac
  apos[i][1] = apos[i][1]/yfac
  apos[i][2] = apos[i][2]/zfac

######### Here, some conversion in case of non-cubic lattices and regarding internal/external coordinates should be done!
######### Skipping that now and testing for cubic lattice and internal coordinates

#########################################################################################################################

##### Find numbers of equivalent atoms, this one is only conceptual.
##### count the number of repetitions of a particular atom type and define them to be the same species
equiv=[]
type=[]
count=0

if (len(atype)-1) == 0:
  equiv.append(1)
  type.append(atype[0])
for entry in range(len(atype)-1):
  count=count+1
  if atype[entry] != atype[entry+1]:
    equiv.append(count)
    type.append(atype[entry])
    if entry == len(atype)-2:
      equiv.append(1)
      type.append(atype[entry+1])
      break
    count=0
  if entry == len(atype)-2:
    equiv.append(count+1)
    type.append(atype[entry+1])
    break

print("Atom species:",str(type))
print("Number of atoms of each species:",str(equiv))

##### Read lattice parameters and scaling factor given in bohr (assumes that the format in the input file is a=number, without space)
for line in range(len(dataline)):
  if 'lattice' in dataline[line]:
    entries=dataline[line].split()
    for entry in range(len(entries)):
      if 'a=' in entries[entry]:
        latcon=entries[entry][2:]
      if 'a0=' in entries[entry]:
        latscl=entries[entry][3:]
##### TODO: noncubic lattices need parameters b, c, alpha, beta, and gamma to be read out

##### Calculate volume and Bravais matrix
BM=[]
if latsys=='sc':
  vol=((float(latcon)*float(latscl)*autoangst)**3)
  BM.append([float(latcon)*float(latscl)*autoangst,0,0])
  BM.append([0,float(latcon)*float(latscl)*autoangst,0])
  BM.append([0,0,float(latcon)*float(latscl)*autoangst])
if latsys=='fcc':
  vol=((float(latcon)*float(latscl)*autoangst)**3)/4
  BM.append([float(latcon)*float(latscl)*autoangst/2,float(latcon)*float(latscl)*autoangst/2,0])
  BM.append([float(latcon)*float(latscl)*autoangst/2,0,float(latcon)*float(latscl)*autoangst/2])
  BM.append([0,float(latcon)*float(latscl)*autoangst/2,float(latcon)*float(latscl)*autoangst/2])
if latsys=='bcc':
  vol=((float(latcon)*float(latscl)*autoangst)**3)*3/4
  BM.append([-float(latcon)*float(latscl)*autoangst/2,float(latcon)*float(latscl)*autoangst/2,float(latcon)*float(latscl)*autoangst/2])
  BM.append([float(latcon)*float(latscl)*autoangst/2,-float(latcon)*float(latscl)*autoangst/2,float(latcon)*float(latscl)*autoangst/2])
  BM.append([float(latcon)*float(latscl)*autoangst/2,float(latcon)*float(latscl)*autoangst/2,-float(latcon)*float(latscl)*autoangst/2])

ensure_dir('phonstart')

##### Prepare lines of POSCAR file
POSCAR=['None']*(7+len(atype))
POSCAR[0]=dataline[0]+' PHON POSCAR'
POSCAR[1]='-'+str(vol)
POSCAR[2]=str(BM[0][0])+' '+str(BM[0][1])+' '+str(BM[0][2])
POSCAR[3]=str(BM[1][0])+' '+str(BM[1][1])+' '+str(BM[1][2])
POSCAR[4]=str(BM[2][0])+' '+str(BM[2][1])+' '+str(BM[2][2])
POSCAR[5]=" ".join(str(entry) for entry in equiv)
POSCAR[6]='Direct'
for line in range(len(atype)):
  POSCAR[7+line]=str(apos[line][0])+' '+str(apos[line][1])+' '+str(apos[line][2])
##### Write POSCAR file in start-subdirectory
wf=open('./phonstart/POSCAR','w')
for line in range(len(POSCAR)):
  wf.write(POSCAR[line]+'\n')
wf.close()


atomic_mass = {  '1':   1.01,   '2':   4.00,   '3':   6.94,   '4':   9.01,   '5':  10.81,   '6':  12.01,
                 '7':  14.01,   '8':  16.00,   '9':  19.00,  '10':  20.18,  '11':  22.99,  '12':  24.31,
                '13':  26.98,  '14':  28.09,  '15':  30.97,  '16':  32.07,  '17':  35.45,  '18':  39.95,
                '19':  39.10,  '20':  40.08,  '21':  44.96,  '22':  47.87,  '23':  50.94,  '24':  52.00,
                '25':  54.94,  '26':  55.85,  '27':  58.93,  '28':  58.69,  '29':  63.55,  '30':  65.39,
                '31':  69.72,  '32':  72.61,  '33':  74.92,  '34':  78.96,  '35':  79.90,  '36':  83.80,
                '37':  85.47,  '38':  87.62,  '39':  88.91,  '40':  91.22,  '41':  92.91,  '42':  95.94,
                '43':  98.00,  '44': 101.07,  '45': 102.91,  '46': 106.42,  '47': 107.87,  '48': 112.41,
                '49': 114.82,  '50': 118.71,  '51': 121.76,  '52': 127.60,  '53': 126.90,  '54': 131.29,
                '55': 132.91,  '56': 137.33,  '57': 138.91,  '58': 140.12,  '59': 140.91,  '60': 144.24,
                '61': 145.00,  '62': 150.36,  '63': 151.96,  '64': 157.25,  '65': 158.93,  '66': 162.50,
                '67': 164.93,  '68': 167.26,  '69': 168.93,  '70': 173.04,  '71': 174.97,  '72': 178.49,
                '73': 180.95,  '74': 183.84,  '75': 186.21,  '76': 190.23,  '77': 192.22,  '78': 195.08,
                '79': 196.97,  '80': 200.59,  '81': 204.38,  '82': 207.20,  '83': 208.98,  '84': 209.00,
                '85': 210.00,  '86': 222.00,  '87': 223.00,  '88': 226.00,  '89': 227.00,  '90': 232.04,
                '91': 231.04,  '92': 238.03,  '93': 237.00,  '94': 244.00,  '95': 243.00,  '96': 247.00,
                '97': 247.00,  '98': 251.00,  '99': 252.00, '100': 257.00, '101': 258.00, '102': 259.00,
               '103': 262.00, '104': 261.00, '105': 262.00, '106': 266.00, '107': 264.00, '108': 269.00,
               '109': 268.00, '110': 271.00, '111': 272.00, '112': 285.00, '113': 284.00, '114': 289.00,
               '115': 288.00, '116': 292.00, '117': 293.00, '118': 294.00 }

##### Prepare lines of INPHON file
INPHON=['none']*9#*10 us_inp: increase to 10 and uncomment corresponding line 
                 # to set particular displacement length (PHON documentation).
                 # Else it will choose for you.
INPHON[0]='MASS = '+" ".join(str(atomic_mass[str(entry)]) for entry in type)
INPHON[1]='LSYMM = .T.'
INPHON[2]='LFORCEOUT = .T.'
INPHON[3]='LSUPER = .T.'
INPHON[4]='NDIM = '+str(scsizex)+' '+str(scsizey)+' '+str(scsizez)
INPHON[5]='NTYPES = '+str(len(type))
INPHON[6]='LCENTRAL = .T.'
INPHON[7]='IPRINT = 3'
INPHON[8]='LEIGEN = .T.'
#INPHON[9]='DISP = 25'
##### Write INPHON file in start-subdirectory
wf=open('./phonstart/INPHON','w')
for line in range(len(INPHON)):
  wf.write(INPHON[line]+'\n')
wf.close()


##### Use PHON code in new directory to generate displacement pattern
os.chdir('./phonstart/')
p=subprocess.Popen(PHON+' > output',shell=True)
p.wait()
os.chdir('../')

########## Could check if everything went well, assuming now that all files were created

##### Read out displacements from DISP file
f=open('./phonstart/DISP')
DISPdata=f.read()
DISP=DISPdata.splitlines()
##### dpat format: (atom, disp in x, disp in y, disp in z)
dpat=[]
for line in range(len(DISP)):
  entries=DISP[line].split()
  dpat.append([entries[1],entries[2],entries[3],entries[4]])

##### Read out coordinates in SPOSCAR file
f=open('./phonstart/SPOSCAR')
SPOSdata=f.read()
SPOS=SPOSdata.splitlines()[7:]
newpos=[]
for line in range(len(SPOS)):
  entries=SPOS[line].split()
  newpos.append([float(entries[0]),float(entries[1]),float(entries[2])])



### Step 2: Prepare displaced calculations and run them
##### Make sure that coordinates are given in internal coordinates

##### Prepare evaluation folder
ensure_dir('./phoneval/')
##### Prepare INPHON for evaluation

nf=open('./phoneval/FORCES','w')
nf.write(str(len(DISP))+'\n')
##### create new input files in different subdirectories
allforces=[len(DISP)]
newPOSCAR=[]
#for i in range(1):
isrunning=['T']*len(DISP)
failsafe=[0]*len(DISP)
for i in range(len(DISP)):
  ensure_dir(str(i)+'/')
  newlines=[]
  for line in range(len(dataline)):
    if '&' in dataline[line] and dataline[line].lstrip()[0] == '&':
      newline=''
      if 'cartesian' in dataline[line]:
        entries=dataline[line].split()
        for entry in entries:
          if 'cartesian' in entry:
            comps=entry.split("=")
            entry=comps[0]+'='+'f'
          newline=newline+entry+' '
        newlines.append(newline.strip())
        continue
      if 'lattice' in dataline[line]:
        entries=dataline[line].split()
        for entry in entries:
          if 'a=' in entry:
            comps=entry.split("=")
            entry=comps[0]+'='+str(scsizex*float(comps[1]))
          newline=newline+entry+' '
        newlines.append(newline.strip())
        continue
      if 'kpt' in dataline[line]:
        entries=dataline[line].split()
        for entry in entries:
          if 'div1' in entry:
            comps=entry.split("=")
            entry=comps[0]+'='+str(int(int(comps[1])/scsizex))
          if 'div2' in entry:
            comps=entry.split("=")
            entry=comps[0]+'='+str(int(int(comps[1])/scsizey))
          if 'div3' in entry:
            comps=entry.split("=")
            entry=comps[0]+'='+str(int(int(comps[1])/scsizez))
          newline=newline+entry+' '
        newlines.append(newline.strip())
        continue
      if not 'factor' in dataline[line]:
        newlines.append(dataline[line])
      continue
##### Block for atoms
    if len(dataline[line].split()) == 1 and int(dataline[line].strip()) == int(sum(int(j) for j in dataline[line].split())):
      newlines.append(str(scsizex*scsizey*scsizez*len(atype)).rjust(4))
      for line2 in range(len(newpos)):
        anum=atype[int(math.floor(line2/(scsizex*scsizey*scsizez)))]
        if line2+1 == int(dpat[i][0]):
          newlines.append(str(anum).rjust(4)+' '+str(float(newpos[line2][0])+float(dpat[i][1]))+' '+str(float(newpos[line2][1])+float(dpat[i][2]))+' '+str(float(newpos[line2][2])+float(dpat[i][3])))
        else:
          newlines.append(str(anum).rjust(4)+' '+str(float(newpos[line2][0]))+' '+str(float(newpos[line2][1]))+' '+str(float(newpos[line2][2])))
      continue
    if line == 0 or len(dataline[line].split()) < 3:
      newlines.append(dataline[line])
##### Write the new input file
  wf=open('./'+str(i)+'/inp.sc','w')
  for line in range(len(newlines)):
    wf.write(str(newlines[line])+'\n')
  wf.close()

##### generate inp file for FLEUR
  os.chdir('./'+str(i)+'/')
  p=subprocess.Popen(inpgen,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  p.wait()
##### When qsubscript is used, set good number of processors, depending on the number of kpts
##### Good number: largest divisor of nkpts up to 12 or largest number i up to 12 that has not nkpts%i == +-1 (12 processors per node)
##### Find nkpts TODO: This should really be checked and improved upon.
  file=open('./kpts.xml')
  f=file.read()
  nkptsline=f.splitlines()[1]
  tetis=nkptsline.split()
  nkpts=int(nkptsline.split()[2][7:-1])
  file.close()
##### Find good number of processors
  nproc=1
  for v in [min(maxproc,nkpts)-w for w in range(min(maxproc,nkpts))]:
    if nkpts % v == 1 or nkpts % v == v-1:
      continue
    if not int(round(float(nkpts)/v)) > int(round(float(nkpts)/nproc)):
      nproc=v
      if nkpts % v == 0:
        break
  print(str(nkpts),"k-points on",str(nproc),"processors in calculation",str(i))
##### construct qsubscript from sample
  ###TODO: this could really be handled more nicely.
  file=open('./scf.sh','w')
  file.write('#!/bin/sh'+'\n')
  file.write('#SBATCH -J '+calc_title+'\n')
  file.write('#SBATCH --nodes='+str(node_count)+'\n')
  #file.write('#SBATCH --ntasks-per-node='+tasks_per_node+'\n')
  file.write('#SBATCH --cpus-per-task='+str(nproc)+'\n')
  file.write('#SBATCH --time=12:00:00'+'\n')
  file.write('#SBATCH -p '+partition+'\n')
  current_dir=os.path.abspath('.')
  file.write('cd '+current_dir+'\n')
  file.write('ulimit -s unlimited'+'\n')
  file.write('mkdir ~/phon_temp'+str(i)+'\n')
  file.write('mkdir '+current_dir+'/'+folder+'\n')
  file.write('cp '+current_dir+'/* ~/phon_temp'+str(i)+'\n')
  file.write('cd ~/phon_temp'+str(i)+'\n') 
  file.write('a=1.0'+'\n')
  file.write('eps='+float_to_str(eps)+'\n')
  file.write('while [ true ]'+'\n')
  file.write('do'+'\n')
  file.write('   srun '+fleursolo+'\n') 
  file.write('   grestr=$(grep "distance of charge densities for spin  1" out | tail -1)'+'\n') 
  file.write('   grestr2=$(echo $grestr | sed "s/.*://")'+'\n') 
  file.write('   a=$(echo $grestr2 | grep -Eo "[0-9]+.([0-9]+)?" | tr "\n" " ")'+'\n')
  file.write('   echo $a'+'\n')
  file.write('   if [ 1 -eq "$(echo "${a} < ${eps}" | bc)" ]'+'\n')
  file.write('   then'+'\n')
  file.write('      break'+'\n')
  file.write('   fi'+'\n')
  file.write('done'+'\n')
  file.write('touch "converged"'+'\n')
  file.write('yes | rm mix* cdn[0234567]* cdn1? Check* coretail* eig* fort* pot* qpw tmat stars wkf2'+'\n') 
  file.write('cp -r ~/phon_temp'+str(i)+'/* '+current_dir+'/'+folder+'\n') 
  file.write(current_dir+'/'+folder+'\n')
  file.write('yes | rm -r ~/phon_temp'+str(i)+'\n')
  file.close()
  
  file=open('./force.sh','w')
  file.write('#!/bin/sh'+'\n')
  file.write('#SBATCH -J '+calc_title+'\n')
  file.write('#SBATCH --nodes='+str(node_count)+'\n')
  #file.write('#SBATCH --ntasks-per-node='+tasks_per_node+'\n')
  file.write('#SBATCH --cpus-per-task='+str(nproc)+'\n')
  file.write('#SBATCH --time=12:00:00'+'\n')
  file.write('#SBATCH -p '+partition+'\n')
  current_dir=os.path.abspath('.')
  file.write('cd '+current_dir+'\n')
  file.write('ulimit -s unlimited'+'\n')
  file.write('mkdir ~/phon_temp'+str(i)+'\n')
  file.write('mkdir '+current_dir+'/'+folder+'\n')
  file.write('cp '+current_dir+'/* ~/phon_temp'+str(i)+'\n')
  file.write('cd ~/phon_temp'+str(i)+'\n') 
  file.write('i=0'+'\n')  
  file.write('j=1'+'\n')  
  file.write('a=1.0'+'\n')
  file.write('eps='+float_to_str(epsforce)+'\n')
  file.write('while [ true ]'+'\n')
  file.write('do'+'\n')
  file.write('   srun '+fleursolo+'\n') 
  file.write('   grestr=$(grep "max force distance" out | tail -1)'+'\n') 
  file.write('   grestr2=$(echo $grestr | sed "s/.*distance=//")'+'\n') 
  file.write('   a=$(echo $grestr2 | grep -Eo "[0-9]+.([0-9]+)?" | tr "\n" " ")'+'\n')
  file.write('   echo $a'+'\n')
  file.write('   if [ 1 -eq "$(echo "${a} < ${eps}" | bc)" ]'+'\n')
  file.write('   then'+'\n')
  file.write('      break'+'\n')
  file.write('   fi'+'\n')
  file.write('done'+'\n')
  file.write('touch "converged"'+'\n')
  file.write('yes | rm mix* cdn[0234567]* cdn1? Check* coretail* eig* fort* pot* qpw tmat stars wkf2'+'\n') 
  file.write('cp -r ~/phon_temp'+str(i)+'/* '+current_dir+'/'+folder+'\n') 
  file.write(current_dir+'/'+folder+'\n')
  file.write('yes | rm -r ~/phon_temp'+str(i)+'\n')
  file.close()
##### Submit jobscript of FLEUR
  p=subprocess.Popen(FLEUR,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  p.wait()
  os.chdir('../')
  
conv = [False for i in range(len(DISP))]

##### Check every minute if all calculations are converged yet. If one is not, restart it.
while (all(conv) != True):
  time.sleep(60)
  for i in range(len(DISP)):
    os.chdir('./'+str(i)+'/')
    if os.path.isfile('./'+folder+'/converged'):
      p=subprocess.Popen('cp ./'+folder+'/* .',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      p.wait()
      p=subprocess.Popen('rm -r ./'+folder+'/',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      p.wait()
      p=subprocess.Popen('rm converged',shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
      p.wait()
      conv[i] = True
    os.chdir('../')

##### Switch on force calculation
for i in range(len(DISP)):
  os.chdir('./'+str(i)+'/')
  f=open('./inp.xml')
  fdata=f.read()
  f.close()
  f=open('./inp.xml','w')
  flines=fdata.splitlines()
  for line in flines:
    if 'l_f="F"' in line:
      comps=line.split('l_f="F"')
      line=comps[0]+'l_f="T" f_level="'+str(level)+'"'+comps[1]
    f.write(line+'\n')
  f.close()
##### Do force calculation
  p=subprocess.Popen(FLEUR,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  p.wait()
  isrunning[i]='T'
  os.chdir('../')

##### Check if force calculation is done for each displacement
while any(i == 'T' for i in isrunning):
  time.sleep(60)
  for i in range(len(DISP)):
    os.chdir('./'+str(i)+'/')
    if os.path.isfile('./'+folder+'/converged'):
      isrunning[i]='F'
    os.chdir('../')

### Step 3: Sort and collect FORCES files in the order of SPOSCAR file generated in Step 1, prepare second PHON folder and run PHON
for i in range(len(DISP)):
  os.chdir('./'+str(i)+'/')
  allforces.append(dpat[i][0]+' '+dpat[i][1]+' '+dpat[i][2]+' '+dpat[i][3])
##### epsilon should be just a little larger than the displacement
  epsilon=float(dpat[i][1])**2+float(dpat[i][2])**2+float(dpat[i][3])**2+0.00001
##### read original positions
  file = open('../phonstart/SPOSCAR')
  pright=file.read()

##### read positions in nth calculation
  file = open('./'+folder+'/POSCAR')
  p1=file.read()

##### read forces from nth calculation
  file = open('./'+folder+'/FORCES')
  f1=file.read()
  f1new=['None']*(len(f1.splitlines()[1:]))
  f1new[0]=dpat[i][0]+' '+dpat[i][1]+' '+dpat[i][2]+' '+dpat[i][3]
  force1=f1.splitlines()[2:]

##### read lines from file data (beginning with first position line)
  if i == 0:
    headright=pright.splitlines()[0:7]
    head1=p1.splitlines()[0:7]
  dataright=pright.splitlines()[7:]
  data1=p1.splitlines()[7:]
##### prepare check for general shift/offset
  [offx,offy,offz]=dataright[0].split()
  [offx,offy,offz]=[float(offx),float(offy),float(offz)]
##### for each a,b,c coordinate in nth calculation file...
  for j in range(len(data1)):
    [a,b,c]=data1[j].split()
    [a,b,c]=[float(a),float(b),float(c)]
##### ...correct any offsets compared to original positions...
    if j == 0:
      if (offx-a)**2+(offy-b)**2+(offz-c)**2 > epsilon:
        [offx,offy,offz]=[offx-a,offy-b,offz-c]
    [a,b,c]=[a+offx,b+offy,c+offz]
##### ...shift the coordinates to be positive...
    while a < -math.sqrt(epsilon):
      a=a+1
    while b < -math.sqrt(epsilon):
      b=b+1
    while c < -math.sqrt(epsilon):
      c=c+1
    found='false'
##### ...and compare them to original file
    for l in range(len(dataright)):
      [x,y,z]=dataright[l].split()
      [x,y,z]=[float(x),float(y),float(z)]
##### if line in original file is found...
      if (x-a)**2+(y-b)**2+(z-c)**2 < epsilon:
        found='true'
##### ...sort forces of nth calculation to that line
        allforces.append(force1[j])
        f1new[l+1]=force1[j]
        break
    if found == 'false':
      print("not found")
      print("DISP",str(i),"pos",str(j))
      print(a,b,c)
##### write force file with new order once all is done
  for line in range(len(f1new)):
    nf.write(str(f1new[line])+'\n')
##### prepare POSCAR file for phonon calculation
  if i == 0:
    for j in range(len(headright)):
      if j > 0 and j < 5:
        newPOSCAR.append(head1[j])
      else:
        newPOSCAR.append(headright[j])
    for j in range(len(dataright)):
      newPOSCAR.append(dataright[j])
    file = open('../phoneval/POSCAR','w')
    for line in newPOSCAR:
      file.write(str(line)+'\n')
    file.close
  os.chdir('../')
nf.close()


BZpaths={ # 'sc': 'ND = 4; NPOINTS = 100\n'+
          #       'QI = 0.0 0.0 0.0   0.0 0.0 0.5   0.5 0.5 0.5   0.0 0.0 0.0\n'+
          #       'QF = 0.0 0.0 0.5   0.5 0.5 0.5   0.0 0.0 0.0   0.0 0.5 0.5',
          # 'sc': 'ND = 5; NPOINTS = 100\n'+
          #       'QI = 0.0 0.0 0.0   0.0 0.0 0.5   0.0 0.5 0.5   0.0 0.0 0.0   0.5 0.5 0.5\n'+
          #       'QF = 0.0 0.0 0.5   0.0 0.5 0.5   0.0 0.0 0.0   0.5 0.5 0.5   0.0 0.5 0.5',
           'sc': 'ND = 7; NPOINTS = 100\n'+
                 'QI = 0.0 0.0 0.0   0.0 0.5 0.0   0.5 0.5 0.0   0.0 0.0 0.0   0.5 0.5 0.5   0.0 0.5 0.0   0.5 0.5 0.0\n'+
                 'QF = 0.0 0.5 0.0   0.5 0.5 0.0   0.0 0.0 0.0   0.5 0.5 0.5   0.0 0.5 0.0   0.5 0.5 0.0   0.5 0.5 0.5',
          'fcc': 'ND = 3; NPOINTS = 100\n'+
                 'QI = 0.0 0.0 0.0   0.5 0.5 1.0   0.0 0.0 0.0   0.5 0.5 0.5   0.00 0.50 0.5   0.25 0.75 0.5\n'+
                 'QF = 0.0 0.5 0.5   0.0 0.0 0.0   0.5 0.5 0.5   0.0 0.5 0.5   0.25 0.75 0.5   0.50 0.50 0.5',
          'bcc': 'ND = 1; NPOINTS = 100\n'+
                 'QI = 0.0 0.0 0.0\n'+
                 'QF = 0.5 0.5 0.5'
        }

##### set up INPHON file for generating phonon bandstructure
INPHON=['none']*11
INPHON[0]='LRECIP = .T.'
INPHON[1]=BZpaths[latsys]
INPHON[2]='MASS = '+" ".join(str(atomic_mass[str(entry)]) for entry in type)
INPHON[3]='LSYMM = .T.'
INPHON[4]='LFORCEOUT = .T.'
INPHON[5]='LSUPER = .F.'
INPHON[6]='NDIM = 1 1 1'
INPHON[7]='NTYPES = '+str(len(type))
INPHON[8]='LCENTRAL = .T.'
INPHON[9]='IPRINT = 3'
INPHON[10]='LEIGEN = .T.'
##### Write INPHON file in start-subdirectory
wf=open('./phoneval/INPHON','w')
for line in range(len(INPHON)):
  wf.write(INPHON[line]+'\n')
wf.close()

##### calculate phonon frequencies
os.chdir('./phoneval')
p=subprocess.Popen(PHON+' > output',shell=True)
p.wait()
##### draw and save phonon band structure
file=open('./FREQ.cm')
f=file.read()
data=f.splitlines()
columns=[]
for line in range(len(data)):
  columns.append(data[line].split())
for column in range(1,len(columns[0])):
  x=[]
  y=[]
  for line in range(len(data)):
    x.append(float(columns[line][0]))
    y.append(float(columns[line][column]))
  mpl.plot(x,y)
mpl.axis([float(columns[0][0]),float(columns[-1][0]),min([float(columns[i][1]) for i in range(len(data))])-10,max([float(columns[i][-1]) for i in range(len(data))])+50])
mpl.ylabel('Phonon frequency in 1/cm') #TODO: Add proper xlabels for the high-symmetry points.
mpl.savefig('./bandstructure.eps')
os.chdir('../')
###TODO: Add plots for thermodynamical properties. They are readily available with PHON (DOS).

print(" done. You should now have a bandstructure.")