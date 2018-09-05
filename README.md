Welcome to the source code of FLEUR
===================================

Please note that the documentation of the
code can be found at the [FLEUR Homepage]
(http://www.flapw.de/pm/index.php?n=User-Documentation.Documentation).

For further instructions on Installation/Usage,
please check the [FLEUR Homepage]
(http://www.flapw.de/pm/index.php?n=User-Documentation.Documentation).

The rest of this document summarizes only the 
structure of the FLEUR source code. Did
we mention to check the Homepage?

The source of FLEUR is organized in several 
subdirectories. Some of them collect code 
specific for particular features, others code
relevant for crutial steps in the code or simply
code that is usually executed together.
Here a short description of the directories:

* main: contains the main-program and several core subroutines
* init: stuff for the initialization (called from fleur_init in main)
* vgen: potential generation (called from vgen in main)
* eigen: setup of the eigenproblem
* diagonalization: various methods to diagonalize the Hamiltonian
* cdn: general code for the generation of charge
* cdn_mt: charge generation in MT-spheres
* force: code related to the evaluation of forces
* mix: charge/potential mixing routines
* ldau: routines needed in case of LDA+U calculations
* inpgen: code for the input generator (seperate executable inpgen)
* fermi: determination of the fermi-level
* eigen_secvar: second variational solution of the Hamiltonian
* eigen_soc: Spin-orbit related code
* core: Core states
* dos: Code for Density of states, bandstructures
* orbdep: Code for quantities depending on orbitals
* optional: code that is used in special cases like inital charge generation
* wannier: wannier related code
* xc-pot: various exchange-correlation potential routines
* mpi: code for parallel execution
* io: subroutines doing IO
* juDFT: timing, error handling, etc
* eels: code for electron-energy loss spectroscopy
* hybrid: code for hybrid functionals 
* math: code providing math functionality
* include: c-type include files
* global: code used everywhere (here you find types.F90 with the data-types)
* cmake: definitions used by cmake

If you modify FLEUR please do so in the develop branch by running
'git checkout -t origin/develop'
after cloning the git repository. For larger changes you might want to
create your own branch.

## Compiling Fleur on Jureca
To compile Fleur on Jureca the following modules need to be loaded:

```bash
module load intel-para CMake HDF5 libxml2/.2.9.7 ELPA/2017.11.001-hybrid
```
which is then configured with
```bash
./configure.sh JURECA_INTEL
```
or
```bash
./configure.sh -external libxc JURECA_INTEL 
```
if libXC should be supported aswell.

## Compilling Fleur on Jureca Booster

To compile on the booster you need to first switch to the KNL environment

```bash
module purge
ml Architecture/KNL
module load intel-para CMake HDF5 libxml2/.2.9.7 ELPA/2017.11.001-hybrid
```
Then you can compile as discrebed above.

## Running Fleur on Jureca Booster

To run on the Boosters you need to switch the architecture in each job script aswell
```bash
#!/bin/bash -x
#SBATCH --nodes=1
##SBATCH --ntasks=10
##SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=68
#SBATCH --output=mpi-%j.out
#SBATCH --error=mpi-%j.err
#SBATCH --time=1:00:00
#SBATCH --partition=develbooster
#SBATCH --gres=mem96
##SBATCH --mail-user=your.name@fz-juelich.de
##SBATCH --mail-type=END
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
module purge
ml Architecture/KNL
module load intel-para CMake HDF5 libxml2/.2.9.7 ELPA/2017.11.001-hybrid
/work/ias-1/s.rost/fleur_booster/fleur/build/fleur
```

## Developing Fleur

We agreed to use a unified indentation-width of 3. 
Hint: [vim](http://vim.wikia.com/wiki/Converting_tabs_to_spaces) [emacs](https://www.gnu.org/software/emacs/manual/html_node/efaq/Changing-the-length-of-a-Tab.html)