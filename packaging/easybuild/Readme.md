To build and use FLEUR with EasyBuild:

```
module load EasyBuild

JSC_OVERRIDE_TOOLCHAIN_CHECK=1 eb --detect-loaded-modules=ignore --use-existing-modules --installpath=$PWD/easybuild --buildpath=$PWD/build --tmpdir=$PWD/  Fleur.eb # change paths to your desired location

module use $PWD/easybuild/modules/

#ml chem/FLEUR/MaX-R8.0-GCCcore-.13.3.0
#ml chem/FLEUR/MaX-R8.1-GCCcore-.13.3.0
ml chem/FLEUR/develop-GCCcore-.13.3.0
```

## FOSS toolchain

`FLEUR-8.1-foss-2025b.eb` builds the MaX-R8.1 release with the community `foss/2025b` toolchain
(GCC, OpenMPI, FlexiBLAS, FFTW, ScaLAPACK) in the standard EasyBuild way. It uses the
MPI compiler wrappers, calls CMake directly and installs `fleur_MPI` and `inpgen` into `bin/`.
All dependencies (HDF5, libxml2, ELPA) are EasyBuild modules from the official easyconfigs:

```
eb FLEUR-8.1-foss-2025b.eb --robot
module load FLEUR/8.1-foss-2025b
```

The build fails before compiling if MPI, HDF5, ScaLAPACK or ELPA were not detected,
so a missing library cannot silently produce a reduced FLEUR. libxc and Wannier90 are not included.
The MKL FFT is switched off explicitly, FFTW from the toolchain is used.
