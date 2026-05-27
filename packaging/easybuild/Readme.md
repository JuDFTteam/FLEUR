To build and use FLEUR with EasyBuild:

```
module load EasyBuild

JSC_OVERRIDE_TOOLCHAIN_CHECK=1 eb --detect-loaded-modules=ignore --use-existing-modules --installpath=$PWD/easybuild --buildpath=$PWD/build --tmpdir=$PWD/  Fleur.eb # change paths to your desired location

module use $PWD/easybuild/modules/

#ml chem/FLEUR/MaX-R8.0-GCCcore-.13.3.0
ml chem/FLEUR/develop-GCCcore-.13.3.0
```
