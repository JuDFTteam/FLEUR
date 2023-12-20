# Changes and Status of different FLEUR versions


Here we try to document the major changes, their status and issues for the different FLEUR releases. Please note that this list is not necessarily complete nor always up to date. Please consult the documentation, commit messages and/or the issues.

Please remember that a lot of bugfixes are not documented here and thus we always recomment to use the newest versions of FLEUR.

## Changes for MaXR7.0

### Features:

- Initial release of DFPT phonon code. (still experimental feature)
- Initial release of LDA+U+V code. (still experimental feature)
- Reactivation of HSE hybrid functional.
- Hybrid functionals in combination with 2nd variation spin-orbit coupling.
- Alternative parameter setup profiles for the input generator.

### Notable bugfixes:

- Fixes for spin-spiral calculations with local orbitals.
- Fixes related to deadlocks in parallelized hybrid functional calculations.
- Fixes related to calculations on noncollinear magnetism.

### Configure/build/usage
- include in-build interface to ELSI
- first version of "fleurist": a python tool to provide parallelization strategies
- possibility to have not fully load-balanced parallelization to use the k-point parallelism more efficiently

