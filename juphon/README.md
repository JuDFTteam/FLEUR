# The FLEUR phonon cookbook

This short document serves as a guide to use FLEUR to generate phonon spectra from scratch. It highlights the necessary preliminaries, how to generate FD spectra in conjunction with the phonopy code and how to use FLEUR's DFPT capability to generate phonon data. A short section on current limitations/future developments is also given.

It will require a workstation, where FLEUR and phonopy are both installed and working. Either have that or run the phonopy specific steps (they need barely any computational effort) on a home PC with modern python.

# Ground-state calculation

First, set up an input file for inpgen. In this case we make one for fcc copper and call it ```inpCu```. It looks as follows:

```
Cu fcc

0.0 1.0 1.0 ! a1
1.0 0.0 1.0 ! a2
1.0 1.0 0.0 ! a3
6.650896828 ! aa
0.5 0.5 0.5 ! scale

1 ! num atoms
29.1 0.0 0.0 0.0

&atom element="Cu" id=29.1 lmax=9 lnonsph=7 jri=981 /
&comp kmax=4.5 gmax=15.0 gmaxxc=15.0 /
&exco xctyp='vwn' /
&kpt div1=16 div2=16 div3=16 tkb=0.0005 /
```

This is a bit more complex than the basic FLEUR input, because it has some comments as flags and a structure, that makes it readable by phonopy. We therefore don't need to reformulate it for the second step.

Now, run inpgen:
```
/path/to/inpgen -f inpCu
```

Ensure several things in the resulting inp.xml:
- set ```ctail="F"```
- reduce ```radius``` by a slight margin
- set the iteration count high enough.
- set ```name``` in ```xcFunctional``` to that of an LDA functional available in libxc; optimally, you directly set (e.g. for vwn):
```
      <xcFunctional name="LibXC" relativisticCorrections="F">
         <LibXCName  exchange="lda_x" correlation="lda_c_vwn"/>
      </xcFunctional>
```

Now you need to optimize your input. Converge it too a minimal energy w.r.t. the lattice constant and cutoff parameters. This is why the MT radius was set lower, so the spheres do not crash into each other when reducing the lattice size.

Once you are satisfied with your optimization, write a new input file (e.g. inpCu_fit) and fill it with the parameters you converged to. This will be your base input for the follow up calculations.

# Finite displacement calculations



# Density-functional perturbation theory

# Current limitations

# Future development