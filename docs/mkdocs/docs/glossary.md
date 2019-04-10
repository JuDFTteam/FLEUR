Glossary
===============
Here we will describe a few terms often used in the context of FLEUR calculations

# atomic units

Almost all input and output in the FLEUR code is given in atomic units, with the
exception of the U and J parameters for the LDA+U method in the input-file and the
 bandstructure and the DOS output-files where the energy unit is eV.

energy units: 1 Hartree (htr) = 2 Rydberg (Ry) = 27.21 electron volt (eV)  
length units: 1 bohr (a.u.) = 0.529177 &Aring;ngstr&ouml;m = 0.0529177 nm   
electron mass, charge and Planks constant h / 2 &pi; (&#x210F;) are unity  
speed of light = e'^2^'/&#x210F; 1/ &alpha; ;
fine-structure constant &alpha;: 1/&alpha; = 137.036

# band gap

The band-gap printed in the output ([[out]] file) of the FLEUR code is the energy separation
between the highest occupied Kohn-Sham eigenvalue and the lowest unoccupied one.
Generally this value differs from the physical band-gap, or the optical band-gap,
due to the fact that Kohn-Sham eigenvalues are in a strict sense Lagrange multipliers
and not quasiparticle energies (see e.g. Perdew & Levy, [PRL 51, 1884 (1983)](http://dx.doi.org/10.1103/PhysRevLett.51.1884)).

# core levels

States, which are localized near the nucleus and show no or negligible dispersion 
can be treated in an atomic-like fashion. These core levels are excluded from the
valence electrons and not described by the FLAPW basisfunctions. 
Nevertheless, their charge is determined at every iteration by solving a Dirac
equation for the actual potential. Either a radially symmetric Dirac equation is solved
(one for spin-up, one for spin-down) or, if @@kcrel=1@@ in the input file,  even a 
magnetic version (cylindrical symmetry) is solved. 

# distance (charge density)

In an iteration of the self consistency cycle, from a
given input charge density, &rho;'^in^', a output density, &rho;'^out^', is calculated.
As a measure, how different these two densities are, the distance of charge densities
(short: distance, d) is calculated. It is defined as the integral over the unit cell: 
{$ d = \int || \rho^{in} - \rho^{out} || d \vec r $}\\
and gives an estimate, whether self-consistency is approached or not. Typically,
values of 0.001 milli-electron per unit volume (a.u.'^3^') are small enough to 
ensure that most properties have converged.
You can find this value in the out-file, e.g. by @@grep dist out@@.
In spin-polarized calculations, distances for the charge- and spin-density are
provided, for non-Collinear magnetism calculations even three 
components exists. Likewise, in an  LDA+U calculation a distance of the
density matrices is given. 

# energy parameters

To construct the FLAPW basisfunctions such, that only the relevant (valence) electrons
are included (and not, e.g. 1s, 2s, 2p for a 3d-metal) we need to specify the energy
range of interest. Depending slightly on the shape of the potential and the muffin-tin radius,
each energy corresponds to a certain principal quantum number "n" for a given "l". E.g.
if for a 3d transition metal all energy parameters are set to the  Fermi-level, the basis
functions should describe the valence electrons 4s, 4p, and 3d. 
Also for the vacuum region we define  energy parameters,
if more than one principal quantum number per "l" is needed, local orbitals can be
specified.

# Fermi level

In a calculation, this is the energy of the highest occupied eigenvalue (or, 
sometimes it can also be the lowest unoccupied eigenvalue, depending on the
"thermal broadening", i.e. numerical issues). In a bulk calculation, this energy
is given relative to the average value of the interstitial potential; in a
film or wire calculation, it is relative to the vacuum zero.
 
# interstitial region

Every part of the unit cell  that does not belong to the  
 muffin-tin spheres and not to the vacuum region. Here, the basis (charge density,
potential) is described as 3D planewaves.

# lattice harmonics

Symmetrized spherical harmonics. According to the point group of the atom, only 
certain linear combinations of spherical harmonics are possible. A list of these
combinations  can be found at the initial section of the out-file.

# local orbitals

To describe states outside the valence energy window,  it is
recommended to use local orbitals. This can be useful for 
lower-lying semicore-states, as well as unoccupied states (note, however, that this just 
enlarges the basis-set and does not cure DFT problems with unoccupied states).

# magnetic moment

The magnetic (spin) moment can be defined as difference between "spin-up" and "spin-down" charge,
either in the entire unit cell, or in the  muffin-tin spheres.
Both quantities can be found in the out-file, the latter one explicitly marked by " --> mm",
the former has to be calculated from the charge analysis (at the end of this file). \\
The orbital moments are found next to the spin-moments, when SOC is included in the calculation.
They are only well defined in the muffin-tin spheres as
{$ m_{orb} = \mu_B \sum_i < \phi_i | r \times v | \phi_i > $}.\\
The in a collinear calculation, the spin-direction without SOC is arbitrary, but assumed to be 
in z-direction. With SOC, it is in the direction of the specified spin-quantization axis. The
orbital moment is projected on this axis. In a non-collinear calculation, the spin-directions
are given explicitely in the input-file.

# muffin-tin sphere

Spherical region around an atom. The muffin-tin radius is an important input parameter.
The basis inside the muffin-tin sphere is described in spherical harmonics times a
radial function. This radial function is given numerically on a logarithmic grid. The
charge density and potential here are also described by a radial function times a
the lattice harmonics.
