# ElectricFields

# Implementation of Electric Fields in FLEUR 

Here, we describe how to put a metallic system in an electric field. For the actual application you might want to jump to see the [input file settings][1]. In case of nonmetallic systems, you might be interested to find the description [how to place an asymmetric field][2] (available from version 0.25g on) or in an electric field with [Dirichlet boundary conditions][3] (available in version 0.26c). 



## The screening charge

The ability to include applied electric fields in our electronic structure calculations is a very powerful tool. We can then begin to study magnetic tunnel junctions as these systems are a junction under the influence of an applied E-field or we can study field induced reconstructions at surfaces. It opens the door for us to a whole new class of systems and effects. 

OK so how do we do it then? 

In the FLEUR code we have a film, so the simplest method of introducing a field is to place a sheet of charge on either side of the film in the vacuum, this then provides an electric field perpendicular to the surface of the film. We need to ensure charge neutrality of the system so we must subtract the same amount of charge from the slab as we have placed on the charged sheets; that's where our screening charge comes from! The figure below shows schematically our film, the charge sheet and the induced screening charge. To answer the question that everyone asks, yes we could just put equal and opposite charges on either side of the film and then not have to subtract any charge from the film, we just have not implemented it yet! [Note: That has changed meanwhile, see below.] 



![][4]

Can we test that it works OK? 



![][5]

Yes, we can test it. We can show that there is a relationship between the change in electronic occupation and the the chemical potential: Remembering that applying an electric field corresponds to changing the occupancy we now have a test. We have two ways of calculating the change in total energy as a function of the change in occupancy (which corresponds to the change in field). We can calculate it directly from the total energy of two calculations (one with an applied field, one without) and we can calculate it by integrating the curve that represents the chemical potential as a function of applied field. The results of the test are shown below. 



![][6]

We see from the above figures that the open squares that represent the results obtained from integrating the chemical potential curve sit almost perfectly centered upon the results calculated directly from the total energy calculations, so we can be confident that the method works well. 



### Example results 

The first results show the screening charge induced in the Ag(001)c(2x2)-Xe structure when a field is applied. We can see that the screening still occurs at the metal surface and not where the adsorbate sits and that the adsorbate atom becomes polarised. 



![][7]

The second example shows the spin-resolved screening charge in Fe(001). Here we see that even though the screening charge in the separate spin channels penetrates right to the center of the slab, the total screening charge only penetrates as deep as the first layer - the field is screened out. 



![][8]



* * *

 



## Input file settings 

The values of the sheets of charge for external electric fields is set by requiring charge neutrality. Thus, you can add or subtracting a fractional charge by e.g. changing the number of electrons in the energy window. The resulting field is shown in the `external electric field` section of the `out` file. 

The position of the sheets of charge relative to the vacuum boundary is set by default to 10 a.u. (= 5.291772 Å), but can be set to a different value in the file [`apwefl`][9]. In the second line of that file, you can specify an additional sigma (`sig_b`) for the upper and one for the lower plate; if the line is not present, 0.0 is assumed for both. 



## Schematic representation of a two-dimensional film in Fleur

    ---------------------------------           -.
                                                  \
    ================================= -.           \ delz*nmz  = 25 a.u.
                                        \          /
    .................................    > zsigma /       -.
                                        /        /          \
    ********************************* \          \           \
    *********************************  > z1 = D/2 \           \
    ********************************* /            > D = Dvac  > Dtilde
    *********************************             /           /
    *********************************            /           /
                                                            /
    ................................                      -'
    
    
    ================================
    
    --------------------------------
    



*   `***` denotes the film (= interstitial + muffin tins); the width Dvac (or short: D) is given in "inp" 
*   `...` denotes a slightly larger region of space to allow for more variational freedom and prevents nodes in the wavefunction at z = z1 (== Dvac/2), where Psiinterstitial and Psivacuum are matched. Dtilde is also specified in "inp" 
*   `===` denotes in electric-field calculations the position of charged plates. Default width of zsigma (= zsheet - z1) is 10 a.u. (= 5.291772 Å) but it can be adjusted via the "apwefl" file 
*   `---` denotes the end of the vacuum; an exponential decay to zero is assumed in the integrals when going from z = z1 + delz*nmz to infinity. The value delz*nmz is hard-wired to 25 a.u. (= 13.229430 Å) in "inped.F"; the number of steps 'nmz' can be set in "fl7para" (default: 250, i.e. 0.1 a.u. steps). Up to "nmzxy" g|| =/= are taken into account; this zone extends to delz*nmzxy, which is by default 10 a.u. (nmzxy can be set in "fl7para"; default: nmzxy = 100). 

 



## Placing asymmetric fields

Sometimes, it might be useful to place a thin film in an antisymmetrc or even asymmetric field. Since this requires to place charges of opposite sign on both sides of the film, it is necessary to provide an additional file, the above-mentioned [`apwefl`][9] file.It allows to place the film into a condenser that creates the electric field like that: 



![][10]

In the image above, we omitted the film for clarity and just show the effect of the plates put at +/-10 a.u. in the vacuum. For polar films this is the correct setup, since it allows to have a flat potential in the vacuum and a vanishing electric field inside the film. Notice, however, that for large applied fields the potential in the vacuum can drop below the Fermi-level and you can get field-emission (i.e. the program stops with a message like `vacuz`). 



## Placing fancy inhomogeneous fields

Since version 0.26b inhomogeneous fields can be generated using the [`apwefl`][9]. The syntax is as follows: Empty lines and lines starting with a hash (`#`) or exclamation mark (`!`) are ignored. In order to use the new syntax, the first line needs to be either empty, start with a comment character (\`!\` or \`#\`) or with the (otherwise optionally) namelist `&efield`. In the next line the variables and default values of the namelist `&efield` are shown. 



    &efield zsigma = 10, sig_b(:) = 0.0, plot_charge = .false., plot_rho = .false., autocomp = .true. /
    

Hereby, `zsigma` is the distance in a.u. of the charge sheets from interstitial-vacuum boundary, `sig_b(1)` and `sig_b(2)` contain the additional (homogeneous) charge for the top and bottom sheet. If `plot_charge` is true and the in-plane charge is inhomogeneous, the files `efield-1.dat` (and, if `nvac=2`, `efield-2.dat`) contain the two-dimensional charge distribution for plotting; if `plot_rho` is true, the result of the real-space-to-G-space-to-real-space transformation is saved in `efield-fft-1.dat` (and `efield-fft-2.dat`). By default, excess (positive/negative) charge of the film is compensated by equally charging the charge sheets; if `autocomp` is true, the user has to do this manually. Note: Fleur requires an overall (film plus top plus bottom sheet) charge neutrality. 

The inhomogeneous charge can be places using the tags `rect`, `circ`, `rectSinX`, `polygon`, and `datafile`. Their syntax is: 



    rect <sheet> <x>,<y> <width>,<height> <charge> [options]
      circ <sheet> <x>,<y> <radius> <charge> [options]
      rectSinX <sheet> <x>,<y> <width>,<height>  <amplitude> <n> <delta> [options]
      polygon <sheet> <n_points> <x1>,<x2> ... <x_n>,<y_n> <charge> [options]
      datafile <filename> [zero_avg] [options]
    

Note that all positions and lengths are currently relative values (i.e. between 0 and 1). The sheet to be used can be set using `<sheet>`, which can be either `top`, `bot` or `topbot`/`bottop`. Options are: `add` (default) to add the charge to the charge of previous tags or `replace` to use the new charge instead; `zero` to place the charge only to areas which were before zero, `nonzero` to place it at areas which where before nonzero or `all` (default) to place it in the whole area covered by the tag. 

Note: The regions can exceed the unit cell plane and then cut off, e.g. 

    circ top 0,0  0.25  0.5
    

places half an electron in a quarter circle with origin (0,0). Note, however, that `circ` creates a perfect circle only on the grid; this generates a circle and not an ellipse only if the *k1d*/*k2d* ratio matches the crystal's *a*/*b* ratio. 

`rectSinX` creates a sinodal potential in *x* direction (constant in *y* direction for any *x* value), i.e. V(x,y) = A \sin(2\pi nx + 2\pi\delta), where *A* is the amplitude; however, the argument in `apwefl` is not *A* directly but A' = A L_z, where L_z is the number of points in *z* direction. Contrary to `circ` and `rect`, charges are mask out without being redistributed to non-masked positions. It is \int_0^{2\pi} A|\sin(x)| {\rm d}x = 4A', `n` is the order and `delta` (\delta) the offset. 

`polygon` creates a polygon-shaped charge distribution; note that the currently used algorithm does not always give the perfectly shaped polygon - and the edge points are not always included in the polygon. 

The `datafile` reads the data from a file; if a `zero_avg` has been given, the charge is averaged to zero, i.e. only the inhomogeneous contributions are taken into account. The option `replace`/`add` is supported, but `zero`/`not_zero` is not. The syntax of the data file itself is as follows. First line: number of *x* and *y* points; second line: charge for point (x=1,y=2), third: (x=1,y=2) etc. The number of points must be `3*k1d` and `3*k2d`. 

**Example 1:** To have two top plates (segments): 

    rect top 0,  0  0.5,1.0  0.2
      rect top 0.5,0  0.5,1.0 -0.2
    

**Example 2:** To have a charged ring with 0.2e and -0.2e of charge evenly distributed outside this ring: 

    circ topBot 0.5,0.5 0.2  1           ! Create temporary an inner ring
      circ topBot 0.5,0.5 0.3  0.2 zero    ! Create outer ring
      circ topBot 0.5,0.5 0.2  0   replace ! Delete inner ring
      rect topBot 0,0     1,1 -0.2 zero    ! Place smeared opposite charge
    

 



## Using Dirichlet boundary conditions

In the scheme above, a (locally) constant surface charge density has been assumed, which matches a Neumann boundary condition. Another common boundary condition is a constant potential (Dirichlet boundary condition) -- such as provides by a metal plate, which can be either grounded or kept at a certain potential. In order to use Dirichlet boundary condition, add `Dirichlet=T` flag to the namelist `&efield`, which has been introduced above. The value given for `sig_b` and in the different shapes are regarded as potential in Hartree or, if `eV=T` is set, as potential in (electron) volts.

 [1]: #file_settings
 [2]: #asymmetricfield
 [3]: #Dirichlet
 [4]: http:/images/newgeom.gif ""
 [5]: http:/images/equ.png ""
 [6]: http:/images/muvse.gif ""
 [7]: http:/images/webfig.gif ""
 [8]: http:/images/webfe.gif ""
 [9]: https://www.flapw.de/pm/index.php?n=User-Documentation.SmallUseful#apwefl
 [10]: http:/images/condenser.png ""
