# PlottingAndSlicing

## Plotting 

If you just want to plot the charge density in a plane or make an isosurface in a 3D-plot, it is sufficient to provide cdn1, inp and sym.out and set `iplot=T`. Then the program will stop, complaining that it wants a plot_inp file and generate a template for you. This has to be modified (see the description of the plot_inp file) and then run `fleur.x`. You will get a two or three-dimensional array of density values, which can be visualized with your favorite plotting program (output for [XCrySDen][1] is provided, a very flexible software is e.g. [OpenDX][2], formerly data-explorer). 



### spin-polarized case

In case of a spin-polarized calculation, the output (e.g. `plot.xsf`) will contain two data-sets for the spin-up and the spin-down part of the density. With most plotting programs (e.g. XCrySDen) you can then sum up the two parts or subtract them to get the full charge density of the spin-density, respectively. 



### potentials

Notice, that you can in this way also plot potentials: Set `plpot=T` and perform a single iteration. Then, you get two files, potcoul\_pl and pottot\_pl, which can be processed further. Copy them to cdn1, then proceed as described above (with `iplot=T` and `plpot=T`). 



### valence density

Normally, the total charge density is not so interesting. If you want to plot it, use a logarithmic scale. Nicer to plot is the valence charge density, since the peak at the nucleus is not so pronounced. For this plot, you also need a cdnc-file, containing the core charge to be subtracted. Set `score=T` in the inp-file and proceed as described above. 



## Slicing 

For plotting charge densities of states in selected energy regions one first calculates the density by choosing `slice=T ` (and iplot=F). You have several possiblilites to make this so-called slices: 

In the inp-file you can find a line: 

    0  0.000000  0.000000,nnne=  0,pallst=F
    



*   Suppose, you want to make a charge density of all occupied states within 1eV from the fermi level (at -0.15 htr). Then this line would read: 0 -0.186751 -0.150000,nnne= 0,pallst=F 

i.e. the second and third number define the e\_low and e\_up in this line of the input file. 



*   If states from the unoccupied part of the spectrum are reqired, set `pallst=T` 0 -0.150000 -0.113249,nnne= 0,pallst=T 

which plots all states (pallst) up to 1eV above the Fermi level (of -0.15 htr). 



*   Alternatively you can also select special k-points (with the first number) and/or eigenvalue (nnne). If these numbers are zero, no further restriction is made. If you want to plot the 16th state of the first k-point, set `nnne= 16` and specify a "1" at the beginning of the line. Use pallst to occupy eigenvalues above E_fermi. Note: For a selected k-point (non-zero first number) and "nnne/=0" only one state is considered and the energy range is ignored. In earlier versions of the code, "nnne" is ignored unless both boundaries of the energy range are zero. 

The result is in any case a file called cdn_slice which can then be taken to make plots. To do this one chooses slice=T and iplot=T and sets up the plotin-file which contains the information of the geometry of the desired planes. Newer versions of the code support the plot_inp file described above.

 [1]: http://www.xcrysden.org/
 [2]: http://www.opendx.org/
