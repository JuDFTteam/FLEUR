# SymmetryOfBands

When doing band structure plots, FLEUR identifies the symmetry of the eigenvectors. Though not implemented for all point-group symmetries, the common cases are present. This shall be illustrated with a band structure of fcc Nickel. 



# Band structure 

fcc Nickel has space group [225][2] or Fm-3m. Its band structure of majority spin along the Delta line (i.e. Gamma-X) is discussed here. Starting from the file dosinp given below, the files bnd1.dat and bnd2.dat were obtained with the command 

    band3xmgr 15 0.314191
    

and bnd1.dat is shown below. The k-points Gamma and X have point group symmetry [0h][3] and [D4h][4] comprising 10 characters each, while the points on the Gamma line have point group symmetry [C4v][5] with five characters. 



![][6]

The Fermi energy is marked with a dashed red line. 



# FLEUR output 

The symmetry data is placed in two files: fort.444 and dosinp. 

When performing a band structure calculation, the file fort.444 is generated (visit the link at the bottom of the page) which contains the symmetry informations of each k-point, but only if it differs from the previous k-point. Therefore we get three tables, for Gamma, for the first point on the Delta line, and for X. Each section begins with the text "Character table for k: ". You can compare this to the tables linked above. These table also contains the degeneracy of eigenvalues. 

In the file dosinp, you have have for each k-point a section listing all eigenvalues followed by its symmetry kind. (About the other information in that line, you can read here: [Calculating the Density of States][6].) For instance for the dosinp file given below starting at line 823 (refering to point number 20), 

    0.377277849596E+00  0.000000000000E+00  0.000000000000E+00    20  0.196078431373E-01
         0.07324 1 0     0.00000
         0.43875     0.08257     0.02693     0.00002
         0.21743 4 0     0.00000
         0.00000     0.00000     0.92980     0.00284
         0.26096 5 0     0.00000
         0.00000     0.00567     0.95576     0.00226
         0.26096 5 0     0.00000
         0.00000     0.00567     0.95576     0.00226
         0.26498 1 0     0.00000
         0.01314     0.00646     0.92652     0.00301
         0.29293 3 0     0.00000
         0.00000     0.00000     0.98334     0.00085
    

you get the symmetry kinds 1, 4, 5, 5, 1, 3 referring to the character listed in tables in file fort.444. 

With this data, you can add some information to the plot. Below the same band structure is shown, with multiplicity and symmetry info for the eigenvalues at the Gamma point and an arbitrarily chosen k-point on the Delta line marked with the vertical red line (to which the cited data a few lines above refer to). The kind of band was read from the dosinp file, name and multiplicity from fort.444. 



![][7]

This provides the proper way to detect band crossings, though it has not been done in the plot above. 



# Footnote 

Example files used: 

*   [![Symbol - externer Link][9]dosinp][9] 
*   [![Symbol - externer Link][10]fort.444][10] 

Comments: 

*   The contents in fort.444 are printed twice. (For each spin, but of course the information is the same. You may ignore this.) 
*   In the given fort.444, the character table is not printed again for Gamma point with spin down (small bug). 
*   In the given fort.444, the band names at the X point are not yet implemented but named 'unknown'.

 [2]: http://www.webelements.com/webelements/elements/text/Ni/xtal.html
 [3]: http://www.cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg=221
 [4]: http://www.cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg=123
 [5]: http://www.cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg=99
 # SymmetryOfBands

When doing band structure plots, FLEUR identifies the symmetry of the eigenvectors. Though not implemented for all point-group symmetries, the common cases are present. This shall be illustrated with a band structure of fcc Nickel. 



# Band structure 

fcc Nickel has space group [![Symbol - externer Link][2]225][2] or Fm-3m. Its band structure of majority spin along the Delta line (i.e. Gamma-X) is discussed here. Starting from the file dosinp given below, the files bnd1.dat and bnd2.dat were obtained with the command 

    band3xmgr 15 0.314191
    

and bnd1.dat is shown below. The k-points Gamma and X have point group symmetry [![Symbol - externer Link][3]0h][3] and [![Symbol - externer Link][4]D4h][4] comprising 10 characters each, while the points on the Gamma line have point group symmetry [![Symbol - externer Link][5]C4v][5] with five characters. 



![][5]

The Fermi energy is marked with a dashed red line. 



# FLEUR output 

The symmetry data is placed in two files: fort.444 and dosinp. 

When performing a band structure calculation, the file fort.444 is generated (visit the link at the bottom of the page) which contains the symmetry informations of each k-point, but only if it differs from the previous k-point. Therefore we get three tables, for Gamma, for the first point on the Delta line, and for X. Each section begins with the text "Character table for k: ". You can compare this to the tables linked above. These table also contains the degeneracy of eigenvalues. 

In the file dosinp, you have have for each k-point a section listing all eigenvalues followed by its symmetry kind. (About the other information in that line, you can read here: [Calculating the Density of States][6].) For instance for the dosinp file given below starting at line 823 (refering to point number 20), 

    0.377277849596E+00  0.000000000000E+00  0.000000000000E+00    20  0.196078431373E-01
         0.07324 1 0     0.00000
         0.43875     0.08257     0.02693     0.00002
         0.21743 4 0     0.00000
         0.00000     0.00000     0.92980     0.00284
         0.26096 5 0     0.00000
         0.00000     0.00567     0.95576     0.00226
         0.26096 5 0     0.00000
         0.00000     0.00567     0.95576     0.00226
         0.26498 1 0     0.00000
         0.01314     0.00646     0.92652     0.00301
         0.29293 3 0     0.00000
         0.00000     0.00000     0.98334     0.00085
    

you get the symmetry kinds 1, 4, 5, 5, 1, 3 referring to the character listed in tables in file fort.444. 

With this data, you can add some information to the plot. Below the same band structure is shown, with multiplicity and symmetry info for the eigenvalues at the Gamma point and an arbitrarily chosen k-point on the Delta line marked with the vertical red line (to which the cited data a few lines above refer to). The kind of band was read from the dosinp file, name and multiplicity from fort.444. 



![][7]

This provides the proper way to detect band crossings, though it has not been done in the plot above. 




Comments: 

*   The contents in fort.444 are printed twice. (For each spin, but of course the information is the same. You may ignore this.) 
*   In the given fort.444, the character table is not printed again for Gamma point with spin down (small bug). 
*   In the given fort.444, the band names at the X point are not yet implemented but named 'unknown'.

 [2]: http://www.webelements.com/webelements/elements/text/Ni/xtal.html
 [3]: http://www.cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg=221
 [4]: http://www.cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg=123
 [5]: http://www.cryst.ehu.es/cgi-bin/rep/programs/sam/point.py?sg=99
 [6]: http://www.flapw.de/pm/datapool/sympsi/fleursym_ni_up_gamma_x_up_1.png ""
 [7]: http://www.flapw.de/pm/datapool/sympsi/fleursym_ni_up_gamma_x_up_2.png ""
