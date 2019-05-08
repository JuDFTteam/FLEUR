# One-dimensionalStructures

## One-dimensional structures

In case of one-dimensional (1D) calculations, the unit cell geometry is different from that for the bulk and film cases. For the details of the geometrical setup and representation of the charge, potential and basis functions in 1D-case, please refer to the following paper: 

*   [Y. Mokrousov, G. Bihlmayer and S. Blügel, Phys. Rev. B, 72, 045402 (2005)][2] 

Briefly, the space is divided into the muffin-tin, interstitial and the vacuum regions. The considered structure is periodic along the z-direction with the period T. The system is embedded in to a cylinder of diameter Dvac with the z-axis as a symmtery axis. The region outside this cylinder is the vacuum region. The region in between the cylinder and the muffin-tin spheres around the atoms is the interstitial region (Fig.1 in the [reference paper][3] and its caption). 

In the muffin-tin region the representation of the charge, potential and basis functions goes in conventional spherical terms. In the interstitial region the quantities are described in terms of plane-waves. To generate three-dimensional reciprocal vectors, an auxiliary in-plane square lattice with the side tilde-D (> Dvac) is used. 

In the vacuum region, cylindrical representation of the quantities is employed, with the plane-wave expansion along the z-axis, decaying on infinity radial expansion and in terms of the angular harmonics around the z-axis. For the charge and potential the angular expansion goes up to an integer MM, while for the basis functions at least a twice smaller parameter vM is used (both are reffered to as mmax in the [reference paper][4]). In order to speedup the process of the Hamiltonian matrix construction an additional angular parameter m_cyl, which is even smaller than vM, is used. 

An extensive use of the rotational, as well as chiral symmetries is possible within the 1D version of FLEUR. Every given chiral symmetry can be uniquely specified by two integers: chi and rot (in the [reference paper][5], M and N), in a way, explained in the 

[reference paper][6]. 

For 1D calculations, certain modifications are required in the inp-file. Below we discuss these modification taking as an example a monoatomic Al chain with the separation between the atoms of 4.573 a.u. and one atom in the unit cell. 



### Al monoatomic chain example input file 

    1|strho=f,film=t,dos=f,isec1=99,ndir= 0,secvar=F
      2|Al MW
      3|squ p1  ,invs=F,zrfs=F,invs2=F,jspins=1,l_noco=F,l_J=F
      4|   7.00000
      5|  5.000000  4.573000  1.000000
      6|rpbe   non-tivistic
      7|igrd=1,lwb=F,ndvgrd=6,idsprs=0,chng= -.100D-11
      8|&odim d1=t,MM=24,vM=12,m_cyl= 8,chi= 1,rot= 4,invs1=t,zrfs1=F /
      9| 1
     10|**********************************
     11|Al  13    4    8  451  2.260000   .023000
     12|
     13| 1,force =f,nlo= 0,llo= 0
     14|  0.000000   .000000  0.000000  1.000000
     15|**********************************
     16| 12.000000 11.000000
     17|vchk=f,cdinf=F,pot8=f
     18|lpr=0,form66=F,l_f=f,eonly=F,eig66=f
     19|  6  6
     20|  1  0
     21|Window # 1
     22|  -0.60000   0.20000   3.00000
     23|   4.00000
     24|gauss=f    .00150tria=F
     25|    .00000    .00000,l_soc=f
     26|frcor=F,slice=F,ctail=F,disp=F,kcrel=0,u2f=F,f2u=F
     27|itmax= 6,maxiter= 19,imix= 7,alpha=  .040,spinf=  1.00
     28|swsp=f  1.50
     29|lflip=F  1  1
     30|vacdos=F,layers= 1,integ=F,star=F,nstars= 0      .00      .00      .00      .00,nstm=0,tworkf=   .000000
     31|
     32|iplot=F,score=F,plpot=F
     33|  0  -.400000   .200000,nnne=  0,pallst=F
     34|xa=   2.00000,thetad= 300.00000,epsdisp=    .00010,epsforce=    .00010
     35|relax 111 111
     36|emin_dos=  -0.60000,emax_dos=   0.10000,sig_dos=    .00450
     37|nkpt= 18
    



*   Line 1: film option is always set to true, film=t. 
*   Line 3: latnam=squ, spgrp=p1, invs=f, zrfs=f and invs2=f, ALWAYS. 
*   Line 4: parameter tilde-D is given. 
*   Line 5: parameters Dvac and T are given as first two numbers, respectively. 
*   Line 8 is the "1D"-line, listing most of the parameters, needed for 1D calculations, which are described in the beginning. For this particular case only c4-rotational symmetry is used (chi=1,rot=4). Parameters invs and zrfs are for inversion and z-reflection symmetry in the system. Be aware, that chi>1 is incompatible with the inversion symmetry. 
*   The coordinates of the atoms are ALWAYS given in rectangular (x,y,z) cartesian coordinates. The atoms are grouped into atom types according to the symmetry of the system. 
*   Line 24: tria is always switched off, tria=f. 
*   The k-points are generated automatically in a half of the 1D Brillouin zone (along the z-axis), the number of k-points is specified by nkpts in line 37. 

The rest of the parameters in the inp-file is of the same meaning as for the bulk and film setups. 



### Important Remarks 

*   To create a band structure, the switch "dos" in the inp-file should be set to true and ndir to 1. By running fleur.x a dosinp file is generated. After that a small executable bnstr.x, obtained by compiling this small fortran program bnstr.F, is to be run in the same directory, which produces file(s) bands1 (and bands2 for jspins=2), suitable for further processing with xmgrace or any other program of the kind. When running bnstr.x the fermi energy (which you can 'grep fermi out'), number of spins, and the name of the file (dosinp) will be asked from you. 
*   For generating a 1D inp-file, a special 1D-input-file-generator is required 
*   The value of the parameter MM is in principle defined from the condition that 2\*MM=kmax\*Dvac. In reality, however, much smaller values of MM are enough to achieve consistency in the values of calculated quantities. 
*   The difference between Dvac and tilde-D values should lie in the range of 2-4 a.u. For larger systems with large diameter Dvac the value of this difference can be reduced, nevertheless, it is always desirable to make sure that your results are consistent with respect to the chosen difference of Dvac and tilde-D. 
*   For 1D calculations the sym.out file is not needed!

 [2]: http://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=PRBMDO000072000004045402000001&idtype=cvips&gifs=yes
 [3]: http://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=PRBMDO000072000004045402000001&idtype=cvips&gifs=yes
 [4]: http://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=PRBMDO000072000004045402000001&idtype=cvips&gifs=yes
 [5]: http://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=PRBMDO000072000004045402000001&idtype=cvips&gifs=yes
 [6]: http://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=PRBMDO000072000004045402000001&idtype=cvips&gifs=yes
