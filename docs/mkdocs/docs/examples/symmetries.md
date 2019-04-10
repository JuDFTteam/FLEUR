#Examples by Lattice

#CrInp
```
bcc Cr

&input film=f /

&lattice latsys=cP a0=1.8897269 a=2.87100 /

2
24.0   0.0 0.0 0.0
24.1   0.5 0.5 0.5

&atom element="Cr" id=24.0 rmt=2.1 jri=981 lmax=12 lnonsph=6 lo="3s 3p" bmu=1.5 /
&atom element="Cr" id=24.1 rmt=2.1 jri=981 lmax=12 lnonsph=6 lo="3s 3p" bmu=-1.5 /
&comp kmax=5.2 gmaxxc=12.5 gmax=15.0 /
&kpt div1=24 div2=24 div3=24 tkb=0.0005 /
```
#AgInp
```
fcc silver

&input film=f /

&lattice latsys=cF a0=1.8897269 a=4.16424 /

1
47 0.0 0.0 0.0

&atom element="Ag" rmt=2.3 jri=981 lmax=12 lnonsph=6 lo="4s 4p" /
&comp kmax=5.0 gmaxxc=12.5 gmax=15.0 /
&kpt div1=25 div2=25 div3=25 tkb=0.0005 /
```
#FeInp
```
bcc Fe

&input film=f /

&lattice latsys=cI a0=1.8897269 a=2.83351 /

1
26     0.0 0.0 0.0

&atom element="Fe" rmt=2.1 jri=981 lmax=12 lnonsph=6 lo="3s 3p" /
&comp kmax=5.2 gmaxxc=12.5 gmax=15.0 /
&kpt div1=27 div2=27 div3=27 tkb=0.0005 /
```
#MnInp
```
tet manganese

&input film=f /

&lattice latsys=tP a0=1.8897269 a=2.5424166 c=3.59552 /

2
25.1   0.0 0.0 0.0
25.2   0.5 0.5 0.5

&atom element="Mn" id=25.1 rmt=2.1 jri=981 lmax=12 lnonsph=6 lo="3s 3p" bmu=2.3 /
&atom element="Mn" id=25.2 rmt=2.1 jri=981 lmax=12 lnonsph=6 lo="3s 3p" bmu=-2.3 /
&comp kmax=5.2 gmaxxc=13.0 gmax=15.5 /
&kpt div1=28 div2=28 div3=20 tkb=0.0005 /
```
#InInp
```
bct In

&input film=f /

&lattice latsys=tI a0=1.8897269 a=3.29841 c=5.06256 /

1
49     0.0 0.0 0.0

&atom element="In" rmt=2.4 jri=981 lmax=12 lnonsph=6  lo="4d" /
&comp kmax=5.0 gmaxxc=12.5 gmax=15.0 /
&kpt div1=17 div2=17 div3=13 tkb=0.0005 /
```
#BrInp
```
Bromium(o)

&input film=f cartesian=T /

&lattice latsys=oC a0=1.8897269 a=8.22860 b=4.22731 c=9.03323 /

 4
   35      0            0.37574      0.88284
   35      0            0.12426      0.38284
   35      0            0.87574      0.61716
   35      0            0.62426      0.11716

&atom element="Br" rmt=2.1 jri=981 lmax=12 lnonsph=6  lo="3d" /
&comp kmax=5.0 gmaxxc=12.5 gmax=15.0 /
&kpt div1= 4 div2= 8 div3= 4 tkb=0.0005 /
```
#HBrInp
```
HBr oF (ICSD # 28842)

&input film=f cartesian=t /

&lattice latsys='oF', a0=1.889727, a=5.5550, b=5.6400, c=6.0630  /

-2
 1 0.0 0.0 0.0
35 0.0 0.0 0.5

&gen 3                          ! in this example the generators
 -1  0  0   0.0                 ! are not necessary, the number of
  0  1  0   0.0                 ! atoms should then be "2"
  0  0  1   0.0

  1  0  0   0.0                 ! as in the int. tables
  0 -1  0   0.0                 ! F m m m(69)
  0  0  1   0.0                 ! http://www.cryst.ehu.es

  1  0  0   0.0
  0  1  0   0.0
  0  0 -1   0.0 /               ! we leave out the centering translation
```
#HgOInp
```
HgO   from ICSD database (http://icsd.fiz-karlsruhe.de/)  #16627

&input film=f cartesian=t /

! 3.311 5.526 3.526 90. 90. 90.
&lattice latsys='oI', a0=1.8897269 a=3.311 b=5.526 c=3.526 /

 -2
80 0.000  0.00    0.00
 8 0.0    0.5     0.17

&gen 2                          ! (generators not really necessary)
 -1  0  0   0.0
  0  1  0   0.0                 ! as in the int. tables
  0  0  1   0.0                 ! for I m m 2(44)

 -1  0  0   0.0
  0 -1  0   0.0
  0  0  1   0.0 /               ! we leave out the translation

&atom element="Hg"  lo="5s 5p" /
&comp kmax=5.0 gmaxxc=12.5 gmax=15.0 /
&exco xctyp="pz" /
```
#Ta2HInp
```
Ta2D  from ICSD database (http://icsd.fiz-karlsruhe.de/)  #61486

&input film=f cartesian=t /

&lattice latsys='oA', a0=1.8897269 a=4.67 b=4.67 c=3.303 /

 -2
73 0.2624 0.25 0.25             ! A-centered setting like in #28022
 1 0.000  0.00 0.00             ! use H for D

&gen 2                          ! as in the int. tables
 -1  0  0   0.0                 ! for C 2 2 2
  0 -1  0   0.0                 ! http://www.cryst.ehu.es
  0  0  1   0.0

 -1  0  0   0.0
  0  1  0   0.0
  0  0 -1   0.0 /               ! we leave out the translation

&atom element="Ta" rmt=2.3 jri=981 lmax=12 lnonsph=6 lo="5s 5p" /
&comp kmax=5.0 gmaxxc=12.5 gmax=15.0 /
&exco xctyp="pz" /
```
#Ta2HBInp
```
Ta2D  from ICSD database #61486 but with oB setting 

&input film=f cartesian=t /

&lattice latsys='oB', a0=1.8897269 a=4.67 b=4.67 c=3.303 /

 -2
73 0.25   0.2624 0.25             ! B centered setting
 1 0.000  0.00   0.00             ! use H for D

&gen 2                          ! as in the int. tables
 -1  0  0   0.0                 ! for C 2 2 2
  0 -1  0   0.0                 ! http://www.cryst.ehu.es
  0  0  1   0.0

 -1  0  0   0.0
  0  1  0   0.0
  0  0 -1   0.0 /               ! we leave out the translation

&atom element="Ta" rmt=2.3 jri=981 lmax=12 lnonsph=6 lo="5s 5p" /
&comp kmax=5.0 gmaxxc=12.5 gmax=15.0 /
&exco xctyp="pz" /
```
#Ta2HCInp
```
Ta2D  from ICSD database #61486 but with oC setting

&input film=f cartesian=t /

&lattice latsys='oC', a0=1.8897269 a=4.67 b=3.303 c=4.67 /

 -2
73 0.25   0.25   0.2624         ! C centered setting
 1 0.000  0.00   0.00           ! use H for D

&gen 2                          ! as in the int. tables
 -1  0  0   0.0                 ! for C 2 2 2
  0 -1  0   0.0                 ! http://www.cryst.ehu.es
  0  0  1   0.0

 -1  0  0   0.0
  0  1  0   0.0
  0  0 -1   0.0 /               ! we leave out the translation

&atom element="Ta" rmt=2.3 jri=981 lmax=12 lnonsph=6 lo="5s 5p" /
&comp kmax=5.0 gmaxxc=12.5 gmax=15.0 /
&exco xctyp="pz" /
```
#PdP2PInp
```
PdP2 (mP) from ICSD database # 48163 (conventional unit cell with all atoms)

! C 1 2/c 1 (15): 6.7787 5.8570 5.8740 90.000 121.769 90.000
! Pd    1       +0.00   4c      0.2500  0.2500  0.0000  0.3220
! P     1       +0.00   8f      0.1886  0.1237  0.3349  0.1620

&input film=f cartesian=f /

&lattice latsys='mP' a0=1.889727, a=5.8740, b=6.7787, c=5.8570, gamma=121.769 /

12
46 0.0000 0.2500  0.2500
46 0.0000 -.2500  -.2500
46 0.5000 0.7500  0.2500
46 0.5000 0.2500  -.2500
15 0.3349 0.1886  0.1237
15 -.3349 -.1886  -.1237
15 0.1651 -.1886  0.1237
15 -.1651 0.1886  -.1237
15 0.3349 0.6886  0.6237
15 -.3349 -.6886  -.6237
15 0.1651 -.6886  0.6237
15 -.1651 0.6886  -.6237
```
#PdP2IInp
```
PdP2 (mI) from ICSD database # 48163 (published data)

! I 1 2/c 1 setting: 6.207 5.857 5.874 90. 111.80 90.
! Pd    1       +0.00   4d      0.25    0.75   0.25
! P     1       +0.00   8f      0.1886  0.1237 0.3537

&input film=f cartesian=t /

&lattice latsys='mI' a0=1.889727, a=5.874, b=6.207, c=5.857, gamma=111.8 /

-2
46 0.2500 0.2500  0.7500
15 0.3537 0.1886  0.1237

&gen 2                          ! as in the int. tables
 -1  0  0   0.5                 ! for   I 1 2/c 1(15)
  0 -1  0   0.0                 ! but now I 1 1 2/a
  0  0  1   0.0

 -1  0  0   0.0
  0 -1  0   0.0
  0  0 -1   0.0 /               ! inversion
```
#PdP2AInp
```
PdP2 (mA)  from ICSD database # 48163 (standardized data)

! C 1 2/c 1 (15): 6.7787 5.8570 5.8740 90.000 121.769 90.000
! Pd    1       +0.00   4c      0.2500  0.2500  0.0000  0.3220
! P     1       +0.00   8f      0.1886  0.1237  0.3349  0.1620

&input film=f cartesian=t /

&lattice latsys='mA' a0=1.889727, a=5.8740, b=6.7787, c=5.8570, gamma=121.769 /

-2
46 0.0000 0.2500  0.2500
15 0.3349 0.1886  0.1237

&gen 2                          ! as in the int. tables
 -1  0  0   0.5                 ! for   A 1 1 2/a (15)
  0 -1  0   0.0                 ! http://www.cryst.ehu.es
  0  0  1   0.0

 -1  0  0   0.0
  0 -1  0   0.5
  0  0 -1   0.5 /               ! inversion
```
#PdP2BInp
```
PdP2 (mB)  from ICSD database # 48163 (standardized data)

! C 1 2/c 1 (15): 6.7787 5.8570 5.8740 90.000 121.769 90.000
! Pd    1       +0.00   4c      0.2500  0.2500  0.0000  0.3220
! P     1       +0.00   8f      0.1886  0.1237  0.3349  0.1620

&input film=f cartesian=t /

&lattice latsys='mB' a0=1.889727, a=6.7787, b=5.8740, c=5.8570, gamma=121.769 /

-2
46  0.2500 0.0000 0.2500
15  0.1886 0.3349 0.1237

&gen 2                          ! as in the int. tables
 -1  0  0   0.0                 ! for   B 1 1 2/b (15)
  0 -1  0   0.5                 ! http://www.cryst.ehu.es
  0  0  1   0.0

 -1  0  0   0.5
  0 -1  0   0.0
  0  0 -1   0.5 /               ! inversion
```
#Binp
```
B (tricl) (APW+lo used for actual calculation)

&input film=f /

&lattice latsys=aP a0=1.8897269 a=4.90067 b=4.90067 c=5.05098 alpha=119.02035 beta=60.97965 gamma=120 /

12
5.0    0.98969      0.01031      0.67465
5.0    0.01031      0.98969      0.32535
5.0    0.98969      0.65404      0.67465
5.0    0.01031      0.34596      0.32535
5.0    0.34596      0.01031      0.67465
5.0    0.65404      0.98969      0.32535
5.0    0.77885      0.22115      0.07274
5.0    0.22115      0.77885      0.92726
5.0    0.77885      0.63043      0.07274
5.0    0.22115      0.36957      0.92726
5.0    0.36957      0.22115      0.07274
5.0    0.63043      0.77885      0.92726

&atom element="B" rmt=1.5 jri=981 lmax=12 lnonsph=6 /
&comp kmax=4.5 gmaxxc=12.0 gmax=13.5 /
&kpt div1=10 div2= 10 div3= 6 tkb=0.0005 /
```
#Cinp
```
grahite C (APW+lo used for actual calculation)

&input film=f /

&lattice latsys=hP a0=1.8897269 a=2.46857 c=8.84079 /

4
 6  0.0  0.0  0.25
 6  0.0  0.0 -0.25
 6  1.0  2.0  0.25
 6  2.0  1.0 -0.25
&factor 3.0 3.0 1.0 /

&atom element="C" rmt=1.3 jri=981 lmax=12 lnonsph=6 /
&comp kmax=5.2 gmaxxc=13.0 gmax=15.5 /
&kpt div1=16 div2=16 div3= 6 tkb=0.0005 /
```
#AsInp
```
As

&input film=f /

&lattice latsys='rho' a0=7.9729138 alpha=53.843117 /

2
33   0.22657  0.22657  0.22657
33  -0.22657 -0.22657 -0.22657

&atom element="As" rmt=2.3 jri=981 lmax=12 lnonsph=6  lo="3d" /
&comp kmax=5.0 gmaxxc=12.5 gmax=15.0 /
&kpt div1=13 div2=13 div3=13 tkb=0.0005 /
```
#S6Inp
```
S6 (ICSD # 40021)

&input film=f cartesian=t /

&lattice latsys='hR2', a0=1.889727, a=10.8180, c=4.2800 /

-1
16 0.0428 0.1882 0.1055

&gen 2                          ! as in the int. tables
  0 -1  0   0.0                 ! for R -3 H(148)
  1 -1  0   0.0                 ! http://www.cryst.ehu.es
  0  0  1   0.0

 -1  0  0   0.0                 ! we leave out the translation
  0 -1  0   0.0
  0  0 -1   0.0 /               ! (and identity op)
```
