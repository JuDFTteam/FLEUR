!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_loccoeff

  ! coefficients used in locrectify:

  COMPLEX coeffi(1:28,1:28,0:13)
  LOGICAL coffkind(1:28,0:13)

  COMPLEX, PARAMETER ::  one = ( 1.0, 0.0)
  COMPLEX, PARAMETER :: mone = (-1.0, 0.0)
  COMPLEX, PARAMETER ::  zro = ( 0.0, 0.0)
  COMPLEX, PARAMETER ::  ci  = ( 0.0, 1.0)
  COMPLEX, PARAMETER :: mci  = ( 0.0,-1.0)

  !***************************************************************************
  !***************************************************************************
  !***************************************************************************
  ! ANGULAR MOMENTUM: 0
  !***************************************************************************
  !***************************************************************************
  !***************************************************************************

  !***************************************************************************
  !     linear combinations of locs that are eigenfunctions
  !     of z-reflection operation for l=0 (two atoms NOT in xy-plane)
  !***************************************************************************
  ! even parity => leave as is
  ! z-refl: even
  DATA coffkind(1,0) /.true./
  DATA coeffi(1:2,1,0) /one,one/
  ! odd parity => multiply with i
  ! z-refl: odd
  DATA coffkind(2,0) /.false./
  DATA coeffi(1:2,2,0) /ci,mci/


  !***************************************************************************
  !     linear combinations of locs that are eigenfunctions
  !     of z-reflection operation for l=0 (four atoms)
  !***************************************************************************
  ! even parity => leave as is
  ! z-refl: even
  DATA coffkind(1,6) /.true./
  DATA coeffi(1:4,1,6) /one,one,one,one/

  ! even parity => leave as is
  ! z-refl: odd
  DATA coffkind(2,6) /.false./
  DATA coeffi(1:4,2,6) /one,one,mone,mone/

  ! odd parity => multiply with i
  ! z-refl: even
  DATA coffkind(3,6) /.true./
  DATA coeffi(1:4,3,6) /ci,mci,ci,mci/

  ! odd parity => multiply with i
  ! z-refl: odd
  DATA coffkind(4,6) /.false./
  DATA coeffi(1:4,4,6) /ci,mci,mci,ci/

  !***************************************************************************
  !***************************************************************************
  !***************************************************************************
  ! ANGULAR MOMENTUM: 1
  !***************************************************************************
  !***************************************************************************
  !***************************************************************************


  !********************************************************
  !     linear combinations of locs that are eigenfunctions
  !     of z-reflection operation for l=1 (single atom)
  !********************************************************
  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//odd parity => multiply with i
  !z-refl: even
  DATA coffkind(1,1) /.true./
  DATA coeffi(1:3,1,1) /one,zro,one/

  !cos(theta) <==> Y10// odd parity => multiply with i
  !z-refl: odd
  DATA coffkind(2,1) /.false./
  DATA coeffi(1:3,2,1) /zro,ci,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// odd parity => multiply with i
  !z-refl: even
  DATA coffkind(3,1) /.true./
  DATA coeffi(1:3,3,1) /ci,zro,mci/

  !********************************************************
  !     linear combinations of locs for l=1 (two atoms; both atoms in xy-plane)
  !********************************************************
  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//even parity => leave as is
  !z-refl: even
  DATA coffkind(1,3) /.true./
  DATA coeffi(1:6,1,3) /ci,zro,ci,mci,zro,mci/

  !cos(theta) <==> Y10//odd parity => multiply with i
  !z-refl: odd
  DATA coffkind(2,3) /.false./
  DATA coeffi(1:6,2,3) /zro,ci,zro,zro,ci,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// odd parity => multiply with i
  !z-refl: even
  DATA coffkind(3,3) /.true./
  DATA coeffi(1:6,3,3) /ci,zro,mci,ci,zro,mci/

  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//odd parity => multiply with i
  !z-refl: even
  DATA coffkind(4,3) /.true./
  DATA coeffi(1:6,4,3) /one,zro,one,one,zro,one/

  !cos(theta) <==> Y10//even parity => leave as is
  !z-refl: odd
  DATA coffkind(5,3) /.false./
  DATA coeffi(1:6,5,3) /zro,one,zro,zro,mone,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// even parity => leave as is
  !z-refl: even
  DATA coffkind(6,3) /.true./
  DATA coeffi(1:6,6,3) /one,zro,mone,mone,zro,one/

  !********************************************************
  !     linear combinations of locs for l=1 (two atoms; atoms NOT in xy-plane)
  !********************************************************
  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(1,4) /.false./
  DATA coeffi(1:6,1,4) /ci,zro,ci,mci,zro,mci/

  !cos(theta) <==> Y10//odd parity => multiply with i
  !z-refl: odd
  DATA coffkind(2,4) /.false./
  DATA coeffi(1:6,2,4) /zro,ci,zro,zro,ci,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// odd parity => multiply with i
  !z-refl: even
  DATA coffkind(3,4) /.true./
  DATA coeffi(1:6,3,4) /ci,zro,mci,ci,zro,mci/

  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//odd parity => multiply with i
  !z-refl: even
  DATA coffkind(4,4) /.true./
  DATA coeffi(1:6,4,4) /one,zro,one,one,zro,one/

  !cos(theta) <==> Y10//even parity => leave as is
  !z-refl: even
  DATA coffkind(5,4) /.true./
  DATA coeffi(1:6,5,4) /zro,one,zro,zro,mone,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// even parity => leave as is
  !z-refl: odd
  DATA coffkind(6,4) /.false./
  DATA coeffi(1:6,6,4) /one,zro,mone,mone,zro,one/


  !********************************************************
  !     linear combinations of locs for l=1 (four atoms)
  !********************************************************
  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//even parity => leave as is
  !z-refl: even
  DATA coffkind(1,7) /.false./
  DATA coeffi(1:12,1,7) /ci,zro,ci,mci,zro,mci, ci,zro,ci,mci,zro,mci/

  !cos(theta) <==> Y10//odd parity => multiply with i
  !z-refl: odd
  DATA coffkind(2,7) /.false./
  DATA coeffi(1:12,2,7) /zro,ci,zro,zro,ci,zro, zro,ci,zro,zro,ci,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// odd parity => multiply with i
  !z-refl: even
  DATA coffkind(3,7) /.true./
  DATA coeffi(1:12,3,7) /ci,zro,mci,ci,zro,mci, ci,zro,mci,ci,zro,mci/

  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//odd parity => multiply with i
  !z-refl: odd
  DATA coffkind(4,7) /.false./
  DATA coeffi(1:12,4,7) / one,zro, one, one,zro, one, mone,zro,mone,mone,zro,mone/

  !cos(theta) <==> Y10//even parity => leave as is
  !z-refl: odd
  DATA coffkind(5,7) /.false./
  DATA coeffi(1:12,5,7) /zro, one,zro,zro,mone,zro, zro, one,zro,zro,mone,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// even parity => leave as is
  !z-refl: even
  DATA coffkind(6,7) /.true./
  DATA coeffi(1:12,6,7) /one,zro,mone,mone,zro,one, one,zro,mone,mone,zro,one/

  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(7,7) /.false./
  DATA coeffi(1:12,7,7) / ci,zro, ci,mci,zro,mci, mci,zro,mci, ci,zro, ci/

  !cos(theta) <==> Y10//odd parity => multiply with i
  !z-refl: even
  DATA coffkind(8,7) /.true./
  DATA coeffi(1:12,8,7) /zro, ci,zro,zro, ci,zro, zro,mci,zro,zro,mci,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// odd parity => multiply with i
  !z-refl: odd
  DATA coffkind(9,7) /.false./
  DATA coeffi(1:12,9,7) / ci,zro,mci, ci,zro,mci, mci,zro, ci,mci,zro, ci/

  !sin(phi)sin(theta) <==> i(Y11  +  Y1-1)//odd parity => multiply with i
  !z-refl: true
  DATA coffkind(10,7) /.true./
  DATA coeffi(1:12,10,7) /one,zro,one,one,zro,one, one,zro,one,one,zro,one/

  !cos(theta) <==> Y10//even parity => leave as is
  !z-refl: even
  DATA coffkind(11,7) /.true./
  DATA coeffi(1:12,11,7) /zro, one,zro,zro,mone,zro, zro,mone,zro,zro, one,zro/

  !cos(phi)sin(theta) <==>  (Y11  -  Y1-1)// even parity => leave as is
  !z-refl: odd
  DATA coffkind(12,7) /.false./
  DATA coeffi(1:12,12,7) / one,zro,mone,mone,zro, one, mone,zro, one, one,zro,mone/


  !***************************************************************************
  !***************************************************************************
  !***************************************************************************
  ! ANGULAR MOMENTUM: 2
  !***************************************************************************
  !***************************************************************************
  !***************************************************************************


  !********************************************************
  !     linear combinations of locs that are eigenfunctions
  !     of z-reflection operation for l=2 (single atom)
  !********************************************************
  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)// even parity => take as is
  !z-refl: even
  DATA coffkind(1,2) /.true./
  DATA coeffi(1:5,1,2) /one,zro,zro,zro,one/

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// even parity => take as is
  !z-refl: odd
  DATA coffkind(2,2) /.false./ 
  DATA coeffi(1:5,2,2) / zro ,one ,zro ,mone ,zro/

  !(3*cos^2(theta)-1) <==> Y20// even parity => take as is
  !z-refl: even
  DATA coffkind(3,2) /.true./ 
  DATA coeffi(1:5,3,2) / zro ,zro ,one ,zro ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// even parity => take as is
  !z-refl: odd
  DATA coffkind(4,2) /.false./ 
  DATA coeffi(1:5,4,2) / zro , ci ,zro , ci ,zro /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// even parity => take as is
  !z-refl: even
  DATA coffkind(5,2) /.true./ 
  DATA coeffi(1:5,5,2) /  ci ,zro ,zro ,zro ,mci /

  !******************************************************************************
  !     linear combinations of locs that are eigenfunctions
  !     of z-reflection operation for l=2; two atoms NOT in xy-plane 
  !******************************************************************************
  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)// even parity => take as is
  !z-refl: even
  DATA coffkind(1,8) /.true./ 
  DATA coeffi(1:10,1,8) / one ,zro ,zro ,zro ,one ,one ,zro ,zro ,zro ,one /

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// even parity => take as is
  !z-refl: odd
  DATA coffkind(2,8) /.false./ 
  DATA coeffi(1:10,2,8) / zro ,one ,zro ,mone ,zro ,zro ,one ,zro ,mone ,zro /

  !(3*cos^2(theta)-1) <==> Y20// even parity => take as is
  !z-refl: even
  DATA coffkind(3,8) /.true./ 
  DATA coeffi(1:10,3,8) / zro ,zro ,one ,zro ,zro ,zro ,zro ,one ,zro ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// even parity => take as is
  !z-refl: odd
  DATA coffkind(4,8) /.false./ 
  DATA coeffi(1:10,4,8) / zro , ci ,zro , ci ,zro ,zro , ci ,zro , ci ,zro /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// even parity => take as is
  !z-refl: even
  DATA coffkind(5,8) /.true./ 
  DATA coeffi(1:10,5,8) /  ci ,zro ,zro ,zro ,mci , ci ,zro ,zro ,zro ,mci /

  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)//odd parity => factor i
  !z-refl: odd
  DATA coffkind(6,8) /.false./ 
  DATA coeffi(1:10,6,8) /  ci ,zro ,zro ,zro , ci ,mci ,zro ,zro ,zro ,mci /

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// odd parity => factor i
  !z-refl: even
  DATA coffkind(7,8) /.true./ 
  DATA coeffi(1:10,7,8) / zro , ci ,zro ,mci ,zro ,zro ,mci ,zro , ci ,zro /

  !(3*cos^2(theta)-1) <==> Y20// odd parity => factor i
  !z-refl: odd
  DATA coffkind(8,8) /.false./ 
  DATA coeffi(1:10,8,8) / zro ,zro , ci ,zro ,zro ,zro ,zro ,mci ,zro ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// odd parity => factor i
  !z-refl: even
  DATA coffkind(9,8) /.true./ 
  DATA coeffi(1:10,9,8) / zro ,one ,zro ,one ,zro ,zro ,mone ,zro ,mone ,zro /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// odd parity => factor i
  !z-refl: odd
  DATA coffkind(10,8) /.false./ 
  DATA coeffi(1:10,10,8) / one ,zro ,zro ,zro ,mone ,mone ,zro ,zro ,zro ,one /

  !******************************************************************************
  !     linear combinations of locs that are eigenfunctions
  !     of z-reflection operation for l=2; two atoms IN xy-plane 
  !******************************************************************************
  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)// odd parity => factor i
  !z-refl: even
  DATA coffkind(1,9) /.true./ 
  DATA coeffi(1:10,1,9) /  ci ,zro ,zro ,zro , ci ,mci ,zro ,zro ,zro ,mci /

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// even parity => take as is
  !z-refl: odd
  DATA coffkind(2,9) /.false./ 
  DATA coeffi(1:10,2,9) / zro ,one ,zro ,mone ,zro ,zro ,one ,zro ,mone ,zro /

  !(3*cos^2(theta)-1) <==> Y20// even parity => take as is
  !z-refl: even
  DATA coffkind(3,9) /.true./ 
  DATA coeffi(1:10,3,9) / zro ,zro ,one ,zro ,zro ,zro ,zro ,one ,zro ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// even parity => take as is
  !z-refl: odd
  DATA coffkind(4,9) /.false./ 
  DATA coeffi(1:10,4,9) / zro , ci ,zro , ci ,zro ,zro , ci ,zro , ci ,zro /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// even parity => take as is
  !z-refl: even
  DATA coffkind(5,9) /.true./ 
  DATA coeffi(1:10,5,9) /  ci ,zro ,zro ,zro ,mci , ci ,zro ,zro ,zro ,mci /

  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)//even parity => leave as is
  !z-refl: even
  DATA coffkind(6,9) /.true./ 
  DATA coeffi(1:10,6,9) / one ,zro ,zro ,zro ,one ,one ,zro ,zro ,zro ,one /

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// odd parity => factor i
  !z-refl: odd
  DATA coffkind(7,9) /.false./ 
  DATA coeffi(1:10,7,9) / zro , ci ,zro ,mci ,zro ,zro ,mci ,zro , ci ,zro /

  !(3*cos^2(theta)-1) <==> Y20// odd parity => factor i
  !z-refl: even
  DATA coffkind(8,9) /.true./ 
  DATA coeffi(1:10,8,9) / zro ,zro , ci ,zro ,zro ,zro ,zro ,mci ,zro ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// odd parity => factor i
  !z-refl: odd
  DATA coffkind(9,9) /.false./ 
  DATA coeffi(1:10,9,9) / zro ,one ,zro ,one ,zro ,zro ,mone ,zro ,mone ,zro /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// odd parity => factor i
  !z-refl: even
  DATA coffkind(10,9) /.true./ 
  DATA coeffi(1:10,10,9) / one ,zro ,zro ,zro ,mone ,mone ,zro ,zro ,zro ,one /


  !******************************************************************************
  !     linear combinations of locs that are eigenfunctions
  !     of z-reflection operation for l=2; four entangled atoms
  !******************************************************************************
  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)// even parity => take as is
  !z-refl: even
  DATA coffkind(1,10) /.true./ 
  DATA coeffi(1:20,1,10) / one ,zro ,zro ,zro ,one ,one ,zro ,zro ,zro ,one&
       ,one ,zro ,zro ,zro ,one ,one ,zro ,zro ,zro ,one /

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// even parity => take as is
  !z-refl: odd
  DATA coffkind(2,10) /.false./ 
  DATA coeffi(1:20,2,10) / zro ,one ,zro ,mone ,zro ,zro ,one ,zro ,mone ,zro&
       ,zro ,one ,zro ,mone ,zro ,zro ,one ,zro ,mone ,zro /

  !(3*cos^2(theta)-1) <==> Y20// even parity => take as is
  !z-refl: even
  DATA coffkind(3,10) /.true./ 
  DATA coeffi(1:20,3,10) / zro ,zro ,one ,zro ,zro ,zro ,zro ,one ,zro ,zro&
       ,zro ,zro ,one ,zro ,zro ,zro ,zro ,one ,zro ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// even parity => take as is
  !z-refl: odd
  DATA coffkind(4,10) /.false./ 
  DATA coeffi(1:20,4,10) / zro , ci ,zro , ci ,zro ,zro , ci ,zro , ci ,zro&
       ,zro , ci ,zro , ci ,zro ,zro , ci ,zro , ci ,zro /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// even parity => take as is
  !z-refl: even
  DATA coffkind(5,10) /.true./ 
  DATA coeffi(1:20,5,10) /  ci ,zro ,zro ,zro ,mci , ci ,zro ,zro ,zro ,mci&
       , ci ,zro ,zro ,zro ,mci , ci ,zro ,zro ,zro ,mci /

  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)//odd parity => factor i
  !z-refl: odd
  DATA coffkind(6,10) /.false./ 
  DATA coeffi(1:20,6,10) /  ci ,zro ,zro ,zro , ci ,mci ,zro ,zro ,zro ,mci&
       ,mci ,zro ,zro ,zro ,mci , ci ,zro ,zro ,zro , ci /

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// odd parity => factor i
  !z-refl: odd
  DATA coffkind(7,10) /.false./ 
  DATA coeffi(1:20,7,10) / zro , ci ,zro ,mci ,zro ,zro ,mci ,zro , ci ,zro&
       ,zro , ci ,zro ,mci ,zro ,zro ,mci ,zro , ci ,zro /

  !(3*cos^2(theta)-1) <==> Y20// odd parity => factor i
  !z-refl: odd
  DATA coffkind(8,10) /.false./ 
  DATA coeffi(1:20,8,10) / zro ,zro , ci ,zro ,zro ,zro ,zro ,mci ,zro ,zro&
       ,zro ,zro ,mci ,zro ,zro ,zro ,zro , ci ,zro ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// odd parity => factor i
  !z-refl: even
  DATA coffkind(9,10) /.true./ 
  DATA coeffi(1:20,9,10) / zro ,one ,zro ,one ,zro ,zro ,mone ,zro ,mone ,zro&
       ,zro ,mone ,zro ,mone ,zro ,zro ,one ,zro ,one ,zro /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// odd parity => factor i
  !z-refl: odd
  DATA coffkind(10,10) /.false./ 
  DATA coeffi(1:20,10,10) / one ,zro ,zro ,zro ,mone ,mone ,zro ,zro ,zro ,one&
       ,mone ,zro ,zro ,zro ,one ,one ,zro ,zro ,zro ,mone /

  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)// even parity => take as is
  !z-refl: odd
  DATA coffkind(11,10) /.false./ 
  DATA coeffi(1:20,11,10) / one ,zro ,zro ,zro ,one ,one ,zro ,zro ,zro ,one&
       ,mone ,zro ,zro ,zro ,mone ,mone ,zro ,zro ,zro ,mone /

  !sin^2(theta)cos(2*phi) <==> (Y22 + Y2-2)//odd parity => factor i
  !z-refl: odd
  DATA coffkind(12,10) /.false./ 
  DATA coeffi(1:20,12,10) /  ci ,zro ,zro ,zro , ci ,mci ,zro ,zro ,zro ,mci&
       , ci ,zro ,zro ,zro , ci ,mci ,zro ,zro ,zro ,mci /

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// even parity => take as is
  !z-refl: even
  DATA coffkind(13,10) /.true./ 
  DATA coeffi(1:20,13,10) / zro ,one ,zro ,mone ,zro ,zro ,one ,zro ,mone ,zro&
       ,zro ,mone ,zro ,one ,zro ,zro ,mone ,zro ,one ,zro /

  !sin(theta)cos(theta)cos(phi) <==> (Y21 - Y2-1)// odd parity => factor i
  !z-refl: even
  DATA coffkind(14,10) /.true./ 
  DATA coeffi(1:20,14,10) / zro , ci ,zro ,mci ,zro ,zro ,mci ,zro , ci ,zro&
       ,zro ,mci ,zro , ci ,zro ,zro , ci ,zro ,mci ,zro /

  !(3*cos^2(theta)-1) <==> Y20// even parity => take as is
  !z-refl: odd
  DATA coffkind(15,10) /.false./ 
  DATA coeffi(1:20,15,10) / zro ,zro ,one ,zro ,zro ,zro ,zro ,one ,zro ,zro&
       ,zro ,zro ,mone ,zro ,zro ,zro ,zro ,mone ,zro ,zro /

  !(3*cos^2(theta)-1) <==> Y20// odd parity => factor i
  !z-refl: even
  DATA coffkind(16,10) /.true./ 
  DATA coeffi(1:20,16,10) / zro ,zro , ci ,zro ,zro ,zro ,zro ,mci ,zro ,zro&
       ,zro ,zro , ci ,zro ,zro ,zro ,zro ,mci ,zro ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// even parity => take as is
  !z-refl: even
  DATA coffkind(17,10) /.true./ 
  DATA coeffi(1:20,17,10) / zro , ci ,zro , ci ,zro ,zro , ci ,zro , ci ,zro&
       ,zro ,mci ,zro ,mci ,zro ,zro ,mci ,zro ,mci ,zro /

  !sin(theta)cos(theta)sin(phi) <==> i(Y21 + Y2-1)// odd parity => factor i
  !z-refl: odd
  DATA coffkind(18,10) /.false./ 
  DATA coeffi(1:20,18,10) / zro ,one ,zro ,one ,zro ,zro ,mone ,zro ,mone ,zro&
       ,zro ,one ,zro ,one ,zro ,zro ,mone ,zro ,mone ,zro /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// even parity => take as is
  !z-refl: odd
  DATA coffkind(19,10) /.false./ 
  DATA coeffi(1:20,19,10) /  ci ,zro ,zro ,zro ,mci , ci ,zro ,zro ,zro ,mci&
       ,mci ,zro ,zro ,zro , ci ,mci ,zro ,zro ,zro , ci /

  !sin^2(theta)sin(2*phi) <==> i(Y22 - Y2-2)// odd parity => factor i
  !z-refl: even
  DATA coffkind(20,10) /.true./ 
  DATA coeffi(1:20,20,10) / one ,zro ,zro ,zro ,mone ,mone ,zro ,zro ,zro ,one&
       ,one ,zro ,zro ,zro ,mone ,mone ,zro ,zro ,zro ,one /


  !***************************************************************************
  !***************************************************************************
  !***************************************************************************
  ! ANGULAR MOMENTUM: 3
  !***************************************************************************
  !***************************************************************************
  !***************************************************************************


  !*************************************************************
  !     linear combinations of locs that are eigenfunctions of
  !     z-reflection operation for l=3 (single locs)
  !*************************************************************
  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(1,5) /.true./ 
  DATA coeffi(1:7,1,5) / one ,zro ,zro ,zro ,zro ,zro ,one /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(2,5) /.false./ 
  DATA coeffi(1:7,2,5) / zro , ci ,zro ,zro ,zro , ci ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(3,5) /.true./ 
  DATA coeffi(1:7,3,5) / zro ,zro ,one ,zro ,one ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//odd parity=>multiply by i
  !z-refl: odd
  DATA coffkind(4,5) /.false./ 
  DATA coeffi(1:7,4,5) / zro ,zro ,zro , ci ,zro ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(5,5) /.true./ 
  DATA coeffi(1:7,5,5) / zro ,zro , ci ,zro ,mci ,zro ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(6,5) /.false./ 
  DATA coeffi(1:7,6,5) / zro ,one ,zro ,zro ,zro ,mone ,zro /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(7,5) /.true./ 
  DATA coeffi(1:7,7,5) /  ci ,zro ,zro ,zro ,zro ,zro ,mci /



  !***************************************************************************
  !     linear combinations of locs that are eigenfunctions of
  !     z-reflection operation for l=3; two locs NOT in xy-plane (double locs)
  !***************************************************************************
  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(1,11) /.true./ 
  DATA coeffi(1:14,1,11) / one ,zro ,zro ,zro ,zro ,zro ,one ,  one ,zro ,zro ,zro ,zro ,zro ,one /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(2,11) /.false./ 
  DATA coeffi(1:14,2,11) / zro , ci ,zro ,zro ,zro , ci ,zro ,  zro , ci ,zro ,zro ,zro , ci ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(3,11) /.true./ 
  DATA coeffi(1:14,3,11) / zro ,zro ,one ,zro ,one ,zro ,zro ,  zro ,zro ,one ,zro ,one ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//odd parity=>multiply by i
  !z-refl: odd
  DATA coffkind(4,11) /.false./ 
  DATA coeffi(1:14,4,11) / zro ,zro ,zro , ci ,zro ,zro ,zro ,  zro ,zro ,zro , ci ,zro ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(5,11) /.true./ 
  DATA coeffi(1:14,5,11) / zro ,zro , ci ,zro ,mci ,zro ,zro ,  zro ,zro , ci ,zro ,mci ,zro ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(6,11) /.false./ 
  DATA coeffi(1:14,6,11) / zro ,one ,zro ,zro ,zro ,mone ,zro ,  zro ,one ,zro ,zro ,zro ,mone ,zro /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(7,11) /.true./ 
  DATA coeffi(1:14,7,11) /  ci ,zro ,zro ,zro ,zro ,zro ,mci ,   ci ,zro ,zro ,zro ,zro ,zro ,mci /

  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(8,11) /.false./ 
  DATA coeffi(1:14,8,11) /  ci ,zro ,zro ,zro ,zro ,zro , ci ,  mci ,zro ,zro ,zro ,zro ,zro ,mci /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)//even parity => leave as is
  !z-refl: even
  DATA coffkind(9,11) /.true./ 
  DATA coeffi(1:14,9,11) / zro ,one ,zro ,zro ,zro ,one ,zro ,  zro ,mone ,zro ,zro ,zro ,mone ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(10,11) /.false./ 
  DATA coeffi(1:14,10,11) / zro ,zro , ci ,zro , ci ,zro ,zro ,  zro ,zro ,mci ,zro ,mci ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//even parity => leave as is
  !z-refl: even
  DATA coffkind(11,11) /.true./ 
  DATA coeffi(1:14,11,11) / zro ,zro ,zro ,one ,zro ,zro ,zro ,  zro ,zro ,zro ,mone ,zro ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(12,11) /.false./ 
  DATA coeffi(1:14,12,11) / zro ,zro ,one ,zro ,mone ,zro ,zro ,  zro ,zro ,mone ,zro ,one ,zro ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)//even parity => leave as is
  !z-refl: true
  DATA coffkind(13,11) /.true./ 
  DATA coeffi(1:14,13,11) / zro , ci ,zro ,zro ,zro ,mci ,zro ,  zro ,mci ,zro ,zro ,zro , ci ,zro /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(14,11) /.false./ 
  DATA coeffi(1:14,14,11) / one ,zro ,zro ,zro ,zro ,zro ,mone ,  mone ,zro ,zro ,zro ,zro ,zro ,one /

  !***************************************************************************
  !     linear combinations of locs that are eigenfunctions of
  !     z-reflection operation for l=3; two locs IN xy-plane (double locs)
  !***************************************************************************
  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(1,12) /.true./ 
  DATA coeffi(1:14,1,12) / one ,zro ,zro ,zro ,zro ,zro ,one ,  one ,zro ,zro ,zro ,zro ,zro ,one /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(2,12) /.false./ 
  DATA coeffi(1:14,2,12) / zro , ci ,zro ,zro ,zro , ci ,zro ,  zro , ci ,zro ,zro ,zro , ci ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(3,12) /.true./ 
  DATA coeffi(1:14,3,12) / zro ,zro ,one ,zro ,one ,zro ,zro ,  zro ,zro ,one ,zro ,one ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//odd parity=>multiply by i
  !z-refl: odd
  DATA coffkind(4,12) /.false./ 
  DATA coeffi(1:14,4,12) / zro ,zro ,zro , ci ,zro ,zro ,zro ,  zro ,zro ,zro , ci ,zro ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(5,12) /.true./ 
  DATA coeffi(1:14,5,12) / zro ,zro , ci ,zro ,mci ,zro ,zro ,  zro ,zro , ci ,zro ,mci ,zro ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(6,12) /.false./ 
  DATA coeffi(1:14,6,12) / zro ,one ,zro ,zro ,zro ,mone ,zro ,  zro ,one ,zro ,zro ,zro ,mone ,zro /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(7,12) /.true./ 
  DATA coeffi(1:14,7,12) /  ci ,zro ,zro ,zro ,zro ,zro ,mci ,   ci ,zro ,zro ,zro ,zro ,zro ,mci /

  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)//even parity => leave as is
  !z-refl: even
  DATA coffkind(8,12) /.true./ 
  DATA coeffi(1:14,8,12) /  ci ,zro ,zro ,zro ,zro ,zro , ci ,  mci ,zro ,zro ,zro ,zro ,zro ,mci /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(9,12) /.false./ 
  DATA coeffi(1:14,9,12) / zro ,one ,zro ,zro ,zro ,one ,zro ,  zro ,mone ,zro ,zro ,zro ,mone ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//even parity => leave as is
  !z-refl: even
  DATA coffkind(10,12) /.true./ 
  DATA coeffi(1:14,10,12) / zro ,zro , ci ,zro , ci ,zro ,zro ,  zro ,zro ,mci ,zro ,mci ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//even parity => leave as is
  !z-refl: odd
  DATA coffkind(11,12) /.false./ 
  DATA coeffi(1:14,11,12) / zro ,zro ,zro ,one ,zro ,zro ,zro ,  zro ,zro ,zro ,mone ,zro ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//even parity => leave as is
  !z-refl: even
  DATA coffkind(12,12) /.true./ 
  DATA coeffi(1:14,12,12) / zro ,zro ,one ,zro ,mone ,zro ,zro ,  zro ,zro ,mone ,zro ,one ,zro ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(13,12) /.false./ 
  DATA coeffi(1:14,13,12) / zro , ci ,zro ,zro ,zro ,mci ,zro ,  zro ,mci ,zro ,zro ,zro , ci ,zro /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)//even parity => leave as is
  !z-refl: even
  DATA coffkind(14,12) /.true./ 
  DATA coeffi(1:14,14,12) / one ,zro ,zro ,zro ,zro ,zro ,mone ,  mone ,zro ,zro ,zro ,zro ,zro ,one /

  !***************************************************************************
  !     linear combinations of locs that are eigenfunctions of
  !     z-reflection operation for l=3; four entangled atoms
  !***************************************************************************
  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(1,13) /.true./ 
  DATA coeffi(1:28,1,13) / one ,zro ,zro ,zro ,zro ,zro ,one ,  one ,zro ,zro ,zro ,zro ,zro ,one&
       , one ,zro ,zro ,zro ,zro ,zro ,one , one ,zro ,zro ,zro ,zro ,zro ,one /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(2,13) /.false./ 
  DATA coeffi(1:28,2,13) / zro , ci ,zro ,zro ,zro , ci ,zro ,  zro , ci ,zro ,zro ,zro , ci ,zro&
       , zro , ci ,zro ,zro ,zro , ci ,zro , zro , ci ,zro ,zro ,zro , ci ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(3,13) /.true./ 
  DATA coeffi(1:28,3,13) / zro ,zro ,one ,zro ,one ,zro ,zro ,  zro ,zro ,one ,zro ,one ,zro ,zro&
       , zro ,zro ,one ,zro ,one ,zro ,zro , zro ,zro ,one ,zro ,one ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//odd parity=>multiply by i
  !z-refl: odd
  DATA coffkind(4,13) /.false./ 
  DATA coeffi(1:28,4,13) / zro ,zro ,zro , ci ,zro ,zro ,zro ,  zro ,zro ,zro , ci ,zro ,zro ,zro&
       , zro ,zro ,zro , ci ,zro ,zro ,zro , zro ,zro ,zro , ci ,zro ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(5,13) /.true./ 
  DATA coeffi(1:28,5,13) / zro ,zro , ci ,zro ,mci ,zro ,zro ,  zro ,zro , ci ,zro ,mci ,zro ,zro&
       , zro ,zro , ci ,zro ,mci ,zro ,zro , zro ,zro , ci ,zro ,mci ,zro ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(6,13) /.false./ 
  DATA coeffi(1:28,6,13) / zro ,one ,zro ,zro ,zro ,mone ,zro ,  zro ,one ,zro ,zro ,zro ,mone ,zro&
       , zro ,one ,zro ,zro ,zro ,mone ,zro , zro ,one ,zro ,zro ,zro ,mone ,zro /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(7,13) /.true./ 
  DATA coeffi(1:28,7,13) /  ci ,zro ,zro ,zro ,zro ,zro ,mci ,   ci ,zro ,zro ,zro ,zro ,zro ,mci&
       ,  ci ,zro ,zro ,zro ,zro ,zro ,mci ,  ci ,zro ,zro ,zro ,zro ,zro ,mci /

  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)//even parity => leave as is
  !z-refl: even
  DATA coffkind(8,13) /.true./ 
  DATA coeffi(1:28,8,13) /  ci ,zro ,zro ,zro ,zro ,zro , ci ,  mci ,zro ,zro ,zro ,zro ,zro ,mci&
       ,  ci ,zro ,zro ,zro ,zro ,zro , ci , mci ,zro ,zro ,zro ,zro ,zro ,mci /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(9,13) /.false./ 
  DATA coeffi(1:28,9,13) / zro ,one ,zro ,zro ,zro ,one ,zro ,  zro ,mone ,zro ,zro ,zro ,mone ,zro&
       , zro ,one ,zro ,zro ,zro ,one ,zro , zro ,mone ,zro ,zro ,zro ,mone ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//even parity => leave as is
  !z-refl: even
  DATA coffkind(10,13) /.true./ 
  DATA coeffi(1:28,10,13) / zro ,zro , ci ,zro , ci ,zro ,zro ,  zro ,zro ,mci ,zro ,mci ,zro ,zro&
       , zro ,zro , ci ,zro , ci ,zro ,zro , zro ,zro ,mci ,zro ,mci ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//even parity => leave as is
  !z-refl: odd
  DATA coffkind(11,13) /.false./ 
  DATA coeffi(1:28,11,13) / zro ,zro ,zro ,one ,zro ,zro ,zro ,  zro ,zro ,zro ,mone ,zro ,zro ,zro&
       , zro ,zro ,zro ,one ,zro ,zro ,zro , zro ,zro ,zro ,mone ,zro ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//even parity => leave as is
  !z-refl: even
  DATA coffkind(12,13) /.true./ 
  DATA coeffi(1:28,12,13) / zro ,zro ,one ,zro ,mone ,zro ,zro ,  zro ,zro ,mone ,zro ,one ,zro ,zro&
       , zro ,zro ,one ,zro ,mone ,zro ,zro , zro ,zro ,mone ,zro ,one ,zro ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(13,13) /.false./ 
  DATA coeffi(1:28,13,13) / zro , ci ,zro ,zro ,zro ,mci ,zro ,  zro ,mci ,zro ,zro ,zro , ci ,zro&
       , zro , ci ,zro ,zro ,zro ,mci ,zro , zro ,mci ,zro ,zro ,zro , ci ,zro /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)//even parity => leave as is
  !z-refl: even
  DATA coffkind(14,13) /.true./ 
  DATA coeffi(1:28,14,13) / one ,zro ,zro ,zro ,zro ,zro ,mone ,  mone ,zro ,zro ,zro ,zro ,zro ,one&
       , one ,zro ,zro ,zro ,zro ,zro ,mone , mone ,zro ,zro ,zro ,zro ,zro ,one /

  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(15,13) /.false./ 
  DATA coeffi(1:28,15,13) / one ,zro ,zro ,zro ,zro ,zro ,one ,  one ,zro ,zro ,zro ,zro ,zro ,one&
       , mone ,zro ,zro ,zro ,zro ,zro ,mone , mone ,zro ,zro ,zro ,zro ,zro ,mone /

  !sin^3(theta)sin(3*phi) <==> i(Y33+Y3-3)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(16,13) /.false./ 
  DATA coeffi(1:28,16,13) /  ci ,zro ,zro ,zro ,zro ,zro , ci ,  mci ,zro ,zro ,zro ,zro ,zro ,mci&
       , mci ,zro ,zro ,zro ,zro ,zro ,mci ,  ci ,zro ,zro ,zro ,zro ,zro , ci /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(17,13) /.true./ 
  DATA coeffi(1:28,17,13) / zro , ci ,zro ,zro ,zro , ci ,zro ,  zro , ci ,zro ,zro ,zro , ci ,zro&
       , zro ,mci ,zro ,zro ,zro ,mci ,zro , zro ,mci ,zro ,zro ,zro ,mci ,zro /

  !sin^2(theta)cos(theta)cos(2*phi)<==>(Y32+Y3-2)//even parity => leave as is
  !z-refl: even
  DATA coffkind(18,13) /.true./ 
  DATA coeffi(1:28,18,13) / zro ,one ,zro ,zro ,zro ,one ,zro ,  zro ,mone ,zro ,zro ,zro ,mone ,zro&
       , zro ,mone ,zro ,zro ,zro ,mone ,zro , zro ,one ,zro ,zro ,zro ,one ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//odd parity=>multiply by i
  !z-refl: odd
  DATA coffkind(19,13) /.false./ 
  DATA coeffi(1:28,19,13) / zro ,zro ,one ,zro ,one ,zro ,zro ,  zro ,zro ,one ,zro ,one ,zro ,zro&
       , zro ,zro ,mone ,zro ,mone ,zro ,zro , zro ,zro ,mone ,zro ,mone ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]sin(phi)<==>i(Y31+Y3-1)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(20,13) /.false./ 
  DATA coeffi(1:28,20,13) / zro ,zro , ci ,zro , ci ,zro ,zro ,  zro ,zro ,mci ,zro ,mci ,zro ,zro&
       , zro ,zro ,mci ,zro ,mci ,zro ,zro , zro ,zro , ci ,zro , ci ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//odd parity=>multiply by i
  !z-refl: even
  DATA coffkind(21,13) /.true./ 
  DATA coeffi(1:28,21,13) / zro ,zro ,zro , ci ,zro ,zro ,zro ,  zro ,zro ,zro , ci ,zro ,zro ,zro&
       , zro ,zro ,zro ,mci ,zro ,zro ,zro , zro ,zro ,zro ,mci ,zro ,zro ,zro /

  ![5*cos^3(theta)-3*cos(theta)]<==>Y30//even parity => leave as is
  !z-refl: even
  DATA coffkind(22,13) /.true./ 
  DATA coeffi(1:28,22,13) / zro ,zro ,zro ,one ,zro ,zro ,zro ,  zro ,zro ,zro ,mone ,zro ,zro ,zro&
       , zro ,zro ,zro ,mone ,zro ,zro ,zro , zro ,zro ,zro ,one ,zro ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//odd parity=>multiply by i
  !z-refl: odd
  DATA coffkind(23,13) /.false./ 
  DATA coeffi(1:28,23,13) / zro ,zro , ci ,zro ,mci ,zro ,zro ,  zro ,zro , ci ,zro ,mci ,zro ,zro&
       , zro ,zro ,mci ,zro , ci ,zro ,zro , zro ,zro ,mci ,zro , ci ,zro ,zro /

  !sin(theta)[5*cos^2(theta)-1]cos(phi)<==>(Y31-Y3-1)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(24,13) /.false./ 
  DATA coeffi(1:28,24,13) / zro ,zro ,one ,zro ,mone ,zro ,zro ,  zro ,zro ,mone ,zro ,one ,zro ,zro&
       , zro ,zro ,mone ,zro ,one ,zro ,zro , zro ,zro ,one ,zro ,mone ,zro ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)// odd parity => multiply by i
  !z-refl: even
  DATA coffkind(25,13) /.true./ 
  DATA coeffi(1:28,25,13) / zro ,one ,zro ,zro ,zro ,mone ,zro ,  zro ,one ,zro ,zro ,zro ,mone ,zro&
       , zro ,mone ,zro ,zro ,zro ,one ,zro , zro ,mone ,zro ,zro ,zro ,one ,zro /

  !sin^2(theta)cos(theta)sin(2*phi)<==>i(Y32-Y3-2)//even parity => leave as is
  !z-refl: even
  DATA coffkind(26,13) /.true./ 
  DATA coeffi(1:28,26,13) / zro , ci ,zro ,zro ,zro ,mci ,zro ,  zro ,mci ,zro ,zro ,zro , ci ,zro&
       , zro ,mci ,zro ,zro ,zro , ci ,zro , zro , ci ,zro ,zro ,zro ,mci ,zro /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)// odd parity => multiply by i
  !z-refl: odd
  DATA coffkind(27,13) /.false./ 
  DATA coeffi(1:28,27,13) /  ci ,zro ,zro ,zro ,zro ,zro ,mci ,   ci ,zro ,zro ,zro ,zro ,zro ,mci&
       , mci ,zro ,zro ,zro ,zro ,zro , ci , mci ,zro ,zro ,zro ,zro ,zro , ci /

  !sin^3(theta)cos(3*phi) <==> (Y33 - Y3-3)//even parity => leave as is
  !z-refl: odd
  DATA coffkind(28,13) /.false./ 
  DATA coeffi(1:28,28,13) / one ,zro ,zro ,zro ,zro ,zro ,mone ,  mone ,zro ,zro ,zro ,zro ,zro ,one&
       , mone ,zro ,zro ,zro ,zro ,zro ,one , one ,zro ,zro ,zro ,zro ,zro ,mone /

END MODULE m_loccoeff

