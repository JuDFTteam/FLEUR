!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_vacuum

   TYPE t_vacuum
      !Stuff for the vacuum
      INTEGER ::nmz=250
      INTEGER ::nmzd=250
      INTEGER ::nmzxy=100
      INTEGER ::nmzxyd=100
      INTEGER :: layerd=1
      INTEGER :: layers=0
      INTEGER :: nvac=2
      INTEGER :: nvacd=2
      REAL :: delz=0.1
      REAL :: dvac=0.0
      INTEGER::nstars=0
      INTEGER:: nstm=0
      REAL :: tworkf=0.0
      REAL :: locx(2)=[0.,0.]
      REAL :: locy(2)=[0.,0.]
      LOGICAL ::starcoeff=.false.
      INTEGER, ALLOCATABLE :: izlay(:, :)= reshape([0.,0.],[1,2])
   END TYPE t_vacuum
 end MODULE m_types_vacuum
