!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_sphhar
   !Data for the spherical harmonics
  TYPE t_sphhar
     !No of symmetry types (must
     !equal maxval(atoms%ntypsy)
     INTEGER ::ntypsd
     !Max no of members of sphhar
     INTEGER ::memd
     !max of nlh
     INTEGER ::nlhd
     !No of sphhar (ntypsd)
     INTEGER,ALLOCATABLE ::nlh(:)
     !l's of sphhar (0:nlhd,ntypsd)
     INTEGER,ALLOCATABLE ::llh(:,:)
     !No of members in sphhar (0:nlh
     INTEGER,ALLOCATABLE ::nmem(:,:)
     !lm's of of members (max(nmem),
     INTEGER,ALLOCATABLE ::mlh(:,:,:)
     !phasefactors (max(nmem),0:nlhd
     COMPLEX,ALLOCATABLE ::clnu(:,:,:)
  END TYPE t_sphhar
END MODULE m_types_sphhar
