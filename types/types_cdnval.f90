!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_cdnval

IMPLICIT NONE

PRIVATE

   TYPE t_orb
      REAL, ALLOCATABLE    :: uu(:,:,:,:)
      REAL, ALLOCATABLE    :: dd(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uup(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uum(:,:,:,:)
      COMPLEX, ALLOCATABLE :: ddp(:,:,:,:)
      COMPLEX, ALLOCATABLE :: ddm(:,:,:,:)

      REAL, ALLOCATABLE    :: uulo(:,:,:,:)
      REAL, ALLOCATABLE    :: dulo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uulop(:,:,:,:)
      COMPLEX, ALLOCATABLE :: uulom(:,:,:,:)
      COMPLEX, ALLOCATABLE :: dulop(:,:,:,:)
      COMPLEX, ALLOCATABLE :: dulom(:,:,:,:)

      REAL, ALLOCATABLE    :: z(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: p(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: m(:,:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init => orb_init
   END TYPE t_orb

   TYPE t_denCoeffs
      ! spherical
      REAL, ALLOCATABLE    :: uu(:,:,:)
      REAL, ALLOCATABLE    :: dd(:,:,:)
      REAL, ALLOCATABLE    :: du(:,:,:)

      ! nonspherical
      REAL, ALLOCATABLE    :: uunmt(:,:,:,:)
      REAL, ALLOCATABLE    :: ddnmt(:,:,:,:)
      REAL, ALLOCATABLE    :: dunmt(:,:,:,:)
      REAL, ALLOCATABLE    :: udnmt(:,:,:,:)

      ! spherical - LOs
      REAL, ALLOCATABLE    :: aclo(:,:,:)
      REAL, ALLOCATABLE    :: bclo(:,:,:)
      REAL, ALLOCATABLE    :: cclo(:,:,:,:)

      ! nonspherical - LOs
      REAL, ALLOCATABLE    :: acnmt(:,:,:,:,:)
      REAL, ALLOCATABLE    :: bcnmt(:,:,:,:,:)
      REAL, ALLOCATABLE    :: ccnmt(:,:,:,:,:)


      CONTAINS
      PROCEDURE,PASS :: init => denCoeffs_init
   END TYPE t_denCoeffs

   TYPE t_denCoeffsOffdiag
      ! spherical
      COMPLEX, ALLOCATABLE :: uu21(:,:)
      COMPLEX, ALLOCATABLE :: dd21(:,:)
      COMPLEX, ALLOCATABLE :: du21(:,:)
      COMPLEX, ALLOCATABLE :: ud21(:,:)

      ! nonspherical
      COMPLEX, ALLOCATABLE :: uunmt21(:,:,:)
      COMPLEX, ALLOCATABLE :: ddnmt21(:,:,:)
      COMPLEX, ALLOCATABLE :: dunmt21(:,:,:)
      COMPLEX, ALLOCATABLE :: udnmt21(:,:,:)

      ! spherical - LOs
      COMPLEX, ALLOCATABLE :: uulo21(:,:)
      COMPLEX, ALLOCATABLE :: dulo21(:,:)
      COMPLEX, ALLOCATABLE :: ulou21(:,:)
      COMPLEX, ALLOCATABLE :: ulod21(:,:)

      COMPLEX, ALLOCATABLE :: uloulop21(:,:,:)

      ! norms
      REAL, ALLOCATABLE     :: uu21n(:,:)
      REAL, ALLOCATABLE     :: ud21n(:,:)
      REAL, ALLOCATABLE     :: du21n(:,:)
      REAL, ALLOCATABLE     :: dd21n(:,:)

      REAL, ALLOCATABLE     :: uulo21n(:,:)
      REAL, ALLOCATABLE     :: dulo21n(:,:)
      REAL, ALLOCATABLE     :: ulou21n(:,:)
      REAL, ALLOCATABLE     :: ulod21n(:,:)

      REAL, ALLOCATABLE     :: uloulop21n(:,:,:)

      CONTAINS
      PROCEDURE,PASS :: init => denCoeffsOffdiag_init
   END TYPE t_denCoeffsOffdiag

   TYPE t_force
      COMPLEX, ALLOCATABLE :: f_a12(:,:)
      COMPLEX, ALLOCATABLE :: f_a21(:,:)
      COMPLEX, ALLOCATABLE :: f_b4(:,:)
      COMPLEX, ALLOCATABLE :: f_b8(:,:)

      COMPLEX, ALLOCATABLE :: e1cof(:,:,:)
      COMPLEX, ALLOCATABLE :: e2cof(:,:,:)
      COMPLEX, ALLOCATABLE :: aveccof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: bveccof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: cveccof(:,:,:,:,:)

      COMPLEX, ALLOCATABLE :: acoflo(:,:,:,:)
      COMPLEX, ALLOCATABLE :: bcoflo(:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init1 => force_init1
         PROCEDURE,PASS :: init2 => force_init2
   END TYPE t_force

   TYPE t_slab
      INTEGER              :: nsld, nsl

      INTEGER, ALLOCATABLE :: nmtsl(:,:)
      INTEGER, ALLOCATABLE :: nslat(:,:)
      REAL,    ALLOCATABLE :: zsl(:,:)
      REAL,    ALLOCATABLE :: volsl(:)
      REAL,    ALLOCATABLE :: volintsl(:)
      REAL,    ALLOCATABLE :: qintsl(:,:)
      REAL,    ALLOCATABLE :: qmtsl(:,:)

      CONTAINS
         PROCEDURE,PASS :: init => slab_init
   END TYPE t_slab

   TYPE t_eigVecCoeffs
      COMPLEX, ALLOCATABLE :: acof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: bcof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: ccof(:,:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init => eigVecCoeffs_init
   END TYPE t_eigVecCoeffs

   TYPE t_mcd
      REAL                 :: emcd_lo, emcd_up

      INTEGER, ALLOCATABLE :: ncore(:)
      REAL,    ALLOCATABLE :: e_mcd(:,:,:)
      REAL,    ALLOCATABLE :: mcd(:,:,:)
      COMPLEX, ALLOCATABLE :: m_mcd(:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init1 => mcd_init1
   END TYPE t_mcd

PUBLIC t_orb, t_denCoeffs, t_denCoeffsOffdiag, t_force, t_slab, t_eigVecCoeffs, t_mcd

CONTAINS

SUBROUTINE orb_init(thisOrb, atoms, noco, jsp_start, jsp_end)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_orb), INTENT(INOUT)    :: thisOrb
   TYPE(t_atoms), INTENT(IN)      :: atoms
   TYPE(t_noco), INTENT(IN)       :: noco
   INTEGER, INTENT(IN)            :: jsp_start
   INTEGER, INTENT(IN)            :: jsp_end

   INTEGER                        :: dim1, dim2, dim3

   IF(ALLOCATED(thisOrb%uu)) DEALLOCATE(thisOrb%uu)
   IF(ALLOCATED(thisOrb%dd)) DEALLOCATE(thisOrb%dd)
   IF(ALLOCATED(thisOrb%uup)) DEALLOCATE(thisOrb%uup)
   IF(ALLOCATED(thisOrb%uum)) DEALLOCATE(thisOrb%uum)
   IF(ALLOCATED(thisOrb%ddp)) DEALLOCATE(thisOrb%ddp)
   IF(ALLOCATED(thisOrb%ddm)) DEALLOCATE(thisOrb%ddm)

   IF(ALLOCATED(thisOrb%uulo)) DEALLOCATE(thisOrb%uulo)
   IF(ALLOCATED(thisOrb%dulo)) DEALLOCATE(thisOrb%dulo)
   IF(ALLOCATED(thisOrb%uulop)) DEALLOCATE(thisOrb%uulop)
   IF(ALLOCATED(thisOrb%uulom)) DEALLOCATE(thisOrb%uulom)
   IF(ALLOCATED(thisOrb%dulop)) DEALLOCATE(thisOrb%dulop)
   IF(ALLOCATED(thisOrb%dulom)) DEALLOCATE(thisOrb%dulom)

   IF(ALLOCATED(thisOrb%z)) DEALLOCATE(thisOrb%z)
   IF(ALLOCATED(thisOrb%p)) DEALLOCATE(thisOrb%p)
   IF(ALLOCATED(thisOrb%m)) DEALLOCATE(thisOrb%m)

   dim1 = 0
   dim2 = 1
   dim3 = 1
   IF (noco%l_soc) THEN
      dim1 = atoms%lmaxd
      dim2 = atoms%ntype
      dim3 = atoms%nlod
   END IF

   ALLOCATE(thisOrb%uu(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%dd(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%uup(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%uum(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%ddp(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%ddm(0:dim1,-atoms%lmaxd:atoms%lmaxd,dim2,jsp_start:jsp_end))

   ALLOCATE(thisOrb%uulo(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%dulo(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%uulop(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%uulom(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%dulop(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%dulom(dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))

   ALLOCATE(thisOrb%z(dim3,dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%p(dim3,dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))
   ALLOCATE(thisOrb%m(dim3,dim3,-atoms%llod:atoms%llod,dim2,jsp_start:jsp_end))

   thisOrb%uu = 0.0
   thisOrb%dd = 0.0
   thisOrb%uup = CMPLX(0.0,0.0)
   thisOrb%uum = CMPLX(0.0,0.0)
   thisOrb%ddp = CMPLX(0.0,0.0)
   thisOrb%ddm = CMPLX(0.0,0.0)

   thisOrb%uulo = 0.0
   thisOrb%dulo = 0.0
   thisOrb%uulop = CMPLX(0.0,0.0)
   thisOrb%uulom = CMPLX(0.0,0.0)
   thisOrb%dulop = CMPLX(0.0,0.0)
   thisOrb%dulom = CMPLX(0.0,0.0)

   thisOrb%z = 0.0
   thisOrb%p = CMPLX(0.0,0.0)
   thisOrb%m = CMPLX(0.0,0.0)

END SUBROUTINE orb_init

SUBROUTINE denCoeffs_init(thisDenCoeffs, atoms, sphhar, jsp_start, jsp_end)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_denCoeffs), INTENT(INOUT) :: thisDenCoeffs
   TYPE(t_atoms),      INTENT(IN)    :: atoms
   TYPE(t_sphhar),     INTENT(IN)    :: sphhar
   INTEGER,            INTENT(IN)    :: jsp_start
   INTEGER,            INTENT(IN)    :: jsp_end

   INTEGER                           :: llpd

   llpd = (atoms%lmaxd*(atoms%lmaxd+3)) / 2

   ALLOCATE (thisDenCoeffs%uu(0:atoms%lmaxd,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%dd(0:atoms%lmaxd,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%du(0:atoms%lmaxd,atoms%ntype,jsp_start:jsp_end))

   ALLOCATE (thisDenCoeffs%uunmt(0:llpd,sphhar%nlhd,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%ddnmt(0:llpd,sphhar%nlhd,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%dunmt(0:llpd,sphhar%nlhd,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%udnmt(0:llpd,sphhar%nlhd,atoms%ntype,jsp_start:jsp_end))

   ALLOCATE (thisDenCoeffs%aclo(atoms%nlod,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%bclo(atoms%nlod,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%cclo(atoms%nlod,atoms%nlod,atoms%ntype,jsp_start:jsp_end))

   ALLOCATE (thisDenCoeffs%acnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%bcnmt(0:atoms%lmaxd,atoms%nlod,sphhar%nlhd,atoms%ntype,jsp_start:jsp_end))
   ALLOCATE (thisDenCoeffs%ccnmt(atoms%nlod,atoms%nlod,sphhar%nlhd,atoms%ntype,jsp_start:jsp_end))

   thisDenCoeffs%uu = 0.0
   thisDenCoeffs%dd = 0.0
   thisDenCoeffs%du = 0.0

   thisDenCoeffs%uunmt = 0.0
   thisDenCoeffs%ddnmt = 0.0
   thisDenCoeffs%dunmt = 0.0
   thisDenCoeffs%udnmt = 0.0

   thisDenCoeffs%aclo = 0.0
   thisDenCoeffs%bclo = 0.0
   thisDenCoeffs%cclo = 0.0

   thisDenCoeffs%acnmt = 0.0
   thisDenCoeffs%bcnmt = 0.0
   thisDenCoeffs%ccnmt = 0.0

END SUBROUTINE denCoeffs_init

SUBROUTINE denCoeffsOffdiag_init(thisDenCoeffsOffdiag, atoms, noco, sphhar, l_fmpl)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_denCoeffsOffdiag), INTENT(INOUT) :: thisDenCoeffsOffdiag
   TYPE(t_atoms),      INTENT(IN)    :: atoms
   TYPE(t_noco),       INTENT(IN)    :: noco
   TYPE(t_sphhar),     INTENT(IN)    :: sphhar
   LOGICAL,            INTENT(IN)    :: l_fmpl

   IF (noco%l_mperp) THEN
      ALLOCATE (thisDenCoeffsOffdiag%uu21(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ud21(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%du21(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dd21(0:atoms%lmaxd,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uulo21(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dulo21(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ulou21(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ulod21(atoms%nlod,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uloulop21(atoms%nlod,atoms%nlod,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uu21n(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ud21n(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%du21n(0:atoms%lmaxd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dd21n(0:atoms%lmaxd,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uulo21n(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dulo21n(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ulou21n(atoms%nlod,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ulod21n(atoms%nlod,atoms%ntype))

      ALLOCATE (thisDenCoeffsOffdiag%uloulop21n(atoms%nlod,atoms%nlod,atoms%ntype))
   ELSE
      ALLOCATE (thisDenCoeffsOffdiag%uu21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ud21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%du21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dd21(1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uulo21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dulo21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ulou21(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ulod21(1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uloulop21(1,1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uu21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ud21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%du21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dd21n(1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uulo21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dulo21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ulou21n(1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ulod21n(1,1))

      ALLOCATE (thisDenCoeffsOffdiag%uloulop21n(1,1,1))
   END IF

   IF (noco%l_mperp.AND.l_fmpl) THEN
      ALLOCATE (thisDenCoeffsOffdiag%uunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%udnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%dunmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype))
      ALLOCATE (thisDenCoeffsOffdiag%ddnmt21((atoms%lmaxd+1)**2,sphhar%nlhd,atoms%ntype))
   ELSE
      ALLOCATE (thisDenCoeffsOffdiag%uunmt21(1,1,1))
      ALLOCATE (thisDenCoeffsOffdiag%udnmt21(1,1,1))
      ALLOCATE (thisDenCoeffsOffdiag%dunmt21(1,1,1))
      ALLOCATE (thisDenCoeffsOffdiag%ddnmt21(1,1,1))
   END IF

   thisDenCoeffsOffdiag%uu21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%ud21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%du21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%dd21 = CMPLX(0.0,0.0)

   thisDenCoeffsOffdiag%uulo21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%dulo21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%ulou21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%ulod21 = CMPLX(0.0,0.0)

   thisDenCoeffsOffdiag%uloulop21 = CMPLX(0.0,0.0)

   thisDenCoeffsOffdiag%uu21n = 0.0
   thisDenCoeffsOffdiag%ud21n = 0.0
   thisDenCoeffsOffdiag%du21n = 0.0
   thisDenCoeffsOffdiag%dd21n = 0.0

   thisDenCoeffsOffdiag%uulo21n = 0.0
   thisDenCoeffsOffdiag%dulo21n = 0.0
   thisDenCoeffsOffdiag%ulou21n = 0.0
   thisDenCoeffsOffdiag%ulod21n = 0.0

   thisDenCoeffsOffdiag%uloulop21n = 0.0

   thisDenCoeffsOffdiag%uunmt21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%udnmt21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%dunmt21 = CMPLX(0.0,0.0)
   thisDenCoeffsOffdiag%ddnmt21 = CMPLX(0.0,0.0)

END SUBROUTINE denCoeffsOffdiag_init

SUBROUTINE force_init1(thisForce,input,atoms)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_force),     INTENT(INOUT) :: thisForce
   TYPE(t_input),      INTENT(IN)    :: input
   TYPE(t_atoms),      INTENT(IN)    :: atoms

   IF (input%l_f) THEN
      ALLOCATE (thisForce%f_a12(3,atoms%ntype))
      ALLOCATE (thisForce%f_a21(3,atoms%ntype))
      ALLOCATE (thisForce%f_b4(3,atoms%ntype))
      ALLOCATE (thisForce%f_b8(3,atoms%ntype))
   ELSE
      ALLOCATE (thisForce%f_a12(1,1))
      ALLOCATE (thisForce%f_a21(1,1))
      ALLOCATE (thisForce%f_b4(1,1))
      ALLOCATE (thisForce%f_b8(1,1))
   END IF

   thisForce%f_a12 = CMPLX(0.0,0.0)
   thisForce%f_a21 = CMPLX(0.0,0.0)
   thisForce%f_b4 = CMPLX(0.0,0.0)
   thisForce%f_b8 = CMPLX(0.0,0.0)

END SUBROUTINE force_init1

SUBROUTINE force_init2(thisForce,noccbd,input,atoms)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_force),     INTENT(INOUT) :: thisForce
   TYPE(t_input),      INTENT(IN)    :: input
   TYPE(t_atoms),      INTENT(IN)    :: atoms
   INTEGER,            INTENT(IN)    :: noccbd

   IF (ALLOCATED(thisForce%e1cof)) DEALLOCATE(thisForce%e1cof)
   IF (ALLOCATED(thisForce%e2cof)) DEALLOCATE(thisForce%e2cof)
   IF (ALLOCATED(thisForce%acoflo)) DEALLOCATE(thisForce%acoflo)
   IF (ALLOCATED(thisForce%bcoflo)) DEALLOCATE(thisForce%bcoflo)
   IF (ALLOCATED(thisForce%aveccof)) DEALLOCATE(thisForce%aveccof)
   IF (ALLOCATED(thisForce%bveccof)) DEALLOCATE(thisForce%bveccof)
   IF (ALLOCATED(thisForce%cveccof)) DEALLOCATE(thisForce%cveccof)

   IF (input%l_f) THEN
      ALLOCATE (thisForce%e1cof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat))
      ALLOCATE (thisForce%e2cof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat))
      ALLOCATE (thisForce%acoflo(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat))
      ALLOCATE (thisForce%bcoflo(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat))
      ALLOCATE (thisForce%aveccof(3,noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat))
      ALLOCATE (thisForce%bveccof(3,noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat))
      ALLOCATE (thisForce%cveccof(3,-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat))
   ELSE
      ALLOCATE (thisForce%e1cof(1,1,1))
      ALLOCATE (thisForce%e2cof(1,1,1))
      ALLOCATE (thisForce%acoflo(1,1,1,1))
      ALLOCATE (thisForce%bcoflo(1,1,1,1))
      ALLOCATE (thisForce%aveccof(1,1,1,1))
      ALLOCATE (thisForce%bveccof(1,1,1,1))
      ALLOCATE (thisForce%cveccof(1,1,1,1,1))
   END IF

   thisForce%e1cof = CMPLX(0.0,0.0)
   thisForce%e2cof = CMPLX(0.0,0.0)
   thisForce%acoflo = CMPLX(0.0,0.0)
   thisForce%bcoflo = CMPLX(0.0,0.0)
   thisForce%aveccof = CMPLX(0.0,0.0)
   thisForce%bveccof = CMPLX(0.0,0.0)
   thisForce%cveccof = CMPLX(0.0,0.0)

END SUBROUTINE force_init2

SUBROUTINE slab_init(thisSlab,banddos,dimension,atoms,cell)

   USE m_types_setup
   USE m_slabdim
   USE m_slabgeom

   IMPLICIT NONE

   CLASS(t_slab),      INTENT(INOUT) :: thisSlab
   TYPE(t_banddos),    INTENT(IN)    :: banddos
   TYPE(t_dimension),  INTENT(IN)    :: dimension
   TYPE(t_atoms),      INTENT(IN)    :: atoms
   TYPE(t_cell),       INTENT(IN)    :: cell

   INTEGER :: nsld

   nsld=1

   IF ((banddos%ndir.EQ.-3).AND.banddos%dos) THEN
      CALL slab_dim(atoms, nsld)
      ALLOCATE (thisSlab%nmtsl(atoms%ntype,nsld))
      ALLOCATE (thisSlab%nslat(atoms%nat,nsld))
      ALLOCATE (thisSlab%zsl(2,nsld))
      ALLOCATE (thisSlab%volsl(nsld))
      ALLOCATE (thisSlab%volintsl(nsld))
      ALLOCATE (thisSlab%qintsl(nsld,dimension%neigd))
      ALLOCATE (thisSlab%qmtsl(nsld,dimension%neigd))
      CALL slabgeom(atoms,cell,nsld,thisSlab%nsl,thisSlab%zsl,thisSlab%nmtsl,&
                    thisSlab%nslat,thisSlab%volsl,thisSlab%volintsl)
   ELSE
      ALLOCATE (thisSlab%nmtsl(1,1))
      ALLOCATE (thisSlab%nslat(1,1))
      ALLOCATE (thisSlab%zsl(1,1))
      ALLOCATE (thisSlab%volsl(1))
      ALLOCATE (thisSlab%volintsl(1))
      ALLOCATE (thisSlab%qintsl(1,1))
      ALLOCATE (thisSlab%qmtsl(1,1))
   END IF
   thisSlab%nsld = nsld

END SUBROUTINE slab_init


SUBROUTINE eigVecCoeffs_init(thisEigVecCoeffs,dimension,atoms,noco,jspin,noccbd)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_eigVecCoeffs), INTENT(INOUT) :: thisEigVecCoeffs
   TYPE(t_dimension),     INTENT(IN)    :: dimension
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_noco),          INTENT(IN)    :: noco

   INTEGER,               INTENT(IN)    :: jspin, noccbd

   IF(ALLOCATED(thisEigVecCoeffs%acof)) DEALLOCATE(thisEigVecCoeffs%acof)
   IF(ALLOCATED(thisEigVecCoeffs%bcof)) DEALLOCATE(thisEigVecCoeffs%bcof)
   IF(ALLOCATED(thisEigVecCoeffs%ccof)) DEALLOCATE(thisEigVecCoeffs%ccof)

   IF (noco%l_mperp) THEN
      ALLOCATE (thisEigVecCoeffs%acof(noccbd,0:dimension%lmd,atoms%nat,dimension%jspd))
      ALLOCATE (thisEigVecCoeffs%bcof(noccbd,0:dimension%lmd,atoms%nat,dimension%jspd))
      ALLOCATE (thisEigVecCoeffs%ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat,dimension%jspd))
   ELSE
      ALLOCATE (thisEigVecCoeffs%acof(noccbd,0:dimension%lmd,atoms%nat,jspin:jspin))
      ALLOCATE (thisEigVecCoeffs%bcof(noccbd,0:dimension%lmd,atoms%nat,jspin:jspin))
      ALLOCATE (thisEigVecCoeffs%ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat,jspin:jspin))
   END IF

   thisEigVecCoeffs%acof = CMPLX(0.0,0.0)
   thisEigVecCoeffs%bcof = CMPLX(0.0,0.0)
   thisEigVecCoeffs%ccof = CMPLX(0.0,0.0)

END SUBROUTINE eigVecCoeffs_init

SUBROUTINE mcd_init1(thisMCD,banddos,dimension,input,atoms)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_mcd),          INTENT(INOUT) :: thisMCD
   TYPE(t_banddos),       INTENT(IN)    :: banddos
   TYPE(t_dimension),     INTENT(IN)    :: dimension
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_atoms),         INTENT(IN)    :: atoms

   ALLOCATE (thisMCD%ncore(atoms%ntype))
   ALLOCATE (thisMCD%e_mcd(atoms%ntype,input%jspins,dimension%nstd))
   IF (banddos%l_mcd) THEN
      thisMCD%emcd_lo = banddos%e_mcd_lo
      thisMCD%emcd_up = banddos%e_mcd_up
      ALLOCATE (thisMCD%m_mcd(dimension%nstd,(3+1)**2,3*atoms%ntype,2))
      ALLOCATE (thisMCD%mcd(3*atoms%ntype,dimension%nstd,dimension%neigd) )
      IF (.NOT.banddos%dos) WRITE (*,*) 'For mcd-spectra set banddos%dos=T!'
   ELSE
      ALLOCATE (thisMCD%m_mcd(1,1,1,1))
      ALLOCATE (thisMCD%mcd(1,1,1))
   ENDIF

   thisMCD%ncore = 0
   thisMCD%e_mcd = 0.0
   thisMCD%mcd = 0.0
   thisMCD%m_mcd = CMPLX(0.0,0.0)

END SUBROUTINE mcd_init1


END MODULE m_types_cdnval
