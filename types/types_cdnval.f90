!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_cdnval
use m_types_jDOS
use m_types_orbcomp
use m_types_mcd
use m_types_slab

IMPLICIT NONE

PRIVATE

PUBLIC :: t_slab,t_orbcomp,t_mcd,t_jDOS
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



   TYPE t_eigVecCoeffs
      COMPLEX, ALLOCATABLE :: acof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: bcof(:,:,:,:)
      COMPLEX, ALLOCATABLE :: ccof(:,:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init => eigVecCoeffs_init
   END TYPE t_eigVecCoeffs



   TYPE t_moments

      REAL, ALLOCATABLE    :: chmom(:,:)
      REAL, ALLOCATABLE    :: clmom(:,:,:)
      COMPLEX, ALLOCATABLE :: qa21(:)

      REAL, ALLOCATABLE    :: stdn(:,:)
      REAL, ALLOCATABLE    :: svdn(:,:)

      REAL, ALLOCATABLE    :: rhoLRes(:,:,:,:,:)

      CONTAINS
         PROCEDURE,PASS :: init => moments_init
   END TYPE t_moments





   TYPE t_cdnvalJob
      LOGICAL              :: l_evp
      INTEGER, ALLOCATABLE :: k_list(:)
      INTEGER, ALLOCATABLE :: ev_list(:)
      INTEGER, ALLOCATABLE :: noccbd(:)    ! Attention, these are for all k-points and all states
      REAL,    ALLOCATABLE :: weights(:,:) !


      CONTAINS
         PROCEDURE,PASS :: init => cdnvalJob_init
         PROCEDURE      :: select_slice
         PROCEDURE      :: compact_ev_list
   END TYPE t_cdnvalJob

   TYPE t_gVacMap

      INTEGER, ALLOCATABLE    :: gvac1d(:)
      INTEGER, ALLOCATABLE    :: gvac2d(:)

      CONTAINS
         PROCEDURE,PASS :: init => gVacMap_init
   END TYPE t_gVacMap

PUBLIC t_orb, t_denCoeffs,  t_eigVecCoeffs
PUBLIC  t_moments,  t_cdnvalJob, t_gVacMap

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




SUBROUTINE eigVecCoeffs_init(thisEigVecCoeffs,input,atoms,jspin,noccbd,l_bothSpins)

   USE m_types_setup

   IMPLICIT NONE

   CLASS(t_eigVecCoeffs), INTENT(INOUT) :: thisEigVecCoeffs

   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_input),         INTENT(IN)    :: input

   INTEGER,               INTENT(IN)    :: jspin, noccbd
   LOGICAL,               INTENT(IN)    :: l_bothSpins

   IF(ALLOCATED(thisEigVecCoeffs%acof)) DEALLOCATE(thisEigVecCoeffs%acof)
   IF(ALLOCATED(thisEigVecCoeffs%bcof)) DEALLOCATE(thisEigVecCoeffs%bcof)
   IF(ALLOCATED(thisEigVecCoeffs%ccof)) DEALLOCATE(thisEigVecCoeffs%ccof)

   IF (l_bothSpins) THEN
      ALLOCATE (thisEigVecCoeffs%acof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat,input%jspins))
      ALLOCATE (thisEigVecCoeffs%bcof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat,input%jspins))
      ALLOCATE (thisEigVecCoeffs%ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat,input%jspins))
   ELSE
      ALLOCATE (thisEigVecCoeffs%acof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat,jspin:jspin))
      ALLOCATE (thisEigVecCoeffs%bcof(noccbd,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat,jspin:jspin))
      ALLOCATE (thisEigVecCoeffs%ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat,jspin:jspin))
   END IF

   thisEigVecCoeffs%acof = CMPLX(0.0,0.0)
   thisEigVecCoeffs%bcof = CMPLX(0.0,0.0)
   thisEigVecCoeffs%ccof = CMPLX(0.0,0.0)

END SUBROUTINE eigVecCoeffs_init


SUBROUTINE moments_init(thisMoments,mpi,input,sphhar,atoms)

   USE m_types_setup
   USE m_types_mpi

   IMPLICIT NONE

   CLASS(t_moments),      INTENT(INOUT) :: thisMoments
   TYPE(t_mpi),           INTENT(IN)    :: mpi
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_sphhar),        INTENT(IN)    :: sphhar
   TYPE(t_atoms),         INTENT(IN)    :: atoms

   ALLOCATE(thisMoments%chmom(atoms%ntype,input%jspins))
   ALLOCATE(thisMoments%clmom(3,atoms%ntype,input%jspins))
   ALLOCATE(thisMoments%qa21(atoms%ntype))

   ALLOCATE(thisMoments%stdn(atoms%ntype,input%jspins))
   ALLOCATE(thisMoments%svdn(atoms%ntype,input%jspins))

   thisMoments%chmom = 0.0
   thisMoments%clmom = 0.0
   thisMoments%qa21 = CMPLX(0.0,0.0)

   thisMoments%stdn = 0.0
   thisMoments%svdn = 0.0

   IF(mpi%irank.EQ.0.and.input%lresMax>0) THEN
      ALLOCATE(thisMoments%rhoLRes(atoms%jmtd,0:sphhar%nlhd,0:(input%lResMax*(input%lResMax+1))/2+input%lResMax,atoms%ntype,4))
      thisMoments%rhoLRes = 0.0
   END IF

END SUBROUTINE moments_init



SUBROUTINE cdnvalJob_init(thisCdnvalJob,mpi,input,kpts,noco,results,jspin)

   USE m_types_mpi
   USE m_types_setup
   USE m_types_kpts
   USE m_types_misc

   IMPLICIT NONE

   CLASS(t_cdnvalJob),             INTENT(OUT)   :: thisCdnvalJob
   TYPE(t_mpi),                    INTENT(IN)    :: mpi
   TYPE(t_input),                  INTENT(IN)    :: input
   TYPE(t_kpts),                   INTENT(IN)    :: kpts
   TYPE(t_noco),                   INTENT(IN)    :: noco
   TYPE(t_results),                INTENT(IN)    :: results


   INTEGER,                        INTENT(IN)    :: jspin

   INTEGER :: jsp, iBand, ikpt, nslibd, noccbd_l, noccbd, ikpt_i

   jsp = MERGE(1,jspin,noco%l_noco)

   thisCdnvalJob%l_evp=mpi%n_size>1
   thisCdnvalJob%k_list=mpi%k_list !includes allocate
   thisCdnvalJob%ev_list=mpi%ev_list

   thisCdnvalJob%weights = results%w_iks(:,:,jsp)*2.0/input%jspins

   ALLOCATE(thisCdnvalJob%noccbd(kpts%nkpt))
   thisCdnvalJob%noccbd = 0

   ! determine bands to be used for each k point, MPI process
   DO ikpt_i = 1,SIZE(thisCdnvalJob%k_list)
      ikpt=thisCdnvalJob%k_list(ikpt_i)
      !Max number of bands
      thisCdnvalJob%noccbd(ikpt)= COUNT(thiscdnvaljob%ev_list<=results%neig(ikpt,jsp))
   ENDDO
 END SUBROUTINE cdnvalJob_init

 SUBROUTINE select_slice(thiscdnvalJob,sliceplot,results,input,kpts,noco,jspin)
   USE m_types_setup
   USE m_types_misc
   USE m_types_kpts
   IMPLICIT NONE
   CLASS(t_cdnvalJob),INTENT(INOUT)  :: thisCdnvalJob
   TYPE(t_sliceplot), INTENT(IN)     :: sliceplot
   TYPE(t_results),    INTENT(IN)    :: results
   TYPE(t_input),INTENT(IN)          :: input
   TYPE(t_kpts),INTENT(IN)           :: kpts
   TYPE(t_noco),INTENT(IN)           :: noco
   INTEGER,INTENT(IN)                :: jspin

   INTEGER :: iband,iband_i,ikpt,ikpt_i,jsp
   jsp = MERGE(1,jspin,noco%l_noco)

   DO ikpt_i=1,SIZE(thiscdnvalJob%k_list)
      ikpt=thiscdnvalJob%k_list(ikpt_i)
      !--->    if slice, only certain bands are taken into account
      IF (sliceplot%slice.AND.input%pallst) thisCdnvalJob%weights(:,ikpt) = kpts%wtkpt(ikpt)*2.0/input%jspins
      IF (sliceplot%slice.AND.thisCdnvalJob%noccbd(ikpt).GT.0) THEN
         IF (sliceplot%kk.EQ.0) THEN
            DO iband_i=1,thisCdnvalJob%noccbd(ikpt)
               iband=thiscdnvaljob%ev_list(iband_i)
               IF (results%eig(iBand,ikpt,jsp).LT.sliceplot%e1s) thisCdnvalJob%weights(iband,ikpt)=0.0
               IF (results%eig(iBand,ikpt,jsp).GT.sliceplot%e2s) thisCdnvalJob%weights(iband,ikpt)=0.0
            END DO
         ELSE IF (sliceplot%kk.EQ.ikpt) THEN
            IF ((sliceplot%e1s.EQ.0.0) .AND. (sliceplot%e2s.EQ.0.0)) THEN
               DO iband_i=1,thisCdnvalJob%noccbd(ikpt)
                  iband=thiscdnvaljob%ev_list(iband_i)
                  IF (iBand.NE.sliceplot%nnne) thisCdnvalJob%weights(iband,ikpt)=0.0
               ENDDO
            ELSE
               DO iband_i=1,thisCdnvalJob%noccbd(ikpt)
                  iband=thiscdnvaljob%ev_list(iband_i)
                  IF (results%eig(iBand,ikpt,jsp).LT.sliceplot%e1s) thisCdnvalJob%weights(iband,ikpt)=0.0
                  IF (results%eig(iBand,ikpt,jsp).GT.sliceplot%e2s) thisCdnvalJob%weights(iband,ikpt)=0.0
               END DO
            END IF
         ELSE
            thisCdnvalJob%weights(:,ikpt)=0.0
         END IF
      END IF ! sliceplot%slice
   END DO
 END SUBROUTINE select_slice

 FUNCTION compact_ev_list(thiscdnvaljob,ikpt,l_empty)
   IMPLICIT NONE
   CLASS(t_cdnvalJob),INTENT(IN)  :: thisCdnvalJob
   INTEGER,INTENT(IN)             :: ikpt
   LOGICAL,INTENT(IN)             :: l_empty

   INTEGER,ALLOCATABLE :: compact_ev_list(:)
   INTEGER :: nk

   nk=thisCdnvalJob%k_list(ikpt)
   IF (l_empty) THEN
      compact_ev_list=thiscdnvalJob%ev_list(:thisCdnvalJob%noccbd(nk))
   ELSE
      compact_ev_list=PACK(thiscdnvalJob%ev_list(:thisCdnvalJob%noccbd(nk)),&
           thisCdnvalJob%weights(thiscdnvalJob%ev_list(:thisCdnvalJob%noccbd(nk)),nk)>1.e-8)
   END IF
 END FUNCTION compact_ev_list


SUBROUTINE gVacMap_init(thisGVacMap,sym,atoms,vacuum,stars,lapw,input,cell,kpts,enpara,vTot,ikpt,jspin)

   USE m_types_setup
   USE m_types_lapw
   USE m_types_enpara
   USE m_types_potden
   USE m_types_kpts
   USE m_nstm3

   IMPLICIT NONE

   CLASS(t_gVacMap),      INTENT(INOUT) :: thisGVacMap

   TYPE(t_sym),           INTENT(IN)    :: sym
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_vacuum),        INTENT(IN)    :: vacuum
   TYPE(t_stars),         INTENT(IN)    :: stars
   TYPE(t_lapw),          INTENT(IN)    :: lapw
   TYPE(t_input),         INTENT(IN)    :: input
   TYPE(t_cell),          INTENT(IN)    :: cell
   TYPE(t_kpts),          INTENT(IN)    :: kpts
   TYPE(t_enpara),        INTENT(IN)    :: enpara
   TYPE(t_potden),        INTENT(IN)    :: vTot

   INTEGER,               INTENT(IN)    :: ikpt
   INTEGER,               INTENT(IN)    :: jspin

   IF (ALLOCATED(thisGVacMap%gvac1d)) DEALLOCATE(thisGVacMap%gvac1d)
   IF (ALLOCATED(thisGVacMap%gvac2d)) DEALLOCATE(thisGVacMap%gvac2d)

   ALLOCATE(thisGVacMap%gvac1d(lapw%dim_nv2d()))
   ALLOCATE(thisGVacMap%gvac2d(lapw%dim_nv2d()))

   thisGVacMap%gvac1d = 0
   thisGVacMap%gvac2d = 0

   IF (vacuum%nstm.EQ.3.AND.input%film) THEN
      CALL nstm3(sym,atoms,vacuum,stars,lapw,ikpt,input,jspin,kpts,&
                 cell,enpara%evac0(1,jspin),vTot%vacz(:,:,jspin),thisGVacMap%gvac1d,thisGVacMap%gvac2d)
   END IF

END SUBROUTINE gVacMap_init

END MODULE m_types_cdnval
