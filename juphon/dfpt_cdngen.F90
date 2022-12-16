!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_cdngen
#ifdef CPP_MPI
   USE mpi
#endif
CONTAINS

SUBROUTINE dfpt_cdngen(eig_id,dfpt_eig_id,fmpi,input,banddosdummy,vacuum,&
                  kpts,kqpts,atoms,sphhar,starsq,sym,gfinp,hub1inp,&
                  enpara,cell,noco,nococonv,vTot,resultsdummy, resultsdummy1,&
                  archiveType, xcpot,outDen,outDenIm,q_dfpt)

   use m_types_vacdos
   USE m_types
   USE m_constants
   USE m_juDFT
   USE m_dfpt_cdnval
   USE m_cdn_io
   USE m_wrtdop
   USE m_cdncore

   IMPLICIT NONE

   ! Type instance arguments
   TYPE(t_results),INTENT(INOUT)    :: resultsdummy, resultsdummy1
   TYPE(t_mpi),INTENT(IN)           :: fmpi
   TYPE(t_enpara),INTENT(INOUT)     :: enpara
   TYPE(t_banddos),INTENT(IN)       :: banddosdummy
   TYPE(t_input),INTENT(IN)         :: input
   TYPE(t_vacuum),INTENT(IN)        :: vacuum
   TYPE(t_noco),INTENT(IN)          :: noco
   TYPE(t_nococonv),INTENT(INOUT)   :: nococonv
   TYPE(t_sym),INTENT(IN)           :: sym
   TYPE(t_stars),INTENT(IN)         :: starsq
   TYPE(t_cell),INTENT(IN)          :: cell
   TYPE(t_kpts),INTENT(IN)          :: kpts, kqpts
   TYPE(t_sphhar),INTENT(IN)        :: sphhar
   TYPE(t_atoms),INTENT(IN)         :: atoms
   TYPE(t_potden),INTENT(IN)        :: vTot
   TYPE(t_gfinp),INTENT(IN)         :: gfinp
   TYPE(t_hub1inp),INTENT(IN)       :: hub1inp
   CLASS(t_xcpot),INTENT(IN)     :: xcpot
   TYPE(t_potden),INTENT(INOUT)     :: outDen, outDenIm

   !Scalar Arguments
   INTEGER, INTENT (IN)             :: eig_id, dfpt_eig_id, archiveType

   REAL, INTENT(IN) :: q_dfpt(3)

   ! Local type instances
   TYPE(t_regionCharges)          :: regCharges
   TYPE(t_dos),TARGET             :: dosdummy
   TYPE(t_vacdos),TARGET          :: vacdosdummy
   TYPE(t_moments)                :: moments
   TYPE(t_cdnvalJob)       :: cdnvalJob, cdnvalJob1

   !Local Scalars
   REAL                  :: fix, qtot, dummy,eFermiPrev
   INTEGER               :: jspin, ierr
   INTEGER               :: dim_idx
   INTEGER               :: n

   LOGICAL               :: l_error,Perform_metagga

   ! Initialization section
   CALL moments%init(fmpi,input,sphhar,atoms)
   !initalize data for DOS
   if (noco%l_noco) resultsdummy%eig(:,:,2)=resultsdummy%eig(:,:,1)
   if (noco%l_noco) resultsdummy1%eig(:,:,2)=resultsdummy1%eig(:,:,1)
   CALL dosdummy%init(input,atoms,kpts,banddosdummy,resultsdummy%eig)
   CALL vacdosdummy%init(input,atoms,kpts,banddosdummy,resultsdummy%eig)

   CALL outDen%init(starsq, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN,.TRUE.)
   CALL outDenIm%init(starsq, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN)

   DO jspin = 1,merge(1,input%jspins,noco%l_mperp)
      CALL cdnvalJob%init(fmpi,input,kpts,noco,resultsdummy,jspin)
      CALL cdnvalJob1%init(fmpi,input,kqpts,noco,resultsdummy1,jspin)
      CALL dfpt_cdnval(eig_id, dfpt_eig_id,fmpi,kpts,jspin,noco,nococonv,input,banddosdummy,cell,atoms,enpara,starsq,&
                       vacuum,sphhar,sym,vTot,cdnvalJob,outDen,dosdummy,vacdosdummy,&
                        hub1inp,kqpts, cdnvalJob1, resultsdummy, resultsdummy1, q_dfpt, outDenIm)
   END DO

   ! TODO: Implement this appropriately.
   !CALL timestart("cdngen: cdncore")
   !CALL cdncore(fmpi,input,vacuum,noco,nococonv,sym,&
   !             starsq,cell,sphhar,atoms,vTot,outDen,moments,results)
   !CALL timestop("cdngen: cdncore")

   ! These should already be broadcast.
!#ifdef CPP_MPI
!   CALL MPI_BCAST(nococonv%alph,atoms%ntype,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
!   CALL MPI_BCAST(nococonv%beta,atoms%ntype,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
!   CALL MPI_BCAST(nococonv%b_con,atoms%ntype*2,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
!   CALL MPI_BCAST(nococonv%qss,3,MPI_DOUBLE_PRECISION,0,fmpi%mpi_comm,ierr)
!#endif
   CALL outDen%distribute(fmpi%mpi_comm)
   CALL outDenIm%distribute(fmpi%mpi_comm)

END SUBROUTINE dfpt_cdngen

END MODULE m_dfpt_cdngen
