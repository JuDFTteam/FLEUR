!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vmt_xc
  USE m_judft
   !.....------------------------------------------------------------------
   !     Calculate the GGA xc-potential in the MT-spheres
   !.....------------------------------------------------------------------
   !     instead of vmtxcor.f: the different exchange-correlation
   !     potentials defined through the key icorr are called through
   !     the driver subroutine vxcallg.f, subroutines vectorized
   !     ** r.pentcheva 22.01.96
   !     *********************************************************
   !     angular mesh calculated on speacial gauss-legendre points
   !     in order to use orthogonality of lattice harmonics and
   !     avoid a least square fit
   !     ** r.pentcheva 04.03.96
   !     *********************************************************
   !     MPI and OpenMP parallelization
   !             U.Alekseeva, February 2017
   !     *********************************************************

CONTAINS
   SUBROUTINE vmt_xc(DIMENSION,mpi,sphhar,atoms,&
                     den,xcpot,input,sym, obsolete,vxc,vx,exc)

#include"cpp_double.h"
      use m_libxc_postprocess_gga
      USE m_mt_tofrom_grid
      USE m_types_xcpot_inbuild
      USE m_types
      IMPLICIT NONE

      CLASS(t_xcpot),INTENT(IN)      :: xcpot
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_mpi),INTENT(IN)         :: mpi
      TYPE(t_obsolete),INTENT(IN)    :: obsolete
      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_sym),INTENT(IN)         :: sym
      TYPE(t_sphhar),INTENT(IN)      :: sphhar
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_potden),INTENT(IN)      :: den
      TYPE(t_potden),INTENT(INOUT)   :: vxc,vx,exc
#ifdef CPP_MPI
      include "mpif.h"
#endif
      !     ..
      !     .. Local Scalars ..
      TYPE(t_gradients)     :: grad
      TYPE(t_xcpot_inbuild) :: xcpot_tmp
      REAL, ALLOCATABLE     :: ch(:,:)
      INTEGER               :: n,nsp,nt,jr
      REAL                  :: divi

      !     ..

      !locals for mpi
      integer :: ierr
      integer:: n_start,n_stride
      REAL:: v_x((atoms%lmaxd+1+MOD(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)*atoms%jmtd,input%jspins)
      REAL:: v_xc((atoms%lmaxd+1+MOD(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)*atoms%jmtd,input%jspins)
      REAL:: e_xc((atoms%lmaxd+1+MOD(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)*atoms%jmtd,1)
      REAL,ALLOCATABLE:: xcl(:,:)
      LOGICAL :: lda_atom(atoms%ntype),l_libxc
      !.....------------------------------------------------------------------
      lda_atom=.FALSE.; l_libxc=.FALSE.
      SELECT TYPE(xcpot)
      TYPE IS(t_xcpot_inbuild)
         lda_atom=xcpot%lda_atom
         IF (ANY(lda_atom)) THEN
            IF((.NOT.xcpot%is_name("pw91"))) &
               CALL judft_warn("Using locally LDA only possible with pw91 functional")
            CALL xcpot_tmp%init("l91",.FALSE.,atoms%ntype)
            ALLOCATE(xcl(SIZE(v_xc,1),SIZE(v_xc,2)))
         ENDIF
      CLASS DEFAULT
         l_libxc=.true. !libxc!!
      END SELECT

      nsp=atoms%nsp()
      ALLOCATE(ch(nsp*atoms%jmtd,input%jspins))
      IF (xcpot%is_gga()) CALL xcpot%alloc_gradients(SIZE(ch,1),input%jspins,grad)

      CALL init_mt_grid(input%jspins,atoms,sphhar,xcpot,sym)

#ifdef CPP_MPI
      n_start=mpi%irank+1
      n_stride=mpi%isize
      IF (mpi%irank>0) THEN
         vxc%mt=0.0
         vx%mt=0.0
         exc%mt=0.0
      ENDIF
#else
      n_start=1
      n_stride=1
#endif

      DO n = n_start,atoms%ntype,n_stride
         CALL mt_to_grid(xcpot, input%jspins, atoms,sphhar,den%mt(:,0:,n,:),n,grad,ch)
         !
         !         calculate the ex.-cor. potential
         CALL xcpot%get_vxc(input%jspins,ch(:nsp*atoms%jri(n),:),v_xc(:nsp*atoms%jri(n),:),v_x(:nsp*atoms%jri(n),:),grad)
         IF (lda_atom(n)) THEN
            ! Use local part of pw91 for this atom
            CALL xcpot_tmp%get_vxc(input%jspins,ch(:nsp*atoms%jri(n),:),xcl(:nsp*atoms%jri(n),:),v_x(:nsp*atoms%jri(n),:),grad)
            !Mix the potentials
            divi = 1.0 / (atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(1,n))
            nt=0
            DO jr=1,atoms%jri(n)
               v_xc(nt+1:nt+nsp,:) = ( xcl(nt+1:nt+nsp,:) * ( atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(jr,n) ) +&
                                      v_xc(nt+1:nt+nsp,:) * ( atoms%rmsh(jr,n) - atoms%rmsh(1,n) ) ) * divi
               nt=nt+nsp
            ENDDO
         ENDIF

         !Add postprocessing for libxc
         IF (l_libxc.AND.xcpot%is_gga()) CALL libxc_postprocess_gga_mt(xcpot,atoms,sphhar,n,v_xc,grad)

         CALL mt_from_grid(atoms,sphhar,n,input%jspins,v_xc,vxc%mt(:,0:,n,:))
         CALL mt_from_grid(atoms,sphhar,n,input%jspins,v_x,vx%mt(:,0:,n,:))

         IF (ALLOCATED(exc%mt)) THEN
            !
            !           calculate the ex.-cor energy density
            !
            CALL xcpot%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),e_xc(:nsp*atoms%jri(n),1),grad)
            IF (lda_atom(n)) THEN
               ! Use local part of pw91 for this atom
               CALL xcpot_tmp%get_exc(input%jspins,ch(:nsp*atoms%jri(n),:),xcl(:nsp*atoms%jri(n),1),grad)
               !Mix the potentials
               nt=0
               DO jr=1,atoms%jri(n)
                  e_xc(nt+1:nt+nsp,1) = ( xcl(nt+1:nt+nsp,1) * ( atoms%rmsh(atoms%jri(n),n) - atoms%rmsh(jr,n) ) +&
                                         e_xc(nt+1:nt+nsp,1) * ( atoms%rmsh(jr,n) - atoms%rmsh(1,n) ) ) * divi
                  nt=nt+nsp
               END DO
            ENDIF
            CALL mt_from_grid(atoms,sphhar,n,1,e_xc,exc%mt(:,0:,n,:))
         ENDIF
      ENDDO

      CALL finish_mt_grid()
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,vx%mt,SIZE(vx%mt),CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,vxc%mt,SIZE(vxc%mt),CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,exc%mt,SIZE(exc%mt),CPP_MPI_REAL,MPI_SUM,mpi%mpi_comm,ierr)
#endif
      !
      RETURN
   END SUBROUTINE vmt_xc
END MODULE m_vmt_xc
