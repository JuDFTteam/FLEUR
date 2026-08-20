!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vmt_xc
#ifdef CPP_MPI
use mpi
#endif
USE m_judft

CONTAINS
   SUBROUTINE dfpt_vmt_xc(fmpi,sphhar,atoms,den,den1,den1im,xcpot,input,sym,noco,vTot,dfptvTotimag)
      use m_libxc_postprocess_gga
      USE m_mt_tofrom_grid
      USE m_types_xcpot_inbuild
      USE m_types
      USE m_metagga
      IMPLICIT NONE

      CLASS(t_xcpot),INTENT(IN)      :: xcpot
      TYPE(t_mpi),INTENT(IN)         :: fmpi
      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_sym),INTENT(IN)         :: sym
      TYPE(t_sphhar),INTENT(IN)      :: sphhar
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_potden),INTENT(IN)      :: den, den1, den1im
      TYPE(t_noco), INTENT(IN)       :: noco
      TYPE(t_potden),INTENT(INOUT)   :: vTot, dfptvTotimag
      !     ..
      !     .. Local Scalars ..
      TYPE(t_gradients)     :: grad
      TYPE(t_xcpot_inbuild) :: xcpot_tmp
      TYPE(t_potden)        :: vTot_tmp
      TYPE(t_noco)          :: noco_loco
      REAL, ALLOCATABLE     :: ch(:,:),chre(:,:),chim(:,:),f_xc(:,:),v_xc1re(:,:),v_xc1im(:,:)
      INTEGER               :: n,nsp,nt,jr
      INTEGER               :: i, j, idx, cnt, iSpin, jSpin, fxcSpin
      REAL                  :: divi

      !locals for fmpi
      integer :: ierr, nfxc
      integer:: n_start,n_stride
      LOGICAL :: lda_atom(atoms%ntype),l_libxc, perform_MetaGGA

      noco_loco = noco
      noco_loco%l_unrestrictMT = .FALSE.

      nfxc = 2 * input%jspins - 1

      l_libxc=.FALSE.
      SELECT TYPE(xcpot)
      TYPE IS(t_xcpot_inbuild)
         lda_atom=atoms%lda_atom
         IF (ANY(lda_atom)) THEN
            CALL judft_error("Using locally LDA not possible with DFPT.")
         ENDIF
      CLASS DEFAULT
         l_libxc=.true.
      END SELECT

      nsp=atoms%nsp()

      CALL init_mt_grid(input%jspins,atoms,sphhar,.FALSE.,sym)

#ifdef CPP_MPI
      n_start=fmpi%irank+1
      n_stride=fmpi%isize
      IF (fmpi%irank>0) THEN
         vTot%mt=0.0
         dfptvTotimag%mt=0.0
      ENDIF
#else
      n_start=1
      n_stride=1
#endif
      DO n = n_start,atoms%ntype,n_stride
         ALLOCATE(ch(nsp*atoms%jri(n),input%jspins),f_xc(nsp*atoms%jri(n),nfxc))
         ALLOCATE(chre(nsp*atoms%jri(n),input%jspins),chim(nsp*atoms%jri(n),input%jspins))
         ALLOCATE(v_xc1re(nsp*atoms%jri(n),input%jspins),v_xc1im(nsp*atoms%jri(n),input%jspins))

         CALL mt_to_grid(.FALSE., input%jspins, atoms,sym,sphhar,.FALSE.,den%mt(:,0:,n,:),n,noco_loco,grad,ch)
         CALL mt_to_grid(.FALSE., input%jspins, atoms,sym,sphhar,.FALSE.,den1%mt(:,0:,n,:),n,noco_loco,grad,chre)
         CALL mt_to_grid(.FALSE., input%jspins, atoms,sym,sphhar,.FALSE.,den1im%mt(:,0:,n,:),n,noco_loco,grad,chim)

#ifdef CPP_LIBXC
        CALL xcpot%get_fxc(input%jspins, ch, f_xc)
#else
        CALL judft_error("You compiled Fleur without libxc but want to use DFPT. Please fix that.")
        !CALL xcpot%get_vxc(input%jspins,ch,v_xc,v_x,grad)
        !TODO: Maybe place the old way with x-Alpha here for fun.
#endif

        v_xc1re = 0.0
        v_xc1im = 0.0
        DO iSpin = 1, input%jspins
            DO jSpin = 1, input%jspins
                fxcSpin = iSpin + jSpin - 1
                v_xc1re(:, iSpin) = v_xc1re(:, iSpin) + f_xc(:, fxcSpin) * chre(:, jSpin)
                v_xc1im(:, iSpin) = v_xc1im(:, iSpin) + f_xc(:, fxcSpin) * chim(:, jSpin)
            END DO
        END DO

         CALL mt_from_grid(atoms,sym,sphhar,n,input%jspins,v_xc1re,vTot%mt(:,0:,n,:))
         CALL mt_from_grid(atoms,sym,sphhar,n,input%jspins,v_xc1im,dfptvTotimag%mt(:,0:,n,:))

         DEALLOCATE (ch,chre,chim,f_xc,v_xc1re,v_xc1im)
      ENDDO

      CALL finish_mt_grid()
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,vTot%mt,SIZE(vTot%mt),MPI_DOUBLE_PRECISION,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,dfptvTotimag%mt,SIZE(dfptvTotimag%mt),MPI_DOUBLE_PRECISION,MPI_SUM,fmpi%mpi_comm,ierr)
#endif
      !
      RETURN
  END SUBROUTINE dfpt_vmt_xc
END MODULE m_dfpt_vmt_xc
