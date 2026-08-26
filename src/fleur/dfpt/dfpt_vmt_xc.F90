!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
      USE m_mt_tofrom_grid
      USE m_types_xcpot_inbuild
      USE m_types_xcpot_libxc
      USE m_types
      USE m_dfpt_gga_kernel
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
      TYPE(t_gradients)     :: gradRho, gradRho1, gradRho1Im
      TYPE(t_gradients)     :: gradDrivsigma, gradDrivsigma1, gradDrivsigma1Im
      TYPE(t_noco)          :: noco_loco
      REAL, ALLOCATABLE     :: ch(:,:),ch1(:,:),ch1Im(:,:),f_xc(:,:),v_xc1(:,:),v_xc1Im(:,:)
      REAL, ALLOCATABLE     :: drivsigma(:,:),driv2rho2(:,:),driv2rhosigma(:,:),driv2sigma2(:,:)
      REAL, ALLOCATABLE     :: drivsigmaT(:,:),drivsigma1(:,:),drivsigma1Im(:,:),drivsigmaMT(:,:,:)
      INTEGER               :: n,nsp,npoints
      INTEGER               :: iSpin, jSpin, fxcSpin

      !locals for fmpi
      integer :: ierr, nfxc, n_sigma
      integer:: n_start,n_stride
      LOGICAL :: lda_atom(atoms%ntype)

      noco_loco = noco
      noco_loco%l_unrestrictMT = .FALSE.

      nfxc = 2 * input%jspins - 1
      n_sigma = MERGE(1,3,input%jspins==1)

      SELECT TYPE(xcpot)
      TYPE IS(t_xcpot_inbuild)
         lda_atom=atoms%lda_atom
         IF (ANY(lda_atom)) THEN
            CALL judft_error("Using locally LDA not possible with DFPT.")
         ENDIF
      END SELECT

      nsp=atoms%nsp()

      CALL init_mt_grid(input%jspins,atoms,sphhar,xcpot%needs_grad(),sym)

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
         npoints = nsp*atoms%jri(n)
         ALLOCATE(ch(npoints,input%jspins),ch1(npoints,input%jspins),ch1Im(npoints,input%jspins))
         ALLOCATE(v_xc1(npoints,input%jspins),v_xc1Im(npoints,input%jspins))

         IF (xcpot%needs_grad()) THEN
            CALL xcpot%alloc_gradients(npoints,input%jspins,gradRho)
            CALL xcpot%alloc_gradients(npoints,input%jspins,gradRho1)
            CALL xcpot%alloc_gradients(npoints,input%jspins,gradRho1Im)
         END IF

         CALL mt_to_grid(xcpot%needs_grad(),input%jspins,atoms,sym,sphhar,.FALSE.,den%mt(:,0:,n,:),n,noco_loco,gradRho,ch)
         CALL mt_to_grid(xcpot%needs_grad(),input%jspins,atoms,sym,sphhar,.FALSE.,den1%mt(:,0:,n,:),n,noco_loco,gradRho1,ch1)
         CALL mt_to_grid(xcpot%needs_grad(),input%jspins,atoms,sym,sphhar,.FALSE.,den1im%mt(:,0:,n,:),n,noco_loco,gradRho1Im,ch1Im)

         IF (xcpot%needs_grad()) THEN
            ALLOCATE(drivsigma(n_sigma,npoints),driv2rho2(nfxc,npoints))
            ALLOCATE(driv2rhosigma(input%jspins*n_sigma,npoints),driv2sigma2(MERGE(1,6,input%jspins==1),npoints))
            SELECT TYPE(xcpot)
            TYPE IS (t_xcpot_libxc)
               CALL xcpot%get_fxc_gga(input%jspins,ch,gradRho%sigma,drivsigma,driv2rho2,driv2rhosigma,driv2sigma2)
            END SELECT

            ALLOCATE(drivsigma1(npoints,n_sigma),drivsigma1Im(npoints,n_sigma))
            CALL dfpt_gga_assemble(drivsigma,driv2rho2,driv2rhosigma,driv2sigma2,gradRho,gradRho1,ch1,v_xc1,drivsigma1)
            CALL dfpt_gga_assemble(drivsigma,driv2rho2,driv2rhosigma,driv2sigma2,gradRho,gradRho1Im,ch1Im,v_xc1Im,drivsigma1Im)

            !The MT expansion carries the exp(iqr) phase already, so no q-shift is
            !needed and the two channels are independent real round trips.
            ALLOCATE(drivsigmaT(npoints,n_sigma),drivsigmaMT(atoms%jri(n),0:sphhar%nlhd,n_sigma))
            drivsigmaT = TRANSPOSE(drivsigma)
            CALL mt_gradient_ftgrid(xcpot,atoms,sym,sphhar,noco_loco,n,drivsigmaT,drivsigmaMT,gradDrivsigma)
            CALL mt_gradient_ftgrid(xcpot,atoms,sym,sphhar,noco_loco,n,drivsigma1,drivsigmaMT,gradDrivsigma1)
            CALL mt_gradient_ftgrid(xcpot,atoms,sym,sphhar,noco_loco,n,drivsigma1Im,drivsigmaMT,gradDrivsigma1Im)

            CALL dfpt_gga_grdotgr(gradRho,gradRho1,gradDrivsigma,gradDrivsigma1,v_xc1)
            CALL dfpt_gga_grdotgr(gradRho,gradRho1Im,gradDrivsigma,gradDrivsigma1Im,v_xc1Im)

            DEALLOCATE(drivsigma,driv2rho2,driv2rhosigma,driv2sigma2)
            DEALLOCATE(drivsigmaT,drivsigma1,drivsigma1Im,drivsigmaMT)
         ELSE
            ALLOCATE(f_xc(npoints,nfxc))
            CALL xcpot%get_fxc(input%jspins, ch, f_xc)

            v_xc1 = 0.0
            v_xc1Im = 0.0
            DO iSpin = 1, input%jspins
                DO jSpin = 1, input%jspins
                    fxcSpin = iSpin + jSpin - 1
                    v_xc1(:, iSpin) = v_xc1(:, iSpin) + f_xc(:, fxcSpin) * ch1(:, jSpin)
                    v_xc1Im(:, iSpin) = v_xc1Im(:, iSpin) + f_xc(:, fxcSpin) * ch1Im(:, jSpin)
                END DO
            END DO
            DEALLOCATE(f_xc)
         END IF

         CALL mt_from_grid(atoms,sym,sphhar,n,input%jspins,v_xc1,vTot%mt(:,0:,n,:))
         CALL mt_from_grid(atoms,sym,sphhar,n,input%jspins,v_xc1Im,dfptvTotimag%mt(:,0:,n,:))

         DEALLOCATE (ch,ch1,ch1Im,v_xc1,v_xc1Im)
      ENDDO

      CALL finish_mt_grid()
#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,vTot%mt,SIZE(vTot%mt),MPI_DOUBLE_PRECISION,MPI_SUM,fmpi%mpi_comm,ierr)
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,dfptvTotimag%mt,SIZE(dfptvTotimag%mt),MPI_DOUBLE_PRECISION,MPI_SUM,fmpi%mpi_comm,ierr)
#endif
      !
      RETURN
  END SUBROUTINE dfpt_vmt_xc

   SUBROUTINE mt_gradient_ftgrid(xcpot,atoms,sym,sphhar,noco,n,f,f_mt,gradF)
      !! Gradient of a real field on the MT grid, obtained by projecting it back
      !! onto the lattice harmonics and differentiating there. f_mt is scratch.
      USE m_mt_tofrom_grid
      USE m_types
      IMPLICIT NONE
      CLASS(t_xcpot),INTENT(IN)       :: xcpot
      TYPE(t_atoms),INTENT(IN)        :: atoms
      TYPE(t_sym),INTENT(IN)          :: sym
      TYPE(t_sphhar),INTENT(IN)       :: sphhar
      TYPE(t_noco),INTENT(IN)         :: noco
      INTEGER,INTENT(IN)              :: n
      REAL,INTENT(IN)                 :: f(:,:)
      REAL,INTENT(INOUT)              :: f_mt(:,0:,:)
      TYPE(t_gradients),INTENT(INOUT) :: gradF

      INTEGER :: jr, nfields

      nfields = SIZE(f,2)
      !mt_from_grid accumulates
      f_mt = 0.0
      CALL mt_from_grid(atoms,sym,sphhar,n,nfields,f,f_mt)
      DO jr = 1, atoms%jri(n)
         f_mt(jr,:,:) = f_mt(jr,:,:)*atoms%rmsh(jr,n)**2
      END DO
      CALL xcpot%alloc_gradients(SIZE(f,1),nfields,gradF)
      CALL mt_to_grid(.TRUE.,nfields,atoms,sym,sphhar,.FALSE.,f_mt,n,noco,grad=gradF)
   END SUBROUTINE mt_gradient_ftgrid
END MODULE m_dfpt_vmt_xc
