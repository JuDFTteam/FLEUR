!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_vgen_xcpot

   USE m_juDFT

#ifdef CPP_MPI
   USE mpi
#endif

CONTAINS

    SUBROUTINE dfpt_vgen_xcpot(atoms, jspins, stars, fmpi, den, denRot, ylm, gWghts, fxcMT, fxcIR, iDatom, ngqdp, gqdp, ngdp, gdp, vTot)

        !---------------------------------------------------------------------
        ! FLAPW potential perturbation generator
        !---------------------------------------------------------------------
        ! Generates the exchange-correlation potential perturbation.
        !---------------------------------------------------------------------

        USE m_types
        USE m_constants

        IMPLICIT NONE

        TYPE(t_mpi), INTENT(IN)              :: fmpi
        !TYPE(t_noco), INTENT(IN)              :: noco
        TYPE(t_stars), INTENT(IN)              :: stars
        TYPE(t_atoms), INTENT(IN)              :: atoms
        TYPE(t_jpPotden), INTENT(IN)              :: den, denRot
        TYPE(t_jpPotden), INTENT(INOUT)           :: vTot
        complex,                       intent(in)  :: ylm(:, :)
        real,                          intent(in)  :: gWghts(:)
        real,                          intent(in)  :: fxcMT(:, :, :, :)
        complex,              intent(in)  :: fxcIR(:, :)
       INTEGER, INTENT(IN) :: jspins, iDatom, ngqdp, ngdp
        INTEGER, INTENT(IN) :: gqdp(:, :), gdp(:, :)

        INTEGER :: i, iType

#ifdef CPP_MPI
        integer:: ierr
#endif

        IF (fmpi%irank == 0) THEN
            CALL timestart("Vxc perturbation in interstitial")
            CALL dfpt_vis_xc(stars, fxcIR, jspins, iDatom, ngqdp, gqdp, ngdp, gdp, den, vTot)
            CALL timestop("Vxc perturbation in interstitial")

            CALL timestart("Vxc perturbation in MT")
        END IF

        CALL dfpt_vmt_xc(fmpi, atoms, ylm, gWghts, fxcMT, jspins, iDatom, den, vTot)

        IF (fmpi%irank == 0) THEN
            CALL timestop("Vxc perturbation in MT")
        END IF ! fmpi%irank == 0

    END SUBROUTINE dfpt_vgen_xcpot

    SUBROUTINE dfpt_vis_xc(stars, fxcIR, jspins, iDatom, ngqdp, gqdp, ngdp, gdp, den, vTot)

       USE m_types

       IMPLICIT NONE

       TYPE(t_stars),INTENT(IN)       :: stars
       TYPE(t_jpPotden),INTENT(IN)    :: den
       TYPE(t_jpPotden),INTENT(INOUT) :: vTot
       complex,              intent(in)  :: fxcIR(:, :)
      INTEGER, INTENT(IN) :: jspins, iDatom, ngqdp, ngdp
       INTEGER, INTENT(IN) :: gqdp(:, :), gdp(:, :)

       REAL, ALLOCATABLE :: rho(:,:)
       integer,        allocatable                 :: igfftg(:)
       integer,        allocatable                 :: igfftgq(:)
       INTEGER           :: jspin, i, js, iDir, iG
       integer :: nfftx, nffty, nfftz, nfftxy, GxFFT, GyFFT, GzFFT

       ! Construct mapping array for mapping a G-vector to its respective mesh point on the FFT mesh.
       ! This constitutes the parallel to init_pw_grid().
       CALL timestart("dfpt_init_pw_grid")
       nfftx  = 3 * stars%mx1
       nffty  = 3 * stars%mx2
       nfftz  = 3 * stars%mx3
       nfftxy = 9 * stars%mx1 * stars%mx2

       ! Mapping array for fxc
       allocate(igfftg(ngdp))
       igfftg(:) = 0
       DO iG = 1, ngdp
         GxFFT = gdp(1, iG)
         GyFFT = gdp(2, iG)
         GzFFT = gdp(3, iG)
         if (GxFFT < 0) GxFFT = GxFFT + nfftx
         if (GyFFT < 0) GyFFT = GyFFT + nffty
         if (GzFFT < 0) GzFFT = GzFFT + nfftz
         igfftg(iG) = GxFFT + GyFFT * nfftx + GzFFT * nfftxy
       END DO

       ! Mapping array for rho1/grad_Rho
       allocate(igfftgq(ngqdp))
       igfftgq = 0
       DO iG = 1, ngqdp
         GxFFT = gqdp(1, iG)
         GyFFT = gqdp(2, iG)
         GzFFT = gqdp(3, iG)
         if (GxFFT < 0) GxFFT = GxFFT + nfftx
         if (GyFFT < 0) GyFFT = GyFFT + nffty
         if (GzFFT < 0) GzFFT = GzFFT + nfftz
         igfftgq(iG) = GxFFT + GyFFT * nfftx + GzFFT * nfftxy
       END DO
       CALL timestop("dfpt_init_pw_grid")

       ! Convolute rho1IR with functional derivative of kernel to get Vxc within the IR.
       DO iDir = 1, 3
         CALL dfpt_pw_tofrom_grid(stars, jspins, ngdp, ngqdp, iDir, den%pw(:, :, iDatom, iDir), fxcIR, igfftg, igfftgq, vTot%pw(:, :, iDatom, iDir))
       END DO

   END SUBROUTINE dfpt_vis_xc

   SUBROUTINE dfpt_vmt_xc(fmpi,atoms,ylm,gWghts,fxcMT,jspins, iDatom, den, vTot)
#include"cpp_double.h"

      USE m_types

      IMPLICIT NONE

      TYPE(t_mpi),INTENT(IN)         :: fmpi
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_jpPotden),INTENT(IN)      :: den
      TYPE(t_jpPotden),INTENT(INOUT)   :: vTot

      complex,                       intent(in)  :: ylm(:, :)
      real,                          intent(in)  :: gWghts(:)
      real,                          intent(in)  :: fxcMT(:, :, :, :)
      INTEGER, INTENT(IN) :: jspins, iDatom
      !     ..
      !     .. Local Scalars ..
      TYPE(t_gradients)     :: grad
      TYPE(t_jpPotden)        :: vTot_tmp
      INTEGER               :: n,nsp,nt,jr, loc_n
      INTEGER               :: i, j, idx, cnt

      integer :: ierr, iAtom
      integer:: n_start,n_stride

      nsp=atoms%nsp()

      ! We can skip this, because we provide it as input.
      !CALL init_mt_grid(input%jspins,atoms,sphhar,xcpot%needs_grad(),sym)

#ifdef CPP_MPI
      n_start=fmpi%irank+1
      n_stride=fmpi%isize
      IF (fmpi%irank>0) THEN
        vTot%mt = cmplx(0.0, 0.0)
      ENDIF
#else
      n_start=1
      n_stride=1
#endif

      DO iAtom = n_start, atoms%nat, n_stride
         CALL dfpt_mt_tofrom_grid(atoms, jspins, atoms%itype(iAtom), den%mt(:, :, iAtom, :, iDatom, :), gWghts, fxcMT(:, :, iAtom, :), ylm, vTot%mt(:, :, iAtom, :, iDatom, :))
      END DO

#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE,vTot%mt,SIZE(vTot%mt),CPP_MPI_COMPLEX,MPI_SUM,fmpi%mpi_comm,ierr)
#endif
      RETURN
  END SUBROUTINE dfpt_vmt_xc

  subroutine dfpt_pw_tofrom_grid(stars, jspins, ngdp, ngqdp, iDir, rho1IR, fxcIR, igfftg, igfftgq, vTot1IR)

    USE m_fft_interface
    use m_types

    IMPLICIT NONE

    type(t_stars),        intent(in)  :: stars
    integer,              intent(in)  :: ngdp
    integer,              intent(in)  :: ngqdp
    integer,              intent(in)  :: iDir
    complex,              intent(in)  :: rho1IR(:, :)
    complex,              intent(in)  :: fxcIR(:, :)
    integer,              intent(in)  :: igfftg(:)
    integer,              intent(in)  :: igfftgq(:)
    complex,              intent(inout) :: vTot1IR(:, :)
    INTEGER, INTENT(IN) :: jspins

    integer                           :: ifftd
    integer                           :: iG
    real                              :: scaling
    INTEGER :: length_zfft(3), iSpin, jSpin, fxcSpin

    COMPLEX, ALLOCATABLE :: zfftfxc1(:), zfftfxc2(:), zfftfxc3(:), zfftfxc(:, :)
    COMPLEX, ALLOCATABLE :: zfftrho1(:), zfftvxc1(:)

    ifftd = 27 * stars%mx1 * stars%mx2 * stars%mx3

    ALLOCATE(zfftfxc1(0:ifftd - 1), zfftfxc2(0:ifftd - 1), zfftfxc3(0:ifftd - 1))
    ALLOCATE(zfftfxc(0:ifftd - 1, 3))
    ALLOCATE(zfftrho1(0:ifftd - 1), zfftvxc1(0:ifftd - 1))

    length_zfft(1) = 3*stars%mx1
    length_zfft(2) = 3*stars%mx2
    length_zfft(3) = 3*stars%mx3
    scaling = 1. / ifftd

    !fxcSpin = 2 * jspins - 1 ! Number of fxc spin-spin indices

    CALL timestart("fxc to grid")
    zfftfxc1 = cmplx(0.0, 0.0)
    zfftfxc2 = cmplx(0.0, 0.0)
    zfftfxc3 = cmplx(0.0, 0.0)
    zfftfxc  = cmplx(0.0, 0.0)

    do iG = 1, ngdp
      zfftfxc1(igfftg(iG)) = fxcIR(iG, 1)
      IF (jspins.EQ.2) zfftfxc2(igfftg(iG)) = fxcIR(iG, 2)
      IF (jspins.EQ.2) zfftfxc3(igfftg(iG)) = fxcIR(iG, 3)
    end do

    ! Complex FFT of fxc into real space; c.f. fft3d
    call fft_interface(3, length_zfft, zfftfxc1, .FALSE., igfftg)
    IF (jspins.EQ.2) call fft_interface(3, length_zfft, zfftfxc2, .FALSE., igfftg)
    IF (jspins.EQ.2) call fft_interface(3, length_zfft, zfftfxc3, .FALSE., igfftg)

    zfftfxc(:, 1)  = zfftfxc1
    zfftfxc(:, 2)  = zfftfxc2
    zfftfxc(:, 3)  = zfftfxc3

    CALL timestop("fxc to grid")

    DO iSpin = 1, jspins
        zfftrho1 = cmplx(0.0, 0.0)
        do iG = 1, ngqdp
          zfftrho1(igfftgq(iG)) = rho1IR(iG, iSpin)
        end do

        ! Complex FFT of rho1 into real space; c.f. fft3d
        call fft_interface(3, length_zfft, zfftrho1, .FALSE., igfftgq)

        DO jSpin = 1, jspins
            zfftvxc1 = cmplx(0.0, 0.0)
            fxcSpin = iSpin + jSpin - 1
            zfftvxc1 = zfftfxc(:, fxcSpin) * zfftrho1

            ! Reintroduce q=0 --> iqpt
            !IF (ANY(AIMAG(zfftVxc1).GT.1e-7)) THEN
            !    write(*, *) 'Warning: Real space pw Vxc1 has imaginary components!'
            !    write(*, *) 'Maximum absolute value:', MAXVAL(ABS(AIMAG(zfft3)))
            !END IF

            ! Complex FFT of Vxc1 into reciprocal space
            CALL fft_interface(3, length_zfft, zfftVxc1, .TRUE., igfftgq)

            ! We have to care for the artifacts of this FFT
            ! Map convoluted quantity from FFT mesh to plane-wave expansion coefficient representation.
            do iG = 1, ngqdp
                vTot1IR(iG, jSpin) = vTot1IR(iG, jSpin) + zfftVxc1(igfftgq(iG)) * scaling
            end do

        END DO
    END DO
    DEALLOCATE(zfftfxc1, zfftfxc2, zfftfxc3, zfftfxc, zfftrho1, zfftvxc1)

  end subroutine dfpt_pw_tofrom_grid

  subroutine dfpt_mt_tofrom_grid(atoms, jspins, iType, rho1MT, gWghts, fxcMT, ylm, vTot1MT)
      
    use m_types

    implicit none

    ! Type parameter
    type(t_atoms),                 intent(in)  :: atoms

    INTEGER, INTENT(IN) :: iType
    ! Array parameter
    complex,                       intent(in)  :: rho1MT(:, :, :, :)
    complex,                       intent(in)  :: ylm(:, :)
    real,                          intent(in)  :: gWghts(:)
    real,                          intent(in)  :: fxcMT(:, :, :)
    complex,                       intent(inout) :: vTot1MT(:, :, :, :)
    INTEGER, INTENT(IN) :: jspins

    ! Local scalar variables
    integer                                    :: idir, iSpin, jSpin, fxcSpin
    integer                                    :: ieqat
    integer                                    :: igmesh ! Loop variable over sampling points of spherical Gauss mesh
    integer                                    :: irmesh ! Loop variable over radial MT mesh
    integer                                    :: oqn_l
    integer                                    :: ll1
    integer                                    :: mqn_m
    integer                                    :: lm
    complex                                    :: vxcMT1KernAdd

    ! Local allocatable variables
    complex,           allocatable             :: rhoMT1Gpts(:, :) !grRhoMT on Gauss mesh
    complex,           allocatable             :: vxcMT1KernGPts(:, :)

    allocate( rhoMT1Gpts(atoms%nsp(), atoms%jmtd), vxcMT1KernGPts(atoms%nsp(), atoms%jmtd) )

        DO iDir = 1, 3
            DO iSpin = 1, jspins
                rhoMT1Gpts(:, :) = cmplx(0., 0.)
                DO oqn_l = 0, atoms%lmax(itype)
                    ll1 = oqn_l * ( oqn_l + 1 ) + 1
                    DO mqn_m = -oqn_l, oqn_l
                        lm = ll1 + mqn_m
                        DO irmesh = 1, atoms%jri(itype)
                            DO igmesh = 1, atoms%nsp()
                                rhoMT1Gpts(igmesh, irmesh) = rhoMT1Gpts(igmesh, irmesh) &
                                                         & + rho1MT(irmesh, lm, idir, iSpin) * ylm(igmesh, lm)
                            END DO ! igmesh
                        END DO ! irmesh
                    END DO ! mqn_m
                END DO ! oqn_l

          ! Reintroduce q=0 --> iqpt
          !IF (ANY(AIMAG(rhoMT1Gpts).GT.1e-7)) THEN
        !    write(*, *) 'Warning: Real space mt Vxc1 has imaginary components!'
        !    write(*, *) 'Maximum absolute value:', MAXVAL(ABS(AIMAG(rhoMT1Gpts)))
         ! END IF
                DO jSpin = 1, jspins
                    vxcMT1KernGPts(:, :) = cmplx(0., 0.)
                    fxcSpin = iSpin + jSpin - 1
                    ! On the spherical Gauss mesh the integral reduces to a weighted (gWghts) sum (over all sampling points on the Gauss mesh)
                    ! of the MT exchange-correlation kernel and either the density's gradient or the first variation of the gradient.
                    DO irmesh = 1, atoms%jri(itype)
                        DO igmesh = 1, atoms%nsp()
                            vxcMT1KernGPts(igmesh, irmesh) = vxcMT1KernGPts(igmesh, irmesh) &
                                                         & + rhoMT1Gpts(igmesh, irmesh) &
                                                         & * fxcMT(igmesh, irmesh, fxcSpin) * gWghts(igmesh)
                        END DO
                    END DO

                    DO oqn_l = 0, atoms%lmax(itype)
                        ll1 = oqn_l * ( oqn_l + 1 ) + 1
                        DO mqn_m = -oqn_l, oqn_l
                            lm = ll1 + mqn_m
                            DO irmesh = 1, atoms%jri(itype)
                                ! Back-transformation of the MT coefficients. Now they are expansion coefficients of the MT grid.
                                vxcMT1KernAdd = dot_product( ylm(:atoms%nsp(), lm), vxcMT1KernGPts(:atoms%nsp(), irmesh) )
                                ! Add this contribution to MT exchange-correlation contribution to the potential
                                vTot1MT(irmesh, lm, idir, jSpin) = vTot1MT(irmesh, lm, idir, jSpin) + vxcMT1KernAdd
                            END DO ! iR
                        END DO ! m
                    END DO ! l
                END DO ! jSpin
            END DO ! iSpin
        END DO ! iDir

    END SUBROUTINE dfpt_mt_tofrom_grid

END MODULE m_dfpt_vgen_xcpot
