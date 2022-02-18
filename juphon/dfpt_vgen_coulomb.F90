!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dfpt_vgen_coulomb

  USE m_juDFT

#ifdef CPP_MPI
  USE mpi
#endif

CONTAINS

    SUBROUTINE dfpt_vgen_coulomb(iSpin, fmpi, stars, cell, atoms, rho0, den, vTot, iDatom, iDtype, gqdp, ngqdp, gdp, ngdp, qvec, l_ext)
        !---------------------------------------------------------------------
        ! FLAPW potential perturbation generator
        !---------------------------------------------------------------------
        ! Generates the Coulomb potential perturbation
        !
        ! It takes a spin variable to indicate in which spin-channel
        ! the charge resides.
        !---------------------------------------------------------------------

        USE m_constants
        USE m_types
        USE m_intnv
        USE m_checkdopall
        USE m_convol
        USE m_cfft

        IMPLICIT NONE

        INTEGER,          INTENT(IN)    :: iSpin
        TYPE(t_mpi),      INTENT(IN)    :: fmpi
        TYPE(t_stars),    INTENT(IN)    :: stars
        TYPE(t_cell),     INTENT(IN)    :: cell
        TYPE(t_atoms),    INTENT(IN)    :: atoms
        TYPE(t_jpPotden), INTENT(IN)    :: rho0
        TYPE(t_jpPotden), INTENT(IN)    :: den
        TYPE(t_jpPotden), INTENT(INOUT) :: vTot
        INTEGER,          INTENT(IN)    :: iDatom, iDtype
        INTEGER,          INTENT(IN)    :: gqdp(:, :)
        INTEGER,          INTENT(IN)    :: ngqdp
        INTEGER,          INTENT(IN)    :: gdp(:, :)
        INTEGER,          INTENT(IN)    :: ngdp
        REAL,             INTENT(IN)    :: qvec(3)
        LOGICAL,          INTENT(IN)    :: l_ext

        INTEGER :: i, j, k, l, iDir, iG
        REAL    :: gqext(3)
        REAL    :: gqnorm

        COMPLEX, ALLOCATABLE :: psq(:, :)

#ifdef CPP_MPI
        INTEGER :: ierr
#endif

        ALLOCATE ( psq(3, ngqdp)  )
        !vCoul%iter = den%iter

        CALL timestart( "dfpt_psqpw" )
        CALL dfpt_psqpw( fmpi, atoms, cell, &
                       & den%pw(:, iSpin, iDatom, :), den%mt(:,:,:,iSpin,iDatom,:), &
                       & rho0%mt(:, :, :, iSpin, 1, 1), rho0%pw(:, iSpin, 1, 1), psq, l_ext, iDatom, iDtype, ngqdp, gqdp, ngdp, gdp, qvec)
        CALL timestop( "dfpt_psqpw" )

        IF (fmpi%irank == 0) THEN
            ! Interstitial potential perturbation
            CALL timestart( "interstitial" )
            WRITE(oUnit, fmt=8010 )
8010  FORMAT (/,5x,'coulomb potential in the interstitial region:')
            DO iG = 1, ngqdp
                gqext = matmul(cell%bmat, gqdp(:, iG) + qvec)
                gqnorm = norm2(gqext)
                IF (gqnorm .LE. 1e-12) CYCLE
                DO iDir = 1, 3
                    vTot%pw(iG, iSpin, iDatom, idir) = fpi_const * psq(idir, iG) / gqnorm**2
                END DO
            END DO
            CALL timestop("interstitial")
        END IF ! fmpi%irank == 0

        ! Muffin Tin potential perturbation
        CALL timestart("MT-spheres")

#ifdef CPP_MPI
        CALL MPI_BARRIER(fmpi%mpi_comm,ierr) !should be totally useless, but needed anyway????
        CALL MPI_BCAST( vTot%pw, size(vTot%pw), MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr )
        CALL MPI_BARRIER(fmpi%mpi_comm,ierr) !should be totally useless, but ...
#endif

        DO iDir = 1, 3
            CALL dfpt_vmts(fmpi, atoms, cell, &
                         & vTot%pw(:, iSpin, iDatom, :), den%mt(:,:,:,iSpin,iDatom,:), vTot%mt(:,:,:,iSpin,iDatom,iDir), &
                         & iDatom, iDtype, iDir, gqdp, ngqdp, qvec, l_ext)
        END DO

        CALL timestop("MT-spheres")

    END SUBROUTINE dfpt_vgen_coulomb

    SUBROUTINE dfpt_vmts( fmpi, atoms, cell, vpw, rho, vr, iDatom, iDtype, iDir, gqdp, ngqdp, qvec, l_ext)

#include"cpp_double.h"

        USE m_constants
        USE m_types
        USE m_intgr, only : intgr2LinIntp
        USE m_sphbes
!        !$ use omp_lib

        IMPLICIT NONE

        TYPE(t_mpi),   INTENT(IN)  :: fmpi
        TYPE(t_atoms), INTENT(IN)  :: atoms
        TYPE(t_cell),  INTENT(IN)  :: cell
        COMPLEX,       INTENT(IN)  :: vpw(:, :)       ! (ngqdp, 3)
        COMPLEX,       INTENT(IN)  :: rho(:, :, :, :) ! (atoms%jmtd,(atoms%lmaxd + 1)**2, atoms%nat, 3)
        COMPLEX,       INTENT(OUT) :: vr(:, :, :)     ! (atoms%jmtd,(atoms%lmaxd + 1)**2, atoms%nat)
        INTEGER,       INTENT(IN)  :: iDatom, iDtype, iDir
        INTEGER,       INTENT(IN)  :: gqdp(:, :)
        INTEGER,       INTENT(IN)  :: ngqdp
        REAL,          INTENT(IN)  :: qvec(3)
        LOGICAL,       INTENT(IN)  :: l_ext

        COMPLEX                           :: cp
        INTEGER                           :: i, jm, k, l, lh, n, nd, lm, n1, m, imax, lmax, iAtom, iEq, iType, iG
        COMPLEX                           :: vtl((atoms%lmaxd + 1)**2, atoms%nat)
        COMPLEX                           :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%nat)
        REAL                              :: green_factor, gqnorm
        COMPLEX                           :: termsR
        REAL                              :: green_1    (1:atoms%jmtd), green_2    (1:atoms%jmtd)
        REAL                              :: integrand_r1(1:atoms%jmtd), integrand_r2(1:atoms%jmtd)
        REAL                              :: integral_r1 (1:atoms%jmtd), integral_r2 (1:atoms%jmtd)
        REAL                              :: integrand_i1(1:atoms%jmtd), integrand_i2(1:atoms%jmtd)
        REAL                              :: integral_i1 (1:atoms%jmtd), integral_i2 (1:atoms%jmtd)
        COMPLEX                           :: integral_1 (1:atoms%jmtd), integral_2 (1:atoms%jmtd)
        REAL                              :: sbf(0:atoms%lmaxd), gqext(3), gqint(3)

      !$ complex, allocatable :: vtl_loc(:,:)
#ifdef CPP_MPI
        INTEGER :: ierr

        COMPLEX, ALLOCATABLE :: c_b(:)
#endif

        vtl(:,:) = cmplx( 0.0, 0.0 )

         ! TODO: Be careful with this; should be ok though once adding iAtom etc.
!        !$omp parallel default( none ) &
!        !$omp& shared( fmpi, vpw, atoms, cell, vtl ) &
!        !$omp& private( k, cp, pylm, n, sbf, nd, lh, jm, m, lm, l ) &
!        !$omp& private( vtl_loc )
!        !$ allocate(vtl_loc((atoms%lmaxd + 1)**2, atoms%nat))
!        !$ vtl_loc(:,:) = cmplx(0.0,0.0)
!        !$omp do
        DO iG = fmpi%irank+1, ngqdp, fmpi%isize
            cp = vpw(iG, iDir)
            gqint = gqdp(:, iG) + qvec(:)
            gqext = matmul(cell%bmat, gqint)
            gqnorm = norm2(gqext)
            IF (norm2(gqext) .LT. 1e-12) CYCLE
            CALL dfpt_phasy1(atoms, cell, gqint, pylm)

            iAtom = 0
            DO iType = 1, atoms%ntype
                CALL sphbes( atoms%lmax(iType), gqnorm * atoms%rmt(iType), sbf )
                DO iEq = 1, atoms%neq(iType)
                    iAtom = iAtom + 1
                    DO l = 0, atoms%lmax(iType)
                        DO m = -l, l
                            lm = l * ( l + 1 ) + m + 1
!                            !$ if (.false.) then
                            vtl(lm, iAtom) = vtl(lm, iAtom) + cp * sbf(l) * pylm(lm, iAtom)
!                            !$ end if
!                            !$ vtl_loc(lh,n) = vtl_loc(lh,n) + cp * sbf(l) * pylm(lm, iatom)
                        END DO
                    END DO
                END DO
            END DO
        END DO

!        !$omp end do
!        !$omp critical
!        !$ vtl = vtl + vtl_loc
!        !$omp end critical
!        !$ deallocate(vtl_loc)
!        !$omp end parallel

#ifdef CPP_MPI
        n1 = (atoms%lmaxd + 1)**2 * atoms%nat
        ALLOCATE( c_b(n1) )
        CALL MPI_REDUCE( vtl, c_b, n1, CPP_MPI_COMPLEX, MPI_SUM, 0, fmpi%mpi_comm, ierr )
        IF ( fmpi%irank == 0 ) vtl = reshape( c_b, (/(atoms%lmaxd + 1)**2,atoms%nat/) )
        DEALLOCATE( c_b )
#endif

        IF ( fmpi%irank == 0 ) THEN
            iAtom = 0
            DO iType = 1, atoms%ntype
                DO iEq = 1, atoms%neq(iType)
                    iAtom = iAtom + 1
                    imax = atoms%jri(iType)
                    do l = 0, atoms%lmax(iType)
                        green_1(1:imax) = atoms%rmsh(1:imax,iType) ** l
                        green_2(1:imax) = 1.0 / ( green_1(1:imax) * atoms%rmsh(1:imax,iType) )
                        green_factor    = fpi_const / ( 2 * l + 1 )
                        do m = -l, l
                          lm = l * (l + 1) + m + 1
                        integrand_r1(1:imax) = green_1(1:imax) * real(rho(1:imax,lm,iAtom, iDir)) * atoms%rmsh(1:imax, iType)**2
                        integrand_r2(1:imax) = green_2(1:imax) * real(rho(1:imax,lm,iAtom, iDir)) * atoms%rmsh(1:imax, iType)**2
                        integrand_i1(1:imax) = green_1(1:imax) * aimag(rho(1:imax,lm,iAtom, iDir)) * atoms%rmsh(1:imax, iType)**2
                        integrand_i2(1:imax) = green_2(1:imax) * aimag(rho(1:imax,lm,iAtom, iDir)) * atoms%rmsh(1:imax, iType)**2

                        CALL intgr2LinIntp( integrand_r1(1:imax), atoms%rmsh(1,iType), atoms%dx(iType), imax, integral_r1(1:imax) )
                        CALL intgr2LinIntp( integrand_r2(1:imax), atoms%rmsh(1,iType), atoms%dx(iType), imax, integral_r2(1:imax) )
                        CALL intgr2LinIntp( integrand_i1(1:imax), atoms%rmsh(1,iType), atoms%dx(iType), imax, integral_i1(1:imax) )
                        CALL intgr2LinIntp( integrand_i2(1:imax), atoms%rmsh(1,iType), atoms%dx(iType), imax, integral_i2(1:imax) )

                        integral_1(1:imax) = cmplx(integral_r1(1:imax), integral_i1(1:imax))
                        integral_2(1:imax) = cmplx(integral_r2(1:imax), integral_i2(1:imax))

                        ! TODO: Old switch made everything except the vtl-term vanish; removed.
                        termsR = integral_2(imax) + ( vtl(lm, iAtom) / green_factor - integral_1(imax) * green_2(imax) ) / green_1(imax)
                        vr(1:imax, lm, iAtom) = green_factor * (   green_1(1:imax) * ( termsR - integral_2(1:imax) ) &
                                                                 + green_2(1:imax) *            integral_1(1:imax)   )
                        END DO
                    END DO
                END DO
            END DO

            IF (l_ext) THEN
                DO lm = 2, 4
                    vr(:, lm, iDatom) = vr(:, lm, iDatom) &
                                    & - atoms%zatom(iDtype) / atoms%rmsh(:, iDtype)**2 &
                                    & * ( 1 - (atoms%rmsh(:, iDtype) / atoms%rmt(iDtype))**3 ) * c_im(iDir, lm - 1)
                END DO
            END IF
        END IF

    END SUBROUTINE dfpt_vmts

    SUBROUTINE dfpt_psqpw(fmpi, atoms, cell, qpw, rho, rho0MTsh, rho0IRpw, psq, l_ext, iDatom, iDtype, ngqdp, gqdp, ngdp, gdp, qvec)

#include"cpp_double.h"

#ifdef CPP_MPI
        USE mpi
#endif

        USE m_constants
        USE m_sphbes
        USE m_qsf
        USE m_types
        USE m_DoubleFactorial

        IMPLICIT NONE

        TYPE(t_mpi),   INTENT(IN)  :: fmpi
        TYPE(t_atoms), INTENT(IN)  :: atoms
        TYPE(t_cell),  INTENT(IN)  :: cell
        COMPLEX,       INTENT(IN)  :: qpw(ngqdp, 3)
        COMPLEX,       INTENT(IN)  :: rho(atoms%jmtd, ( atoms%lmaxd + 1 ) ** 2, atoms%nat, 3)
        COMPLEX,       INTENT(IN)  :: rho0MTsh(:, :, :) ! (atoms%jmtd, lmd, atoms%nat)
        COMPLEX,       INTENT(IN)  :: rho0IRpw(:)       ! (ngdp)
        COMPLEX,       INTENT(OUT) :: psq(3, ngqdp)
        LOGICAL,       INTENT(IN)  :: l_ext
        INTEGER,       INTENT(IN)  :: iDatom
        INTEGER,       INTENT(IN)  :: iDtype
        INTEGER,       INTENT(IN)  :: ngqdp
        INTEGER,       INTENT(IN)  :: gqdp(:, :)
        INTEGER,         INTENT(IN)   :: ngdp
        INTEGER,         INTENT(IN)   :: gdp(:, :)
        REAL,          INTENT(IN)  :: qvec(3)

        COMPLEX :: sa, sl, sm
        REAL    :: fpo, gqnorm
        INTEGER :: l, n1, ncvn, lm, ll1, nd, m, iAtom, iDir, iEq, iG, iType

        REAL    :: pn(0:atoms%lmaxd,atoms%ntype), gqint(3), gqext(3)
        REAL    :: aj(0:maxval(atoms%ncv)+1)
        COMPLEX :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%nat)
        COMPLEX :: qlm((atoms%lmaxd + 1)**2,atoms%nat,3)

#ifdef CPP_MPI
        INTEGER :: ierr

        COMPLEX, ALLOCATABLE :: c_b(:)
#endif

        ! Calculate vectorial multipole moments
        CALL timestart("dfpt_mpmom")
        CALL dfpt_mpmom( fmpi, atoms, cell, l_ext, iDatom, iDtype, ngqdp, gqdp, ngdp, gdp, qvec, qpw, rho, rho0MTsh, rho0IRpw, qlm )
        CALL timestop("dfpt_mpmom")

#ifdef CPP_MPI
        psq(:, :) = cmplx(0.0, 0.0)
        CALL MPI_BCAST( qpw, size(qpw), CPP_MPI_COMPLEX, 0, fmpi%mpi_comm, ierr )
        nd = (atoms%lmaxd + 1)**2 * atoms%nat * 3
        CALL MPI_BCAST( qlm, nd, CPP_MPI_COMPLEX, 0, fmpi%MPI_COMM, ierr )
#endif

        pn = 0.0
        DO iType = 1, atoms%ntype
            DO l = 0, min( atoms%ncv(iType) - 1, atoms%lmax(iType) )
                pn(l, iType) = DoubleFactorial( atoms%ncv(iType) + 1, l ) / ( atoms%rmt(iType) ** ( atoms%ncv(iType) + 1 ) )
            END DO
        END DO

        fpo = 1. / cell%omtil

        CALL timestart("loop")
        ! TODO: This loop may need to be put further inwards, due to the iAtom counter.
!        !$omp parallel do default( shared ) private( pylm, sa, ncvn, aj, sl, l, n1, ll1, sm, m, lm, iDir)
        DO iG = fmpi%irank+1, ngqdp, fmpi%isize
            gqint = gqdp(:, iG) + qvec(:)
            gqext = matmul(cell%bmat, gqint)
            gqnorm = norm2(gqext)
            IF (norm2(gqext) .LT. 1e-12) CYCLE
            CALL dfpt_phasy1(atoms, cell, gqint, pylm)
            DO iDir =1, 3
                sa = cmplx(0.0, 0.0)
                iAtom = 0
                DO iType = 1, atoms%ntype
                    ncvn = atoms%ncv(iType)
                    CALL sphbes( ncvn + 1 , gqnorm * atoms%rmt(iType), aj )
                    DO iEq = 1, atoms%neq(iType)
                        iAtom = iAtom + 1
                        sl = cmplx(0.0, 0.0)
                        DO l = 0, atoms%lmax(iType)
                            IF ( l >= ncvn ) CYCLE
                            n1 = ncvn - l + 1
                            ll1 = l * ( l + 1 ) + 1
                            sm = cmplx(0.0, 0.0)
                            DO m = -l, l
                                lm = ll1 + m
                                sm = sm + qlm(lm, iAtom, iDir) * conjg(pylm(lm, iAtom))
                            END DO
                            sl = sl + pn(l,iType) / ( (gqnorm) ** n1 ) * aj( ncvn + 1 ) * sm
                        END DO
                        sa = sa + sl
                    END DO
                END DO
                psq(iDir, iG) = qpw(iG, iDir) + fpo * sa
            END DO
        END DO
!        !$omp end parallel do

        CALL timestop("loop")

#ifdef CPP_MPI
        ALLOCATE(c_b(3 * ngqdp))
        CALL MPI_REDUCE( psq, c_b, 3 * ngqdp, CPP_MPI_COMPLEX, MPI_SUM, 0, fmpi%MPI_COMM, ierr )
        IF ( fmpi%irank == 0 ) THEN
            psq = reshape(c_b, (/3,ngqdp/))
        END IF
        DEALLOCATE(c_b)
#endif

    END SUBROUTINE dfpt_psqpw

    SUBROUTINE dfpt_mpmom(fmpi, atoms, cell, l_ext, iDatom, iDtype, ngqdp, gqdp, ngdp, gdp, qvec, qpw, rho, rho0MTsh, rho0IRpw, qlm)

        USE m_types
        USE m_constants

        IMPLICIT NONE

        TYPE(t_mpi),     INTENT(IN)   :: fmpi
        TYPE(t_cell),    INTENT(IN)   :: cell
        TYPE(t_atoms),   INTENT(IN)   :: atoms

        LOGICAL,         INTENT(IN)   :: l_ext
        INTEGER,         INTENT(IN)   :: iDatom
        INTEGER,         INTENT(IN)   :: iDtype
        INTEGER,         INTENT(IN)   :: ngqdp
        INTEGER,         INTENT(IN)   :: gqdp(:, :)
        INTEGER,         INTENT(IN)   :: ngdp
        INTEGER,         INTENT(IN)   :: gdp(:, :)
        REAL,            INTENT(IN)   :: qvec(3)
        COMPLEX,         INTENT(IN)   :: rho(:,:,:,:)      ! (atoms%jmtd, lmd, atoms%nat, 3)
        COMPLEX,         INTENT(IN)   :: qpw(:, :)         ! (ngqdp, 3)
        COMPLEX,         INTENT(IN)   :: rho0MTsh(:, :, :) ! (atoms%jmtd, lmd, atoms%nat)
        COMPLEX,         INTENT(IN)   :: rho0IRpw(:)       ! (ngdp)
        COMPLEX,         INTENT(OUT)  :: qlm((atoms%lmaxd + 1)**2,atoms%nat,3)

        INTEGER :: iDir, iType, iEq, iAtom, l, m, ll1, lm
        COMPLEX :: qlmo((atoms%lmaxd + 1)**2,atoms%nat,3)
        COMPLEX :: qlmp((atoms%lmaxd + 1)**2,atoms%nat,3)

        ! Vectorial multipole moments of original charge density perturbation and atomic charge
        IF ( fmpi%irank == 0 ) then
            CALL dfpt_mt_moments(atoms, l_ext, iDatom, rho, qlmo)
            CALL dfpt_mt_moments_SF(atoms, iDtype, iDatom, rho0MTsh, qlmo)
        END IF

        ! Vectorial multipole moments of the interstitial charge density in the spheres
        CALL dfpt_pw_moments( fmpi, atoms, cell, ngqdp, qvec, gqdp, qpw, qlmp )
        CALL dfpt_pw_moments_SF( atoms, cell, ngdp, iDtype, iDatom, gdp, rho0IRpw, qlmp )

        IF ( fmpi%irank == 0 ) then
            qlm = qlmo - qlmp
            DO iDir = 1, 3
                iAtom = 0
                DO iType = 1, atoms%ntype
                    DO iEq =1, atoms%neq(iType)
                        iAtom = iAtom + 1
                        WRITE(oUnit, fmt=8000) iAtom, iDir
                        DO l = 0, atoms%lmax(iType)
                            DO m = -l, l
                                lm = l * (l + 1) + 1 + m
                                IF ((qlmo(lm, iAtom, iDir)/=CMPLX(0.0)).OR.(qlmp(lm, iAtom, iDir)/=CMPLX(0.0))) THEN
                                    WRITE(oUnit, fmt=8010 ) l, m, qlmo(lm, iAtom, iDir), qlmp(lm, iAtom, iDir)
                                END IF
                            END DO
                        END DO
                    END DO
                END DO
            END DO

8000 FORMAT (/,10x,'Multipole moments for atom',i5,'and direction',i1,/,/,t3,'l',t7,'m',t27,'original',t57,'plane wave')
8010 FORMAT (1x,i2,2x,i3,2x,2 (5x,2e15.5))

        END IF ! fmpi%irank == 0

    END SUBROUTINE dfpt_mpmom

    SUBROUTINE dfpt_mt_moments(atoms, l_ext, iDatom, rho, qlmo)

        USE m_intgr,     only: intgr3LinIntp ! TODO: Once fed with behaving quantities, use intgr3.
        USE m_constants, only: sfp_const
        USE m_types
        USE m_juDFT

        IMPLICIT NONE

        TYPE(t_atoms),  INTENT(IN)  :: atoms
        LOGICAL,        INTENT(IN)  :: l_ext
        INTEGER,        INTENT(IN)  :: iDatom
        COMPLEX,        INTENT(IN)  :: rho(: , :, :, :)
        COMPLEX,        INTENT(OUT) :: qlmo(:, :, :)

        INTEGER :: iDir, iType, iEq, iAtom, l, m, lm, iR
        REAL    :: fintReal, fintImag
        REAL    :: fReal(atoms%jmtd), fImag(atoms%jmtd)

        qlmo = cmplx(0.0, 0.0)

        DO iDir = 1, 3
            iAtom = 0
            DO iType = 1, atoms%ntype
                DO iEq =  1, atoms%neq(iType)
                    iAtom = iAtom + 1
                    DO l = 0, atoms%lmax(iType)
                        DO m = -l, l
                            lm = l * (l + 1) + 1 + m
                            DO iR = 1, atoms%jri(iType)
                                fReal(iR) = atoms%rmsh(iR, iType)**(l + 2) *  real(rho(iR, lm, iAtom, iDir))
                                fImag(iR) = atoms%rmsh(iR, iType)**(l + 2) * aimag(rho(iR, lm, iAtom, iDir))
                            END DO
                            CALL intgr3LinIntp(fReal, atoms%rmsh(:, iType), atoms%dx(iType), atoms%jri(iType), fintReal, 1)
                            CALL intgr3LinIntp(fImag, atoms%rmsh(:, iType), atoms%dx(iType), atoms%jri(iType), fintImag, 1)
                            qlmo(lm, iAtom, iDir) = cmplx(fintReal, fintImag)
                        END DO
                    END DO ! DO l = 0, atoms%lmax(iType)
                    IF (l_ext.AND.(iAtom.EQ.iDatom)) THEN
                        qlmo(2:4, iDatom, iDir) = qlmo(2:4, iDatom, iDir) - 3.0 / fpi_const * atoms%zatom(iType) * c_im(iDir, :)
                    END IF
                END DO ! iEq =  1, atoms%neq(iType)
            END DO ! iType = 1, atoms%ntype
        END DO ! iDir = 1, 3

    END SUBROUTINE dfpt_mt_moments

    SUBROUTINE dfpt_pw_moments(fmpi, atoms, cell, ngqdp, qvec, gqdp, qpw_in, qlmp_out)

#ifdef CPP_MPI
        USE mpi
#endif

        USE m_sphbes
        USE m_types

        IMPLICIT NONE

        TYPE(t_mpi),   INTENT(IN)  :: fmpi
        TYPE(t_atoms), INTENT(IN)  :: atoms
        TYPE(t_cell),  INTENT(IN)  :: cell
        INTEGER,       INTENT(IN)  :: ngqdp
        REAL,          INTENT(IN)  :: qvec(3)
        INTEGER,       INTENT(IN)  :: gqdp(:, :)
        COMPLEX,       INTENT(IN)  :: qpw_in(:, :)
        COMPLEX,       INTENT(OUT) :: qlmp_out(:, :, :)

        INTEGER :: iG, iDir, iType, iEq, iAtom, l, m, ll1, lm, ierr
        REAL    :: sk3r, rl2, gqnorm
        COMPLEX :: sk3i, cil

        REAL    :: aj(0:atoms%lmaxd + 1 ), gqint(3), gqext(3)
        COMPLEX :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%nat)
        COMPLEX :: qlmp(( atoms%lmaxd + 1 ) ** 2, atoms%nat, 3)

        qlmp = cmplx(0.0, 0.0)

#ifdef CPP_MPI
        CALL MPI_BCAST( qpw_in, size(qpw_in), MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr )
#endif

        DO iG = fmpi%irank+1, ngqdp, fmpi%isize
            gqint = gqdp(:, iG) + qvec(:)
            gqext = matmul(cell%bmat, gqint)
            gqnorm = norm2(gqext)
            IF (gqnorm .LT. 1e-12) CYCLE

            CALL dfpt_phasy1(atoms, cell, gqint, pylm)

            DO iDir = 1, 3
                iAtom = 0
                DO iType = 1, atoms%ntype
                    sk3r = gqnorm * atoms%rmt(iType)
                    CALL sphbes( atoms%lmax(iType) + 1, sk3r, aj )
                    sk3i = qpw_in(iG, iDir) / gqnorm
                    rl2 = atoms%rmt(iType) ** 2
                    DO iEq = 1, atoms%neq(iType)
                        iAtom = iAtom + 1
                        DO l = 0, atoms%lmax(iType)
                            cil = aj(l+1) * sk3i * rl2
                            ll1 = l * ( l + 1 ) + 1
                            DO m = -l, l
                                lm = ll1 + m
                                qlmp(lm, iAtom, iDir) = qlmp(lm, iAtom, iDir) + cil * pylm(lm, iAtom)
                            END DO
                            rl2 = rl2 * atoms%rmt(iType)
                        END DO ! l = 0, atoms%lmax(n)
                    END DO ! iEq = 1, atoms%neq(iType)
                END DO ! iType = 1, atoms%ntype
            END DO ! iDir = 1, 3
        END DO ! iG = 1, ngqdp

#ifdef CPP_MPI
        CALL MPI_REDUCE( qlmp, qlmp_out, SIZE(qlmp), MPI_DOUBLE_COMPLEX, MPI_SUM,0, fmpi%mpi_comm, ierr )
#else
        qlmp_out = qlmp
#endif

    END SUBROUTINE dfpt_pw_moments

    SUBROUTINE dfpt_mt_moments_SF(atoms, iDtype, iDatom, rho0MTsh, qlmo)

        USE m_types, only : t_atoms
        USE m_gaunt, only : Gaunt1
        USE m_constants

        IMPLICIT NONE

        TYPE(t_atoms), INTENT(IN)    :: atoms
        INTEGER,       INTENT(IN)    :: iDtype
        INTEGER,       INTENT(IN)    :: iDatom
        COMPLEX,       INTENT(IN)    :: rho0MTsh(:, :, :)
        COMPLEX,       INTENT(INOUT) :: qlmo(:, :, :)

        INTEGER :: l, lp, m, mp, m2p, lm, lm1p, iDir
        REAL    :: gauntFactor

        DO l = 0, atoms%lmax(iDtype)
            DO m = -l, l
                lm = l * (l + 1) + m + 1
                DO lp = 0, atoms%lmax(iDtype)
                    DO mp = -lp, lp
                        lm1p = lp * (lp + 1) + mp + 1
                        DO m2p = -1, 1
                            gauntFactor = Gaunt1( l, lp, 1, m, mp, m2p, atoms%lmax(iDtype))
                            DO iDir = 1, 3
                                qlmo(lm, iDatom, iDir) = qlmo(lm, iDatom, iDir) + &
                                                       & c_im(iDir, m2p + 2) * atoms%rmt(iDtype)**(l + 2) * &
                                                       & rho0MTsh(atoms%jri(iDtype), lm1p, iDatom) * gauntFactor
                            END DO ! iDir
                        END DO ! m2p
                    END DO ! mp
                END DO ! lp
            END DO ! m
        END DO ! l

    END SUBROUTINE dfpt_mt_moments_SF

    SUBROUTINE dfpt_pw_moments_SF( atoms, cell, ngdp, iDtype, iDatom, gdp, rho0IRpw, qlmp )

        USE m_types,  only : t_atoms, t_cell
        USE m_ylm,    only : ylm4
        USE m_sphbes, only : sphbes
        USE m_gaunt,  only : Gaunt1
        USE m_constants

        IMPLICIT NONE

        TYPE(t_atoms), INTENT(IN)    :: atoms
        TYPE(t_cell),  INTENT(IN)    :: cell
        INTEGER,       INTENT(IN)    :: ngdp
        INTEGER,       INTENT(IN)    :: iDtype
        INTEGER,       INTENT(IN)    :: iDatom
        INTEGER,       INTENT(IN)    :: gdp(:, :)
        COMPLEX,       INTENT(IN)    :: rho0IRpw(:)
        COMPLEX,       INTENT(INOUT) :: qlmp(:, :,:)

        INTEGER :: iG, l, m, lp, mp, m2p, lm, lmp, iDir
        COMPLEX :: pref, phaseFac, temp1, temp2, temp3
        REAL    :: gnorm, gauntFactor

        REAL    :: gext(3)
        REAL    :: sbes(0:atoms%lmax(iDtype))
        COMPLEX :: ylm((atoms%lmax(iDtype) + 1)**2)


        pref = fpi_const * atoms%rmt(iDtype) * atoms%rmt(iDtype)
        DO iG = 1, ngdp
            gext = matmul(cell%bmat, real(gdp(:, iG)))
            gnorm = norm2(gExt)

            call ylm4( atoms%lmax(iDtype), gExt(1:3), ylm )
            call sphbes(atoms%lmax(iDtype), gnorm * atoms%rmt(iDtype), sbes)

            phaseFac = exp( ImagUnit * tpi_const * dot_product(gdp(:, iG), atoms%taual(:, iDatom)))

            DO l = 0, atoms%lmax(iDtype)
                temp1 = pref * phaseFac * atoms%rmt(iDtype)**l * rho0IRpw(iG)
                DO m = -l, l
                    lm = l * (l + 1) + m + 1
                    DO lp = 0, atoms%lmax(iDtype)
                        temp2 = temp1 * sbes(lp) * ImagUnit**lp
                        DO mp = -lp, lp
                            lmp = lp * (lp + 1) + mp + 1
                            temp3 = temp2 * conjg(ylm(lmp))
                            DO m2p = -1, 1
                                gauntFactor = Gaunt1( l, lp, 1, m, mp, m2p, atoms%lmax(iDtype))
                                DO iDir = 1, 3
                                    qlmp(lm, iDatom, iDir) = qlmp(lm, iDatom, iDir) + c_im(iDir, m2p + 2) * gauntFactor * temp3
                                END DO ! iDir
                            END DO ! m2p
                        END DO ! mp
                    END DO ! lp
                END DO ! m
            END DO ! l
        END DO ! iG

    END SUBROUTINE dfpt_pw_moments_SF

    SUBROUTINE dfpt_phasy1(atoms, cell, gq, pylm)

        USE m_constants
        USE m_types
        USE m_ylm

        TYPE(t_atoms), INTENT(IN)  :: atoms
        TYPE(t_cell),  INTENT(IN)  :: cell
        REAL,          INTENT(IN)  :: gq(3)
        COMPLEX, 	INTENT(OUT) :: pylm(:,:)

        REAL    :: x
        COMPLEX :: sf,csf

        INTEGER :: iType, iEq, iAtom, l, m, lm, ll1

        COMPLEX :: ciall(0:atoms%lmaxd)
        REAL    :: Gqext(3)

        COMPLEX, ALLOCATABLE :: ylm(:)

        ciall(0) = fpi_const
        DO l = 1,atoms%lmaxd
            ciall(l) = ciall(0)*ImagUnit**l
        END DO

        ALLOCATE (ylm((atoms%lmaxd + 1)**2))
        ylm   = cmplx(0.0, 0.0)
        gqext = matmul(gq, cell%bmat)
        CALL ylm4(atoms%lmaxd, Gqext, ylm)
        ylm   = conjg( ylm )

        iAtom = 0
        pylm(:, :) = cmplx(0.0, 0.0)
        DO iType = 1, atoms%ntype
            DO iEq = 1, atoms%neq(iType)
                iAtom = iAtom + 1
                x = tpi_const * dot_product(gq, atoms%taual(:, iAtom))
                sf = exp(ImagUnit *  x)
                DO l = 0, atoms%lmax(iType)
                    ll1 = l * (l + 1) + 1
                    csf = ciall(l) * sf
                    DO m = -l, l
                        lm = ll1 + m
                        pylm(lm, iatom) = csf * ylm(lm)
                    END DO ! m
                END DO ! l
            END DO ! iEq
        END DO ! iType

        DEALLOCATE ( ylm )

    END SUBROUTINE dfpt_phasy1

END MODULE m_dfpt_vgen_coulomb
