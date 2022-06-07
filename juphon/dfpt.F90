!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt
   USE m_juDFT_time
   USE m_constants
   USE m_types
   USE m_dfpt_check
   USE m_dfpt_test
   USE m_dfpt_init
   USE m_dfpt_sternheimer
   USE m_dfpt_dynmat
   USE m_jpSternheimer,     only : solveSternheimerSCC
   USE m_jp2ndOrdQuant,     only : CalcIIEnerg2
   USE m_jpSetupDynMat,     only : SetupDynamicMatrix
   USE m_jpProcessDynMat!,   only : DiagonalizeDynMat, CalculateFrequencies
   USE m_juDFT_stop, only : juDFT_error

   IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt(fi, sphhar, stars, nococonv, qpts, fmpi, results, enpara, &
                 & rho, vTot, vCoul, vxc, exc, eig_id, nvfull, GbasVec_eig, oldmode, xcpot, hybdat, mpdata, forcetheo)

      TYPE(t_fleurinput), INTENT(IN)  :: fi
      TYPE(t_sphhar),   INTENT(IN)  :: sphhar
      TYPE(t_stars),    INTENT(IN)  :: stars
      TYPE(t_nococonv), INTENT(IN)  :: nococonv
      TYPE(t_kpts),     INTENT(IN)  :: qpts
      TYPE(t_mpi),      INTENT(IN)  :: fmpi
      TYPE(t_results),  INTENT(INOUT)  :: results
      TYPE(t_enpara),   INTENT(INOUT)  :: enpara
      TYPE(t_potden),   INTENT(IN)  :: rho, vTot, vCoul, vxc, exc
      INTEGER,          INTENT(IN)  :: eig_id
      INTEGER,          INTENT(IN)  :: nvfull(:, :), GbasVec_eig(:, :, :, :)

      ! New input:
      LOGICAL, INTENT(IN)        :: oldmode
      CLASS(t_xcpot), INTENT(IN) :: xcpot
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat
      TYPE(t_mpdata),     INTENT(INOUT) :: mpdata
      CLASS(t_forcetheo), INTENT(INOUT) :: forcetheo

      TYPE(t_usdus)                 :: usdus
      TYPE(t_potden)                :: vTotclean, rhoclean, grRho
      TYPE(t_potden)                :: grRho3(3), grVtot3(3)
      TYPE(t_jpPotden)              :: rho0, grRho0, vTot0, grVTot0
      TYPE(t_tlmplm)                :: tdHS0
      TYPE(t_results)               :: results1

        integer                       :: logUnit = 100
        integer                       :: ngpqdp

        COMPLEX, ALLOCATABLE          :: loosetd(:, :, :, :)
        REAL,             ALLOCATABLE :: El(:, :, :, :)
        INTEGER,          ALLOCATABLE :: recG(:, :), GbasVec(:, :), ilst(:, :, :)
        INTEGER                       :: ngdp2km
        INTEGER,          ALLOCATABLE :: gdp2Ind(:, :, :)
        INTEGER                       :: gdp2iLim(2, 3)
        INTEGER,          ALLOCATABLE :: nRadFun(:, :), iloTable(:, :, :), ilo2p(:, :)
        REAL,             ALLOCATABLE :: uuilon(:, :)
        REAL,             ALLOCATABLE :: duilon(:, :)
        REAL,             ALLOCATABLE :: ulouilopn(:, :, :)
        INTEGER,          ALLOCATABLE :: kveclo(:,:)
        REAL,             ALLOCATABLE :: rbas1(:, :, :, :, :)
        REAL,             ALLOCATABLE :: rbas2(:, :, :, :, :)
        REAL,             ALLOCATABLE :: gridf(:, :)
        COMPLEX,          ALLOCATABLE :: z0(:, :, :, :)
        complex,           allocatable :: grVxcIRKern(:)
        real,              allocatable :: dKernMTGPts(:, :, :)
        real,              allocatable :: gausWts(:)
        complex,           allocatable :: ylm(:, :)
        complex,           allocatable :: qpwcG(:, :)
        complex,           allocatable :: rho1MTCoreDispAt(:, :, :, :)
        complex,           allocatable :: grVeff0MT_init(:, :, :, :)
        complex,           allocatable :: grVeff0MT_main(:, :, :, :)
        complex,           allocatable :: grVext0IR_DM(:, :)
        complex,           allocatable :: grVext0MT_DM(:, :, :, :)
        complex,           allocatable :: grVCoul0IR_DM_SF(:, :)
        complex,           allocatable :: grVCoul0MT_DM_SF(:, :, :, :)
        complex,           allocatable :: grVeff0IR_DM(:, :)
        complex,           allocatable :: grVeff0MT_DM(:, :, :, :)
        complex,           allocatable :: dynMat(:, :)
        complex,           allocatable :: E2ndOrdII(:, :)
        complex,           allocatable :: eigenFreqs(:)
        real,              allocatable :: eigenVals(:)
        complex,           allocatable :: eigenVecs(:, :)
        integer,           allocatable :: gpqdp(:, :)
        INTEGER,           ALLOCATABLE :: nocc(:, :)
        complex,            allocatable :: rho1IR(:, :, :)
        complex,            allocatable :: rho1MT(:, :, :, :, :)
        complex,            allocatable :: rho1MTDelta(:, :, :, :, :)
        complex,            allocatable :: rho1MTz0(:, :, :, :)
        complex,            allocatable :: vCoul1IRtempNoVol(:, :)
        complex,            allocatable :: vCoul1MTtempNoVol(:, :, :, :)
        complex,            allocatable :: vEff1IR(:, :, :)
        complex,            allocatable :: vEff1MT(:, :, :, :, :)
        complex,            allocatable :: vEff1MTnoVol(:, :, :, :, :)
        complex,            allocatable :: vExt1IR_final(:, :, :)
        complex,            allocatable :: vExt1MT(:, :, :, :, :)
        complex,            allocatable :: vExt1MTDelta(:, :, :, :, :)
        complex,            allocatable :: vExt1MTq0(:, :, :, :, :)
        complex,            allocatable :: vExt1noqIR_final(:, :, :)
        complex,            allocatable :: vHar1IR_final(:, :, :)
        complex,            allocatable :: vHar1MTDelta(:, :, :, :, :)
        complex,            allocatable :: vHar1MTq0(:, :, :, :, :)
        complex,            allocatable :: vHar1MT_final(:, :, :, :, :)
        complex,            allocatable :: vXc1MTDelta(:, :, :, :, :)
        complex,            allocatable :: vXc1MTq0(:, :, :, :, :)
        integer,      allocatable :: mapKpq2K(:, :)
        integer,      allocatable :: kpq2kPrVec(:, :, :)

        REAL, ALLOCATABLE :: dynmatrow(:)

        INTEGER :: ngdp, iSpin, iType, iR, ilh, iQ, iDir, iDtype

        CHARACTER(len=20) :: dfpt_tag

        INTEGER, ALLOCATABLE :: q_list(:)

#ifndef CPP_FFTW
        call juDFT_error('juPhon is only usable with fftw support.', calledby='dfpt')
#endif

      CALL results1%init(fi%input, fi%atoms, fi%kpts, fi%noco)

        WRITE (oUnit,*) '------------------------------------------------------'
        WRITE (oUnit,*) 'This output is generated by juPhon, FLEURs DFPT addon.'
        !WRITE (oUnit,*) 'l_dfpt = ', fi%juPhon%l_dfpt
        !WRITE (oUnit,*) 'l_jpCheck = ', fi%juPhon%l_jpCheck
        !WRITE (oUnit,*) 'l_jpTest = ', fi%juPhon%l_jpTest
        !WRITE (oUnit,*) 'l_potout = ', fi%juPhon%l_potout
        !WRITE (oUnit,*) 'l_eigout = ', fi%juPhon%l_eigout
        !WRITE (oUnit,*) 'l_symTsh = ', fi%juPhon%l_symTsh
        !WRITE (oUnit,*) 'l_symTdm = ', fi%juPhon%l_symTdm
        !WRITE (oUnit,*) 'l_bfkq = ', fi%juPhon%l_bfkq
        !WRITE (oUnit,*) 'jplPlus = ', fi%juPhon%jplPlus
        !WRITE (oUnit,*) 'kgqmax = ', fi%juPhon%kgqmax
        !WRITE (oUnit,*) 'gqmax = ', fi%juPhon%gqmax
        !WRITE (oUnit,*) 'eps_pert = ', fi%juPhon%eps_pert
        !WRITE (oUnit,*) 'eDiffCut = ', fi%juPhon%eDiffCut
        !WRITE (oUnit,*) 'qpt_ph = ', fi%juPhon%qpt_ph

      ! TODO: This is a test set of qpoints for a fixed fcc system.
      !       We need to read out actual q-points at some point.
      ALLOCATE(q_list(5))
      q_list = [1, 10, 19, 28, 37]! 512 k-points: \Gamma to X

      ALLOCATE(dynmatrow(3*fi%atoms%ntype))
      DO iQ = 1, SIZE(q_list)
         DO iDtype = 1, fi%atoms%ntype
            DO iDir = 1, 3
               dfpt_tag = ''
               WRITE(dfpt_tag,'(a1,i0,a2,i0,a2,i0)') 'q', q_list(iQ), '_b', iDtype, '_j', iDir
               !WRITE(8001,*) dfpt_tag
               ! This is where the magic happens. The Sternheimer equation is solved
               ! iteratively, providing the scf part of dfpt calculations.
               CALL dfpt_sternheimer(fi, xcpot, sphhar, stars, nococonv, qpts, fmpi, results, enpara, hybdat, mpdata, forcetheo, &
                                     rho, vTot, grRho3(iDir), grVtot3(iDir), q_list(iQ), iDtype, iDir, &
                                     dfpt_tag, eig_id, results1, dynmatrow)
               ! Once the first order quantities are converged, we can construct all
               ! additional quantities necessary and from that the dynamical matrix.
               CALL dfpt_dynmat()
            END DO
         END DO
      END DO
        IF (fi%juPhon%l_jpCheck) THEN
            ! This function will be used to check the validity of juPhon's
            ! input. I.e. check, whether all prohibited switches are off and,
            ! once there is more expertise considering this topic, check whether
            ! the cutoffs are chosen appropriately.
            CALL dfpt_check(fi, xcpot)
        END IF

        ! Construct potential without the l=0 prefactor.
        CALL vTotclean%copyPotDen(vTot)
        CALL rhoclean%copyPotDen(rho)

        DO iSpin = 1, fi%input%jspins
            DO iType = 1, fi%atoms%ntype
                DO ilh = 0, sphhar%nlhd
                    DO iR = 1, fi%atoms%jri(iType)
                        IF (ilh.EQ.0) THEN
                            vTotclean%mt(iR, 0, iType, iSpin) &
                        & = vTotclean%mt(iR, 0, iType, iSpin) * sqrt(fpi_const) / fi%atoms%rmsh(iR, iType)
                        END IF
                        rhoclean%mt(iR, ilh, iType, iSpin) &
                    & = rhoclean%mt(iR, ilh, iType, iSpin) / fi%atoms%rmsh(iR, iType) / fi%atoms%rmsh(iR, iType)
                    END DO
                END DO
            END DO
        END DO

        ! This routine will initialize everything for juPhon that isn't already
        ! provided by the FLEUR code/must be provieded in a modified way.
        ! This includes for example the de-symmetrized MT and pw quantities and
        ! their gradients. Notably, q-dependent quantities are initialized and
        ! constructed elsewhere, within the q-loop.
        ! TODO: I ignored the actual significance of clnu_atom etc. They are not type-dependent, but actually
        ! refer to each atom respectively. So this will explode for iatom > 1. This is easily fixed.
        CALL timestart("juPhon DFPT initialization")
        CALL dfpt_init(fi%juPhon, fi%sym, fi%input, fi%atoms, sphhar, stars, fi%cell, fi%noco, nococonv, fi%kpts, &
                     & fmpi, results, enpara, rho, vTot, eig_id, nvfull, GbasVec_eig, usdus, rho0, grRho0, vTot0, grVTot0, &
                     & ngdp, El, recG, ngdp2km, gdp2Ind, gdp2iLim, GbasVec, ilst, nRadFun, iloTable, ilo2p, &
                     & uuilon, duilon, ulouilopn, kveclo, rbas1, rbas2, gridf, z0, grVxcIRKern, dKernMTGPts, &
                     & gausWts, ylm, qpwcG, rho1MTCoreDispAt, grVeff0MT_init, grVeff0MT_main, grVext0IR_DM, grVext0MT_DM, &
                     & grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, grVeff0IR_DM, grVeff0MT_DM, tdHS0, loosetd, nocc, rhoclean, oldmode, xcpot, grRho)
        CALL timestop("juPhon DFPT initialization")

        IF (fi%juPhon%l_jpTest) THEN
            ! This function will be used to run (parts of) the test suite for
            ! OG juPhon, as provided by CRG.
            CALL dfpt_test(fi, sphhar, stars, fmpi, rho, grRho, rho0, grRho0, xcpot, ngdp, recG, grVxcIRKern, ylm, dKernMTGPts, gausWts, hybdat)
        END IF
        ! < Imagine starting a q-grid-loop here. >
        ! < For now we just select one q-point from the input. >

        call createkqMapArrays( fi%kpts, qpts, 0, fi%kpts%nkpt3, [0], mapKpq2K, kpq2kPrVec )

        CALL timestart("juPhon DFPT scf loop")
        call solveSternheimerSCC( fmpi,  fi%atoms, fi%sym, stars, sphhar, fi%cell, enpara, usdus, fi%input, fi%kpts, qpts, results, usdus,      &
          & logUnit, ngdp, rbas1, rbas2, kveclo, uuilon, duilon, ulouilopn, &
          & recG, mapKpq2K, results%neig(:, 1), results%eig, GbasVec, ilst, z0, nvfull, El, nradFun, iloTable, nocc, ilo2p, gdp2Ind,     &
          & gdp2iLim, kpq2kPrVec, qpwcG, fi%juPhon%singleQpt, tdHS0, loosetd, ylm, grRho0%pw(:, 1, 1, :), grRho0%mt(:, :, :, 1, 1, :), grVeff0MT_init, grVeff0MT_main, dKernMTGPts,       &! main --> DM for Alex
          & grVxcIRKern, rho1MTCoreDispAt, gausWts, rho1IR, rho1MT, vExt1MT, vEff1IR, vEff1MT, fi%juPhon%oneSternhCycle, ngpqdp, gpqdp,&
          & vExt1IR_final, vHar1IR_final, vHar1MT_final, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, vXc1MTDelta, &
          & vXc1MTq0, rho0%pw(:, :, 1, 1), rho0%mt(:, :, :, :, 1, 1), vTot0%pw(:, :, 1, 1), fi%juPhon%noPtsCon, vEff1MTnoVol, vExt1noqIR_final, rho1MTz0, vCoul1IRtempNoVol, vCoul1MTtempNoVol )
        CALL timestop("juPhon DFPT scf loop")

        CALL timestart("juPhon DFPT Eii2")
        CALL CalcIIEnerg2(fi%atoms, fi%cell, qpts, stars, fi%input, fi%juPhon%singleQpt, ngdp, recG, E2ndOrdII)
        CALL timestop("juPhon DFPT Eii2")

        CALL timestart("juPhon DFPT dynmat setup")
        CALL SetupDynamicMatrix( fmpi, fi%noco, nococonv,  fi%atoms, fi%input, fi%sym, fi%cell, sphhar, stars, fi%kpts, qpts, usdus, results, vTotclean, fi%juPhon%singleQpt, ngdp, ngpqdp, recG, sphhar%mlh, sphhar%nmem,&
            & sphhar%clnu, rho%pw, rho1IR, rho1MT, vExt1MT, vEff1IR, vEff1MT, vTot%pw_w, vTotclean%mt(:, 0:, :, 1),&
            & rhoclean%mt, E2ndOrdII, El, results%eig, rbas1, rbas2, iloTable, nvfull, nocc, ilst, GbasVec, z0, kveclo, nRadFun, mapKpq2K, kpq2kPrVec,       &
            & gpqdp, sphhar%memd, logUnit, vxc%pw, exc%pw(:, 1), vxc%mt, exc%mt(:, 0:, :, 1), vExt1IR_final, vHar1IR_final, vHar1MT_final, grRho0%pw(:, 1, 1, :), grRho0%mt(:, :, :, 1, 1, :), &
            & grVext0IR_DM, grVext0MT_DM, grVeff0IR_DM, grVeff0MT_DM, dynMat, rho1MTDelta, vExt1MTDelta, vExt1MTq0, vHar1MTDelta, vHar1MTq0, &
            & vXc1MTDelta, vXc1MTq0, vEff1MTnoVol, vExt1noqIR_final, rho1MTz0, &
            & grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, vCoul1IRtempNoVol, vCoul1MTtempNoVol)
        CALL timestop("juPhon DFPT dynmat setup")

        CALL timestart("juPhon DFPT dynmat diagonalization")
        CALL DiagonalizeDynMat(fi%atoms, qpts, fi%juPhon%calcEigenVec, dynMat, eigenVals, eigenVecs, fi%juPhon%singleQpt)
        CALL timestop("juPhon DFPT dynmat diagonalization")

        CALL timestart("juPhon DFPT frequency calculation")
        CALL CalculateFrequencies(fi%atoms, fi%juPhon%singleQpt, eigenVals, eigenFreqs)
        CALL timestop("juPhon DFPT frequency calculation")

        ! < Imagine ending a q-grid-loop here. >

        WRITE (oUnit,*) '------------------------------------------------------'

        CALL juDFT_end("Phonon calculation finished.")

    END SUBROUTINE dfpt

    subroutine createkqMapArrays( kpts, qpts, nrAddQs, kSetDim, addQnkptis, mapKpq2K, kpq2kPrVec )

      use m_types
      use m_juDFT_stop, only : juDFT_error

      implicit none

      ! Type parameters
      type(t_kpts),              intent(in)  :: kpts
      type(t_kpts),              intent(in)  :: qpts

      ! Array parameter
      integer,                   intent(in)  :: addQnkptis(:)
      integer,                   intent(in)  :: kSetDim(:)
      integer,      allocatable, intent(out) :: mapKpq2K(:, :)
      integer,      allocatable, intent(out) :: kpq2kPrVec(:, :, :)

      ! Scalar parameter
      integer,                   intent(in)  :: nrAddQs

      ! Array variable
      integer,      allocatable              :: mapK2Ind(:, :, :) ! takes components of kpts%bk * lcm and gives its index
      integer,      allocatable              :: mapK2mK(:)
      integer                                :: kpqNomin(3)       ! helps to find k + q mapped back to Brillouin zone
      real                                   :: kpqTemp(3)        ! helps to find k + q mapped back to Brillouin zone
      character(len=1024)                    :: errorMessage      ! stores error message for error output

      ! Scalar variable
      integer                                :: ikpt              ! loop variable
      integer                                :: maxKcomp(3)       ! stores maximal k-point component * lcm
      integer                                :: idir              ! loop variable
      integer                                :: iqpt              ! loop variable
      integer                                :: nkptShift         ! stores shift of shifted k-point set
      integer                                :: ikptSh            ! loop variable
      logical                                :: matchFound        ! is true if index for k + q is found
      real :: lcm

      lcm = real( kgv(kSetDim, 3) )

      ! Determine maximal value of k-vector per direction in internal representation for allocation of mapKpq2K array and allocate it
      maxKcomp = 0
      do idir = 1, 3
        maxKcomp(idir) = maxval( kpts%bk(idir, :kpts%nkpt) * lcm )
      end do
      allocate( mapK2Ind(0:maxKcomp(1), 0:maxKcomp(2), 0:maxKcomp(3)) )
      mapK2Ind = 0
      allocate( mapKpq2K(kpts%nkpt, qpts%nkptf + nrAddQs) )

      ! Fill up array which stores the index of a given k-point
      do ikpt = 1, kpts%nkpt
        mapK2Ind(nint(kpts%bk(1, ikpt) * lcm), nint(kpts%bk(2, ikpt) * lcm), nint(kpts%bk(3, ikpt) * lcm)) = ikpt
      end do

      ! Determine k-point on which k + q can be folded back and determine the respective reciprocal lattice vector.
      ! The absolute value of every coordinate of the reciprocal lattice vector can be 1 maximally.
      allocate( kpq2kPrVec(3, kpts%nkpt, qpts%nkpt) )
      kpq2kPrVec = 0
      do iqpt = 1, qpts%nkpt
        do ikpt = 1, kpts%nkpt
          kpqNomin = 0
          do idir = 1, 3
            kpqNomin(idir) = nint( mod( kpts%bk(idir, ikpt) + qpts%bk(idir, iqpt), 1. ) * lcm )
            !kpqNomin(idir) = nint(lcm*kpts%bk(idir, ikpt) + lcm*qpts%bk(idir, iqpt))
            ! Is in 1st Brillouin zone
            if (abs(real(kpqNomin(idir)) / real(lcm) - (kpts%bk(idir, ikpt) + qpts%bk(idir, iqpt))) < 1e-5) then
            !if (kpqNomin(idir).lt.int(lcm)) then
              kpq2kPrVec(idir, ikpt, iqpt) = 0
            ! Has to be backfolded
            else
              kpq2kPrVec(idir, ikpt, iqpt) = -1
            end if
          end do
          if(.false.) then
            write(1005, '(i5, i5, 3(i5))') iqpt, ikpt, kpq2kPrVec(:, ikpt, iqpt)
          end if
          mapKpq2K(ikpt, iqpt) = mapK2Ind( kpqNomin(1), kpqNomin(2), kpqNomin(3) )
          !mapKpq2K(ikpt, iqpt) = mapK2Ind( mod(kpqNomin(1),int(lcm)), mod(kpqNomin(2),int(lcm)), mod(kpqNomin(3),int(lcm)) )
        end do
      end do

      ! For this found k-vector equals to k + q, determine the index and fill up array which connects the kpts indices of the k and q
      ! qpoint with the index of the k-vector equals to k + q.
      nkptShift = 0
      matchFound = .false.
      do iqpt = qpts%nkpt + 1, qpts%nkpt + nrAddQs
        do ikpt = 1, kpts%nkpt
          kpqTemp(:) = modulo1r( kpts%bk(:, ikpt) + qpts%bk(:, iqpt) )
          do ikptSh = 1, addQnkptis(iqpt - qpts%nkpt)
            if ( norm2( kpts%bk(:, kpts%nkpt + ikptSh + nkptShift) - kpqTemp(:) ) < 1e-7 ) then
              mapKpq2K(ikpt, iqpt) = kpts%nkpt + nkptShift + ikptSh
              matchFound = .true.
              exit
            end if
          end do
          if ( .not.matchFound ) then
            write (errorMessage, '(a,1x,i3)') 'no match for k+q-point', kpqTemp
            call juDFT_error( errorMessage, calledby='createkqMapArrays', hint='Please check k-points, q-points and mapKpq2K array!' )
          else
            matchFound = .false.
          end if
        end do
        nkptShift = nkptShift + addQnkptis(iqpt - qpts%nkpt)
      end do

      if (.false.) then
        ! Finds out which k' results when sending k to -k and then backfolding into 1st Brillouin zone.
        do ikpt = 1, kpts%nkpt
          do idir = 1, 3
            if ( kpts%bk(idir, ikpt)==0 ) then
              kpqTemp(idir) = kpts%bk(idir, ikpt)
            else
              kpqTemp(idir) = -kpts%bk(idir, ikpt) + 1
            end if
          end do ! idir
          mapK2mK(ikpt) = mapK2Ind(int(kpqTemp(1) * lcm), int(kpqTemp(2) * lcm), int(kpqTemp(3) * lcm))
        end do ! ikpt
      end if

    end subroutine createkqMapArrays

    function modulo1r(kpoint)
      implicit none
      real(8), intent(in) :: kpoint(3)
      real(8)             :: modulo1r(3)
      integer             :: i

      modulo1r = modulo (kpoint , 1d0)

      do i = 1,3
        if(abs(1-abs(modulo1r(i))).lt.1d-13) modulo1r(i) = 0d0
      enddo
    end function modulo1r

    function kgv(iarr,n)
      implicit none
      integer              :: kgv
      integer, intent(in)  :: n,iarr(n)
      logical              :: lprim(2:maxval(iarr))
      integer, allocatable :: prim(:),expo(:)
      integer              :: nprim,marr
      integer              :: i,j,ia,k
      ! Determine prime numbers
      marr  = maxval(iarr)
      lprim = .true.
      do i = 2,marr
        j = 2
        do while (i*j.le.marr)
          lprim(i*j) = .false.
          j          = j + 1
        enddo
      enddo
      nprim = count(lprim)
      allocate ( prim(nprim),expo(nprim) )
      j = 0
      do i = 2,marr
        if(lprim(i)) then
          j       = j + 1
          prim(j) = i
        endif
      enddo
      ! Determine least common multiple
      expo = 0
      do i = 1,n
        ia = iarr(i)
        if(ia.eq.0) cycle
        do j = 1,nprim
          k = 0
          do while(ia/prim(j)*prim(j).eq.ia)
            k  = k + 1
            ia = ia / prim(j)
          enddo
          expo(j) = max(expo(j),k)
        enddo
      enddo
      kgv = 1
      do j = 1,nprim
        kgv = kgv * prim(j)**expo(j)
      enddo
      deallocate ( prim,expo )
    end function kgv
END MODULE m_dfpt
