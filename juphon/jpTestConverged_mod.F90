module m_jpTestConverged

  implicit none

  contains

    subroutine TestConvergedQuantities( atoms, sym, lathar, dimens, stars, input, cell, kpts, qpts, enpara, usdus, Veff0, results, ngdp, ngdp2km, logUnit, memd_atom, numAddQs, gdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, rho0MT, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, kveclo, mapKpq2K, ne, eig, gbas, ilst, z, nv, El, nRadFun, iloTable, nobd, kpq2kPrVec, gdp2Ind, rho0IRpw, rho0MTsh, vEff0MTsh, vEff0IRpwUw, noPotsCon )

#include "cppmacro.h"

      use m_types

      implicit none

      type(t_atoms),                  intent(in) :: atoms
      type(t_sym),                    intent(in) :: sym
      type(t_sphhar),                 intent(in) :: lathar
      type(t_dimension),              intent(in) :: dimens
      type(t_stars),                  intent(in) :: stars
      type(t_input),                  intent(in) :: input
      type(t_cell),                   intent(in) :: cell
      type(t_kpts),                   intent(in) :: kpts
      type(t_kpts),                   intent(in) :: qpts
      type(t_enpara),                 intent(in) :: enpara
      type(t_usdus),                  intent(in) :: usdus
      type(t_potential),              intent(in) :: Veff0
      type(t_results),                intent(in) :: results

      ! Scalar parameter
      integer,                        intent(in) :: ngdp
      integer,                        intent(in) :: ngdp2km
      integer,                        intent(in) :: logUnit
      integer,                        intent(in) :: memd_atom
      integer,                        intent(in) :: numAddQs
      integer,                        intent(in) :: noPotsCon

      !Array parameter
      integer,                        intent(in) :: gdp(:, :)
      integer,                        intent(in) :: mlh_atom(:,0:,:)
      integer,                        intent(in) :: nmem_atom(0:, :)
      complex,                        intent(in) :: clnu_atom(:,0:,:)
      complex,                        intent(in) :: rho0IR(:, :) !n3d !todo  with spin
      real,                           intent(in) :: rho0MT(:, 0:, :, :) ! todo with spin
      real,                           intent(in) :: rbas1(:, :, 0:, :, :)
      real,                           intent(in) :: rbas2(:, :, 0:, :, :)
      real,                           intent(in) :: uuilon(:, :)
      real,                           intent(in) :: duilon(:, :)
      real,                           intent(in) :: ulouilopn(:, :, :)
      integer,                        intent(in) :: ilo2p(:, :)
      integer,                        intent(in) :: kveclo(:,:)
      integer,                        intent(in) :: mapKpq2K(:, :) ! todo possible error here
      integer,                        intent(in) :: ne(:)
      real,                           intent(in) :: eig(:,:,:)
      integer,                        intent(in) :: gbas(:, :)
      integer,                        intent(in) :: ilst(:, :, :)
      MCOMPLEX,                       intent(in) :: z(:, :, :, :) ! Attention: see zBar
      integer,                        intent(in) :: nv(:, :)
      real,                           intent(in) :: El(:, 0:, :, :)
      integer,                        intent(in) :: nRadFun(0:, :)
      integer,                        intent(in) :: iloTable(:, 0:, :)
      integer,                        intent(in) :: nobd(:, :)
      integer,                        intent(in) :: kpq2kPrVec(:, :, :)
      integer,                        intent(in) :: gdp2Ind(:, :, :)
      complex,                        intent(in) :: rho0IRpw(:, :)
      complex,                        intent(in) :: rho0MTsh(:, :, :, :)
      complex,                        intent(in) :: vEff0MTsh(:, :, :, :)
      complex,                        intent(in) :: vEff0IRpwUw(:, :)

      !call TestRhoConverged(atoms, sym, qpts, lathar, dimens, stars, input, cell, vacuum,   ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho1IR, rho1MT)
      if (.false.) then
        call testGassignmWithQSternheimer( atoms, sym, cell, dimens, input, lathar, enpara, usdus, Veff0, stars, qpts, kpts, results, &
          & logUnit, ngdp, memd_atom, ngdp2km, numAddQs, gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, &
          & duilon, ulouilopn, ilo2p, kveclo, mapKpq2K, ne, eig, gbas, ilst, z, nv, El, nRadFun, iloTable, nobd, kpq2kPrVec, gdp2Ind, rho0IRpw, rho0MTsh, vEff0MTsh, vEff0IRpwUw, noPotsCon )
      end if

    end subroutine TestConvergedQuantities

    ! This routine tests whether the quantities coming from the Sternheimer routine are consistent with their G assignment of two
    ! related q-points.

    subroutine testGassignmWithQSternheimer(atoms, sym, cell, dimens, input, lathar, enpara, usdus, Veff0, stars, qpts, kpts, results, &
        & logUnit, ngdp, memd_atom, ngdp2km, numAddQs, gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, kveclo, mapKpq2K, ne, eig, gbas, ilst, z, nv, El, nRadFun, iloTable, nobd, kpq2kPrVec, gdp2Ind, rho0IRpw, rho0MTsh, vEff0MTsh , vEff0IRpwUw, noPotsCon)

#include "cppmacro.h"

      use m_types
      use m_jpSternheimer, only : solveSternheimerSCC, initSternheimerSCC
      use m_jpPotDensHelper, only : genPertPotDensGvecs

      implicit none

      ! Type parameters
      type(t_atoms),                  intent(in) :: atoms
      type(t_sym),                    intent(in) :: sym
      type(t_cell),                   intent(in) :: cell
      type(t_dimension),              intent(in) :: dimens
      type(t_input),                  intent(in) :: input
      type(t_sphhar),                 intent(in) :: lathar
      type(t_enpara),                 intent(in) :: enpara
      type(t_usdus),                  intent(in) :: usdus
      type(t_potential),              intent(in) :: Veff0
      type(t_stars),                  intent(in) :: stars
      type(t_kpts),                   intent(in) :: qpts
      type(t_kpts),                   intent(in) :: kpts
      type(t_results),                intent(in) :: results


      ! Scalar parameters
      integer,                        intent(in) :: logUnit
      integer,                        intent(in) :: ngdp
      integer,                        intent(in) :: memd_atom
      integer,                        intent(in) :: ngdp2km
      integer,                        intent(in) :: numAddQs
      integer,                        intent(in) :: noPotsCon

      ! Array parameters
      integer,                        intent(in) :: gdp(:, :)
      complex,                        intent(in) :: rho0IR(:, :) !n3d !todo  with spin
      real,                           intent(in) :: rho0MT(:, 0:, :, :) ! todo with spin
      integer,                        intent(in) :: mlh_atom(:,0:,:)
      integer,                        intent(in) :: nmem_atom(0:, :)
      complex,                        intent(in) :: clnu_atom(:,0:,:)
      real,                           intent(in) :: rbas1(:, :, 0:, :, :)
      real,                           intent(in) :: rbas2(:, :, 0:, :, :)
      real,                           intent(in) :: uuilon(:, :)
      real,                           intent(in) :: duilon(:, :)
      real,                           intent(in) :: ulouilopn(:, :, :)
      integer,                        intent(in) :: ilo2p(:, :)
      integer,                        intent(in) :: kveclo(:,:)
      integer,                        intent(in) :: mapKpq2K(:, :) ! todo possible error here
      integer,                        intent(in) :: ne(:)
      real,                           intent(in) :: eig(:,:,:)
      integer,                        intent(in) :: gbas(:, :)
      integer,                        intent(in) :: ilst(:, :, :)
      MCOMPLEX,                       intent(in) :: z(:, :, :, :) ! Attention: see zBar
      integer,                        intent(in) :: nv(:, :)
      real,                           intent(in) :: El(:, 0:, :, :)
      integer,                        intent(in) :: nRadFun(0:, :)
      integer,                        intent(in) :: iloTable(:, 0:, :)
      integer,                        intent(in) :: nobd(:, :)
      integer,                        intent(in) :: kpq2kPrVec(:, :, :)
      integer,                        intent(in) :: gdp2Ind(:, :, :)
      complex,                        intent(in) :: rho0IRpw(:, :)
      complex,                        intent(in) :: rho0MTsh(:, :, :, :)
      complex,                        intent(in) :: vEff0MTsh(:, :, :, :)
      complex,                        intent(in) :: vEff0IRpwUw(:, :)

      ! Type variables
      type(t_tlmplm)                             :: tdHS0

      ! Scalar variables
      integer                                    :: lcm
      integer                                    :: iqpt
      integer                                    :: subIqpt
      integer                                    :: idir
      integer                                    :: ngpqdp
      integer                                    :: ngpqdp2
      integer                                    :: ngpqdp2km
      logical                                    :: oneSternhCycle
      real                                       :: paPoX
      real                                       :: paPoY
      real                                       :: paPoZ
      integer                                    :: iatom
      integer                                    :: itype
      integer                                    :: ieqat
      logical                                    :: passed
      real                                       :: minQcomp
      integer                                    :: iG
      logical                                    :: l_dummy

      ! Array variables
      integer                                    :: gShift(3)
      real                                       :: qpoint(3)
      integer                                    :: gpqdp2iLim(2, 3)
      integer,           allocatable             :: qIndex(:, :, :)
      integer,           allocatable             :: gpqdp2(:, :)
    integer,             allocatable             :: gpqdp(:, :)
      integer,           allocatable             :: gdpIndex(:, :, :)
      integer,           allocatable             :: gpqdp2Ind(:, :, :)
      real,              allocatable             :: acoff(:)
      real,              allocatable             :: alpha(:)
      complex,           allocatable             :: qpwcG(:, :)
      complex,           allocatable             :: rho1MTCoreDispAt(:, :, :, :)
      complex,           allocatable             :: grVxcIRKern(:)
      real,              allocatable             :: dKernMTGPts(:, :, :)
      complex,           allocatable             :: grVeff0MT_init(:, :, :, :)
      complex,           allocatable             :: grVeff0MT_main(:, :, :, :)
      complex,           allocatable             :: grRho0MT(:, :, :, :)
      real,              allocatable             :: gausWts(:) ! gaussian weights belonging to gausPts
      complex,           allocatable             :: ylm(:, :)
      complex,           allocatable             :: grVext0IR_DM(:, :)
      complex,           allocatable             :: grVext0MT_DM(:, :, :, :)
      complex,           allocatable             :: grVeff0IR_DM(:, :)
      complex,           allocatable             :: grVeff0MT_DM(:, :, :, :)
      complex,           allocatable             :: grVeff0MT_DMhxc(:, :, :, :)
      complex,           allocatable             :: rho1IR(:, :, :)
      complex,           allocatable             :: rho1IR2(:, :, :)
      complex,           allocatable             :: rho1MT(:, :, :, :, :)
      complex,           allocatable             :: vExt1MT(:, :, :, :, :)
      complex,           allocatable             :: vExt1MTnoVol(:, :, :, :, :)
      complex,           allocatable             :: w_vEff1IR(:, :, :)
      complex,           allocatable             :: w_vEff1IR2(:, :, :)
      complex,           allocatable             :: vEff1MT(:, :, :, :, :)
      complex,           allocatable             :: vEff1MTnoVol(:, :, :, :, :)
      complex,           allocatable             :: vExt1IR_final(:, :, :)
      complex,           allocatable             :: vExt1IR_final2(:, :, :)
      complex,           allocatable             :: vHar1IR_final(:, :, :)
      complex,           allocatable             :: vCoul1IRtempNoVol(:, :)
      complex,           allocatable             :: vCoul1MTtempNoVol(:, :, :, :)
      complex,           allocatable             :: grVCoul0IR_DM_SF(:, :)
      complex,           allocatable             :: grVCoul0MT_DM_SF(:, :, :, :)
      complex,           allocatable             :: vHar1IR_final2(:, :, :)
      complex,           allocatable             :: vHar1MT_final(:, :, :, :, :)
      complex,           allocatable             :: grRho0IR(:, :)
      complex,           allocatable             :: dummy1(:, :, :, :, :),dummy2(:, :, :, :, :),dummy3(:, :, :, :),dummy4(:, :, :)
      complex,           allocatable             :: dummy5(:, :, :, :)
      integer                                    :: gdp2iLim(2, 3)

      oneSternhCycle = .false.
      paPoX = 1.
      paPoY = 1.
      paPoZ = 1.

      minQcomp = 2.
      do iqpt = 1, qpts%nkpt
        do idir = 1, 3
          if ( (qpts%bk(idir, iqpt) > 1e-12) .and. (qpts%bk(idir, iqpt) < minQcomp)) minQcomp = qpts%bk(idir, iqpt)
        end do ! idir
      end do ! iqpt
      lcm = int(1. / minQcomp)
      allocate(qIndex(0:lcm, 0:lcm, 0:lcm))
      qIndex = 0
      do iqpt = 1, qpts%nkpt
        qIndex(nint(qpts%bk(1, iqpt) * lcm), nint(qpts%bk(2, iqpt) * lcm), nint(qpts%bk(3, iqpt) * lcm)) = iqpt
      end do ! iqpt


      do iqpt = 40, 40!qpts%nkpt

        do idir = 1, 3
          if (qpts%bk(idir, iqpt) > 1e-12) then
            gShift(idir) = 1
          else
            gshift(idir) = 0
          end if
        end do ! idir

        subIqpt = iqpt
        qpoint(1:3) = qpts%bk(1:3, subIqpt)
        call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpoint, gpqdp, gpqdp2Ind, gpqdp2iLim )
!        do iG = 1, ngpqdp
!          write(5002, '(4i8)') iG, gpqdp(:, iG)
!        end do ! iG

        ! We determine all quantities not dependent on the q-points
        call initSternheimerSCC( atoms, sym, cell, stars, dimens, lathar, enpara, usdus, input, tdHS0, qpts, logUnit, ngdp, subiqpt,&
          & gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, rho0IRpw,         &
          & rho0MTsh, vEff0MTsh, qpwcG, rho1MTCoreDispAt, grVxcIRKern, dKernMTGPts, grVeff0MT_init, grVeff0MT_main, grRho0IR,      &
          & grRho0MT, gausWts, ylm, grVext0IR_DM, grVext0MT_DM, grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, grVeff0IR_DM, grVeff0MT_DM,grVeff0MT_DMhxc  )
        write(*, '(a)') 'Sternheimer self-consistency cycle initialized!'

        ! Solve Sternheimer self-consistency cycle and obtain first order densities and first order potentials for all atoms
        write(*, '(a,1x,i3,a)') 'Sternheimer self-consistency cycle for q-point', iqpt, '...'
        call solveSternheimerSCC( atoms, sym, stars, lathar, dimens, cell, enpara, usdus, input, kpts, qpts, results, usdus,      &
          & logUnit, ngdp,  rbas1, rbas2, kveclo, uuilon, duilon, ulouilopn, &
          & gdp, mapKpq2K, ne, eig, gbas, ilst, z, nv, El, nradFun, iloTable, nobd, ilo2p, gdp2Ind,    &
          & gdp2iLim, kpq2kPrVec, qpwcG, subiqpt, tdHS0, ylm, grRho0IR, grRho0MT, grVeff0MT_init, grVeff0MT_main, dKernMTGPts,            &
          & grVxcIRKern, rho1MTCoreDispAt, gausWts, rho1IR, rho1MT, vExt1MT, w_vEff1IR, vEff1MT, oneSternhCycle, ngpqdp, gpqdp,&
          & vExt1IR_final, vHar1IR_final, vHar1MT_final, rho1MT, vHar1MT_final, vHar1MT_final, vHar1MT_final, vHar1MT_final, &
          & vHar1MT_final, vHar1MT_final, rho0IRpw, rho0MTsh, vEff0IRpwUw, noPotsCon, vExt1MTnoVol, vEff1MTnoVol, dummy1, dummy2, dummy3, dummy4, dummy5, vCoul1IRtempNoVol, vCoul1MTtempNoVol)
        write(*, '(a)') '1st-order quantities determined!'

        deallocate(grVeff0MT_DM, grVeff0MT_DMhxc, grVeff0IR_DM, grVext0MT_DM, grVext0IR_DM, ylm, gausWts, grRho0MT, grRho0IR, grVeff0MT_main, grVeff0MT_init, dKernMTGPts, grVxcIRKern, rho1MTCoreDispAt, qpwcG, alpha, acoff, gpqdp2Ind)

        qPoint(1:3) = nint((-qpts%bk(1:3, iqpt) + gShift(1:3)) * lcm)
        subIqpt = qIndex(qPoint(1), qPoint(2), qPoint(3))
        write(*, *) 'subiqpt', subIqpt
        call genPertPotDensGvecs( stars, cell, input, ngpqdp2, ngpqdp2km, qpoint, gpqdp2, gpqdp2Ind, gpqdp2iLim )

        ! We determine all quantities not dependent on the q-points
        call initSternheimerSCC( atoms, sym, cell, stars, dimens, lathar, enpara, usdus, input, tdHS0, qpts, logUnit, ngdp, subiqpt,&
          & gdp, rho0IR, rho0MT, mlh_atom, nmem_atom, clnu_atom, rbas1, rbas2, uuilon, duilon, ulouilopn, ilo2p, rho0IRpw,         &
          & rho0MTsh, vEff0MTsh, qpwcG, rho1MTCoreDispAt, grVxcIRKern, dKernMTGPts, grVeff0MT_init, grVeff0MT_main, grRho0IR,      &
          & grRho0MT, gausWts, ylm, grVext0IR_DM, grVext0MT_DM, grVCoul0IR_DM_SF, grVCoul0MT_DM_SF, grVeff0IR_DM, grVeff0MT_DM, grVeff0MT_DMhxc )
        write(*, '(a)') 'Sternheimer self-consistency cycle initialized!'

        ! Solve Sternheimer self-consistency cycle and obtain first order densities and first order potentials for all atoms
        write(*, '(a,1x,i3,a)') 'Sternheimer self-consistency cycle for q-point', subiqpt, '...'
        call solveSternheimerSCC( atoms, sym, stars, lathar, dimens, cell, enpara, usdus, input, kpts, qpts, results, usdus,      &
          & logUnit, ngdp, rbas1, rbas2, kveclo, uuilon, duilon, ulouilopn, &
          & gdp, mapKpq2K, ne, eig, gbas, ilst, z, nv, El, nradFun, iloTable, nobd, ilo2p, gdp2Ind,     &
          & gdp2iLim, kpq2kPrVec, qpwcG, subIqpt, tdHS0, ylm, grRho0IR, grRho0MT, grVeff0MT_init, grVeff0MT_main, dKernMTGPts,            &
          & grVxcIRKern, rho1MTCoreDispAt, gausWts, rho1IR2, rho1MT, vExt1MT, w_vEff1IR2, vEff1MT, oneSternhCycle, ngpqdp2, gpqdp2,&
          & vExt1IR_final2, vHar1IR_final2, vHar1MT_final, rho1MT, vHar1MT_final, vHar1MT_final, vHar1MT_final, vHar1MT_final, vHar1MT_final, &
          & vHar1MT_final, rho0IRpw, rho0MTsh, vEff0IRpwUw, noPotsCon, vExt1MTnoVol, vEff1MTnoVol, dummy1, dummy2, dummy3, dummy4, dummy5, vCoul1IRtempNoVol, vCoul1MTtempNoVol)
        write(*, '(a)') '1st-order quantities determined!'

!        do iG = 1, ngpqdp2
!          write(5003, '(4i8)') iG, gpqdp2(:, iG)
!        end do ! iG
!        NOstopNO
        ! Test rho1
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            write(*, *) 'rho'
            call testGassignmWithQ( ngpqdp, ngpqdp2, gpqdp, gpqdp2, gShift, rho1IR(:, :, iatom), rho1IR2(:, :, iatom), passed )
            write(*, *) 'w_veff1'
            call testGassignmWithQ( ngpqdp, ngpqdp2, gpqdp, gpqdp2, gShift, w_vEff1IR(:, :, iatom), w_vEff1IR2(:, :, iatom), passed )
            write(*, *) 'vext1'
            call testGassignmWithQ( ngpqdp, ngpqdp2, gpqdp, gpqdp2, gShift, vExt1IR_final(:, :, iatom), vExt1IR_final2(:, :, iatom), passed )
            write(*, *) 'vhar1'
            call testGassignmWithQ( ngpqdp, ngpqdp2, gpqdp, gpqdp2, gShift, vHar1IR_final(:, :, iatom), vHar1IR_final2(:, :, iatom), passed )
          end do ! ieqat
        end do ! itype

        deallocate( acoff, alpha, qpwcG, rho1MTCoreDispAt, grVxcIRKern, dKernMTGPts, grVeff0MT_init, grVeff0MT_main, grRho0IR, grRho0MT, gausWts, ylm, grVext0IR_DM, grVext0MT_DM, grVeff0IR_DM, grVeff0MT_DM, rho1IR, rho1IR2, rho1MT, vExt1MT, vExt1IR_final, vExt1IR_final2, vHar1IR_final, vHar1IR_final2, vHar1MT_final )

      end do ! iqpt

    end subroutine testGassignmWithQSternheimer
    !todo now there is a old test from Sternheimer loop that also checks the symmetry for the muffin-tin for q = 0 0 1/4 and 0 0 3/4
    !if (.false.) then
    !  if ( iqpt == 2 ) then
    !    write(*, *) 'was here mixed'
    !    iDatom = 0
    !    do iDtype = 1, atoms%ntype
    !      do iDeqat = 1, atoms%neq(iDtype)
    !        iDatom = iDatom + 1
    !        do idir = 1, 3
    !          do iG = 1, ngpqdp
    !            write(2201, '(5i8,2f15.8)') idir, iG, gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG), rho1IRDS(iG, idir, iDatom)
    !            write(2203, '(5i8,2f15.8)') idir, iG, gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG), w_vExt1IR(iG, idir, iDatom)
    !            write(2205, '(5i8,2f15.8)') idir, iG, gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG), w_vEff1IR_final(iG, idir, iDatom)
    !          end do ! iG
    !        end do ! idir
    !      end do ! iDeqat
    !    end do ! iDtype
    !    iDatom = 0
    !    do iDtype = 1, atoms%ntype
    !      do iDeqat = 1, atoms%neq(iDtype)
    !        iDatom = iDatom + 1
    !        do idir = 1, 3
    !          iatom = 0
    !          do itype = 1, atoms%ntype
    !            do ieqat = 1, atoms%neq(itype)
    !              iatom = iatom + 1
    !              do oqn_l = 0, atoms%lmax(itype)
    !                do mqn_m = -oqn_l, oqn_l
    !                  lm = oqn_l * (oqn_l + 1) + mqn_m + 1! todo does everythink start at 1?
    !                  do imesh = 1, atoms%jri(itype)
    !                    write(2207, '(5(i8),2(f15.8))') iDatom, idir, iatom, lm, imesh, rho1MT(imesh, lm, iatom, idir, iDatom)
    !                    write(2209, '(5(i8),2(f15.8))') iDatom, idir, iatom, lm, imesh, vExt1MT(imesh, lm, iatom, idir, iDatom)
    !                    write(2211, '(5(i8),2(f15.8))') iDatom, idir, iatom, lm, imesh, vEff1MT_final(imesh, lm, iatom, idir, iDatom)
    !                  end do ! imesh
    !                end do ! mqn_m
    !              end do ! oqn_l
    !            end do ! ieqat
    !          end do ! itype
    !        end do ! idir
    !      end do ! iDeqat
    !    end do ! iDtype
    !    write(2186) gpqdp
    !  else if ( iqpt == 4 ) then
    !    write(*, *) 'and here'
    !    allocate(gpqdp2(3, ngpqdp))
    !    gpqdp2(:, :) = 0
    !    read(2186) gpqdp2(:, :)
    !    allocate(gdpIndex(minval(gpqdp(1, :)) : maxval(gpqdp(1, :)), minval(gpqdp(2, :)) : maxval(gpqdp(2, :)), &
    !                                                                                 & minval(gpqdp(3, :))-1 : maxval(gpqdp(3, :))))
    !    gdpIndex(:, :, :) = 0

    !    do iG = 1, ngpqdp
    !      gdpIndex(gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG)) = iG
    !    end do ! iG
    !    iDatom = 0
    !    do iDtype = 1, atoms%ntype
    !      do iDeqat = 1, atoms%neq(iDtype)
    !        iDatom = iDatom + 1
    !        do idir = 1, 3
    !          do iG = 1, ngpqdp
    !            iGvar = gdpIndex(-gpqdp2(1, iG), -gpqdp2(2, iG), -gpqdp2(3, iG)-1)
    !            ! If I want to access 0.75 from 0.25, I have to go to -0.25 which gives a -G and then go + G_T to 0.75 so that the Gs
    !            ! have to be subtracted by - G_T. If I want to make the new G-Vectors equal, I first go back to -0.25 and than to 0.25
    !            ! so that it gives me a -(Gz + 1) = - Gz - 1 in analogy to the previous test.
    !            write(2202, '(5i8,2f15.8)') idir, iG, -gpqdp(1, iGvar), -gpqdp(2, iGvar), -gpqdp(3, iGvar)-1, &
    !                                                                                          & conjg(rho1IRDS(iGvar, idir, iDatom))
    !            write(2204, '(5i8,2f15.8)') idir, iG, -gpqdp(1, iGvar), -gpqdp(2, iGvar), -gpqdp(3, iGvar)-1, &
    !                                                                                         & conjg(w_vExt1IR(iGvar, idir, iDatom))
    !            write(2206, '(5i8,2f15.8)') idir, iG, -gpqdp(1, iGvar), -gpqdp(2, iGvar), -gpqdp(3, iGvar)-1, &
    !                                                                                   & conjg(w_vEff1IR_final(iGvar, idir, iDatom))
    !          end do ! iG
    !        end do ! idir
    !      end do ! iDeqat
    !    end do ! iDtype
    !    iDatom = 0
    !    do iDtype = 1, atoms%ntype
    !      do iDeqat = 1, atoms%neq(iDtype)
    !        iDatom = iDatom + 1
    !        do idir = 1, 3
    !          iatom = 0
    !          do itype = 1, atoms%ntype
    !            do ieqat = 1, atoms%neq(itype)
    !              iatom = iatom + 1
    !              do oqn_l = 0, atoms%lmax(itype)
    !                do mqn_m = -oqn_l, oqn_l
    !                  lm = oqn_l * (oqn_l + 1) + mqn_m + 1! todo does everythink start at 1?
    !                  do imesh = 1, atoms%jri(itype)
    !                    write(2208, '(5i8,2f15.8)') iDatom, idir, iatom, lm, imesh, (-1)**(oqn_l + 1) * &
    !                                                                                        & rho1MT(imesh, lm, iatom, idir, iDatom)
    !                    write(2210, '(5(i8),2(f15.8))') iDatom, idir, iatom, lm, imesh, (-1)**(oqn_l + 1) * &
    !                                                                                       & vExt1MT(imesh, lm, iatom, idir, iDatom)
    !                    write(2212, '(5(i8),2(f15.8))') iDatom, idir, iatom, lm, imesh, (-1)**(oqn_l + 1) * &
    !                                                                                 & vEff1MT_final(imesh, lm, iatom, idir, iDatom)
    !                  end do ! imesh
    !                end do ! mqn_m
    !              end do ! oqn_l
    !            end do ! ieqat
    !          end do ! itype
    !        end do ! idir
    !      end do ! iDeqat
    !    end do ! iDtype
    !  end if ! q-point selection
    !end if ! end of maintenance

  ! This test is a fork of TestQLatticePeriodicity in jpTestPotential. It takes two quantities that are expanded in Gs and checks
  ! whether they fulfill the relation of related q-points where a q is made to a -q and then shifted by a vector back to the first
  ! Brillouin zone.
  subroutine testGassignmWithQ( ngpqdp, ngpqdp2, gpqdp, gpqdp2, gShift, f1IR, f2IR, passed )

    implicit none

    ! Scalar parameters
    integer,         intent(in)              :: ngpqdp
    integer,         intent(in)              :: ngpqdp2
    logical,         intent(out)             :: passed

    ! Array parameters
    integer,         intent(in)              :: gpqdp(:, :)
    integer,         intent(in)              :: gpqdp2(:, :)
    integer,         intent(in)              :: gShift(:)
    complex,         intent(in)              :: f1IR(:, :)
    complex,         intent(in)              :: f2IR(:, :)

    ! Scalar variables
    integer                                  :: idir
    integer                                  :: iG
    integer                                  :: iGvar

    ! Array variables
    integer,           allocatable              :: gdpIndex(:, :, :)

    passed = .true.

    ! Output of the reference quantity
    if (.true.) then
      do idir = 1, 3
        do iG = 1, ngpqdp
          write(4100, '(5i8,2f15.8)') idir, iG, gpqdp(1, iG), gpqdp(2, iG), gpqdp(3, iG), f1IR(iG, idir)
        end do
      end do
    end if

    ! Gset for -q projected back to the shifted first Brillouin zone (0 to 1)
    allocate(gdpIndex(minval(gpqdp2(1, :)) : maxval(gpqdp2(1, :)), minval(gpqdp2(2, :)) : maxval(gpqdp2(2, :)), minval(gpqdp2(3, :))-1 : maxval(gpqdp2(3, :))))
    gdpIndex(:, :, :) = -1

    ! Create mapping array from G-vector to G-vector index
    do iG = 1, ngpqdp2
      gdpIndex(gpqdp2(1, iG), gpqdp2(2, iG), gpqdp2(3, iG)) = iG
    end do ! iG

    do idir = 1, 3
      do iG = 1, ngpqdp2
        ! Find index in second G-set that corresponds to a certain G in the first G-set
        iGvar = gdpIndex(-(gpqdp(1, iG) + gShift(1)), -(gpqdp(2, iG) + gShift(2)), -(gpqdp(3, iG) + gShift(3)))
        if (iGvar == -1) then
          write(*, *) 'Error in Gset at -q - Gf'
        end if
        ! Note: if we start from q, we find the partner in the first Brillouin zone if we additionaly shift -q by the shiftvector Gf
        if ( (abs(conjg(f2IR(iGvar, idir)) - f1IR(iG, idir)) > 1e-8) ) then
          passed = .false.
        end if

        if (.true.) then
          write(4101, '(5i8,2f15.8)') idir, iG, gpqdp2(1, iGvar), gpqdp2(2, iGvar), gpqdp2(3, iGvar), f2IR(iGvar, idir)
          write(4102, '(5i8,2f15.8)') idir, iG, -(gpqdp2(1, iGvar) + gShift(1)), -(gpqdp2(2, iGvar) + gShift(2)), -(gpqdp2(3, iGvar) +gShift(3)), conjg(f2IR(iGvar, idir))
        end if

      end do ! iG
    end do ! idir
    write(*, *) 'fancy', passed

  end subroutine testGassignmWithQ

!  subroutine TestRhoConverged(atoms, sym, qpts, lathar, dimens, stars, input, cell, vacuum,   ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho1IR, rho1MT)
!
!    use m_jpPotDensHelper, only : genPertPotDensGvecs
!    use m_jpPlotObservables, only : plotPathUCScalar
!    use m_jpPlotObservables, only : plotVestPotDens1, plot1stDensVar
!    use mod_juPhonUtils, only : Fopen, Fclose
!    use m_types
!    use m_loddop
!    use mod_juPhonUtils, only : convertStar2G
!
!    implicit none
!
!    ! Type parameters
!    type(t_atoms),                  intent(in) :: atoms
!    type(t_sym),                    intent(in) :: sym
!    type(t_sphhar),                 intent(in) :: lathar
!    type(t_dimension),              intent(in) :: dimens
!    type(t_stars),                  intent(in) :: stars
!    type(t_input),                  intent(in) :: input
!    type(t_cell),                   intent(in) :: cell
!    type(t_vacuum),                 intent(in) :: vacuum
!     
!    type(t_kpts),                   intent(in) :: qpts
!
!    ! Scalar parameter
!    integer,                        intent(in) :: ngdp
!
!    !Array parameter
!    integer,                        intent(in) :: gdp(:, :)
!    integer,                        intent(in) :: mlh_atom(:,0:,:)
!    integer,                        intent(in) :: nmem_atom(0:, :)
!    complex,                        intent(in) :: clnu_atom(:,0:,:)
!    complex,                        intent(in) :: rho1IR(:, :, :)
!    complex,                        intent(in) :: rho1MT(:, :, :, :, :)
!
!    ! Scalar variables
!    integer                                    :: iqpt
!    integer(4)                                 :: iunit = 205        ! unit number for density
!    integer                                    :: iter               ! irrelevant for juPhon
!    integer                                    :: ispin              ! loop variable
!    integer                                    :: itype              ! loop variable
!    integer                                    :: ilh                ! loop variable
!    integer                                    :: imesh              ! loop variable
!    integer                                    :: istar              ! loop variable
!    integer                                    :: iatom
!    integer                                    :: idir
!    integer                                    :: ieqat
!    real                                       :: paPoXCo
!    real                                       :: paPoYCo
!    real                                       :: paPoZCo
!    integer                                    :: ptsym
!    integer                                    :: oqn_l
!    integer                                    :: imem
!    integer                                    :: mqn_m
!    integer                                    :: lm
!    integer                                    :: counter
!    logical                                    :: lmaxp1
!    integer                                    :: nrDX
!    integer                                    :: ngpqdp
!    integer                                    :: ngpqdp2km
!    real                                       :: dispDelta
!
!    ! Array variables
!    character(len=17)                          :: plotname
!    character(len=14)                          :: filenameRD !Reference density
!    character(len=14)                          :: filenameDD !Displaced density
!    character(len=2)                           :: qStrApp
!    character(8)                               :: dop, iop, name(10) ! irrelevant for juPhon calculation
!    complex,           allocatable             :: rho0IRrd(:,:)
!    real,              allocatable             :: rho0MTrd(:,:,:,:)
!    complex,           allocatable             :: rho0IRdd(:,:)
!    real,              allocatable             :: rho0MTdd(:,:,:,:)
!    complex,           allocatable             :: rho1DiffQuotIRst(:)
!    complex,           allocatable             :: rho1DiffQuotIRpw(:)
!    real,              allocatable             :: rho1DiffQuotMTlh(:, :, :, :)
!    complex,           allocatable             :: rho1DiffQuotMTsh(:, :, :, :)
!    real,              allocatable             :: nz(:,:,:)          ! irrelevant for juPhon calculation
!    complex,           allocatable             :: nzxy(:,:,:,:)      ! irrelevant for juPhon calcualtion
!    integer,           allocatable             :: gpqdp(:, :)
!    integer,           allocatable             :: gpqdp2Ind(:, :, :)
!    integer                                    :: gpqdp2iLim(2, 3)
!
!
!    ! Set up control parameters
!    !--------------------------
!    iqpt = 2
!    if (iqpt < 10) then
!      write(qStrApp, '(a1,i1)') '0', iqpt
!    else
!      ! qStrAppendix
!      write(qStrApp, '(i2)') iqpt
!    end if
!    paPoXCo = 1.
!    paPoYCo = 1.
!    paPoZCo = 1.
!
!    dispDelta = 0.0001
!
!    call genPertPotDensGvecs( stars, cell, input, ngpqdp, ngpqdp2km, qpts%bk(1:3, iqpt), gpqdp, gpqdp2Ind, gpqdp2iLim )
!
!    ! Load densities from FLEUR
!    !--------------------------
!    allocate( rho0IRrd(stars%ng3,input%jspins), rho0IRdd(stars%ng3,input%jspins), nzxy(vacuum%nmzxyd, stars%ng2-1,2, input%jspins), nz(vacuum%nmzd, 2, input%jspins), &
!      & rho0MTrd(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins), rho0MTdd(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
!
!    ! Load density for neutral supercell
!    write(filenameRD, '(a9)') 'cdn1SCq01'
!    call Fopen(iunit, name=filenameRD, status='old', action='read', form='unformatted')
!
!    rewind (iunit)
!
!    write(*, *) 'foo'
!    call loddop( input%jspins, stars%ng3, stars%ng2, vacuum%nmzxyd, vacuum%nmzd, atoms%jmtd, lathar%nlhd, atoms%ntype,           &
!      & input%jspins, stars%ng3, stars%ng2, vacuum%nvac, atoms%ntype, sym%invs, sym%invs2, input%film, lathar%nlh, atoms%jri,    &
!      & lathar%ntypsd, atoms%ntypsy, iunit, atoms%natd, atoms%neq, iop, dop, iter, rho0MTrd, rho0IRrd, nz, nzxy, name )
!
!    call Fclose(iunit)
!
!    ! Load density for displaced supercell
!    write(filenameDD, '(a7,a2)') 'cdn1SCq', qStrApp
!    write(*, *) filenameDD
!    call Fopen(iunit, name=filenameRD, status='old', action='read', form='unformatted')
!
!    rewind (iunit)
!
!    call loddop( input%jspins, stars%ng3, stars%ng2, vacuum%nmzxyd, vacuum%nmzd, atoms%jmtd, lathar%nlhd, atoms%ntype,           &
!      & input%jspins, stars%ng3, stars%ng2, vacuum%nvac, atoms%ntype, sym%invs, sym%invs2, input%film, lathar%nlh, atoms%jri,    &
!      & lathar%ntypsd, atoms%ntypsy, iunit, atoms%natd, atoms%neq, iop, dop, iter, rho0MTrd, rho0IRrd, nz, nzxy, name )
!
!    call Fclose(iunit)
!
!    ! Divide out factor r^2 so that rho0MT is only rho0MT and not r^2 rho0MT as it is done in Fleur. The factor sqrt(4pi) for the
!    ! zeroth component was already dividied out when constructing rho0MT in cdnmt routine in Fleur.
!    do ispin = 1, input%jspins
!      do itype = 1, atoms%ntype
!        do ilh = 0, lathar%nlhd
!          do imesh = 1, atoms%jri(itype)
!            rho0MTrd(imesh, ilh, itype, ispin) = rho0MTrd(imesh, ilh, itype, ispin) / atoms%rmsh(imesh, itype)                          &
!              & / atoms%rmsh(imesh, itype)
!            rho0MTdd(imesh, ilh, itype, ispin) = rho0MTdd(imesh, ilh, itype, ispin) / atoms%rmsh(imesh, itype)                          &
!              & / atoms%rmsh(imesh, itype)
!          end do
!        end do
!      end do
!    end do

    ! Calculate differential quotient
    !--------------------------------

!    !todo allocate the quotient quantities
!    allocate( rho1DiffQuotIRst(stars%n3d), rho1DiffQuotMTlh(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
!    do ispin = 1, input%jspins
!      do iStar = 1, stars%n3d
!        rho1DiffQuotIRst(iStar) = (rho0IRrd(iStar, 1) - rho0IRdd(iStar, 1)) / dispDelta
!      end do
!    end do ! ispin
!
!    do ispin = 1, input%jspins
!      do itype = 1, atoms%ntype
!        do ilh = 0, lathar%nlhd
!          do imesh = 1, atoms%jri(itype)
!            rho1DiffQuotMTlh(imesh, ilh, itype, ispin) = (rho0MTrd(imesh, ilh, itype, ispin) - rho0MTdd(imesh, ilh, itype, ispin)) / dispDelta
!          end do
!        end do
!      end do
!    end do

    ! Convert from star to plane waves and from lattice harmonics to spherical harmonics
    !-----------------------------------------------------------------------------------
!    allocate(rho1DIFFQuotIRpw(ngdp), rho1DIFFQuotMTSH(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3))
!    rho1DIFFQuotIRpw = cmplx(0., 0.)
!    call convertStar2G(rho1DiffQuotIRst(:), rho1DIFFQuotIRpw(:), stars, ngdp, gdp)
!
!    rho1DIFFQuotMTSH = cmplx(0., 0.)
!    iatom = 0
!    do itype = 1, atoms%ntype
!      do ieqat = 1, atoms%neq(itype)
!        iatom = iatom + 1
!        ptsym = atoms%ntypsy(iatom)
!        do ilh = 0, lathar%nlh(ptsym)
!          oqn_l = lathar%llh(ilh, ptsym)
!          do imem = 1, nmem_atom(ilh, iatom)
!            mqn_m = mlh_atom(imem, ilh, iatom)
!            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!            do imesh = 1, atoms%jri(itype)
!              rho1DIFFQuotMTSH(imesh, lm, iatom, 1) = rho1DIFFQuotMTSH(imesh, lm, iatom, 1) + rho1DiffQuotMTlh(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
!            end do ! imesh
!          end do ! imem
!        end do ! ilh
!      end do
!    end do
!    NOstopNO

    ! Plot rho1 the differential quotient and the difference along a chosen path
    !---------------------------------------------------------------------------
!    iatom = 0
!    do itype = 1, atoms%ntype
!      do ieqat = 1, atoms%neq(itype)
!        iatom = iatom + 1
!      !todo is ngdp, gdp here okkay?
!        counter = 1
!        iqpt = 1
!        lmaxp1 = .false.
!        nrDX = 600
!        call plotPathUCScalar(cell, atoms, rho1DIFFQuotIRpw, ngdp, gdp, qpts%bk(:, iqpt), rho1DIFFQuotMTSH(:, :, :, 1), paPoXCo, paPoYCo, paPoZCo, counter, iqpt, lmaxp1, nrDX )
!        !call plotPathUCScalar(cell, atoms, rho1IR(:, 3, iatom), ngpqdp, gpqdp, qpts%bk(:, iqpt), rho1MT(:, :, :, 3, 1), paPoXCo, paPoYCo, paPoZCo, counter, iqpt, lmaxp1, nrDX )
!
!        iqpt = 2
!        !Note: rho0MTrd is here only a dummy variable, because the plot1stDensVar routine does still want to have a rho
!        call plot1stDensVar(sym, cell, input, lathar, atoms, rho1IR(:, :, iatom), ngpqdp, gpqdp, mlh_atom, nmem_atom, clnu_atom, rho0MTrd, qpts%bk(:, iqpt), rho1MT(:, :, :, :, iatom), paPoXCo, paPoYCo, paPoZCo, iatom, counter, iqpt) !todo remove counter!
!      !  do idir = 1, 3
!      !    plotname = ''
!      !    write(plotname, '(a,i1,a,i1,a)') 'plotVesta', idir, 'at',iatom, '.xsf'
!      !    ! Create file of density variation for Vesta
!!     !     call plotVestPotDens1(dimens, stars, lathar, atoms, input, sym, cell, gdp, ngdp, qpts%bk(:, iqpt), rho1IRDS(:, idir, iatom), rho1MT(:, :, :, idir, iatom), plotname, iqpt)
!      !  end do
!      end do
!    end do

    ! Expand the difference of rho1 and the difference quotient on a 3D mesh and check that that it vanishes. One might also plot it
    !-------------------------------------------------------------------------------------------------------------------------------
    ! for visualization in VESTA
    !---------------------------
!  end subroutine TestRhoConverged

end module m_jpTestConverged
