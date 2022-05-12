!todo we have a multiple definition of calcsurfintmtdynmat in mt and ir check this
module m_jpTestDynMatDeprecated

  implicit none

  contains


  subroutine TestDynMatDeprecated( atoms, enpara, lathar, sym, cell, kpts, dimens, usdus, input, results, qpts, stars, Veff0,  &
      & logUnit, ngdp, memd_atom, GbasVec, gdp2iLim, kpq2kPrVec, gdp2Ind, gdp, mapKpq2K,    &
      & rho0IR, rho0MT, nv, mapGbas, z0, kveclo, nRadFun, rbas1, rbas2, iloTable, El, eig, nobd, clnu_atom, nmem_atom, mlh_atom,   &
      & uuilon, duilon, ulouilopn, ilo2p,   vacuum, ne )

#include "cppmacro.h"
    use m_types

    implicit none

    ! Scalar parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_enpara),                 intent(in) :: enpara
    type(t_sphhar),                 intent(in) :: lathar
    type(t_sym),                    intent(in) :: sym
    type(t_cell),                   intent(in) :: cell
    type(t_kpts),                   intent(in) :: kpts
    type(t_dimension),              intent(in) :: dimens
    type(t_usdus),                  intent(in) :: usdus
    type(t_input),                  intent(in) :: input
    type(t_results),                intent(in) :: results
    type(t_kpts),                   intent(in) :: qpts
    type(t_stars),                  intent(in) :: stars
    type(t_potential),              intent(in) :: Veff0
    type(t_vacuum),                 intent(in) :: vacuum
     

    ! Scalar parameters
    integer,                        intent(in) :: logUnit
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom

    ! Array parameters
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: gdp2iLim(2, 3)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    integer,                        intent(in) :: gdp2Ind(:, :, :)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: mapKpq2K(:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    complex,                        intent(in) :: rho0IR(:,:)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    MCOMPLEX,                       intent(in) :: z0(:,:,:,:)
    integer,                        intent(in) :: kveclo(:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    real,                           intent(in) :: eig(:, :, :)
    integer,                        intent(in) :: nobd(:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    real,                           intent(in) :: uuilon(:, :)
    real,                           intent(in) :: duilon(:, :)
    real,                           intent(in) :: ulouilopn(:, :, :)
    integer,                        intent(in) :: ilo2p(:, :)
    integer,                        intent(in) :: ne(:)

    if ( .false. ) then
      call testGoldsteinModes( atoms, dimens, kpts, cell, input, results, sym, usdus, lathar, qpts, stars, vacuum,   logUnit, ngdp,          &
        & memd_atom, nv, GbasVec, gdp, mapGbas, nobd, kveclo, mapKpq2K, z0, gdp2Ind, gdp2iLim, kpq2kPrVec, rbas1, rbas2, ilo2p,    &
        & clnu_atom, nmem_atom, mlh_atom, rho0MT, rho0IR )
    end if
    !todo put here all routines which which do not need any Sternheimer quantities

    if (.false.) then
    ! This routine is to test the value and derivative of the first basis variation. It turns out that the first basis variation vanishes
    ! within numerical accuracy. However, the gradient derivative does not. For k = G = 0 the matching coefficients only have a
    ! contribution within the s-channel (l = 0). It is therefore very simple to track certain contributions. The gradient derivative of the
    ! first basis variation collapses to the radial derivative term times Gaunt coefficient and constant prefactor. However, latter lies
    ! in the order of 1e-1. Therefore we need additional surface terms as Aaron has analyzed probably for every system.
    ! Needs to be made more efficient later!!!
    ! Only 64 k-points need over an hour.
      call plotBasVar1( atoms, sym, cell, kpts, dimens, usdus, GbasVec, nv, mapGbas, kveclo, nRadFun, rbas1, rbas2, iloTable )
    end if

    if (.false.) then
      write(*, *) 'this test activated'
      call TesttensorGradMatElems( atoms, dimens, stars, lathar, kpts, qpts, sym, cell, usdus, results, Veff0, input, ngdp, logUnit,       &
      & rho0IR, gdp, rho0MT, nRadFun, rbas1, rbas2, nmem_atom, mlh_atom, clnu_atom, El, nv, mapGbas, GbasVec, kveclo,   &
      & nobd, z0, iloTable, eig, kpq2kPrVec, memd_atom )
    end if

    if (.false.) then
      call TestTensorGradProdKet( atoms, lathar, stars, dimens, cell, sym, kpts, usdus, results, input, Veff0, qpts, ngdp, memd_atom, rho0IR, rho0MT, mlh_atom, nmem_atom,      &
        & clnu_atom, gdp, nRadFun, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nv, mapGbas, gBasVec, kveclo, nobd, z0, iloTable, kpq2kPrVec, El, eig )
    end if

    if (.false.) then
      call testGradBraCorrection( atoms, dimens, lathar, Veff0, kpts, sym, cell, usdus, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nRadFun, El, mlh_atom, nmem_atom, clnu_atom, nv, mapGbas, gBasVec, kveclo, nobd, z0, iloTable )
    end if

    if (.false.) then
      ! Test routine tests for calculating double MT gradient One Channel of the gradient is calculated in a different way to check
      ! the general routine.
      call TestDobGrad(atoms, rbas1(:, :, :, :, 1), rbas2(:, :, :, :, 1), nRadFun)
    end if

    if (.false.) then
      call TestGrHepsGrt(atoms, dimens, lathar, kpts, cell, input, stars, Veff0, sym, usdus, results, ngdp, nRadFun, nobd, rbas1, rbas2, nv, mapGbas, gBasVec, z0, gdp, kpq2kPrVec, El, eig, mlh_atom, Veff0%vr(:, :, :, 1), nmem_atom, clnu_atom, kveclo, iloTable )
    end if

    if (.true.) then
      call Calc2ArgSFIntegrals( atoms, stars, cell, lathar, Veff0, ngdp, rho0IR, gdp, clnu_atom, nmem_atom, mlh_atom, rho0MT )
    end if
  end subroutine TestDynMatDeprecated

  ! This routine is to test the value and derivative of the first basis variation. It turns out that the first basis variation vanishes
  ! within numerical accuracy. However, the gradient derivative does not. For k = G = 0 the matching coefficients only have a
  ! contribution within the s-channel (l = 0). It is therefore very simple to track certain contributions. The gradient derivative of the
  ! first basis variation collapses to the radial derivative term times Gaunt coefficient and constant prefactor. However, latter lies
  ! in the order of 1e-1. Therefore we need additional surface terms as Aaron has analyzed probably for every system.
  ! Needs to be made more efficient later!!!
  ! Only 64 k-points need over an hour.
  subroutine plotBasVar1( atoms, sym, cell, kpts, dimens, usdus, GbasVec, nv, mapGbas, kveclo, nRadFun, rbas1, rbas2, iloTable )

    use m_types, only : t_atoms, t_sym, t_cell, t_kpts, t_dimension, t_usdus
    use m_od_types, only : od_inp, od_sym
    use m_abcof3
    use m_jpConstants, only : iu, fpi
    use m_jpPotDensHelper, only : grFlmpYlmPerType
    use m_ylm_old
    use mod_juPhonUtils, only : Derivative

    implicit none

    ! Type paramter
    type(t_atoms),                  intent(in) :: atoms
    type(t_sym),                    intent(in) :: sym
    type(t_cell),                   intent(in) :: cell
    type(t_kpts),                   intent(in) :: kpts
    type(t_dimension),              intent(in) :: dimens
    type(t_usdus),                  intent(in) :: usdus

    ! Array parameters
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: kveclo(:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:, :, 0:, :, :)
    real,                           intent(in) :: rbas2(:, :, 0:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)

    ! Local type variables
    type (od_inp)                              :: odi_type !Both types, this and the one beneath have to be taken from od_types module
    type (od_sym)                              :: ods_type

    ! Scalar variables
    integer                                    :: nTestDir
    integer                                    :: idir
    integer                                    :: itype
    integer                                    :: iatom
    integer                                    :: lm
    integer                                    :: irdPt
    integer                                    :: maxlmp
    integer                                    :: oqn_l
    integer                                    :: ieqat
    integer                                    :: nmat
    integer                                    :: ispin
    integer                                    :: mqn_m
    integer                                    :: ikpt
    integer                                    :: iBas
    integer                                    :: imesh
    integer                                    :: lm_pre
    integer                                    :: p
    integer                                    :: irel

    ! Array variables
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    complex,              allocatable             :: rbasLM1(:,:,:)
    complex,              allocatable             :: rbasLM2(:,:,:)
    complex,           allocatable             :: grRsYlm1R(:, :, :, :)
    complex,           allocatable             :: grRsYlm1I(:, :, :, :)
    complex,           allocatable             :: grgrRsYlm(:, :, :, :)
    integer,           allocatable             :: lmaxExtended(:)
    complex,           allocatable             :: basCorr(:, :, :, :)
    complex,           allocatable             :: basCorr1(:, :, :)
    complex,           allocatable             :: basCorr2(:, :, :)
    complex                                    :: iuKpG(3)

    real                                    :: test(atoms%jmtd)
    real, allocatable :: dus(:, :)
    real, allocatable :: duds(:, :)
    real              :: test1(atoms%jmtd)
    real              :: test2(atoms%jmtd)
    real :: aQ
    integer :: tP
    real :: aQ1
    integer :: tP1
    integer :: idir2
    integer, allocatable :: nRadFunHack(:, :)

    ! Needed to recycle the gradient routine dependent on p also for a p-independent quantity.
    allocate( nRadFunHack(atoms%lmaxd, atoms%ntype) )
    nRadFunHack = 1

    ! Test quantites to do some statistics
    tP = 0
    aQ = 0.
    tP1 = 0
    aQ1 = 0.

    ! determine random numbers on unit sphere
    !! Different random numbers for every call
    !nTestDir = 1
    !call random_seed ()

    !allocate( rPtsCart(3, nTestDir, atoms%ntype, 3) )
    !! these idir show the displacement directions not the spatial direction
    !do idir = 1, 3
    !  do itype = 1, atoms%ntype
    !    do irdPt = 1, nTestDir
    !      ! Create random number
    !      call random_number(rPtsCart(:, irdPt, itype, idir))
    !      rPtsCart(:, irdPt, itype, idir) = rPtsCart(:, irdPt, itype, idir) / norm2(rPtsCart(:, irdPt, itype, idir)) * atoms%rmt(itype)
    !    end do ! irdPt
    !  end do ! itype
    !end do ! idir
    !write(*, '(3(f15.8,1x))') rPtsCart(:, 1, 1, 1)

    allocate( a( dimens%nvd, 0:(atoms%lmaxd + 1)**2 - 1, atoms%nat), b(dimens%nvd, 0:(atoms%lmaxd + 1)**2 - 1, atoms%nat), &
      ! 1. index alo blo oder clo, 2. index Gbasisvektoren der los, 3. index number of los, 4. index atom index
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat))

    ! matching coefficient times radial basis function
    allocate( rbasLM1( atoms%jmtd, 2 + atoms%nlod, 0:(atoms%lmaxd + 1)**2 - 1 ) )
    allocate( rbasLM2( atoms%jmtd, 2 + atoms%nlod, 0:(atoms%lmaxd + 1)**2 - 1 ) )
    rbasLM1(:, :, :) = cmplx(0.0, 0.0)
    rbasLM2(:, :, :) = cmplx(0.0, 0.0)

    ! Results of gradient routines
    allocate( grRsYlm1R(atoms%jmtd, 2 + atoms%nlod, 3, 0:(atoms%lmaxd + 2)**2 - 1) )
    allocate( grRsYlm1I(atoms%jmtd, 2 + atoms%nlod, 3, 0:(atoms%lmaxd + 2)**2 - 1) )
    allocate( grgrRsYlm(atoms%jmtd, 2 + atoms%nlod, 3, 0:(atoms%lmaxd + 2)**2 - 1 ) )
    grRsYlm1R(:, :, :, :) = cmplx(0.0, 0.0)
    grRsYlm1I(:, :, :, :) = cmplx(0.0, 0.0)
    grgrRsYlm(:, :, :, :) = cmplx(0.0, 0.0)

    !we do not want rotated local systems
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1

    ! Help arrays for calculating the first basis variation
    allocate( basCorr( atoms%jmtd, 1, 0:(atoms%lmaxd + 1)**2 - 1, 3 ) )
    allocate( basCorr1(atoms%jmtd, 3, 0:(atoms%lmaxd + 1)**2 - 1))
    allocate( basCorr2(atoms%jmtd, 3, 0:(atoms%lmaxd + 1)**2 - 1))

    !do ikpt = 1, kpts%nkpt
    do ikpt =  1, 1!43,  43!kpts%nkpt
    !write(*, *) ikpt

      ! dimension of Hamiltonian
      nmat = nv(1, ikpt) + atoms%nlotot
      a = cmplx(0.0, 0.0)
      b = cmplx(0.0, 0.0)
      bascof_lo = cmplx(0.0, 0.0)

      ! Determine the matching coefficients for the basis functions, i.e. these are not the matching coefficients acof, bcof, ccof!
      !todo adjust ylm Norm!
      ! dus duds
      allocate(dus(0:atoms%lmaxd, atoms%ntype))
      allocate(duds(0:atoms%lmaxd, atoms%ntype))
      dus(:, :) = usdus%dus(:, :, 1)
      duds(:, :) = usdus%duds(:, :, 1)
      do oqn_l = 0, atoms%lmax(1)
        test1 = 0.
        test2 = 0.
        ! We determine the derivative by calculating the integral, because if we would use the Fleur value, this test would not work
        ! for numerical reasons. dus and duds are determined as a side-product of solving a differential equation.
        call Derivative( rbas1(:, 1, oqn_l, 1, 1) / atoms%rmsh(:, 1), 1, atoms, test1 )
        call Derivative( rbas1(:, 2, oqn_l, 1, 1) / atoms%rmsh(:, 1), 1, atoms, test2 )
        dus(oqn_l, 1) = test1(atoms%jri(1))
        duds(oqn_l, 1) = test2(atoms%jri(1))
      end do
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), GbasVec(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        & GbasVec(2, mapGbas(:nv(1, ikpt), ikpt, 1)), GbasVec(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        & usdus%us, dus, usdus%uds, duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi_type, ods_type, a, b, bascof_lo )


      do iBas = 1, 1! nv(1, ikpt)
    !  write(*, *) 'foo', gbasVec(:, mapGbas(iBas, ikpt, 1))
        iuKpG(:) = iu * matmul( cell%bmat, kpts%bk(:, ikpt) + GbasVec(:, mapGbas(iBas, ikpt, 1)))
    !    write(*, *) GbasVec(:, mapGbas(iBas, ikpt, 1))
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            rbasLM1(:, :, :) = cmplx(0.0, 0.0)
            rbasLM2(:, :, :) = cmplx(0.0, 0.0)
            grRsYlm1R = cmplx(0.0, 0.0)
            grRsYlm1I = cmplx(0.0, 0.0)
            grgrRsYlm = cmplx(0.0, 0.0)
            do oqn_l = 0, atoms%lmax(itype)
              lm_pre = oqn_l * (oqn_l + 1)
              do mqn_m = -oqn_l, oqn_l
                lm = lm_pre + mqn_m
                !manual test of abcof3
               ! if(lm == 0)  write(1346, *) a(iBas, lm, iatom), b(iBas, lm, iatom), sqrt(fpi) / sqrt(cell%omtil) * atoms%rmt(1)**2 / (-2)  * usdus%duds(0, 1, 1)
               ! if(lm == 0)  write(1346, *) a(iBas, lm, iatom), b(iBas, lm, iatom), sqrt(fpi) / sqrt(cell%omtil) * atoms%rmt(1)**2 / (2)  * usdus%dus(0, 1, 1)
    !           if (lm == 0) write(*, *) a(iBas, lm, iatom) * usdus%dus(0, 1, 1) + b(iBas, lm, iatom) * usdus%duds(0, 1, 1)
                do imesh = 1, atoms%jri(itype)
                ! we only match the large component, the small component is not matched.
                  rbasLM1(imesh, 1, lm) = a(iBas, lm, iatom) * iu**oqn_l * rbas1(imesh, 1, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)

            !      rbasLM2(imesh, 1, lm) = a(iBas, lm, iatom) * rbas2(imesh, 1, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
                end do ! imesh
                do imesh = 1, atoms%jri(itype)
                  rbasLM1(imesh, 2, lm) = b(iBas, lm, iatom) * iu**oqn_l * rbas1(imesh, 2, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            !      rbasLM2(imesh, 2, lm) = b(iBas, lm, iatom) * rbas2(imesh, 2, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
                end do ! imesh
        !        do p = 3, nRadFun(oqn_l, itype)
        !          do imesh = 1, atoms%jri(itype)
        !            rbasLM1(imesh, 1, lm) = rbasLM1(imesh, 1, lm) + bascof_lo(1, mqn_m, iBas, iloTable(p, oqn_l, itype), iatom) &
        !                                   & * rbas1(imesh, 1, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
        !            rbasLM1(imesh, 2, lm) = rbasLM1(imesh, 2, lm) + bascof_lo(2, mqn_m, iBas, iloTable(p, oqn_l, itype), iatom) &
        !                                   & * rbas1(imesh, 2, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
        !            rbasLM1(imesh, p, lm) = bascof_lo(3, mqn_m, iBas, iloTable(p, oqn_l, itype), iatom) &
        !                                   & * rbas1(imesh, p, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            !        rbasLM2(imesh, 1, lm) = bascof_lo(1, mqn_m, iBas, iloTable(p, oqn_l, itype), iatom) &
            !                               & * rbas2(imesh, 1, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            !        rbasLM2(imesh, 2, lm) = bascof_lo(2, mqn_m, iBas, iloTable(p, oqn_l, itype), iatom) &
            !                               & * rbas2(imesh, 2, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            !        rbasLM2(imesh, p, lm) = bascof_lo(3, mqn_m, iBas, iloTable(p, oqn_l, itype), iatom) &
            !                                      & * rbas2(imesh, p, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
               !   end do ! imesh
               ! end do ! p
              end do ! mqn_m
            end do ! oqn_l

            call grFlmpYlmPerType(atoms, itype, atoms%lmax(itype), nRadFun, real(rbasLM1(:, :, :)), grRsYlm1R)
            call grFlmpYlmPerType(atoms, itype, atoms%lmax(itype), nRadFun, aimag(rbasLM1(:, :, :)), grRsYlm1I)
            !call grFlmpYlmPerType(atoms, itype, atoms%lmax(itype), nRadFun, rbasLM2(:, :, :), grRsYlm2)

            ! radial part of the basis function
            basCorr = cmplx(0.0, 0.0)
            basCorr1 = cmplx(0.0, 0.0)
            basCorr2 = cmplx(0.0, 0.0)
            do oqn_l = 0, atoms%lmax(itype)
              lm_pre = oqn_l * (oqn_l + 1)
              do mqn_m = -oqn_l, oqn_l
                lm = lm_pre + mqn_m
                do idir = 1, 3
                  do imesh = 1, atoms%jri(itype)
                    do p = 1, nRadFun(oqn_l, itype)
!                      basCorr1(imesh, idir, lm) = basCorr1(imesh, idir, lm) + iuKpG(idir) * rbasLM1(imesh, p, lm) &
!                                                & - grRsYlm1(imesh, p, idir, lm)
                      basCorr1(imesh, idir, lm) = basCorr1(imesh, idir, lm) + iuKpG(idir) * rbasLM1(imesh, p, lm)
                      basCorr2(imesh, idir, lm) = basCorr2(imesh, idir, lm) - (grRsYlm1R(imesh, p, idir, lm) + iu * grRsYlm1I(imesh, p, idir, lm))
                    end do ! p
                    basCorr(imesh, 1, lm, idir) = basCorr(imesh, 1, lm, idir) + basCorr1(imesh, idir, lm) + basCorr2(imesh, idir, lm)
                  end do ! imesh
                  tP = tP + 1
                  aQ = aQ + abs(basCorr1(atoms%jri(1), idir, lm) + basCorr2(atoms%jri(1), idir, lm))
                end do ! idir
              end do ! mqn_m
            end do ! oqn_l
           ! do idir = 1, 3
           !   !write(*, '(2(f15.8))') basCorr1(atoms%jri(1), idir, 1) + basCorr2(atoms%jri(1), idir, 1)
           !   write(*, '(2(f15.8))') basCorr(atoms%jri(1), 1, 3, idir)
           ! end do
           ! NOstopNO
            grRsYlm1R = cmplx(0.0, 0.0)
            grRsYlm1I = cmplx(0.0, 0.0)
            !basCorr(atoms%jri(1) - 1, 1, :, 1:3) = cmplx(0.0, 0.0)
            write(1020, '(2(i8))', advance='no') 0, 0
            do imesh = 1, atoms%jri(itype) - 1
              write(1020, '(2(f15.8,4x))', advance='no') atoms%rmsh(imesh, 1)
              ! T-matrix produces an overlap between m = -1, and m= 1
              write(1021, '(7(f15.8))') atoms%rmsh(imesh, 1), basCorr(imesh, 1, 1, 1), basCorr(imesh, 1, 1, 2), basCorr(imesh, 1, 1, 3)
            end do
            write(1021, '(7(f15.8))') atoms%rmsh(atoms%jri(1), 1), basCorr(atoms%jri(1), 1, 1, 1), basCorr(atoms%jri(1), 1, 1, 2), basCorr(atoms%jri(1), 1, 1, 3)
            write(1020, '(2(f15.8,4x))', advance='yes') atoms%rmsh(atoms%jri(1), 1)
            do idir = 1, 3
              do oqn_l = 0, atoms%lmax(1)
                lm_pre = oqn_l * (oqn_l + 1)
                do mqn_m = -oqn_l, oqn_l
                  lm = lm_pre + mqn_m
                  write(1020, '(2(i8))', advance='no') idir, lm
                  do imesh = 1, atoms%jri(itype)
                    write(1020, '(2(f15.8))', advance='no') basCorr(imesh, 1, lm, idir)
                  end do ! imesh
                  write(1020, *)
                  !write(1020, '(2(i8),2(f15.8))') oqn_l, mqn_m, basCorr(atoms%jri(1), 1, lm, 1)
                end do !mqn_m
              end do ! oqn_l
            end do ! idir
            NOstopNO
            do idir = 1, 3
              call grFlmpYlmPerType(atoms, itype, atoms%lmax(itype), nRadFunHack, real(basCorr(:, :, :, idir)), grRsYlm1R)
              call grFlmpYlmPerType(atoms, itype, atoms%lmax(itype), nRadFunHack, aimag(basCorr(:, :, :, idir)), grRsYlm1I)
              do idir2 = 1, 3
                do oqn_l = 0, atoms%lmax(1)
                  lm_pre = oqn_l * (oqn_l + 1)
                  do mqn_m = -oqn_l, oqn_l
                    lm = lm_pre + mqn_m
                    tp1 = tp1 + 1
                    aQ1 = aQ1 + abs(grRsYlm1R(atoms%jri(1), 1, idir2, lm) + iu * grRsYlm1I(atoms%jri(1), 1, idir2, lm))
                  end do ! mqn_m
                end do ! oqn_l
              end do ! idir2
            end do ! idir
            grgrRsYlm = grRsYlm1R + iu * grRsYlm1I
            iatom = iatom + atoms%neq(itype)
          end do ! ieqat
        end do ! itype
      end do ! iBas
    end do ! ikpt
    write(*, *) 'test result:', aQ / tP
    write(*, *) 'test result2:', aQ1 / tP1
    do oqn_l = 0, atoms%lmax(1)
      lm_pre = oqn_l * (oqn_l + 1)
      do mqn_m = -oqn_l, oqn_l
        lm = lm_pre + mqn_m
       ! write(1020, '(2(i8),6(f15.8))') oqn_l, mqn_m, real(basCorr1(atoms%jri(1), 1, lm)), aimag(basCorr1(atoms%jri(1), 1, lm)), real(basCorr2(atoms%jri(1), 1, lm)), aimag(basCorr2(atoms%jri(1), 1, lm)), basCorr1(atoms%jri(1), 1, lm) + basCorr2(atoms%jri(1), 1, lm)
        write(1020, '(2(i8),4(f15.8))') oqn_l, mqn_m, basCorr1(atoms%jri(1), 1, lm) + basCorr2(atoms%jri(1), 1, lm), grgrRsYlm(atoms%jri(1), 1, 3, lm)
      end do
    end do
    NOstopNO

  end subroutine plotBasVar1

  ! Test whether for for q = 0 the frequencies are zero. Currently only IR HF works.
  subroutine testGoldsteinModes( atoms, dimens, kpts, cell, input, results, sym, usdus, lathar, qpts, stars, vacuum,   logUnit,&
      & ngdp, memd_atom, nv, gbas, gdp, mapGbas, nobd, kveclo, mapKpq2K, z0, gdp2Ind, gdp2iLim, kpq2kPrVec, rbas1, rbas2, ilo2p,   &
      & clnu_atom, nmem_atom, mlh_atom, rho0MT, rho0IR )

    use m_types
    use m_jpConstants,     only : iu
    use m_jpPotDensHelper, only : calcGrFinLH, warpIRPot, calcGrR2FinLH
    use m_jp2ndOrdQuant,   only : GenVext2, CalcIIEnerg2
    use m_jpDens1stVar, only : calcKdepValRho1MT, calcRho1IRValDS, multRadSolVzcRho1MT
    use mod_juPhonUtils, only : convertStar2G
    use m_jpProcessDynMat, only : DiagonalizeDynMat
    use m_jpVcoul1, only : GenVeff1
    use m_jpSetupDynMat, only : SetupDynMatHF
    use m_jpSternheimer, only : storeZ1nG
    use m_jpSetupDynMatHelper, only : ReadInz1

    use m_jpGrVeff0, only : GenGrVeff0
    use m_jpTestDynMatHF, only : TestPlotVext2

    implicit none

    ! Type parameter
    type(t_atoms),                 intent(in)  :: atoms
    type(t_dimension),             intent(in)  :: dimens
    type(t_kpts),                  intent(in)  :: kpts
    type(t_cell),                  intent(in)  :: cell
    type(t_input),                 intent(in)  :: input
    type(t_results),               intent(in)  :: results
    type(t_sym),                   intent(in)  :: sym
    type(t_usdus),                 intent(in)  :: usdus
    type(t_sphhar),                intent(in)  :: lathar
    type(t_kpts),                  intent(in)  :: qpts
    type(t_stars),                 intent(in)  :: stars
    type(t_vacuum),                intent(in)  :: vacuum
     

    ! Scalar parameter
    integer,                       intent(in)  :: logUnit
    integer,                       intent(in)  :: ngdp
    integer,                       intent(in)  :: memd_atom

    ! Array parameter
    integer,                       intent(in)  :: nv(:, :)
    integer,                       intent(in)  :: gbas(:, :)
    integer,                       intent(in)  :: gdp(:, :)
    integer,                       intent(in)  :: mapGbas(:, :, :)
    integer,                       intent(in)  :: nobd(:, :)
    integer,                       intent(in)  :: kveclo(:,:)
    integer,                       intent(in)  :: mapKpq2K(:, :)
    complex,                       intent(in)  :: z0(:, :, :, :)
    integer,                       intent(in)  :: gdp2Ind(:, :, :)
    integer,                       intent(in)  :: gdp2iLim(2, 3)
    integer,                       intent(in)  :: kpq2kPrVec(:, :, :)
    real,                          intent(in)  :: rbas1(:, :, 0:, :, :)
    real,                          intent(in)  :: rbas2(:, :, 0:, :, :)
    integer,                       intent(in)  :: ilo2p(:, :)
    complex,                       intent(in)  :: clnu_atom(:, 0:, :) !member, 0,lathar, nat
    integer,                       intent(in)  :: nmem_atom(0:, :) !0,lathar, nat
    integer,                       intent(in)  :: mlh_atom(:, 0:, :) ! memd_atom, 0:lathar, nat
    real,                          intent(in)  :: rho0MT(:, 0:, :, :)
    complex,                       intent(in)  :: rho0IR(:,:)


    ! Scalar variables
    integer                                   :: iDatom
    logical                                   :: calcEigenVec
    integer                                   :: iDtype
    integer                                   :: iDeqat
    real                                      :: invsAtNr
    integer                                   :: ikpt
    integer                                   :: idir
    integer                                   :: iBas
    integer                                   :: iband
    integer                                   :: iqpt
    integer                                   :: iatom
    integer                                   :: itype
    integer                                   :: ieqat
    integer                                   :: oqn_l
    integer                                   :: mqn_m
    integer                                   :: lm
    integer                                   :: imesh
    logical                                   :: harSw
    logical                                   :: extSw
    logical                                   :: xcSw
    logical                                   :: vExt1FullSw
    logical                                   :: grRhoTermSw
    integer                                   :: iG
    integer                                   :: idirC
    logical                                   :: testGoldstein = .true.


    ! Array variables
    complex,           allocatable            :: grRho0MT(:, :, :, :)
    complex,           allocatable            :: z1Special(:, :, :)
    complex,           allocatable            :: rho1IR(:, :, :)
    complex,           allocatable            :: rho1MT(:, :, :, :, :)
    complex,           allocatable            :: uu(:,:,:)
    complex,           allocatable            :: du(:,:,:)
    complex,           allocatable            :: dd(:,:,:)
    complex,           allocatable            :: ud(:,:,:)
    complex,           allocatable            :: aclo(:,:,:)
    complex,           allocatable            :: bclo(:,:,:)
    complex,           allocatable            :: cclo(:,:,:,:)
    complex,           allocatable            :: uunmt(:,:,:,:,:)
    complex,           allocatable            :: udnmt(:,:,:,:,:)
    complex,           allocatable            :: dunmt(:,:,:,:,:)
    complex,           allocatable            :: ddnmt(:,:,:,:,:)
    complex,           allocatable            :: acnmt(:,:,:,:,:)
    complex,           allocatable            :: bcnmt(:,:,:,:,:)
    complex,           allocatable            :: ccnmt(:,:,:,:,:)
    complex,           allocatable            :: vExt2IR(:, :, :, :)
    complex,           allocatable            :: vExt2MT(:, :, :, :)
    complex,           allocatable            :: E2ndOrdII(:, :)
    real,              allocatable            :: eigenVals(:)
    complex,           allocatable            :: eigenVecs(:, :)
    complex,           allocatable            :: w_vExt1IR(:, :, :)
    complex,           allocatable            :: vExt1IR(:, :, :)
    complex,           allocatable            :: vExt1MT(:, :, :, :, :)
    complex,           allocatable            :: vExt1IRtemp(:, :)
    complex,           allocatable            :: vExt1MTtemp(:, :, :, :)
    complex,           allocatable            :: vxc1IRKern(:)
    complex,           allocatable            :: rho0IRpw(:)
    complex,           allocatable            :: ylm(:, :)
    real,              allocatable            :: gWghts(:) ! gaussian weights belonging to gausPts
    real,              allocatable            :: dKernMTGPts(:, :, :)
    complex,           allocatable            :: dynMat(:, :)
    real                                      :: kExt(3)
    real                                      :: Gext(3)
    real                                      :: qpoint(3)
    complex,           allocatable            :: w_vExt2IR(:, :, :, :)

    complex,           allocatable            :: grVext0IR(:, :)
    complex,           allocatable            :: grVext0MT(:, :, :, :)
    complex,           allocatable            :: rho0IRncT(:, :)
    complex,           allocatable            :: rho0IRG(:)
    complex,           allocatable            :: gradRhoIR(:, :)
    complex                                   :: surfInt(3, 3)
    complex,           allocatable            :: vEff0IRpw(:, :)
    complex,        allocatable              :: vext1AltMT(:, :, :, :, :)
    complex,        allocatable              :: vext2AltMT(:, :, :, :)
    real,          allocatable                 :: r2Rho0MT(:, :, :, :)
    complex,       allocatable                :: r2GrRho0MT(:, :, :, :)

    complex :: testSum
    integer :: idirR
    integer :: ilh
    integer :: lm_pre
    complex, allocatable :: z1nG(:, :, :, :)


    write(logUnit, '(a)') 'Test of Dynamical Matrix for q = 0 (Goldstein modes).'
    write(logUnit, '(a)') '---------------------------------------------------------'

    ! Allocation of all required quantities
    allocate( z1Special(dimens%nbasfcn, dimens%neigd, 3) )
    allocate( rho1IR( ngdp, 3, atoms%nat ) )
    allocate( rho1MT( atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3, atoms%nat ) )
    allocate( gWghts(dimens%nspd), ylm(dimens%nspd, ( atoms%lmaxd + 1 )**2), dKernMTGPts(dimens%nspd, atoms%jmtd, atoms%nat))
    allocate( vxc1IRKern(ngdp) )
    allocate( vExt1IR( ngdp, 3, atoms%nat ), vExt1MT( atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3, atoms%nat ) )
    allocate( w_vExt1IR( ngdp, 3, atoms%nat ) )
    allocate( uu( 0 : atoms%lmaxd, atoms%nat, 3 ), du( 0 : atoms%lmaxd, atoms%nat, 3 ), dd( 0 : atoms%lmaxd, atoms%nat, 3 ), &
      & ud( 0 : atoms%lmaxd, atoms%nat, 3) )
    allocate( uunmt( 0 : atoms%lmaxd, 0: atoms%lmaxd, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      & ddnmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      & udnmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      & dunmt( 0 : atoms%lmaxd, 0 : atoms%lmaxd, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3) )
    allocate( acnmt( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      &bcnmt( 0 : atoms%lmaxd, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3), &
      &ccnmt(atoms%nlod, atoms%nlod, atoms%lmaxd * (atoms%lmaxd + 2), atoms%nat, 3) )
    allocate( aclo(atoms%nlod, atoms%nat, 3), bclo(atoms%nlod, atoms%nat, 3), cclo(atoms%nlod, atoms%nlod, atoms%nat, 3) )
    allocate( dynMat( 3 * atoms%nat, 3 * atoms%nat ) )

    ! Initialization
    rho1IR (:, :, :) = cmplx(0., 0.)
    rho1MT(:, :, :, :, :) = cmplx(0., 0.)
    vxc1IRKern(:) = cmplx(0., 0.)
    ylm(:, :) = cmplx(0., 0.)
    dKernMTGPts(:, :, :) = 0.
    gWghts(:) = 0.
    vExt1IR(:, :, :) = cmplx(0., 0)
    vExt1MT(:, :, :, :, :) = cmplx(0., 0.)
    w_vExt1IR(:, :, :) = cmplx(0.0, 0.0)
    dynMat(:, :) = cmplx(0., 0.)

    ! Goldstein modes are at a collective displacement of the lattice
    iqpt = 1
    qpoint(:) = 0.


    ! The factor r^2 has beeen divided out so that rho0MT is only rho0MT and not r^2 rho0MT as it is done in Fleur. The factor
    ! sqrt(4pi) for the zeroth component was already dividied out when constructing rho0MT in cdnmt routine in Fleur. Here to
    ! improve stability of the gradient routine we derive r2Rho0MT and divide out the r^2 again later. Doing so avoids the
    ! subtraction of small numbers close to the core.
    allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    r2Rho0MT(:, :, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)

    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          r2Rho0MT(imesh, ilh, itype, 1) = rho0MT(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
        end do
      end do
    end do

    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT(:, :, :, 1), r2GrRho0MT )

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                grRho0MT(imesh, lm, iatom, idir) = r2GrRho0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype) &
                                                                                                        & / atoms%rmsh(imesh, itype)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir


    allocate( z1nG(dimens%nbasfcn, 3, atoms%nat, maxval(nobd(:, :))) )
    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1
        invsAtNr = 1. / real(atoms%nat)


        ! Create external potential.
        ! NOTE: Only working for external potential because xc kernel has not been generated
        harSw = .false.
        extSw = .true.
        xcSw = .false.
        ! There is nothing canceling away as in Sternheimer therefore we need the full contribution of the external potential
        vExt1FullSw = .true.
        ! todo we don't need rho0IR in the argument list
        ! todo attention we do not have a factor r^2 here
        call GenVeff1( stars, cell, lathar, atoms, dimens, sym, harSw, extSw, xcSw, vExt1FullSw, ngdp, rho0IR(:, 1),                  &
          & rho0MT(:,:, :, 1), qpoint, rho1IR(:, :, iDatom), rho1MT(:, :, :, :, iDatom), grRho0MT, mlh_atom, nmem_atom, clnu_atom, &
          & memd_atom, gdp, vExt1IRtemp, vExt1MTtemp, vxc1IRKern, ylm, dKernMTGPts, gWghts, iDatom, iDtype, iqpt, ngdp, gdp ) ! add spin

        ! We are operating with routines expecting a more complicated array dimension structure. Therefore, we have to recast.
        do idir = 1, 3
          do iG = 1, ngdp
            vExt1IR(iG, idir, iDatom) = vExt1IRtemp(iG, idir)
          end do ! iG
        end do ! idir

        write(*, *) 'solve the r*2 vext or vext decision'
        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do oqn_l = 0, atoms%lmax(itype)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itype)
                    vExt1MT(imesh, lm, iatom, idir, iDatom) = vExt1MTtemp(imesh, lm, iatom, idir) !* atoms%rmsh(imesh, itype)**2
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! ieqat
          end do ! itype
        end do ! idir

        deallocate( vExt1IRtemp, vExt1MTtemp )

        ! grRho0MT is created within GenGrVeff0 which will be changed later to be consistent
        ! todo at the moment we calculate r^2 Vext
        ! todo we have to refine the gradient of the density later
!        deallocate(grRho0MT)
        grRhoTermSw = .false.
        !TODO REVIEW VEXTFULL BEFORE COMMIT
        ! we use veff1 instead of genGrVeff0
        call GenGrVeff0( atoms, cell, lathar, stars, dimens, memd_atom, ngdp, .false., .true., .false., gdp, &
          & rho0IR( :, 1 ), rho0MT(:, :, :, 1), nmem_atom, mlh_atom, clnu_atom, grVext0IR, grVext0MT, grRho0MT, gWghts, ylm, ylm, &
          & dKernMTGPts, vxc1IRKern, testGoldstein, grRhoTermSw ) ! add spin

        ! Set this 0 for safety reasons
        rho1IR(:, :, iDatom) = cmplx(0., 0.)
        rho1MT(:, :, :, :, iDatom) = cmplx(0., 0.)

        do idir = 1, 3
          call warpIRPot(stars, ngdp, idir, gdp, vExt1IR(:, :, iDatom), w_vExt1IR(:, idir, iDatom))
        end do

        uu = 0.0
        du = 0.0
        ud = 0.0
        dd = 0.0
        aclo = 0.0
        bclo = 0.0
        cclo = 0.0
        uunmt = 0.0
        ddnmt = 0.0
        udnmt = 0.0
        dunmt = 0.0
        acnmt = 0.0
        bcnmt = 0.0
        ccnmt = 0.0
        ! For q = 0, we assume z1 to be z1 = -i (k + G) z0
        do ikpt = 1, kpts%nkpt
!
          z1nG = cmplx(0., 0.)
          call ReadInz1( atoms, ikpt, 1, ikpt, nobd, nv, z1nG)
          z1Special(:, :, :) = cmplx(0.0, 0.0)
          do iband = 1, nobd(ikpt, 1)
            do idir = 1, 3
              do iBas = 1, nv(1, ikpt) + atoms%nlotot
                z1Special(iBas, iband, idir) = z1nG(iBas, idir, iDatom, iband)
              end do ! iBas
            end do ! idir
          end do ! iband
!          kExt(1:3) = matmul( cell%bmat(1:3, 1:3), kpts%bk(1:3, ikpt) )
          do idir = 1, 3
!            do iBas = 1, nv(1, ikpt)
!              Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(iBas, ikpt, 1)) )
!              do iband = 1, nobd(ikpt, 1) !todo LOs mit berücksichtigen, auch mit Gs?, + or - q?
!                z1Special(iBas, iband, idir) = z1Special(iBas, iband, idir) - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) )&
!                  & * z0(iBas, iband, ikpt, 1)
!
!              end do ! iband
!            end do ! iBas
!            do iBas = nv(1, ikpt) + 1, nv(1, ikpt) + atoms%nlotot
!              Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(kveclo(iBas - nv(1, ikpt), ikpt), ikpt, 1)) )
!              do iband = 1, nobd(ikpt, 1) !todo LOs mit berücksichtigen, auch mit Gs?, + or - q?
!                z1Special(iBas, iband, idir) = z1Special(iBas, iband, idir) - iu * cmplx(invsAtNr, 0) * ( Gext(idir) + kExt(idir) )&
!                  & * z0(iBas, iband, ikpt, 1)
!              end do ! iband
!            end do ! iBas
            ! calculate interstitial density with the special z1
            call calcRho1IRValDS( cell, input, results, nobd(ikpt, 1), nv(1, :), ikpt, iqpt, mapKpq2K(ikpt,iqpt), idir, gbas(:, :),&
              & z0(:, :, ikpt, 1), z1Special, rho1IR(:, :, iDatom), gdp2Ind, mapGbas, gdp2iLim, kpq2kPrVec)
          end do ! idir

          ! calculate k-dependent part of muffin-tin density with the special z1
          call calcKdepValRho1MT( atoms, dimens, sym, cell, kpts, input, usdus, results, ikpt, mapKpq2K(ikpt, iqpt), iDatom, nv,   &
            & mapGbas, gbas, nobd(:, 1), uu, du, ud, dd, aclo, bclo, cclo, uunmt, udnmt, dunmt, ddnmt, acnmt, bcnmt, ccnmt,        &
            & z0(:, :, ikpt, 1), z1Special, kveclo, rbas1, rbas2, ilo2p )
!
!          if (.false.) then
!            call storeZ1nG( atoms, ikpt, iqpt, mapKpq2K, iDatom, nobd, nv, z1Special )
!          end if

        end do ! ikpt


        ! Calculate the k-dependent parts of the density. By looping over the k-points, we perform the sum over the k-points.
        call multRadSolVzcRho1MT( atoms, input, aclo, bclo, cclo, acnmt, bcnmt, ccnmt, rbas1, rbas2, uu, du, ud, dd, uunmt, udnmt, &
          & dunmt, ddnmt, ilo2p, rho1MT(:, :, :, :, iDatom) )

        ! todo We are a bit inconsistent because we do not vary thoroughly the core terms of the density. Sometimes they are are
        ! sometimes not
!        rho1MT = cmplx(0., 0.)
        do idir = 1, 3
          do oqn_l = 0, atoms%lmax(iDtype) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + mqn_m + 1! todo does everythink start at 1?
              do imesh = 1, atoms%jri(iDtype)
                if (abs(rho1MT(imesh, lm, iDatom, idir, iDatom)) > 9e-5) write(*, *) 'control'
                ! If we are at the displaced atom we have to add the gradient of the density
                rho1MT(imesh, lm, iDatom, idir, iDatom) =  rho1MT(imesh, lm, iDatom, idir, iDatom) &
                  &- grRho0MT(imesh, lm, iDatom, idir)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! idir
      end do ! iDeqat
    end do ! iDtype


    ! Import rho0IR before adding coretail corrections to it. This makes no difference for Neon so we use full rho0IR
    ! todo test also for other systems
    if ( .false. ) then
      allocate(rho0IRncT(stars%n3d, 1))
      rho0IRncT(:, :) = cmplx(0.,0.)
      rewind(7800)
      read(7800) rho0IRncT
      allocate(rho0IRG(ngdp))
      rho0IRG = cmplx(0., 0.)
      call convertStar2G(rho0IRncT(:, 1), rho0IRG, stars, ngdp, gdp)
    end if

    ! The trace has to be remained by setting .true. the last argument otherwise this test does not work because Vext2 is not the
    ! gradient
    call GenVext2( atoms, cell, dimens, ngdp, gdp, vExt2IR, vExt2MT, .true. )

    ! Warp interstitial second-order external potential
    allocate(w_vExt2IR(ngdp, 3, 3, atoms%nat))
    w_vExt2IR = cmplx(0.0, 0.0)
    do idirC = 1, 3
      do idir = 1, 3
        call warpIRPot(stars, ngdp, idir, gdp, vExt2IR(:, :, idirC, iDatom), w_vExt2IR(:, idir, idirC, iDatom))
      end do
    end do

    !todo testVexts has to be outsourced and to be run seperately as a test. There is no need to leave it here in the long term
    !call testVexts( atoms, stars, lathar, cell, sym, ngdp, gdp, grVext0IR, clnu_atom, nmem_atom, mlh_atom, grVext0MT, vExt2IR,    &
    !  & vExt2MT )

    !todo this has to be expanded
    ! Checks whether the warping works correctly by inserting the reciprocal representation of the step function
    !call testWarpingHFdynMat(atoms, cell, stars, ngdp, 1, gdp, rho1IR, vExt1IR, w_vExt1IR)

    ! Expand interstitial densiyt from stars to plane-waves
    allocate(rho0IRpw(ngdp))
    rho0IRpw = cmplx(0., 0.)
    call convertStar2G( rho0IR(:, 1), rho0IRpw, stars, ngdp, gdp )

    ! Calculate Eii2 vanishing for q = 0 anyway but we leave it here as a test
    call CalcIIEnerg2(atoms, cell, dimens, qpts, stars, input, 1, ngdp, gdp, E2ndOrdII)

    ! Analytical tests of gradient routine and test of rho with conventional and enhanced method. Also many methods to integrate
    ! the integral vext1 rho1 and working method for MT Goldstone test
    ! Maybe it makes sense to outsource the call later. It was just an auxilliary routine to play around.
!    call testMTIntegrandsHFDynMat( atoms, stars, input, vacuum,   lathar, sym, cell, ngdp, gdp, clnu_atom, nmem_atom,      &
!      & mlh_atom, w_vExt2IR, vExt2MT, rho1IR, w_vExt1IR, -grVext0MT, E2ndOrdII )

    ! This routine was used to derive the correct approach to the HF dynmat integrals
!    call testSimpleRhoVextMTInt( atoms, cell, lathar, ngdp, mlh_atom, nmem_atom, clnu_atom, grVext0MT, vExt1MT, vExt2MT, rho0MT )

    call SetupDynMatHF( atoms, cell, lathar, ngdp, ngdp, mlh_atom, nmem_atom, clnu_atom, rho0IRpw, w_vExt2IR, vExt2MT, rho1IR, rho1MT,   &
      & w_vExt1IR, vExt1MT, rho0MT, E2ndOrdII, dynMat )

!    call SetupDynMatHF( atoms, cell, lathar, ngdp, mlh_atom, nmem_atom, clnu_atom, rho0IRpw, w_vExt2IR, vExt2MT, rho1IR, rho1MT,   &
!      & w_vExt1IR, vExt1MT, r2Rho0MT, E2ndOrdII, dynMat )

    write(*, '(a)') 'Complete Hellmann-Feynman Dynamical Matrix for q = 0'
    write(*, '(3(2(es16.8,1x),3x))') dynMat(1, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMat(2, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMat(3, :)

    if (.false.) then
      ! Compare with surface integral of dynamical matrix
      surfInt(:, :) = cmplx(0.0, 0.0)
      call CalcSurfIntIRDynMat( atoms, cell, ngdp, gdp, rho0IRpw, grVext0IR, surfInt )
      write(*, '(a)') 'IR Surface integral'
      write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)

      !todo find a solution for the correct mt-surface integral
      ! todo Does only work if goldstein test set true in grVext0MT routine, solve that
      surfInt(:, :) = cmplx(0.0, 0.0)
      call CalcSurfIntMTDynMat(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT, grVext0MT, surfInt)
      write(*, '(a)') 'MT Surface integral'
      write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)

      ! Difference
      write(*, '(a)') 'Difference'
      write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :) - dynMat(1, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :) - dynMat(2, :)
      write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :) - dynMat(3, :)
    end if

  end subroutine testGoldsteinModes

  ! For q = 0 the sum of the IR Hellmann-Feynman integrals contributing is this surface integral
  subroutine CalcSurfIntIRDynMat( atoms, cell, ngdp, gdp, rho0IRpw, grVext0IR, surfInt )

    use m_types
    use m_jPConstants, only : iu, tpi, fpi, c_im
    use m_ylm_old
    use m_sphbes

    implicit none

    ! Type parameter
    type(t_atoms),        intent(in)  :: atoms
    type(t_cell),         intent(in)  :: cell

    ! Scalar parameter
    integer                           :: ngdp

    ! Array parameter
    integer,              intent(in)  :: gdp(:, :)
    complex,              intent(in)  :: rho0IRpw(:)
    complex,              intent(in)  :: grVext0IR(:, :)
    complex,              intent(out) :: surfInt(3, 3)

    ! Scalar variables
    integer                           :: idirC
    integer                           :: idirR
    integer                           :: iatom
    integer                           :: itype
    integer                           :: ieqat
    integer                           :: iG
    integer                           :: it
    integer                           :: iGp
    complex                           :: phaseFac
    complex                           :: tSummedCitY1t
    complex                           :: pref
    complex                           :: surfIntNonMat

    ! Array variables
    integer                           :: gSum(3)
    real                              :: gSumCart(3)
    complex                           :: ylm(4)
    real                              :: sbes(0:1)

    surfInt(:, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      pref = fpi * iu * atoms%rmt(itype) * atoms%rmt(itype)
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do iG = 1, ngdp
          do iGp = 1, ngdp

            gSum(1:3) = gdp(1:3, iG) + gdp(1:3, iGp)
            gSumCart(1:3) = matmul( cell%bmat(1:3, 1:3), gSum(1:3) )

            ylm(:) = cmplx(0., 0.)
            call ylm4( 1, gSumCart, ylm )

            sbes(:) = 0
            call sphbes(1, norm2(gSumCart) * atoms%rmt(itype), sbes)

            phaseFac = exp( iu * tpi * dot_product(gSum(1:3), atoms%taual(1:3, iatom)))

            surfIntNonMat = pref * phaseFac * sbes(1) * rho0IRpw(iG)

            do idirC = 1, 3
              do idirR = 1, 3
                ! Corresponds to the magnetic quantum number m for the orbital quantum number l = 1
                tSummedCitY1t = cmplx(0., 0.)
                do it = -1, 1
                  tSummedCitY1t = tSummedCitY1t + c_im(idirR, it + 2) * ylm(it + 3)
                end do ! it

                surfInt(idirR, idirC) = surfInt(idirR, idirC) - surfIntNonMat * tSummedCitY1t * grVext0IR(iGp, idirC)
              end do ! idirR
            end do ! idirC
          end do ! iGp
        end do ! iG
      end do ! ieqat
    end do ! itype

  end subroutine CalcSurfIntIRDynMat

  ! For q = 0 the sum of the MT Hellmann-Feynman integrals contributing is this surface integral
  subroutine CalcSurfIntMTDynMat(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT, grVext0MT, surfInt)

    use m_types, only : t_atoms, t_sphhar
    use m_gaunt, only : Gaunt1
    use m_jpConstants, only : c_im

    implicit none

    ! Type parameters
    type(t_atoms),                     intent(in)  :: atoms
    type(t_sphhar),                    intent(in)  :: lathar

    ! Array parameters
    complex,                           intent(in)  :: clnu_atom(:, 0:, :)
    integer,                           intent(in)  :: nmem_atom(0:, :)
    integer,                           intent(in)  :: mlh_atom(:, 0:, :)
    real,                              intent(in)  :: rho0MT(:, 0:, :, :)
    complex,                           intent(in)  :: grVext0MT(:, :, :, :)
    complex,                           intent(out) :: surfInt(3, 3)

    ! Scalar variables
    integer                                        :: iatom
    integer                                        :: itype
    integer                                        :: ieqat
    integer                                        :: ptsym
    integer                                        :: ilh
    integer                                        :: oqn_l
    integer                                        :: imem
    integer                                        :: mqn_m
    integer                                        :: mqn_mp
    integer                                        :: oqn_lpp
    integer                                        :: mqn_mpp
    integer                                        :: lmpp
    real                                           :: gauntFactor
    integer                                        :: idirC
    integer                                        :: idirR

    surfInt(:, :) = 0
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            do mqn_mp = -1, 1
              ! lmax + 1 is okay for grVext1
             ! do oqn_lpp = abs(oqn_l - 1), oqn_l + 1
              do oqn_lpp = 0, atoms%lmax(itype)
                do mqn_mpp = -oqn_lpp, oqn_lpp
             !   mqn_mpp = mqn_m + mqn_mp
                lmpp = oqn_lpp * (oqn_lpp + 1) + 1 + mqn_mpp
                gauntFactor = Gaunt1( oqn_lpp, oqn_l, 1, mqn_mpp, mqn_m, mqn_mp, atoms%lmax(itype) + 1)
                do idirC = 1, 3
                  do idirR = 1, 3
                    surfInt(idirR, idirC) = surfInt(idirR, idirC) + conjg(c_im(idirR, mqn_mp + 2))           &
                      ! NOTE TODO NOTE TODO : This conj is probably wrong the conjg should be around the grVext0MT this gives
                      ! the correct muffin-tin behavior, i.e. the then the IR and the MT surface integral are equal.
                      & * conjg(rho0MT(atoms%jri(itype), ilh, itype, 1) * clnu_atom(imem, ilh, iatom)) *                           &
                      & grVext0MT(atoms%jri(itype), lmpp, idirC, iatom) * gauntFactor
                  end do ! idirC
                end do ! idirR
                end do
              end do ! oqn_lpp
            end do ! mqn_mp
          end do ! imem
        end do ! ilh
      end do ! ieqat
    end do ! itype

  end subroutine CalcSurfIntMTDynMat

  ! test the cutoffs of the warping function by performing the warping via the reciprocal representation of the step function
  ! this is caled from testGoldsteinModes but deprecated in this place, it is used in the same way within the tests of the pulay contributions
  subroutine testWarpingHFdynMat( atoms, cell, stars, ngdp, iDatom, gdp, rho1IR, vExt1IR, w_vExt1IR )

    use m_types, only : t_atoms, t_cell, t_stars
    use m_jpConstants, only : iu, fpi, tpi
    use m_sphbes
    use mod_juPhonUtils, only : convertStar2G

    implicit none

    type(t_atoms), intent(in) :: atoms
    type(t_cell),  intent(in) :: cell
    type(t_stars), intent(in) :: stars

    integer,       intent(in) :: ngdp
    integer,       intent(in) :: iDatom

    integer,       intent(in) :: gdp(:, :)
    complex,       intent(in) :: rho1IR(:, :, :)
    complex,       intent(in) :: vExt1IR(:, :, :)
    complex,       intent(in) :: w_vExt1IR(:, :, :)

    integer                   :: idirR
    integer                   :: idirC
    real                      :: pref
    integer                   :: iG
    integer                   :: iGp
    real                      :: GdiffCartAbsV
    complex                   :: stepFunc
    integer                   :: itype
    integer                   :: iatom
    integer                   :: ieqat
    complex                   :: dynMatHfFFT(3, 3)
    complex                   :: dynMatHfSUM(3, 3)
    real                      :: Gext(3)
    real                      :: Gdiff(3)
    real                      :: Gpext(3)
    real                      :: GdiffCart(3)
    real, allocatable         :: sbes(:)

    ! If we want to debug stars%ustep in stepf.F, we have to comment out the read in of wkf2 containing the stepfunction


    ! Integrate two vector quantities in the interstitial so that the Bloch character vanishes by warping one of them and
    ! calculating the dot_product
    dynMatHfFFT(:, :) = cmplx(0., 0.)
    do idirC =  1, 3
      do idirR =  1, 3
        dynMatHfFFT(idirR, idirC) = cell%omtil * dot_product( rho1IR(:, idirR, iDatom), w_vExt1IR(:, idirC, iDatom) )
      end do ! idirR
    end do ! idirC

    ! Integrate the same quantities by double-sum a product of these quantities and the step function in reciprocal space
    allocate(sbes(0:atoms%lmax(atoms%ntype)))

    sbes = cmplx(0., 0.)
    dynMatHfSUM(:, :) = cmplx(0., 0.)
    pref = fpi * atoms%rmt(1)**2
    do iG = 1, ngdp
      Gext(:) = matmul(cell%bmat, gdp(:, iG))
      do iGp = 1, ngdp
        Gdiff(1:3) = gdp(1:3, iG) - gdp(1:3, iGp)
        Gpext(1:3) = matmul(cell%bmat, gdp(1:3, iGp))
        GdiffCart(1:3) = Gext(1:3) - Gpext(1:3)
        GdiffCartAbsV = norm2(GdiffCart)
        stepFunc = cmplx(0., 0.)
        if (iG == iGp) then
          ! It's just the interstitial volume
          stepFunc = cell%omtil
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              stepFunc = stepFunc - atoms%volmts(itype)
            end do ! ieqat
          end do ! itype
        else
          iatom = 0
          do itype = 1, atoms%ntype
            call sphbes(atoms%lmax(itype), GdiffCartAbsV * atoms%rmt(itype), sbes)
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              stepFunc =  -pref * exp(-tpi * iu * dot_product(gDiff(1:3), atoms%taual(1:3, itype))) * sbes(1) / GdiffCartAbsV
            end do ! ieqat
          end do ! itype
        end if
        do idirC = 1, 3
          do idirR = 1, 3
            dynMatHfSUM(idirR, idirC) = dynMatHfSUM(idirR, idirC) + conjg(rho1IR(iG, idirR, iDatom)) * vExt1IR(iGp, idirC, iDatom) &
                                                                                                                        & * stepFunc
          end do ! idirR
        end do ! idirC
      end do ! iGp
    end do ! iG

    write(*, *) 'Test of warping'
    write(*, *) 'Max. deviation:', maxval(abs(dynMatHfSUM - dynMatHfFFT))

  end subroutine testWarpingHFdynMat

  ! this is for testGoldstein Modes and irrelevant now
  subroutine testMTIntegrandsHFDynMat( atoms, stars, input, vacuum,   lathar, sym, cell, ngdp, gdp, clnu_atom, nmem_atom,      &
        & mlh_atom, w_vExt2IR, vExt2MT, rho1IR, w_vExt1IR, vExt1MT, E2ndOrdII  )

    use m_types
    use mod_juPhonUtils, only : Fopen, Fclose
    use m_loddop
    use m_jpPlotObservables, only : plot1stDensVar
    use mod_juPhonUtils, only : convertStar2G
    use m_jpPotDensHelper, only : calcGrFinLH, calcGrR2FinLH, calcGrr2FinSH, calcGrFinSH
    use m_intgr, only : intgr3
    use m_jpConstants, only : fpi, tpi, pi
    use m_jpSetupDynMat, only : SetupDynMatHF

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_stars),                  intent(in)  :: stars
    type(t_input),                  intent(in)  :: input
    type(t_vacuum),                 intent(in)  :: vacuum
     
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_sym),                    intent(in)  :: sym
    type(t_cell),                   intent(in)  :: cell

    ! Scalar parameters
    integer,                        intent(in)  :: ngdp

    ! Array parameters
    integer,                        intent(in)  :: gdp(:, :)
    complex,                        intent(in)  :: clnu_atom(:, 0:, :)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    integer,                        intent(in)  :: mlh_atom(:, 0:, :)
    complex,                        intent(in)  :: w_vExt2IR(:, :, :, :) ! keep in memory
    complex,                        intent(in)  :: vExt2MT(:, :, :, :) ! keep in memory

    complex,                        intent(in)  :: vExt1MT(:, :, :, :)
    complex,                        intent(in)  :: rho1IR(:, :, :) ! read from hard disk
    complex,                        intent(in)  :: w_vExt1IR(:, :, :)     ! read from hard disk
    complex,                        intent(in)  :: E2ndOrdII(:, :) ! keep in memory

    ! Scalar variables
    integer(4)                                  :: iunit = 205        ! unit number for density
    integer                                     :: iter               ! irrelevant for juPhon
    integer                                     :: iatom
    integer                                     :: itype
    integer                                     :: ieqat
    integer                                     :: ptsym
    integer                                     :: idir
    integer                                     :: idirR
    integer                                     :: idirC
    integer                                     :: ilh
    integer                                     :: oqn_l
    integer                                     :: imem
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: imesh
    integer                                     :: counter
    real                                        :: intgrReal
    real                                        :: intgrImag

    ! Array variables
    complex                                     :: sumVec(3)
    complex,           allocatable              :: rho0IRAlt(:,:)
    complex,           allocatable              :: rho0IRPW(:, : )
    real,              allocatable              :: rho0MTAlt(:,:,:,:)
    real,              allocatable              :: testFunction(:,:,:,:)
    complex,           allocatable              :: testFunctionSH(:,:,:,:)
    real,              allocatable              :: rho0MTAltPure(:,:,:,:)
    complex,           allocatable              :: rho0MTSH(:,:,:,:)
    complex,           allocatable              :: grRho0MT(:, :, :, :)
    complex,           allocatable              :: grTestFunction(:, :, :, :)
    complex,           allocatable              :: grTestFunctionSH(:, :, :, :)
    complex,           allocatable              :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable              :: benchMarkTestFunction(:, :, :, :)
    complex,           allocatable              :: r2GrTestFunction(:, :, :, :)
    complex,           allocatable              :: r2GrTestFunctionSH(:, :, :, :)
    complex,           allocatable              :: corrGrR2fMT(:, :, :, :)
    character(8)                                :: dop, iop, name(10) ! irrlevant for juPhon
    complex,           allocatable              :: nzxy(:,:,:,:)      ! irrelevant for juPhon calcualtion
    real,              allocatable              :: nz(:,:,:)          ! irrelevant for juPhon calculation
    complex,           allocatable              :: rho1Dummy(:, :, :, :, :)
    complex,           allocatable              :: vext1MTDummy(:, :, :, :, :)
    complex,           allocatable              :: r2vext1(:, :, :)
    complex,           allocatable              :: r2GrVext1(:, :, :, :, :)
    complex,           allocatable              :: vExt2MTDummy(:, :, :, :)
    complex,           allocatable              :: r2Rho0MTsh(:, :, :)
    complex,           allocatable              :: grRho0MTsh(:, :, :, :)
    complex,           allocatable              :: r2GrRho0MTsh(:, :, :, :)

    real :: quotient

    real,          allocatable                  :: fReal(:)
    real,          allocatable                  :: fImag(:)

    complex                                     :: mtIntgrEnhd(3, 3)
    complex                                     :: mtIntgrConv(3, 3)
    complex                                     :: dynMatHF(3, 3)
    complex                                     :: surfInt(3, 3)

    ! Allocate required quantities
    allocate( fReal(atoms%jmtd) )
    allocate( fImag(atoms%jmtd) )

    allocate( rho0IRAlt(stars%ng3,input%jspins), nzxy(vacuum%nmzxyd, stars%nq2-1,2, input%jspins),                              &
      & nz(vacuum%nmzd, 2, input%jspins), rho0MTAlt(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
    allocate( rho0MTAltPure(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins),                                                 &
      & testFunction(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
    allocate( testFunctionSH( atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%ntype, input%jspins))
    allocate( benchMarkTestFunction(atoms%jmtd, ( atoms%lmaxd + 2)**2, atoms%nat, 3) )
    allocate( rho1Dummy(atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3, atoms%nat) )
    allocate( vext1MTDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat))
    allocate(r2GrVext1(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3, 3))
    allocate(vExt2MTDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, 3))
    allocate(r2vext1(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%ntype))
    rho1Dummy = cmplx(0., 0.)
    vExt1MTDummy = cmplx(0., 0.)
    rho0MTAltPure = cmplx(0., 0.)
    testFunction(:, :, :, :) = cmplx(0., 0.)

    ! Test whether it makes a difference when rho is divided and multiplied by r^2, therefore we need the unchanged density
    call Fopen(iunit, name='cdn1', status='old', action='read', form='unformatted')

    rewind (iunit)

    call loddop( input%jspins, stars%ng3, stars%ng2, vacuum%nmzxyd, vacuum%nmzd, atoms%jmtd, lathar%nlhd, atoms%ntype,          &
      & input%jspins, stars%ng3, stars%ng2, vacuum%nvac, atoms%ntype, sym%invs, sym%invs2, input%film, lathar%nlh, atoms%jri,   &
      & lathar%ntypsd, atoms%ntypsy, iunit, atoms%natd, atoms%neq, iop, dop, iter, rho0MTAlt, rho0IRAlt, nz, nzxy, name )

    call Fclose(iunit)


    ! Divide out factor r^2 so that rho0MT is only rho0MT and not r^2 rho0MT as it is done in Fleur. The factor sqrt(4pi) for the
    ! zeroth component was already dividied out when constructing rho0MT in cdnmt routine in Fleur.
    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          rho0MTAltPure(imesh, ilh, itype, 1) = rho0MTAlt(imesh, ilh, itype, 1) / atoms%rmsh(imesh, itype)                          &
            & / atoms%rmsh(imesh, itype)
        end do
      end do
    end do

    ! Unfold rho0MT from spherical harmonic to lattice harmonic representation
    allocate( rho0MTSH(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3), r2Rho0MTsh(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat), &
      & grRho0MTsh(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3), r2grRho0MTsh(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3) )
    rho0MTSH(:, :, :, :) = cmplx(0., 0.)
    r2Rho0MTsh(:, :, :) = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              rho0MTsh(imesh, lm, iatom, 1) = rho0MTsh(imesh, lm, iatom, 1) + rho0MTAltPure(imesh, ilh, itype, 1)                &
                                                                                                   & * clnu_atom(imem, ilh, iatom)
              r2Rho0MTsh(imesh, lm, iatom) = r2Rho0MTsh(imesh, lm, iatom) + rho0MTAlt(imesh, ilh, itype, 1)                      &
                                                                                                   & * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh
      end do
    end do

    ! Different ways of deriving the density
    if (.false.) then
      allocate( grRho0MTsh(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3),                                                        &
                                                                    & r2grRho0MTsh(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3) )
      grRho0MTsh(:, :, :, :) = cmplx(0., 0.)
      r2grRho0MTsh(:, :, :, :) = cmplx(0., 0.)
      call calcGrFinSH(atoms, rho0MTsh(:, :, :, 1), grRho0MTsh)
      call calcGrR2FinSH(atoms, r2Rho0MTsh(:, :, :), r2GrRho0MTsh)
      call calcGrFinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MTAltPure(:, :, :, 1), grRho0MT)
    end if

    call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MTAlt(:, :, :, 1), r2GrRho0MT)
    do idir = 1, 3
      do lm = 1, (atoms%lmaxd + 1)**2
        do imesh = 1, atoms%jri(1)
          rho1Dummy(imesh, lm, 1, idir, 1) = -r2GrRho0MT(imesh, lm, 1, idir) / atoms%rmsh(imesh, 1)**2
          vext1MTDummy(imesh, lm, 1, idir, 1) = vext1MT(imesh, lm, idir, 1)
        end do ! imesh
      end do ! lm
    end do ! idir

    if (.false.) then
      ! Writeout plot data for different ways of deriving the density
      do idir = 1, 3
        do lm = 1, (atoms%lmaxd + 1)**2
          do imesh = 1, atoms%jri(1)
            write(1100, '(3(i8),1x,2(es15.8))') idir, lm, imesh, grRho0MTsh(imesh, lm, 1, idir)
            write(1101, '(3(i8),1x,2(es15.8))') idir, lm, imesh, grRho0MT(imesh, lm, 1, idir)
            write(1102, '(3(i8),1x,2(es15.8))') idir, lm, imesh, r2GrRho0MTsh(imesh, lm, 1, idir)
            write(1103, '(3(i8),1x,2(es15.8))') idir, lm, imesh, r2GrRho0MT(imesh, lm, 1, idir)
          end do ! imesh
        end do ! lm
      end do ! idir
    end if

  ! for l = 1 write out gradrho and gradvext1
  do idir = 1, 3
    do mqn_m = -1, 1
      do imesh = 1, atoms%jri(1)
        lm = 3 + mqn_m
        write(1130, '(i8, i8, i8, 1x, 4(es15.8,1x))') idir, mqn_m, imesh, r2GrRho0MT(imesh, lm, 1, idir), vext1MTDummy(imesh, lm, 1, idir, 1)
      end do ! imesh
    end do ! mqn_m
  end do ! idir

  ! Write out rho and r^2 rho as a function of the mesh
  do imesh = 1, atoms%jri(1)
    write(1010, '(3(es15.8))') atoms%rmsh(imesh, 1), real(rho0MTsh(imesh, 1, 1, 1)), real(r2Rho0MTsh(imesh, 1, 1))
  end do ! imesh

  ! Write out various combinations of integrations to find the optimal configuration
  idir = 1
  lm = 2
  do imesh = 1, atoms%jri(1)
    write(1131, '(8(es15.8,1x))') atoms%rmsh(imesh, 1), real(r2GrRho0MT(imesh, lm, 1, idir)),                                      &
      & real(r2GrRho0MT(imesh, lm, 1, idir) / atoms%rmsh(imesh, 1)**2),  real(vext1MTDummy(imesh, lm, 1, idir, 1)),                &
      & real(vext1MTDummy(imesh, lm, 1, idir, 1) * r2GrRho0MT(imesh, lm, 1, idir) / atoms%rmsh(imesh, 1)**2),                      &
      & real(vext1MTDummy(imesh, lm, 1, idir, 1) * r2GrRho0MT(imesh, lm, 1, idir)),                                                &
      & real(vext1MTDummy(imesh, lm, 1, idir, 1) * atoms%rmsh(imesh, 1)**2), atoms%rmsh(imesh, 1)**2
  end do ! imesh
  write(*, *) 'fit finished'
  idir = 1
  lm = 4
  do imesh = 1, atoms%jri(1)
    write(1132, '(8(es15.8,1x))') atoms%rmsh(imesh, 1), real(r2GrRho0MT(imesh, lm, 1, idir)), real(r2GrRho0MT(imesh, lm, 1, idir)  &
      & / atoms%rmsh(imesh, 1)**2),  real(vext1MTDummy(imesh, lm, 1, idir, 1)), real(vext1MTDummy(imesh, lm, 1, idir, 1) *         &
      & r2GrRho0MT(imesh, lm, 1, idir) / atoms%rmsh(imesh, 1)**2), real(vext1MTDummy(imesh, lm, 1, idir, 1) *                      &
      & r2GrRho0MT(imesh, lm, 1, idir)), real(vext1MTDummy(imesh, lm, 1, idir, 1) * atoms%rmsh(imesh, 1)**2),                      &
      & atoms%rmsh(imesh, 1)**2
  end do ! imesh
  idir = 2
  lm = 2
  do imesh = 1, atoms%jri(1)
    write(1133, '(8(es15.8,1x))') atoms%rmsh(imesh, 1), aimag(r2GrRho0MT(imesh, lm, 1, idir)), aimag(r2GrRho0MT(imesh, lm, 1, idir)&
      & / atoms%rmsh(imesh, 1)**2),  aimag(vext1MTDummy(imesh, lm, 1, idir, 1)), real(vext1MTDummy(imesh, lm, 1, idir, 1) *        &
      & r2GrRho0MT(imesh, lm, 1, idir) / atoms%rmsh(imesh, 1)**2), real(vext1MTDummy(imesh, lm, 1, idir, 1) *                      &
      & r2GrRho0MT(imesh, lm, 1, idir)), aimag(vext1MTDummy(imesh, lm, 1, idir, 1) * atoms%rmsh(imesh, 1)**2),                     &
      & atoms%rmsh(imesh, 1)**2
  end do ! imesh
  idir = 2
  lm = 4
  do imesh = 1, atoms%jri(1)
    write(1134, '(8(es15.8,1x))') atoms%rmsh(imesh, 1), aimag(r2GrRho0MT(imesh, lm, 1, idir)), aimag(r2GrRho0MT(imesh, lm, 1, idir)&
      & / atoms%rmsh(imesh, 1)**2),  aimag(vext1MTDummy(imesh, lm, 1, idir, 1)), real(vext1MTDummy(imesh, lm, 1, idir, 1) *        &
      & r2GrRho0MT(imesh, lm, 1, idir) / atoms%rmsh(imesh, 1)**2), real(vext1MTDummy(imesh, lm, 1, idir, 1) *                      &
      & r2GrRho0MT(imesh, lm, 1, idir)), aimag(vext1MTDummy(imesh, lm, 1, idir, 1) * atoms%rmsh(imesh, 1)**2),                     &
      & atoms%rmsh(imesh, 1)**2
  end do ! imesh
  idir = 3
  lm = 3
  do imesh = 1, atoms%jri(1)
    write(1135, '(8(es15.8,1x))') atoms%rmsh(imesh, 1), real(r2GrRho0MT(imesh, lm, 1, idir)), real(r2GrRho0MT(imesh, lm, 1, idir) /&
      & atoms%rmsh(imesh, 1)**2),  real(vext1MTDummy(imesh, lm, 1, idir, 1)), real(vext1MTDummy(imesh, lm, 1, idir, 1) *           &
      & r2GrRho0MT(imesh, lm, 1, idir) / atoms%rmsh(imesh, 1)**2), real(vext1MTDummy(imesh, lm, 1, idir, 1) *                      &
      & r2GrRho0MT(imesh, lm, 1, idir)), real(vext1MTDummy(imesh, lm, 1, idir, 1) * atoms%rmsh(imesh, 1)**2),                      &
      & atoms%rmsh(imesh, 1)**2
  end do ! imesh

  do imesh = 1, atoms%jri(1)
    write(1136, '(8(es15.8,1x))') atoms%rmsh(imesh, 1), real(rho0MTAlt(imesh, 0, 1, 1)), real(rho0MTAlt(imesh, 0, 1, 1) /          &
      & atoms%rmsh(imesh, 1)),  real(rho0MTAlt(imesh, 0, 1, 1) / atoms%rmsh(imesh, 1)**2), 1. / atoms%rmsh(imesh, 1), 1. /         &
      & atoms%rmsh(imesh, 1)**2, atoms%rmsh(imesh, 1), atoms%rmsh(imesh, 1)**2
  end do ! imesh

  lm = 2
  idir = 1
  do imesh = 1, atoms%jri(1)
    write(1137, '(8(es15.8,1x))') atoms%rmsh(imesh, 1), real(r2GrRho0MT(imesh, lm, 1, idir)), real(r2GrRho0MT(imesh, lm, 1, idir) /&
      & atoms%rmsh(imesh, 1)),  real(r2GrRho0MT(imesh, lm, 1, idir) / atoms%rmsh(imesh, 1)**2), 1. / atoms%rmsh(imesh, 1), 1. /    &
      & atoms%rmsh(imesh, 1)**2, atoms%rmsh(imesh, 1), atoms%rmsh(imesh, 1)**2
  end do ! imesh


  if (.false.) then
    do idirC = 1, 3
      do idirR = 1, 3
        do lm = 1, (atoms%lmaxd + 1)**2
          do imesh = 1, atoms%jri(1)
            vExt2MTDummy(imesh, lm, idirR, idirC) = r2GrVext1(imesh, lm, 1, idirR, idirC) / atoms%rmsh(imesh, 1)**2
          end do ! imesh
        end do ! lm
      end do ! idirR
    end do! idirC

    do idir = 1, 3
      do lm = 1, (atoms%lmaxd + 1)**2
        do imesh = 1, atoms%jri(1)
          write(5922, '(3(i8),1x,2(es15.8))') idir, lm, imesh, 2 * vext2MTDummy(imesh, lm, idir, idir)
        end do ! imesh
      end do ! lm
    end do ! idir
  end if


    ! The interstitial quantities are wrong now, because we are only interested in MT quantities
    call SetupDynMatHF( atoms, cell, lathar, ngdp, ngdp, mlh_atom, nmem_atom, clnu_atom, rho0IRAlt(:, 1), w_vExt2IR, vext2MT, rho1IR,    &
      & rho1Dummy, w_vExt1IR, vext1MTDummy, rho0MTAlt, E2ndOrdII, dynMatHF )

    write(*, '(a)') 'dynmat HF'
    write(*, '(3(2(es16.8,1x),3x))') dynMatHF(1, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMatHF(2, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMatHF(3, :)

    call CalcSurfIntMTDynMat(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MTAltPure, -vext1MT, surfInt)
    write(*, '(a)') 'dynmat surfmt'
    write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)


    ! Test stability of gradient routine when integrating rho1 v1 for q = 0
    mtIntgrEnhd(:, :) = cmplx(0., 0.)
    do idirC = 1, 3
      do idirR = 1, 3
        do oqn_l = 0, atoms%lmax(1)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            ! Muffin-tin part of 7.112b PhDAK
            fReal(:) = 0.0
            fImag(:) = 0.0
            do imesh = 1, atoms%jri(1)
              fReal(imesh) = real ( conjg(-r2GrRho0MT(imesh, lm, 1, idirR)) * vExt1MT(imesh, lm, idirC, 1))
              fImag(imesh) = aimag( conjg(-r2GrRho0MT(imesh, lm, 1, idirR)) * vExt1MT(imesh, lm, idirC, 1))
            end do
            call intgr3(fReal, atoms%rmsh(1,1), atoms%dx(1), atoms%jri(1), intgrReal)
            call intgr3(fImag, atoms%rmsh(1,1), atoms%dx(1), atoms%jri(1), intgrImag)
            mtIntgrEnhd(idirR, idirC) = mtIntgrEnhd(idirR, idirC) + cmplx(intgrReal, intgrImag)
          end do
        end do
      end do ! idirR
    end do !idirC

    mtIntgrConv(:, :) = cmplx(0., 0.)
    do idirC = 1, 3
      do idirR = 1, 3
        do oqn_l = 0, atoms%lmax(1)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            ! Muffin-tin part of 7.112b PhDAK
            fReal(:) = 0.0
            fImag(:) = 0.0
            do imesh = 1, atoms%jri(1)
              fReal(imesh) = real ( conjg(-grRho0MT(imesh, lm, 1, idirR)) * vExt1MT(imesh, lm, idirC, 1)) * atoms%rmsh(imesh, 1)**2
              fImag(imesh) = aimag( conjg(-grRho0MT(imesh, lm, 1, idirR)) * vExt1MT(imesh, lm, idirC, 1)) * atoms%rmsh(imesh, 1)**2
            end do
            call intgr3(fReal, atoms%rmsh(1,1), atoms%dx(1), atoms%jri(1), intgrReal)
            call intgr3(fImag, atoms%rmsh(1,1), atoms%dx(1), atoms%jri(1), intgrImag)
            mtIntgrConv(idirR, idirC) = mtIntgrConv(idirR, idirC) + cmplx(intgrReal, intgrImag)
          end do
        end do
      end do
    end do
    write(*, *) 'difference is'
    do idirC = 1, 3
      do idirR = 1, 3
        write(*, '(2(i8),2(es15.8))') idirR, idirC, mtIntgrConv(idirR, idirC)
        write(*, '(2(i8),2(es15.8))') idirR, idirC, mtIntgrEnhd(idirR, idirC)
        write(*, '(2(i8),2(es15.8))') idirR, idirC, mtIntgrConv(idirR, idirC) - mtIntgrEnhd(idirR, idirC)
      end do ! idirR
    end do ! idirC

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(5790, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, grRho0MT(imesh, lm, iatom, idir)
                write(5793, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, r2GrRho0MT(imesh, lm, iatom, idir) &
                  &  / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir



    ! numerical stable and unstable gradient is integrated to show influence on integral
    sumVec(:) = cmplx(0., 0.)
    do idir = 1, 3
      oqn_l = 1
       do mqn_m = -1, 1
         lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
         call intgr3(real(grRho0MT(:, lm, 1, idir)) * atoms%rmsh(:, 1)**2,atoms%rmsh(1, 1), atoms%dx(1), atoms%jri(1), intgrReal)
         call intgr3(aimag(grRho0MT(:, lm, 1, idir))* atoms%rmsh(:,1)**2, atoms%rmsh(1, 1),atoms%dx(1), atoms%jri(1), intgrImag)
         sumVec(idir) = sumVec(idir) + cmplx(intgrReal, intgrImag)
       end do
    end do
    write(*, *) 'unprecise'
    write(*, '(6(es15.8))') sumVec

    sumVec(:) = cmplx(0., 0.)
    do idir = 1, 3
      oqn_l = 1
       do mqn_m = -1, 1
         lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
   !      call intgr3(real(r2GrRho0MT(:, lm, 1, idir) + corrGrR2fMT(:, lm, 1, idir)), atoms%rmsh(1, 1), atoms%dx(1), atoms%jri(1), intgrReal)
   !      call intgr3(aimag(r2GrRho0MT(:, lm, 1, idir) + corrGrR2fMT(:, lm, 1, idir)),atoms%rmsh(1, 1),  atoms%dx(1), atoms%jri(1), intgrImag)
         call intgr3(real(r2GrRho0MT(:, lm, 1, idir)), atoms%rmsh(1, 1), atoms%dx(1), atoms%jri(1), intgrReal)
         call intgr3(aimag(r2GrRho0MT(:, lm, 1, idir)),atoms%rmsh(1, 1),  atoms%dx(1), atoms%jri(1), intgrImag)
         sumVec(idir) = sumVec(idir) + cmplx(intgrReal, intgrImag)
       end do
    end do
    write(*, *) 'precise'
    write(*, '(6(es15.8))') sumVec

    ! Analytical test of gradient routine where we see a significant failure of the unprecise gradient routine
    write(*, *) 'analytical test'

   ! Testfunction f(r) = r
   testFunction(:, :, :, :) = cmplx(0., 0.)
   testFunctionSH(:, :, :, :) = cmplx(0., 0.)
   do imesh = 1, atoms%jri(1)
     testFunction(imesh, 0, 1, 1) = atoms%rmsh(imesh, 1)**3 * sqrt(fpi)
     testFunctionSH(imesh, 1, 1, 1) = cmplx(atoms%rmsh(imesh, 1)**3 * sqrt(fpi), 0.)
   end do

   allocate( r2GrTestFunctionSH(atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
   r2GrTestFunctionSH = cmplx(0., 0.)
   call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1),  r2GrTestFunction)
   call calcGrR2FinSH(atoms, testFunctionSH(:, :, :, 1), r2GrTestFunctionSH)

   do imesh = 1, atoms%jri(1)
    testFunction(imesh, 0, 1, 1) = testFunction(imesh, 0, 1, 1) / atoms%rmsh(imesh, 1)**2
    testFunctionSH(imesh, 1, 1, 1) = testFunctionSH(imesh, 1, 1, 1) / atoms%rmsh(imesh, 1)**2
   end do

   call calcGrFinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1), grTestFunction)
   call calcGrFinSH(atoms, testFunctionSH(:, :, :, 1), grTestFunctionSH)


   benchMarkTestFunction(:, :, :, :) = cmplx(0., 0)
   benchMarkTestFunction(:, 2, 1, 1) = cmplx( sqrt(tpi / 3), 0.)
   benchMarkTestFunction(:, 4, 1, 1) = cmplx(-sqrt(tpi / 3), 0.)
   benchMarkTestFunction(:, 2, 1, 2) = cmplx(0., sqrt(tpi / 3))
   benchMarkTestFunction(:, 4, 1, 2) = cmplx(0., sqrt(tpi / 3))
   benchMarkTestFunction(:, 3, 1, 3) = cmplx(sqrt(fpi / 3), 0.)

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                r2GrTestFunction(imesh, lm, iatom, idir) = r2GrTestFunction(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
                r2GrTestFunctionSH(imesh, lm, iatom, idir) = r2GrTestFunctionSH(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

   if (.true.) then
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(5794, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, grTestFunction(imesh, lm, iatom, idir)
                write(1110, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, grTestFunctionSH(imesh, lm, iatom, idir)
                write(5795, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, r2GrTestFunction(imesh, lm, iatom, idir)
                write(1111, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, r2GrTestFunctionSH(imesh, lm, iatom, idir)
                write(5796, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, benchMarkTestFunction(imesh, lm, iatom, idir)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir
  end if

    if (any (abs(grTestFunction(:, :, :, :) - benchMarkTestFunction(:, :, :, :)) > 1e-10) &
      & .or. any (abs(r2GrTestFunction(:, :, :, :) - benchMarkTestFunction(:, :, :, :)) > 1e-10) ) then
      write (*, *) 'Test f(r) = r failed!'
    else
      write (*, *) 'Test f(r) = r passed!'
    end if

   ! Testfunction f(r) = 1 / r
   testFunction(:, :, :, :) = cmplx(0., 0.)
   do imesh = 1, atoms%jri(1)
    testFunction(imesh, 0, 1, 1) = atoms%rmsh(imesh, 1) * sqrt(fpi)
    testFunctionSH(imesh, 1, 1, 1) = atoms%rmsh(imesh, 1) * sqrt(fpi)
   end do

   call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1),  r2GrTestFunction)
   call calcGrR2FinSH(atoms, testFunctionSH(:, :, :, 1), r2GrTestFunctionSH)

   do imesh = 1, atoms%jri(1)
    testFunction(imesh, 0, 1, 1) = testFunction(imesh, 0, 1, 1) / atoms%rmsh(imesh, 1)**2
    testFunctionSH(imesh, 1, 1, 1) = testFunctionSH(imesh, 1, 1, 1) / atoms%rmsh(imesh, 1)**2
   end do

   call calcGrFinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1), grTestFunction)
   call calcGrFinSH(atoms, testFunctionSH(:, :, :, 1), grTestFunctionSH)

   benchMarkTestFunction(:, :, :, :) = cmplx(0., 0)
   do imesh = 1, atoms%jri(1)
     quotient = -atoms%rmsh(imesh, 1) * atoms%rmsh(imesh, 1)
     benchMarkTestFunction(imesh, 2, 1, 1) = cmplx( sqrt(tpi / 3), 0.) / quotient
     benchMarkTestFunction(imesh, 4, 1, 1) = cmplx(-sqrt(tpi / 3), 0.) / quotient
     benchMarkTestFunction(imesh, 2, 1, 2) = cmplx(0., sqrt(tpi / 3)) / quotient
     benchMarkTestFunction(imesh, 4, 1, 2) = cmplx(0., sqrt(tpi / 3)) / quotient
     benchMarkTestFunction(imesh, 3, 1, 3) = cmplx(sqrt(fpi / 3), 0.) / quotient
   end do ! imesh

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                r2GrTestFunction(imesh, lm, iatom, idir) = r2GrTestFunction(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
                r2GrTestFunctionSH(imesh, lm, iatom, idir) = r2GrTestFunctionSH(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

  if (.true.) then
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(5797, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, grTestFunction(imesh, lm, iatom, idir)
                write(1112, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, grTestFunctionSH(imesh, lm, iatom, idir)
                write(5798, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, r2GrTestFunction(imesh, lm, iatom, idir)
                write(1113, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, r2GrTestFunctionSH(imesh, lm, iatom, idir)
                write(5799, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, benchMarkTestFunction(imesh, lm, iatom, idir)
                write(5800, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, grTestFunction(imesh, lm, iatom, idir) - benchMarkTestFunction(imesh, lm, iatom, idir)
                write(5801, '(4(i8),2x,2(f24.9))') idir, oqn_l, mqn_m, imesh, r2GrTestFunction(imesh, lm, iatom, idir) - benchMarkTestFunction(imesh, lm, iatom, idir)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir
  end if


  call Fopen(1000, name='grad1orConvtl', status='replace', action='write', form='formatted')
  write(1000, '(a)') '# Difference conventional routine to analytical result'
  write(1000, '(a)') '# Columns (only l = 1, each real and imag part): dir=x m=-1, dir=x m=-1, dir=y m=-1, dir=y m = 1, dir=z m = 0'
  do imesh = 1, atoms%jri(1)
    write(1000, '(6(f18.9))') atoms%rmsh(imesh, 1),                                                   &
                             & abs(real(grTestFunction(imesh, 2, 1, 1) - benchMarkTestFunction(imesh, 2, 1, 1))), &
                             & abs(real(grTestFunction(imesh, 4, 1, 1) - benchMarkTestFunction(imesh, 4, 1, 1))), &
                             & abs(aimag(grTestFunction(imesh, 2, 1, 2) - benchMarkTestFunction(imesh, 2, 1, 2))), &
                             & abs(aimag(grTestFunction(imesh, 4, 1, 2) - benchMarkTestFunction(imesh, 4, 1, 2))), &
                             & abs(real(grTestFunction(imesh, 3, 1, 3) - benchMarkTestFunction(imesh, 3, 1, 3)))
  end do
  call Fclose(1000)

  call Fopen(1000, name='grad1orEnhcd', status='replace', action='write', form='formatted')
  write(1000, '(a)') '# Difference enhanced routine to analytical result'
  write(1000, '(a)') '# Columns (only l = 1, each real and imag part): dir=x m=-1, dir=x m=-1, dir=y m=-1, dir=y m = 1, dir=z m = 0'
  do imesh = 1, atoms%jri(1)
    write(1000, '(6(f18.9))') atoms%rmsh(imesh, 1),                                                     &
                             & abs(real(r2GrTestFunction(imesh, 2, 1, 1) - benchMarkTestFunction(imesh, 2, 1, 1))), &
                             & abs(real(r2GrTestFunction(imesh, 4, 1, 1) - benchMarkTestFunction(imesh, 4, 1, 1))), &
                             & abs(aimag(r2GrTestFunction(imesh, 2, 1, 2) - benchMarkTestFunction(imesh, 2, 1, 2))), &
                             & abs(aimag(r2GrTestFunction(imesh, 4, 1, 2) - benchMarkTestFunction(imesh, 4, 1, 2))), &
                             & abs(real(r2GrTestFunction(imesh, 3, 1, 3) - benchMarkTestFunction(imesh, 3, 1, 3)))
  end do
  call Fclose(1000)

   ! Testfunction f(r) = 1 / r^2
   testFunction(:, :, :, :) = cmplx(0., 0.)
   testFunctionSH(:, :, :, :) = cmplx(0., 0.)
   do imesh = 1, atoms%jri(1)
    testFunction(imesh, 0, 1, 1) = sqrt(fpi)
    testFunctionSH(imesh, 1, 1, 1) = sqrt(fpi)
   end do
   deallocate(r2GrTestFunction, grTestFunction)

   call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1),  r2GrTestFunction)
   call CalcGrR2FinSH(atoms, testFunctionSH(:, :, :, 1), r2GrTestFunctionSH)

   do imesh = 1, atoms%jri(1)
!    testFunction(imesh, 0, 1, 1) = testFunction(imesh, 0, 1, 1) / atoms%rmsh(imesh, 1)**2
!    testFunctionSH(imesh, 1, 1, 1) = testFunctionSH(imesh, 1, 1, 1) / atoms%rmsh(imesh, 1)**2
   end do

   call calcGrFinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, testFunction(:, :, :, 1), grTestFunction)
   write(1115, '(2(f15.8))') grTestFunction
   call calcGrFinSH(atoms, testFunctionSH(:, :, :, 1), grTestFunctionSH)

   benchMarkTestFunction(:, :, :, :) = cmplx(0., 0)
   do imesh = 1, atoms%jri(1)
     quotient = atoms%rmsh(imesh, 1)**3
     benchMarkTestFunction(imesh, 2, 1, 1) = cmplx(-sqrt(8 * pi / 3), 0.) / quotient
     benchMarkTestFunction(imesh, 4, 1, 1) = cmplx( sqrt(8 * pi / 3), 0.) / quotient
     benchMarkTestFunction(imesh, 2, 1, 2) = cmplx(0., -sqrt(8 * pi / 3)) / quotient
     benchMarkTestFunction(imesh, 4, 1, 2) = cmplx(0., -sqrt(8 * pi / 3)) / quotient
     benchMarkTestFunction(imesh, 3, 1, 3) = cmplx(-4 * sqrt(pi / 3), 0.) / quotient
   end do ! imesh

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                r2GrTestFunction(imesh, lm, iatom, idir) = r2GrTestFunction(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
                r2GrTestFunctionSH(imesh, lm, iatom, idir) = r2GrTestFunctionSH(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

   if (.true.) then
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                write(5802, '(4(i8),2x,2(f30.9))') idir, oqn_l, mqn_m, imesh, grTestFunction(imesh, lm, iatom, idir)
                write(1114, '(4(i8),2x,2(f30.9))') idir, oqn_l, mqn_m, imesh, grTestFunctionSH(imesh, lm, iatom, idir)
                write(5803, '(4(i8),2x,2(f30.9))') idir, oqn_l, mqn_m, imesh, r2GrTestFunction(imesh, lm, iatom, idir)
                write(1115, '(4(i8),2x,2(f30.9))') idir, oqn_l, mqn_m, imesh, r2GrTestFunctionSH(imesh, lm, iatom, idir)
                write(5804, '(4(i8),2x,2(f30.9))') idir, oqn_l, mqn_m, imesh, benchMarkTestFunction(imesh, lm, iatom, idir)
                write(5805, '(4(i8),2x,2(f30.9))') idir, oqn_l, mqn_m, imesh, grTestFunction(imesh, lm, iatom, idir) - benchMarkTestFunction(imesh, lm, iatom, idir)
                write(5806, '(4(i8),2x,2(f30.9))') idir, oqn_l, mqn_m, imesh, r2GrTestFunction(imesh, lm, iatom, idir) - benchMarkTestFunction(imesh, lm, iatom, idir)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir
  end if

  call Fopen(1000, name='grad1orSqConvtl', status='replace', action='write', form='formatted')
  write(1000, '(a)') '# Difference conventional routine to analytical result'
  write(1000, '(a)') '# Columns (only l = 1, each real and imag part): dir=x m=-1, dir=x m=-1, dir=y m=-1, dir=y m = 1, dir=z m = 0'
  do imesh = 1, atoms%jri(1)
    write(1000, '(6(f22.9))') atoms%rmsh(imesh, 1),                                                   &
                             & abs(real(grTestFunction(imesh, 2, 1, 1) - benchMarkTestFunction(imesh, 2, 1, 1))), &
                             & abs(real(grTestFunction(imesh, 4, 1, 1) - benchMarkTestFunction(imesh, 4, 1, 1))), &
                             & abs(aimag(grTestFunction(imesh, 2, 1, 2) - benchMarkTestFunction(imesh, 2, 1, 2))), &
                             & abs(aimag(grTestFunction(imesh, 4, 1, 2) - benchMarkTestFunction(imesh, 4, 1, 2))), &
                             & abs(real(grTestFunction(imesh, 3, 1, 3) - benchMarkTestFunction(imesh, 3, 1, 3)))
  end do
  call Fclose(1000)

  call Fopen(1000, name='grad1orSqEnhcd', status='replace', action='write', form='formatted')
  write(1000, '(a)') '# Difference enhanced routine to analytical result'
  write(1000, '(a)') '# Columns (only l = 1, each real and imag part): dir=x m=-1, dir=x m=-1, dir=y m=-1, dir=y m = 1, dir=z m = 0'
  do imesh = 1, atoms%jri(1)
    write(1000, '(6(f22.9))') atoms%rmsh(imesh, 1),                                                     &
                             & abs(real(r2GrTestFunction(imesh, 2, 1, 1) - benchMarkTestFunction(imesh, 2, 1, 1))), &
                             & abs(real(r2GrTestFunction(imesh, 4, 1, 1) - benchMarkTestFunction(imesh, 4, 1, 1))), &
                             & abs(aimag(r2GrTestFunction(imesh, 2, 1, 2) - benchMarkTestFunction(imesh, 2, 1, 2))), &
                             & abs(aimag(r2GrTestFunction(imesh, 4, 1, 2) - benchMarkTestFunction(imesh, 4, 1, 2))), &
                             & abs(real(r2GrTestFunction(imesh, 3, 1, 3) - benchMarkTestFunction(imesh, 3, 1, 3)))
  end do
  call Fclose(1000)

  ! Plot output of density and its gradient
    allocate(rho0IRPW(ngdp, 3))
    rho0IRPW = cmplx(0.0, 0.0)
    rho0MTSH = cmplx(0.0, 0.0)
    call convertStar2G(rho0IRAlt(:, 1), rho0IRPW(:, 1), stars, ngdp, gdp)
      rho0MTsh = 0
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ptsym = atoms%ntypsy(iatom)
          do ilh = 0, lathar%nlh(ptsym)
            oqn_l = lathar%llh(ilh, ptsym)
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                rho0MTsh(imesh, lm, iatom, 1) = rho0MTsh(imesh, lm, iatom, 1) + rho0MTAltPure(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
              end do ! imesh
            end do ! imem
          end do ! ilh
        end do
      end do
    ! Plot rho0 with counter 0
    counter = 0
    call plot1stDensVar(sym, cell, input, lathar, atoms, rho0IRPW(:, :), ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MTAltPure,&
      & [0., 0., 0.], rho0MTsh(:, :, :, :), 1., 1., 1., 1, counter, 1)



    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                grRho0MT(imesh, lm, iatom, idir) = grRho0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir
    counter = 1
    rho0IRPW = cmplx(0., 0.)
    call plot1stDensVar(sym, cell, input, lathar, atoms, rho0IRPW(:, :), ngdp, gdp, mlh_atom, nmem_atom, clnu_atom, rho0MTAltPure,&
      & [0., 0., 0.], grRho0MT(:, :, :, :), 1., 1., 1., 1, counter, 1)

  end subroutine testMTIntegrandsHFDynMat

  ! Just tests different algorithms of generating the external potentials i.e. by deriving it numerically/analytically or by using the Weinert method
  ! is part of the testGoldsteinModes and deprecated
  subroutine testVexts( atoms, stars, lathar, cell, sym, ngdp, gdp, grVext0IR, clnu_atom, nmem_atom, mlh_atom, grVext0MT, vExt2IR, &
      & vExt2MT )

    use m_types, only : t_atoms, t_sphhar, t_stars, t_cell, t_sym
    use mod_juPhonUtils, only : Fopen, Fclose
    use mod_juPhonUtils, only : convertStar2G
    use m_jpConstants, only : fpi, iu
    use m_jpPotDensHelper, only : calcGrR2FinLH, CalcGrR2FinSH, plotCompVnVbench, calcGrFinSH
    use m_jpTestPotential, only : checkjuPhPots

    implicit none

    ! Type parameters
    type(t_atoms),               intent(in)  :: atoms
    type(t_stars),               intent(in)  :: stars
    type(t_sym),                 intent(in)  :: sym
    type(t_sphhar),              intent(in)  :: lathar
    type(t_cell),                intent(in)  :: cell

    ! Scalar parameters
    integer,                     intent(in)  :: ngdp

    ! Array parameters
    integer,                     intent(in)  :: gdp(:, :)
    complex,                     intent(in)  :: grVext0IR(:, :)
    complex,                     intent(in)  :: clnu_atom(:, 0:, :)
    integer,                     intent(in)  :: nmem_atom(0:, :)
    integer,                     intent(in)  :: mlh_atom(:, 0:, :)
    complex,                     intent(in)  :: grVext0MT(:, :, :, :)
    complex,                     intent(in)  :: vExt2IR(:, :, :, :)
    complex,                     intent(in)  :: vExt2MT(:, :, :, :)

    ! Scalar variables
    integer                                  :: idir
    integer                                  :: idirR
    integer                                  :: idirC
    integer                                  :: iG
    integer                                  :: iatom
    integer                                  :: itype
    integer                                  :: ieqat
    integer                                  :: ptsym
    integer                                  :: ilh
    integer                                  :: oqn_l
    integer                                  :: lm_pre
    integer                                  :: imem
    integer                                  :: mqn_m
    integer                                  :: lm
    integer                                  :: imesh
    real                                     :: paPoX
    real                                     :: paPoY
    real                                     :: paPoZ
    real                                     :: resolution

    ! Array variables
    complex,        allocatable              :: vpwStar(:, :)
    complex,        allocatable              :: npsExtG(:)
    complex,        allocatable              :: vpwStarDiff(:, :)
    real,           allocatable              :: vrFleur(:, :, :, :)
    real,           allocatable              :: vrFleurDiff(:, :, :, :)
    real,           allocatable              :: r2vExt0MT(:, :, :)
    complex,        allocatable              :: r2GrvExt0MT(:, :, :, :)
    complex,        allocatable              :: r2GrGrTVext0MT(:, :, :, :, :)
    complex,        allocatable              :: vExt0IRFleurDummy(:, :)
    complex,        allocatable              :: vExt0MTFleurDummy(:, :, :, :)
    complex,        allocatable              :: grvExt0IRFleurDummy(:, :)
    complex,        allocatable              :: grvExt0MTFleurDummy(:, :, :, :)
    complex,        allocatable              :: grGrTvExt0IRFleurDummy(:, :, :)
    complex,        allocatable              :: grGrTvExt0MTTest(:, :, :, :)
    complex,        allocatable              :: vext1AltMT(:, :, :, :, :)
    complex,        allocatable              :: vext2AltMT(:, :, :, :)

    real                                     :: Gext(3)

    ! Allocate required arrays
    allocate( vpwStar(stars%n3d, 1), vpwStarDiff(stars%n3d, 1) )
    allocate( vrFleur(atoms%jmtd, 0:lathar%nlhd, atoms%ntypd, 1), &
            & vrFleurDiff(atoms%jmtd, 0:lathar%nlhd, atoms%ntypd, 1) )
    allocate( r2vExt0MT(atoms%jmtd, 0:lathar%nlhd, atoms%ntypd), r2GrvExt0MT(atoms%jmtd, ( atoms%lmaxd + 2)**2, atoms%nat, 3 ) )
    allocate( npsExtG( ngdp) )
    allocate( vExt0IRFleurDummy( ngdp, 3), vExt0MTFleurDummy(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3) )
    allocate( grVExt0IRFleurDummy( ngdp, 3), grVExt0MTFleurDummy(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3) )

    npsExtG = cmplx(0., 0.)
    vpwStar(:, :) = cmplx(0., 0.)
    vpwStarDiff(:, :) = cmplx(0., 0.)
    vrFleur(:, :, :, :) = cmplx(0., 0.)
    vrFleurDiff(:, :, :, :) = cmplx(0., 0.)
    r2vExt0MT(:, :, :) = cmplx(0., 0.)
    r2GrvExt0MT(:, :, :, :) = cmplx(0., 0.)
    vExt0IRFleurDummy = cmplx(0., 0.)
    vExt0MTFleurDummy = cmplx(0., 0.)
    grVExt0IRFleurDummy = cmplx(0., 0.)
    grVExt0MTFleurDummy = cmplx(0., 0.)

    ! Determine pseudodensity for Vext0IR in plane-wave represenation and Vext0MT in lattice harmonic representation
    ! from a difference of Fleur Coulomb and Hartree potential
    call fopen( 1000, name='psRhoIRSt_Coul', status='old', action='read', form='unformatted' )
      read( 1000 ) vpwStar
    call fclose( 1000 )

    call fopen( 1000, name='psRhoIRSt_Hart', status='old', action='read', form='unformatted' )
      read( 1000 ) vpwStarDiff
    call fclose( 1000 )

    vpwStar = vpwStar - vpwStarDiff
    call convertStar2G( vpwStar(:, 1), npsExtG, stars, ngdp, gdp )


    call fopen(1000, name='v0MTFLEUR_coulTest', status='old', action='read', form='unformatted')
      read(1000) vrFLEUR
    call fclose(1000)

    call fopen(1000, name='v0MTFLEUR_hart', status='old', action='read', form='unformatted')
      read(1000) vrFLEURDiff
    call fclose(1000)

    vrFLEUR = vrFLEUR - vrFLEURDiff


    ! Compare pseudodensities of grVext0 produced with Weinert in Fleur and read in from Fleur and store grVext0IR originating from
    ! vExt0IR pseudodensity read in from FLEUR
    do idir = 1, 3
      do iG = 1, ngdp
        Gext(:) = matmul(cell%bmat, gdp(:, iG))
        write(4187, '(2(i8),1x,2(f15.8))') idir, iG, grVext0IR(iG, idir) * norm2(Gext)**2 / fpi
        write(4188, '(2(i8),1x,2(f15.8))') idir, iG, iu * Gext(idir) * npsExtG(iG)
        if( norm2(Gext) > 1e-10) then
          grVExt0IRFleurDummy(iG, idir) = iu * Gext(idir) * fpi / norm2(Gext)**2 * npsExtG(iG)
        end if
      end do ! iG
    end do ! idirR

    ! Generate interstitial potential of FLEUR
    do iG = 1, ngdp
      Gext(:) = matmul(cell%bmat, gdp(:, iG))
      if ( norm2(Gext) > 1e-10) then
        vExt0IRFleurDummy( iG, 1 ) = fpi * npsExtG(iG) / norm2(Gext)**2
      end if
    end do ! iG

    ! Unfold lattice-harmonic representation of vExt0MT from FLEUR to spherical-harmonic representation and calculate a version of
    ! it multiplied with r^2 for numerical stability
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1) + 1
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            do imesh = 1, atoms%jri(itype)
!              write(4186, '(3(i8),2x,2(f20.8))') oqn_l, mqn_m, imesh, vrFLEUR(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
              vExt0MTFleurDummy(imesh, lm, iatom, 1) = vrFLEUR(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
          do imesh = 1, atoms%jri(itype)
            r2vExt0MT(imesh, ilh, itype) = vrFLEUR(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype)**2
          end do ! imesh
        end do ! ilh
      end do ! ieqat
    end do ! itype

    ! Numerical stable muffin-tin gradient routine
    call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2vExt0MT, r2GrvExt0MT)

    ! The result of the muffin-tin gradient routine has to be divided by r^2
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                grVExt0MTFleurDummy(imesh, lm, iatom, idir) = r2GrvExt0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    ! Plot Vext0 from Fleur
    paPoX = 1.
    paPoY = 1.
    paPoZ = 1.
    resolution = 600
    call plotCompVnVbench( atoms, cell, paPoX, paPoY, paPoZ, 1, ngdp, resolution, gdp, vExt0IRFleurDummy,             &
       & vExt0IRFleurDummy, vExt0MTFleurDummy, vExt0MTFleurDummy(:, :, 1, :) )

    ! Plot grVext0 generated by deriving Vext0 of Fleur
    paPoX = 1.
    paPoY = 1.
    paPoZ = 1.
    resolution = 600
    call plotCompVnVbench( atoms, cell, paPoX, paPoY, paPoZ, 2, ngdp, resolution, gdp, grVExt0IRFleurDummy, grVext0IR,             &
       & grVExt0MTFleurDummy, grVext0MT(:, :, :, 1) )

    ! Section the logfile by a heading
    write(100, '(a)') 'Continuity check of grVext0 made from Fleur Vext0'
    call checkjuPhPots(30, atoms, ngdp, gdp, cell, lathar, sym, grVext0IR, grVext0MT, 100)

    ! Compare quantitatively the muffin-tin potentials of grVext0, i.e. the derived Vext0 from Fleur with the Weinert version of
    ! grVext0 produced in juphon also by calculating the difference
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                ! without dividing by r^2 l = 1 is a constant
                write(4189, '(4(i8),2(f25.8))') idir, oqn_l, mqn_m, imesh, r2GrvExt0MT(imesh, lm, iatom, idir)                     &
                                                                                                     & / atoms%rmsh(imesh, itype)**2
                write(4190, '(4(i8),2(f25.8))') idir, oqn_l, mqn_m, imesh, grvExt0MT(imesh, lm, idir, iatom)
                write(4191, '(4(i8),(es15.8))') idir, oqn_l, mqn_m, imesh, abs(grvExt0MT(imesh, lm, idir, iatom)                   &
                                                             & - r2GrvExt0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2 )
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    ! Plot whether error of gradient routine is proportional to increment of mesh. We compare the derived Vext0 from Fleur with the
    ! Weinert version of grVext0 produced in juphon for all lm
    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              write(4192, '(4(i8),2(f20.8))')  idir, oqn_l, mqn_m, 1, atoms%rmsh(1, itype), abs(grvExt0MT(1, lm, idir, iatom)      &
                & - r2GrvExt0MT(1, lm, iatom, idir) / atoms%rmsh(1, itype)**2) / atoms%rmsh(1, itype)
              do imesh = 2, atoms%jri(1)
                write(4192, '(4(i8),2(f20.8))') idir, oqn_l, mqn_m, imesh, atoms%rmsh(imesh, itype),                               &
                  & abs(grvExt0MT(imesh, lm, idir, iatom) - r2GrvExt0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, 1)**2)         &
                  & / (atoms%rmsh(imesh, 1) - atoms%rmsh(imesh - 1, 1))
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    ! Plot whether error of gradient routine is proportional to increment of mesh. We compare the derived Vext0 from Fleur with the
    ! Weinert version of grVext0 produced in juphon for l = 1 and all relevant ms and directions giving a contribution
    ! x-direction l = 1, m = -1
    write(4193, '(2(f20.8))')  atoms%rmsh(1, 1), abs(grvExt0MT(1, 2, 1, 1) - r2GrvExt0MT(1, 2, 1, 1) / atoms%rmsh(1, 1)**2)        &
                                                                                                                & / atoms%rmsh(1, 1)
    do imesh = 2, atoms%jri(1)
      write(4193, '(2(f20.8))') atoms%rmsh(imesh, 1), abs(grvExt0MT(imesh, 2, 1, 1) - r2GrvExt0MT(imesh, 2, 1, 1)                  &
                                                     / atoms%rmsh(imesh, 1)**2) / (atoms%rmsh(imesh, 1) - atoms%rmsh(imesh - 1, 1) )
    end do ! imesh

    ! x-direction l = 1, m = 1
    write(4194, '(2(f20.8))')  atoms%rmsh(1, 1), abs(grvExt0MT(1, 4, 1, 1) - r2GrvExt0MT(1, 4, 1, 1) / atoms%rmsh(1, 1)**2) /      &
                                                                                                                  & atoms%rmsh(1, 1)
    do imesh = 2, atoms%jri(1)
      write(4194, '(2(f20.8))') atoms%rmsh(imesh, 1), abs(grvExt0MT(imesh, 4, 1, 1) - r2GrvExt0MT(imesh, 4, 1, 1) /                &
                                                     & atoms%rmsh(imesh, 1)**2) / (atoms%rmsh(imesh, 1) - atoms%rmsh(imesh - 1, 1) )
    end do ! imesh

    ! y-direction l = 1, m = 0
    write(4195, '(2(f20.8))')  atoms%rmsh(1, 1), abs(grvExt0MT(1, 2, 2, 1) - r2GrvExt0MT(1, 2, 1, 2) / atoms%rmsh(1, 1)**2)        &
                                                                                                                & / atoms%rmsh(1, 1)
    do imesh = 2, atoms%jri(1)
      write(4195, '(2(f20.8))') atoms%rmsh(imesh, 1), abs(grvExt0MT(imesh, 2, 2, 1) - r2GrvExt0MT(imesh, 2, 1, 2)                  &
                                                   & / atoms%rmsh(imesh, 1)**2) / (atoms%rmsh(imesh, 1) - atoms%rmsh(imesh - 1, 1) )
    end do ! imesh

    ! y-direction l = 1, m = -1
    write(4196, '(2(f20.8))')  atoms%rmsh(1, 1), abs(grvExt0MT(1, 4, 2, 1) - r2GrvExt0MT(1, 4, 1, 2) / atoms%rmsh(1, 1)**2)        &
                                                                                                                  / atoms%rmsh(1, 1)
    do imesh = 2, atoms%jri(1)
      write(4196, '(2(f20.8))') atoms%rmsh(imesh, 1), abs(grvExt0MT(imesh, 4, 2, 1) - r2GrvExt0MT(imesh, 4, 1, 2)                  &
                                                   & / atoms%rmsh(imesh, 1)**2) / (atoms%rmsh(imesh, 1) - atoms%rmsh(imesh - 1, 1) )
    end do ! imesh

    ! y-direction l = 1, m = 1
    write(4197, '(2(f20.8))')  atoms%rmsh(1, 1), abs(grvExt0MT(1, 3, 3, 1) - r2GrvExt0MT(1, 3, 1, 3) / atoms%rmsh(1, 1)**2)        &
                                                                                                                & / atoms%rmsh(1, 1)
    do imesh = 2, atoms%jri(1)
      write(4197, '(2(f20.8))') atoms%rmsh(imesh, 1), abs(grvExt0MT(imesh, 3, 3, 1) - r2GrvExt0MT(imesh, 3, 1, 3)                  &
                                                     / atoms%rmsh(imesh, 1)**2) / (atoms%rmsh(imesh, 1) - atoms%rmsh(imesh - 1, 1) )
    end do ! imesh


    ! Generate interstitial gradGradTvExt0 by deriving the interstitial external Fleur potential Vext0 the second time
    allocate( grGrTvExt0IRFleurDummy(ngdp, 3, 3) )
    grGrTvExt0IRFleurDummy = cmplx(0., 0.)
    do idirC = 1, 3
      do idirR = 1, 3
        do iG = 1, ngdp
          Gext(:) = matmul(cell%bmat, gdp(:, iG))
          grGrTvExt0IRFleurDummy(iG, idirR, idirC) = iu * Gext(idirR) * grVExt0IRFleurDummy(iG, idirC)
        end do ! iG
      end do ! idirR
    end do ! idirC

    ! Divide the twice derived Vext0MT from FLEUR by r^2 to test the unstable numerical MT gradient
    allocate( vext1AltMT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat), vext2AltMT(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, 3) )
    vext1AltMT(:, :, :, :, :) = cmplx(0., 0.)
    do idirR = 1, 3
      do oqn_l = 0, atoms%lmax(1)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          do imesh = 1, atoms%jri(1)
            vext1AltMT(imesh, lm, 1, idirR, 1) = r2GrvExt0MT(imesh, lm, 1, idirR) / atoms%rmsh(imesh, 1)**2
          end do ! imesh
        end do ! mqn_m
      end do ! oqn_l
    end do ! idirR

    ! Compare the twice derived Vext0MT from Fleur with
    allocate( r2GrGrTVext0MT(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3, 3) )
    r2GrGrTVext0MT(:, :, :, :, :) = cmplx(0., 0.)
    vext2AltMT(:, :, :, :) = cmplx(0., 0.)
    do idirC = 1, 1!3
      !todo The following line is wrong for more than one atom because the routine only runs over different types
      call CalcGrR2FinSH(atoms, r2GrVext0MT(:, :, :, idirC), r2GrGrTVext0MT(:, :, :, :, idirC))

      call calcGrFinSH(atoms, vext1AltMT(:, :, :, idirC, 1), grGrTvExt0MTTest)

      do idirR = 1, 3
        do oqn_l = 0, atoms%lmax(1)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
            do imesh = 1, atoms%jri(1)
              vext2AltMT(imesh, lm, idirR, idirC) = r2GrGrTVext0MT(imesh, lm, 1, idirR, idirC) / atoms%rmsh(imesh, 1)**2
              write(4198, '(5(i8),2(es15.8))') idirC, idirR, oqn_l, mqn_m, imesh, vext2AltMT(imesh, lm, idirR, idirC)
              write(4199, '(5(i8),2(es15.8))') idirC, idirR, oqn_l, mqn_m, imesh, grGrTvExt0MTTest(imesh, lm, 1, idirR)
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end do ! idirR
      paPoX = 0.
      paPoY = 0.
      paPoZ = 1.
      resolution = 600
      call plotCompVnVbench( atoms, cell, paPoX, paPoY, paPoZ, idirC + 2, ngdp, resolution, gdp, grGrTvExt0IRFleurDummy(:, :, idirC), &
        & vExt2IR(:, :, idirC, 1), vext2AltMT, vExt2MT(:, :, :, idirC) )
    end do ! idirC

  end subroutine testVexts

  ! is part of testGoldstein modes and deprecated shows that 4pi delta is not part of the integral over Vext2
  subroutine testSimpleRhoVextMTInt( atoms, cell, lathar, ngdp, mlh_atom, nmem_atom, clnu_atom, grVext0MT, vExt1MT, vExt2MT, rho0MT)

    use m_jpConstants, only : fpi, tpi
    use m_types, only : t_atoms, t_cell, t_sphhar
    use m_jpSetupDynMat, only : setupDynMatHF
    use m_jpPotDensHelper, only : calcGrR2FinLH

    implicit none

    ! Type parameters
    type(t_atoms),               intent(in)  :: atoms
    type(t_cell),                intent(in)  :: cell
    type(t_sphhar),              intent(in)  :: lathar

    ! Scalar parameters
    integer,                     intent(in)  :: ngdp

    ! Array parameters
    integer,                     intent(in)  :: mlh_atom(:, 0:, :)
    integer,                     intent(in)  :: nmem_atom(0:, :)
    complex,                     intent(in)  :: clnu_atom(:, 0:, :)
    complex,                     intent(in)  :: grVext0MT(:, :, :, :)
    complex,                     intent(in)  :: vExt1MT(:, :, :, :, :) ! read from hard disk
    complex,                     intent(in)  :: vExt2MT(:, :, :, :) ! keep in memory
    real,                        intent(in)  :: rho0MT(:, 0:, :, :) ! keep in memory

    ! Integer variables
    integer                                  :: imesh
    real                                     :: recR
    integer                                  :: iatom
    integer                                  :: itype
    integer                                  :: ieqat
    integer                                  :: ptsym
    integer                                  :: ilh
    integer                                  :: oqn_l
    integer                                  :: imem
    integer                                  :: mqn_m
    integer                                  :: lm
    integer                                  :: iDtype
    integer                                  :: iDatom
    integer                                  :: iDeqat
    integer                                  :: idir

    ! Array variables
    complex,        allocatable              :: rho0IR(:)
    complex,        allocatable              :: w_vExt2IR(:, :, :, :)
    complex,        allocatable              :: rho1IR(:, :, :)
    complex,        allocatable              :: rho1MT(:, :, :, :, :)
    complex,        allocatable              :: rho1MTSmpl(:, :, :, :, :)
    complex,        allocatable              :: w_vExt1IR(:, :, :)
    real,           allocatable              :: rho0MTSmpl(:, :, :, :)
    complex,        allocatable              :: grVext0MTSmpl(:, :, :, :)
    complex,        allocatable              :: vExt1MTSmpl(:, :, :, :, :)
    complex,        allocatable              :: vExt2MTSmpl(:, :, :, :)
    real,           allocatable              :: r2Rho0MT(:, :, :, :)
    complex,        allocatable              :: r2GrRho0MT(:, :, :, :)
    complex,        allocatable              :: vExt1MTinvGrad(:, :, :, :, :)
    complex                                  :: E2ndOrdII(3, 3)
    complex                                  :: dynMat(3, 3)
    complex                                  :: surfInt(3, 3)

    allocate( rho0MTSmpl(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1), w_vExt1IR(ngdp, 3, atoms%nat),                                      &
      & rho1MTSmpl(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat), rho1IR(ngdp, 3, atoms%nat),                             &
      & w_vExt2IR(ngdp, 3, 3, atoms%nat), rho0IR(ngdp), rho1MT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat) )
    allocate( vExt1MTSmpl(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat),                                              &
      & grVext0MTSmpl(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat), vExt2MTSmpl(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, 3) )
    allocate( r2Rho0MT(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 1) )
    allocate( vExt1MTinvGrad(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat, 3, atoms%nat) )

    rho0MTSmpl(:, 0, 1, 1) = cmplx(sqrt(fpi), 0.)
    w_vExt1IR(:, :, :) = cmplx(0., 0.)
    ! In a first step we derive analytically and the gradient of a constant is zero
    rho1MTSmpl(:, :, :, :, :) = cmplx(0., 0.)
    rho1IR(:, :, :) = cmplx(0., 0.)
    w_vExt2IR(:, :, :, :) = cmplx(0., 0.)
    rho0IR(:) = cmplx(0., 0.)
    E2ndOrdII(:, :) = cmplx(0., 0.)
    vExt1MTSmpl(:, :, :, :, :) = cmplx(0., 0.)
    vExt2MTSmpl(:, :, :, :) = cmplx(0., 0.)
    grVext0MTSmpl(:, :, :, :) = cmplx(0., 0.)

    ! Set vext1 = 1 / r e_r and vext2 1./r^2 in the l = 0 compoenent and the density 0 to check whether the partial works at all
    do imesh = 1, atoms%jri(1)
      recR = 1 / atoms%rmsh(imesh, 1)
      vExt1MTSmpl(imesh, 2, 1, 1, 1) = cmplx(-sqrt(tpi / 3) * recR, 0.)
      vExt1MTSmpl(imesh, 4, 1, 1, 1) = cmplx( sqrt(tpi / 3) * recR, 0.)
      vExt1MTSmpl(imesh, 2, 1, 2, 1) = cmplx(0., -sqrt(tpi / 3) * recR)
      vExt1MTSmpl(imesh, 4, 1, 2, 1) = cmplx(0., -sqrt(tpi / 3) * recR)
      vExt1MTSmpl(imesh, 3, 1, 3, 1) = cmplx(-sqrt(fpi / 3) * recR, 0.)

      grVext0MTSmpl(imesh, 2, 1, 1) = cmplx( sqrt(tpi / 3) * recR, 0.)
      grVext0MTSmpl(imesh, 4, 1, 1) = cmplx(-sqrt(tpi / 3) * recR, 0.)
      grVext0MTSmpl(imesh, 2, 2, 1) = cmplx(0.,  sqrt(tpi / 3) * recR)
      grVext0MTSmpl(imesh, 4, 2, 1) = cmplx(0.,  sqrt(tpi / 3) * recR)
      grVext0MTSmpl(imesh, 3, 3, 1) = cmplx(sqrt(fpi / 3) * recR, 0.)
      recR  = recR / atoms%rmsh(imesh, 1)
      ! We fill only the diagonal and (because the density has no l = 2 component) only the l = 0 component
      vExt2MTSmpl(imesh, 1, 1, 1) = cmplx(sqrt(fpi) / 3 * recR, 0.)
      vExt2MTSmpl(imesh, 1, 2, 2) = cmplx(sqrt(fpi) / 3 * recR, 0.)
      vExt2MTSmpl(imesh, 1, 3, 3) = cmplx(sqrt(fpi) / 3 * recR, 0.)
    end do ! imesh

    call SetupDynMatHF( atoms, cell, lathar, ngdp, ngdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, w_vExt2IR, vExt2MTSmpl, rho1IR, rho1MTSmpl,   &
      & w_vExt1IR, vExt1MTSmpl, rho0MTSmpl, E2ndOrdII, dynMat )

    write(*, '(a)') 'MT Hellmann-Feynman Dynamical Matrix'
    write(*, '(3(2(es16.8,1x),3x))') dynMat(1, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMat(2, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMat(3, :)

    call CalcSurfIntMTDynMat(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MTSmpl, grVext0MTSmpl, surfInt)
    write(*, '(a)') 'MT Surface integral'
    write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)

    write(*, *) 'This test finished'

    ! Now assume additionaly that the density and its gradient is also in the gradient routine
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              r2Rho0MT(imesh, lm, iatom, 1) = rho0MT(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype)**2
            end do ! imesh
          end do ! imem
        end do ! ilh
      end do
    end do
    call calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT(:, :, :, 1), r2GrRho0MT)

    iDatom = 0
    do iDtype = 1, atoms%ntype
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1
        do idir = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do oqn_l = 0, atoms%lmax(itype)
                do mqn_m = -oqn_l, oqn_l
                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                  do imesh = 1, atoms%jri(itype)
                    rho1MT(imesh, lm, iatom, idir, iDatom) = -r2GrRho0MT(imesh, lm, iatom, idir)
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! ieqat
          end do ! itype
        end do ! idir
      end do ! iDeqat
    end do ! iDtype


    call SetupDynMatHF( atoms, cell, lathar, ngdp, ngdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, w_vExt2IR, vExt2MTSmpl, rho1IR, rho1MT, &
      & w_vExt1IR, vExt1MTSmpl, rho0MT, E2ndOrdII, dynMat )

    write(*, '(a)') 'MT Hellmann-Feynman Dynamical Matrix'
    write(*, '(3(2(es16.8,1x),3x))') dynMat(1, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMat(2, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMat(3, :)

    call CalcSurfIntMTDynMat(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT, grVext0MTSmpl, surfInt)
    write(*, '(a)') 'MT Surface integral'
    write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)

    write(*, *) 'This test finished'
    !allocate

    write(*, *) surfInt(1, 1) - dynMat(1, 1)

!    ! Fake vext2 by gradient
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
!                  lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!                  do imesh = 1, atoms%jri(itype)
!                    vExt1MTinvGrad(imesh, lm, iatom, idir, iDatom) = -grVext0(imesh, lm, idir, iatom)
!
!                  end do ! imesh
!                end do ! mqn_m
!              end do ! oqn_l
!            end do ! ieqat
!          end do ! itype
!        end do ! idir
!      end do ! iDeqat
!    end do ! iDtype
!    call calcGrR2FinSH(atoms, r2Rho0MTsh(:, :, :), r2GrRho0MTsh)
!
    call SetupDynMatHF( atoms, cell, lathar, ngdp, ngdp, mlh_atom, nmem_atom, clnu_atom, rho0IR, w_vExt2IR, vExt2MT, rho1IR, rho1MT, &
      & w_vExt1IR, vExt1MT, rho0MT, E2ndOrdII, dynMat )

    write(*, '(a)') 'MT Hellmann-Feynman Dynamical Matrix'
    write(*, '(3(2(es16.8,1x),3x))') dynMat(1, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMat(2, :)
    write(*, '(3(2(es16.8,1x),3x))') dynMat(3, :)

    call CalcSurfIntMTDynMat(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT, grVext0MT, surfInt)
    write(*, '(a)') 'MT Surface integral'
    write(*, '(3(2(es16.8,1x),3x))') surfInt(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfInt(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfInt(3, :)

    ! Surprise, maybe we need a delta function?
    write(*, *) rho0MT(1, 0, 1, 1) * fpi
    !allocate

  end subroutine testSimpleRhoVextMTInt

  ! Test relation with tensor gradient probably somewhere written in the notes
  subroutine TesttensorGradMatElems( atoms, dimens, stars, lathar, kpts, qpts, sym, cell, usdus, results, Veff0, input, ngdp, logUnit,     &
      & rho0IR, gdp, rho0MT, nRadFun, rbas1, rbas2, nmem_atom, mlh_atom, clnu_atom, El, nv, mapGbas, GbasVec, kveclo,   &
      & nobd, z0, iloTable, eig, kpq2kPrVec, memd_atom )

    use m_types, only : t_atoms, t_dimension, t_kpts, t_sphhar, t_cell, t_input, t_stars, t_potential, t_sym, t_usdus, t_results
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat, CalcChannelsGrGrtFlpNat, convertStar2G
    use m_jpSetupDynMat, only : CalcGrVarphiHepsGrtVarphiElem, CalcVecBasfMatElems, calcPsi1HepsPsi1IR
    use m_od_types, only : od_inp, od_sym
    use m_jpConstants, only : iu, Tmatrix, Tmatrix_transposed
    use m_abcof3
    use m_jpSetupDynMatSF, only : CalcSFintIRgradPsiHepsPsi, CalcSFintMTgradPsiHepsPsi
    use m_jpPotDensHelper, only : warpIRPot
    use m_jpSternhHF, only : CalcMEPotIR
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_dimension),              intent(in) :: dimens
    type(t_stars),                  intent(in) :: stars
    type(t_sphhar),                 intent(in) :: lathar
    type(t_kpts),                   intent(in) :: kpts
    type(t_kpts),                   intent(in) :: qpts
    type(t_sym),                    intent(in) :: sym
    type(t_cell),                   intent(in) :: cell
    type(t_usdus),                  intent(in) :: usdus
    type(t_results),                intent(in) :: results
    type(t_potential),              intent(in) :: Veff0
    type(t_input),                  intent(in) :: input

    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: logUnit
    integer,                        intent(in) :: memd_atom

    ! Array paramters
    complex,                        intent(in) :: rho0IR(:, :)
    integer,                        intent(in) :: gdp(:, :)
    real,                           intent(in) :: rho0MT(:, :, :, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:,:,0:,:,:)
    real,                           intent(in) :: rbas2(:,:,0:,:,:)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    complex,                        intent(in) :: clnu_atom(:,0:,:)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: nobd(:, :)
    MCOMPLEX,                       intent(in) :: z0(:, :, :, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    real,                           intent(in) :: eig(:, :, :)

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods

    ! Scalar variables
    integer                                    :: iatom
    integer                                    :: ieqat
    integer                                    :: itype
    integer                                    :: oqn_l
    integer                                    :: iradf
    integer                                    :: imesh
    integer                                    :: nRadFunMax
    integer                                    :: mqn_m2prR
    integer                                    :: lmp
    integer                                    :: mqn_m
    integer                                    :: ichan
    integer                                    :: lmpMax
    integer                                    :: ptsym
    integer                                    :: ilh
    integer                                    :: lm_pre
    integer                                    :: imem
    integer                                    :: lm
    integer                                    :: mqn_m2PrC
    integer                                    :: ikpt
    integer                                    :: nmat
    integer                                    :: iband
    integer                                    :: iDatom
    integer                                    :: iDtype
    integer                                    :: iDeqat
    integer                                    :: pMaxLocal
    integer                                    :: idirR
    integer                                    :: idirC
    integer                                    :: iBas

    ! Array variables
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    real,              allocatable             :: delrGrVarphiCh1(:, :, :, :)
    real,              allocatable             :: delrGrVarphiCh2(:, :, :, :)
    integer,           allocatable             :: grGrtVarphiChLout(:, :)
    integer,           allocatable             :: grGrtVarphiChMout(:, :, :)
    real,              allocatable             :: grGrtVarphiCh1(:, :, :, :, :)
    real,              allocatable             :: grGrtVarphiCh2(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    real,              allocatable             :: r2(:)
    real,              allocatable             :: grVarphiGrtVarphiNatNat(:, :, :, :, :)
    complex,           allocatable             :: grVarphiHpreGrtVarphiNatNat(:, :, :, :, :)
    real,              allocatable             :: grGrtVarphiVarphiNatNat(:, :, :, :, :)
    complex,           allocatable             :: grGrtVarphiHvarphiNatNat(:, :, :, :, :)
    complex,           allocatable             :: grVarphiGrtVeff0SphVarphiNat(:, :, :, :, :)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    complex,           allocatable             :: ab0cofScl(:)
    complex,           allocatable             :: grVarphiGrtVeff0SphPsiNat(:)
    complex,           allocatable             :: grVarphiHepsGrtPsiNatNat(:)
    complex,           allocatable             :: z1nG(:, :)
    complex,           allocatable             :: gradGradTz0(:, :, :)
    complex,           allocatable             :: grGrtVarphiHPsiNatNat(:)
    complex,           allocatable             :: grGrtVarphiPsiNatNat(:)
    complex,           allocatable             :: surfIntIR(:, :)
    complex,           allocatable             :: surfIntMT(:, :)
    complex                                    :: grPsiHepsGrtPsiNatNat(3, 3)
    complex                                    :: grPsiGrtVeff0SphPsiNat(3, 3)
    complex                                    :: grPsiHepsGrtPsiNat(3, 3)
    complex                                    :: grPsiHepsGrtPsi(3, 3)
    complex                                    :: grGrtPsiHepsPsiNat(3, 3)
    complex                                    :: grGrtPsiHepsPsi(3, 3)
    real                                       :: kExt(3)
    real                                       :: gbasExt(3)
    complex                                    :: z1VarphiHepsVarphiz1IR(3, 3)
    complex                                    :: z0gradGradTVarphiHepsVarphiz0IR(3, 3)
    complex                                    :: grGrtPsiHepsPsiNatNat(3, 3)
    complex                                    :: testSum(1:3, 1:3)
    complex                                    :: convIntgrl(1:3, 1:3)

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( grGrtVarphiChLout(4, 0:atoms%lmaxd), grGrtVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1, -1:1), &
            & grGrtVarphiCh1(atoms%jmtd, 4, lmpMax, -1:1, -1:1), grGrtVarphiCh2(atoms%jmtd, 4, lmpMax, -1:1, -1:1) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( z1nG(dimens%nbasfcn, 3), gradGradTz0(dimens%nbasfcn, 3, 3) )
    allocate( grVarphiGrtVarphiNatNat(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grVarphiHpreGrtVarphiNatNat(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( grGrtVarphiVarphiNatNat(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grGrtVarphiHvarphiNatNat(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( grVarphiGrtVeff0SphVarphiNat(lmpMax, lmpMax, -1:1, 1:3, atoms%nat) )
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
            & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( grVarphiHepsGrtPsiNatNat(lmpMax), grVarphiGrtVeff0SphPsiNat(lmpMax) )
    allocate( grGrtVarphiHPsiNatNat(lmpMax), grGrtVarphiPsiNatNat(lmpMax) )
    allocate( r2(atoms%jmtd) )
    allocate( surfIntIR(3, 3), surfIntMT(3, 3) )

    grVarphiGrtVarphiNatNat(:, :, :, :, :) = 0.
    grVarphiHpreGrtVarphiNatNat(:, :, :, :, :) = 0.
    grGrtVarphiVarphiNatNat = cmplx(0., 0.)
    grGrtVarphiHVarphiNatNat = cmplx(0., 0.)
    grVarphiGrtVeff0SphVarphiNat = cmplx(0., 0.)
    testSum(:, :) = cmplx(0., 0.)
    r2(:) = 0.

    convIntgrl(:, :) = cmplx(0., 0.)

    !is the old version within this module
!    call TestDynMatPulInt( atoms, sym, dimens, lathar, stars, input, usdus, Veff0, kpts, cell, results, ngdp, logUnit, rho0IR, rho0MT,      &
!      & Veff0%vr, gdp, clnu_atom, nmem_atom, mlh_atom, nRadFun, rbas1, rbas2, nv, mapGbas, gbasVec, nobd, z0, kpq2kPrVec, El,      &
!      & kveclo, iloTable, ConvIntgrl )

    call CalcGrVeffGrtRhoInt( atoms, cell, lathar, dimens, stars, Veff0, input, ngdp, memd_atom, clnu_atom, nmem_atom, mlh_atom,      &
      & rho0IR, gdp, rho0MT, ConvIntgrl )
    ! todo
    ! NOTE: lmax + 1 is missing but not relevant for Neon

    iatom = 0
    do itype = 1, atoms%ntype

      do imesh = 1, atoms%jri(itype)
        r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      end do ! imesh

      ! Prepare the basis functions, the gradient of the basis function and the tensor product of gradients applied to the basis
      ! functions
      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      delrVarphi1(:, :, :) = 0.
      delrVarphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(itype)
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of
            ! the Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
          end do ! imesh
          ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
          call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
          call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
        end do ! iradf
      end do ! oqn_l

      ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching
      ! coefficients have no spatial dependence) and determing its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout,&
                                                                                                        & grVarphiCh1, grVarphiCh2 )
    write(*, *) 'gr finished'
      ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times.
      delrGrVarphiCh1(:, :, :, :) = 0.
      delrGrVarphiCh2(:, :, :, :) = 0.
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              do ichan = 1, 2
                call Derivative( grVarphiCh1(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
                call Derivative( grVarphiCh2(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
              end do ! ichan
            end do !iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR

      ! Calculate tensorial product of gradients applied onto MT basis function expanded in spherical harmonics
      grGrtVarphiChLout(:, :) = 0
      grGrtVarphiChMout(:, :, :) = 0
      grGrtVarphiCh1(:, :, :, :, :) = 0.
      grGrtVarphiCh2(:, :, :, :, :) = 0.
      call CalcChannelsGrGrtFlpNat( atoms, itype, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, &
                          & delrGrVarphiCh1, delrGrVarphiCh2, grGrtVarphiChLout, grGrtVarphiChMout, grGrtVarphiCh1, grGrtVarphiCh2 )
    write(*, *) 'grGrt finished'

      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1


        call CalcSFintIRgradPsiHepsPsi( atoms, stars, cell, dimens, kpts, results, Veff0, ngdp, iatom, itype, nobd, eig, gBasVec, &
                                                                                                       & mapGbas, nv, gdp, z0, surfIntIR )


        r2grVeff0SphVarphi = cmplx(0., 0.)
        vEff0NsphGrVarphi = cmplx(0., 0.)
        ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
            ! maybe construct a pointer and run only over them to make it memory efficient.
            do imesh = 1, atoms%jri(itype)
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + Veff0%vr(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh

write(*, *) 'calchngrv0varphiINit'
        ! Prepare MT matrix elements for the test
        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1),          &
          & clnu_atom, nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi,  &
          & r2grVeff0SphVarphi )
write(*, *) 'calchngrv0varphi'

        call CalcSFintMTgradPsiHepsPsi( atoms, cell, kpts, sym, dimens, results, usdus, itype, iatom, nRadFun, &
          & lmpMax, nobd, nv, gBasVec, mapGbas, kveclo, z0, eig, hVarphi, iloTable, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, varphi1, varphi2, surfIntMT )

        do mqn_m2PrC = -1, 1
        write(*, *) mqn_m2PrC
          ! Calculate <grad gradT varphi| H - eps | varphi >
          call CalcVecBasfMatElems( atoms, itype, 4, nRadFun, r2, grGrtVarphiChLout, grGrtVarphiChMout(:, :, mqn_m2PrC),  &
            & varphi1, varphi2, grGrtVarPhiCh1(:, :, :, :, mqn_m2PrC), grGrtVarPhiCh2(:, :, :, :, mqn_m2PrC), hVarphi,        &
            & grGrtVarphiVarphiNatNat(:, :, :, mqn_m2PrC, itype), grGrtVarphiHVarphiNatNat(:, :, :, mqn_m2PrC, iatom) )
        write(*, *) mqn_m2PrC

          ! Calculate <grad varphi| H - eps | gradT varphi >
          call CalcGrVarphiHepsGrtVarphiElem( atoms, itype, iatom, mqn_m2PrC, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, r2, &
            & grVarPhiCh1, grVarPhiCh2, El, lmpT, vEff0NsphGrVarphi, r2grVeff0SphVarphi, grVarphiGrtVarphiNatNat,                  &
            & grVarphiHpreGrtVarphiNatNat, grVarphiGrtVeff0SphVarphiNat)
        end do ! mqn_m2PrC
      end do ! ieqat

    end do ! itype

    write(*, *) 'finished preparations'
    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1
    allocate( ab0cofScl(lmpMax) )

    grPsiHepsGrtPsi(:, :) = cmplx(0., 0.)

    do ikpt = 1, kpts%nkpt

      ! Calculate small matching coefficients to be multiplied onto the matrix element
      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbasVec(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        & gbasVec(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbasVec(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

      do iband = 1, nobd(ikpt, 1)
        iDatom = 0
        do iDtype = 1, atoms%ntype
          do iDeqat = 1, atoms%neq(iDtype)
            iDatom = iDatom + 1

            z1nG(:, :) = cmplx(0., 0.)
            kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
            do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
              gbasExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, mapGbas(iBas, ikpt, 1)))
              do idirC = 1, 3
                z1nG(iBas, idirC) = iu * ( kExt(idirC) + gbasExt(idirC) ) * z0(iBas, iband, ikpt, 1)
                do idirR = 1, 3
                  gradGradTz0(iBas, idirR, idirC) = - ( kExt(idirR) + gbasExt(idirR) ) * ( kExt(idirC) + gbasExt(idirC) )          &
                                                                                                          * z0(iBas, iband, ikpt, 1)
                end do ! idirR
              end do ! idir
            end do ! iBas
            ! Calculation of the IR quantities

            z1VarphiHepsVarphiz1IR(:, :) = cmplx(0., 0.)
            z0gradGradTVarphiHepsVarphiz0IR(:, :) = cmplx(0., 0.)
            do idirR = 1, 3
              do idirC = 1, 3
                call calcPsi1HepsPsi1IR( kpts, qpts, stars, dimens, cell, gbasVec(:, mapGbas(:nv(1, ikpt), ikpt, 1)), nv, &
                  & z1nG(:nv(1, ikpt), idirR), z1nG(:nv(1, ikpt), idirC), nmat, ikpt, 1, ikpt, kpq2kPrVec, &
                  & Veff0%vpw, eig, iband, z1VarphiHepsVarphiz1IR(idirR, idirC) )

                !todo do not give the parameters with limits!!!!!
                call calcPsi1HepsPsi1IR( kpts, qpts, stars, dimens, cell, gbasVec(:, mapGbas(:nv(1, ikpt), ikpt, 1)), nv, &
                  & gradGradTz0(:nv(1, ikpt), idirR, idirC), z0(:nv(1, ikpt), iband, ikpt, 1), nmat, ikpt, 1, ikpt, kpq2kPrVec, &
                  & Veff0%vpw, eig, iband, z0gradGradTVarphiHepsVarphiz0IR(idirR, idirC) )
              end do ! idirC
            end do ! idirR

            lm  = 0
            lmp = 0
            ! Calculate large matching coefficients to be multiplied onto the matrix element
            ab0cofScl(:) = cmplx(0., 0.)
            do oqn_l = 0, atoms%lmax(iDtype)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtype)
                ! p = 1
                ab0cofScl(lmp + 1) = iu**oqn_l * dot_product(conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom))
                ! p = 2
                ab0cofScl(lmp + 2) = iu**oqn_l * dot_product(conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom))
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) &
                    & + iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) &
                    & + iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                end do ! iradf
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do ! oqn_l
            ! Multiply large matching coefficients onto the basis function matrix elements
            grPsiHepsGrtPsiNatNat(:, :) = cmplx(0., 0.)
            grVarphiHepsGrtPsiNatNat = cmplx(0., 0.)
            grVarphiGrtVeff0SphPsiNat(:) = cmplx(0., 0.)
            grPsiGrtVeff0SphPsiNat(:, :) = cmplx(0., 0.)
            do idirC = -1, 1
              do idirR = -1, 1
                grVarphiHepsGrtPsiNatNat(1:lmpT(iDtype)) = &
                  &   matmul(grVarphiHpreGrtVarphiNatNat(1:lmpT(iDtype), 1:lmpT(iDtype), idirR, idirC, iDatom),                    &
                                                                                                      & ab0cofScl(1:lmpT(iDtype)) )&
                  & - eig(iband, ikpt, 1) * matmul(grVarphiGrtVarphiNatNat(1:lmpT(iDtype), 1:lmpT(iDtype), idirR, idirC, iDtype)   &
                                                                                                    &  , ab0cofScl(1:lmpT(iDtype)) )

                grPsiHepsGrtPsiNatNat(idirR + 2, idirC + 2) = dot_product(ab0cofScl(1:lmpT(iDtype)),                              &
                                                                                        & grVarphiHepsGrtPsiNatNat(1:lmpT(iDtype)) )

                grVarphiGrtVeff0SphPsiNat(1:lmpT(iDtype)) = &
                  & matmul( grVarphiGrtVeff0SphVarphiNat(1:lmpT(iDtype), 1:lmpT(iDtype), idirR, idirC + 2, iDatom),                &
                                                                                                       & ab0cofScl(1:lmpT(iDtype)) )
                grPsiGrtVeff0SphPsiNat(idirR + 2, idirC + 2) = dot_product(ab0cofScl(1:lmpT(iDtype)), &
                                                                                       & grVarphiGrtVeff0SphPsiNat(1:lmpT(iDtype)) )

                ! Calculation of the last line in A.50 PhD thesis A. Klueppelberg
                grGrtVarphiHPsiNatNat(:) = cmplx(0.0, 0.0)
                grGrtVarphiHPsiNatNat(1:lmpT(iDtype)) =                                                                            &
                    & matmul(grGrtVarphiHvarphiNatNat(:lmpT(iDtype), :lmpT(iDtype), idirR, idirC, iDatom), ab0cofScl(:lmpT(iDtype)))

                ! Precalculation of the last line in A.50 PhD thesis A. Klueppelberg.
                grGrtVarphiPsiNatNat(:) = cmplx(0.0, 0.0)
                grGrtVarphiPsiNatNat(1:lmpT(iDtype)) = &
                    & matmul( grGrtVarphiVarphiNatNat(1:lmpT(iDtype), 1:lmpT(iDtype), idirR, idirC, iDtype), &
                                                                                                       & ab0cofScl(1:lmpT(iDtype)) )
                grGrtPsiHepsPsiNatNat(idirR + 2, idirC + 2) =                                                                      &
                            & dot_product(ab0cofScl(:lmpT(iDtype)), grGrtVarphiHPsiNatNat(1:lmpT(iDtype)))&
                            & - eig(iband, ikpt, 1) * dot_product(ab0cofScl(:lmpT(iDtype)), grGrtVarphiPsiNatNat(1:lmpT(iDtype)) )
              end do ! idirR
            end do ! idirC
            ! Transform from natural into cartesian coordinates. Note, that grPsiGrtVeff0SphPsiNat, already contains the gradient of
            ! the spherical effective potential transformed into cartesian coordinates, hence the respective matrix element only
            ! contains one multiplication with the T-matrix.
            grPsiHepsGrtPsiNat(1:3, 1:3) = cmplx(0., 0.)
            grPsiHepsGrtPsiNat(1:3, 1:3) = matmul(grPsiHepsGrtPsiNatNat(1:3, 1:3), Tmatrix_transposed(1:3, 1:3))
!            grPsiHepsGrtPsiNat(1:3, 1:3) = cmplx(0., 0.)
            grPsiHepsGrtPsiNat(1:3, 1:3) = grPsiHepsGrtPsiNat(1:3, 1:3) - grPsiGrtVeff0SphPsiNat(1:3, 1:3)
            grPsiHepsGrtPsi(1:3, 1:3) = matmul(conjg(Tmatrix(1:3, 1:3)), grPsiHepsGrtPsiNat(1:3, 1:3))

            grGrtPsiHepsPsiNat(1:3, 1:3) = matmul(conjg(Tmatrix(1:3, 1:3)), grGrtPsiHepsPsiNatNat(1:3, 1:3))
            grGrtPsiHepsPsi(1:3, 1:3) = matmul( grGrtPsiHepsPsiNat(1:3, 1:3), conjg(Tmatrix_transposed(1:3, 1:3)) )

            testSum(1:3, 1:3) = testSum(1:3, 1:3) + results%w_iks(iband, ikpt, 1)                                              &
                                                            & * (                                                                  &
                                                            &   z1VarphiHepsVarphiz1IR(1:3, 1:3)                                   &
                                                            & + z0gradGradTVarphiHepsVarphiz0IR(1:3, 1:3)                          &
                                                            & + grPsiHepsGrtPsi(1:3, 1:3)                                          &
                                                            & + grGrtPsiHepsPsi(1:3, 1:3)                                          &
                                                            )
          end do ! iDeqat
        end do ! iDtype
      end do ! iband
    end do ! ikpt


!    write(*, '(a)') 'grVarphiGrtvarphiIR for last k-point'
!    write(*, '(3(2(es16.8,1x),3x))') z1VarphiHepsVarphiz1IR(1, :)
!    write(*, '(3(2(es16.8,1x),3x))') z1VarphiHepsVarphiz1IR(2, :)
!    write(*, '(3(2(es16.8,1x),3x))') z1VarphiHepsVarphiz1IR(3, :)
!
!    write(*, '(a)') 'grGrtVarphiVarphiIR for last k-point'
!    write(*, '(3(2(es16.8,1x),3x))') z0gradGradTVarphiHepsVarphiz0IR(1, :)
!    write(*, '(3(2(es16.8,1x),3x))') z0gradGradTVarphiHepsVarphiz0IR(2, :)
!    write(*, '(3(2(es16.8,1x),3x))') z0gradGradTVarphiHepsVarphiz0IR(3, :)
!
!    write(*, '(a)') 'grVarphiGrtvarphMT for last k-point'
!    write(*, '(3(2(es16.8,1x),3x))') grPsiHepsGrtPsi(1, :)
!    write(*, '(3(2(es16.8,1x),3x))') grPsiHepsGrtPsi(2, :)
!    write(*, '(3(2(es16.8,1x),3x))') grPsiHepsGrtPsi(3, :)
!
!    write(*, '(a)') 'grGrtVarphiVarphi for last k-point'
!    write(*, '(3(2(es16.8,1x),3x))') grGrtPsiHepsPsi(1, :)
!    write(*, '(3(2(es16.8,1x),3x))') grGrtPsiHepsPsi(2, :)
!    write(*, '(3(2(es16.8,1x),3x))') grGrtPsiHepsPsi(3, :)
    ! The gradient of the spherical potential in grad varphi varphi is equal to the integral, but we have to leave out the IR and the nonspherical components
    ! We might test the nonspherical componenet analogously to the core terms recasting!

    write(*, '(a)') 'grGrtPsiHepsPsi + grPsiHepsGrtPsi'
    write(*, '(3(2(es16.8,1x),3x))') testSum(1, :)
    write(*, '(3(2(es16.8,1x),3x))') testSum(2, :)
    write(*, '(3(2(es16.8,1x),3x))') testSum(3, :)
!
    write(*, '(a)') 'rho GrVeff'
    write(*, '(3(2(es16.8,1x),3x))') 0.25 * convIntgrl(1, :)
    write(*, '(3(2(es16.8,1x),3x))') 0.25 * convIntgrl(2, :)
    write(*, '(3(2(es16.8,1x),3x))') 0.25 * convIntgrl(3, :)
    testSum(1:3, 1:3) = testSum(1:3, 1:3) + 0.25 * convIntgrl(1:3, 1:3)

    write(*, '(a)') 'surfIntIR'
    write(*, '(3(2(es16.8,1x),3x))') surfIntIR(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntIR(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntIR(3, :)

    write(*, '(a)') 'surfIntMT'
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(3, :)

    write(*, '(a)') 'surfIntMT + surfIntIR'
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(1, :) + surfIntIR(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(2, :) + surfIntIR(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(3, :) + surfIntIR(3, :)

    write(*, '(a)') 'grGrtPsiHepsPsi + grPsiHepsGrtPsi + rho GrVeff'
    write(*, '(3(2(es16.8,1x),3x))') testSum(1, :)
    write(*, '(3(2(es16.8,1x),3x))') testSum(2, :)
    write(*, '(3(2(es16.8,1x),3x))') testSum(3, :)


  end subroutine TesttensorGradMatElems

  subroutine TestDynMatPulInt( atoms, sym, dimens, lathar, stars, input, usdus, Veff0, kpts, cell, results, ngdp, logUnit, rho0IRst, rho0MTlh, Veff0Mtlh, gdp, clnu_atom, nmem_atom, mlh_atom, nRadFun, rbas1, rbas2, nv, ilst, GbasVec, nobd, z, kpq2kPrVec, El, kveclo, iloTable, ConvIntgrl)

    use m_types, only : t_atoms, t_sym, t_dimension, t_sphhar, t_stars, t_potential, t_kpts, t_cell, t_input, t_usdus, t_results
    use m_jpConstants, only : iu, Tmatrix
    use mod_juPhonUtils, only : convertStar2G, Derivative, CalcChannelsGrFlpNat
    use m_jpSternhHF, only : calcMEPotIR
    use m_jpPotDensHelper, only : warpIRPot, calcGrR2FinLH
    use m_jpSetupDynMat, only : EvalIntRho1Veff1, CalcGrVarphiHepsGrtVarphiElem
    use m_od_types, only : od_inp, od_sym
    use m_abcof3
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi, CalcFnsphVarphi

    use m_jpTestPotential, only : checkjuPhPots

    use m_intgr

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_sym),                    intent(in) :: sym
    type(t_dimension),              intent(in) :: dimens
    type(t_sphhar),                 intent(in) :: lathar
    type(t_stars),                  intent(in) :: stars
    type(t_potential),              intent(in) :: Veff0
    type(t_kpts),                   intent(in) :: kpts
    type(t_cell),                   intent(in) :: cell
    type(t_input),                  intent(in) :: input
    type(t_usdus),                  intent(in) :: usdus
    type(t_results),                intent(in) :: results

    ! Scalar parameters
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: logUnit

    ! Array parameter
    complex,                        intent(out):: convIntgrl(:, :)

    ! Array parameters
    complex,                        intent(in) :: rho0IRst(:, :)
    real,                           intent(in) :: rho0Mtlh(:, 0:, :, :)
    real,                           intent(in) :: Veff0Mtlh(:, 0:, :, :)
    integer,                        intent(in) :: gdp(:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    integer,                        intent(in) :: nRadFun(0:, :)
    real,                           intent(in) :: rbas1(:,:,0:,:,:)
    real,                           intent(in) :: rbas2(:,:,0:,:,:)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    integer,                        intent(in) :: GbasVec(:, :)
    integer,                        intent(in) :: nobd(:, :)
    MCOMPLEX,                       intent(in) :: z(:, :, :, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in)    :: iloTable(:, 0:, :)

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods

    ! Scalar variables
    integer                                    :: ikpt
    integer                                    :: iG
    integer                                    :: idirR
    integer                                    :: idirC
    integer                                    :: itype
    integer                                    :: ilh
    integer                                    :: imesh
    integer                                    :: oqn_l
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: iqpt
    integer                                    :: iatom
    integer                                    :: ieqat
    integer                                    :: lm_pre
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: iradf
    integer                                    :: iBas
    integer                                    :: iband
    integer                                    :: nmat
    integer                                    :: ispin
    integer                                    :: mqn_m2prC
    integer                                    :: lmp
    integer                                    :: pMaxLocal
    integer                                    :: idir
    integer                                    :: mqn_m2PrR


    integer         :: lmp1Pr
    integer         :: lm1Pr
    integer         :: oqn_l1Pr
    integer         :: mqn_m1Pr
    integer         :: mqn_m3Pr
    integer         :: iradf1Pr
    integer         :: ichanPr
    integer         :: oqn_l3Pr
    integer         :: lm3Pr
    real, allocatable            :: intgrdR(:)
    real, allocatable            :: intgrdI(:)
    real :: integralR
    real :: integralI

    ! Array variables
    integer,           allocatable             :: lmpT(:)
    complex,           allocatable             :: grVeff0IR(:, :)
    complex,           allocatable             :: grRho0IR(:, :)
    complex,           allocatable             :: vEff0Pw(:)
    complex,           allocatable             :: rho0IRpw(:)
    complex,           allocatable             :: w_grVeff0IR(:, :)
    real,              allocatable             :: r2Rho0MTlh(:, :, :, :)
    real,              allocatable             :: r2Veff0MTlh(:, :, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    complex,           allocatable             :: rho1IRContainer(:, :, :)
    complex,           allocatable             :: w_vEff1IRContainer(:, :, :)
    complex,           allocatable             :: rho1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: vEff1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: r2GrRho0MTsh(:, :, :, :)
    complex,           allocatable             :: r2GrVeff0Mtsh(:, :, :, :)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: r2(:)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    complex,           allocatable             :: ikpGz0(:, :, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable             :: grVarphiGrtVeff0SphVarphi(:, :, :, :, :)
    real,              allocatable             :: vEff0MtLhDummy(:, :, :)
    real,              allocatable             :: grVarphiGrtVarphi(:, :, :, :, :)
    complex,           allocatable             :: grVarphiHpreGrtVarphi(:, :, :, :, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    complex,           allocatable             :: altIntgrlIRband(:, :)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: ab0cofScl(:)
    complex,           allocatable             :: grVarphiGrtVeff0PsiNat(:)
    complex                                    :: grPsiGrtVeff0PsiNat(-1:1, 3)
    complex                                    :: grPsiGrtVeff0PsiMt(3, 3)
    real                                       :: Gext(3)
    real                                       :: kExt(3)
    complex                                    :: altIntgrlIR(3, 3)
    complex                                    :: altIntgrl(3, 3)
    integer :: ii, jj
!!!!!!!!!!!!!!!!!!!!!!
complex, allocatable :: grVarphiGrtVeff0SphVarphiMeshNat(:, :, :, :, :, :)
complex, allocatable :: grPsiGrtVeff0PsiMtMesh(:, :, :)
complex, allocatable :: grPsiGrtVeff0PsiNatMesh(:, :, :)
complex, allocatable :: grVarphiGrtVeff0PsiNatMesh(:)
complex, allocatable :: integrandInt(:, :, :)
real, allocatable :: integralRMesh(:), integralIMesh(:)
complex, allocatable :: grVeff0MTCont(:, :, :, :)

ConvIntgrl = cmplx(0., 0.)

allocate(intgrdR(atoms%jmtd), intgrdI(atoms%jmtd))
allocate(integralRMesh(atoms%jmtd), integralIMesh(atoms%jmtd))
    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
allocate(grPsiGrtVeff0PsiMtMesh(atoms%jmtd, 3, 3), grPsiGrtVeff0PsiNatMesh(atoms%jmtd, -1:1, 3), grVarphiGrtVeff0PsiNatMesh(lmpMax))
intgrdR(:) = 0.
intgrdI(:) = 0.
integralRMesh(:) = 0.
integralIMesh(:) = 0.
grPsiGrtVeff0PsiMtMesh(:, :, :) = cmplx(0., 0.)
grPsiGrtVeff0PsiNatMesh(:, :, :) = cmplx(0., 0.)
grVarphiGrtVeff0PsiNatMesh(:) = cmplx(0., 0.)

    altIntgrl(:, :) = cmplx(0., 0.)

!!!!!!
    ! This test is only applicable for q = 0
    iqpt = 1

    ! Convert from star expansion coefficients to plane-wave expansion coefficients
    allocate( vEff0Pw(ngdp), rho0IRpw(ngdp) )
    vEff0Pw(:) = cmplx(0., 0.)
    rho0IRpw(:) = cmplx(0., 0.)
    call convertStar2G( Veff0%vpw_uw(:, 1), vEff0Pw(:), stars, ngdp, gdp )
    call convertStar2G( rho0IRst(:, 1), rho0IRpw(:), stars, ngdp, gdp )

    ! Perform analytical gradient of unperturbed density and effective potential in the interstitial region
    allocate( grVeff0IR(ngdp, 3), grRho0IR(ngdp, 3) )
    grVeff0IR(:, :) = cmplx(0., 0.)
    grRho0IR(:, :) = cmplx(0., 0.)
    do iG = 1, ngdp
      Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
      grVeff0IR(iG, 1:3) = iu * Gext(1:3) * vEff0Pw(iG)
      grRho0IR(iG, 1:3)  = iu * Gext(1:3) * rho0IRpw(iG)
    end do ! iG

    ! Warp interstitial second-order external potential
    allocate(w_grVeff0IR(ngdp, 3))
    w_grVeff0IR = cmplx(0.0, 0.0)
    do idirC = 1, 3
      call warpIRPot(stars, ngdp, idirC, gdp, grVeff0IR, w_grVeff0IR(:, idirC))
    end do ! idirC

    ! Gradient of density and effective potential in MT
    allocate( r2Rho0MTlh( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
    allocate( r2Veff0MTlh( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    allocate( grVeff0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    r2Rho0MTlh(:, :, :, :) = 0.
    r2Veff0MTlh(:, :, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)
    grVeff0MT(:, :, :, :) = cmplx(0., 0.)

    ! Read in valence density form FLEUR
    read(1040) r2Rho0MtLh

    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          r2Veff0MTlh(imesh, ilh, itype, 1) = vEff0MTlh(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
        end do
      end do
    end do

    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MTlh(:, :, :, 1), r2GrRho0MTsh )
    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0MTlh(:, :, :, 1), r2GrVeff0MTsh )

    do idirC = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                ! For reasons of performance, this calculation is done within the idirC loop
                grRho0MT(imesh, lm, iatom, idirC) = r2GrRho0MTsh(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
                grVeff0MT(imesh, lm, iatom, idirC) = r2GrVeff0MTsh(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idirC


    ! Fill container arrays with gradient coefficients so that they can be processed by EvelIntRho1Veff1
    allocate( rho1IRContainer(ngdp, 3, atoms%nat), w_vEff1IRContainer(ngdp, 3, atoms%nat),  &
      & rho1MTContainer( atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3, atoms%nat ), &
      & vEff1MTContainer( atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3, atoms%nat ) )

    rho1IRContainer(:, :, :) = cmplx(0., 0.)
    w_vEff1IRContainer(:, :, :) = cmplx(0., 0.)
    rho1MTContainer(:, :, :, :, :) = cmplx(0., 0.)
    vEff1MTContainer(:, :, :, :, :) = cmplx(0., 0.)

!    ! todo think about this for more than one atom


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOTE: NO IR FOR THIS TEST TO WORK. THIS WAS ADDED LATER TO TEST SOMETHONG DIFFERENT

    do idirC = 1, 3
      do iG = 1, ngdp
        rho1IRContainer(iG, idirC, 1) = -grRho0IR(iG, idirC)
        write(4004, '(2i8,2f15.8)') iG, idirC, -rho1IRContainer(iG, idirC, 1)
        w_vEff1IRContainer(iG, idirC, 1) = -w_grVeff0IR(iG, idirC)
        write(4005, '(2i8,2f15.8)') iG, idirC, -w_vEff1IRContainer(iG, idirC, 1)
      end do ! iG
    end do ! idirC

    allocate(grVeff0MTCont(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat))
    grVeff0MTCont(:, :, :, :) = cmplx(0., 0.)
    !todo does not work for more than one atom
    do idirC = 1, 3
      ! There should be an itype
      do oqn_l = 0, atoms%lmax(1)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          do imesh = 1, atoms%jri(1)
            rho1MTContainer(imesh, lm, 1, idirC, 1) = -grRho0MT(imesh, lm, 1, idirC)
            write(4006, '(3i8,2f15.8)') imesh, lm, idirC, -rho1MTContainer(imesh, lm, 1, idirC, 1)
            vEff1MTContainer(imesh, lm, 1, idirC, 1) = -grVeff0MT(imesh, lm, 1, idirC)
            grVeff0MTCont(imesh, lm, idirC, 1) = grVeff0MT(imesh, lm, 1, idirC)
            write(4007, '(3i8,2f15.8)') imesh, lm, idirC, -vEff1MTContainer(imesh, lm, 1, idirC, 1)
          end do ! imesh
        end do ! mqn_m
      end do ! oqn_l
    end do ! idirC

    !todo continious version of grVeff instead of two gradients
    ! todo when calculating continious version is there the same integral used
    call checkjuPhPots(30, atoms, ngdp, gdp, cell, lathar, sym, grVeff0IR, grVeff0MTCont, logUnit)
    call checkjuPhPots(30, atoms, ngdp, gdp, cell, lathar, sym, w_grVeff0IR, grVeff0MTCont, logUnit)

    ! Calculate the integral in the conventional way
    call EvalIntRho1Veff1(atoms, cell, rho1IRContainer, rho1MTContainer, w_vEff1IRContainer, vEff1MTContainer, ConvIntgrl)
    write(*, '(a)') 'Matrix of conventional integrals'
    write(*, '(3(2(es16.8,1x),3x))') 0.25* convIntgrl(1, :)
    write(*, '(3(2(es16.8,1x),3x))') 0.25* convIntgrl(2, :)
    write(*, '(3(2(es16.8,1x),3x))') 0.25* convIntgrl(3, :)

    return

    ! Setup zBra and zKet for evaluation of respective benchmark braket

    ! Quantities for initialization
    !allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( ikpGz0(dimens%nbasfcn, dimens%neigd, 3, kpts%nkpt) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( grVarphiGrtVeff0SphVarphi(lmpMax, lmpMax, -1:1, 1:3, atoms%nat), vEff0MtLhDummy( atoms%jmtd, 0:lathar%nlhd, atoms%ntype ) )
    allocate( r2(atoms%jmtd) )
    allocate( grVarphiGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grVarphiHpreGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( altIntgrlIRband(maxval(nobd), maxval(nobd)) )
    allocate( ab0cofScl(lmpMax) )
    allocate(grVarphiGrtVeff0PsiNat(lmpMax))

    grVarphiCh1(:, :, :, :) = 0.
    grVarphiCh2(:, :, :, :) = 0.
    grVarphiChLout(:, :) = 0
    grVarphiChMout(:, :) = 0
    ikpGz0(:, :, :, :) = cmplx(0., 0.)
    vEff0MtSpH = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiGrtVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    vEff0MtLhDummy(:, :, :) = cmplx(0., 0.)
    r2(:) = 0.
    grVarphiGrtVarphi(:, :, :, :, :) = 0.
    grVarphiHpreGrtVarphi(:, :, :, :, :) = cmplx(0., 0.)
    altIntgrlIRband(:, :) = cmplx(0., 0.)
    grVarphiGrtVeff0PsiNat(:) = cmplx(0., 0.)
    ab0cofScl(:) = cmplx(0., 0.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(grVarphiGrtVeff0SphVarphiMeshNat(atoms%jmtd, lmpMax, lmpMax, -1:1, 1:3, atoms%nat))
grVarphiGrtVeff0SphVarphiMeshNat = 0.
!!!!!!!!!!!!!!!!!!!!!!!
    iatom = 0
    do itype = 1, atoms%ntype

      ! These arrays are set to zero here, because
      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      delrVarphi1(:, :, :) = 0.
      delrVarphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(itype)
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
          end do ! imesh
          ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
          call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
          call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
        end do ! iradf
      end do ! oqn_l

      ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching coefficients
      ! have no spatial dependence) and determing its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
        do idirC = 1, 3
          call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, r2GrVeff0MTSh(:, :, iatom, idirC), &
                                                                                          & r2grVeff0SphVarphi(:, :, :, :, idirC) )
        end do

!        do idirR = 1, 3
!          do ilh = 1, lmpMax
!            do lm = 0, (atoms%lmax(itype) + 3)**2 -1
!              do imesh = 1, atoms%jri(1)
!                write(3009, '(5i8,2f30.8)') 1, idirR, lm, ilh, imesh, r2grVeff0SphVarphi(1, imesh, lm, ilh, idirR)
!                write(3009, '(5i8,2f30.8)') 2, idirR, lm, ilh, imesh, r2grVeff0SphVarphi(2, imesh, lm, ilh, idirR)
!              end do ! imesh
!            end do ! lm
!          end do ! ilh
!        end do ! idir
!        write(*, *) 'written1'
!        return
!
        ! Calculate the integral
        do mqn_m2PrC = -1, 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    do mqn_m2PrR = -1, 1
!      lmp = 0
!      lm = -1
!      do oqn_l = 0, atoms%lmax(itype)
!        do mqn_m = -oqn_l, oqn_l
!          lm = lm + 1
!          do iradf = 1, nRadFun(oqn_l, itype)
!            lmp = lmp + 1
!            lmp1Pr = 0
!            lm1Pr = -1
!            do oqn_l1Pr = 0, atoms%lmax(itype)
!              do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
!                mqn_m3Pr = grVarphiChMout(mqn_m1Pr, mqn_m2PrR)
!                lm1Pr = lm1Pr  + 1
!                do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
!                  lmp1Pr = lmp1Pr + 1
!                  ! overlap of grVarphi grVarphi
!                  do ichanPr = 1, 2
!                    oqn_l3Pr = grVarphiChLout(ichanPr, oqn_l1Pr)
!                    ! We have to catch channels where l < 0 and where m"' does not contribute, lmax + 1 is allowed due to the gradient
!                    if ( ( abs(mqn_m3Pr) > oqn_l3Pr )  .or. ( oqn_l3Pr < 0 )  )  cycle
!                    lm3Pr = oqn_l3Pr * (oqn_l3Pr + 1) + mqn_m3Pr
!
!                   ! intgrdR(:) = 0.
!                   ! intgrdI(:) = 0.
!                    do imesh = 1, atoms%jri(itype)
!      !                intgrdR(imesh) = real ( &
!      !                  &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC + 2) &
!      !                  &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC + 2) ) )
!      !                intgrdI(imesh) = aimag( &
!      !                  !todo why without iatom
!      !                  &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC + 2) &
!      !                  &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC + 2) ) )
!                      grVarphiGrtVeff0SphVarphiMeshNat(imesh, lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) =  grVarphiGrtVeff0SphVarphiMeshNat(imesh, lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) + ( &
!                        &   (grVarPhiCh1(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(1, imesh, lm3Pr, lmp, mqn_m2PrC + 2) &
!                        &  + grVarPhiCh2(imesh, ichanPr, lmp1Pr, mqn_m2PrR) * r2grVeff0SphVarphi(2, imesh, lm3Pr, lmp, mqn_m2PrC + 2) ) )
!                    end do ! imesh
!                   ! call intgr3(intgrdR(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralR)
!                   ! call intgr3(intgrdI(:), atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralI)
!                   ! ! < grVarphi| grTveff0Sph |varphi >
!                   ! ! idir = mqn_m2PrC + 2 due to performance reasons
!                   ! grVarphiGrtVeff0SphVarphi(lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) = &
!                   !          & grVarphiGrtVeff0SphVarphi(lmp1Pr, lmp, mqn_m2PrR, mqn_m2PrC + 2, iatom) + cmplx(integralR, integralI)
!                  end do ! ichanPr
!                end do ! iradf1Pr
!              end do ! mqn_m1Pr
!            end do ! oqn_l1Pr
!          end do ! iradf
!        end do ! mqn_m
!      end do ! oqn_l
!    end do ! mqn_m2PrR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          call CalcGrVarphiHepsGrtVarphiElem( atoms, itype, iatom, mqn_m2PrC, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, r2, &
            & grVarPhiCh1, grVarPhiCh2, El, lmpT, vEff0NsphGrVarphi, r2grVeff0SphVarphi, grVarphiGrtVarphi, grVarphiHpreGrtVarphi, &
            & grVarphiGrtVeff0SphVarphi)
        end do ! mqn_m2PrC
      end do ! ieqat
    end do ! itype

!    do itype = 1, atoms%ntype
!      do mqn_m2PrR  = -1, 1
!        do mqn_m2PrC = -1, 1
!          do ii = 1, lmpMax
!            do jj = 1, lmpMax
!              write(2292, '(4i8,f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiGrtVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!              write(2293, '(4i8,2f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiHpreGrtVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!              write(2294, '(4i8,2f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiGrtVeff0SphVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!            end do ! jj
!          end do ! ii
!        end do ! mqn_m2PrC
!      end do ! mqn_m2PrR
!    end do ! itype
!
!
    ispin = 1
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( ngoprI(atoms%nat) )
    a(:, :, :)               = cmplx(0., 0.)
    b(:, :, :)               = cmplx(0., 0.)
    bascof_lo(:, :, :, :, :) = cmplx(0., 0.)
    ! We do not want the local coordinate systems to be rotated for non-representative atoms constructing the matching coefficients.
    ngoprI(:) = 1
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        grPsiGrtVeff0PsiMt(:, :) = cmplx(0., 0.)
!!        altIntgrlIR(:, :) = cmplx(0., 0.)
        do ikpt = 1, kpts%nkpt
!!write(*, *) 'ikpt = ', ikpt
          nmat = nv(1, ikpt) + atoms%nlotot
          a(:, :, :) = cmplx(0.0, 0.0)
          b(:, :, :) = cmplx(0.0, 0.0)
          bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
          call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
            & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
            & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)), &
            & gbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), gbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
            & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
            & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

          grPsiGrtVeff0PsiNat(:, :) = cmplx(0., 0.)
          !!!!!!
!          grPsiGrtVeff0PsiNatMesh(:, :, :) = cmplx(0., 0.)
          !!!!!!!
          do iband = 1, nobd(ikpt, 1)
!write(*, *) 'iband = ', iband
            ab0cofScl(:) = cmplx(0., 0.)
            lmp = 0
            lm  = 0
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, itype)
                ! p = 1
                ab0cofScl(lmp + 1) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iatom) )
                ! p = 2
                ab0cofScl(lmp + 2) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iatom) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + &
                    & iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + &
                    & iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                end do ! iradf


                ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
                ! least the p as lm1Pr = lm. But for sake of performance we place it here.
                ! sake of performance
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l
            do idir = 1, 3
              do mqn_m2PrR = -1, 1
                grVarphiGrtVeff0PsiNat(1:lmpT(itype)) = cmplx(0., 0.)
                grVarphiGrtVeff0PsiNat(1:lmpT(itype)) = matmul(grVarphiGrtVeff0SphVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, idir, iatom), ab0cofScl(1:lmpT(itype)))
                grPsiGrtVeff0PsiNat(mqn_m2PrR, idir) = grPsiGrtVeff0PsiNat(mqn_m2PrR, idir) + 2 * results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), grVarphiGrtVeff0PsiNat(1:lmpT(itype)))


                !!!!!!!!!!!!!!!!!!!!!!!!

!                do imesh = 1, atoms%jri(itype)
!                  grVarphiGrtVeff0PsiNatMesh(1:lmpT(itype)) = cmplx(0., 0.)
!                  grVarphiGrtVeff0PsiNatMesh(1:lmpT(itype)) = matmul(grVarphiGrtVeff0SphVarphiMeshNat(imesh, 1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, idir, iatom), ab0cofScl(1:lmpT(itype)))
!                  grPsiGrtVeff0PsiNatMesh(imesh, mqn_m2PrR, idir) = grPsiGrtVeff0PsiNatMesh(imesh, mqn_m2PrR, idir) + 2 * results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), grVarphiGrtVeff0PsiNatMesh(1:lmpT(itype)))
!                end do ! imesh

                !!!!!!!!
              end do ! mqn_m2PrR
            end do ! idir
          end do ! iband
          grPsiGrtVeff0PsiMt(1:3, 1:3) = grPsiGrtVeff0PsiMt(1:3, 1:3) + matmul(conjg(Tmatrix(1:3, 1:3)), grPsiGrtVeff0PsiNat(-1:1, 1:3))
          !!!!!!!!
!          do imesh = 1, atoms%jri(itype)
!            grPsiGrtVeff0PsiMtMesh(imesh, 1:3, 1:3) = grPsiGrtVeff0PsiMtMesh(imesh, 1:3, 1:3) + matmul(conjg(Tmatrix(1:3, 1:3)), grPsiGrtVeff0PsiNatMesh(imesh, -1:1, 1:3))
!          end do ! imesh
          !!!!!!!

          ikpGz0 = cmplx(0.0, 0.0)
          !kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
          !gExt(:) = 0.
          !do idirR = 1, 3
          !  do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
          !    gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, ilst(iBas, ikpt, 1)))
          !    do iband = 1, nobd(ikpt, 1)
          !      ikpGz0(iBas, iband, idirR, ikpt) = iu * ( kExt(idirR) + gExt(idirR) ) * z(iBas, iband, ikpt, 1)
          !    end do ! iband
          !  end do ! iBas
          !end do ! idirR

          !nmat = nv(1, ikpt) + atoms%nlotot
          !do idirC = 1, 3
          !  do idirR = 1, 3
          !    altIntgrlIRband(:, :) = cmplx(0., 0.)
          !    call calcMEPotIR( input, stars, dimens, GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)),                   &
          !      & GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, w_vEff1IRContainer(:, idirC, 1),                  &
          !      & ikpGz0(:, :, idirR, ikpt), z(:, :, ikpt, 1), gdp, nmat, nobd(ikpt, 1), nobd(ikpt, 1), ispin, ikpt, iqpt, ikpt, ngdp, &
          !      & altIntgrlIRband(:, :), kpq2kPrVec )
          !    do iband = 1, nobd(ikpt, 1)
          !      altIntgrlIR(idirR, idirC) = altIntgrlIR(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * altIntgrlIRband(iband, iband)
          !    end do ! iband
          !  end do ! idirR
          !end do ! idirC
!
!
!          !do idirC = 1, 3
!          !  do idirR = 1, 3
!          !    do iband = 1, nobd(ikpt, 1)
!          !      altIntgrlIR = altIntgrlIR + altIntgrlIRband(iband, iband, idirR, idirC)
!          !    end do ! iband
!          !    altIntgrl(idirR, idirC) = altIntgrl(idirR, idirC) + altIntgrlIR(idirR, idirC) + altIntgrlMT(idirR, idirC)
!          !  end do ! idirR
!          !end do ! idirC
!

       end do !ikptR
        !altIntgrl(1:3, 1:3) = 2 * (altIntgrlIR(1:3, 1:3) + grPsiGrtVeff0PsiMt(1:3, 1:3))
        altIntgrl(1:3, 1:3) = 2 * (grPsiGrtVeff0PsiMt(1:3, 1:3))

!        allocate(integrandInt(atoms%jmtd, 3, 3))
!        integrandInt(:, :, :) = cmplx(0., 0.)
!        rewind(2284)
!        read(2284) integrandInt
!        do idirC = 1, 3
!          do idirR = 1, 3
!            do imesh = 1, atoms%jri(itype)
!!              intgrdR(imesh) = real ( integrandInt(imesh, idirR, idirC) - 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC) )
!!              intgrdI(imesh) = aimag( integrandInt(imesh, idirR, idirC) - 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC) )
!              intgrdR(imesh) = real ( 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC) )
!              intgrdI(imesh) = aimag( 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC) )
!            end do ! imesh
!            call intgr2(intgrdR, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), integralRMesh)
!            call intgr2(intgrdI, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), integralIMesh)
!            do imesh = 1, atoms%jri(itype)
!              write(2285, '(3i8,2f22.8)') idirC, idirR, imesh, integralRMesh(imesh), integralIMesh(imesh)
!            end do ! imesh
!          end do ! idirR
!        end do ! idirC

!        do idirC = 1, 3
!          do idirR = 1, 3
!            do imesh = 1, atoms%jri(itype)
!              write(2283, '(3i8,2f22.8)') idirC, idirR, imesh, 2 * grPsiGrtVeff0PsiMtMesh(imesh, idirR, idirC)
!            end do ! imesh
!          end do ! idirR
!        end do ! idirC


        write(*, '(a)') 'Matrix of brakets'
        write(*, '(3(2(es16.8,1x),3x))') altIntgrl(1, :)
        write(*, '(3(2(es16.8,1x),3x))') altIntgrl(2, :)
        write(*, '(3(2(es16.8,1x),3x))') altIntgrl(3, :)
      end do ! ieqat
    end do ! itype
    deallocate(a, b, bascof_lo)

  end subroutine TestDynMatPulInt

  ! Part of TesttensorGradMatElems
  subroutine CalcGrVeffGrtRhoInt( atoms, cell, lathar, dimens, stars, Veff0, input, ngdp, memd_atom, clnu_atom, nmem_atom, mlh_atom,      &
      & rho0IRst, gdp, rho0MT, convIntgrl )

    use m_types
    use m_jpConstants, only : iu
    use m_jpSetupDynMat, only : EvalIntRho1Veff1
    use mod_juPhonUtils, only : convertStar2G
    use m_jpPotDensHelper, only : calcGrR2FinLH, CalcIRdVxcKern, CalcMTdVxcKern, warpIRPot
    use m_jpGrVeff0, only : GenGrVeff0

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in) :: atoms
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_dimension),              intent(in) :: dimens
    type(t_stars),                  intent(in) :: stars
    type(t_potential),              intent(in) :: Veff0
    type(t_input),                  intent(in) :: input

    ! Scalar parameter
    integer,                        intent(in) :: ngdp
    integer,                        intent(in) :: memd_atom

    ! Array parameters
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    complex,                        intent(in) :: rho0IRst(:, :)
    integer,                        intent(in) :: gdp(:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)
    complex,                        intent(out):: convIntgrl(:, :)

    ! Scalar variables
    integer                                    :: iG
    integer                                    :: idirC
    integer                                    :: iatom
    integer                                    :: itype
    integer                                    :: ieqat
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: ilh
    integer                                    :: imesh
    logical                                    :: numericalGradient
    logical                                    :: harSw =.true.
    logical                                    :: extSw =.true.
    logical                                    :: xcSw =.true.
    logical                                    :: testGoldstein
    logical                                    :: fullGrVeff0

    ! Array variables
    complex,           allocatable             :: grRho0IR(:, :)
    real,              allocatable             :: r2Veff0MT(:, :, :)
    complex,           allocatable             :: grVeff0IR(:, :)
    complex,           allocatable             :: vEff0IRpw(:)
    complex,           allocatable             :: rho0IRpw(:)
    real,              allocatable             :: r2Rho0MTlh(:, :, :, :)
    real,              allocatable             :: gaussWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable             :: ylm2(:, :)
    complex,           allocatable             :: ylm1(:, : )
    real,              allocatable             :: dKernMTGPts(:, :, :)
    complex,           allocatable             :: grVxcIRKern(:)
    complex,           allocatable             :: rho1IRContainer(:, :, :)
    complex,           allocatable             :: rho1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: w_grVeff0IR(:, :)
    complex,           allocatable             :: w_vEff1IRContainer(:, :, :)
    complex,           allocatable             :: vEff1MTContainer(:, :, :, :, :)
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    complex,           allocatable             :: r2GrVeff0MT(:, :, :, :)
    complex,           allocatable             :: r2GrRho0MTsh(:, :, :, :)
    complex,           allocatable             :: grRho0MT(:, :, :, :)
    real,              allocatable             :: r2Rho0MT(:, :, :)
    complex,           allocatable             :: r2GrRho0MT(:, :, :, :)
    real                                       :: Gext(3)

    convIntgrl(:, :) = cmplx(0., 0.)


    !todo this code is copied and can be outsourced to a routine
    ! One can choose whether to calculate the numerical gradient of the unperturbed potential or to use the Weinert method for
    ! calculating it. The Weinert method has the advantage that continuity is enforced by construction at the MT boundary.
    ! It has been tested that both methods deliver the same numbers except for the just named behavior.
    numericalGradient = .false.
    if (numericalGradient) then

      ! Calculate the numerical gradient of the unperturbed effective potential.
      ! ------------------------------------------------------------------------

      ! Calculate the IR gradient.
      !   Convert the unperturbed effective potential from star expansion coefficients to plane-wave expansion coefficients
      allocate( vEff0IRpw(ngdp) )
      vEff0IRpw(:) = cmplx(0., 0.)
      ! SPIN MISSING
      call ConvertStar2G( Veff0%vpw_uw(:, 1), vEff0IRpw(:), stars, ngdp, gdp )

      ! Perform analytical gradient of effective potential in the interstitial region
      allocate( grVeff0IR(ngdp, 3) )
      grVeff0IR(:, :) = cmplx(0., 0.)
      do iG = 1, ngdp
        Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
        grVeff0IR(iG, 1:3) = iu * Gext(1:3) * vEff0IRpw(iG)
      end do ! iG

      ! Calculate the MT gradient.
      !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
      !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
      !   avoided.
      allocate( r2Veff0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
      allocate( grVeff0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
      r2Veff0MT(:, :, :) = 0.
      grVeff0MT(:, :, :, :) = cmplx(0., 0.)

      do itype = 1, atoms%ntype
        do ilh = 0, lathar%nlhd
          do imesh = 1, atoms%jri(itype)
            ! SPIN MISSING
            r2Veff0MT(imesh, ilh, itype) = Veff0%vr(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
          end do ! imesh
        end do ! ilh
      end do ! itype

      call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0MT, r2GrVeff0MT )

      do idirC = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              lm_pre = oqn_l * (oqn_l + 1) + 1
              do mqn_m = -oqn_l, oqn_l
                lm = lm_pre + mqn_m
                do imesh = 1, atoms%jri(itype)
                  grVeff0MT(imesh, lm, iatom, idirC) = r2GrVeff0MT(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! idirC

    else

      ! Use the Weinert method to calculate the gradient of the unperturbed effective potential
      ! ---------------------------------------------------------------------------------------

      ! Calculate the gradient of the unperturbed density.
      !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
      !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
      !   avoided.
      allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
      allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
      r2Rho0MT(:, :, :) = 0.
      grRho0MT(:, :, :, :) = cmplx(0., 0.)

      do itype = 1, atoms%ntype
        do ilh = 0, lathar%nlhd
          do imesh = 1, atoms%jri(itype)
            ! SPIN MISSING
            r2Rho0MT(imesh, ilh, itype) = rho0MT(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
          end do ! imesh
        end do ! ilh
      end do ! itype

      call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT, r2GrRho0MT )

      do idirC = 1, 3
        iatom = 0
        do itype = 1, atoms%ntype
          do ieqat = 1, atoms%neq(itype)
            iatom = iatom + 1
            do oqn_l = 0, atoms%lmax(itype)
              lm_pre = oqn_l * (oqn_l + 1) + 1
              do mqn_m = -oqn_l, oqn_l
                lm = lm_pre + mqn_m
                do imesh = 1, atoms%jri(itype)
                  grRho0MT(imesh, lm, iatom, idirC) = r2GrRho0MT(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
                end do ! imesh
              end do ! mqn_m
            end do ! oqn_l
          end do ! ieqat
        end do ! itype
      end do ! idirC

      ! Precalculate quantities for the kernel of the unperturbed effective potential's gradient.
      ! SPIN MISSING
      call CalcIRdVxcKern(stars, gdp, ngdp, rho0IRst(:, 1), grVxcIRKern)
      ! SPIN MISSING
      call CalcMTdVxcKern(atoms, dimens, lathar, rho0MT(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gaussWghts, ylm1, ylm2, dKernMTGPts)

      ! Generates gradient of the unperturbed effective potential by using the Weinert method.
      harSw = .true.
      extSw = .true.
      xcSw = .true.
      fullGrVeff0 = .true.
      testGoldstein = .false.
      ! SPIN MISSING
      call GenGrVeff0(atoms, cell, lathar, stars, dimens, memd_atom, ngdp, harSw, extSw, xcSw, gdp, &
        & rho0IRst( :, 1 ), rho0MT(:, :, :, 1), nmem_atom, mlh_atom, clnu_atom, grVeff0IR, grVeff0MT, grRho0MT, gaussWghts, ylm1, ylm2, &
        & dKernMTGPts, grVxcIRKern, testGoldstein, fullGrVeff0 )
      deallocate(grRho0MT)

    end if ! numericalGradient

    allocate(w_grVeff0IR(ngdp, 3))
    w_grVeff0IR = cmplx(0.0, 0.0)
    do idirC = 1, 3
      call warpIRPot(stars, ngdp, idirC, gdp, grVeff0IR, w_grVeff0IR(:, idirC))
    end do ! idirC

    ! Convert from star expansion coefficients to plane-wave expansion coefficients
    allocate( rho0IRpw(ngdp) )
    rho0IRpw(:) = cmplx(0., 0.)
    call convertStar2G( rho0IRst(:, 1), rho0IRpw(:), stars, ngdp, gdp )

    ! Perform analytical gradient of unperturbed density and effective potential in the interstitial region
    allocate( grRho0IR(ngdp, 3) )
    grRho0IR(:, :) = cmplx(0., 0.)
    do iG = 1, ngdp
      Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
      grRho0IR(iG, 1:3)  = iu * Gext(1:3) * rho0IRpw(iG)
    end do ! iG

    ! Gradient of density and effective potential in MT
    allocate( r2Rho0MTlh( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    r2Rho0MTlh(:, :, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)

    ! Read in valence density form FLEUR multiplied with r2
    read(1040) r2Rho0MtLh

    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MTlh(:, :, :, 1), r2GrRho0MTsh )

    do idirC = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                ! For reasons of performance, this calculation is done within the idirC loop
                grRho0MT(imesh, lm, iatom, idirC) = r2GrRho0MTsh(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idirC


    ! Fill container arrays with gradient coefficients so that they can be processed by EvelIntRho1Veff1
    allocate( rho1IRContainer(ngdp, 3, atoms%nat), w_vEff1IRContainer(ngdp, 3, atoms%nat), &
      & rho1MTContainer( atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3, atoms%nat ), &
      & vEff1MTContainer( atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3, atoms%nat ) )

    w_vEff1IRContainer(:, :, :) = cmplx(0., 0.)
    rho1IRContainer(:, :, :) = cmplx(0., 0.)
    rho1MTContainer(:, :, :, :, :) = cmplx(0., 0.)
    vEff1MTContainer(:, :, :, :, :) = cmplx(0., 0.)

    ! todo think about this for more than one atom

    do idirC = 1, 3
      do iG = 1, ngdp
        rho1IRContainer(iG, idirC, 1) = grRho0IR(iG, idirC)
        write(4000, '(2i8,2f15.8)') iG, idirC, rho1IRContainer(iG, idirC, 1)
        w_vEff1IRContainer(iG, idirC, 1) = w_grVeff0IR(iG, idirC)
        write(4001, '(2i8,2f15.8)') iG, idirC, w_vEff1IRContainer(iG, idirC, 1)
      end do ! iG
    end do ! idirC

    !todo does not work for more than one atom
    do idirC = 1, 3
      ! There should be an itype
      do oqn_l = 0, atoms%lmax(1)
        do mqn_m = -oqn_l, oqn_l
          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
          do imesh = 1, atoms%jri(1)
            rho1MTContainer(imesh, lm, 1, idirC, 1) = grRho0MT(imesh, lm, 1, idirC)
            write(4002, '(3i8,2f15.8)') imesh, lm, idirC, rho1MTContainer(imesh, lm, 1, idirC, 1)
            vEff1MTContainer(imesh, lm, 1, idirC, 1) = grVeff0MT(imesh, lm, idirC, 1)
            write(4003, '(3i8,2f15.8)') imesh, lm, idirC, vEff1MTContainer(imesh, lm, 1, idirC, 1)
          end do ! imesh
        end do ! mqn_m
      end do ! oqn_l
    end do ! idirC

    call EvalIntRho1Veff1(atoms, cell, rho1IRContainer, rho1MTContainer, w_vEff1IRContainer, vEff1MTContainer, ConvIntgrl)

  end subroutine CalcGrVeffGrtRhoInt

  ! One can apply the product rule and then has some terms with the tensor gradient which have to be fulfilled or be the same than
  ! surface integrals. This test does not work so well and is deprecated when deciding not to use tensor gradients.
  subroutine TestTensorGradProdKet( atoms, lathar, stars, dimens, cell, sym, kpts, usdus, results, input, Veff0, qpts, ngdp, memd_atom, rho0IRst, rho0MTlh, mlh_atom, nmem_atom,      &
      & clnu_atom, gdp, nRadFun, rbas1, rbas2, nv, mapGbas, gbas, kveclo, nobd, z0, iloTable, kpq2kPrVec, El, eig )

    use m_types, only : t_atoms, t_sphhar, t_stars, t_dimension, t_cell, t_sym, t_kpts, t_usdus, t_results, t_input, t_potential
    use m_jPConstants, only : iu, fpi, Tmatrix_transposed, Tmatrix
    use m_jpPotDensHelper, only : CalcGrR2FinLh, CalcIRdVxcKern, CalcMTdVxcKern, CalcGrR2FinSH, WarpIRPot
    use m_jpGrVeff0, only : GenGrVeff0
    use m_jpTestPotential, only : checkjuPhPots
    use m_jpSetupDynMatHelper, only : CalcFnsphVarphi, CalcFnsphGrVarphi, CalcHnGrV0Varphi
    use m_jpSetupDynMat, only : CalcScalBasfMatElems, CalcSelfAdjCorrection, CalcVecBasfMatElems, calcPsi1HepsPsi1IR, CalcGrVarphiHepsGrtVarphiElem, SetupDynMatHF
    use m_od_types, only : od_inp, od_sym
    use m_abcof3
    use m_jpSternhHF, only : CalcMEPotIR
    use mod_juPhonUtils, only : convertStar2G, Derivative, CalcChannelsGrFlpNat, CalcChannelsGrGrtFlpNat
    use m_jpSetupDynMatSF, only : CalcSFintIRPsiHepsGradPsi, CalcSFintMTPsiHepsGradPsi


    use mod_juPhonUtils, only : fopen, fclose
    use m_intgr, only : intgr3NoIntp, intgr3LinIntp

    implicit none

    ! Type parameters
    type(t_atoms),                  intent(in)  :: atoms
    type(t_sphhar),                 intent(in)  :: lathar
    type(t_stars),                  intent(in)  :: stars
    type(t_dimension),              intent(in)  :: dimens
    type(t_cell),                   intent(in)  :: cell
    type(t_sym),                    intent(in)  :: sym
    type(t_kpts),                   intent(in)  :: kpts
    type(t_kpts),                   intent(in)  :: qpts
    type(t_usdus),                  intent(in) :: usdus
    type(t_results),                intent(in) :: results
    type(t_input),                  intent(in) :: input
    type(t_potential),              intent(in) :: Veff0

    ! Scalar parameters
    integer,                        intent(in)  :: ngdp
    integer,                        intent(in)  :: memd_atom

    ! Array parameters
    complex,                        intent(in)  :: rho0IRst(:, :)
    real,                           intent(in)  :: rho0MTlh(:, 0:, :, :)
    integer,                        intent(in)  :: mlh_atom(:,0:,:)
    integer,                        intent(in)  :: nmem_atom(0:, :)
    complex,                        intent(in)  :: clnu_atom(:,0:,:)
    integer,                        intent(in)  :: gdp(:, :)
    integer,                        intent(in)  :: nRadFun(0:, :)
    real,                           intent(in)  :: rbas1(:, :, 0:, :)
    real,                           intent(in)  :: rbas2(:, :, 0:, :)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: gbas(:, :)
    integer,                        intent(in) :: mapGbas(:, :, :)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in) :: nobd(:, :)
    MCOMPLEX,                       intent(in) :: z0(:, :, :, :)
    integer,                        intent(in) :: iloTable(:, 0:, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    real,                           intent(in)  :: El(:, 0:, :, :)
    real,                           intent(in)    :: eig(:, :, :)

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods

    ! Scalar variables
    integer                                     :: ilh
    integer                                     :: itype
    integer                                     :: idirC
    integer                                     :: idirR
    integer                                    :: iband
    integer                                     :: iG
    integer                                     :: iatom
    integer                                     :: ieqat
    integer                                     :: lmp
    integer                                     :: oqn_l
    integer                                     :: lm_pre
    integer                                     :: lmpMax
    integer                                     :: nRadFunMax
    integer                                     :: nmat
    integer                                     :: iradf
    integer                                     :: mqn_m
    integer                                     :: lm
    integer                                     :: imesh
    integer                                     :: ikpt
    logical                                     :: harSw
    integer                                    :: pMaxLocal
    logical                                     :: extSw
    logical                                     :: xcSw
    logical                                     :: testGoldstein
    logical                                     :: fullGrVeff0
    integer                                     :: nobdMax
    integer                                     :: iqpt
    integer                                     :: ptsym
    integer                                     :: imem
    integer                                     :: mqn_m2PrC
    integer                                     :: iBas
    integer                                     :: mqn_m2PrR
    integer                                     :: ichan
    complex                                     :: psiHepsgrGrtPsiIR
    complex                                     :: grPsiHepsGrtPsiIR
    complex                                     :: traceLeft
    complex                                     :: singleTrace
    integer                                     :: ii

    ! Array variables
    real,              allocatable              :: r2Rho0MT(:, :, :)
    complex,           allocatable              :: grRho0MT(:, :, :, :) ! Dummy quantity at the moment
    complex,           allocatable              :: r2GrRho0MT(:, :, :, :)
    complex,           allocatable              :: grVxcIRKern(:)
    real,              allocatable              :: gaussWghts(:) ! gaussian weights belonging to gausPts
    complex,           allocatable              :: ylm2(:, :)
    complex,           allocatable              :: ylm1(:, : )
    real,              allocatable              :: dKernMTGPts(:, :, :)
    complex,           allocatable              :: grVxc0IR(:, :)
    complex,           allocatable              :: grVeff0IR(:, :)
    complex,           allocatable             :: ab0cofScl(:)
    complex,           allocatable              :: grVxc0MT(:, :, :, :)
    complex,           allocatable              :: grVeff0MT(:, :, :, :)
    complex,           allocatable              :: r2GrVxc0MT(:, :, :, :)
    complex,           allocatable              :: grGrtVxc0IR(:, :, :)
    complex,           allocatable              :: r2grGrtVxc0MTtemp(:, :, :, :)
    complex,           allocatable              :: grGrtVXc0MT(:, :, :, :, :)
    real,              allocatable              :: varphi1(:, :, :)
    real,              allocatable              :: varphi2(:, :, :)
    complex,           allocatable              :: grGrtVxc0MtSphVarphi(:, :, :, :)
    real,              allocatable              :: varphiVarphiDummy(:, :, :)
    complex,           allocatable              :: varphigradGradVxcvarphi(:, :, :, :, :)
    integer,           allocatable              :: lmpT(:)
    real,              allocatable              :: r2(:)
    integer,           allocatable              :: ngoprI(:)
    complex,           allocatable              :: a(:, :, :)
    complex,           allocatable              :: b(:, :, :)
    complex,           allocatable              :: bascof_lo(:, :, :, :, :)
    complex,           allocatable              :: varphigradGradVxcPsi(:)
    complex,           allocatable              :: gradz0(:, :, :)
    complex,           allocatable              :: psigradGradVxcPsi(:, :, :)
    complex,           allocatable              :: w_grGrtVxcIR(:, :, :)
    complex,           allocatable              :: psiGrGrtVxcPsiIR(:, :)
    complex,           allocatable              :: psiGrGrtVxcPsiIRSum(:, :, :)
    complex,           allocatable              :: trPsiGrGrtVharPsiIR(:, :)
    complex,           allocatable              :: trPsiGrGrtVharPsiIRSum(:)
    complex,           allocatable              :: trGrGrtVharVarphi(:, :, :, :)
    real,              allocatable              :: varphiVarphi(:, :, :)
    complex,           allocatable              :: trvarphiGrGrtVharVarphi(:, :, :)
    complex,           allocatable              :: trPsiGrGrtVharPsiMTSum(:)
    complex,           allocatable              :: trVarphiGrGrtVharPsi(:)
    complex,           allocatable              :: varphiPsi(:)
    complex,           allocatable              :: rhoValMt(:)
    complex,           allocatable              :: rho0IRpwContainer(:, :)
    complex,           allocatable              :: w_rho0IRpw(:)
    complex,           allocatable              :: rho0MTSpH(:, :, :)
    real,              allocatable              :: delrVarphi1(:, :, :)
    real,              allocatable              :: delrVarphi2(:, :, :)
    real,              allocatable              :: grVarphiCh1(:, :, :, :)
    real,              allocatable              :: grVarphiCh2(:, :, :, :)
    integer,           allocatable              :: grVarphiChLout(:, :)
    integer,           allocatable              :: grVarphiChMout(:, :)
    complex,           allocatable              :: grVeff0GrtVarphi(:, :, :, :)
    complex,           allocatable              :: varphiGrVeff0GrtVarphi(:, :, :, :, :)
    complex,           allocatable              :: psiGrVeffGrtPsiIR(:, :)
    complex,           allocatable              :: psiGrVeff0GrtPsiIRSum(:, :, :)
    complex,           allocatable              :: varphiGrVeff0GrtPsi(:)
    complex,           allocatable              :: psiGrVeff0GrtPsiMT(:, :, :)
    complex,           allocatable              :: w_grVeff0IR(:, :)
    complex,           allocatable              :: vEff0MtSpH(:, :)
    real,              allocatable              :: delrGrVarphiCh1(:, :, :, :)
    real,              allocatable              :: delrGrVarphiCh2(:, :, :, :)
    integer,           allocatable              :: grGrtVarphiChLout(:, :)
    integer,           allocatable              :: grGrtVarphiChMout(:, :, :)
    real,              allocatable              :: grGrtVarphiCh1(:, :, :, :, :)
    real,              allocatable              :: grGrtVarphiCh2(:, :, :, :, :)
    real,              allocatable              :: sAdjCorrVarphiHkinGrGrtVarphi(:, :, :, :)
    real,              allocatable              :: delrGrGrtVarphiCh1(:, :, :, :)
    real,              allocatable              :: delrGrGrtVarphiCh2(:, :, :, :)
    complex,           allocatable              :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable              :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable              :: hVarphi(:, :, :, :)
    complex,           allocatable              :: grGrtVarphiHvarphi(:, :, :, :, :)
    real,              allocatable              :: grGrtVarphiVarphi(:, :, :, :, :)
    complex,           allocatable              :: grGrtPsiHepsPsi(:, :, :)
    complex,           allocatable              :: sAdjCorrPsiHkinGrGrtPsi(:, :, :)
    complex,           allocatable              :: sAdjCorrPsiGrGrtVarphiNat(:)
    complex,           allocatable              :: grGrtVarphiPsiNatNat(:)
    complex,           allocatable              :: grGrtVarphiHPsiNatNat(:)
    complex,           allocatable              :: psiHepsgrGrtPsiIRSum(:, :, :)
    complex,           allocatable              :: grGrTz0(:, :, :)
    complex,           allocatable              :: grPsiGrtVeffPsiIR(:, :)
    complex,           allocatable              :: grPsiGrtVeff0PsiIRSum(:, :, :)
    complex,           allocatable              :: grtVeff0MtSphVarphi(:, :, :, :)
    complex,           allocatable              :: grVarphiGrtVeff0Varphi(:, :, :, :, :)
    complex,           allocatable              :: grVarphiGrtVeff0PsiNat(:)
    complex,           allocatable              :: grPsiGrtVeff0Psi(:, :, :)
    complex,           allocatable              :: grVarphiGrtVeff0SphVarphi(:, :, :, :, :)
    real,              allocatable              :: grVarphiGrtVarphi(:, :, :, :, :)
    complex,           allocatable              :: grVarphiHpreGrtVarphi(:, :, :, :, :)
    complex,           allocatable              :: grPsiHepsGrtPsi(:, :, :)
    complex,           allocatable              :: grVarphiGrtVeff0SphPsiNat(:)
    complex,           allocatable              :: grVarphiHepsGrtPsiNatNat(:)
    complex,           allocatable              :: grPsiHepsGrtPsiIRSum(:, :, :)
    complex,           allocatable              :: surfInts(:, :, :)
    real,              allocatable              :: r2Rho0MTlhVal(:, :, :, :)
    complex,        allocatable       :: r2grVeff0SphSh( :, :, :, : )
    complex                                     :: grPsiHepsGrtPsiNatNat(-1:1, -1:1)
    complex                                     :: grPsiHepsGrtPsiNat(-1:1, 3)
    complex                                     :: grPsiGrtVeff0SphPsiNat(-1:1, 3)
    complex                                     :: grPsiGrtVeff0PsiNat(-1:1, 1:3)
    complex                                     :: grGrtPsiHepsPsiNat(3, -1:1)
    complex                                     :: grGrtPsiHepsPsiNatNat(-1:1, -1:1)
    complex                                     :: sAdjCorrPsiHkinGrGrtPsiNat(3, -1:1)
    complex                                     :: sAdjCorrPsiHkinGrGrtPsiNatNat(-1:1, -1:1)
    complex                                     :: psiGrVeff0GrtPsiMTNat(3, -1:1)
    complex                                     :: surfInt(3, 3)
    real                                        :: Gext(3)
    real                                        :: kExt(3)


    real,           allocatable       :: r2Veff0SphLh(:, :, :)
    real,           allocatable       :: vxc0mt(:, :, :, :)
    real, allocatable :: r2Unit(:)
    complex, allocatable :: grVeff0MTtest(:, :, :)
    integer :: oqn_l1Pr
    integer :: mqn_m1Pr
    integer :: lmp1Pr
    integer :: iradf1Pr

    complex, allocatable :: r2Rho0MTsphVal(:, :)
    real, allocatable :: fReal(:)
    real, allocatable :: fImag(:)
    complex :: integral
    real :: intgrReal
    real :: intgrImag

    complex,          allocatable               :: w_vExt2IRDummy(:, :, :, :)
    complex,          allocatable               :: vExt2MTDummy(:, :, :, :)
    complex,          allocatable               :: rho1IRDummy(:, :, :)
    complex,          allocatable               :: rho1MTDummy(:, :, :, :, :)
    complex,          allocatable               :: w_vExt1IRDummy(:, :, :)
    complex,          allocatable               :: vExt1MTDummy(:, :, :, :, :)
    complex :: E2ndOrdII(3, 3)
    complex :: dynMatHF(3, 3)
    complex, allocatable :: r2GrRho0MTVal(:, :, :, :)
    complex, allocatable :: grRho0MTVal(:, :, :, :)
    real, allocatable :: rho0MTlhVal(:, :, :, :)

    complex, allocatable :: grGrtRho0MTval(:, :, :, :, :)
    complex, allocatable :: r2grGrtRho0MTvalTemp(:, :, :, :)
    complex, allocatable :: grGrtRho0IR(:, :, :)


    logical :: hfMTOn = .false.
    logical :: hfIRON = .false.
    logical :: tensGradOvlOn = .true.

    complex :: integralMat(3, 3)


    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )
    nobdMax    = maxval( nobd(:, :) )

    allocate( r2Rho0MTlhVal( atoms%jmtd, 0:lathar%nlhd, atoms%ntype, input%jspins) )
    r2Rho0MtlhVal(:, :, :, :) = cmplx(0., 0.)
    read(1040) r2Rho0MtLhVal
    ! Calculate the gradient of the xc-potential
    ! ---------------------------------------------------------------------------------------

    ! Calculate the gradient of the unperturbed density.
    !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
    !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
    !   avoided.
    allocate( r2Rho0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
    allocate( grRho0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    r2Rho0MT(:, :, :) = 0.
    grRho0MT(:, :, :, :) = cmplx(0., 0.)
    allocate( grRho0MTVal( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    grRho0MTVal(:, :, :, :) = cmplx(0., 0.)

    allocate( rho0MtlhVal(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1))
    rho0MtlhVal(:, :, :, :) = cmplx(0., 0.)
    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          ! SPIN MISSING
          r2Rho0MT(imesh, ilh, itype) = rho0MTlh(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
          rho0MtlhVal(imesh, ilh, itype, 1) = r2Rho0MtlhVal(imesh, ilh, itype, 1) / atoms%rmsh(imesh, itype) / atoms%rmsh(imesh, itype)
        end do ! imesh
      end do ! ilh
    end do ! itype

    call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MT, r2GrRho0MT )
    call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Rho0MTlhVal(:, :, :, 1), r2GrRho0MTVal )

    do idirC = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                grRho0MT(imesh, lm, iatom, idirC) = r2GrRho0MT(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
                grRho0MTVal(imesh, lm, iatom, idirC) = r2GrRho0MTVal(imesh, lm, iatom, idirC) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idirC

    ! Precalculate quantities for the kernel of the unperturbed effective potential's gradient.
    ! SPIN MISSING
    call CalcIRdVxcKern(stars, gdp, ngdp, rho0IRst(:, 1), grVxcIRKern)
    ! SPIN MISSING
    call CalcMTdVxcKern(atoms, dimens, lathar, rho0MTlh(:, :, :, 1), nmem_atom, clnu_atom, mlh_atom, gaussWghts, ylm1, ylm2, dKernMTGPts)

    ! Generates gradient of the unperturbed effective potential by using the Weinert method.
    harSw = .true.
    extSw = .true.
    xcSw = .true.
    fullGrVeff0 = .true.
    testGoldstein = .false.
    ! SPIN MISSING
    call GenGrVeff0(atoms, cell, lathar, stars, dimens, memd_atom, ngdp, harSw, extSw, xcSw, gdp, &
      & rho0IRst( :, 1 ), rho0MTlh(:, :, :, 1), nmem_atom, mlh_atom, clnu_atom, grVeff0IR, grVeff0MT, grRho0MT, gaussWghts, ylm1, ylm2, &
      & dKernMTGPts, grVxcIRKern, testGoldstein, fullGrVeff0 )


    allocate( w_vExt2IRDummy(ngdp, 3, 3, 1), vExt2MTDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, 3, 3), rho1IRDummy(ngdp, 3, 1), &
      & rho1MTDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, 1, 3, 1), w_vExt1IRDummy(ngdp, 3, 1), &
      & vExt1MTDummy(atoms%jmtd, (atoms%lmaxd + 1)**2, 1, 3, 1) )

    w_vExt2IRDummy(:, :, :, :) = cmplx(0., 0.)
    vExt2MTDummy(:, :, :, :) = cmplx(0., 0.)
    rho1IRDummy(:, :, :) = cmplx(0., 0.)
    rho1MTDummy(:, :, :, :, :) = cmplx(0., 0.)
    w_vExt1IRDummy(:, :, :) = cmplx(0., 0.)
    vExt1MTDummy(:, :, :, :, :) = cmplx(0., 0.)
    E2ndOrdII(:, :) = cmplx(0., 0.)
    dynMatHF(:, :) = cmplx(0., 0.)
    if (.false.) then
      do idirC = 1, 3
        do oqn_l = 0, atoms%lmax(1)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atoms%jri(1)
              rho1MTDummy(imesh, lm, 1, idirC, 1) = -grRho0MTVal(imesh, lm, 1, idirC)
              vExt1MTDummy(imesh, lm, 1, idirC, 1) = -grVeff0MT(imesh, lm, idirC, 1)
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end do !
    end if
!    write(4000, *) rho1MTDummy
!    write(4001, *) vExt1MTDummy
if (.false.) then
    call SetupDynMatHF(atoms, cell, lathar, ngdp, ngdp, mlh_atom, nmem_atom, clnu_atom, rho0IRst(:, 1), w_vExt2IRDummy, vExt2MTDummy, rho1IRDummy, rho1MTDummy, w_vExt1IRDummy,     &
      & vExt1MTDummy, rho0MtlhVal, E2ndOrdII, dynMatHF)
end if

    if (.false.) then
    ! Generates gradient of the unperturbed effective potential by using the Weinert method.
    harSw = .false.
    extSw = .false.
    xcSw = .true.
    fullGrVeff0 = .true.
    testGoldstein = .false.
    ! SPIN MISSING
    call GenGrVeff0(atoms, cell, lathar, stars, dimens, memd_atom, ngdp, harSw, extSw, xcSw, gdp, &
      & rho0IRst( :, 1 ), rho0MTlh(:, :, :, 1), nmem_atom, mlh_atom, clnu_atom, grVxc0IR, grVxc0MT, grRho0MT, gaussWghts, ylm1, ylm2, &
      & dKernMTGPts, grVxcIRKern, testGoldstein, fullGrVeff0 )
    deallocate(grRho0MT)

    !grVeff0IR(1:ngdp, 1:3) = grVeff0IR(1:ngdp, 1:3)
    grVeff0IR(1:ngdp, 1:3) = grVeff0IR(1:ngdp, 1:3) + grVxc0IR(1:ngdp, 1:3)
    grVeff0MT(1:atoms%jmtd, 1:(atoms%lmaxd + 2)**2, 1:3, 1:atoms%nat) = grVeff0MT(1:atoms%jmtd, 1:(atoms%lmaxd + 2)**2, 1:3, 1:atoms%nat) + grVxc0MT(1:atoms%jmtd, 1:(atoms%lmaxd + 2)**2, 1:3, 1:atoms%nat)
!    grVeff0MT(1:atoms%jmtd, 1:(atoms%lmaxd + 2)**2, 1:3, 1:atoms%nat) = grVxc0MT(1:atoms%jmtd, 1:(atoms%lmaxd + 2)**2, 1:3, 1:atoms%nat)
    end if


    ! Calculate the tensor gradient of the xc-potential
    allocate( r2grVxc0MT( atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%nat, 3 ) )
    r2grVxc0MT = cmplx(0., 0.)
    do idirC = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, 2!atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                r2grVxc0MT(imesh, lm, iatom, idirC) = grVeff0MT(imesh, lm, idirC, iatom) * atoms%rmsh(imesh, itype)**2
!                r2grVxc0MT(imesh, lm, iatom, idirC) = grVxc0MT(imesh, lm, idirC, iatom) * atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_mj
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idirC

    allocate( grGrtVxc0IR(ngdp, 3,  3), r2grGrtVxc0MTtemp(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%ntype, 3), &
                                                                & grGrtVXc0MT(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%ntype, 3) )
    allocate( grGrtRho0IR(ngdp, 3, 3), r2grGrtRho0MTvalTemp(atoms%jmtd, (atoms%lmaxd + 2)**2, atoms%ntype, 3), grGrtRho0MTval(atoms%jmtd, (atoms%lmaxd + 2)**2, 3, atoms%nat, 3) )
    grGrtVxc0IR(:, :, :) = cmplx(0., 0.)
    grGrtVxc0MT(:, :, :, :, :) = cmplx(0., 0.)
    grGrtRho0MTval(:, :, :, :, :) = cmplx(0., 0.)
    grGrtRho0IR(:, :, :) = cmplx(0., 0.)
    allocate( rho0IRpwContainer(ngdp, 1) )
    rho0IRpwContainer(:, :) = cmplx(0., 0.)
    call convertStar2G( rho0IRst(:, 1), rho0IRpwContainer(:, 1), stars, ngdp, gdp )

    do idirC = 1, 3
      do idirR = 1, 3
        do iG = 1, ngdp
          Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
!          grGrtVxc0IR(iG, idirR, idirC) = iu * Gext(idirR) * grVxc0IR(iG, idirC)
          grGrtVxc0IR(iG, idirR, idirC) = iu * Gext(idirR) * grVeff0IR(iG, idirC)
          grGrtRho0IR(iG, idirR, idirC) = - Gext(idirR) * Gext(idirC) * rho0IRpwContainer(iG, 1)
        end do ! iG
      end do ! idirR
      if (hfMTOn .or. tensGradOvlOn) then
        !todo be careful it is only the type here, we have to change that!
        r2grGrtVxc0MTtemp = cmplx(0., 0.)
        r2grGrtRho0MTvaltemp = cmplx(0., 0.)
        call CalcGrR2FinSH( atoms, r2grVxc0MT(:, :, :, idirC), r2grGrtVxc0MTtemp)
        call CalcGrR2FinSH( atoms, r2GrRho0MTVal(:, :, :, idirC), r2grGrtRho0MTvalTemp )
        do idirR = 1, 3
          iatom = 0
          do itype = 1, atoms%ntype
            do ieqat = 1, atoms%neq(itype)
              iatom = iatom + 1
              do oqn_l = 0, atoms%lmax(itype)
                lm_pre = oqn_l * (oqn_l + 1) + 1
                do mqn_m = -oqn_l, oqn_l
                  lm = lm_pre + mqn_m
                  do imesh = 1, atoms%jri(itype)
                    grGrtVxc0MT(imesh, lm, idirR, itype, idirC) = r2grGrtVxc0MTtemp(imesh, lm, itype, idirR) / atoms%rmsh(imesh, itype)**2
                    grGrtRho0MTval(imesh, lm, idirR, iatom, idirC) = r2grGrtRho0MTvalTemp(imesh, lm, itype, idirR) / atoms%rmsh(imesh, itype)**2
                  end do ! imesh
                end do ! mqn_m
              end do ! oqn_l
            end do ! ieqat
          end do ! itype
        end do ! idirR
        if (.false.) then
          write(100, '(a)') 'foo', idirC
          call checkjuPhPots(30, atoms, ngdp, gdp, cell, lathar, sym, grGrtVxc0IR(:, :, idirC), grGrtVxc0MT(:, :, :, :, idirC), 100)
          call checkjuPhPots(30, atoms, ngdp, gdp, cell, lathar, sym, grGrtRho0IR(:, :, idirC), grGrtRho0MTval(:, :, :, :, idirC), 100)
        end if
      end if
    end do ! idirC


    allocate( w_grGrtVxcIR(ngdp, 3, 3))
    allocate( w_rho0IRpw(ngdp))
    allocate( w_grVeff0IR(ngdp, 3))
    w_grGrtVxcIR(:, :, :) = cmplx(0., 0.)
    w_rho0IRpw(:) = cmplx(0., 0.)


    do idirR = 1, 3
      call warpIRPot( stars, ngdp, idirR, gdp, grVeff0IR, w_grVeff0IR(:, idirR) )
    end do ! idirR

    call warpIRPot( stars, ngdp, 1, gdp, rho0IRpwContainer, w_rho0IRpw )
    do idirC = 1, 3
      do idirR = 1, 3
        call warpIRPot(stars, ngdp, idirR, gdp, grGrtVxc0IR(:, :, idirC), w_grGrtVxcIR(:, idirR, idirC))
      end do ! idirR
    end do !idirC

    if (.false.) then
      allocate( rho0MtSpH(atoms%jmtd, (atoms%lmaxd + 1)**2, atoms%nat) )
      rho0MtSpH(:, :, :) = cmplx(0.0, 0.0)
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          ptsym = atoms%ntypsy(iatom)
          do ilh = 0, lathar%nlh(ptsym)
            oqn_l = lathar%llh(ilh, ptsym)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)
              lm = lm_pre + mqn_m
              !todo one could only evaluate the rho0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
              ! maybe construct a pointer and run only over them to make it memory efficient.
              do imesh = 1, atoms%jri(itype)
                rho0MtSpH(imesh, lm, iatom) = rho0MtSpH(imesh, lm, iatom) + rho0MTLh(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
              end do ! imesh
            end do ! imem
          end do ! ilh
        end do ! ieqat
      end do ! itype
    end if

    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( r2(atoms%jmtd) )

    allocate( grGrtVxc0MtSphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )
    allocate( varphiVarphiDummy(lmpMax, lmpMax, atoms%ntype), varphigradGradVxcvarphi(lmpMax, lmpMax, atoms%nat, 3, 3))
    allocate( psigradGradVxcPsi(atoms%ntype, 3, 3))
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
            & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate(ab0cofScl(lmpMax))
    allocate( varphigradGradVxcPsi(lmpMax))
    allocate( ngoprI(atoms%nat) )
    allocate( psiGrGrtVxcPsiIR(nobdMax, nobdMax) )
    allocate( psiGrGrtVxcPsiIRSum(atoms%ntype, 3, 3) )
    allocate( trPsiGrGrtVharPsiIR(nobdMax, nobdMax) )
    allocate( trPsiGrGrtVharPsiIRSum(atoms%nat) )
    allocate( trGrGrtVharVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax) )
    allocate( grVeff0GrtVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax ) )
    allocate( grtVeff0MtSphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax ) )
    allocate( varphiVarphi(lmpMax, lmpMax, atoms%ntype), trvarphiGrGrtVharVarphi(lmpMax, lmpMax, atoms%nat) )
    allocate( trPsiGrGrtVharPsiMTSum(atoms%nat), trvarphiGrGrtVharPsi(lmpMax) )
    allocate( varphiPsi(lmpMax), rhoValMt(atoms%nat) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( varphiGrVeff0GrtVarphi(lmpMax, lmpMax, atoms%nat, 3, -1:1) )
    allocate( grVarphiGrtVeff0Varphi(lmpMax, lmpMax, -1:1, 3, atoms%nat) )
    allocate( psiGrVeffGrtPsiIR(nobdMax, nobdMax), psiGrVeff0GrtPsiIRSum(3, 3, atoms%nat) )
    allocate( varphiGrVeff0GrtPsi(lmpMax) )
    allocate( psiGrVeff0GrtPsiMT(3, 3, atoms%nat) )
    allocate( gradz0( dimens%nbasfcn, nobdMax, 3) )
    allocate( grGrtVarphiVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grGrtVarphiHvarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( grGrtVarphiChLout(4, 0:atoms%lmaxd), grGrtVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1, -1:1), &
            & grGrtVarphiCh1(atoms%jmtd, 4, lmpMax, -1:1, -1:1), grGrtVarphiCh2(atoms%jmtd, 4, lmpMax, -1:1, -1:1) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( sAdjCorrVarphiHkinGrGrtVarphi(lmpMax, lmpMax, -1:1, -1:1) )
    allocate( delrGrGrtVarphiCh1(atoms%jmtd, 4, lmpMax, -1:1), delrGrGrtVarphiCh2(atoms%jmtd, 4, lmpMax, -1:1) )
    allocate( delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( grGrtPsiHepsPsi(3, 3, atoms%nat) )
    allocate( sAdjCorrPsiHkinGrGrtPsi(3, 3, atoms%ntype) )
    allocate( sAdjCorrPsiGrGrtVarphiNat(lmpMax), grGrTz0(dimens%nbasfcn, 3, 3) )
    allocate( psiHepsgrGrtPsiIRSum(3, 3, atoms%nat))
    allocate( grGrtVarphiHPsiNatNat(lmpMax), grGrtVarphiPsiNatNat(lmpMax) )
    allocate( grPsiGrtVeff0PsiIRSum(3, 3, atoms%nat), grPsiGrtVeffPsiIR(nobdMax, nobdMax) )
    allocate( grVarphiGrtVeff0PsiNat(lmpMax) )
    allocate( grPsiGrtVeff0Psi(3, 3, atoms%nat) )
    allocate( grVarphiGrtVeff0SphVarphi(lmpMax, lmpMax, -1:1, 1:3, atoms%nat) )
    allocate( grVarphiGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grVarphiHpreGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( grPsiHepsGrtPsi(3, 3, atoms%nat))
    allocate( grVarphiGrtVeff0SphPsiNat(lmpMax), grVarphiHepsGrtPsiNatNat(lmpMax))
    allocate( grPsiHepsGrtPsiIRSum(3, 3, atoms%nat) )
    allocate( surfInts(3, 3, atoms%nat) )

    varphiVarphiDummy(:, :, :) = cmplx(0., 0.)
    varphiVarphi(:, :, :) = cmplx(0., 0.)
    trvarphiGrGrtVharVarphi(:, :, :) = cmplx(0., 0.)
    varphiGrVeff0GrtVarphi(:, :, :, :, :) = cmplx(0., 0.)
    psiGrVeffGrtPsiIR(:, :) = cmplx(0., 0.)
    psiGrVeff0GrtPsiIRSum(:, :, :) = cmplx(0., 0.)
    grGrtVarphiHvarphi(:, :, :, :, :)   = cmplx(0., 0.)
    grGrtVarphiVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grGrtPsiHepsPsi(:, :, :) = cmplx(0., 0.)
    psiHepsgrGrtPsiIRSum(:, :, :) = cmplx(0., 0.)
    sAdjCorrPsiHkinGrGrtPsi(:, :, :) = cmplx(0., 0.)
    grPsiGrtVeff0PsiIRSum(:, :, :) = cmplx(0., 0.)
    grVarphiGrtVeff0Varphi(:, :, :, :, :) = cmplx(0., 0.)
    grPsiGrtVeff0Psi(:, :, :) = cmplx(0., 0.)
    grVarphiGrtVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiGrtVarphi(:, :, :, :, :) = 0.
    grVarphiHpreGrtVarphi(:, :, :, :, :)   = cmplx(0., 0.)
    grPsiHepsGrtPsiIRSum(:, :, :) = cmplx(0., 0.)
    surfInts(:, :, :) = cmplx(0., 0.)

    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1

        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.
        delrVarphi1(:, :, :) = 0.
        delrVarphi2(:, :, :) = 0.
        do oqn_l = 0, atoms%lmax(itype)
          do iradf = 1, nRadFun(oqn_l, itype)
            do imesh = 1, atoms%jri(itype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype) / atoms%rmsh(imesh, itype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype) / atoms%rmsh(imesh, itype)
            end do ! imesh
            call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
            call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
          end do ! iradf
        end do ! oqn_l

       ! Precalculate radial Jacobi determinant for later integrals
       r2(:) = 0.
       do imesh = 1, atoms%jri(itype)
         r2(imesh) = atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
       end do ! imesh

      ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching coefficients
      ! have no spatial dependence) and determing its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

      ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times.
      delrGrVarphiCh1(:, :, :, :) = 0.
      delrGrVarphiCh2(:, :, :, :) = 0.
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              do ichan = 1, 2
                call Derivative( grVarphiCh1(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
                call Derivative( grVarphiCh2(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
              !  write(4002, '(4i8,3f15.8)') mqn_m2PrR, lmp, ichan, imesh, atoms%rmsh(imesh, 1), grVarphiCh1(imesh, ichan, lmp, mqn_m2PrR), grVarphiCh2(imesh, ichan, lmp, mqn_m2PrR)
              end do ! ichan
            end do !iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR

      grGrtVarphiChLout(:, :) = 0
      grGrtVarphiChMout(:, :, :) = 0
      grGrtVarphiCh1(:, :, :, :, :) = 0.
      grGrtVarphiCh2(:, :, :, :, :) = 0.
      ! Not tested
      call CalcChannelsGrGrtFlpNat( atoms, itype, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, &
                          & delrGrVarphiCh1, delrGrVarphiCh2, grGrtVarphiChLout, grGrtVarphiChMout, grGrtVarphiCh1, grGrtVarphiCh2 )

 !     ! Recycling <Psi2|H|Psi0> for the calculation of <Psi0|H|Psi2> requires additional terms correcting that the kinetic energy is not
 !     ! self-adjoint in LAPW.
 !     sAdjCorrVarphiHkinGrGrtVarphi(:, :, :, :) = 0.
 !     delrGrGrtVarphiCh1(:, :, :, :) = 0.
 !     delrGrGrtVarphiCh2(:, :, :, :) = 0.
 !     do mqn_m2PrC = -1, 1
 !       do mqn_m2PrR = -1, 1
 !         lmp = 0
 !         do oqn_l = 0, atoms%lmax(itype)
 !           do mqn_m = -oqn_l, oqn_l
 !             do iradf = 1, nRadFun(oqn_l, itype)
 !               lmp = lmp + 1
 !               do ichan = 1, 4
 !                 call Derivative( grGrtVarphiCh1(:, ichan, lmp, mqn_m2PrR, mqn_m2PrC), itype, atoms, delrGrGrtVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
 !                 call Derivative( grGrtVarphiCh2(:, ichan, lmp, mqn_m2PrR, mqn_m2PrC), itype, atoms, delrGrGrtVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
 !               end do ! ichan
 !             end do ! iradf
 !           end do ! mqn_m
 !         end do ! oqn_l
 !         call CalcSelfAdjCorrection( atoms, dimens, 4, nRadFunMax, itype, mqn_m2PrR, nRadFun, grGrtVarphiChLout,                  &
 !           & grGrtVarphiChMout(:, :, mqn_m2PrC), varphi1, varphi2, delrVarphi1, delrVarphi2, grGrtVarphiCh1(:, :, :, :, mqn_m2PrC), &
 !           & grGrtVarphiCh2(:, :, :, :, mqn_m2PrC), delrGrGrtVarphiCh1, delrGrGrtVarphiCh2, &
 !           & sAdjCorrVarphiHkinGrGrtVarphi(:, :, mqn_m2PrR, mqn_m2PrC) )
 !       end do ! mqn_m2PrR
 !     end do ! mqn_m2PrC

      vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
      ptsym = atoms%ntypsy(iatom)
      do ilh = 0, lathar%nlh(ptsym)
        oqn_l = lathar%llh(ilh, ptsym)
        lm_pre = oqn_l * (oqn_l + 1)
        do imem = 1, nmem_atom(ilh, iatom)
          mqn_m = mlh_atom(imem, ilh, iatom)
          lm = lm_pre + mqn_m
          !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
          ! maybe construct a pointer and run only over them to make it memory efficient.
          do imesh = 1, atoms%jri(itype)
            vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + Veff0%vr(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
          end do ! imesh
        end do ! imem
      end do ! ilh
      hVarphi = cmplx(0.0, 0.0)
      vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
      r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
      call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1), clnu_atom, &
        & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi )

      if (tensGradOvlOn) then
        do mqn_m2PrC = -1, 1
          call CalcVecBasfMatElems( atoms, itype, 4, nRadFun, r2, grGrtVarphiChLout, grGrtVarphiChMout(:, :, mqn_m2PrC),  &
            & varphi1, varphi2, grGrtVarPhiCh1(:, :, :, :, mqn_m2PrC), grGrtVarPhiCh2(:, :, :, :, mqn_m2PrC), hVarphi,        &
            & grGrtVarphiVarphi(:, :, :, mqn_m2PrC, itype), grGrtVarphiHVarphi(:, :, :, mqn_m2PrC, iatom) )

          call CalcGrVarphiHepsGrtVarphiElem( atoms, itype, iatom, mqn_m2PrC, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, r2, &
            & grVarPhiCh1, grVarPhiCh2, El, lmpT, vEff0NsphGrVarphi, r2grVeff0SphVarphi, grVarphiGrtVarphi, grVarphiHpreGrtVarphi, &
            & grVarphiGrtVeff0SphVarphi)
        end do ! mqn_m2PrC
      end if


!
!    allocate(r2Veff0SphLh(atoms%jmtd, 0:lathar%nlhd, atoms%ntype))
!    r2Veff0SphLh(:, :, :) = cmplx(0.0, 0.0)
!
!    allocate( vXC0MT(atoms%jmtd, 0:lathar%nlhd, atoms%ntype, 1) )
!    vxc0mt(:, :, :, :) = cmplx(0., 0.)
!    ! XC potential in the MT read in from FLEUR
!    call fopen(1000, name='v0MTFLEUR_xc', status='old', action='read', form='unformatted')
!    read(1000) vXC0MT(:, :, :, :)
!    call fclose(1000)
!    do ilh = 0, lathar%nlhd
!    do imesh = 1, atoms%jri(itype)
!      r2Veff0SphLh(imesh, 0, itype) = vEff0Mtlh(imesh, 0, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
      !r2Veff0SphLh(imesh, ilh, itype) = vEff0Mtlh(imesh, ilh, itype) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!      r2Veff0SphLh(imesh, ilh, itype) = vxc0mt(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!      r2Veff0SphLh(imesh, ilh, itype) = (Veff0%vr(imesh, ilh, itype, 1) - vxc0mt(imesh, ilh, itype, 1)) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
!    end do
!    end do

!    allocate( r2grVeff0SphSh( atoms%jmtd, (atoms%lmax(itype) + 2)**2, atoms%nat, 3 ) )
!    r2grVeff0SphSh = cmplx(0.0, 0.0)
!
!    call calcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0SphLh, r2grVeff0SphSh )
!    allocate(r2Unit(atoms%jmtd))
!    r2Unit(:) = 1.

!    do idirR = 1, 3
!      do oqn_l = 0, atoms%lmax(itype) + 1
!        do mqn_m = -oqn_l, oqn_l
!          lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
!          do imesh = 1, atoms%jmtd
!            write(4001, '(4i8,2f15.8)') idirR, oqn_l, mqn_m, imesh, r2grVeff0SphSh(imesh, lm, 1, idirR)
!          end do ! imesh
!        end do ! mqn_m
!      end do ! oqn_l
!    end do ! idir
!
!    allocate(grVeff0MTtest(atoms%jmtd,(atoms%lmax(itype) + 2)**2, 3 ))
!    grVeff0MTtest(:, :, :) = cmplx(0., 0.)
!    do idirR = 1, 3
!            do oqn_l = 0, atoms%lmax(itype)
!              lm_pre = oqn_l * (oqn_l + 1) + 1
!              do mqn_m = -oqn_l, oqn_l
!                lm = lm_pre + mqn_m
!                do imesh = 1, atoms%jri(itype)
!                  !grVeff0MTtest(imesh, lm, idirR) = r2grVeff0SphSh(imesh, lm, iatom, idirR) / atoms%rmsh(imesh, itype)**2
!                  grVeff0MTtest(imesh, lm, idirR) = grVeff0MT(imesh, lm, idirR, iatom)
!                end do ! imesh
!              end do ! mqn_m
!            end do ! oqn_l
!      end do ! idirR
!    !todo before commnit: is it critical that lmPr is running until lmax + 2
      do idirC = 1, 3
        grtVeff0MtSphVarphi(:, :, :, :) = cmplx(0., 0.)
        call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, grVeff0MT(:, :, idirC, iatom), &
!        call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, r2grVeff0SphSh(:, :, iatom, idirC), &
!        call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, grVeff0MTtest(:, :, idirC), &
                                                                                                & grtVeff0MtSphVarphi(:, :, :, :) )

!    do lmp = 1, lmpT(itype)
!      write(*, *) lmp
!      do oqn_l = 0, atoms%lmax(itype) + 2
!        lm_pre = oqn_l * (oqn_l + 1)
!        do mqn_m = -oqn_l, oqn_l
!          lm = lm_pre + mqn_m
!          do imesh = 1, atoms%jri(itype)
!            do ii = 1, 2
!              if (abs(r2grVeff0SphVarphi(ii, imesh, lm, lmp, 3)) > 1e-8 .or. abs(grtVeff0MtSphVarphi(ii, imesh, lm, lmp)) > 1e-8 ) then
!              write(4000, '(4i8,2f15.8)') ii, imesh, lm, lmp, r2grVeff0SphVarphi(ii, imesh, lm, lmp, 3)
!              write(4001, '(4i8,2f15.8)') ii, imesh, lm, lmp, grtVeff0MtSphVarphi(ii, imesh, lm, lmp)
!            end if
!            end do ! ii
!          end do ! imesh
!        end do ! mqn_m
!      end do ! oqn_l
!    end do ! lmp
!NOstopNO
        call CalcVecBasfMatElems( atoms, itype, 2, nRadFun, r2, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
!        call CalcVecBasfMatElems( atoms, itype, 2, nRadFun, r2Unit, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
!        call CalcVecBasfMatElems( atoms, itype, 2, nRadFun, r2Unit, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
          & grVarPhiCh1, grVarphiCh2, grtVeff0MtSphVarphi, varphiVarphiDummy, grVarphiGrtVeff0Varphi(:, :, :, idirC, iatom) )
      end do ! idirC

!      do idirC = 3, 3!3
!        do mqn_m2PrR = -1, 1
!      lmp = 0
!          do oqn_l = 0, atoms%lmax(itype)
!            write(*, *) oqn_l
!            do mqn_m = -oqn_l, oqn_l
!              do iradf = 1, nRadFun(oqn_l, itype)
!                lmp = lmp + 1
!                lmp1Pr = 0
!                do oqn_l1Pr = 0, atoms%lmax(itype)
!                  do mqn_m1Pr = -oqn_l1Pr, oqn_l1Pr
!                    do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
!                      lmp1Pr = lmp1Pr + 1
!                      if (abs(grVarphiGrtVeff0Varphi(lmp1Pr, lmp, mqn_m2PrR, idirC, 1)) > 1e-8 .or. abs(grVarphiGrtVeff0SphVarphi(lmp1Pr, lmp, mqn_m2PrR, idirC, 1)) > 1e-8 ) then
!                        write(4000, '(4i8,2f15.8)') lmp1Pr, lmp, mqn_m2PrR, idirC, grVarphiGrtVeff0Varphi(lmp1Pr, lmp, mqn_m2PrR, idirC, 1)
!                        write(4001, '(4i8,2f15.8)') lmp1Pr, lmp, mqn_m2PrR, idirC, grVarphiGrtVeff0SphVarphi(lmp1Pr, lmp, mqn_m2PrR, idirC, 1)
!                      end if
!                    end do ! iradf1Pr
!                  end do ! mqn_m1Pr
!                end do ! oqn_l1Pr
!              end do ! iradf
!            end do ! mqn_m
!          end do ! oqn_l
!        end do ! mqn_m2PrR
!      end do ! idirC
!      NOstopNO

       if (hfMTOn) then
         do mqn_m2PrC = -1, 1
           do idirR = 1, 3
             grVeff0GrtVarphi(:, :, :, :) = cmplx(0., 0.)
             call calcFnsphGrVarphi( atoms, itype, 0, mqn_m2PrC, nRadFun, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, &
               & grVeff0MT(:, :, idirR, iatom), grVeff0GrtVarphi(:, :, :, :) )
                                                                                                   !todo set output zero for every k for every atom
             call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, &
               & grVeff0GrtVarphi(:, :, :, :), varphiVarphiDummy, varphiGrVeff0GrtVarphi(:, :, :, idirR, mqn_m2PrC) )
           end do
         end do ! mqn_m2PrC
       end if

       if (.false.) then
         trGrGrtVharVarphi(:, :, :, :) = cmplx(0., 0.)
         call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, rho0MtSpH(:, :, iatom), &
                                                                                                  & trGrGrtVharVarphi(:, :, :, :) )

                                                                                                   !todo set output zero for every k for every atom
         call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, trGrGrtVharVarphi, varphiVarphi,    &
                                                                                                  & trvarphiGrGrtVharVarphi(:, :, :) )
       end if

       if (hfMTOn) then
         do idirC = 1, 3
           do idirR = 1, 3
             grGrtVxc0MtSphVarphi(:, :, :, :) = cmplx(0., 0.)
!             !todo this vxc0 should be over atoms but the gradient routine is not able to
             call CalcFnsphVarphi( atoms, itype, 0, nRadFun, varphi1, varphi2, grGrtVxc0MT(:, :, idirR, itype, idirC), &
                                                                                                     & grGrtVxc0MtSphVarphi(:, :, :, :) )
                                                                                                   !todo set output zero for every k for every atom
             call CalcScalBasfMatElems( atoms, itype, iatom, nRadFun, r2, varphi1, varphi2, grGrtVxc0MtSphVarphi, varphiVarphiDummy,    &
                                                                                       & varphigradGradVxcvarphi(:, :, :, idirR, idirC) )
           end do ! idirR
         end do ! idirC
       end if

!       surfInt(:, :) = cmplx(0., 0.)
!       call CalcSFintIRPsiHepsGradPsi( atoms, stars, cell, dimens, kpts, results, Veff0, ngdp, itype, iatom, nobd, eig, gBas, &
!         mapGbas, nv, gdp, z0, surfInt )
!       surfInts(1:3, 1:3, iatom) = surfInts(1:3, 1:3, iatom) + surfInt(1:3, 1:3)
!
       if (hfMTOn) then
         surfInt(:, :) = cmplx(0., 0.)
         call CalcSFintMTPsiHepsGradPsi( atoms, kpts, dimens, sym, cell, usdus, results, itype, iatom, varphi1, varphi2, &
           & nv, El, gBas, eig, lmpMax, mapGbas, nRadFun, kveclo, nobd, z0, iloTable, grVarphiChLout, grVarphiChMout, &
           & r2grVeff0SphVarphi, vEff0NsphGrVarphi, lmpT, grVarphiCh1, grVarphiCh2, surfInt )
       end if

       ! Surface integral with gradient in Ket
       surfInts(1:3, 1:3, iatom) = surfInts(1:3, 1:3, iatom) + surfInt(1:3, 1:3)


       ngoprI(:) = 1
       iqpt = 1

       psiGrGrtVxcPsiIRSum(:, :, :) = cmplx(0., 0.)
       trPsiGrGrtVharPsiIRSum(:) = cmplx(0., 0.)
       psigradGradVxcPsi(:, :, :) = cmplx(0., 0.)
       trPsiGrGrtVharPsiMTSum(:) = cmplx(0., 0.)
       rhoValMt(:) = cmplx(0., 0.)
       psiGrVeff0GrtPsiIRSum(:, :, :) = cmplx(0., 0.)
       psiGrVeff0GrtPsiMTNat(:, :) = cmplx(0., 0.)
       sAdjCorrPsiHkinGrGrtPsiNatNat(:, :) = cmplx(0., 0.)
       grPsiGrtVeff0PsiNat(:, :) = cmplx(0., 0.)
       grPsiHepsGrtPsiNatNat(:, :) = cmplx(0., 0.)
       grPsiGrtVeff0SphPsiNat(:, :) = cmplx(0., 0.)
       grPsiHepsGrtPsiNat(:, :) = cmplx(0., 0.)
       grPsiHepsGrtPsi(:, :, :) = cmplx(0., 0.)
       grGrtPsiHepsPsiNatNat(:, :) = cmplx(0., 0.)

       do ikpt = 1, kpts%nkpt

         nmat = nv(1, ikpt) + atoms%nlotot

!         gradz0(:, :, :) = cmplx(0., 0.)
!         kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
!         gExt(:) = 0.
!         do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
!           gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
!           do iband = 1, nobd(ikpt, 1)
!             do idirR = 1, 3
!               gradz0(iBas, iband, idirR) = iu * ( kExt(idirR) + gExt(idirR) ) * z0(iBas, iband, ikpt, 1)
!             end do ! idirR
!           end do ! iband
!         end do ! iBas
!
!         do idirC = 1, 3
!           do idirR = 1, 3
!             grPsiGrtVeffPsiIR(:, :) = cmplx(0., 0.)
!             call calcMEPotIR( input, stars, dimens, gBas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
!               & gBas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, w_grVeff0IR(:, idirC), z0(:, :, ikpt, 1), &
!               & gradz0(:, :, idirR), gdp, nmat, nobd(ikpt, 1), nobd(ikpt, 1), 1, ikpt, iqpt, ikpt, ngdp, grPsiGrtVeffPsiIR, kpq2kPrVec ) !todo spin-relevant
!             do iband = 1, nobd(ikpt, 1)
!               grPsiGrtVeff0PsiIRSum(idirR, idirC, iatom) = grPsiGrtVeff0PsiIRSum(idirR, idirC, iatom) + results%w_iks(iband, ikpt, 1) * grPsiGrtVeffPsiIR(iband, iband)
!             end do ! iband
!           end do ! idirR
!         end do ! idirC
!
!         do iband = 1, nobd(ikpt, 1)
!           do idirC = 1, 3
!             do idirR = 1, 3
!               grPsiHepsgrtPsiIR = cmplx(0., 0.)
!               !todo is the kinetic energy correct?
!               call calcPsi1HepsPsi1IR( kpts, qpts, stars, dimens, cell, gbas(:, mapGbas(:nv(1, ikpt), ikpt, 1)), nv, &
!                 & gradz0(:nv(1, ikpt), iband, idirR), gradz0(:nv(1, ikpt), iband, idirC), nmat, ikpt, 1, ikpt, kpq2kPrVec, &
!                 & Veff0%vpw, eig, iband, grPsiHepsGrtPsiIR )
!               grPsiHepsGrtPsiIRSum(idirR, idirC, iatom) = grPsiHepsGrtPsiIRSum(idirR, idirC, iatom) + results%w_iks(iband, ikpt, 1) * grPsiHepsGrtPsiIR
!             end do ! idirR
!           end do ! idirC
!         end do ! iband
!
!         do idirC = 1, 3
!           do idirR = 1, 3
!             psiGrVeffGrtPsiIR(:, :) = cmplx(0., 0.)
!             call calcMEPotIR( input, stars, dimens, gBas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
!               & gBas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, w_grVeff0IR(:, idirR), z0(:, :, ikpt, 1), &
!               & gradz0(:, :, idirC), gdp, nmat, nobd(ikpt, 1), nobd(ikpt, 1), 1, ikpt, iqpt, ikpt, ngdp, psiGrVeffGrtPsiIR, kpq2kPrVec ) !todo spin-relevant
!             do iband = 1, nobd(ikpt, 1)
!               psiGrVeff0GrtPsiIRSum(idirR, idirC, iatom) = psiGrVeff0GrtPsiIRSum(idirR, idirC, iatom) + results%w_iks(iband, ikpt, 1) * psiGrVeffGrtPsiIR(iband, iband)
!             end do ! iband
!           end do ! idirR
!         end do ! idirC
!
!         trPsiGrGrtVharPsiIR(:, :) = cmplx(0., 0.)
!         call calcMEPotIR( input, stars, dimens, gBas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
!           & gBas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, w_rho0IRpw, z0(:, :, ikpt, 1), &
!           & z0(:, :, ikpt, 1), gdp, nmat, nobd(ikpt, 1), nobd(ikpt, 1), 1, ikpt, iqpt, ikpt, ngdp, trPsiGrGrtVharPsiIR, kpq2kPrVec ) !todo spin-relevant
!         do iband = 1, nobd(ikpt, 1)
!           !todo that should be atom, but the potential is not atom
!           trPsiGrGrtVharPsiIRSum(iatom) = trPsiGrGrtVharPsiIRSum(iatom) + results%w_iks(iband, ikpt, 1) * trPsiGrGrtVharPsiIR(iband, iband)
!         end do ! iband
!
!         do idirC = 1, 3
!          do idirR = 1, 3
!             psiGrGrtVxcPsiIR(:, :) = cmplx(0., 0.)
!             call calcMEPotIR( input, stars, dimens, gBas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), &
!               & gBas(:, mapGbas(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, w_grGrtVxcIR(:,idirR, idirC), z0(:, :, ikpt, 1), &
!               & z0(:, :, ikpt, 1), gdp, nmat, nobd(ikpt, 1), nobd(ikpt, 1), 1, ikpt, iqpt, ikpt, ngdp, psiGrGrtVxcPsiIR, kpq2kPrVec ) !todo spin-relevant
!
!             do iband = 1, nobd(ikpt, 1)
!               !todo that should be atom, but the potential is not atom
!               psiGrGrtVxcPsiIRSum(itype, idirR, idirC) = psiGrGrtVxcPsiIRSum(itype, idirR, idirC) &
!                    & + results%w_iks(iband, ikpt, 1) * psiGrGrtVxcPsiIR(iband, iband)
!             end do ! iband
!           end do ! idirR
!         end do ! idirC
!
!         grGrTz0(:, :, :) = cmplx(0., 0.)
!         kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
!         do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
!           Gext(1:3) = matmul( cell%bmat(1:3, 1:3), gbas(1:3, mapGbas(iBas, ikpt, 1)))
!           do idirC = 1, 3
!             do idirR = 1, 3
!               grGrTz0(iBas, idirR, idirC) = - ( kExt(idirR) + gExt(idirR) ) * ( kExt(idirC) + gExt(idirC) )          &
!                                                                                                       * z0(iBas, iband, ikpt, 1)
!             end do ! idirR
!           end do ! idir
!         end do ! iBas
!         ! Calculation of the IR quantities
!
!         do iband = 1, nobd(ikpt, 1)
!           do idirC = 1, 3
!             do idirR = 1, 3
!               psiHepsgrGrtPsiIR = cmplx(0., 0.)
!               !todo is the kinetic energy correct?
!               call calcPsi1HepsPsi1IR( kpts, qpts, stars, dimens, cell, gbas(:, mapGbas(:nv(1, ikpt), ikpt, 1)), nv, &
!                 & z0(:nv(1, ikpt), iband, ikpt, 1), grGrTz0(:nv(1, ikpt), idirR, idirC), nmat, ikpt, 1, ikpt, kpq2kPrVec, &
!                 & Veff0%vpw, eig, iband, psiHepsgrGrtPsiIR )
!               psiHepsgrGrtPsiIRSum(idirR, idirC, iatom) = psiHepsgrGrtPsiIRSum(idirR, idirC, iatom) + results%w_iks(iband, ikpt, 1) * psiHepsgrGrtPsiIR
!             end do ! idirR
!           end do ! idirC
!         end do ! iband

         a(:, :, :) = cmplx(0.0, 0.0)
         b(:, :, :) = cmplx(0.0, 0.0)
         bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
         call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
           & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
           & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
           & gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
           & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
           & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

         do iband = 1, nobd(ikpt, 1)
           ab0cofScl(:) = cmplx(0., 0.)
           lmp = 0
           lm  = 0
           do oqn_l = 0, atoms%lmax(itype)
             do mqn_m = - oqn_l, oqn_l
               pMaxLocal = nRadFun(oqn_l, itype)
               ! p = 1
               ab0cofScl(lmp + 1) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iatom) )
               ! p = 2
               ab0cofScl(lmp + 2) = iu**oqn_l * dot_product( conjg(z0(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iatom) )
               ! Add LO contributions
               do iradf = 3, pMaxLocal
                 ! p = 1
                 ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + iu**oqn_l * &
                   & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                 ! p = 2
                 ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + iu**oqn_l * &
                   & dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
                 ! 2 < p < LOs for that l and that atom type
                 ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z0(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, itype), iatom) )
               end do ! iradf


               ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
               ! least the p as lm1Pr = lm. But for sake of performance we place it here.
               ! sake of performance
               lm = lm + 1
               lmp = lmp + pMaxLocal
             end do ! mqn_m
           end do !oqn_l


           do idirC = 1, 3
             do mqn_m2PrR = -1, 1
               grVarphiGrtVeff0PsiNat(1:lmpT(itype)) = cmplx(0., 0.)
               grVarphiGrtVeff0PsiNat(1:lmpT(itype)) = matmul( grVarphiGrtVeff0Varphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, idirC, iatom), ab0cofScl(1:lmpT(itype)) )
               grPsiGrtVeff0PsiNat(mqn_m2PrR, idirC) = grPsiGrtVeff0PsiNat(mqn_m2PrR, idirC) + results%w_iks(iband, ikpt, 1) * dot_product( ab0cofScl(1:lmpT(itype)), grVarphiGrtVeff0PsiNat(1:lmpT(itype)) )

            ! write(4001, '(4i8,2f15.8)') ikpt, iband, mqn_m2PrR, idirC, grPsiGrtVeff0PsiNat(mqn_m2PrR, idirC)
             end do ! mqn_m2PrR
           end do ! idirC

           do mqn_m2PrC = -1, 1
             do mqn_m2PrR = -1, 1
                grVarphiHepsGrtPsiNatNat(:) = cmplx(0., 0.)
                grVarphiHepsGrtPsiNatNat(1:lmpT(itype)) = &
         !         &   matmul(grVarphiHpreGrtVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC, iatom),          &
         !                                                                                   & ab0cofScl(1:lmpT(itype)) )&
         !         & - eig(iband, ikpt, 1) * matmul(grVarphiGrtVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC, itype)&
         !                                                                                 &  , ab0cofScl(1:lmpT(itype)) )
                  & matmul(grVarphiGrtVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC, itype)&
                                                                                          &  , ab0cofScl(1:lmpT(itype)) )

                grPsiHepsGrtPsiNatNat(mqn_m2PrR, mqn_m2PrC) = grPsiHepsGrtPsiNatNat(mqn_m2PrR, mqn_m2PrC) +                       &
                  & results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), grVarphiHepsGrtPsiNatNat(1:lmpT(itype)) )

                grVarphiGrtVeff0SphPsiNat(:) = cmplx(0., 0.)
                grVarphiGrtVeff0SphPsiNat(1:lmpT(itype)) = &
                  & matmul( grVarphiGrtVeff0SphVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC + 2, iatom), &
                                                                                             & ab0cofScl(1:lmpT(itype)) )

                grPsiGrtVeff0SphPsiNat(mqn_m2PrR, mqn_m2PrC + 2) = grPsiGrtVeff0SphPsiNat(mqn_m2PrR, mqn_m2PrC+2) + &
                  & results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), grVarphiGrtVeff0SphPsiNat(1:lmpT(itype)) )
!                write(4000, '(4i8,2f15.8)') ikpt, iband, mqn_m2PrR, mqn_m2PrC + 2, grPsiGrtVeff0SphPsiNat(mqn_m2PrR, mqn_m2PrC + 2)
             end do ! mqn_m2PrR
           end do ! mqn_m2PrC

           do mqn_m2PrC = -1, 1
             do mqn_m2PrR = -1, 1
!               ! Calculation of the last line in A.50 PhD thesis A. Klueppelberg
!               grGrtVarphiHPsiNatNat(:) = cmplx(0.0, 0.0)
!               grGrtVarphiHPsiNatNat(1:lmpT(itype)) = matmul(grGrtVarphiHvarphi(:lmpT(itype), :lmpT(itype), mqn_m2PrR, mqn_m2PrC, iatom), ab0cofScl(:lmpT(itype)))
!               ! Precalculation of the last line in A.50 PhD thesis A. Klueppelberg.
               grGrtVarphiPsiNatNat(:) = cmplx(0.0, 0.0)
               grGrtVarphiPsiNatNat(1:lmpT(itype)) = &
                     & matmul( grGrtVarphiVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC, itype), &
                                                                                                     & ab0cofScl(1:lmpT(itype)) )
               grGrtPsiHepsPsiNatNat(mqn_m2PrR, mqn_m2PrC) = grGrtPsiHepsPsiNatNat(mqn_m2PrR, mqn_m2PrC) &
!                 & + results%w_iks(iband, ikpt, 1) * ( dot_product(ab0cofScl(:lmpT(itype)), grGrtVarphiHPsiNatNat(1:lmpT(itype))) &
!                           & - eig(iband, ikpt, 1) * dot_product(ab0cofScl(:lmpT(itype)), grGrtVarphiPsiNatNat(1:lmpT(itype)) ) )
                           & + results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(:lmpT(itype)), grGrtVarphiPsiNatNat(1:lmpT(itype)) )
!               sAdjCorrPsiGrGrtVarphiNat(:) = cmplx(0., 0.)
!               sAdjCorrPsiGrGrtVarphiNat(1:lmpT(itype)) = &
!                 & matmul( sAdjCorrVarphiHkinGrGrtVarphi(1:lmpT(itype), 1:lmpT(itype), mqn_m2PrR, mqn_m2PrC), &
!                                                                                              & conjg(ab0cofScl(1:lmpT(itype))) )
!               ! Precalculation for the correction of <varphi|H_kin|grGrtvarphi> if we complex conjugate it for calculating A.51 PhD
!               ! thesis A. Klueppelberg.
!               sAdjCorrPsiHkinGrGrtPsiNatNat(mqn_m2PrR, mqn_m2PrC) = sAdjCorrPsiHkinGrGrtPsiNatNat(mqn_m2PrR, mqn_m2PrC) + &
!            & results%w_iks(iband, ikpt, 1) * dot_product( conjg(ab0cofScl(:lmpT(itype))), sAdjCorrPsiGrGrtVarphiNat(:lmpT(itype)) )
             end do ! mqn_m2PrR
           end do ! mqn_m2PrC
           if (hfMTOn) then
             do mqn_m2PrC = -1, 1
               do idirR = 1, 3
                 varphiGrVeff0GrtPsi(:) = cmplx(0., 0.)
                 varphiGrVeff0GrtPsi(1:lmpT(itype)) = matmul( varphiGrVeff0GrtVarphi(1:lmpT(itype), 1:lmpT(itype), iatom, idirR, mqn_m2PrC), ab0cofScl(1:lmpT(itype)) )
                 PsiGrVeff0GrtPsiMTNat(idirR, mqn_m2PrC) = PsiGrVeff0GrtPsiMTNat(idirR, mqn_m2PrC) + results%w_iks(iband, ikpt, 1) * dot_product(ab0cofScl(1:lmpT(itype)), varphiGrVeff0GrtPsi(1:lmpT(itype)))
               end do ! idirR
             end do ! mqn_m2PrC
           end if

           if (.false.) then
             trvarphiGrGrtVharPsi(:) = cmplx(0., 0.)
             trvarphiGrGrtVharPsi(1:lmpT(itype)) = matmul(trvarphiGrGrtVharVarphi(1:lmpT(itype), 1:lmpT(itype), iatom), ab0cofScl(1:lmpT(itype)) )
             trPsiGrGrtVharPsiMTSum(iatom) = trPsiGrGrtVharPsiMTSum(iatom) + dot_product( ab0cofScl(1:lmpT(itype)), trvarphiGrGrtVharPsi(1:lmpT(itype)) )
             varphiPsi(:) = cmplx(0., 0.)
             varphiPsi(1:lmpT(itype)) = matmul( varphiVarphi(1:lmpT(itype), 1:lmpT(itype), itype), ab0cofScl(1:lmpT(itype)) )
             rhoValMt(iatom) = rhoValMt(iatom) + results%w_iks(iband, ikpt, 1) * dot_product( ab0cofScl(1:lmpT(itype)), varphiPsi(1:lmpT(itype)) )
             !todo same for varphiPsi then multiply the external ptoential tensor gradient to it
           end if

           if (hfMTOn) then
             do idirC = 1, 3
               do idirR = 1, 3
                 varphigradGradVxcPsi(:) = cmplx(0., 0.)
                 varphigradGradVxcPsi(1:lmpT(itype)) = &
                      & matmul( varphigradGradVxcvarphi(1:lmpT(itype), 1:lmpT(itype), itype, idirR, idirC), ab0cofScl(1:lmpT(itype)) )
                 psigradGradVxcPsi(itype, idirR, idirC) = psigradGradVxcPsi(itype, idirR, idirC) + results%w_iks(iband, ikpt, 1) &
                                                      & * dot_product( ab0cofScl(1:lmpT(itype)), varphigradGradVxcPsi(1:lmpT(itype)) )
               end do ! idirR
             end do ! idirC
           end if
         end do ! iband
       end do ! ikpt

       grPsiGrtVeff0Psi(1:3, 1:3, iatom) = matmul( conjg(Tmatrix(1:3, 1:3)), grPsiGrtVeff0PsiNat(-1:1, 1:3) )

       grPsiHepsGrtPsiNat(-1:1, 1:3) = matmul(grPsiHepsGrtPsiNatNat(-1:1, -1:1), Tmatrix_transposed(1:3, 1:3))
       ! grPsiHepsGrtPsiNat(:, :) = cmplx(0., 0.)
       ! grPsiHepsGrtPsiNat(-1:1, 1:3) = grPsiHepsGrtPsiNat(-1:1, 1:3) - grPsiGrtVeff0SphPsiNat(-1:1, 1:3)
       grPsiHepsGrtPsi(1:3, 1:3, iatom) = matmul(conjg(Tmatrix(1:3, 1:3)), grPsiHepsGrtPsiNat(-1:1, 1:3))

       if (hfMTOn) then
         PsiGrVeff0GrtPsiMt(1:3, 1:3, iatom) = matmul(PsiGrVeff0GrtPsiMTNat(1:3, -1:1), Tmatrix_transposed(1:3, 1:3))
         psiGrVeff0GrtPsiMTNat(:, :) = cmplx(0., 0.)
       end if

!       sAdjCorrPsiHkinGrGrtPsiNat(:, :) = cmplx(0., 0.)
!       ! This Tmatrix is not complex conjugated because it is not extracted from the bra for the correction terms
!       sAdjCorrPsiHkinGrGrtPsiNat(1:3, -1:1) = matmul( Tmatrix(1:3, 1:3), sAdjCorrPsiHkinGrGrtPsiNatNat(-1:1, -1:1) )
!       sAdjCorrPsiHkinGrGrtPsi(1:3, 1:3, itype) = matmul( sAdjCorrPsiHkinGrGrtPsiNat(1:3, -1:1), Tmatrix_transposed(1:3, 1:3) )
!
       grGrtPsiHepsPsiNat(:, :) = cmplx(0.0, 0.0)
       grGrtPsiHepsPsiNat(1:3, -1:1) = matmul(conjg(Tmatrix(1:3, 1:3)), grGrtPsiHepsPsiNatNat(-1:1, -1:1))
       grGrtPsiHepsPsi(1:3, 1:3, iatom) = matmul( grGrtPsiHepsPsiNat(1:3, -1:1), conjg(Tmatrix_transposed(1:3, 1:3)) )
      end do ! ieqat
    end do ! itype


    if (.false.) then
      allocate( r2Rho0MTsphVal(atoms%jmtd, (atoms%lmaxd + 1)**2), fReal(atoms%jmtd), fImag(atoms%jmtd) )
    end if

    traceLeft = cmplx(0., 0.)
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1

        if (.false.) then
          r2Rho0MTsphVal(:, :) = cmplx(0.0, 0.0)
          ptsym = atoms%ntypsy(iatom)
          do ilh = 0, lathar%nlh(ptsym)
            oqn_l = lathar%llh(ilh, ptsym)
            lm_pre = oqn_l * (oqn_l + 1)
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)
              lm = lm_pre + mqn_m + 1
              !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
              ! maybe construct a pointer and run only over them to make it memory efficient.
              do imesh = 1, atoms%jri(itype)
                r2Rho0MTsphVal(imesh, lm) = r2Rho0MTsphVal(imesh, lm) + r2Rho0MtlhVal(imesh, ilh, itype, 1) * clnu_atom(imem, ilh, iatom)
              end do ! imesh
            end do ! imem
          end do ! ilh

          !rho val * rho
          integral = cmplx(0., 0.)
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              fReal(:) = 0.0
              fImag(:) = 0.0
              do imesh = 1, atoms%jri(itype)
                fReal(imesh) = real( conjg(r2Rho0MTsphVal(imesh, lm)) * rho0MtSpH(imesh, lm, iatom) )
                fImag(imesh) = aimag( conjg(r2Rho0MTsphVal(imesh, lm)) * rho0MtSpH(imesh, lm, iatom) )
              end do ! imesh
              call intgr3LinIntp(fReal, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
              call intgr3LinIntp(fImag, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrImag, 1)
              integral = integral + cmplx(intgrReal, intgrImag)
            end do ! mqn_m
          end do ! oqn_l

          do idirR = 1, 3
            do oqn_l = 0, atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                fReal(:) = 0.0
                fImag(:) = 0.0
                do imesh = 1, atoms%jri(itype)
                  fReal(imesh) = real( conjg(r2Rho0MTsphVal(imesh, lm)) * grGrtVxc0MT(imesh, lm, idirR, 1, idirR) )
                  fImag(imesh) = aimag( conjg(r2Rho0MTsphVal(imesh, lm)) * grGrtVxc0MT(imesh, lm, idirR, 1, idirR) )
                end do ! imesh
                call intgr3LinIntp(fReal, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
                call intgr3LinIntp(fImag, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrImag, 1)
                integral = integral + cmplx(intgrReal, intgrImag)
              end do ! mqn_m
            end do ! oqn_l
          end do ! idirR
        end if ! end if


        allocate(fReal(atoms%jmtd), fImag(atoms%jmtd))
        integralMat(:, :) = cmplx(0., 0.)
        do idirC = 1, 3
          do idirR = 1, 3
            do oqn_l = 0, 0!atoms%lmax(itype)
              do mqn_m = -oqn_l, oqn_l
                lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
                fReal(:) = 0.0
                fImag(:) = 0.0
                do imesh = 1, atoms%jri(itype)
                  write(4003, '(5i8,2f15.8)') idirC, idirR, oqn_l, mqn_m, imesh, atoms%rmsh(imesh, itype)**2 * grGrtRho0MTval(imesh, lm, idirR, 1, idirC )
                  fReal(imesh) = real( grGrtRho0MTval(imesh, lm, idirR, 1, idirC ) * atoms%rmsh(imesh, itype)**2 )
                  fImag(imesh) = aimag( grGrtRho0MTval(imesh, lm, idirR, 1, idirC) * atoms%rmsh(imesh, itype)**2 )
                end do ! imesh
                call intgr3LinIntp(fReal, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrReal, 1)
                call intgr3LinIntp(fImag, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), intgrImag, 1)
                integralMat(idirR, idirC) = integralMat(idirR, idirC) + 0.25 * sqrt(fpi) * cmplx(intgrReal, intgrImag)
              end do ! mqn_m
            end do ! oqn_l
          end do ! idirR
        end do ! idirC

! todo times 2 because of spin
! todo core?
!        write(*, '(a)') 'grPsiGrtVeff0PsiIR'
!        write(*, '(3(2(es16.8,1x),3x))') transpose(grPsiGrtVeff0PsiIRSum(:, :, iatom))
!        singleTrace = cmplx(0., 0.)
!        do ii = 1, 3
!          singleTrace = singleTrace + grPsiGrtVeff0PsiIRSum(ii, ii, iatom)
!        end do ! ii
!        write(*, '(a)') 'TRACE grPsiGrtVeff0PsiIR'
!        write(*, '(2f15.8)') singleTrace
!        traceLeft = traceLeft + singleTrace
!        write(*, *)
!
!        write(*, '(a)') 'grPsiGrtVeff0PsiMT'
!        write(*, '(3(2(es16.8,1x),3x))') transpose(grPsiGrtVeff0Psi(:, :, iatom))
!        singleTrace = cmplx(0., 0.)
!        do ii = 1, 3
!          singleTrace = singleTrace + grPsiGrtVeff0Psi(ii, ii, iatom)
!        end do ! ii
!        write(*, '(a)') 'TRACE grPsiGrtVeff0PsiMT'
!        write(*, '(2f15.8)') singleTrace
!        traceLeft = traceLeft + singleTrace
!        write(*, *)

!        write(*, '(a)') 'grPsiHepsGrtPsiIR'
!        write(*, '(3(2(es16.8,1x),3x))') transpose(grPsiHepsGrtPsiIRSum(:, :, iatom))
!        singleTrace = cmplx(0., 0.)
!        do ii = 1, 3
!          singleTrace = singleTrace + grPsiHepsGrtPsiIRSum(ii, ii, iatom)
!        end do ! ii
!        write(*, '(a)') 'TRACE grPsiHepsGrtPsiIR'
!        write(*, '(2f15.8)') singleTrace
!        traceLeft = traceLeft + singleTrace
!        write(*, *)
!
        write(*, *)
        write(*, '(a)') 'grPsiGrtPsiMT'
        write(*, '(3(2(es16.8,1x),3x))') transpose(grPsiHepsGrtPsi(:, :, iatom))
        write(*, *)
        singleTrace = cmplx(0., 0.)
        do ii = 1, 3
          singleTrace = singleTrace + grPsiHepsGrtPsi(ii, ii, iatom)
        end do ! ii
!        write(*, '(a)') 'TRACE grPsiHepsGrtPsiMT'
!        write(*, '(2f15.8)') singleTrace
        traceLeft = traceLeft + singleTrace
!        write(*, *)
!
!        write(*, '(a)') 'psiGrGrtVxcPsiIR'
!        write(*, '(3(2(es16.8,1x),3x))') 2 * transpose(psiGrGrtVxcPsiIRSum(itype, :, :))
!        singleTrace = cmplx(0., 0.)
!        do ii = 1, 3
!          singleTrace = singleTrace + 2 * psiGrGrtVxcPsiIRSum(itype, ii, ii)
!        end do ! ii
!        write(*, '(a)') 'TRACE psiGrGrtVxcPsiIR'
!        write(*, '(2f15.8)') singleTrace
!        traceLeft = traceLeft + singleTrace
!        write(*, *)
         write(*, *) 'Integral over r^2 grad gradT rho'
         write(*, '(3(2(es16.8,1x),3x))') integralMat(:, :)
          singleTrace = cmplx(0., 0.)
          do ii = 1, 3
            singleTrace = singleTrace + integralMat(ii, ii)
          end do ! ii
         write(*, '(a)') 'trace of upper matrix'
         write(*, '(2f15.8)') singleTrace
        write(*, *)
!
        if (hfMTOn) then
          write(*, '(a)') 'psiGrGrtVxcPsiMT'
          write(*, '(3(2(es16.8,1x),3x))') 2 * transpose(psigradGradVxcPsi(itype, :, :))
          singleTrace = cmplx(0., 0.)
          do ii = 1, 3
            ! 2 for spin degeneracy
            singleTrace = singleTrace + 2 * psigradGradVxcPsi(itype, ii, ii)
          end do ! ii
          write(*, '(a)') 'TRACE psiGrGrtVxcPsiMT'
          write(*, '(2f15.8)') singleTrace
          traceLeft = traceLeft + singleTrace
          write(*, *)
        end if

!        write(*, '(a)') 'trpsiGrGrtVharPsiIR'
!        write(*, '(2f15.8)') trPsiGrGrtVharPsiIRSum(iatom)
!        write(*, '(a)') 'TRACE trpsiGrGrtVharPsiIR'
!        write(*, '(2f15.8)') trPsiGrGrtVharPsiIRSum(iatom)
!        traceLeft = traceLeft - fpi * trPsiGrGrtVharPsiIRSum(iatom)
!        write(*, *)

        ! I comment this out because, one gets better results if we add the hartree weinert potential to the grad gradT xc algorithm
        if (.false.) then
          write(*, '(a)') 'TRACE trpsiGrGrtVharPsiMT'
          !write(*, '(2f15.8)') - fpi * trPsiGrGrtVharPsiMTSum(iatom)
          !write(*, '(2f15.8)') - fpi * integral
          ! the upper is if we use rho0Full, the lower is if we use grad gradT Vext
          write(*, '(2f15.8)') integral
          write(*, *)
        end if

        if (hfMTOn) then
          ! Note: This together with the psiGrGrtVxcPsiMT if there the complete effective potential is a gradient applied to, gives
          ! a quiet good result. The question is in how far a uncontinious quantity is okay for the comparison.
          ! An idea might also to put this into HF and compensate the discontinuity with a surface integral, still the integration
          ! va;ues are a bit different because one funcitions strives for continuity the other not.
          write(*, '(a)') 'Trace trpsiGrGrtVextPsi'
          write(*, '(2f15.8)') sqrt(fpi) * atoms%zatom(itype) * r2Rho0MtLhVal(1, 0, itype, 1) * clnu_atom(1, 0, iatom) / atoms%rmsh(1, itype)**2
          traceLeft = traceLeft + sqrt(fpi) * atoms%zatom(itype) * r2Rho0MtLhVal(1, 0, itype, 1) * clnu_atom(1, 0, iatom) / atoms%rmsh(1, itype)**2
          write(*, *)
        end if

!        write(*, '(a)') 'psiGrVeff0GrtPsiIR'
!        write(*, '(3(2(es16.8,1x),3x))') 4 * transpose(psiGrVeff0GrtPsiIRSum(:, :, iatom))
!        singleTrace = cmplx(0., 0.)
!        do ii = 1, 3
!          singleTrace = singleTrace + 4 * psiGrVeff0GrtPsiIRSum(ii, ii, iatom)
!        end do ! ii
!        write(*, '(a)') 'TRACE psiGrVeff0GrtPsiIR'
!        write(*, '(2f15.8)') singleTrace
!        traceLeft = traceLeft + singleTrace
!        write(*, *)

        if (hfMTOn) then
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Note!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!:
          ! this integral is only the same as in HF, when we leave out the first sixth mesh points. consistently as in HF. It
          ! does not depend on the integration method i.e. the interpolation done. This due to the fact that here, the lmp lmp1Pr
          ! are integrated step by step, which does not result in that large integrals. After that the abcofs are multiplied and there
          ! is a sum over the k-points and the bands. This is a difference to directly evaluating the integral of rho1 times vext1,
          ! where a sum over k-points and mesh points and the spin factor is already within. So Integration of the sum is something
          ! different than integration of every summand. We should note however, that leaving away the first 6 mesh points causes
          ! a difference of 1e-2
          write(*, '(a)') 'psiGrVeff0GrtPsiMt'
          ! Note: We need a factor of 4 to compensate the spin factor of the valence density from FLEUR.
          write(*, '(3(2(es16.8,1x),3x))') 4 * transpose(PsiGrVeff0GrtPsiMt(:, :, iatom))
          singleTrace = cmplx(0., 0.)
          do ii = 1, 3
            singleTrace = singleTrace + 4 * PsiGrVeff0GrtPsiMt(ii, ii, iatom)
          end do ! ii
          write(*, '(a)') 'TRACE psiGrVeff0GrtPsiMt'
          write(*, '(2f15.8)')  singleTrace
          traceLeft = traceLeft +  singleTrace
          write(*, *)
        end if

!
!        write(*, '(a)') 'psiHepsGrGrtPsiIR'
!        write(*, '(3(2(es16.8,1x),3x))') transpose( psiHepsgrGrtPsiIRSum(:, :, iatom) )
!        singleTrace = cmplx(0., 0.)
!        do ii = 1, 3
!          singleTrace = singleTrace + psiHepsgrGrtPsiIRSum(ii, ii, iatom)
!        end do ! ii
!        write(*, '(a)') 'TRACE psiHepsGrGrtPsiIR'
!        write(*, '(2f15.8)') singleTrace
!        traceLeft = traceLeft + singleTrace
!        write(*, *)
!
          write(*, *)
        write(*, '(a)') 'psiGrGrtPsiMT'
!        write(*, '(3(2(es16.8,1x),3x))') transpose( conjg(grGrtPsiHepsPsi(:, :, iatom)) + sAdjCorrPsiHkinGrGrtPsi(:, :, iatom) )
        write(*, '(3(2(es16.8,1x),3x))') transpose( conjg(grGrtPsiHepsPsi(:, :, iatom)))
          write(*, *)
        singleTrace = cmplx(0., 0.)
        do ii = 1, 3
          singleTrace = singleTrace + conjg(grGrtPsiHepsPsi(ii, ii, iatom)) + sAdjCorrPsiHkinGrGrtPsi(ii, ii, iatom)
        end do ! ii
!        write(*, '(a)') 'TRACE psiHepsGrGrtPsiMT'
!        write(*, '(2f15.8)') singleTrace
        traceLeft = traceLeft + singleTrace
!        write(*, *)
!
!        write(*, '(a)') 'Surfints'
!        write(*, '(3(2(es16.8,1x),3x))') transpose( surfInts(:, :, iatom) )
!        singleTrace = cmplx(0., 0.)
!        do ii = 1, 3
!          singleTrace = singleTrace + surfInts(ii, ii, iatom)
!        end do ! ii

        write(*, *)
        write(*, *)
        write(*, '(a)') 'LeftSide'
        write(*, '(2f15.8)') traceLeft

        write(*, '(a)') 'Surfints'
        write(*, '(2f15.8)') singleTrace
      end do ! ieqat
    end do ! itype

  end subroutine TestTensorGradProdKet

  ! Test the correction term that Aaron suggests to compensate the non-hermiticiy of the Hamiltonian due to the kinetic energy.
  ! The problem with that is that we have decided to use an unsymmetric version of the kientic energy in the IR applied to the
  ! ket. Therefore, we have to also apply the kinetic energy in the MT to the ket, i.e. the spherical Hamiltonian in total.
  ! Doing so leads to some side effects so that, to use Aarons equation, one would have to compensate the difference between
  ! <grPsi|Hsph|Psi> and <Psi|Hsph|GrPsi> and then can also directly calculate <Psi|Hsph|GrPsi> according to the trick in A.54 and
  ! A.55.
  ! TODO Make this test nice, in principle it works but has the problem of the unsymmetric kinetic energy
  subroutine testGradBraCorrection(atoms, dimens, lathar, Veff0, kpts, sym, cell, usdus, rbas1, rbas2, nRadFun, El, mlh_atom, nmem_atom, clnu_atom, nv, mapGbas, gBas, kveclo, nobd, z, iloTable)

    use m_types,         only : t_atoms, t_dimension, t_sphhar, t_potential, t_kpts, t_sym, t_cell, t_usdus
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat
    use m_jpSetupDynMat, only : CalcGrVarphiHepsGrtVarphiElem, CalcVecBasfMatElems, CalcScalBasfMatElems, CalcSelfAdjCorrection
    use m_jpSetupDynMatSF, only : CalcHGrVarphi
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi
    use m_jpConstants, only : iu, Tmatrix_transposed
    use m_od_types, only : od_inp, od_sym
    use m_abcof3


    implicit none

    ! Type parameter
    type(t_atoms),              intent(in)  :: atoms
    type(t_dimension),          intent(in)  :: dimens
    type(t_sphhar),             intent(in)  :: lathar
    type(t_potential),          intent(in)  :: Veff0
    type(t_kpts),                   intent(in)    :: kpts
    type(t_sym),                intent(in) :: sym
    type(t_cell),                   intent(in)    :: cell
    type(t_usdus),                  intent(in)    :: usdus

    ! Array parameters
    real,                       intent(in)  :: rbas1(:, :, 0:, :)
    real,                       intent(in)  :: rbas2(:, :, 0:, :)
    integer,                    intent(in)  :: nRadFun(0:, :)
    real,                       intent(in)  :: El(:, 0:, :, :)
    integer,                    intent(in)  :: mlh_atom(:, 0:, :)
    integer,                    intent(in)  :: nmem_atom(0:, :)
    complex,                    intent(in)  :: clnu_atom(:, 0:, :)
    integer,                        intent(in)    :: nv(:, :)
    integer,                        intent(in)    :: mapGbas(:, :, :)
    integer,                        intent(in)    :: gBas(:, :)
    integer,                        intent(in)    :: kveclo(:,:)
    integer,                        intent(in)  :: nobd(:, :)
    MCOMPLEX,                       intent(in)    :: z(:,:,:,:)
    integer,                        intent(in)    :: iloTable(:, 0:, :)

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods

    ! Scalar variables
    integer                                       :: nmat
    integer                                 :: oqn_l
    integer                                 :: iradf
    integer                                 :: imesh
    integer                                 :: lmpMax
    integer                                 :: nRadFunMax
    integer                                 :: lmaxBra
    integer                                 :: iDeqat
    integer                                 :: iDtype
    integer                                 :: iDatom
    integer                                 :: ptsym
    integer                                 :: ilh
    integer                                 :: lm_pre
    integer                                 :: imem
    integer                                 :: mqn_m
    integer                                 :: lm
    integer                                 :: mqn_m2PrC
    integer                                 :: mqn_m2PrR
    integer                                 :: ii
    integer                                 :: jj
    integer                                 :: idir
    integer                                 :: lmp
    integer                                 :: ichan
    integer                                 :: ikpt
    integer                                 :: iband
    integer                                 :: iDatomA
    integer                                 :: iDtypeA
    integer                                 :: iDeqatA
    integer                                       :: pMaxLocal
    integer                                 :: idirR
    integer                                 :: idirC
    integer                                :: iBas

    ! Array variables
    integer,       allocatable              :: lmpT(:)
    real,          allocatable              :: varphi1(:, :, :)
    real,          allocatable              :: varphi2(:, :, :)
    real,          allocatable              :: delrVarphi1(:, :, :)
    real,          allocatable              :: delrVarphi2(:, :, :)
    integer,       allocatable              :: grVarphiChLout(:, :)
    integer,       allocatable              :: grVarphiChMout(:, :)
    real,          allocatable              :: grVarphiCh1(:, :, :, :)
    real,          allocatable              :: grVarphiCh2(:, :, :, :)
    complex,       allocatable              :: vEff0MtSpH(:, :)
    real,          allocatable              :: grVarphiGrtVarphi(:, :, :, :, :)
    complex,       allocatable              :: grVarphiHpreGrtVarphi(:, :, :, :, :)
    complex,       allocatable              :: grVarphiGrtVeff0SphVarphibench(:, :, :, :, :)
    complex,       allocatable              :: grVarphiGrtVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable              :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,       allocatable              :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,       allocatable              :: hVarphi(:, :, :, :)
    complex,       allocatable              :: hGrVarphi(:, :, :, :, :)
    complex,       allocatable              :: grVarphiHGrVarphi(:, :, :, :, :)
    real,          allocatable              :: varphiGrVarphi(:, :, :, :)
    real,          allocatable              :: grVarphiVarphi(:, :, :, :)
    complex,       allocatable              :: grVarphiHVarphiDummy(:, :, :, :)
    real,          allocatable              :: varphiVarphiDummy(:, :, :)
    complex,       allocatable              :: varphiHGrvarphi(:, :, :, :)
    complex,       allocatable              :: varphiGrVeffvarphi(:, :, :, :)
    real,          allocatable              :: delrGrVarphiCh1(:, :, :, :)
    real,          allocatable              :: delrGrVarphiCh2(:, :, :, :)
    real,          allocatable              :: sAdjCorrVarphiHkinGrVarphi(:, :, :, :)
    complex,       allocatable              :: grVarphiHVarphi(:, :, :, :)
    complex,           allocatable                :: ab0cofScl(:)
    real,          allocatable              :: r2(:)
    complex,           allocatable          :: a(:, :, :)
    complex,           allocatable          :: b(:, :, :)
    complex,           allocatable          :: bascof_lo(:, :, :, :, :)
    integer,           allocatable          :: ngoprI(:)
    complex,           allocatable                :: ab1cofVec(:, :, :)
    complex,       allocatable              :: z1nG(:, :, :, :)
    complex,           allocatable                :: grVarphiPsi(:)
    complex, allocatable                     :: varphiGrPsi(:)
    complex, allocatable :: grVarphiHPsi(:)
    complex, allocatable :: sAdjCorrVarphiHkinGrPsi(:)
    complex, allocatable :: varphiHGrPsi(:)
    complex, allocatable :: varphiGrVeffPsi(:)

    complex                                 :: psiGrPsiNat(3,-1:1)
    complex                                 :: grPsiPsiNat(3,-1:1)

    complex                                 :: psiGrPsi(3,3)
    complex                                 :: grPsiPsi(3,3)
    real :: kExt(3)
    real :: gExt(3)
    complex :: psiGrVeffPsi(3, 3)
    complex :: psiHGrPsiNat(3, -1:1)
    complex :: grPsiHPsiNat(3, -1:1)
    complex :: sAdjCorrPsiHkinGrPsiNat(3, -1:1)
    complex :: psiHGrPsi(3, 3)
    complex :: grPsiHPsi(3, 3)
    complex :: sAdjCorrPsiHkinGrPsi(3, 3)
    write(*, *) 'myfoo'

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do iDtype = 1, atoms%ntype
      lmpT(iDtype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, iDtype), oqn_l = 0, atoms%lmax(iDtype) ) ] )
    end do ! iDtype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    ! Allocation of required arrays
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( hGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1))
    allocate( grVarphiHGrVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( r2(atoms%jmtd) )
    allocate( grVarphiGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grVarphiHpreGrtVarphi(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( grVarphiGrtVeff0SphVarphibench(lmpMax, lmpMax, -1:1, 1:3, atoms%nat) )
    allocate( grVarphiGrtVeff0SphVarphi(lmpMax, lmpMax, -1:1, 1:3, atoms%nat) )
    allocate( varphiGrVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), grVarphiHVarphiDummy(lmpMax, lmpMax, -1:1, atoms%nat ) )
    allocate( grVarphiVarphi(lmpMax, lmpMax, -1:1, atoms%ntype), grVarphiHVarphi(lmpMax, lmpMax, -1:1, atoms%nat ) )
    allocate( varphiVarphiDummy(lmpMax, lmpMax, atoms%ntype), varphiHGrvarphi(lmpMax, lmpMax, atoms%nat, -1:1))
    allocate( varphigrVeffvarphi(lmpMax, lmpMax, atoms%nat, 3))
    allocate( delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( sAdjCorrVarphiHkinGrVarphi(lmpMax, lmpMax, -1:1, atoms%ntype) )
    allocate( ab0cofScl(lmpMax) )
    allocate( grVarphiHPsi(lmpMax), sAdjCorrVarphiHkinGrPsi(lmpMax), varphiHGrPsi(lmpMax), varphiGrVeffPsi(lmpMax))

    ! Initialization
    grVarphiGrtVeff0SphVarphibench(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiGrtVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiHpreGrtVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiGrtVarphi(:, :, :, :, :) = cmplx(0., 0.)
    grVarphiHGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    hGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    r2grVeff0SphVarphi(:, :, :, :, :) = cmplx(0., 0.)
    vEff0NsphGrVarphi(:, :, :, :, :) = cmplx(0., 0.)
    varphiGrVarphi(:, :, :, :)  = 0.
    grVarphiHVarphiDummy(:, :, :, :) = cmplx(0., 0.)
    grVarphiVarphi(:, :, :, :)  = 0.
    varphiVarphiDummy(:, :, :) = 0.
    varphiHGrvarphi(:, :, :, :) = cmplx(0., 0.)
    varphigrVeffvarphi(:, :, :, :) = cmplx(0., 0.)
    grVarphiHPsi(:) = cmplx(0., 0.)
    sAdjCorrVarphiHkinGrPsi(:) = cmplx(0., 0.)
    varphiHGrPsi(:) = cmplx(0., 0.)
    varphiGrVeffPsi(:) = cmplx(0., 0.)

    iDatom = 0
    do iDtype = 1, atoms%ntype
      r2(:) = 0.
      ! Jacobi determinant of radial part
      do imesh = 1, atoms%jri(iDtype)
        r2(imesh) = atoms%rmsh(imesh, iDtype) * atoms%rmsh(imesh, iDtype)
      end do ! imesh
      do iDeqat = 1, atoms%neq(iDtype)
        iDatom = iDatom + 1

        ! Calculate the gradient of varphi
        varphi1(:, :, :) = 0.
        varphi2(:, :, :) = 0.
        delrVarphi1(:, :, :) = 0.
        delrVarphi2(:, :, :) = 0.
        do oqn_l = 0, atoms%lmax(iDtype)
          do iradf = 1, nRadFun(oqn_l, iDtype)
            do imesh = 1, atoms%jri(iDtype)
              ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
              ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
              varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, iDtype) / atoms%rmsh(imesh, iDtype)
              varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, iDtype) / atoms%rmsh(imesh, iDtype)
            end do ! imesh
            ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
            call Derivative( varphi1(1:atoms%jri(iDtype), iradf, oqn_l), iDtype, atoms, delrVarphi1(1:atoms%jri(iDtype), iradf, oqn_l) )
            call Derivative( varphi2(1:atoms%jri(iDtype), iradf, oqn_l), iDtype, atoms, delrVarphi2(1:atoms%jri(iDtype), iradf, oqn_l) )
          end do ! iradf
        end do ! oqn_l

        call CalcChannelsGrFlpNat( atoms, iDtype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )

        write(*, *) 'here 1'
        ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = atoms%ntypsy(iDatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iDatom)
            mqn_m = mlh_atom(imem, ilh, iDatom)
            lm = lm_pre + mqn_m
            !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
            ! maybe construct a pointer and run only over them to make it memory efficient.
            do imesh = 1, atoms%jri(iDtype)
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + Veff0%vr(imesh, ilh, iDtype, 1) * clnu_atom(imem, ilh, iDatom)
            end do ! imesh
          end do ! imem
        end do ! ilh
        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, iDtype, iDatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, Veff0%vr(:, :, :, 1), clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi )
         write(*, *) 'here2'

        ! Call routine to calculate H grVarphi
        ! For testing reasons we want this cutoff because we project on the gradient of varphi
        lmaxBra = atoms%lmax(iDtype) + 1
         write(*, *) 'here3'
        do mqn_m2PrC = -1, 1
          call CalcHGrVarphi( atoms, iDtype, mqn_m2PrC, lmpMax, lmaxBra, grVarphiChMout, nRadFun, grVarPhiCh1, grVarPhiCh2,        &
                                                      & grVarphiChLout, vEff0NsphGrVarphi, El, lmpT, hGrVarphi )
        end do ! mqn_m2PrC

        call CalcVecBasfMatElems( atoms, iDtype, 2, nRadFun, r2, grVarphiChLout, grVarphiChMout, varphi1, varphi2, &
          & grVarPhiCh1, grVarphiCh2, hVarphi, grVarphiVarphi(:, :, :, iDtype), grVarphiHVarphi(:, :, :, iDatom) )


        do mqn_m2PrC = -1, 1
          do jj = 1, lmpMax
            do ii = 1, lmpMax
              varphiGrVarphi(ii, jj, mqn_m2PrC, iDtype) = grVarphiVarphi(jj, ii, mqn_m2PrC, iDtype)
            end do ! ii
          end do ! jj

          call CalcScalBasfMatElems( atoms, iDtype, iDatom, nRadFun, r2, varphi1, varphi2, hGrVarphi(:, :, :, :, mqn_m2PrC), varphiVarphiDummy, varphiHGrvarphi(:, :, :, mqn_m2PrC) )
        end do ! mqn_m2PrC

        r2(:) = 1.
        do idir = 1, 3
          call CalcScalBasfMatElems( atoms, iDtype, iDatom, nRadFun, r2, varphi1, varphi2, r2grVeff0SphVarphi(:, :, :, :, idir), varphiVarphiDummy, varphigrVeffvarphi(:, :, :, idir) )
        end do ! idir

        ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times.
        delrGrVarphiCh1(:, :, :, :) = 0.
        delrGrVarphiCh2(:, :, :, :) = 0.
        do mqn_m2PrR = -1, 1
          lmp = 0
          do oqn_l = 0, atoms%lmax(iDtype)
            do mqn_m = -oqn_l, oqn_l
              do iradf = 1, nRadFun(oqn_l, iDtype)
                lmp = lmp + 1
                do ichan = 1, 2
                  call Derivative( grVarphiCh1(:, ichan, lmp, mqn_m2PrR), iDtype, atoms, delrGrVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
                  call Derivative( grVarphiCh2(:, ichan, lmp, mqn_m2PrR), iDtype, atoms, delrGrVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
                end do ! ichan
              end do !iradf
            end do ! mqn_m
          end do ! oqn_l
        end do ! mqn_m2PrR

        ! Recycling <Psi2|H|Psi0> for the calculation of <Psi0|H|Psi2> requires additional terms correcting that the kinetic energy is not
        ! self-adjoint in LAPW.
        sAdjCorrVarphiHkinGrVarphi(:, :, :, :) = 0.
        do mqn_m2PrC = -1, 1
          call CalcSelfAdjCorrection( atoms, dimens, 2, nRadFunMax, iDtype, mqn_m2PrC, nRadFun, grVarphiChLout, grVarphiChMout,     &
            & varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiCh1, grVarphiCh2, delrGrVarphiCh1, delrGrVarphiCh2,              &
            sAdjCorrVarphiHkinGrVarphi(:, :, mqn_m2PrC, iDtype) )
        end do ! mqn_m2PrC

         write(*, *) 'here4'



      end do ! iDeqat
    end do ! iDtype

    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
      & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( grVarphiPsi(lmpMax), varphiGrPsi(lmpMax) )
    allocate( z1nG(dimens%nbasfcn, 3, atoms%nat, maxval(nobd(:, :))) )
    allocate( ab1cofVec( lmpMax, atoms%nat, 3) )
    allocate( ngoprI(atoms%nat) )
    ngoprI(:) = 1
    psiGrPsiNat(:, :) = cmplx(0., 0.)
    grPsiPsiNat(:, :) = cmplx(0., 0.)
    psiGrVeffPsi(:, :) = cmplx(0., 0.)
    psiHGrPsiNat(:, :) = cmplx(0., 0.)
    grPsiHPsiNat(:, :) = cmplx(0., 0.)
    sAdjCorrPsiHkinGrPsiNat(:, :) = cmplx(0., 0.)

    do ikpt = 1, kpts%nkpt
      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbas(1, mapGbas(:nv(1, ikpt), ikpt, 1)), &
        & gbas(2, mapGbas(:nv(1, ikpt), ikpt, 1)), gbas(3, mapGbas(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

      z1nG = cmplx(0.0, 0.0)
      kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
      gExt(:) = 0.
      do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
        gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBas(1:3, mapGbas(iBas, ikpt, 1)))
        do iband = 1, nobd(ikpt, 1)
          do idir = 1, 3
            z1nG(iBas, idir, 1, iband) = -iu * ( kExt(idir) + gExt(idir) ) * z(iBas, iband, ikpt, 1)
          end do ! iBas
        end do ! idir
      end do ! iband

      do iband = 1, nobd(ikpt, 1)
        iDatomA = 0
        do iDtypeA = 1, atoms%ntype
          do iDeqatA = 1, atoms%neq(iDtypeA)
            iDatomA = iDatomA + 1
            ab0cofScl(:) = cmplx(0., 0.)
            lmp = 0
            lm  = 0
            do oqn_l = 0, atoms%lmax(iDtypeA)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtypeA)
                ! p = 1
                ab0cofScl(lmp + 1) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatomA) )
                ! p = 2
                ab0cofScl(lmp + 2) = iu**oqn_l * dot_product( conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatomA) )
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) + iu**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) + iu**oqn_l * &
                    & dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                         & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                end do ! iradf


                ! This is a precalculation of the 1st and 3rd line of A.50. Actually the resulting lmp should have been primed or at
                ! least the p as lm1Pr = lm. But for sake of performance we place it here.
                ! sake of performance
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do !oqn_l

            ! Calculate the vector-like large matching coefficients using the first-order wave-function expansion coefficients gained
            ! from solving the Sternheimer equation.
            ! todo does not work for more than one atom
            ab1cofVec(:, :, :) = cmplx(0.0, 0.0)
            do idirR = 1, 3
              lm  = 0
              lmp = 0
              do oqn_l = 0, atoms%lmax(iDtypeA)
                do mqn_m = -oqn_l, oqn_l
                  pMaxLocal = nRadFun(oqn_l, iDtypeA)
                  ! p = 1
                  ab1cofVec(lmp + 1, iDatomA, idirR) = iu**oqn_l * &
                                 & dot_product( conjg( z1nG(:nv(1, ikpt), idirR, iDatomA, iband) ), a(:nv(1, ikpt), lm, iDatomA) )
                  ! p = 2
                  ab1cofVec(lmp + 2, iDatomA, idirR) = iu**oqn_l * &
                                 & dot_product( conjg( z1nG(:nv(1, ikpt), idirR, iDatomA, iband) ), b(:nv(1, ikpt), lm, iDatomA) )
                  do iradf = 3, pMaxLocal
                    ! p = 1
                    ab1cofVec(lmp + 1, iDatomA, idirR) = ab1cofVec(lmp + 1, iDatomA, idirR) &
                      & + iu**oqn_l * dot_product( conjg( z1nG(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR, iDatomA, iband) ), &
                                                  & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                    ! p = 2
                    ab1cofVec(lmp + 2, iDatomA, idirR) = ab1cofVec(lmp + 2, iDatomA, idirR) &
                      & + iu**oqn_l * dot_product(conjg(z1nG(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR, iDatomA, iband)), &
                                                  & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                    ! 2 < p < LOs for that l and that atom type
                    ab1cofVec(lmp + iradf, iDatomA, idirR) = &
                      & iu**oqn_l * dot_product( conjg( z1nG(nv(1, ikpt) + 1:nv(1, ikpt) + atoms%nlotot, idirR, iDatomA, iband) ),&
                                                  & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtypeA), iDatomA) )
                  end do ! iradf
                  lm = lm + 1
                  lmp = lmp + nRadFun(oqn_l, iDtypeA)
                end do ! mqn_m
              end do ! oqn_l
            end do ! idirR

            do idirC = 1, 3
              do idirR = 1, 3
                grVarphiPsi(:) = cmplx(0., 0.)
                grVarphiPsi(1:lmpT(iDtypeA)) = matmul(grVarphiVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), idirC - 2, iDtypeA), ab0cofScl(1:lmpT(iDtypeA)))
                grPsiPsiNat(idirR, idirC - 2) = grPsiPsiNat(idirR, idirC - 2) + dot_product(ab1cofVec(1:lmpT(iDtypeA), iDatomA, idirR), grVarphiPsi(1:lmpT(iDtypeA)) )
                varphiGrPsi(:) = cmplx(0., 0.)
                varphiGrPsi(1:lmpT(iDtypeA)) = matmul(varphiGrVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), idirC - 2, iDtypeA), ab1cofVec(1:lmpT(iDtypeA), iDatomA, idirR))
                psiGrPsiNat(idirR, idirC - 2) = psiGrPsiNat(idirR, idirC - 2) + dot_product(ab0cofScl(1:lmpT(iDtypeA)), varphiGrPsi(1:lmpT(iDtypeA)))

                grVarphiHPsi(:) = cmplx(0., 0.)
                grVarphiHPsi(1:lmpT(iDtypeA)) = matmul(grVarphiHVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), idirC - 2, iDatom), ab0cofScl(1:lmpT(iDtypeA)))
                grPsiHPsiNat(idirR, idirC - 2) = grPsiHPsiNat(idirR, idirC - 2) + dot_product(ab1cofVec(1:lmpT(iDatomA), iDatomA, idirR), grVarphiHPsi(1:lmpT(iDtypeA)))

                sAdjCorrVarphiHkinGrPsi(:) = cmplx(0., 0.)
                sAdjCorrVarphiHkinGrPsi(1:lmpT(iDtypeA)) = matmul(sAdjCorrVarphiHkinGrVarphi(1:lmpT(iDatomA), 1:lmpT(iDatomA), idirC - 2, iDtypeA), ab1cofVec(1:lmpT(iDatomA), iDatomA, idirR))
                sAdjCorrPsiHkinGrPsiNat(idirR, idirC - 2) = sAdjCorrPsiHkinGrPsiNat(idirR, idirC - 2) + dot_product(ab0cofScl(1:lmpT(iDtypeA)), sAdjCorrVarphiHkinGrPsi(1:lmpT(iDtypeA)))

                varphiHGrPsi(:) = cmplx(0., 0.)
                varphiHGrPsi(1:lmpT(iDtypeA)) = matmul(varphiHGrvarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDatomA, idirC - 2), ab1cofVec(1:lmpT(iDatomA), iDatomA, idirR))
                psiHGrPsiNat(idirR, idirC - 2) = psiHGrPsiNat(idirR, idirC - 2) + dot_product(ab0cofScl(1:lmpT(iDtypeA)), varphiHGrPsi(1:lmpT(iDtypeA)))

                varphiGrVeffPsi(:) = cmplx(0., 0.)
                varphiGrVeffPsi(1:lmpT(iDtypeA)) = matmul(varphiGrVeffVarphi(1:lmpT(iDtypeA), 1:lmpT(iDtypeA), iDatomA, idirC), ab1cofVec(1:lmpT(iDatomA), iDatomA, idirR))
                psiGrVeffPsi(idirR, idirC) = psiGrVeffPsi(idirR, idirC) + dot_product(ab0cofScl(1:lmpT(iDtypeA)), varphiGrVeffPsi(1:lmpT(iDtypeA)))

              end do ! idirR
            end do ! idirC
          end do ! iDeqatA
        end do ! iDtypeA
      end do ! iband
    end do ! ikpt

    psiGrPsi(1:3, 1:3) = matmul(psiGrPsiNat(1:3, -1:1), Tmatrix_transposed(1:3, 1:3))
    grPsiPsi(1:3, 1:3) = matmul(grPsiPsiNat(1:3, -1:1), conjg(Tmatrix_transposed(1:3, 1:3)))

    psiHGrPsi(1:3, 1:3) = matmul(psiHGrPsiNat(1:3, -1:1), Tmatrix_transposed(1:3, 1:3))

    grPsiHPsi(1:3, 1:3) = matmul(grPsiHPsiNat(1:3, -1:1), conjg(Tmatrix_transposed(1:3, 1:3)))

    sAdjCorrPsiHkinGrPsi(1:3, 1:3) = matmul(sAdjCorrPsiHkinGrPsiNat(1:3, -1:1), Tmatrix_transposed(1:3, 1:3))

    write(*, '(a)') 'grPsiPsi'
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiPsi(1, :))
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiPsi(2, :))
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiPsi(3, :))

    write(*, '(a)') 'psiGrPsi'
    write(*, '(3(2(es16.8,1x),3x))') psiGrPsi(1, :)
    write(*, '(3(2(es16.8,1x),3x))') psiGrPsi(2, :)
    write(*, '(3(2(es16.8,1x),3x))') psiGrPsi(3, :)

    write(*, '(a)') 'psiHGrPsiBench'
    write(*, '(3(2(es16.8,1x),3x))') psiHGrPsi(1, :)
    write(*, '(3(2(es16.8,1x),3x))') psiHGrPsi(2, :)
    write(*, '(3(2(es16.8,1x),3x))') psiHGrPsi(3, :)

    write(*, '(a)') 'psiHGrPsiBenchgrveff'
    write(*, '(3(2(es16.8,1x),3x))') psiGrVeffPsi(1, :)
    write(*, '(3(2(es16.8,1x),3x))') psiGrVeffPsi(2, :)
    write(*, '(3(2(es16.8,1x),3x))') psiGrVeffPsi(3, :)

    write(*, '(a)') 'psiHGrPsiBenchSum'
    write(*, '(3(2(es16.8,1x),3x))') psiHGrPsi(1, :) !- psiGrVeffPsi(1, :)
    write(*, '(3(2(es16.8,1x),3x))') psiHGrPsi(2, :) !- psiGrVeffPsi(2, :)
    write(*, '(3(2(es16.8,1x),3x))') psiHGrPsi(3, :) !- psiGrVeffPsi(3, :)

    write(*, '(a)') 'psiHGrPsiOrig'
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiHPsi(1, :))
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiHPsi(2, :))
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiHPsi(3, :))

    write(*, '(a)') 'psiHGrPsiOrigcorr'
    write(*, '(3(2(es16.8,1x),3x))') sAdjCorrPsiHkinGrPsi(1, :)
    write(*, '(3(2(es16.8,1x),3x))') sAdjCorrPsiHkinGrPsi(2, :)
    write(*, '(3(2(es16.8,1x),3x))') sAdjCorrPsiHkinGrPsi(3, :)

    write(*, '(a)') 'psiHGrPsiOrigSum'
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiHPsi(1, :)) !- sAdjCorrPsiHkinGrPsi(1, :)
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiHPsi(2, :)) !- sAdjCorrPsiHkinGrPsi(2, :)
    write(*, '(3(2(es16.8,1x),3x))') conjg(grPsiHPsi(3, :)) !- sAdjCorrPsiHkinGrPsi(3, :)

    ! The kinetic energy is not sufficient to correct but we have to correct the spherical Hamiltonian, especially the mixing terms are relevant. It gives only the same here if we concentrate everything to the potential, so add the spherical component of the potential to the non-spherical component and switch off everything with El, probably the El terms is still working but not the mixing terms.
    ! Call HgrVarphi

    ! Call ProjOnGrVarphi with varphi

    ! Generate gradvarphi H - eps varphi

    ! Generate correction term

    ! Compare
  end subroutine testGradBraCorrection

  subroutine TestDobGrad(atoms, rbas1, rbas2, nRadFun)

    use mod_juPhonUtils, only : Derivative
    use m_types, only : t_atoms
    use mod_juPhonUtils, only : CalcChannelsGrFlpNat, CalcChannelsGrGrtFlpNat
    use m_jpConstants, only : Tmatrix, Tmatrix_transposed

    implicit none

    ! Type parameters
    type(t_atoms),              intent(in)  :: atoms

    ! Array parameters
    real,                       intent(in)  :: rbas1(:, :, 0:, :)
    real,                       intent(in)  :: rbas2(:, :, 0:, :)
    integer,                    intent(in)  :: nRadFun(0:, :)

    ! Scalar variables
    integer                                 :: oqn_l
    integer                                 :: iradf
    integer                                 :: imesh
    integer                                 :: itype
    integer                                 :: lmpMax
    integer                                 :: mqn_m
    integer                                 :: nRadFunMax
    integer                                 :: mqn_m2Pr
    integer                                 :: mqn_m1Pr
    integer                                 :: lmp
    integer                                 :: ichan
    integer                                 :: oqn_l1Pr
    integer                                 :: lmp1Pr
    integer                                 :: oqn_l3Pr
    integer                                 :: mqn_m3Pr
    integer                                 :: iradf1Pr
    integer                                 :: mqn_m2PrC
    integer                                 :: mqn_m2PrR



    ! Array variables
    real,          allocatable              :: varphi1(:, :, :)
    real,          allocatable              :: varphi2(:, :, :)
    real,          allocatable              :: delrVarphi1(:, :, :)
    real,          allocatable              :: delrVarphi2(:, :, :)
    integer,       allocatable              :: grVarphiChLout(:, :)
    integer,       allocatable              :: grVarphiChMout(:, :)
    integer,       allocatable              :: grGrtVarphiChLout(:, :)
    integer,       allocatable              :: grGrtVarphiChMout(:, :, :)
    real,          allocatable              :: grVarphiCh1(:, :, :, :)
    real,          allocatable              :: grVarphiCh2(:, :, :, :)
    integer,       allocatable              :: lmpT(:)
    real,          allocatable              :: grVarphiCh1Res(:, :, :)
    real,          allocatable              :: grVarphiCh2Res(:, :, :)
    real,          allocatable              :: delGrVarphiCh1Res(:, :, :)
    real,          allocatable              :: delGrVarphiCh2Res(:, :, :)
    integer,       allocatable              :: testgrVarphiChLout(:, :)
    integer,       allocatable              :: testgrVarphiChMout(:, :)
    real,          allocatable              :: testgrVarphiCh1(:, :, :, :, :)
    real,          allocatable              :: testgrVarphiCh2(:, :, :, :, :)
    real,          allocatable              :: delrGrVarphiCh1(:, :, :, :)
    real,          allocatable              :: delrGrVarphiCh2(:, :, :, :)
    real,          allocatable              :: grGrtVarphiCh1(:, :, :, :, :)
    real,          allocatable              :: grGrtVarphiCh2(:, :, :, :, :)


    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1Res( atoms%jmtd, lmpMax, -1:1), grVarphiCh2Res( atoms%jmtd, lmpMax, -1:1) )
    allocate( delGrVarphiCh1Res(atoms%jmtd, lmpMax, -1:1), delGrVarphiCh2Res(atoms%jmtd, lmpMax, -1:1) )
    allocate( delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( grGrtVarphiChLout(4, 0:atoms%lmaxd), grGrtVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1, -1:1), &
            & grGrtVarphiCh1(atoms%jmtd, 4, lmpMax, -1:1, -1:1), grGrtVarphiCh2(atoms%jmtd, 4, lmpMax, -1:1, -1:1) )

    ! Fill only l = 0 channel, so that a double application of the gradient only shows contributions in l' = 0 und l' = 2
    itype = 1
    varphi1(:, :, :) = 0.
    varphi2(:, :, :) = 0.
    delrVarphi1(:, :, :) = 0.
    delrVarphi2(:, :, :) = 0.
    do oqn_l = 0, 0
      do iradf = 1, nRadFun(oqn_l, itype)
        do imesh = 1, atoms%jri(itype)
          ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
          ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
          varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype) / atoms%rmsh(imesh, itype)
          varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype) / atoms%rmsh(imesh, itype)
        end do ! imesh
        ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
        call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
        call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
      end do ! iradf
    end do ! oqn_l


    ! Calculate the application of the gradient and onto the MT basis functions (matching coefficients
    ! have no spatial dependence) and determing its scattering channels.
    grVarphiChLout(:, :) = 0
    grVarphiChMout(:, :) = 0
    grVarphiCh1(:, :, :, :) = 0.
    grVarphiCh2(:, :, :, :) = 0.
    call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )


    ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times.
    delrGrVarphiCh1(:, :, :, :) = 0.
    delrGrVarphiCh2(:, :, :, :) = 0.
    do mqn_m2PrR = -1, 1
      lmp = 0
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          do iradf = 1, nRadFun(oqn_l, itype)
            lmp = lmp + 1
            do ichan = 1, 2
              call Derivative( grVarphiCh1(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
              call Derivative( grVarphiCh2(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
            end do ! ichan
          end do !iradf
        end do ! mqn_m
      end do ! oqn_l
    end do ! mqn_m2PrR

    grGrtVarphiChLout(:, :) = 0
    grGrtVarphiChMout(:, :, :) = 0
    grGrtVarphiCh1(:, :, :, :, :) = 0.
    grGrtVarphiCh2(:, :, :, :, :) = 0.
    ! Not tested
    call CalcChannelsGrGrtFlpNat( atoms, itype, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, &
                        & delrGrVarphiCh1, delrGrVarphiCh2, grGrtVarphiChLout, grGrtVarphiChMout, grGrtVarphiCh1, grGrtVarphiCh2 )


    ! Put the gradient together. Does only work without LOs
    itype = 1
    grVarphiCh1Res(:, :, :) = 0.
    grVarphiCh2Res(:, :, :) = 0.
    do mqn_m2Pr = -1, 1
      lmp = 0
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          mqn_m1Pr = grVarphiChMout(mqn_m, mqn_m2Pr)
          do iradf = 1, nRadFun(oqn_l, itype)
            lmp = lmp + 1
            do ichan = 1, 2
              oqn_l1Pr = grVarphiChLout(ichan, oqn_l)
              if ( ( abs(mqn_m1Pr) > oqn_l1Pr ) .or. ( oqn_l1Pr < 0 ) .or. ( oqn_l1Pr > atoms%lmax(itype) ) ) cycle
              lmp1Pr = 0
              ! Approach oqn_l1Pr
              do oqn_l3Pr = 0, oqn_l1Pr - 1
                do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
                  do iradf1Pr = 1, nRadFun(oqn_l3Pr, itype)
                    lmp1Pr = lmp1Pr + 1
                  end do ! iradf1Pr
                end do ! mqn_m2Pr
              end do ! oqn_l2Pr
              ! Approach mqn_m1Pr
              do mqn_m3Pr = -oqn_l1Pr, mqn_m1Pr - 1
                do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                  lmp1Pr = lmp1Pr + 1
                end do ! iradf1Pr
              end do ! mqn_m2Pr
              ! Approach correct radial function
              do iradf1Pr = 1, iradf
                lmp1Pr = lmp1Pr + 1
              end do ! iradf1Pr
              do imesh = 1, atoms%jri(itype)
                grVarphiCh1Res(imesh, lmp1Pr, mqn_m2Pr) = grVarphiCh1Res(imesh, lmp1Pr, mqn_m2Pr) + grVarphiCh1(imesh, ichan, lmp, mqn_m2Pr)
                grVarphiCh2Res(imesh, lmp1Pr, mqn_m2Pr) = grVarphiCh2Res(imesh, lmp1Pr, mqn_m2Pr) + grVarphiCh2(imesh, ichan, lmp, mqn_m2Pr)
              end do ! imesh
            end do ! ichan
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l
    end do ! mqn_m2Pr



    ! Derive these gradients
    delGrVarphiCh1Res(:, :, :) = cmplx(0., 0.)
    delGrVarphiCh2Res(:, :, :) = cmplx(0., 0.)
    itype = 1
    do mqn_m2Pr = -1, 1
      lmp = 0
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          do iradf = 1, nRadFun(oqn_l, itype)
            lmp = lmp + 1
            call Derivative( grVarphiCh1Res(1:atoms%jri(itype), lmp, mqn_m2Pr), itype, atoms, delGrVarphiCh1Res(1:atoms%jri(itype), lmp, mqn_m2Pr) )
            call Derivative( grVarphiCh2Res(1:atoms%jri(itype), lmp, mqn_m2Pr), itype, atoms, delGrVarphiCh2Res(1:atoms%jri(itype), lmp, mqn_m2Pr) )
          end do ! iradf
        end do ! mqn_m
      end do ! oqn_l
    end do ! mqn_m2Pr

    allocate( testgrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1, -1:1), testgrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1, -1:1), &
            & testgrVarphiChLout(2, 0:atoms%lmaxd), testgrVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    do mqn_m2PrC = -1, 1
      call TestcalcChannelsGrFlpNat( atoms, itype, nRadFun, grVarphiCh1Res(:, :, mqn_m2PrC), grVarphiCh2Res(:, :, mqn_m2PrC),      &
        & delGrVarphiCh1Res(:, :, mqn_m2PrC), delGrVarphiCh2Res(:, :, mqn_m2PrC), testgrVarphiChLout, testgrVarphiChMout,          &
        & testgrVarphiCh1(:, :, :, :, mqn_m2PrC), testgrVarphiCh2(:, :, :, :, mqn_m2PrC ) )
    end do ! mqn_m2PrC


    do mqn_m2PrC = -1, 1
      do mqn_m2PrR = -1, 1
        do iradf = 1, 2
          do ichan = 1, 1
          do imesh = 1, atoms%jri(itype)
            write(2288, '(5i8,2f20.8)') mqn_m2PrR, mqn_m2PrC, iradf, ichan, imesh, grGrtVarphiCh1(imesh, ichan, iradf, mqn_m2PrR, mqn_m2PrC)
          end do ! imesh
          end do
        end do ! iradf
      end do ! mqn_m2PrR
    end do ! mqn_m2PrC

    do mqn_m2PrC = -1, 1
      do mqn_m2PrR = -1, 1
        do iradf = 1, 2
         do ichan = 1, 1
          do imesh = 1, atoms%jri(itype)
            write(2289, '(5i8,2f20.8)') mqn_m2PrR, mqn_m2PrC, iradf, ichan, imesh, testgrVarphiCh1(imesh, ichan, iradf, mqn_m2PrR, mqn_m2PrC)
          end do ! imesh
         end do ! ichan
        end do ! iradf
      end do ! mqn_m2PrR
    end do ! mqn_m2PrC

    ! works until here!


    if (.false.) then
      ! Fill only l = 0 channel, so that a double application of the gradient only shows contributions in l' = 0 und l' = 2
      itype = 1
      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      delrVarphi1(:, :, :) = 0.
      delrVarphi2(:, :, :) = 0.
      do oqn_l = 2, 2
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype) / atoms%rmsh(imesh, itype)
          end do ! imesh
          ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
          call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
          call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
        end do ! iradf
      end do ! oqn_l


      ! Calculate the application of the gradient and onto the MT basis functions (matching coefficients
      ! have no spatial dependence) and determing its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                          & grVarphiCh1, grVarphiCh2 )


      ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times.
      delrGrVarphiCh1(:, :, :, :) = 0.
      delrGrVarphiCh2(:, :, :, :) = 0.
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              do ichan = 1, 2
                call Derivative( grVarphiCh1(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
                call Derivative( grVarphiCh2(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
              end do ! ichan
            end do !iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR

      grGrtVarphiChLout(:, :) = 0
      grGrtVarphiChMout(:, :, :) = 0
      grGrtVarphiCh1(:, :, :, :, :) = 0.
      grGrtVarphiCh2(:, :, :, :, :) = 0.
      ! Not tested
      call CalcChannelsGrGrtFlpNat( atoms, itype, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, &
                          & delrGrVarphiCh1, delrGrVarphiCh2, grGrtVarphiChLout, grGrtVarphiChMout, grGrtVarphiCh1, grGrtVarphiCh2 )


      ! Put the gradient together. Does only work without LOs
      itype = 1
      grVarphiCh1Res(:, :, :) = 0.
      grVarphiCh2Res(:, :, :) = 0.
      do mqn_m2Pr = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            mqn_m1Pr = grVarphiChMout(mqn_m, mqn_m2Pr)
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              do ichan = 1, 2
                oqn_l1Pr = grVarphiChLout(ichan, oqn_l)
                if ( ( abs(mqn_m1Pr) > oqn_l1Pr ) .or. ( oqn_l1Pr < 0 ) .or. ( oqn_l1Pr > atoms%lmax(itype) ) ) cycle
                lmp1Pr = 0
                ! Approach oqn_l1Pr
                do oqn_l3Pr = 0, oqn_l1Pr - 1
                  do mqn_m3Pr = -oqn_l3Pr, oqn_l3Pr
                    do iradf1Pr = 1, nRadFun(oqn_l3Pr, itype)
                      lmp1Pr = lmp1Pr + 1
                    end do ! iradf1Pr
                  end do ! mqn_m2Pr
                end do ! oqn_l2Pr
                ! Approach mqn_m1Pr
                do mqn_m3Pr = -oqn_l1Pr, mqn_m1Pr - 1
                  do iradf1Pr = 1, nRadFun(oqn_l1Pr, itype)
                    lmp1Pr = lmp1Pr + 1
                  end do ! iradf1Pr
                end do ! mqn_m2Pr
                ! Approach correct radial function
                do iradf1Pr = 1, iradf
                  lmp1Pr = lmp1Pr + 1
                end do ! iradf1Pr
                do imesh = 1, atoms%jri(itype)
                  grVarphiCh1Res(imesh, lmp1Pr, mqn_m2Pr) = grVarphiCh1Res(imesh, lmp1Pr, mqn_m2Pr) + grVarphiCh1(imesh, ichan, lmp, mqn_m2Pr)
                  grVarphiCh2Res(imesh, lmp1Pr, mqn_m2Pr) = grVarphiCh2Res(imesh, lmp1Pr, mqn_m2Pr) + grVarphiCh2(imesh, ichan, lmp, mqn_m2Pr)
                end do ! imesh
              end do ! ichan
            end do ! iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2Pr



      ! Derive these gradients
      delGrVarphiCh1Res(:, :, :) = cmplx(0., 0.)
      delGrVarphiCh2Res(:, :, :) = cmplx(0., 0.)
      itype = 1
      do mqn_m2Pr = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              call Derivative( grVarphiCh1Res(1:atoms%jri(itype), lmp, mqn_m2Pr), itype, atoms, delGrVarphiCh1Res(1:atoms%jri(itype), lmp, mqn_m2Pr) )
              call Derivative( grVarphiCh2Res(1:atoms%jri(itype), lmp, mqn_m2Pr), itype, atoms, delGrVarphiCh2Res(1:atoms%jri(itype), lmp, mqn_m2Pr) )
            end do ! iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2Pr

      allocate( testgrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1, -1:1), testgrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1, -1:1), &
              & testgrVarphiChLout(2, 0:atoms%lmaxd), testgrVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
            !todo write this more general!
      do mqn_m2PrC = -1, 1
        call TestcalcChannelsGrFlpNat( atoms, itype, nRadFun, grVarphiCh1Res(:, :, mqn_m2PrC), grVarphiCh2Res(:, :, mqn_m2PrC),      &
          & delGrVarphiCh1Res(:, :, mqn_m2PrC), delGrVarphiCh2Res(:, :, mqn_m2PrC), testgrVarphiChLout, testgrVarphiChMout,          &
          & testgrVarphiCh1(:, :, :, :, mqn_m2PrC), testgrVarphiCh2(:, :, :, :, mqn_m2PrC ) )
      end do ! mqn_m2PrC


      do mqn_m2PrC = -1, 1
        do mqn_m2PrR = -1, 1
          do iradf = 1, 2
            do ichan = 1, 1
            do imesh = 1, 1!atoms%jri(itype)
              write(2288, '(5i8,2f20.8)') mqn_m2PrR, mqn_m2PrC, iradf, ichan, imesh, grGrtVarphiCh1(imesh, ichan, iradf, mqn_m2PrR, mqn_m2PrC)
            end do ! imesh
            end do
          end do ! iradf
        end do ! mqn_m2PrR
      end do ! mqn_m2PrC

      do mqn_m2PrC = -1, 1
        do mqn_m2PrR = -1, 1
          do iradf = 1, 2
           do ichan = 1, 1
            do imesh = 1, 1!atoms%jri(itype)
              write(2289, '(5i8,2f20.8)') mqn_m2PrR, mqn_m2PrC, iradf, ichan, imesh, testgrVarphiCh1(imesh, ichan, iradf, mqn_m2PrR, mqn_m2PrC)
            end do ! imesh
           end do ! ichan
          end do ! iradf
        end do ! mqn_m2PrR
      end do ! mqn_m2PrC
    end if

    ! Recalculate it with CalcChannelsGrFlpNat



  end subroutine TestDobGrad

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the spherical-harmonics expansion-coefficients coordinates of a general MT function's gradient resolved in its
  !> scattering channels within the natural coordinate system.
  !>
  !> @details
  !> The MT gradient of a function expanded spherical harmonics scatters into two neighbouring l channels and respective outgoing
  !> m-channels. This routine gives the resulting output l and m quantum numbers and the respective gradient channels. Furthermore,
  !> this routine remains within the natural coordinate system.
  !>
  !> @note
  !> The input function should be dependent on the orbital quantum number l and the radial function index p and can contain a large
  !> and a small relativistic component.
  !>
  !> @todo
  !> review documentation
  !> move whole routine two spaces to the left
  !--------------------------------------------------------------------------------------------------------------------------------------
  subroutine TestcalcChannelsGrFlpNat( atoms, itype, nRadFun, flp1, flp2, delrFlp1, delrFlp2, grFlpChLout, grFlpChMout, grFlpCh1,    &
      & grFlpCh2 )

    use m_JPConstants, only : fpi
    use m_types, only : t_atoms
    use m_gaunt, only : gaunt1

    implicit none

    ! Type parameter
    type(t_atoms),               intent(in)  :: atoms

    ! Scalar parameter
    integer,                     intent(in)  :: itype

    ! Array parameter
    ! nRadFun     : number of radial functions per l and type
    ! flp1        : 1st(large) component of function with l and p(p loops over the radial function types)
    ! flp2        : 2nd(small) component of function with l and p(p loops over the radial function types)
    ! delrFlp1    : Radial derivative of the 1st component of function with l and p(p loops over the radial function types)
    ! delrFlp2    : Radial derivative of the 2nd component of function with l and p(p loops over the radial function types)
    ! grFlpChLout : Contains the resulting outgoing l of a respective scattering channel of the gradient
    ! grFlpChMout : Contains the resulting outgoing m of a respective scattering channel of the gradient
    ! grFlpCh1    : Contains the 1st (large) component of the gradient scattering channels of the input function
    ! grFlpCh2    : Contains the 2nd (small) component of the gradient scattering channels of the input function
    integer,                     intent(in)  :: nRadFun(0:, :)
    real,                        intent(in)  :: flp1(:, : )
    real,                        intent(in)  :: flp2(:, :)
    real,                        intent(in)  :: delrFlp1(:, :)
    real,                        intent(in)  :: delrFlp2(:, :)
    integer,                     intent(out) :: grFlpChLout(:, 0:)
    integer,                     intent(out) :: grFlpChMout(-atoms%lmaxd:, -1:)
    real,                        intent(out) :: grFlpCh1(:, :, :, -1:)
    real,                        intent(out) :: grFlpCh2(:, :, :, -1:)

    ! Local Variables:
    !
    ! pfac     : contains sqrt(4 pi / 3)
    ! lm          : encodes oqn_l and mqn_m
    ! idirec      : runs over 3 directions the atom can be displaced to
    ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
    ! imesh       : runs over mesh points of current grid
    ! mqn_m       : magnetic quantum number m
    ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
    ! oqn_l       : orbital quantum number l

    ! Local Scalar Variables
    real                                     :: pfac
    integer                                  :: oqn_l
    integer                                  :: mqn_m
    integer                                  :: iradf
    integer                                  :: lmp
    real                                     :: gauntCoeff
    integer                                  :: imesh
    integer                                  :: mqn_m2Pr


    ! Local Array Variables

    grFlpCh1(:, :, :, :) = 0.
    grFlpCh2(:, :, :, :) = 0.


    pfac = sqrt( fpi / 3. )

    ! Precalculate L-output channels
    do oqn_l = 0, atoms%lmax(itype)
      grFlpChLout(1, oqn_l) = oqn_l + 1
      grFlpChLout(2, oqn_l) = oqn_l - 1
    end do ! oqn_l
    ! Set this lout = 0 - 1 to abnormal value as not allowed
    grFlpChLout(2, 0) = -9999

    ! Precalculate M-output channels
    do mqn_m2Pr = -1, 1
      do mqn_m = -atoms%lmax(itype), atoms%lmax(itype)
        grFlpChMout(mqn_m, mqn_m2Pr) = mqn_m - mqn_m2Pr
      end do ! mqn_m
    end do ! mqn_m2Pr

    do mqn_m2Pr = -1, 1
      lmp = 0
      do oqn_l = 0, atoms%lmax(itype)
        do mqn_m = -oqn_l, oqn_l
          do iradf = 1, nRadFun(oqn_l, itype)
            lmp = lmp + 1

            ! scattering channel (l + 1)
            if ( abs(mqn_m - mqn_m2Pr) <= oqn_l + 1 ) then
              gauntCoeff = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmax(itype) + 1 )
              do imesh = 1, atoms%jri(itype)
                ! Consider large and small relativistic components of radial solution
                grFlpCh1(imesh, 1, 0 + iradf, mqn_m2Pr) = grFlpCh1(imesh, 1, 0 + iradf, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                               & * (delrFlp1(imesh, lmp) -  oqn_l      * flp1(imesh, lmp) / atoms%rmsh(imesh, itype))
                grFlpCh2(imesh, 1, 0 + iradf, mqn_m2Pr) = grFlpCh2(imesh, 1, 0 + iradf, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                               & * (delrFlp2(imesh, lmp) -  oqn_l      * flp2(imesh, lmp) / atoms%rmsh(imesh, itype))
              enddo ! imesh
            endif ! scattering channel (l + 1)

            ! scattering channel (l - 1)
            ! This condition ensures that oqn_l = 0 is not accepted due to the emerging false condition 0 or 1 <= -1
            if ( ( abs(mqn_m - mqn_m2Pr) <= oqn_l - 1 ) ) then
              gauntCoeff = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmax(itype) + 1)
              do imesh = 1, atoms%jri(itype)
                ! Consider large and small relativistic components of radial solution
                grFlpCh1(imesh, 2, 0 + iradf, mqn_m2Pr) = grFlpCh1(imesh, 2, 0+ iradf, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                               & * (delrFlp1(imesh, lmp) + (oqn_l + 1) * flp1(imesh, lmp) / atoms%rmsh(imesh, itype))
                grFlpCh2(imesh, 2, 0 + iradf, mqn_m2Pr) = grFlpCh2(imesh, 2, 0 + iradf, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                               & * (delrFlp2(imesh, lmp) + (oqn_l + 1) * flp2(imesh, lmp) / atoms%rmsh(imesh, itype))
              enddo ! imesh
            endif ! scattering channel (l - 1)

          end do ! p
        end do ! mqn_m
      end do ! oqn_l
    end do ! mqn_m2Pr

  end subroutine TestcalcChannelsGrFlpNat

  ! Just gives out the matrix elmeent <grPsi|H - eps|grPsi>
  subroutine TestGrHepsGrt(atoms, dimens, lathar, kpts, cell, input, stars, Veff0, sym, usdus, results, ngdp, nRadFun, nobd, rbas1, rbas2, nv, ilst, GbasVec, z, gdp, kpq2kPrVec, El, eig, mlh_atom, vEff0MtLh, nmem_atom, clnu_atom, kveclo, iloTable )

    use m_types, only : t_atoms, t_dimension, t_kpts, t_sphhar, t_cell, t_input, t_stars, t_potential, t_sym, t_usdus, t_results
    use m_jpConstants, only : iu, Tmatrix, Tmatrix_transposed
    use mod_juPhonUtils, only : Derivative, CalcChannelsGrFlpNat, CalcChannelsGrGrtFlpNat
    use m_jpSternhHF, only : calcMEPotIR
    use m_jpSetupDynMat, only : CalcGrVarphiHepsGrtVarphiElem, CalcVecBasfMatElems, calcPsi1HepsPsi1IR
    use m_od_types, only : od_inp, od_sym
    use m_abcof3
    use m_jpSetupDynMatHelper, only : CalcHnGrV0Varphi

    implicit none

    type(t_atoms),                  intent(in) :: atoms
    type(t_dimension),              intent(in) :: dimens
    type(t_kpts),                   intent(in) :: kpts
    type(t_sphhar),                 intent(in) :: lathar
    type(t_cell),                   intent(in) :: cell
    type(t_input),                  intent(in) :: input
    type(t_stars),                  intent(in) :: stars
    type(t_potential),              intent(in) :: Veff0
    type(t_sym),                    intent(in) :: sym
    type(t_usdus),                  intent(in)    :: usdus
    type(t_results),                intent(in) :: results

    integer,                        intent(in) :: ngdp

    integer,                        intent(in) :: nRadFun(0:, :)
    integer,                        intent(in) :: nobd(:, :)
    real,                           intent(in) :: rbas1(:,:,0:,:,:)
    real,                           intent(in) :: rbas2(:,:,0:,:,:)
    integer,                        intent(in) :: nv(:, :)
    integer,                        intent(in) :: ilst(:, :, :)
    integer,                        intent(in) :: GbasVec(:, :)
    MCOMPLEX,                       intent(in) :: z(:, :, :, :)
    integer,                        intent(in) :: gdp(:, :)
    integer,                        intent(in) :: kpq2kPrVec(:, :, :)
    real,                           intent(in) :: El(:, 0:, :, :)
    real,                           intent(in) :: eig(:, :, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    real,                        intent(in) :: vEff0MtLh(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    complex,                        intent(in) :: clnu_atom(:,0:,:)
    integer,                        intent(in) :: kveclo(:,:)
    integer,                        intent(in)    :: iloTable(:, 0:, :)

    integer                                    :: iqpt
    integer                                    :: itype
    integer                                    :: lmpMax
    integer                                    :: nRadFunMax
    integer                                    :: oqn_l
    integer                                    :: iatom
    integer                                    :: iradf
    integer                                    :: imesh
    integer                                    :: ieqat
    integer                                    :: ikpt
    integer                                    :: idirR
    integer                                    :: iBas
    integer                                    :: iband
    integer                                    :: nmat
    integer                                    :: idirC
    integer                                    :: ispin
    integer                                    :: mqn_m2prC
    integer                                    :: mqn_m2prR
    integer                                    :: lmp
    integer                                    :: mqn_m
    integer                                    :: ichan
    integer                                    :: ptsym
    integer                                    :: ilh
    integer                                    :: lm_pre
    integer                                    :: imem
    integer                                    :: lm
    integer                                    :: iDatom
    integer                                    :: iDtype
    integer                                    :: iDeqat
    integer                                    :: pMaxLocal

    ! Type variable
    type(od_inp)                                  :: odi
    type(od_sym)                                  :: ods

    integer,           allocatable             :: ngoprI(:)
    complex,           allocatable             :: a(:, :, :)
    complex,           allocatable             :: b(:, :, :)
    complex,           allocatable             :: bascof_lo(:, :, :, :, :)
    integer,           allocatable             :: lmpT(:)
    real,              allocatable             :: varphi1(:, :, :)
    real,              allocatable             :: varphi2(:, :, :)
    real,              allocatable             :: grVarphiCh1(:, :, :, :)
    real,              allocatable             :: grVarphiCh2(:, :, :, :)
    integer,           allocatable             :: grVarphiChLout(:, :)
    integer,           allocatable             :: grVarphiChMout(:, :)
    real,              allocatable             :: delrVarphi1(:, :, :)
    real,              allocatable             :: delrVarphi2(:, :, :)
    complex,           allocatable             :: ikpGz0(:, :, :, :)
    complex,           allocatable             :: vEff0MtSpH(:, :)
    complex,           allocatable             :: vEff0NsphGrVarphi(:, :, :, :, :)
    complex,           allocatable             :: r2grVeff0SphVarphi(:, :, :, :, :)
    complex,           allocatable             :: grVarphiGrtVeff0SphVarphiNat(:, :, :, :, :)
    real,              allocatable             :: vEff0MtLhDummy(:, :, :)
    real,              allocatable             :: r2(:)
    real,              allocatable             :: grVarphiGrtVarphiNatNat(:, :, :, :, :)
    complex,           allocatable             :: grVarphiHpreGrtVarphiNatNat(:, :, :, :, :)
    complex,           allocatable             :: altIntgrlIRband(:, :, :, :)
    complex,           allocatable             :: w_vEff1IRContainer(:, :, :)
    real,          allocatable              :: delrGrVarphiCh1(:, :, :, :)
    real,          allocatable              :: delrGrVarphiCh2(:, :, :, :)
    integer,       allocatable              :: grGrtVarphiChLout(:, :)
    integer,       allocatable              :: grGrtVarphiChMout(:, :, :)
    real,          allocatable              :: grGrtVarphiCh1(:, :, :, :, :)
    real,          allocatable              :: grGrtVarphiCh2(:, :, :, :, :)
    complex,           allocatable             :: hVarphi(:, :, :, :)
    real,              allocatable             :: grGrtVarphiVarphiNatNat(:, :, :, :, :)
    complex,           allocatable             :: grGrtVarphiHvarphiNatNat(:, :, :, :, :)
    complex,        allocatable              :: kpGKpGTz0(:, :, :, :, :)
    complex,           allocatable                :: grVarphiGrtVeff0SphPsiNat(:)
    complex,           allocatable                :: grVarphiHepsGrtPsiNatNat(:)
    complex, allocatable :: ab0cofScl(:)

    real                                       :: kExt(3)
    real                                       :: Gext(3)
    complex                                    :: altIntgrl(3, 3)
    complex                                    :: altIntgrlIR(3, 3)
    complex                                    :: altIntgrlMT(3, 3)
    complex                                    :: grGrtPsiHepsPsiBand
    complex                                    :: grGrtPsiHepsPsi(3, 3)
    complex                                       :: grPsiHepsGrtPsiNatNat(3, 3)
    complex                                       :: grPsiGrtVeff0SphPsiNat(3, 3)
    complex                                       :: grPsiHepsGrtPsiNat(3, 3)
    complex                                       :: grPsiHepsGrtPsi(3, 3)

    iqpt = 1

    ! Quantities for initialization
    allocate( lmpT(atoms%ntype) )
    lmpT(:) = 0
    do itype = 1, atoms%ntype
      lmpT(itype) = sum( [ ( (2 * oqn_l + 1)* nRadFun(oqn_l, itype), oqn_l = 0, atoms%lmax(itype) ) ] )
    end do ! itype
    lmpMax     = maxval( lmpT(:) )
    nRadFunMax = maxval( nRadFun(:, :) )

    allocate( varphi1(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), varphi2(atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( grVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), grVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1), &
            & grVarphiChLout(2, 0:atoms%lmaxd), grVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1) )
    allocate( delrVarphi1( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd), delrVarphi2( atoms%jmtd, nRadFunMax, 0:atoms%lmaxd) )
    allocate( ikpGz0(dimens%nbasfcn, dimens%neigd, 3, kpts%nkpt) )
    allocate( vEff0MtSpH( atoms%jmtd, 0:dimens%lmd), vEff0NsphGrVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, -1:1), &
            & r2grVeff0SphVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax, 3) )
    allocate( grVarphiGrtVeff0SphVarphiNat(lmpMax, lmpMax, -1:1, 1:3, atoms%nat), vEff0MtLhDummy( atoms%jmtd, 0:lathar%nlhd, atoms%ntype ) )
    allocate( r2(atoms%jmtd) )
    allocate( hVarphi(2, atoms%jmtd, 0:(atoms%lmaxd + 3)**2 - 1, lmpMax)  )
    allocate( grVarphiGrtVarphiNatNat(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grVarphiHpreGrtVarphiNatNat(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( altIntgrlIRband(maxval(nobd), maxval(nobd), 3, 3) )
    allocate( delrGrVarphiCh1(atoms%jmtd, 2, lmpMax, -1:1), delrGrVarphiCh2(atoms%jmtd, 2, lmpMax, -1:1) )
    allocate( grGrtVarphiChLout(4, 0:atoms%lmaxd), grGrtVarphiChMout(-atoms%lmaxd:atoms%lmaxd, -1:1, -1:1), &
            & grGrtVarphiCh1(atoms%jmtd, 4, lmpMax, -1:1, -1:1), grGrtVarphiCh2(atoms%jmtd, 4, lmpMax, -1:1, -1:1) )
    allocate( grGrtVarphiVarphiNatNat(lmpMax, lmpMax, -1:1, -1:1, atoms%ntype), &
            & grGrtVarphiHvarphiNatNat(lmpMax, lmpMax, -1:1, -1:1, atoms%nat) )
    allocate( a( dimens%nvd, 0:dimens%lmd, atoms%nat), b(dimens%nvd, 0:dimens%lmd, atoms%nat), &
            & bascof_lo(3, -atoms%llod:atoms%llod, 4 * atoms%llod + 2, atoms%nlod, atoms%nat) )
    allocate( grVarphiGrtVeff0SphPsiNat(lmpMax))
    vEff0MtSpH = cmplx(0., 0.)
    grVarphiGrtVarphiNatNat(:, :, :, :, :) = 0.
    grVarphiHpreGrtVarphiNatNat(:, :, :, :, :) = 0.
    grGrtVarphiVarphiNatNat = cmplx(0., 0.)
    grGrtVarphiHVarphiNatNat = cmplx(0., 0.)
    grVarphiGrtVeff0SphVarphiNat = cmplx(0., 0.)
    !todo we have to reset stuff

    iatom = 0
    do itype = 1, atoms%ntype

      varphi1(:, :, :) = 0.
      varphi2(:, :, :) = 0.
      delrVarphi1(:, :, :) = 0.
      delrVarphi2(:, :, :) = 0.
      do oqn_l = 0, atoms%lmax(itype)
        do iradf = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            ! In Fleur the radial solutions u_lp are multiplied by a factor r by default to avoid an additional multiplication of the
            ! Jacobi determinant r^2 in radial integrals given spherical coordinates.
            varphi1(imesh, iradf, oqn_l) = rbas1(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
            varphi2(imesh, iradf, oqn_l) = rbas2(imesh, iradf, oqn_l, itype, 1) / atoms%rmsh(imesh, itype)
          end do ! imesh
          ! Precalculate partial derivatives of varphis in r-direction since it is needed twice
          call Derivative( varphi1(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi1(1:atoms%jri(itype), iradf, oqn_l) )
          call Derivative( varphi2(1:atoms%jri(itype), iradf, oqn_l), itype, atoms, delrVarphi2(1:atoms%jri(itype), iradf, oqn_l) )
        end do ! iradf
      end do ! oqn_l

      ! Calculate the application of the gradient and the gradient's dyadic product onto the MT basis functions (matching coefficients
      ! have no spatial dependence) and determing its scattering channels.
      grVarphiChLout(:, :) = 0
      grVarphiChMout(:, :) = 0
      grVarphiCh1(:, :, :, :) = 0.
      grVarphiCh2(:, :, :, :) = 0.
      call CalcChannelsGrFlpNat( atoms, itype, nRadFun, varphi1, varphi2, delrVarphi1, delrVarphi2, grVarphiChLout, grVarphiChMout, &
                                                                                                        & grVarphiCh1, grVarphiCh2 )
      ! Precalculate partial derivatives of the varphis' gradients in r-direction since it is required two times.
      delrGrVarphiCh1(:, :, :, :) = 0.
      delrGrVarphiCh2(:, :, :, :) = 0.
      do mqn_m2PrR = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1
              do ichan = 1, 2
                call Derivative( grVarphiCh1(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh1(:, ichan, lmp, mqn_m2PrR) )
                call Derivative( grVarphiCh2(:, ichan, lmp, mqn_m2PrR), itype, atoms, delrGrVarphiCh2(:, ichan, lmp, mqn_m2PrR) )
              end do ! ichan
            end do !iradf
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2PrR

      grGrtVarphiChLout(:, :) = 0
      grGrtVarphiChMout(:, :, :) = 0
      grGrtVarphiCh1(:, :, :, :, :) = 0.
      grGrtVarphiCh2(:, :, :, :, :) = 0.
      call CalcChannelsGrGrtFlpNat( atoms, itype, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, grVarphiCh1, grVarphiCh2, &
                          & delrGrVarphiCh1, delrGrVarphiCh2, grGrtVarphiChLout, grGrtVarphiChMout, grGrtVarphiCh1, grGrtVarphiCh2 )

      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1

        r2grVeff0SphVarphi = cmplx(0., 0.)
        vEff0NsphGrVarphi = cmplx(0., 0.)
        ! Expand the coefficients of the lattice-harmonic given potential into spherical-harmonic coefficients for the given atom.
        vEff0MtSpH(:, :) = cmplx(0.0, 0.0)
        ptsym = atoms%ntypsy(iatom)
        do ilh = 0, lathar%nlh(ptsym)
          oqn_l = lathar%llh(ilh, ptsym)
          lm_pre = oqn_l * (oqn_l + 1)
          do imem = 1, nmem_atom(ilh, iatom)
            mqn_m = mlh_atom(imem, ilh, iatom)
            lm = lm_pre + mqn_m
            !todo one could only evaluate the vEff0MtSpH which have a contribution, i.e. the oqn_l and mqn_m which are in llh and mlh_atom
            ! maybe construct a pointer and run only over them to make it memory efficient.
            do imesh = 1, atoms%jri(itype)
              vEff0MtSpH(imesh, lm) = vEff0MtSpH(imesh, lm) + vEff0MtLh(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom)
            end do ! imesh
          end do ! imem
        end do ! ilh

        hVarphi = cmplx(0.0, 0.0)
        call CalcHnGrV0Varphi( atoms, lathar, itype, iatom, lmpMax, El, varphi1, varphi2, nRadFun, vEff0MtSpH, vEff0MtLh, clnu_atom, &
          & nmem_atom, mlh_atom, grVarphiCh1, grVarphiCh2, grVarphiChLout, grVarphiChMout, hVarphi, vEff0NsphGrVarphi, r2grVeff0SphVarphi )


        do mqn_m2PrC = -1, 1
          call CalcVecBasfMatElems( atoms, itype, iatom, nRadFun, r2, grGrtVarphiChLout, grGrtVarphiChMout(:, :, mqn_m2PrC),  &
            & varphi1, varphi2, grGrtVarPhiCh1(:, :, :, :, mqn_m2PrC), grGrtVarPhiCh2(:, :, :, :, mqn_m2PrC), hVarphi,        &
            & grGrtVarphiVarphiNatNat(:, :, :, mqn_m2PrC, itype), grGrtVarphiHVarphiNatNat(:, :, :, mqn_m2PrC, iatom) )

          call CalcGrVarphiHepsGrtVarphiElem( atoms, itype, iatom, mqn_m2PrC, lmpMax, nRadFun, grVarphiChLout, grVarphiChMout, r2, &
            & grVarPhiCh1, grVarPhiCh2, El, lmpT, vEff0NsphGrVarphi, r2grVeff0SphVarphi, grVarphiGrtVarphiNatNat, grVarphiHpreGrtVarphiNatNat, &
            & grVarphiGrtVeff0SphVarphiNat)
        end do ! mqn_m2PrC
      end do ! ieqat
    end do ! itype
!    do itype = 1, atoms%ntype
!      do mqn_m2PrR  = -1, 1
!        do mqn_m2PrC = -1, 1
!          do ii = 1, lmpMax
!            do jj = 1, lmpMax
!              write(2295, '(4i8,f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiGrtVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!              write(2296, '(4i8,2f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiHpreGrtVarphi(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!              write(2297, '(4i8,2f20.8)') ii, jj, mqn_m2PrR, mqn_m2PrC, grVarphiGrtVeff0SphVarphiNat(ii, jj, mqn_m2PrR, mqn_m2PrC, 1)
!            end do ! jj
!          end do ! ii
!        end do ! mqn_m2PrC
!      end do ! mqn_m2PrR
!    end do ! itype
    grGrtVarphiVarphiNatNat = cmplx(0., 0.)
    grGrtVarphiHVarphiNatNat = cmplx(0., 0.)

    allocate(ngoprI(atoms%nat))
    ngoprI(:) = 1
    allocate( ab0cofScl(lmpMax) )
    allocate( grVarphiHepsGrtPsiNatNat(lmpMax) )
    grPsiHepsGrtPsi(:, :) = cmplx(0., 0.)
    do ikpt = 1, kpts%nkpt
      nmat = nv(1, ikpt) + atoms%nlotot
      a(:, :, :) = cmplx(0.0, 0.0)
      b(:, :, :) = cmplx(0.0, 0.0)
      bascof_lo(:, :, :, :, :) = cmplx(0.0, 0.0)
      call abcof3( atoms%lmaxd, atoms%ntype, atoms%nat, sym%nop, dimens%nvd, dimens%jspd, 1, dimens%lmd, dimens%nbasfcn, &
        & atoms%llod, atoms%nlod, atoms%nlotot, sym%invtab, atoms%ntype, sym%mrot, ngoprI, atoms%taual, atoms%neq, atoms%lmax, &
        & atoms%rmt, cell%omtil, cell%bmat, cell%bbmat, kpts%bk(:, ikpt), gbasVec(1, ilst(:nv(1, ikpt), ikpt, 1)), &
        & gbasVec(2, ilst(:nv(1, ikpt), ikpt, 1)), gbasVec(3, ilst(:nv(1, ikpt), ikpt, 1)), nv(:, ikpt), nmat, &
        & usdus%us, usdus%dus, usdus%uds, usdus%duds, usdus%ddn, atoms%invsat, sym%invsatnr, usdus%ulos, usdus%uulon, usdus%dulon, &
        & usdus%dulos, atoms%llo, atoms%nlo, atoms%l_dulo, atoms%lapw_l, kveclo(:,ikpt), odi, ods, a, b, bascof_lo )

      do iband = 1, nobd(ikpt, 1)
        iDatom = 0
        do iDtype = 1, atoms%ntype
          do iDeqat = 1, atoms%neq(iDtype)
            iDatom = iDatom + 1
            lm  = 0
            lmp = 0
            ab0cofScl(:) = cmplx(0., 0.)
            do oqn_l = 0, atoms%lmax(iDtype)
              do mqn_m = - oqn_l, oqn_l
                pMaxLocal = nRadFun(oqn_l, iDtype)
                ! p = 1
                ab0cofScl(lmp + 1) = iu**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), a(:nv(1, ikpt), lm, iDatom))
                ! p = 2
                ab0cofScl(lmp + 2) = iu**oqn_l * dot_product(conjg(z(:nv(1, ikpt), iband, ikpt, 1)), b(:nv(1, ikpt), lm, iDatom))
                ! Add LO contributions
                do iradf = 3, pMaxLocal
                  ! p = 1
                  ab0cofScl(lmp + 1) = ab0cofScl(lmp + 1) &
                    & + iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(1, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! p = 2
                  ab0cofScl(lmp + 2) = ab0cofScl(lmp + 2) &
                    & + iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(2, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                  ! 2 < p < LOs for that l and that atom type
                  ab0cofScl(lmp + iradf) = iu**oqn_l * dot_product( conjg(z(nv(1, ikpt) + 1 :nv(1, ikpt) + atoms%nlotot, iband, ikpt, 1)), &
                                                        & bascof_lo(3, mqn_m, :atoms%nlotot, iloTable(iradf, oqn_l, iDtype), iDatom) )
                end do ! iradf
                lm = lm + 1
                lmp = lmp + pMaxLocal
              end do ! mqn_m
            end do ! oqn_l
            grPsiHepsGrtPsiNatNat(:, :) = cmplx(0., 0.)
            grVarphiHepsGrtPsiNatNat = cmplx(0., 0.)
            grVarphiGrtVeff0SphPsiNat(:) = cmplx(0., 0.)
            grPsiGrtVeff0SphPsiNat(:, :) = cmplx(0., 0.)
            do idirC = -1, 1
              do idirR = -1, 1
                grVarphiHepsGrtPsiNatNat(1:lmpT(iDtype)) = &
                  &   matmul(grVarphiHpreGrtVarphiNatNat(1:lmpT(iDtype), 1:lmpT(iDtype), idirR, idirC, iDatom),                 &
                                                                                            & ab0cofScl(1:lmpT(iDtype)) )&
                  & - eig(iband, ikpt, 1) * matmul(grVarphiGrtVarphiNatNat(1:lmpT(iDtype), 1:lmpT(iDtype), idirR, idirC, iDtype)&
                                                                                          &  , ab0cofScl(1:lmpT(iDtype)) )

                grPsiHepsGrtPsiNatNat(idirR + 2, idirC + 2) = dot_product(ab0cofScl(1:lmpT(iDtype)),                              &
                                                                                       & grVarphiHepsGrtPsiNatNat(1:lmpT(iDtype)) )

                grVarphiGrtVeff0SphPsiNat(1:lmpT(iDtype)) = &
                  & matmul( grVarphiGrtVeff0SphVarphiNat(1:lmpT(iDtype), 1:lmpT(iDtype), idirR, idirC + 2, iDatom), &
                                                                                             & ab0cofScl(1:lmpT(iDtype)) )

                grPsiGrtVeff0SphPsiNat(idirR + 2, idirC + 2) = dot_product(ab0cofScl(1:lmpT(iDtype)), &
                                                                                      & grVarphiGrtVeff0SphPsiNat(1:lmpT(iDtype)) )

              end do ! idirR
            end do ! idirC
            grPsiHepsGrtPsiNat(1:3, 1:3) = cmplx(0., 0.)
            grPsiHepsGrtPsiNat(1:3, 1:3) = matmul(grPsiHepsGrtPsiNatNat(1:3, 1:3), Tmatrix_transposed(1:3, 1:3))
            grPsiHepsGrtPsiNat(1:3, 1:3) = grPsiHepsGrtPsiNat(1:3, 1:3) + grPsiGrtVeff0SphPsiNat(1:3, 1:3)
!            grPsiHepsGrtPsiNat(1:3, 1:3) = grPsiGrtVeff0SphPsiNat(1:3, 1:3)
            grPsiHepsGrtPsi(1:3, 1:3) = grPsiHepsGrtPsi(1:3, 1:3) + results%w_iks(iband, ikpt, 1) * matmul(conjg(Tmatrix(1:3, 1:3)), grPsiHepsGrtPsiNat(1:3, 1:3))
          end do ! iDeqat
        end do ! iDtype
!        write(*, *) iband, ikpt
      end do ! iband
    end do ! ikpt

    write(*, '(a)') 'Matrix of brakets'
    write(*, '(3(2(es16.8,1x),3x))') 4 * grPsiHepsGrtPsi(1, :)
    write(*, '(3(2(es16.8,1x),3x))') 4 * grPsiHepsGrtPsi(2, :)
    write(*, '(3(2(es16.8,1x),3x))') 4 * grPsiHepsGrtPsi(3, :)

    write(*, *) 'foo'

!    allocate( kpGKpGTz0(dimens%nvd, maxval(nobd), 3, 3, kpts%nkpt))
!    allocate( w_vEff1IRContainer(ngdp, 3, atoms%nat) )
!    w_vEff1IRContainer = cmplx(0., 0.)
!    iatom = 0
!    do itype = 1, atoms%ntype
!      do ieqat = 1, atoms%neq(itype)
!        iatom = iatom + 1
!        grGrtPsiHepsPsi(:, :) = cmplx(0., 0.)
!        do ikpt = 1, kpts%nkpt
!          ikpGz0 = cmplx(0.0, 0.0)
!          kExt(1:3) = matmul(cell%bmat, kpts%bk(1:3, ikpt))
!          gExt(:) = 0.
!          do idirR = 1, 3
!            do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
!              gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, ilst(iBas, ikpt, 1)))
!              do iband = 1, nobd(ikpt, 1)
!                ikpGz0(iBas, iband, idirR, ikpt) = iu * ( kExt(idirR) + gExt(idirR) ) * z(iBas, iband, ikpt, 1)
!              end do ! iband
!            end do ! iBas
!          end do ! idirR
!          do idirC = 1, 3
!            do idirR = 1, 3
!              do iBas = 1, nv(1, ikpt) !+ atoms%nlotot
!                gExt(1:3) = matmul( cell%bmat(1:3, 1:3), gBasVec(1:3, ilst(iBas, ikpt, 1)))
!                do iband = 1, nobd(ikpt, 1)
!                  kpGKpGTz0(iBas, iband, idirR, idirC, ikpt) = ( kExt(idirR) + gExt(idirR) ) * ( kExt(idirC) + gExt(idirC) ) &
!                                                                                                         & * z(iBas, iband, ikpt, 1)
!                end do ! iband
!              end do ! iBas
!            end do ! idirR
!          end do ! idirC
!
!          ispin = 1
!          nmat = nv(1, ikpt) + atoms%nlotot
!          altIntgrlIRband(:, :, :, :) = cmplx(0., 0.)
!          do idirC = 1, 3
!            do idirR = 1, 3
!              call calcMEPotIR( input, stars, dimens, GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)),                   &
!                & GbasVec(:, ilst(:nv(input%jspins, ikpt), ikpt, input%jspins)), nv, w_vEff1IRContainer(:, idirC, 1),                  &
!                & ikpGz0(:, :, idirR, ikpt), z(:, :, ikpt, 1), gdp, nmat, nobd(ikpt, 1), nobd(ikpt, 1), ispin, ikpt, iqpt, ikpt, ngdp, &
!                & altIntgrlIRband(:, :, idirR, idirC), kpq2kPrVec )
!              do iband = 1, nobd(ikpt, 1)
!                altIntgrlIR(idirR, idirC) = altIntgrlIR(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * altIntgrlIRband(iband, iband, idirR, idirC)
!                grGrtPsiHepsPsiBand = cmplx(0., 0.)
!                call calcPsi1HepsPsi1IR( kpts, stars, dimens, cell, gBasVec(:, ilst(:nv(1, ikpt), ikpt, 1)), nv, &
!                  & kpGKpGTz0(:, iband, idirR, idirC, ikpt), z(:, idirC, iatom, iband), nmat, ikpt, iqpt, ikpt, kpq2kPrVec, &
!                  & Veff0%vpw, eig, iband, grGrtPsiHepsPsiBand )
!                  grGrtPsiHepsPsi(idirR, idirC) = grGrtPsiHepsPsi(idirR, idirC) + 2 * results%w_iks(iband, ikpt, 1) * grGrtPsiHepsPsiBand
!              end do ! iband
!            end do ! idirR
!          end do ! idirC
!        end do !ikpt
!      end do ! ieqat
!    end do ! itype

  end subroutine TestGrHepsGrt

  ! Calculates the surface integral of rho0MT and the grVeff, while grVeff is derived nuemrically and not exactly the same on the
  ! boundary than calculated with the Weinert method
  subroutine Calc2ArgSFIntegrals( atoms, stars, cell, lathar, Veff0, ngdp, rho0IRst, gdp, clnu_atom, nmem_atom, mlh_atom, rho0MT )

    use m_types
    use mod_juPhonUtils, only : ConvertStar2G
    use m_jpConstants, only : iu
    use m_jpPotDensHelper, only : calcGrR2FinLH, WarpIRPot
    use m_jpSetupDynMatSF, only : CalcSurfIntIRDynMat, CalcSurfIntMTDynMat

    implicit none

    ! Type parameter
    type(t_atoms),                  intent(in) :: atoms
    type(t_stars),                  intent(in) :: stars
    type(t_cell),                   intent(in) :: cell
    type(t_sphhar),                 intent(in) :: lathar
    type(t_potential),              intent(in) :: Veff0

    ! Scalar parameter
    integer,                        intent(in) :: ngdp

    ! Array parameter
    complex,                        intent(in) :: rho0IRst(:, :)
    integer,                        intent(in) :: gdp(:, :)
    complex,                        intent(in) :: clnu_atom(:, 0:, :)
    integer,                        intent(in) :: nmem_atom(0:, :)
    integer,                        intent(in) :: mlh_atom(:, 0:, :)
    real,                           intent(in) :: rho0MT(:, 0:, :, :)

    ! Scalar variables
    integer                                    :: iG
    integer                                    :: itype
    integer                                    :: ilh
    integer                                    :: imesh
    integer                                    :: idir
    integer                                    :: iatom
    integer                                    :: oqn_l
    integer                                    :: lm_pre
    integer                                    :: mqn_m
    integer                                    :: lm
    integer                                    :: ieqat

    ! Array variables
    complex,           allocatable             :: vEff0IRpw(:)
    complex,           allocatable             :: rho0IRpw(:)
    complex,           allocatable             :: grVeff0IR(:, :)
    real,              allocatable             :: r2Veff0MT(:, :, :)
    complex,           allocatable             :: grVeff0MT(:, :, :, :)
    complex,           allocatable             :: r2GrVeff0MT(:, :, :, :)
    complex,           allocatable             :: w_grVeff0IR(:, :)
    complex                                    :: surfIntIR(3, 3)
    complex                                    :: surfIntMT(3, 3)
    real                                       :: Gext(3)

    allocate( vEff0IRpw(ngdp), rho0IRpw(ngdp) )
    vEff0IRpw(:) = cmplx(0., 0.)
    call ConvertStar2G( Veff0%vpw_uw(:, 1), vEff0IRpw(:), stars, ngdp, gdp )
    call ConvertStar2G( rho0IRst(:, 1), rho0IRpw(:), stars, ngdp, gdp )

    ! Perform analytical gradient of effective potential in the interstitial region
    allocate( grVeff0IR(ngdp, 3) )
    grVeff0IR(:, :) = cmplx(0., 0.)
    do iG = 1, ngdp
      Gext(1:3) = matmul(cell%bmat(1:3, 1:3), gdp(1:3, iG))
      grVeff0IR(iG, 1:3) = iu * Gext(1:3) * vEff0IRpw(iG)
    end do ! iG

    ! Calculate the MT gradient.
    !   Any factors as in Fleur have been removed. In order to improve the numerical accuracy of the gradient r^2 is multiplied
    !   before and removed after the calculation of the gradient. Thereyby, subtraction of small numbers close to the core are
    !   avoided.
    allocate( r2Veff0MT( atoms%jmtd, 0:lathar%nlhd, atoms%ntype) )
    allocate( grVeff0MT( atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    r2Veff0MT(:, :, :) = 0.
    grVeff0MT(:, :, :, :) = cmplx(0., 0.)

    do itype = 1, atoms%ntype
      do ilh = 0, lathar%nlhd
        do imesh = 1, atoms%jri(itype)
          ! SPIN MISSING
          r2Veff0MT(imesh, ilh, itype) = Veff0%vr(imesh, ilh, itype, 1) * atoms%rmsh(imesh, itype) * atoms%rmsh(imesh, itype)
        end do ! imesh
      end do ! ilh
    end do ! itype

    call CalcGrR2FinLH( atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2Veff0MT, r2GrVeff0MT )

    do idir = 1, 3
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            lm_pre = oqn_l * (oqn_l + 1) + 1
            do mqn_m = -oqn_l, oqn_l
              lm = lm_pre + mqn_m
              do imesh = 1, atoms%jri(itype)
                grVeff0MT(imesh, lm, iatom, idir) = r2GrVeff0MT(imesh, lm, iatom, idir) / atoms%rmsh(imesh, itype)**2
              end do ! imesh
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! idir

    surfIntIR(:, :) = cmplx(0.0, 0.0)
    call CalcSurfIntIRDynMat( atoms, cell, ngdp, gdp, rho0IRpw, grVeff0IR, surfIntIR )
    write(*, '(a)') 'IR Surface integral'
    write(*, '(3(2(es16.8,1x),3x))') surfIntIR(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntIR(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntIR(3, :)

    surfIntMT(:, :) = cmplx(0.0, 0.0)
    call CalcSurfIntMTDynMat(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, rho0MT(:, :, :, 1), grVeff0MT, surfIntMT)
    write(*, '(a)') 'MT Surface integral'
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(1, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(2, :)
    write(*, '(3(2(es16.8,1x),3x))') surfIntMT(3, :)

    write(*, '(a)') 'MT + IR Surface integral'
    write(*, '(3(2(es16.8,1x),3x))') (surfIntMT(1, :) + surfIntIR(1, :))
    write(*, '(3(2(es16.8,1x),3x))') (surfIntMT(2, :) + surfIntIR(2, :))
    write(*, '(3(2(es16.8,1x),3x))') (surfIntMT(3, :) + surfIntIR(3, :))

  end subroutine Calc2ArgSFIntegrals

  subroutine testNplotGradGradTrlYlm(atoms)

    use m_jpPotDensHelper, only : calcGrFinSH

    implicit none

    ! Type parameters
    type(t_atoms),              intent(in) :: atoms

    ! Scalar variables
    integer                                :: idirR
    integer                                :: idirC
    integer                                :: oqn_l
    integer                                :: mqn_m
    integer                                :: imesh
    integer                                :: lm

    ! Array variables
    complex,       allocatable             :: testrl(:, :, :)
    complex,       allocatable             :: gradTestrl(:, :, :, :)
    complex,       allocatable             :: gradGradTtestrl(:, :, :, :)

    allocate( testrl(atoms%jmtd, (atoms%lmaxd + 1)**2, 1), gradTestrl(atoms%jmtd, (atoms%lmaxd + 2)**2, 1, 3), &
                                                                         & gradGradTtestrl(atoms%jmtd, (atoms%lmaxd + 2)**2, 1, 3) )

    gradTestrl(:, :, :, :) = cmplx(0., 0.)
    testrl(:, :, :) = cmplx(0., 0.)
    gradGradTtestrl(:, :, :, :) = cmplx(0., 0.)

    oqn_l = 2
    idirR = 3
    idirC = 3

    do mqn_m = -oqn_l, oqn_l
      lm = oqn_l * (oqn_l + 1) + mqn_m + 1
      do imesh = 1, atoms%jri(1)
        testrl(imesh, lm, 1) = atoms%rmsh(imesh, 1)**oqn_l
      end do ! imesh
    end do ! mqn_m

    call calcGrFinSH( atoms, testrl, gradTestrl )

    call calcGrFinSH( atoms, gradTestrl(:, :(atoms%lmaxd + 1)**2, :, idirC), gradGradTtestrl )

    oqn_l = 0
    mqn_m = 0
!    do oqn_l = 0, atoms%lmax(1)
!      do mqn_m = -oqn_l, oqn_l
        lm = oqn_l * (oqn_l + 1) + mqn_m + 1
        do imesh = 1, atoms%jri(1)
!          write(1978, '(i8, i8, i8,f15.8, f15.8)') oqn_l, mqn_m, imesh, testrl(imesh, lm, 1)
!          write(1570, '(i8, i8, i8,f15.8, f15.8)') oqn_l, mqn_m, imesh, gradGradTtestrl(imesh, lm, 1, idirR)
          write(1978, '(f15.8,f15.8, f15.8)') atoms%rmsh(imesh, 1), testrl(imesh, lm, 1)
          write(1570, '(f15.8,f15.8, f15.8)') atoms%rmsh(imesh, 1), gradGradTtestrl(imesh, lm, 1, idirR)
        end do ! imesh
!      end do ! mqn_m
!    end do ! oqn_l

  end subroutine testNplotGradGradTrlYlm

 subroutine testLegendrePoly(atoms, cell, dimens, ngdp, gdp)

   use m_jp2ndOrdQuant, only : GenPsDens2ndOrd
   use m_jpConstants, only : fpi, iu, tpi
   use m_sphbes
   use m_ylm_old

   implicit none

   ! Type parameter
   type(t_atoms),                  intent(in) :: atoms
   type(t_cell),                   intent(in) :: cell
   type(t_dimension),              intent(in) :: dimens

   ! Scalar parameter
   integer,                        intent(in) :: ngdp

   ! Array parameter
   integer,                        intent(in) :: gdp(:, :)

   ! Scalar variable
   integer                                    :: G0index
   integer                                    :: itype
   integer                                    :: iatomTemp
   integer                                    :: iatom
   integer                                    :: ieqat
   integer                                    :: iG
   integer                                    :: ii
   integer                                    :: idirC
   integer                                    :: idirR
   integer                                    :: mqn_m
   integer                                    :: lm
   complex                                    :: mSum
   real                                       :: prefactor
   logical                                    :: testMode

   ! Array variables
   real                                       :: qpt(3)
   complex                                    :: qlm(-2:2,3, 3)
   complex                                    :: ylm(9)

   complex,           allocatable             :: psDens2ndOrd(:, :, :, :)
   complex,           allocatable             :: psDens2ndOrdPre(:, :, :, :)
   real,              allocatable             :: Gext(:, :)
   real,              allocatable             :: GextRmt(:)
   real,              allocatable             :: sbes(:, :, :)
   complex,           allocatable             :: phaseFactor(:, :)

   write(*, *)' testmode'
   qpt = [0., 0., 0.]
   testMode = .false.

   call GenPsDens2ndOrd( atoms, cell, dimens, ngdp, G0index, gdp, qpt, psDens2ndOrd, testMode )

   qlm(:, :, :) = cmplx(0., 0.)

   qlm(-2, 1, 1) = cmplx(sqrt(30. / 4.), 0.)
   qlm( 2, 1, 1) = cmplx(sqrt(30. / 4.), 0.)
   qlm( 0, 1, 1) = cmplx(-sqrt(5.), 0.)

   qlm(-2, 2, 2) = cmplx(-sqrt(30. / 4.), 0.)
   qlm( 2, 2, 2) = cmplx(-sqrt(30. / 4.), 0.)
   qlm( 0, 2, 2) = cmplx(-sqrt(5.), 0.)

   qlm( 0, 3, 3) = cmplx(sqrt(20.), 0.)


   allocate( GextRmt(ngdp), Gext(3, ngdp) )
   allocate( sbes( 0:dimens%ncvd + 1, ngdp, atoms%ntype ) )
   allocate( phaseFactor(ngdp, atoms%nat) )

   iatom = 0
   do itype = 1, atoms%ntype

     ! Prefactor is initialized with neutral element of multiplication for every atom type
     prefactor = 1.

     ! Calculate double factorial (2 N + 7)!! in 7.78 (PhD thesis Klüppelberg). Note, atoms%ncv(itype) is already the value taken from
     ! Table 1 in the Weinert paper. For the orbital quantum number l = 2, this leads to
     ! (2 * atoms%ncv(itype) + 3)!! = (2 * N + 2 * l + 3)!! = (2 * N + 7)!!.
     ! The right hand side of the recent equation is consistent with 7.78 (PhD thesis Klüppelberg), where N is the Weinert parameter.
     ! todo MB prefactor should epend on itype
     do ii = 1, 2 * atoms%ncv(itype) + 3, 2
       prefactor = prefactor * ii
     end do ! ii

     ! Complete prefactor. (We include also here the fpi from the Poisson equation leading from the pseudo density to the second order external potential.
     prefactor = prefactor * atoms%zatom(itype) / cell%omtil * fpi / 5. / 3. / atoms%rmt(itype) / atoms%rmt(itype)

     iatomTemp = iatom
     ! does this work?
     do iG = 1, ngdp

       ! No pseudo density contribution for G = 0
       if (iG == 483) cycle

       !todo MB GpqextRmt auch von itype abhängig!
       Gext(:, iG) = matmul(cell%bmat, gdp(:, iG))
       GextRmt(iG) = norm2(Gext(:, iG)) * atoms%rmt(itype)
       call sphbes(atoms%ncv(itype) + 1, GextRmt(iG), sbes(:, iG, itype))
       iatom = iatomTemp
       do ieqat = 1, atoms%neq(itype)
         iatom = iatom + 1
         phaseFactor(iG, iatom) = exp(-tpi * iu * dot_product(gdp(:, iG), atoms%taual(:, iatom)))
       end do ! ieqat
     end do ! iG
   end do ! itype

   allocate( psDens2ndOrdPre(ngdp, 3, atoms%nat, 3) )
   psDens2ndOrdPre = cmplx(0., 0.)
   do idirC = 1, 3
     iatom = 0
     do itype = 1, atoms%ntype
       do ieqat = 1, atoms%neq(itype)
         iatom = iatom + 1
         do idirR = 1, 3
           do iG = 1, ngdp
             if (iG == 483) cycle
             call Ylm4( 2, Gext(:, iG), ylm )
             mSum = 0
             do mqn_m = -2, 2
               lm = 7 + mqn_m
               mSum = mSum + qlm(mqn_m, idirR, idirC) * ylm(lm)
             end do ! mqn_m
             ! Note that atoms%ncv = N + l, so for l = 2, we need only an exponent or the spherical Bessel function, respectively,
             ! of ncv + 1.
             psDens2ndOrdPre(iG, idirR, iatom, idirC) = prefactor / GextRmt(iG)**(atoms%ncv(itype) - 1) &
                                                          & * sbes(atoms%ncv(itype) + 1, iG, itype) * phaseFactor(iG, iatom) * mSum
           end do ! iG
         end do ! idirC
       end do ! ieqat
     end do ! itype
   end do ! idirC

   do idirC = 1, 3
     iatom = 0
     do itype = 1, atoms%ntype
       do ieqat = 1, atoms%neq(itype)
         iatom = iatom + 1
         !do idirR = 1, 3
           do iG = 1, ngdp
             write(4285, '(i8, i8, 2(es15.8))') idirC, iG, psDens2ndOrd(iG, idirC, iatom, idirC)
             write(4286, '(i8, i8, 2(es15.8))') idirC, iG, psDens2ndOrdPre(iG, idirC, iatom, idirC)
             if (real(psDens2ndOrdPre(iG, idirC, iatom, idirC)) /= 0.) then
             !write(4287, '(i8, i8, 2(es15.8))') idirC, iG, real(psDens2ndOrdPre(iG, idirC, iatom, idirC)) / (psDens2ndOrd(iG, idirC, iatom, idirC))
               write(4288, '(i8, i8, 2(es15.8))') idirC, iG, real(psDens2ndOrd(iG, idirC, iatom, idirC)) / real(psDens2ndOrdPre(iG, idirC, iatom, idirC))
              else
               write(4288, '(i8, i8, 2(es15.8))') idirC, iG, -9.
             end if
             write(4289, '(i8, i8, 2(es15.8))') idirC, iG, psDens2ndOrd(iG, idirC, iatom, idirC) - psDens2ndOrdPre(iG, idirC, iatom, idirC)
           end do ! iG
         !end do ! idirR
       end do ! ieqat
     end do ! itype
   end do ! idirC

 end subroutine testLegendrePoly

end module m_jpTestDynMatDeprecated
