module m_mpmom
  !     ***********************************************************
  !     calculation of the multipole moments of the original charge
  !     density minus the interstitial charge density
  !     for each atom type
  !
  !     For yukawa_residual = .true. the subroutines calculate the
  !     multipole moments for the Yukawa potential instead of the
  !     Coulomb potential. This is used in the preconditioning of
  !     the SCF iteration for metallic systems.
  !
  !     qlmo(m,l,n) : mult.mom. of the muffin-tin charge density
  !     qlmp(m,l,n) : mult.mom. of the plane-wave charge density
  !     qlm (m,l,n) : (output) difference of the former quantities
  !
  !     references:
  !     for both the Coulomb and the Yukawa pseudo charge density:
  !     F. Tran, P. Blaha: Phys. Rev. B 83, 235118 (2011)
  !     or see the original paper for the normal Coulomb case only:
  !     M. Weinert: J.Math.Phys. 22(11) (1981) p.2434 eq. (10)-(15)
  !     ***********************************************************

contains

  subroutine mpmom( input, fmpi, atoms, sphhar, stars, sym, cell,   qpw, rho, potdenType, qlm,l_coreCharge,&
                  & rhoimag, stars2, iDtype, iDir, rho0, qpw0, iDir2, mat2ord )

    use m_types
    USE m_constants
    implicit none

    type(t_input),   intent(in)   :: input
    type(t_mpi),     intent(in)   :: fmpi

    type(t_sym),     intent(in)   :: sym
    type(t_stars),   intent(in)   :: stars
    type(t_cell),    intent(in)   :: cell
    type(t_sphhar),  intent(in)   :: sphhar
    type(t_atoms),   intent(in)   :: atoms
    real,            intent(in)   :: rho(:,0:,:) !(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    complex,         intent(in)   :: qpw(:)      !(stars%ng3)
    integer,         intent(in)   :: potdenType
    complex,         intent(out)  :: qlm(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    LOGICAL, OPTIONAL, INTENT(IN) :: l_coreCharge

    REAL, OPTIONAL, INTENT(IN)          :: rhoimag(:,0:,:), rho0(:,0:,:)
    INTEGER, OPTIONAL, INTENT(IN)       :: iDtype, iDir ! DFPT: Type and direction of displaced atom
    COMPLEX, OPTIONAL, INTENT(IN)       :: qpw0(:)
    TYPE(t_stars), OPTIONAL, INTENT(IN) :: stars2
    INTEGER, OPTIONAL, INTENT(IN)       :: iDir2
    COMPLEX, OPTIONAL, INTENT(IN)       :: mat2ord(5,3,3)

    integer                       :: j, jm, lh, mb, mem, mems, n, nd, l, nat, m
    complex                       :: qlmo(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    complex                       :: qlmp(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    complex                       :: qlmp_SF(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)

    LOGICAL :: l_dfptvgen ! If this is true, we handle things differently!

    l_dfptvgen = PRESENT(stars2)

    ! multipole moments of original charge density
    if ( fmpi%irank == 0 ) then
      qlmo = 0.0
!      call mt_moments( input, atoms, sphhar, rho(:,:,:), potdenType, qlmo )
      IF (.NOT.l_dfptvgen) THEN
          call mt_moments( input, atoms, sym,sphhar, rho(:,:,:), potdenType,qlmo,l_coreCharge)
      ELSE IF (.NOT.PRESENT(iDir2)) THEN
          ! qlmo for the real part of rho1:
          call mt_moments( input, atoms, sym,sphhar, rho(:,:,:), potdenType,qlmo,l_coreCharge=.FALSE.)
          ! qlmo for the imaginary part of rho1 and the perturbation of vExt in the displaced atom:
          call mt_moments( input, atoms, sym,sphhar, rhoimag(:,:,:), potdenType,qlmo,l_coreCharge=.TRUE.,l_rhoimag=.TRUE.,iDtype=iDtype,iDir=iDir)
          CALL dfpt_mt_moments_SF(atoms, sym, sphhar, iDtype, iDir, rho0(:,:,:), qlmo)
      ELSE
          call mt_moments( input, atoms, sym,sphhar, rho(:,:,:), potdenType,qlmo,l_coreCharge=.TRUE.,l_rhoimag=.FALSE.,iDtype=iDtype,iDir=iDir,iDir2=iDir2,mat2ord=mat2ord)
      END IF
    end if

    ! multipole moments of the interstitial charge density in the spheres
    call pw_moments( input, fmpi, stars, atoms, cell, sym,   qpw(:), potdenType, qlmp , l_dfptvgen)

    IF (l_dfptvgen.AND..NOT.PRESENT(iDir2)) THEN
      CALL dfpt_pw_moments_SF( fmpi, stars2, atoms, cell, sym, iDtype, iDir, qpw0(:), qlmp_SF )
      qlmp = qlmp + qlmp_SF
    END IF

    if ( fmpi%irank == 0 ) then
      ! see (A14)
      qlm = qlmo - qlmp
      ! output section
      do n = 1, atoms%ntype
         nat = atoms%firstAtom(n)
         write(oUnit, fmt=8000 ) n
         do l = 0, atoms%lmax(n)
            do m = -l, l
               if ( qlmo(m,l,n)/=CMPLX(0.0) .or. qlmp(m,l,n)/=CMPLX(0.0) ) then
                  write(oUnit, fmt=8010 ) l, m, qlmo(m,l,n), qlmp(m,l,n)
               end if
            end do
         end do
      end do

8000   FORMAT (/,10x,'multipole moments for atom type=',i5,/,/,t3,'l',t7,&
            &       'm',t27,'original',t57,'plane wave')
8010   FORMAT (1x,i2,2x,i3,2x,2 (5x,2e15.5))
       !
    end if ! fmpi%irank == 0

  end subroutine mpmom


!  subroutine mt_moments( input, atoms, sphhar, rho, potdenType, qlmo )
  subroutine mt_moments( input, atoms, sym,sphhar, rho, potdenType,qlmo,l_coreCharge,l_rhoimag,iDtype,iDir,iDir2,mat2ord)
    ! multipole moments of original charge density
    ! see (A15) (Coulomb case) or (A17) (Yukawa case)

    use m_intgr,     only: intgr3
    use m_constants, only: sfp_const, POTDEN_TYPE_POTYUK, POTDEN_TYPE_CRYSTALFIELD
    use m_types
    use m_DoubleFactorial
    use m_SphBessel
    use m_juDFT
    implicit none

    type(t_input),  intent(in)        :: input
    type(t_sphhar), intent(in)        :: sphhar
    type(t_atoms),  intent(in)        :: atoms
    type(t_sym),    intent(in)        :: sym
    real,           intent(in)        :: rho(: ,0:, :)
    integer,        intent(in)        :: potdenType
    complex,        intent(inout)     :: qlmo(-atoms%lmaxd:,0:,:)
    LOGICAL, OPTIONAL, INTENT(IN)     :: l_coreCharge,l_rhoimag
    INTEGER, OPTIONAL, INTENT(IN)     :: iDtype, iDir ! DFPT: Type and direction of displaced atom
    INTEGER, OPTIONAL, INTENT(IN)     :: iDir2
    COMPLEX, OPTIONAL, INTENT(IN)     :: mat2ord(5,3,3)

    integer                           :: n, ns, jm, nl, l, j, mb, m, nat, i, imax, lmax
    real                              :: fint
    real                              :: f( maxval( atoms%jri ) )
    real, allocatable, dimension(:,:) :: il, kl
    LOGICAL                           :: l_subtractCoreCharge

    LOGICAL :: l_dfptvgen ! If this is true, we handle things differently!

    l_dfptvgen = PRESENT(iDtype)

    if ( potdenType == POTDEN_TYPE_POTYUK ) then
      allocate( il(0:atoms%lmaxd, 1:atoms%jmtd), kl(0:atoms%lmaxd, 1:atoms%jmtd) )
    end if

    l_subtractCoreCharge = .TRUE.
    if ( potdenType == POTDEN_TYPE_POTYUK ) l_subtractCoreCharge = .FALSE.
    IF(PRESENT(l_coreCharge)) l_subtractCoreCharge = l_coreCharge

    nat = 1
    do n = 1, atoms%ntype
      ns = sym%ntypsy(nat)
      jm = atoms%jri(n)
      imax = atoms%jri(n)
      lmax = sphhar%llh(sphhar%nlh(ns), ns)
      if ( potdenType == POTDEN_TYPE_POTYUK ) then
        !do concurrent (i = 1:imax)
        do i = 1,imax
          call ModSphBessel( il(:, i), kl(:, i), input%preconditioning_param * atoms%rmsh(i, n), lmax )
        end do
      end if
      do nl = 0, sphhar%nlh(ns)
        l = sphhar%llh(nl,ns)
        if(jm < 2) call juDFT_error("This would be uninit in intgr3.")
        do j = 1, jm
          if ( potdenType /= POTDEN_TYPE_POTYUK ) then
            f(j) = atoms%rmsh(j,n) ** l * rho(j,nl,n)
          else
            f(j) = il(l, j) * rho(j,nl,n)
          end if
        end do
        call intgr3( f, atoms%rmsh(:,n), atoms%dx(n), jm, fint )
        if ( potdenType == POTDEN_TYPE_POTYUK ) then
          fint = fint * DoubleFactorial( l ) / input%preconditioning_param ** l
        end if
        do mb = 1, sphhar%nmem(nl,ns)
          m = sphhar%mlh(mb,nl,ns)
          IF (.NOT.PRESENT(l_rhoimag)) THEN
              qlmo(m,l,n) = qlmo(m,l,n) + sphhar%clnu(mb,nl,ns) * fint
          ELSE
              qlmo(m,l,n) = qlmo(m,l,n) + ImagUnit*sphhar%clnu(mb,nl,ns) * fint
          END IF
        end do
      end do
!      if ( potdenType /= POTDEN_TYPE_POTYUK ) then
      if (l_subtractCoreCharge) then
        IF (.NOT.l_dfptvgen) THEN
            qlmo(0,0,n) = qlmo(0,0,n) - atoms%zatom(n) / sfp_const
        ELSE IF (.NOT.PRESENT(iDir2)) THEN
            IF ((n.EQ.iDtype)) THEN
                qlmo(-1:1,1,n) = qlmo(-1:1,1,n) - 3.0 / fpi_const * atoms%zatom(n) * c_im(iDir, :)
            ELSE IF ((0.EQ.iDtype)) THEN
                qlmo(-1:1,1,n) = qlmo(-1:1,1,n) + 3.0 / fpi_const * atoms%zatom(n) * c_im(iDir, :)
            END IF
         ELSE
            IF ((n.EQ.iDtype).OR.(0.EQ.iDtype)) qlmo(-2:2,2,n) = qlmo(-2:2,2,n) - 5.0 / fpi_const * atoms%zatom(n) * mat2ord(:,iDir2,iDir)
        END IF
      end if
      nat = nat + atoms%neq(n)
    end do

  end subroutine mt_moments


!  subroutine pw_moments( input, fmpi, stars, atoms, cell, sym,   qpw, potdenType, qlmp_out )
  subroutine pw_moments( input, fmpi, stars, atoms, cell, sym,   qpw_in, potdenType, qlmp_out, l_dfptvgen )
    ! multipole moments of the interstitial charge in the spheres

    use m_mpi_bc_tool
    use m_mpi_reduce_tool
    use m_phasy1
    use m_sphbes

    use m_constants, only: sfp_const, POTDEN_TYPE_POTYUK
    use m_types
    use m_DoubleFactorial
    use m_SphBessel
    implicit none

    type(t_input),    intent(in)   :: input
    type(t_mpi),      intent(in)   :: fmpi

    type(t_sym),      intent(in)   :: sym
    type(t_stars),    intent(in)   :: stars
    type(t_cell),     intent(in)   :: cell
    type(t_atoms),    intent(in)   :: atoms
    complex,          intent(in)   :: qpw_in(:)
    integer,          intent(in)   :: potdenType
    complex,          intent(out)  :: qlmp_out(-atoms%lmaxd:,0:,:)
    LOGICAL,          INTENT(IN)   :: l_dfptvgen

    integer                        :: n, k, l, ll1, lm, ierr, m
    integer                        :: maxBunchSize, numBunches, iBunch, firstStar, iTempArray
    complex                        :: sk3i, cil, nqpw, temp
    complex                        :: pylm(( maxval( atoms%lmax ) + 1 ) ** 2, atoms%ntype)
    real                           :: sk3r, rl2
    real                           :: aj(0:maxval( atoms%lmax ) + 1 )
    complex, ALLOCATABLE           :: qpw(:)
    real                           :: il(0:maxval( atoms%lmax ) + 1 )
    real                           :: kl(0:maxval( atoms%lmax ) + 1 )
    complex                        :: qlmp(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    complex, ALLOCATABLE           :: qlmpStars(:,:,:,:)
    TYPE(t_parallelLoop)           :: mpiLoop, ompLoop

    ALLOCATE(qpw(stars%ng3))
    qpw = qpw_in(:stars%ng3)
    qlmp = 0.0

    call mpi_bc(qpw,0,fmpi%mpi_comm)

    firstStar = MERGE(2,1,norm2(stars%center)<=1e-8)
    maxBunchSize = 2*getNumberOfThreads() ! This bunch size is kind of a magic number detrmined from some
                                          ! naive performance tests for a 64 atom unit cell
    CALL calcNumberComputationBunches(firstStar, stars%ng3, maxBunchSize, numBunches)

    ALLOCATE(qlmpStars(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype,maxBunchSize))
    qlmpStars = CMPLX(0.0,0.0)

    CALL mpiLoop%init(fmpi%irank,fmpi%isize,0,numBunches-1)
    DO iBunch = mpiLoop%bunchMinIndex, mpiLoop%bunchMaxIndex
       CALL ompLoop%init(iBunch,numBunches,firstStar,stars%ng3)
       !$OMP parallel do default( NONE ) &
       !$OMP SHARED(input,atoms,stars,sym,cell,ompLoop,qpw,qlmpStars,potdenType) &
       !$OMP private(iTempArray,pylm,nqpw,n,sk3r,aj,rl2,il,kl,sk3i,l,cil,ll1,m,lm)
       do k = ompLoop%bunchMinIndex, ompLoop%bunchMaxIndex
          iTempArray = k - ompLoop%bunchMinIndex + 1
          call phasy1( atoms, stars, sym, cell, k, pylm )
          nqpw = qpw(k) * stars%nstr(k)
          do n = 1, atoms%ntype
             sk3r = stars%sk3(k) * atoms%rmt(n)
             call sphbes( atoms%lmax(n) + 1, sk3r, aj )
             rl2 = atoms%rmt(n) ** 2
             if ( potdenType == POTDEN_TYPE_POTYUK ) then
                call ModSphBessel( il(0:atoms%lmax(n)+1), kl(0:atoms%lmax(n)+1), input%preconditioning_param * atoms%rmt(n), atoms%lmax(n) + 1 )
                sk3i = nqpw / ( stars%sk3(k) ** 2 + input%preconditioning_param ** 2 ) * rl2
             else
                sk3i = nqpw / stars%sk3(k)
             end if
             do l = 0, atoms%lmax(n)
                if ( potdenType == POTDEN_TYPE_POTYUK ) then
                   cil = ( stars%sk3(k) * il(l) * aj(l+1) + input%preconditioning_param * il(l+1) * aj(l) ) * ( DoubleFactorial( l ) / input%preconditioning_param ** l ) * sk3i
                else
                   cil = aj(l+1) * sk3i * rl2
                   rl2 = rl2 * atoms%rmt(n)
                end if
                ll1 = l * ( l + 1 ) + 1
                do m = -l, l
                   lm = ll1 + m
                   qlmpStars(m,l,n,iTempArray) = qlmpStars(m,l,n,iTempArray) + cil * pylm(lm,n)
                end do
             end do                  ! l = 0, atoms%lmax(n)
          end do                    ! n = 1, atoms%ntype
       end do
       !$OMP END PARALLEL DO
    END DO
    !$OMP parallel do default( NONE ) &
    !$OMP SHARED(atoms,maxBunchSize,qlmpStars,qlmp) &
    !$OMP private(l,m,temp,iTempArray)
    DO n = 1, atoms%ntype
       DO l = 0, atoms%lmax(n)
          DO m = -l, l
             temp = CMPLX(0.0,0.0)
             DO iTempArray = 1, maxBunchSize
                temp = temp + qlmpStars(m,l,n,iTempArray)
             END DO
             qlmp(m,l,n) = temp
          END DO
       END DO
    END DO
    !$OMP END PARALLEL DO

    if ( fmpi%irank == 0 .AND. .NOT. l_dfptvgen) then
      ! q=0 term: see (A19) (Coulomb case) or (A20) (Yukawa case)
      do n = 1, atoms%ntype
        if ( potdenType /= POTDEN_TYPE_POTYUK ) then
          qlmp(0,0,n) = qlmp(0,0,n) + qpw(1) * stars%nstr(1) * atoms%volmts(n) / sfp_const
        else
          call ModSphBessel( il(0:1), kl(0:1), input%preconditioning_param * atoms%rmt(n), 1 )
          qlmp(0,0,n) = qlmp(0,0,n) + qpw(1) * stars%nstr(1) * sfp_const * atoms%rmt(n) ** 2 * il(1) / input%preconditioning_param
        end if
      end do
    end if

    CALL mpi_sum_reduce(qlmp, qlmp_out, fmpi%mpi_comm)

  end subroutine pw_moments

   SUBROUTINE dfpt_mt_moments_SF(atoms, sym, sphhar, iDtype, iDir, rho0, qlmo)
      USE m_types
      USE m_gaunt, only : Gaunt1
      USE m_constants

      IMPLICIT NONE

      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_sym),    INTENT(IN)    :: sym
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      INTEGER,        INTENT(IN)    :: iDtype, iDir
      REAL,           INTENT(IN)    :: rho0(:, 0:, :)
      COMPLEX,        INTENT(INOUT) :: qlmo(-atoms%lmaxd:,0:, :)

      INTEGER :: mb, n, nat, nl, ns, jm, l, lp, m, mp, mVec, pref
      REAL    :: fint, gauntFactor

      pref = -1

      IF (iDtype.NE.0) THEN
         nat = iDtype
         pref = 1
      END IF

      DO n = MERGE(1,iDtype,iDtype.EQ.0), MERGE(atoms%ntype,iDtype,iDtype.EQ.0)
         nat = atoms%firstAtom(n)
         ns = sym%ntypsy(nat)
         jm = atoms%jri(n)
         DO nl = 0, sphhar%nlh(ns)
            lp = sphhar%llh(nl,ns)
            DO l = MERGE(1, lp - 1, lp.EQ.0), MERGE(1, lp + 1, lp.EQ.0), 2 ! Gaunt selection
               IF (l.GT.atoms%lmax(n)) CYCLE
               fint = atoms%rmt(n)**l * rho0(jm,nl,n)
               DO mb = 1, sphhar%nmem(nl,ns)
                  mp = sphhar%mlh(mb,nl,ns)
                  DO mVec = -1, 1
                     m = mVec + mp ! Gaunt selection
                     IF (ABS(m).GT.l) CYCLE
                     gauntFactor = Gaunt1(l, 1, lp, m, mVec, mp, atoms%lmax(n))
                     qlmo(m, l, n) = qlmo(m, l, n) + c_im(iDir, mVec + 2) * gauntFactor * &
                                                   & sphhar%clnu(mb,nl,ns) * fint * pref
                  END DO ! mVec
               END DO ! mb
            END DO ! l
         END DO ! nl
      END DO ! n

   END SUBROUTINE dfpt_mt_moments_SF

   SUBROUTINE dfpt_pw_moments_SF( fmpi, stars, atoms, cell, sym, iDtype, iDir, qpw_in, qlmp_SF )

      use m_mpi_bc_tool
      use m_mpi_reduce_tool
      use m_phasy1
      use m_sphbes
      use m_constants
      use m_types
      USE m_gaunt, only : Gaunt1

      implicit none

      type(t_mpi),      intent(in)   :: fmpi
      type(t_sym),      intent(in)   :: sym
      type(t_stars),    intent(in)   :: stars
      type(t_cell),     intent(in)   :: cell
      type(t_atoms),    intent(in)   :: atoms
      INTEGER,       INTENT(IN)    :: iDtype, iDir
      complex,          intent(in)   :: qpw_in(:)
      complex,          intent(out)  :: qlmp_SF(-atoms%lmaxd:,0:,:)

      integer                        :: n, k, l, ll1p, lmp, ierr, m, lp, mp, mVec, pref
      complex                        :: cil, nqpw
      complex                        :: pylm(( maxval( atoms%lmax ) + 1 ) ** 2, atoms%ntype)
      real                           :: sk3r, rl2
      real                           :: aj(0:maxval( atoms%lmax ) + 1 )
      complex, ALLOCATABLE           :: qpw(:)
      complex                        :: qlmp(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)

!      TYPE(t_atoms), INTENT(IN)    :: atoms
!      TYPE(t_cell),  INTENT(IN)    :: cell
!      INTEGER,       INTENT(IN)    :: ngdp

!      INTEGER,       INTENT(IN)    :: gdp(:, :)
!      COMPLEX,       INTENT(IN)    :: rho0IRpw(:)
!      COMPLEX,       INTENT(INOUT) :: qlmp(:, :,:)

!      INTEGER :: iG, l, m, lp, mp, m2p, lm, lmp, iDir
!      COMPLEX :: pref, phaseFac, temp1, temp2, temp3
      REAL    :: gauntFactor

      ALLOCATE (qpw(stars%ng3))
      qpw = qpw_in(:stars%ng3)
      qlmp = 0.0

      pref = -1
      IF (iDtype.NE.0) pref = 1

      if ( fmpi%irank == 0 ) then
         do n = MERGE(1,iDtype,iDtype.EQ.0), MERGE(atoms%ntype,iDtype,iDtype.EQ.0)
            DO mVec = -1, 1
               qlmp(mVec,1,n) = pref * c_im(iDir, mVec + 2) * qpw(1) * stars%nstr(1) * atoms%rmt(n)**3
            END DO
         end do
      end if

      call mpi_bc(qpw,0,fmpi%mpi_comm)

      do k = fmpi%irank+2, stars%ng3, fmpi%isize
         call phasy1( atoms, stars, sym, cell, k, pylm )

         nqpw = qpw(k) * stars%nstr(k)

         do n = MERGE(1,iDtype,iDtype.EQ.0), MERGE(atoms%ntype,iDtype,iDtype.EQ.0)

            sk3r = stars%sk3(k) * atoms%rmt(n)
            call sphbes( atoms%lmax(n), sk3r, aj )
            rl2 = atoms%rmt(n) ** 2

            DO lp = 0, atoms%lmax(n)
               cil = aj(lp) * nqpw * rl2
               ll1p = lp * ( lp + 1 ) + 1
               DO l = MERGE(1, lp - 1, lp.EQ.0), MERGE(1, lp + 1, lp.EQ.0), 2 ! Gaunt selection
                  IF (l.GT.atoms%lmax(n)) CYCLE
                  DO mp = -lp, lp
                     lmp = ll1p + mp
                     DO mVec = -1, 1
                        m = mVec + mp ! Gaunt selection
                        IF (ABS(m).GT.l) CYCLE
                        gauntFactor = Gaunt1(l, 1, lp, m, mVec, mp, atoms%lmax(n))
                        qlmp(m,l,n) = qlmp(m,l,n) + c_im(iDir, mVec + 2) * gauntFactor * &
                                                  & cil * atoms%rmt(n)**l * pylm(lmp,n) * pref
                     END DO ! mVec
                  END DO ! mp
               END DO ! l
            END DO ! lp
         END DO ! n = 1, atoms%ntype
      END DO ! k = 2, stars%ng3

      CALL mpi_sum_reduce(qlmp, qlmp_SF, fmpi%mpi_comm)

!      pref = fpi_const * atoms%rmt(iDtype) * atoms%rmt(iDtype)
!      DO iG = 1, ngdp
!          gext = matmul(cell%bmat, real(gdp(:, iG)))
!          gnorm = norm2(gExt)

!          call ylm4( atoms%lmax(iDtype), gExt(1:3), ylm )
!          call sphbes(atoms%lmax(iDtype), gnorm * atoms%rmt(iDtype), sbes)

!          phaseFac = exp( ImagUnit * tpi_const * dot_product(gdp(:, iG), atoms%taual(:, iDatom)))

!          DO l = 0, atoms%lmax(iDtype)
!             temp1 = pref * phaseFac * atoms%rmt(iDtype)**l * rho0IRpw(iG)
!             DO m = -l, l
!                  lm = l * (l + 1) + m + 1
!                  DO lp = 0, atoms%lmax(iDtype)
!                      temp2 = temp1 * sbes(lp) * ImagUnit**lp
!                      DO mp = -lp, lp
!                          lmp = lp * (lp + 1) + mp + 1
!                          temp3 = temp2 * conjg(ylm(lmp))
!                          DO m2p = -1, 1
!                              gauntFactor = Gaunt1( l, lp, 1, m, mp, m2p, atoms%lmax(iDtype))
!                              DO iDir = 1, 3
!                                  qlmp(lm, iDatom, iDir) = qlmp(lm, iDatom, iDir) + c_im(iDir, m2p + 2) * gauntFactor * temp3
!                              END DO ! iDir
!                          END DO ! m2p
!                      END DO ! mp
!                  END DO ! lp
!             END DO ! m
!          END DO ! l
!      END DO ! iG

   END SUBROUTINE dfpt_pw_moments_SF

end module m_mpmom
