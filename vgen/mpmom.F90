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

  subroutine mpmom( input, mpi, atoms, sphhar, stars, sym, cell, oneD, qpw, rho, potdenType, qlm,l_coreCharge )

    use m_types
    implicit none

    type(t_input),   intent(in)   :: input
    type(t_mpi),     intent(in)   :: mpi
    type(t_oneD),    intent(in)   :: oneD
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

    integer                       :: j, jm, lh, mb, mem, mems, n, nd, l, nat, m
    complex                       :: qlmo(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    complex                       :: qlmp(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)

    ! multipole moments of original charge density
    if ( mpi%irank == 0 ) then
!      call mt_moments( input, atoms, sphhar, rho(:,:,:), potdenType, qlmo )
      call mt_moments( input, atoms, sphhar, rho(:,:,:), potdenType,qlmo,l_coreCharge)
    end if

    ! multipole moments of the interstitial charge density in the spheres
    call pw_moments( input, mpi, stars, atoms, cell, sym, oneD, qpw(:), potdenType, qlmp )

    if ( mpi%irank == 0 ) then
      ! see (A14)
      qlm = qlmo - qlmp
      ! output section
      nat = 1
      do n = 1, atoms%ntype
        write( 6, fmt=8000 ) n
        nd = atoms%ntypsy(nat)
        do lh = 0, sphhar%nlh(nd)
          l = sphhar%llh(lh,nd)
          mems = sphhar%nmem(lh,nd)
          do mem = 1, mems
            m = sphhar%mlh(mem,lh,nd)
            write( 6, fmt=8010 ) l, m, qlmo(m,l,n), qlmp(m,l,n)
          end do
        end do
        nat = nat + atoms%neq(n)
      end do

8000   FORMAT (/,10x,'multipole moments for atom type=',i5,/,/,t3,'l',t7,&
            &       'm',t27,'original',t57,'plane wave')
8010   FORMAT (1x,i2,2x,i3,2x,2 (5x,2e15.5))
       !
    end if ! mpi%irank == 0

  end subroutine mpmom


!  subroutine mt_moments( input, atoms, sphhar, rho, potdenType, qlmo )
  subroutine mt_moments( input, atoms, sphhar, rho, potdenType,qlmo,l_coreCharge)
    ! multipole moments of original charge density
    ! see (A15) (Coulomb case) or (A17) (Yukawa case)

    use m_intgr,     only: intgr3
    use m_constants, only: sfp_const, POTDEN_TYPE_POTYUK
    use m_types
    use m_DoubleFactorial
    use m_SphBessel
    implicit none

    type(t_input),  intent(in)        :: input
    type(t_sphhar), intent(in)        :: sphhar
    type(t_atoms),  intent(in)        :: atoms
    real,           intent(in)        :: rho(: ,0:, :)
    integer,        intent(in)        :: potdenType
    complex,        intent(out)       :: qlmo(-atoms%lmaxd:,0:,:)
    LOGICAL, OPTIONAL, INTENT(IN)     :: l_coreCharge

    integer                           :: n, ns, jm, nl, l, j, mb, m, nat, i, imax, lmax
    real                              :: fint
    real                              :: f( maxval( atoms%jri ) )
    real, allocatable, dimension(:,:) :: il, kl
    LOGICAL                           :: l_subtractCoreCharge

    if ( potdenType == POTDEN_TYPE_POTYUK ) then
      allocate( il(0:atoms%lmaxd, 1:atoms%jmtd), kl(0:atoms%lmaxd, 1:atoms%jmtd) )
    end if

    l_subtractCoreCharge = .TRUE.
    if ( potdenType == POTDEN_TYPE_POTYUK ) l_subtractCoreCharge = .FALSE.
    IF(PRESENT(l_coreCharge)) l_subtractCoreCharge = l_coreCharge

    qlmo = 0.0
    nat = 1
    do n = 1, atoms%ntype
      ns = atoms%ntypsy(nat)
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
          qlmo(m,l,n) = qlmo(m,l,n) + sphhar%clnu(mb,nl,ns) * fint
        end do
      end do
!      if ( potdenType /= POTDEN_TYPE_POTYUK ) then
      if (l_subtractCoreCharge) then
        qlmo(0,0,n) = qlmo(0,0,n) - atoms%zatom(n) / sfp_const
      end if
      nat = nat + atoms%neq(n)
    end do

  end subroutine mt_moments


!  subroutine pw_moments( input, mpi, stars, atoms, cell, sym, oneD, qpw, potdenType, qlmp_out )
  subroutine pw_moments( input, mpi, stars, atoms, cell, sym, oneD, qpw_in, potdenType, qlmp_out )
    ! multipole moments of the interstitial charge in the spheres

    use m_phasy1
    use m_sphbes
    use m_od_phasy
    use m_constants, only: sfp_const, POTDEN_TYPE_POTYUK
    use m_types
    use m_DoubleFactorial
    use m_SphBessel
    implicit none

    type(t_input),    intent(in)   :: input
    type(t_mpi),      intent(in)   :: mpi
    type(t_oneD),     intent(in)   :: oneD
    type(t_sym),      intent(in)   :: sym
    type(t_stars),    intent(in)   :: stars
    type(t_cell),     intent(in)   :: cell
    type(t_atoms),    intent(in)   :: atoms
    complex,          intent(in)   :: qpw_in(:)
    integer,          intent(in)   :: potdenType
    complex,          intent(out)  :: qlmp_out(-atoms%lmaxd:,0:,:)

    integer                        :: n, k, l, ll1, lm, ierr(3), m
    complex                        :: sk3i, cil, nqpw
    complex                        :: pylm(( maxval( atoms%lmax ) + 1 ) ** 2, atoms%ntype)
    real                           :: sk3r, rl2
    real                           :: aj(0:maxval( atoms%lmax ) + 1 )
    complex                        :: qpw(stars%ng3)
    logical                        :: od
    real                           :: il(0:maxval( atoms%lmax ) + 1 )
    real                           :: kl(0:maxval( atoms%lmax ) + 1 )
#ifdef CPP_MPI
    include 'mpif.h'
#endif
    complex                        :: qlmp(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)

    qpw = qpw_in(:stars%ng3)
    qlmp = 0.0
    if ( mpi%irank == 0 ) then
      ! q=0 term: see (A19) (Coulomb case) or (A20) (Yukawa case)
      do n = 1, atoms%ntype
        if ( potdenType /= POTDEN_TYPE_POTYUK ) then  
          qlmp(0,0,n) = qpw(1) * stars%nstr(1) * atoms%volmts(n) / sfp_const
        else
          call ModSphBessel( il(0:1), kl(0:1), input%preconditioning_param * atoms%rmt(n), 1 )
          qlmp(0,0,n) = qpw(1) * stars%nstr(1) * sfp_const * atoms%rmt(n) ** 2 * il(1) / input%preconditioning_param
        end if
      end do
    end if

#ifdef CPP_MPI
    call MPI_BCAST( qpw, size(qpw), MPI_DOUBLE_COMPLEX, 0, mpi%mpi_comm, ierr )
#endif

    ! q/=0 terms: see (A16) (Coulomb case) or (A18) (Yukawa case)
    od = oneD%odi%d1
!    !$omp parallel do default( shared ) private( pylm, nqpw, n, sk3r, aj, rl2, sk3i, &
!    !$omp& l, cil, ll1, m, lm, k ) reduction( +:qlmp )
    do k = mpi%irank+2, stars%ng3, mpi%isize
      if ( od ) then
        call od_phasy( atoms%ntype, stars%ng3, atoms%nat, atoms%lmaxd, atoms%ntype, &
             atoms%neq, atoms%lmax, atoms%taual, cell%bmat, stars%kv3, k, oneD%odi, oneD%ods, pylm)
      else
        call phasy1( atoms, stars, sym, cell, k, pylm )
      end if
     
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
            qlmp(m,l,n) = qlmp(m,l,n) + cil * pylm(lm,n)
          end do
        end do                  ! l = 0, atoms%lmax(n)
      end do                    ! n = 1, atoms%ntype
    end do                      ! k = 2, stars%ng3
!    !$omp end parallel do
#ifdef CPP_MPI
    CALL MPI_REDUCE( qlmp, qlmp_out, SIZE(qlmp), MPI_DOUBLE_COMPLEX, MPI_SUM,0, mpi%mpi_comm, ierr )
#else
    qlmp_out = qlmp
#endif

  end subroutine pw_moments

end module m_mpmom
