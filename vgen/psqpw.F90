module m_psqpw
  !     ***********************************************************
  !     generates the fourier coefficients of pseudo charge density
  !
  !     For yukawa_residual = .true. the subroutines calculate the
  !     pseudo charge density for the generation of the Yukawa
  !     potential instead of the Coulomb potential. This is used in
  !     the preconditioning of the SCF iteration for metallic systems.
  !
  !     references:
  !     for both the Coulomb and Yukawa cases:
  !     F. Tran, P. Blaha: Phys. Rev. B 83, 235118 (2011)
  !     or see the original paper for the normal Coulomb case only:
  !     M. Weinert: J. Math. Phys. 22(11) (1981) p.2434 eq. (10)-(15)
  !     ***********************************************************

contains

  subroutine psqpw( fmpi, atoms, sphhar, stars, vacuum,  cell, input, sym,   &
       &     den, ispin, l_xyav, potdenType, psq, rhoimag, stars2, iDtype, iDir, rho0, qpw0 )

#ifdef CPP_MPI
    use mpi
#endif
    use m_constants
    use m_phasy1
    use m_mpmom
    use m_sphbes
    use m_qsf
    USE m_mpi_reduce_tool
     
     
    use m_types
    use m_DoubleFactorial
    use m_SphBessel
    implicit none

    type(t_mpi),        intent(in)  :: fmpi
    type(t_atoms),      intent(in)  :: atoms
    type(t_sphhar),     intent(in)  :: sphhar
    type(t_stars),      intent(in)  :: stars
    type(t_vacuum),     intent(in)  :: vacuum
    type(t_cell),       intent(in)  :: cell
    type(t_input),      intent(in)  :: input
    type(t_sym),        intent(in)  :: sym
    type(t_potden),     intent(in)  :: den
     
    logical,            intent(in)  :: l_xyav
    !complex,            intent(in)  :: qpw(stars%ng3)
    !real,               intent(in)  :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    !real,               intent(in)  :: rht(vacuum%nmzd,2)
    integer,            intent(in)  :: potdenType, ispin
    complex,            intent(out) :: psq(stars%ng3)

    REAL, OPTIONAL, INTENT(IN)      :: rhoimag(atoms%jmtd,0:sphhar%nlhd,atoms%ntype), rho0(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)

    TYPE(t_stars), OPTIONAL, INTENT(IN) :: stars2
    COMPLEX, OPTIONAL, INTENT(IN)   :: qpw0(:)

    INTEGER, OPTIONAL, INTENT(IN)     :: iDtype, iDir ! DFPT: Type and direction of displaced atom

    complex                         :: psint, sa, sl, sm, qvac, fact
    real                            :: f, fpo, gz, p, rmtl, s, fJ, gr, g
    integer                         :: ivac, k, l, n, n1, nc, ncvn, lm, ll1, nd, m, nz, kStart, kEnd
    complex                         :: psq_local(stars%ng3)
    complex                         :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%ntype)
    complex                         :: qlm(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    complex                         :: q2(vacuum%nmzd)
    real                            :: q2r(vacuum%nmzd),q2i(vacuum%nmzd)
    real                            :: pn(0:atoms%lmaxd,atoms%ntype)
    real                            :: aj(0:atoms%lmaxd+maxval(atoms%ncv)+1)
    real, allocatable, dimension(:) :: il, kl
    real                            :: g0(atoms%ntype)
    complex                         :: qpw(stars%ng3)
    real                            :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    complex                         :: rht(vacuum%nmzd,2)
    LOGICAL :: l_dfptvgen ! If this is true, we handle things differently!

#ifdef CPP_MPI
    integer                         :: ierr
#endif

    l_dfptvgen = PRESENT(stars2)
    qpw = den%pw(:,ispin)
    rho = den%mt(:,:,:,ispin)
    IF (input%film) rht = den%vac(:,1,:,ispin)

    ! Calculate multipole moments
    call timestart("mpmom")
    IF (.NOT.l_dfptvgen) THEN
        call mpmom( input, fmpi, atoms, sphhar, stars, sym, cell,   qpw, rho, potdenType, qlm )
    ELSE
        call mpmom( input, fmpi, atoms, sphhar, stars, sym, cell,   qpw, rho, potdenType, qlm, &
                  & rhoimag=rhoimag, stars2=stars2, iDtype=iDtype, iDir=iDir, rho0=rho0, qpw0=qpw0 )
    END IF
    call timestop("mpmom")

    psq(:) = cmplx( 0.0, 0.0 )
    psq_local(:) = cmplx( 0.0, 0.0 )
#ifdef CPP_MPI
    call MPI_BCAST( qpw, size(qpw), MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr )
    nd = ( 2 * atoms%lmaxd + 1 ) * ( atoms%lmaxd + 1 ) * atoms%ntype
    call MPI_BCAST( qlm, nd, MPI_DOUBLE_COMPLEX, 0, fmpi%MPI_COMM, ierr )
#endif

    ! prefactor in (A10) (Coulomb case) or (A11) (Yukawa case)
    ! nc(n) is the integer p in the paper; ncv(n) is l + p
    ! Coulomb case: pn(l,n) = (2 * l + 2 * p + 3)!! / ( (2 * l + 1)!! * R ** (ncv(n) + 1) ),
    ! Yukawa case: pn(l,n) = lambda ** (l + p + 1) / ( i_{l+p+1}(lambda * R) * (2 * l + 1)!! )
    ! g0 is the prefactor for the q=0 component in (A13)
    pn = 0.
    do n = 1, atoms%ntype
      if ( potdenType /= POTDEN_TYPE_POTYUK ) then
        do l = 0, min( atoms%ncv(n) - 1, atoms%lmax(n) )
          pn(l, n) = DoubleFactorial( atoms%ncv(n) + 1, l ) / ( atoms%rmt(n) ** ( atoms%ncv(n) + 1 ) )
        end do
      else
        allocate( il(0:atoms%ncv(n)+1), kl(0:atoms%ncv(n)+1) )
        call ModSphBessel( il(0:), kl(0:), input%preconditioning_param * atoms%rmt(n), atoms%ncv(n) + 1 )
        g0(n) = ( input%preconditioning_param * atoms%rmt(n) ) ** ( atoms%ncv(n) + 1 ) / DoubleFactorial( atoms%ncv(n) + 1 ) / il( atoms%ncv(n) + 1 ) !p=ncv(n)
        do l = 0, min( atoms%ncv(n) - 1, atoms%lmax(n) )
          pn(l, n) = input%preconditioning_param ** ( atoms%ncv(n) + 1 ) / il( atoms%ncv(n) + 1 ) / DoubleFactorial( l )
        end do
        deallocate( il, kl )
      end if
    end do

    ! q=0 term: see (A12) (Coulomb case) or (A13) (Yukawa case)
    if( fmpi%irank == 0 .AND. norm2(stars%center)<=1e-8) then
       s = 0.
       do n = 1, atoms%ntype
         if ( potdenType /= POTDEN_TYPE_POTYUK ) then
           s = s + atoms%neq(n) * real( qlm(0,0,n) )
         else
           s = s + atoms%neq(n) * real( qlm(0,0,n) ) * g0(n)
         end if
       end do
    !if( fmpi%irank == 0 ) then
       psq_local(1) = qpw(1) + ( sfp_const / cell%omtil ) * s
    end if

    ! q/=0 term: see (A10) (Coulomb case) or (A11) (Yukawa case)
    fpo = 1. / cell%omtil

    call timestart("loop in psqpw")

    CALL calcIndexBounds(fmpi, MERGE(2,1,norm2(stars%center)<=1e-8), stars%ng3, kStart, kEnd)
    !$OMP parallel do default( NONE ) &
    !$OMP SHARED(atoms,stars,sym,cell,kStart,kEnd,psq_local,qpw,qlm,pn,fpo) &
    !$OMP private( pylm, sa, n, ncvn, aj, sl, l, n1, ll1, sm, m, lm )
    do k = kStart, kEnd
      call phasy1( atoms, stars, sym, cell, k, pylm )
      sa = 0.0
      do n = 1, atoms%ntype
        ncvn = atoms%ncv(n)
        call sphbes( ncvn + 1 , stars%sk3(k) * atoms%rmt(n), aj )
        sl = 0.
        do l = 0, atoms%lmax(n)
          IF(l.LT.ncvn) THEN
             n1 = ncvn - l + 1
             ll1 = l * ( l + 1 ) + 1
             sm = 0.
             do m = -l, l
               lm = ll1 + m
               sm = sm + qlm(m,l,n) * conjg( pylm(lm,n) )
             end do
          END IF
        sl = sl + pn(l,n) / ( stars%sk3(k) ** n1 ) * aj( ncvn + 1 ) * sm
        end do
        sa = sa + atoms%neq(n) * sl
      end do
      psq_local(k) = qpw(k) + fpo * sa
    end do
    !$omp end parallel do

    CALL mpi_sum_reduce(psq_local,psq,fmpi%MPI_COMM)

    call timestop("loop in psqpw")

    IF (l_dfptvgen) RETURN ! TODO: Change this!

    if ( fmpi%irank == 0 ) then
      if ( potdenType == POTDEN_TYPE_POTYUK ) return
      ! Check: integral of the pseudo charge density within the slab
      if ( input%film ) then
        psint = psq(1) * stars%nstr(1) * vacuum%dvac
        do k = 2, stars%ng3
          if ( stars%ig2(k) == 1 ) then
            gz = stars%kv3(3,k) * cell%bmat(3,3)
            f = 2. * sin( gz * cell%z1 ) / gz
            psint = psint + stars%nstr(k) * psq(k) * f
          end if
        end do
        psint = cell%area * psint
      else if ( .not. input%film ) then
        psint = psq(1) * stars%nstr(1) * cell%omtil
      end if
      write(oUnit, fmt=8000 ) psint
8000  format (/,10x,'integral of pseudo charge density inside the slab='&
            &       ,5x,2f11.6)
      if ( .not. input%film .or. potdenType == POTDEN_TYPE_POTYUK ) return

      ! Normalized pseudo density
        qvac = cmplx(0.0,0.0)
        do ivac = 1, vacuum%nvac
          if (.not.l_dfptvgen) then
            call qsf( vacuum%delz, real(rht(:,ivac)), q2r, vacuum%nmz, 0 )
            q2 = q2r
          else
            call qsf( vacuum%delz, real(rht(:,ivac)), q2r, vacuum%nmz, 0 )
            call qsf( vacuum%delz,aimag(rht(:,ivac)), q2i, vacuum%nmz, 0 )
            q2 = q2r + ImagUnit*q2i
          end if
          q2(1) = q2(1) * cell%area
          qvac = qvac + q2(1) * 2. / real( vacuum%nvac )
        end do
        !TODO: reactivate electric fields
        !qvac = qvac - 2 * input%sigma
      if ( l_xyav ) return
      fact = ( qvac + psint ) / ( stars%nstr(1) * cell%vol )
      psq(1) = psq(1) - fact
      if (.not.l_dfptvgen) write(oUnit, fmt=8010 ) fact * 1000
8010  format (/,10x,'                     1000 * normalization const. ='&
            &       ,5x,2f11.6)
    end if ! fmpi%irank == 0

  end subroutine psqpw

end module m_psqpw
