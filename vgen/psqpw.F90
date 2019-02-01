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

  subroutine psqpw( mpi, atoms, sphhar, stars, vacuum, dimension, cell, input, sym, oneD, &
       &     qpw, rho, rht, l_xyav, potdenType, psq )

#include"cpp_double.h"
    use m_constants
    use m_phasy1
    use m_mpmom 
    use m_sphbes
    use m_qsf
    use m_od_phasy
    use m_od_cylbes
    use m_types
    use m_DoubleFactorial
    use m_SphBessel
    implicit none

    type(t_mpi),        intent(in)  :: mpi
    type(t_atoms),      intent(in)  :: atoms
    type(t_sphhar),     intent(in)  :: sphhar
    type(t_stars),      intent(in)  :: stars
    type(t_vacuum),     intent(in)  :: vacuum
    type(t_dimension),  intent(in)  :: dimension
    type(t_cell),       intent(in)  :: cell
    type(t_input),      intent(in)  :: input
    type(t_sym),        intent(in)  :: sym
    type(t_oneD),       intent(in)  :: oneD
    logical,            intent(in)  :: l_xyav
    complex,            intent(in)  :: qpw(stars%ng3) 
    real,               intent(in)  :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype) 
    real,               intent(in)  :: rht(vacuum%nmzd,2)
    integer,            intent(in)  :: potdenType
    complex,            intent(out) :: psq(stars%ng3)

    complex                         :: psint, sa, sl, sm
    real                            :: f, fact, fpo, gz, p, qvac, rmtl, s, fJ, gr, g
    integer                         :: ivac, k, l, n, n1, nc, ncvn, lm, ll1, nd, m, nz
    complex                         :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%ntype)
    complex                         :: qlm(-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype)
    real                            :: q2(vacuum%nmzd)
    real                            :: pn(0:atoms%lmaxd,atoms%ntype)
    real                            :: aj(0:atoms%lmaxd+DIMENSION%ncvd+1)
    real                            :: rht1(vacuum%nmz)
    real, allocatable, dimension(:) :: il, kl
    real                            :: g0(atoms%ntype)
#ifdef CPP_MPI
    include 'mpif.h'
    integer                         :: ierr(3)
    complex, allocatable            :: c_b(:)
#endif

    ! Calculate multipole moments
    call mpmom( input, mpi, atoms, sphhar, stars, sym, cell, oneD, qpw, rho, potdenType, qlm )
#ifdef CPP_MPI
    psq(:) = cmplx( 0.0, 0.0 )
    call MPI_BCAST( qpw, size(qpw), CPP_MPI_COMPLEX, 0, mpi%mpi_comm, ierr )
    nd = ( 2 * atoms%lmaxd + 1 ) * ( atoms%lmaxd + 1 ) * atoms%ntype
    call MPI_BCAST( qlm, nd, CPP_MPI_COMPLEX, 0, mpi%MPI_COMM, ierr )
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
    if( mpi%irank == 0 ) then
    s = 0.
    do n = 1, atoms%ntype
      if ( potdenType /= POTDEN_TYPE_POTYUK ) then
        s = s + atoms%neq(n) * real( qlm(0,0,n) )
      else
        s = s + atoms%neq(n) * real( qlm(0,0,n) ) * g0(n)
      end if
    end do
    !if( mpi%irank == 0 ) then
      psq(1) = qpw(1) + ( sfp_const / cell%omtil ) * s
    end if

    ! q/=0 term: see (A10) (Coulomb case) or (A11) (Yukawa case)
    fpo = 1. / cell%omtil
    !$omp parallel do default( shared ) private( pylm, sa, n, ncvn, aj, sl, l, n1, ll1, sm, m, lm )
    do k = mpi%irank+2, stars%ng3, mpi%isize
      if ( .not. oneD%odi%d1 ) then
        call phasy1( atoms, stars, sym, cell, k, pylm )
      else
        call od_phasy( atoms%ntype, stars%ng3, atoms%nat, atoms%lmaxd, atoms%ntype, atoms%neq, &
             atoms%lmax, atoms%taual, cell%bmat, stars%kv3, k, oneD%odi, oneD%ods, pylm )
      end if
      sa = 0.
      do n = 1, atoms%ntype
        ncvn = atoms%ncv(n)
        call sphbes( ncvn + 1 , stars%sk3(k) * atoms%rmt(n), aj )
        sl = 0.
        do l = 0, atoms%lmax(n)
          if ( l >= ncvn ) go to 60
          n1 = ncvn - l + 1
          ll1 = l * ( l + 1 ) + 1
          sm = 0.
          do m = -l, l
            lm = ll1 + m 
            sm = sm + qlm(m,l,n) * conjg( pylm(lm,n) )
          end do
60        sl = sl + pn(l,n) / ( stars%sk3(k) ** n1 ) * aj( ncvn + 1 ) * sm
        end do
        sa = sa + atoms%neq(n) * sl
      end do
      psq(k) = qpw(k) + fpo * sa
    end do
    !$omp end parallel do
#ifdef CPP_MPI
    allocate( c_b(stars%ng3) )
    call MPI_REDUCE( psq, c_b, stars%ng3, CPP_MPI_COMPLEX, MPI_SUM, 0, mpi%MPI_COMM, ierr )
    if ( mpi%irank == 0 ) then
       psq(:stars%ng3) = c_b(:stars%ng3)
    end if
    deallocate( c_b )
#endif

    if ( mpi%irank == 0 ) then
      if ( potdenType == POTDEN_TYPE_POTYUK ) return
      ! Check: integral of the pseudo charge density within the slab
      if ( input%film .and. .not. oneD%odi%d1 ) then
        psint = psq(1) * stars%nstr(1) * vacuum%dvac
        do k = 2, stars%ng3
          if ( stars%ig2(k) == 1 ) then
            gz = stars%kv3(3,k) * cell%bmat(3,3)
            f = 2. * sin( gz * cell%z1 ) / gz
            psint = psint + stars%nstr(k) * psq(k) * f
          end if
        end do
        psint = cell%area * psint
      else if ( input%film .and. oneD%odi%d1 ) then
        psint = (0.0, 0.0)
        do k = 2, stars%ng3
          if ( stars%kv3(3,k) == 0 ) then
            g = ( stars%kv3(1,k) * cell%bmat(1,1) + stars%kv3(2,k) * cell%bmat(2,1) ) ** 2 + &
              & ( stars%kv3(1,k) * cell%bmat(1,2) + stars%kv3(2,k) * cell%bmat(2,2) ) ** 2
            gr = sqrt( g )
            call od_cylbes( 1, gr * cell%z1, fJ )
            f = 2 * cell%vol * fJ / ( gr * cell%z1 )
            psint = psint + stars%nstr(k) * psq(k) * f
          end if
        end do
        psint = psint + psq(1) * stars%nstr(1) * cell%vol
      else if ( .not. input%film ) then
        psint = psq(1) * stars%nstr(1) * cell%omtil
      end if
      write( 6, fmt=8000 ) psint
8000  format (/,10x,'integral of pseudo charge density inside the slab='&
            &       ,5x,2f11.6)
      if ( .not. input%film .or. potdenType == POTDEN_TYPE_POTYUK ) return

      ! Normalized pseudo density
      if ( .not. oneD%odi%d1 ) then
        qvac = 0.0
        do ivac = 1, vacuum%nvac
          call qsf( vacuum%delz, rht(1,ivac), q2, vacuum%nmz, 0 )
          q2(1) = q2(1) * cell%area
          qvac = qvac + q2(1) * 2. / real( vacuum%nvac )
        end do
        qvac = qvac - 2 * input%sigma
      else
        qvac = 0.0
        do nz = 1, vacuum%nmz
          rht1(nz) = ( cell%z1 + ( nz - 1 ) * vacuum%delz ) * rht(nz,vacuum%nvac)
        end do
        call qsf( vacuum%delz, rht1(1), q2, vacuum%nmz, 0 )
        qvac = cell%area * q2(1)
      end if
      if ( l_xyav ) return
      fact = ( qvac + psint ) / ( stars%nstr(1) * cell%vol )
      psq(1) = psq(1) - fact
      write( 6, fmt=8010 ) fact * 1000
8010  format (/,10x,'                     1000 * normalization const. ='&
            &       ,5x,2f11.6)
    end if ! mpi%irank == 0 

  end subroutine psqpw

end module m_psqpw
