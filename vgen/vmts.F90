module m_vmts

contains

  subroutine vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, vpw, rho, potdenType, vr )

  !-------------------------------------------------------------------------
  ! This subroutine calculates the lattice harmonics expansion coefficients 
  ! of the Coulomb / Yukawa potential for all atom types.                     
  !
  ! In more detail:
  ! the radius-dependent coefficient of the potential is
  ! V_{lm}(r) = \int_0^R G_l(r,r') \rho_{lm}(r') r'^2 dr' + V_{lm}(R) P_l(r)
  ! where 
  ! [sphere boundary contribution - first part of code:]
  ! V_{lm}(R) is the spherical harmonics expansion of the interstitial 
  ! potential,
  ! P_l(r) is the derivative of the spherical Green's function on the sphere 
  ! boundary (times a constant),
  ! [sphere interior contribution - second part of code:]
  ! G_l(r,r') ) = green_factor u_1(r_<) u_2(r_>) is the spherical Green's 
  ! function with homogeneous solutions u_1 and u_2 and 
  ! r_< = min(r,r') and r_> = max(r,r')
  ! The integral is split in a part where r'=r_< and a part where r'=r_> and
  ! the integral from r to R is split in \int_0^R - \int_0^r. Resulting 
  ! terms depending solely on R (not on r) are summarised in the variable 
  ! termsR.
  !
  ! More information and equations can be found in
  ! F. Tran, P. Blaha: Phys. Rev. B 83, 235118(2011) 
  !-------------------------------------------------------------------------

#include"cpp_double.h"
    use m_constants
    use m_types
    use m_intgr, only : intgr2
    use m_phasy1
    use m_sphbes
    use m_od_phasy
    use m_SphBessel
    implicit none

    type(t_input),  intent(in)        :: input
    type(t_mpi),    intent(in)        :: mpi
    type(t_stars),  intent(in)        :: stars
    type(t_sphhar), intent(in)        :: sphhar
    type(t_atoms),  intent(in)        :: atoms
    type(t_sym),    intent(in)        :: sym
    type(t_cell),   intent(in)        :: cell
    type(t_oneD),   intent(in)        :: oneD
    complex,        intent(in)        :: vpw(:)!(stars%ng3,input%jspins)
    real,           intent(in)        :: rho(:,0:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    integer,        intent(in)        :: potdenType
    real,           intent(out)       :: vr(:,0:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)

    complex                           :: cp, sm
    integer                           :: i, jm, k, l, lh, n, nd, lm, n1, nat, m, imax, lmax
    complex                           :: vtl(0:sphhar%nlhd, atoms%ntype)
    complex                           :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%ntype)
    real                              :: green_factor, termsR
    real                              :: green_1    (1:atoms%jmtd), green_2    (1:atoms%jmtd)
    real                              :: integrand_1(1:atoms%jmtd), integrand_2(1:atoms%jmtd)
    real                              :: integral_1 (1:atoms%jmtd), integral_2 (1:atoms%jmtd)
    real                              :: sbf(0:atoms%lmaxd)
    real, allocatable, dimension(:,:) :: il, kl
    
    !$ complex, allocatable :: vtl_loc(:,:)
#ifdef CPP_MPI
    include 'mpif.h'
    integer                       :: ierr(3)
    complex, allocatable          :: c_b(:)

    external MPI_REDUCE
#endif
    integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM



    ! SPHERE BOUNDARY CONTRIBUTION to the coefficients calculated from the values
    ! of the interstitial Coulomb / Yukawa potential on the sphere boundary

    vtl(:,:) = cmplx( 0.0, 0.0 )

    
    ! q=0 component
    if ( mpi%irank == 0 ) then
      vtl(0,1:atoms%ntype) = sfp_const * vpw(1)
    end if

    ! q/=0 components
    !$omp parallel default( none ) &
    !$omp& shared( mpi, stars, vpw, oneD, atoms, sym, cell, sphhar, vtl ) &
    !$omp& private( k, cp, pylm, nat, n, sbf, nd, lh, sm, jm, m, lm, l ) &
    !$omp& private( vtl_loc )
    !$ allocate(vtl_loc(0:sphhar%nlhd,atoms%ntype)) 
    !$ vtl_loc(:,:) = cmplx(0.0,0.0)
    !$omp do
    do k = mpi%irank+2, stars%ng3, mpi%isize
      cp = vpw(k) * stars%nstr(k)
      if ( .not. oneD%odi%d1 ) then
        call phasy1( atoms, stars, sym, cell, k, pylm )
      else
        call od_phasy( atoms%ntype, stars%ng3, atoms%nat, atoms%lmaxd, atoms%ntype, &
            atoms%neq, atoms%lmax, atoms%taual, cell%bmat, stars%kv3, k, oneD%odi, oneD%ods, pylm )
      end if
      nat = 1
      do n = 1, atoms%ntype
        call sphbes( atoms%lmax(n), stars%sk3(k) * atoms%rmt(n), sbf )
        nd = atoms%ntypsy(nat)
        do lh = 0, sphhar%nlh(nd)
          l = sphhar%llh(lh,nd)
          sm = (0.0,0.0)
          do jm = 1, sphhar%nmem(lh,nd)
            m = sphhar%mlh(jm,lh,nd)
            lm = l * ( l + 1 ) + m + 1 
            sm = sm + conjg( sphhar%clnu(jm,lh,nd) ) * pylm(lm,n)
          end do
          !$ if (.false.) then
          vtl(lh,n) = vtl(lh,n) + cp * sbf(l) * sm
          !$ end if
          !$ vtl_loc(lh,n) = vtl_loc(lh,n) + cp * sbf(l) * sm
        end do
        nat = nat + atoms%neq(n)
      end do
    end do
    !$omp end do
    !$omp critical
    !$ vtl = vtl + vtl_loc
    !$omp end critical
    !$ deallocate(vtl_loc) 
    !$omp end parallel
#ifdef CPP_MPI
    n1 = ( sphhar%nlhd + 1 ) * atoms%ntype
    allocate( c_b(n1) )
    call MPI_REDUCE( vtl, c_b, n1, CPP_MPI_COMPLEX, MPI_SUM, 0, mpi%mpi_comm, ierr )
    if ( mpi%irank == 0 ) vtl = reshape( c_b, (/sphhar%nlhd+1,atoms%ntype/) )
    deallocate( c_b )
#endif



    ! SPHERE INTERIOR CONTRIBUTION to the coefficients calculated from the 
    ! values of the sphere Coulomb/Yukawa potential on the sphere boundary

    if( mpi%irank == 0 ) then
    if ( potdenType == POTDEN_TYPE_POTYUK ) then
      allocate( il(0:atoms%lmaxd, 1:atoms%jmtd), kl(0:atoms%lmaxd, 1:atoms%jmtd) )
    end if

    nat = 1
    do n = 1, atoms%ntype
      nd = atoms%ntypsy(nat)
      imax = atoms%jri(n)
      lmax = sphhar%llh(sphhar%nlh(nd), nd)
      if ( potdenType == POTDEN_TYPE_POTYUK ) then
        !do concurrent (i = 1:imax)
        do i = 1,imax
          call ModSphBessel( il(0:,i), kl(0:,i), input%preconditioning_param * atoms%rmsh(i,n), lmax )
        end do
      end if
      do lh = 0, sphhar%nlh(nd)
        l = sphhar%llh(lh,nd)
        if ( potdenType == POTDEN_TYPE_POTYUK ) then
          green_1(1:imax) = il(l,1:imax)
          green_2(1:imax) = kl(l,1:imax)
          green_factor    = fpi_const * input%preconditioning_param
        else
          green_1(1:imax) = atoms%rmsh(1:imax,n) ** l
          green_2(1:imax) = 1.0 / ( green_1(1:imax) * atoms%rmsh(1:imax,n) )
          green_factor    = fpi_const / ( 2 * l + 1 )
        end if
        integrand_1(1:imax) = green_1(1:imax) * rho(1:imax,lh,n)
        integrand_2(1:imax) = green_2(1:imax) * rho(1:imax,lh,n)
        call intgr2( integrand_1(1:imax), atoms%rmsh(1,n), atoms%dx(n), imax, integral_1(1:imax) )
        call intgr2( integrand_2(1:imax), atoms%rmsh(1,n), atoms%dx(n), imax, integral_2(1:imax) )
        termsR = integral_2(imax) + ( vtl(lh,n) / green_factor - integral_1(imax) * green_2(imax) ) / green_1(imax)
        vr(1:imax,lh,n) = green_factor * (   green_1(1:imax) * ( termsR - integral_2(1:imax) ) &
                                           + green_2(1:imax) *            integral_1(1:imax)   )
      end do
      nat = nat + atoms%neq(n)
    end do
    if ( potdenType == POTDEN_TYPE_POTYUK ) then
      deallocate( il, kl )
    end if

    

    if ( potdenType /= POTDEN_TYPE_POTYUK ) then
      do n = 1, atoms%ntype
        vr(1:atoms%jri(n),0,n) = vr(1:atoms%jri(n),0,n) - sfp_const * ( 1.0 / atoms%rmsh(1:atoms%jri(n),n) - 1.0 / atoms%rmt(n) ) * atoms%zatom(n)
      end do
    end if
    end if

  end subroutine vmts

end module m_vmts
