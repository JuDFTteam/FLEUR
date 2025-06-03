module m_vmts
#ifdef CPP_MPI
  use mpi
#endif
contains

  subroutine vmts( input, fmpi, stars, sphhar, atoms, sym, cell, juphon, dosf, vpw, rho, potdenType, vr, rhoIm, vrIm, iDtype, iDir, iDir2, mat2ord )

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

    use m_constants
    use m_types
    use m_mpi_reduce_tool
    use m_intgr, only : intgr2, sfint
    use m_phasy1
    use m_sphbes
     
    use m_SphBessel
    !$ use omp_lib
    implicit none

    type(t_input),  intent(in)        :: input
    type(t_mpi),    intent(in)        :: fmpi
    type(t_stars),  intent(in)        :: stars
    type(t_sphhar), intent(in)        :: sphhar
    type(t_atoms),  intent(in)        :: atoms
    type(t_sym),    intent(in)        :: sym
    type(t_cell),   intent(in)        :: cell
    type(t_juphon), intent(in)        :: juphon
     
    LOGICAL,        INTENT(IN)        :: dosf
    complex,        intent(in)        :: vpw(:)!(stars%ng3,input%jspins)
    real,           intent(in)        :: rho(:,0:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    integer,        intent(in)        :: potdenType
    real,           intent(out)       :: vr(:,0:,:)!(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    REAL,    OPTIONAL, INTENT(IN)     :: rhoIm(:,0:,:)
    REAL,    OPTIONAL, INTENT(OUT)    :: vrIm(:,0:,:)
    INTEGER, OPTIONAL, INTENT(IN)     :: iDtype, iDir
    INTEGER, OPTIONAL, INTENT(IN)     :: iDir2
    COMPLEX, OPTIONAL, INTENT(IN)     :: mat2ord(5,3,3)

    complex                           :: cp, sm
    integer                           :: i, jm, k, l, lh, n, nd, lm, m, imax, lmax, iMem, ptsym
    integer                           :: maxBunchSize, numBunches, iBunch, firstStar, iTempArray
    complex                           :: temp
    complex                           :: vtl(0:sphhar%nlhd, atoms%ntype)
    complex                           :: pylm(( atoms%lmaxd + 1 ) ** 2, atoms%ntype)
    real                              :: green_factor, termsR, pref
    real                              :: green_1    (1:atoms%jmtd), green_2    (1:atoms%jmtd)
    real                              :: integrand_1(1:atoms%jmtd), integrand_2(1:atoms%jmtd)
    real                              :: integral_1 (1:atoms%jmtd), integral_2 (1:atoms%jmtd)!, integral_3 (1:atoms%jmtd)
    real                              :: sbf(0:atoms%lmaxd)
    real, allocatable, dimension(:,:) :: il, kl
    LOGICAL                           :: l_dfptvgen
    COMPLEX, ALLOCATABLE              :: vtlStars(:,:,:), vtlLocal(:,:)
    TYPE(t_parallelLoop)              :: mpiLoop, ompLoop

    l_dfptvgen = PRESENT(rhoIm)

    ! SPHERE BOUNDARY CONTRIBUTION to the coefficients calculated from the values
    ! of the interstitial Coulomb / Yukawa potential on the sphere boundary

    ALLOCATE (vtlLocal(0:sphhar%nlhd,atoms%ntype))
    vtlLocal(:,:) = cmplx(0.0,0.0)

    firstStar = MERGE(2,1,norm2(stars%center)<=1e-8)
    maxBunchSize = 2*getNumberOfThreads() ! This bunch size is kind of a magic number determined from some
                                          ! naive performance tests for a 64 atom unit cell
    CALL calcNumberComputationBunches(firstStar, stars%ng3, maxBunchSize, numBunches)

    CALL mpiLoop%init(fmpi%irank,fmpi%isize,0,numBunches-1)

    ALLOCATE(vtlStars(0:sphhar%nlhd,atoms%ntype,maxBunchSize))
    vtlStars = CMPLX(0.0,0.0)

    ! q/=0 components
    DO iBunch = mpiLoop%bunchMinIndex, mpiLoop%bunchMaxIndex
       CALL ompLoop%init(iBunch,numBunches,firstStar,stars%ng3)
       !$OMP parallel do default( NONE ) &
       !$OMP SHARED(ompLoop,atoms,stars,sym,cell,sphhar,vpw,vtlStars,potdenType) &
       !$OMP private(iTempArray,cp,pylm,n,sbf,nd,lh,l,sm,jm,m,lm)
       do k = ompLoop%bunchMinIndex, ompLoop%bunchMaxIndex
          iTempArray = k - ompLoop%bunchMinIndex + 1
          cp = vpw(k) * stars%nstr(k)
          call phasy1( atoms, stars, sym, cell, k, pylm )
          do n = 1, atoms%ntype
             call sphbes( atoms%lmax(n), stars%sk3(k) * atoms%rmt(n), sbf )
             nd = sym%ntypsy(atoms%firstAtom(n))
             do lh = 0, sphhar%nlh(nd)
                l = sphhar%llh(lh,nd)
                sm = (0.0,0.0)
                do jm = 1, sphhar%nmem(lh,nd)
                   m = sphhar%mlh(jm,lh,nd)
                   lm = l * ( l + 1 ) + m + 1
                   sm = sm + conjg( sphhar%clnu(jm,lh,nd) ) * pylm(lm,n)
                end do
                vtlStars(lh,n,iTempArray) = vtlStars(lh,n,iTempArray) + cp * sbf(l) * sm
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    END DO
    
    !$OMP parallel do default( NONE ) &
    !$OMP SHARED(atoms,sym,sphhar,maxBunchSize,vtlStars,vtlLocal) &
    !$OMP private(nd,lh,temp,iTempArray)
    DO n = 1, atoms%ntype
       nd = sym%ntypsy(atoms%firstAtom(n))
       do lh = 0, sphhar%nlh(nd)
          temp = CMPLX(0.0,0.0)
          DO iTempArray = 1, maxBunchSize
             temp = temp + vtlStars(lh,n,iTempArray)
          END DO
          vtlLocal(lh,n) = temp
       END DO
    END DO
    !$OMP END PARALLEL DO
    
    ! q=0 component
    if ( fmpi%irank == 0 .AND. norm2(stars%center)<=1e-8 ) then
       DO n = 1, atoms%ntype
          vtlLocal(0,n) = vtlLocal(0,n) + sfp_const * vpw(1)
       END DO
    end if

    CALL mpi_sum_reduce(vtlLocal, vtl,fmpi%mpi_comm)

    ! SPHERE INTERIOR CONTRIBUTION to the coefficients calculated from the
    ! values of the sphere Coulomb/Yukawa potential on the sphere boundary

    if( fmpi%irank == 0 ) then
    if ( potdenType == POTDEN_TYPE_POTYUK ) then
      allocate( il(0:atoms%lmaxd, 1:atoms%jmtd), kl(0:atoms%lmaxd, 1:atoms%jmtd) )
    end if

    do n = 1, atoms%ntype
      nd = sym%ntypsy(atoms%firstAtom(n))
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
        if (.not.dosf) THEN
         call intgr2( integrand_1(1:imax), atoms%rmsh(1,n), atoms%dx(n), imax, integral_1(1:imax) )
         call intgr2( integrand_2(1:imax), atoms%rmsh(1,n), atoms%dx(n), imax, integral_2(1:imax) )
        else
           call sfint(integrand_1(1:imax),atoms%rmsh(:,n),atoms%dx(n),imax,integral_1(1:imax))
           call sfint(integrand_2(1:imax),atoms%rmsh(:,n),atoms%dx(n),imax,integral_2(1:imax))
        end if
        termsR = integral_2(imax) + ( vtl(lh,n) / green_factor - integral_1(imax) * green_2(imax) ) / green_1(imax)
        vr(1:imax,lh,n) = green_factor * (   green_1(1:imax) * ( termsR - integral_2(1:imax) ) &
                                           + green_2(1:imax) *            integral_1(1:imax)   )
        IF (l_dfptvgen) THEN
           ! Integrate the imaginary part of the density perturbation as well.
           integrand_1(1:imax) = green_1(1:imax) * rhoIm(1:imax,lh,n)
           integrand_2(1:imax) = green_2(1:imax) * rhoIm(1:imax,lh,n)
           call intgr2( integrand_1(1:imax), atoms%rmsh(1,n), atoms%dx(n), imax, integral_1(1:imax) )
           call intgr2( integrand_2(1:imax), atoms%rmsh(1,n), atoms%dx(n), imax, integral_2(1:imax) )
           termsR = integral_2(imax) + ( AIMAG(vtl(lh,n)) / green_factor - integral_1(imax) * green_2(imax) ) / green_1(imax)
           vrIm(1:imax,lh,n) = green_factor * (   green_1(1:imax) * ( termsR - integral_2(1:imax) ) &
                                                + green_2(1:imax) *            integral_1(1:imax)   )
        END IF
      end do
    end do
    if ( potdenType == POTDEN_TYPE_POTYUK ) then
      deallocate( il, kl )
    end if

    if ( potdenType /= POTDEN_TYPE_POTYUK .AND. potdenType /= POTDEN_TYPE_CRYSTALFIELD) then
      IF (.NOT.l_dfptvgen) THEN
         do n = 1, atoms%ntype
         vr(1:atoms%jri(n),0,n) = vr(1:atoms%jri(n),0,n) - sfp_const * ( 1.0 / atoms%rmsh(1:atoms%jri(n),n) - 1.0 / atoms%rmt(n) ) * atoms%zatom(n)
         end do
      ELSE IF (juphon%l_phonon) THEN
         IF (.NOT.PRESENT(iDir2)) THEN
            ! DFPT(-phonon) case:
            ! l=1 contributions from the Coulomb singularity instead of l=0 (1/r -> 1/r^2)
            DO n = MERGE(1,iDtype,iDtype==0), MERGE(atoms%ntype,iDtype,iDtype==0)
               ptsym = sym%ntypsy(atoms%firstAtom(n))
               pref = MERGE(atoms%zatom(n),-atoms%zatom(n),iDtype==0)
               DO lh = 1, 3
                  l = sphhar%llh(lh, ptsym)
                  DO iMem = 1, sphhar%nmem(lh, ptsym)
                     m = sphhar%mlh(iMem, lh, ptsym)
                     lm = l*(l+1) + m + 1
                     vr(1:atoms%jri(n),lh,n) = vr(1:atoms%jri(n),lh,n) + &
                                                conjg(sphhar%clnu(iMem, lh, ptsym)) * c_im(iDir, lm - 1) * pref * &
                                                ( 1 - (atoms%rmsh(1:atoms%jri(n), n) / atoms%rmt(n))**3) / atoms%rmsh(1:atoms%jri(n),n)**2
                  END DO
               END DO
            END DO
         ELSE
            ! DFPT 2nd order case:
            ! l=2 contributions from the Coulomb singularity instead of l=0 (1/r -> 1/r^3)
            DO n = 1, atoms%ntype!MERGE(1,iDtype,iDtype==0), MERGE(atoms%ntype,iDtype,iDtype==0)
               ptsym = sym%ntypsy(atoms%firstAtom(n))
               pref = -atoms%zatom(n)!MERGE(atoms%zatom(n),-atoms%zatom(n),iDtype==0)
               DO lh = 4, 8
                  l = sphhar%llh(lh, ptsym)
                  DO iMem = 1, sphhar%nmem(lh, ptsym)
                     m = sphhar%mlh(iMem, lh, ptsym)
                     lm = l*(l+1) + m + 1
                     IF ((n.EQ.iDtype).OR.(0.EQ.iDtype)) vr(1:atoms%jri(n),lh,n) = vr(1:atoms%jri(n),lh,n) + &
                                                         conjg(sphhar%clnu(iMem, lh, ptsym)) * mat2ord(lm-4,iDir2,iDir) * pref * &
                                                         ( 1 - (atoms%rmsh(1:atoms%jri(n), n) / atoms%rmt(n))**5) / atoms%rmsh(1:atoms%jri(n),n)**3
                  END DO
               END DO
            END DO
         END IF
      END IF
    end if
    end if

  end subroutine vmts

end module m_vmts
