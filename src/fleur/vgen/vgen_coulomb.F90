!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_vgen_coulomb

  use m_juDFT
#ifdef CPP_MPI
  use mpi
#endif
contains

  subroutine vgen_coulomb( ispin, fmpi,    input, field, vacuum, sym, juphon, stars, &
             cell, sphhar, atoms, dosf, den, vCoul, sigma_disc, results, dfptdenimag, dfptvCoulimag, dfptden0, stars2, iDtype, iDir, iDir2, sigma_disc2 )
    !----------------------------------------------------------------------------
    ! FLAPW potential generator
    !----------------------------------------------------------------------------
    ! Generates the Coulomb or Yukawa potential and optionally the
    ! density-potential integrals
    ! vCoul%potdenType = POTDEN_TYPE_POTYUK -> Yukawa case
    ! It takes a spin variable to indicate in which spin-channel the charge
    ! resides.
    !----------------------------------------------------------------------------

    use m_constants
    use m_types
    use m_vmts
    use m_intnv
    use m_vvac
    use m_vvacis
    use m_vvacxy
    use m_vintcz
    use m_checkdopall
    use m_convol
    use m_psqpw
    use m_cfft
    use m_dfpt_vefield_mt
    
    implicit none

    integer,            intent(in)               :: ispin
    type(t_mpi),        intent(in)               :: fmpi

     
    type(t_input),      intent(in)               :: input
    type(t_field),      intent(in)               :: field
    type(t_vacuum),     intent(in)               :: vacuum
    type(t_sym),        intent(in)               :: sym
    type(t_juphon),     intent(in)               :: juphon
    type(t_stars),      intent(in)               :: stars
    type(t_cell),       intent(in)               :: cell
    type(t_sphhar),     intent(in)               :: sphhar
    type(t_atoms),      intent(in)               :: atoms
    LOGICAL,            INTENT(IN)               :: dosf
    type(t_potden),     intent(in)               :: den
    type(t_potden),     intent(inout)            :: vCoul
    COMPLEX,            INTENT(INOUT)            :: sigma_disc(2)
    type(t_results),    intent(inout), optional  :: results

    TYPE(t_potden),     OPTIONAL, INTENT(IN)     :: dfptdenimag,  dfptden0
    TYPE(t_potden),     OPTIONAL, INTENT(INOUT)  :: dfptvCoulimag
    TYPE(t_stars),      OPTIONAL, INTENT(IN)     :: stars2
    INTEGER, OPTIONAL, INTENT(IN)                :: iDtype, iDir ! DFPT: Type and direction of displaced atom
    INTEGER, OPTIONAL, INTENT(IN)                :: iDir2 ! DFPT: 2nd direction for 2nd order VC
    COMPLEX, OPTIONAL, INTENT(IN)                :: sigma_disc2(2)

    complex                                      :: vintcza, xint, rhobar,vslope
    integer                                      :: i, i3, irec2, irec3, ivac, j, js, k, k3
    integer                                      :: lh, n, nzst1, first_star
    integer                                      :: imz, imzxy, ichsmrg, ivfft
    integer                                      :: l, nat
    real                                         :: ani, g3, z
    complex                                      :: sig1dh, vz1dh, vmz1dh, vmz1dh_is
    complex                                      :: mat2ord(5,3,3), sigma_loc(2), sigma_loc2(2)
    complex, allocatable                         :: alphm(:,:), psq(:)
    real,    allocatable                         :: af1(:), bf1(:)
    LOGICAL :: l_dfptvgen ! If this is true, we handle things differently!
    LOGICAL :: l_2ndord, l_corr

    ! For non analytic corrections:
    !REAL :: q_lim
#ifdef CPP_MPI
    integer:: ierr
#endif
    l_dfptvgen = PRESENT(stars2)
    l_2ndord = PRESENT(iDir2)
    l_corr = ALL(ABS(den%vac)<1e-12)
    vmz1dh_is = cmplx(0.0,0.0)
    sigma_loc = sigma_disc
    sigma_loc2 = MERGE(sigma_disc,cmplx(0.0,0.0),PRESENT(sigma_disc2))

    vintcza = cmplx(0.0,0.0)
    sig1dh = cmplx(0.0,0.0)
    vz1dh = cmplx(0.0,0.0)
    vmz1dh = cmplx(0.0,0.0)
    vslope = cmplx(0.0,0.0)

    allocate ( alphm(stars%ng2,2), af1(3*stars%mx3), bf1(3*stars%mx3), psq(stars%ng3)  )
    vCoul%iter = den%iter


    ! PSEUDO-CHARGE DENSITY COEFFICIENTS
    call timestart( "psqpw" )
    if (.not.l_dfptvgen) then
        call psqpw( fmpi, atoms, sphhar, stars, vacuum,  cell, input, sym,   &
            & juphon, den, ispin, .false., vCoul%potdenType, psq, sigma_loc )
    else if (.not.l_2ndord) then
        ! If we do DFPT, the MT density perturbation has an imaginary part that needs to be explicitly carried
        ! as another variable dfptdenimag%mt and results in the same component for the Coulomb potential later on.
        ! Also, the ionic qlm behave differently.
        call psqpw( fmpi, atoms, sphhar, stars, vacuum,  cell, input, sym,   &
            & juphon, den, ispin, .false., vCoul%potdenType, psq, sigma_loc,&
            & dfptdenimag%mt(:,:,:,ispin), stars2, iDtype, iDir, dfptden0%mt(:,:,:,ispin), dfptden0%pw(:,ispin) )
    else
        call make_mat_2nd(mat2ord)
        call psqpw( fmpi, atoms, sphhar, stars, vacuum,  cell, input, sym,   &
            & juphon, den, ispin, .false., vCoul%potdenType, psq, sigma_loc,&
            & dfptdenimag%mt(:,:,:,ispin), stars2, iDtype, iDir, dfptden0%mt(:,:,:,ispin), dfptden0%pw(:,ispin), iDir2, mat2ord )
    end if
    call timestop( "psqpw" )

    ! VACUUM POTENTIAL
    if ( fmpi%irank == 0 ) then
      if ( input%film ) then
        !     ----> potential in the  vacuum  region
        call timestart( "Vacuum" )
        if ((.not.l_dfptvgen).or.norm2(stars%center)<1e-8) then
          ! If we do DPFT AND q/=0, there is no G_||+q_||=0 part! So all components are
          ! handled by the G_||/=0 parts in vvacis/vvacxy, that are told to explicitly
          ! start at star 1 instead of 2 for this!
          if (.not.l_2ndord) then
            !call vvac( vacuum, stars, cell,  input, field, psq, den%vac(:,1,:,ispin), vCoul%vac(:,1,:,ispin), rhobar, sig1dh, vz1dh,vslope,.FALSE..AND.l_dfptvgen,vmz1dh,sigma_disc )
            call vvac( vacuum, stars, cell,  input, field, psq, den%vac(:,1,:,ispin), vCoul%vac(:,1,:,ispin), rhobar, sig1dh, vz1dh,vslope,l_corr,vmz1dh,sigma_disc )
          else
            !call vvac( vacuum, stars, cell,  input, field, psq, den%vac(:,1,:,ispin), vCoul%vac(:,1,:,ispin), rhobar, sig1dh, vz1dh,vslope,.TRUE.,vmz1dh,sigma_disc,sigma_disc2 )
            call vvac( vacuum, stars, cell,  input, field, psq, den%vac(:,1,:,ispin), vCoul%vac(:,1,:,ispin), rhobar, sig1dh, vz1dh,vslope,l_corr,vmz1dh,sigma_disc,sigma_disc2 )
          end if
        end if
        call vvacis( stars, vacuum, cell, psq, input, field, vCoul%vac(:vacuum%nmzxyd,:,:,ispin), l_dfptvgen, l_corr )
        call vvacxy( stars, vacuum, cell, sym, input, field, den%vac(:vacuum%nmzxyd,:,:,ispin), vCoul%vac(:vacuum%nmzxyd,:,:,ispin), alphm, l_dfptvgen )
        call timestop( "Vacuum" )
      end if

      ! INTERSTITIAL POTENTIAL
      call timestart( "interstitial" )
      write(oUnit, fmt=8010 )
8010  format (/,5x,'coulomb potential in the interstitial region:')
      ! in case of a film:
      if ( input%film ) then
        ! create v(z) for each 2-d reciprocal vector
        ivfft = 3 * stars%mx3
        ani = 1.0 / real( ivfft )
        do irec2 = 1, stars%ng2
          ! If we do DFPT, we want to fix the second vacuum at infinity to 0. This is WIP,
          ! as to how we want to do it eventually. Here, we calculate the necessary offset.
          !IF (l_dfptvgen.AND.irec2 == 1) vmz1dh_is = vintcz( stars, vacuum, cell,  input, field, -cell%z1+1e-15, irec2, psq, &
          !                    vCoul%vac(:,:,:,ispin), &
          !                    rhobar, sig1dh, vz1dh, alphm, vslope, .TRUE., CMPLX(0.0,0.0) )
          i = 0
          do i3 = 0, ivfft - 1
            i = i + 1
            z = cell%amat(3,3) * i3 * ani
            if ( z > cell%amat(3,3) / 2. ) z = z - cell%amat(3,3)
            if (.not.l_2ndord) then
              vintcza = vintcz( stars, vacuum, cell,  input, field, z, irec2, psq, &
                                vCoul%vac(:,:,:,ispin), &
                                rhobar, sig1dh, vz1dh, alphm, vslope, sigma_disc, l_dfptvgen, l_corr, vmz1dh-vmz1dh_is )
            else
            vintcza = vintcz( stars, vacuum, cell,  input, field, z, irec2, psq, &
                              vCoul%vac(:,:,:,ispin), &
                                rhobar, sig1dh, vz1dh, alphm, vslope, sigma_disc, l_dfptvgen, l_corr, vmz1dh-vmz1dh_is, sigma_disc2 )
            end if
            af1(i) = real( vintcza )
            bf1(i) = aimag( vintcza )
          end do
          !                z = (i_sm-1)*ani
          !                IF (z > 0.5) z = z - 1.0
          !                af1(i_sm) = af1(i_sm) + z * delta
          !                bf1(i_sm) = bf1(i_sm) + z * deltb
          !              ENDDO
          !            ENDIF


          !        --> 1-d fourier transform and store the coefficients in vTot%pw( ,1)
          call cfft( af1, bf1, ivfft, ivfft, ivfft, -1 )
          !            delta = ivfft * delta * 2 / fpi ! * amat(3,3)**2 * ani
          i = 0
          do i3 = 0, ivfft - 1
            k3 = i3
            if ( k3 > floor( ivfft / 2. ) ) k3 = k3 - ivfft
            i = i + 1
            if ( ( k3 >= -stars%mx3 ) .and. ( k3 <= stars%mx3 ) ) then
              irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),k3)
              !                 IF ( (irec2 == 1).AND.(i3 > 0) ) THEN                 ! smooth potential
              !                   corr = 2.0*mod(abs(k3),2) - 1.0
              !                   bf1(i) = bf1(i) + delta * corr / k3
              !                 ENDIF
              !       ----> only stars within g_max sphere (shz oct.97)
              if ( irec3 /= 0 ) then
                xint = cmplx( af1(i), bf1(i) ) * ani
                nzst1 = stars%nstr(irec3) / stars%nstr2(irec2)
                vCoul%pw(irec3,1) = vCoul%pw(irec3,1) + xint / nzst1
              end if
            end if
          end do
        end do
        sigma_disc = sigma_loc
      ! in case of a bulk system:
      else if ( .not. input%film ) then
        if ( vCoul%potdenType == POTDEN_TYPE_POTYUK ) then
          vCoul%pw(1:stars%ng3,ispin) = fpi_const * psq(1:stars%ng3) &
            / ( stars%sk3(1:stars%ng3) ** 2 + input%preconditioning_param ** 2 )
          ! if( abs( real( psq(1) ) ) * cell%omtil < 0.01 ) vCoul%pw(1,ispin) = 0.0
          ! there is a better option now using qfix in mix
        else
          vCoul%pw(1,ispin) = cmplx(0.0,0.0)
          first_star = MERGE(2,1,stars%sk3(1)< 1E-9)
          vCoul%pw(first_star:stars%ng3,ispin) = fpi_const * psq(first_star:stars%ng3) / stars%sk3(first_star:stars%ng3) ** 2
        end if
      end if
    call timestop("interstitial")
    end if ! fmpi%irank == 0

    ! MUFFIN-TIN POTENTIAL
    !print*,"dfptvCoulimag%mt(:,0:,:,ispin)"
    !print*, dfptvCoulimag%mt(:,0:,:,ispin)
    call timestart( "MT-spheres" )
#ifdef CPP_MPI
    CALL MPI_BARRIER(fmpi%mpi_comm,ierr) !should be totally useless, but needed anyway????
    call MPI_BCAST( vcoul%pw, size(vcoul%pw), MPI_DOUBLE_COMPLEX, 0, fmpi%mpi_comm, ierr )
    CALL MPI_BARRIER(fmpi%mpi_comm,ierr) !should be totally useless, but ...
#endif

    IF (.NOT.l_dfptvgen) THEN
      call vmts( input, fmpi, stars, sphhar, atoms, sym, cell, juphon,  dosf, vCoul%pw(:,ispin), &
                 den%mt(:,0:,:,ispin), vCoul%potdenType, vCoul%mt(:,0:,:,ispin) )
    ELSE IF (.NOT.l_2ndord) THEN
      !write(1500,*) "mache normal dfpt" ! fort.1500 
      ! For DFPT there is a) an imaginary part to the potential and b) a different treatment
      ! for the ionic 1/r (now 1/r^2) contribution.
      call vmts( input, fmpi, stars, sphhar, atoms, sym, cell, juphon,  dosf, vCoul%pw(:,ispin), &
                 den%mt(:,0:,:,ispin), vCoul%potdenType, vCoul%mt(:,0:,:,ispin), &
                 dfptdenimag%mt(:,0:,:,ispin), dfptvCoulimag%mt(:,0:,:,ispin), iDtype, iDir )
      !write(1000,*) vCoul%pw(:,ispin)
      !call save_npy()
    ELSE
      call vmts( input, fmpi, stars, sphhar, atoms, sym, cell, juphon,  dosf, vCoul%pw(:,ispin), &
                 den%mt(:,0:,:,ispin), vCoul%potdenType, vCoul%mt(:,0:,:,ispin), &
                 dfptdenimag%mt(:,0:,:,ispin), dfptvCoulimag%mt(:,0:,:,ispin), iDtype, iDir, iDir2, mat2ord )
    END IF
    call timestop( "MT-spheres" )

    if( vCoul%potdenType == POTDEN_TYPE_POTYUK ) return

    if ( fmpi%irank == 0 ) then
      CHECK_CONTINUITY: if ( input%vchk ) then ! TODO: We could use this for DFPT as well if we
        call timestart( "checking" )           !       passed an optional to checkDOPAll and modded
        call checkDOPAll( input,  sphhar, stars, atoms, sym, vacuum,   & ! slightly.
                          cell, vCoul, ispin )
        call timestop( "checking" )
      end if CHECK_CONTINUITY

      CALCULATE_DENSITY_POTENTIAL_INTEGRAL: if ( present( results ) ) then
          call timestart( "den-pot integrals" )
          !     CALCULATE THE INTEGRAL OF n*Vcoulomb
          write(oUnit, fmt=8020 )
8020      format (/,10x,'density-coulomb potential integrals',/)
          !       interstitial first
          !       convolute ufft and pot: F(G) = \sum_(G') U(G - G') V(G')
          call convol( stars, vCoul%pw_w(:,ispin), vCoul%pw(:,ispin))
          results%te_vcoul = 0.0
          call int_nv( ispin, stars, vacuum, atoms, sphhar, cell, sym, input,   &
                       vCoul, den, results%te_vcoul )

          write(oUnit, fmt=8030 ) results%te_vcoul
8030      format (/,10x,'total density-coulomb potential integral :', t40,f20.10)

          call timestop( "den-pot integrals" )
      end if CALCULATE_DENSITY_POTENTIAL_INTEGRAL
    end if !irank==0
  end subroutine vgen_coulomb

  subroutine make_mat_2nd(mat2ord)
     use m_constants
     complex, intent(out) :: mat2ord(5,3,3)
     mat2ord = cmplx(0.0,0.0)

     mat2ord(1,1,1) =  sqrt(3.0/2.0)
     mat2ord(1,1,2) =  sqrt(3.0/2.0)*ImagUnit
     mat2ord(1,2,1) =  sqrt(3.0/2.0)*ImagUnit
     mat2ord(1,2,2) = -sqrt(3.0/2.0)

     mat2ord(2,1,3) =  sqrt(3.0/2.0)
     mat2ord(2,2,3) =  sqrt(3.0/2.0)*ImagUnit   
     mat2ord(2,3,1) =  sqrt(3.0/2.0)
     mat2ord(2,3,2) =  sqrt(3.0/2.0)*ImagUnit  

     mat2ord(3,1,1) = -1.0! - 0.5 ! TODO: What the hell is this value???
     mat2ord(3,2,2) = -1.0! - 0.5 ! TODO: What the hell is this value???
     mat2ord(3,3,3) =  2.0! - 0.5 ! TODO: What the hell is this value???

     mat2ord(4,1,3) = -sqrt(3.0/2.0)
     mat2ord(4,2,3) =  sqrt(3.0/2.0)*ImagUnit   
     mat2ord(4,3,1) = -sqrt(3.0/2.0)
     mat2ord(4,3,2) =  sqrt(3.0/2.0)*ImagUnit  

     mat2ord(5,1,1) =  sqrt(3.0/2.0)
     mat2ord(5,1,2) = -sqrt(3.0/2.0)*ImagUnit
     mat2ord(5,2,1) = -sqrt(3.0/2.0)*ImagUnit
     mat2ord(5,2,2) = -sqrt(3.0/2.0)

     mat2ord = sqrt(4.0*pi_const/5.0) * mat2ord
  end subroutine

end module m_vgen_coulomb
