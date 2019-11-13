!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_vgen_coulomb

  use m_juDFT

contains

  subroutine vgen_coulomb( ispin, mpi, dimension, oneD, input, field, vacuum, sym, stars, &
             cell, sphhar, atoms, dosf, den, vCoul, results )
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
    use m_od_vvac
    use m_od_vvacis
    use m_convol
    use m_psqpw
    use m_cfft
    implicit none

    integer,            intent(in)               :: ispin
    type(t_mpi),        intent(in)               :: mpi
    type(t_dimension),  intent(in)               :: dimension
    type(t_oneD),       intent(in)               :: oneD
    type(t_input),      intent(in)               :: input
    type(t_field),      intent(inout)            :: field
    type(t_vacuum),     intent(in)               :: vacuum
    type(t_sym),        intent(in)               :: sym
    type(t_stars),      intent(in)               :: stars
    type(t_cell),       intent(in)               :: cell
    type(t_sphhar),     intent(in)               :: sphhar
    type(t_atoms),      intent(in)               :: atoms 
    LOGICAL,            INTENT(IN)               :: dosf
    type(t_potden),     intent(in)               :: den
    type(t_potden),     intent(inout)            :: vCoul
    type(t_results),    intent(inout), optional  :: results

    complex                                      :: vintcza, xint, rhobar
    integer                                      :: i, i3, irec2, irec3, ivac, j, js, k, k3
    integer                                      :: lh, n, nzst1
    integer                                      :: imz, imzxy, ichsmrg, ivfft
    integer                                      :: l, nat
    real                                         :: ani, g3, z, sig1dh, vz1dh
    complex, allocatable                         :: alphm(:,:), psq(:)
    real,    allocatable                         :: af1(:), bf1(:)
#ifdef CPP_MPI
    include 'mpif.h'
    integer:: ierr
#endif
 
    
    allocate ( alphm(stars%ng2,2), af1(3*stars%mx3), bf1(3*stars%mx3), psq(stars%ng3)  )
    vCoul%iter = den%iter



    ! PSEUDO-CHARGE DENSITY COEFFICIENTS
    call timestart( "psqpw" )      
    call psqpw( mpi, atoms, sphhar, stars, vacuum,  cell, input, sym, oneD, &
         den%pw(:,ispin), den%mt(:,:,:,ispin), den%vacz(:,:,ispin), .false., vCoul%potdenType, psq )
    call timestop( "psqpw" )



    ! VACUUM POTENTIAL
    if ( mpi%irank == 0 ) then
      if ( oneD%odi%d1 ) then
        call timestart( "Vacuum" )
        !---> generates the m=0,gz=0 component of the vacuum potential
        call od_vvac( stars, vacuum, cell, psq, den%vacz(:,:,ispin), vCoul%vacz(:,:,ispin) )
        !---> generation of the vacuum warped potential components and
        !---> interstitial pw potential
        !---> vvacxy_5.F is a symmetrized potential generator
        call od_vvacis( oneD%odi%n2d, input, vacuum, oneD%odi%nq2, &
             oneD%odi%kv, cell, oneD%odi%M, stars, oneD%odi%nst2, &
             oneD, den%vacz(:,:,ispin), den%vacxy(:,:,:,ispin), psq, &
             vCoul%vacz(:,:,ispin), sym, vCoul%vacxy(:,:,:,ispin), vCoul%pw(:,ispin) )
        call timestop( "Vacuum" )
        !---> generation of the vacuum warped potential components and       
      elseif ( input%film .and. .not. oneD%odi%d1 ) then
        !     ----> potential in the  vacuum  region
        call timestart( "Vacuum" ) 
        call vvac( vacuum, stars, cell, sym, input, field, psq, den%vacz(:,:,ispin), vCoul%vacz(:,:,ispin), rhobar, sig1dh, vz1dh )
        call vvacis( stars, vacuum, sym, cell, psq, input, field, vCoul%vacxy(:,:,:,ispin) )
        call vvacxy( stars, vacuum, cell, sym, input, field, den%vacxy(:,:,:,ispin), vCoul%vacxy(:,:,:,ispin), alphm )
        call timestop( "Vacuum" )
      end if



      ! INTERSTITIAL POTENTIAL
      call timestart( "interstitial" )
      write( 6, fmt=8010 )
8010  format (/,5x,'coulomb potential in the interstitial region:')
      ! in case of a film:
      if ( input%film .and. .not.oneD%odi%d1 ) then
        ! create v(z) for each 2-d reciprocal vector
        ivfft = 3 * stars%mx3 
        ani = 1.0 / real( ivfft )
        do irec2 = 1, stars%ng2
          i = 0
          do i3 = 0, ivfft - 1
            i = i + 1
            z = cell%amat(3,3) * i3 * ani
            if ( z > cell%amat(3,3) / 2. ) z = z - cell%amat(3,3)
            vintcza = vintcz( stars, vacuum, cell, sym, input, field, z, irec2, psq, &
                              vCoul%vacxy(:,:,:,ispin), vCoul%vacz(:,:,ispin), &
                              rhobar, sig1dh, vz1dh, alphm )
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
      ! in case of a bulk system:
      else if ( .not. input%film ) then
        if ( vCoul%potdenType == POTDEN_TYPE_POTYUK ) then
          vCoul%pw(1:stars%ng3,ispin) = fpi_const * psq(1:stars%ng3) &
            / ( stars%sk3(1:stars%ng3) ** 2 + input%preconditioning_param ** 2 )
          ! if( abs( real( psq(1) ) ) * cell%omtil < 0.01 ) vCoul%pw(1,ispin) = 0.0
          ! there is a better option now using qfix in mix
        else
          vCoul%pw(1,ispin) = cmplx(0.0,0.0) 
          vCoul%pw(2:stars%ng3,ispin) = fpi_const * psq(2:stars%ng3) / stars%sk3(2:stars%ng3) ** 2
        end if
      end if
    call timestop("interstitial")
    end if ! mpi%irank == 0

    if (dosf) then
       vCoul%pw(:,:)=0.0
    end do

    ! MUFFIN-TIN POTENTIAL
    call timestart( "MT-spheres" )
#ifdef CPP_MPI
    CALL MPI_BARRIER(mpi%mpi_comm,ierr) !should be totally useless, but needed anyway????
    call MPI_BCAST( vcoul%pw, size(vcoul%pw), MPI_DOUBLE_COMPLEX, 0, mpi%mpi_comm, ierr )
    CALL MPI_BARRIER(mpi%mpi_comm,ierr) !should be totally useless, but ...
#endif
    call vmts( input, mpi, stars, sphhar, atoms, sym, cell, oneD, vCoul%pw(:,ispin), &
               den%mt(:,0:,:,ispin), vCoul%potdenType, vCoul%mt(:,0:,:,ispin) )
    call timestop( "MT-spheres" )



    if( vCoul%potdenType == POTDEN_TYPE_POTYUK ) return



    if ( mpi%irank == 0 ) then
      CHECK_CONTINUITY: if ( input%vchk ) then
        call timestart( "checking" )
        call checkDOPAll( input, dimension, sphhar, stars, atoms, sym, vacuum, oneD, &
                          cell, vCoul, ispin )
        call timestop( "checking" )
      end if CHECK_CONTINUITY

      CALCULATE_DENSITY_POTENTIAL_INTEGRAL: if ( present( results ) ) then
        if ( input%total ) then
          call timestart( "den-pot integrals" )
          !     CALCULATE THE INTEGRAL OF n*Vcoulomb
          write( 6, fmt=8020 )
8020      format (/,10x,'density-coulomb potential integrals',/)
          !       interstitial first
          !       convolute ufft and pot: F(G) = \sum_(G') U(G - G') V(G')
          call convol( stars, vCoul%pw_w(:,ispin), vCoul%pw(:,ispin), stars%ufft )
          results%te_vcoul = 0.0
          call int_nv( ispin, stars, vacuum, atoms, sphhar, cell, sym, input, oneD, &
                       vCoul, den, results%te_vcoul )

          write( 6, fmt=8030 ) results%te_vcoul
8030      format (/,10x,'total density-coulomb potential integral :', t40,f20.10)

          call timestop( "den-pot integrals" )
        end if
      end if CALCULATE_DENSITY_POTENTIAL_INTEGRAL
    end if !irank==0

  end subroutine vgen_coulomb

end module m_vgen_coulomb
