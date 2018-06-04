!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_vgen_coulomb

  USE m_juDFT

CONTAINS

  SUBROUTINE vgen_coulomb( ispin, mpi, DIMENSION, oneD, input, field, vacuum, sym, stars, &
             cell, sphhar, atoms, den, yukawa_residual, vCoul, results )
    !     ***********************************************************
    !     FLAPW potential generator                           
    !     ***********************************************************
    !     Generates the Coulomb potential and optionally the density-potential integrals
    !     is takes a spin variable to indicate in which spin-channel the charge resides...
    !     
    !     ***********************************************************
    USE m_constants
    USE m_vmts
    USE m_intnv
    USE m_vvac
    USE m_vvacis
    USE m_vvacxy
    USE m_vintcz
    USE m_checkdopall
    USE m_types
    USE m_od_vvac
    USE m_od_vvacis
    USE m_convol
    USE m_psqpw
    USE m_cfft
    IMPLICIT NONE

    INTEGER,            INTENT(IN)               :: ispin
    TYPE(t_mpi),        INTENT(IN)               :: mpi
    TYPE(t_dimension),  INTENT(IN)               :: dimension
    TYPE(t_oneD),       INTENT(IN)               :: oneD
    TYPE(t_input),      INTENT(IN)               :: input
    TYPE(t_field),      INTENT(INOUT)            :: field
    TYPE(t_vacuum),     INTENT(IN)               :: vacuum
    TYPE(t_sym),        INTENT(IN)               :: sym
    TYPE(t_stars),      INTENT(IN)               :: stars
    TYPE(t_cell),       INTENT(IN)               :: cell
    TYPE(t_sphhar),     INTENT(IN)               :: sphhar
    TYPE(t_atoms),      INTENT(IN)               :: atoms 
    TYPE(t_potden),     INTENT(IN)               :: den
    logical,            intent(in)               :: yukawa_residual
    TYPE(t_potden),     INTENT(INOUT)            :: vCoul
    TYPE(t_results),    INTENT(INOUT), OPTIONAL  :: results

    COMPLEX                                      :: vintcza, xint, rhobar
    INTEGER                                      :: i, i3, irec2, irec3, ivac, j, js, k, k3
    integer                                      :: lh, n, nzst1
    INTEGER                                      :: imz, imzxy, ichsmrg, ivfft
    INTEGER                                      :: l, nat
    REAL                                         :: ani, g3, z, sig1dh, vz1dh
    COMPLEX, ALLOCATABLE                         :: alphm(:,:), psq(:)
    REAL,    ALLOCATABLE                         :: af1(:), bf1(:)
#ifdef CPP_MPI
    include 'mpif.h'
    integer:: ierr
#endif

    ALLOCATE ( alphm(stars%ng2,2), af1(3*stars%mx3), bf1(3*stars%mx3), psq(stars%ng3)  )


    vCoul%iter = den%iter

    !     ************** coulomb potential ***********************
    !     ----> create pesudo-charge density coefficients


    CALL timestart( "psqpw" )      
    CALL psqpw( mpi, atoms, sphhar, stars, vacuum, DIMENSION, cell, input, sym, oneD, &
         den%pw(:,ispin), den%mt(:,:,:,ispin), den%vacz(:,:,ispin), .FALSE., yukawa_residual, psq )
    CALL timestop( "psqpw" )
    IF (mpi%irank == 0) THEN
       !First the vacuum
       !     ------------------------------------------
       IF (oneD%odi%d1) THEN
          !-odim
          CALL timestart( "Vacuum" )

          !---> generates the m=0,gz=0 component of the vacuum potential
          CALL od_vvac( stars, vacuum, cell, psq, den%vacz(:,:,ispin), vCoul%vacz(:,:,ispin) )

          !---> generation of the vacuum warped potential components and
          !---> interstitial pw potential
          !---> vvacxy_5.F is a symmetrized potential generator

          CALL od_vvacis( oneD%odi%n2d, dimension, vacuum, oneD%odi%nq2, &
               oneD%odi%kv, cell, oneD%odi%M, stars, oneD%odi%nst2, &
               oneD, den%vacz(:,:,ispin), den%vacxy(:,:,:,ispin), psq, &
               vCoul%vacz(:,:,ispin), sym, vCoul%vacxy(:,:,:,ispin), vCoul%pw(:,ispin) )
          CALL timestop( "Vacuum" )

          !+odim
       ELSEIF (input%film .AND. .NOT.oneD%odi%d1) THEN
          !     ----> potential in the  vacuum  region
          !       
          CALL timestart( "Vacuum" ) 
          CALL vvac( vacuum, stars, cell, sym, input, field, psq, den%vacz(:,:,ispin), vCoul%vacz(:,:,ispin), rhobar, sig1dh, vz1dh )
          CALL vvacis( stars, vacuum, sym, cell, psq, input, field, vCoul%vacxy(:,:,:,ispin) )

          CALL vvacxy( stars, vacuum, cell, sym, input, field, den%vacxy(:,:,:,ispin), vCoul%vacxy(:,:,:,ispin), alphm )
          CALL timestop( "Vacuum" )
       END IF
       !     ------------------------------------------
       !     ----> potential in the  interstitial  region
       CALL timestart( "interstitial" )
       WRITE (6,FMT=8010)
8010   FORMAT (/,5x,'coulomb potential in the interstitial region:')
       IF (input%film .AND. .NOT.oneD%odi%d1) THEN
          !           -----> create v(z) for each 2-d reciprocal vector
          ivfft =  3*stars%mx3 
          !         ivfft = 2*mx3 - 1
          ani = 1.0/REAL(ivfft)
          DO  irec2 = 1,stars%ng2
             i = 0
             DO i3 = 0,ivfft - 1
                i = i + 1
                z = cell%amat(3,3)*i3*ani
                IF (z.GT.cell%amat(3,3)/2.) z = z - cell%amat(3,3)
                vintcza = vintcz( stars, vacuum, cell, sym, input, field, &
                     z, irec2, psq, vCoul%vacxy(:,:,:,ispin), vCoul%vacz(:,:,ispin), &
                     rhobar, sig1dh, vz1dh, alphm )
                af1(i) = REAL( vintcza )
                bf1(i) = AIMAG( vintcza )
             ENDDO
             !                z = (i_sm-1)*ani
             !                IF (z > 0.5) z = z - 1.0
             !                af1(i_sm) = af1(i_sm) + z * delta
             !                bf1(i_sm) = bf1(i_sm) + z * deltb
             !              ENDDO
             !            ENDIF
             !        --> 1-d fourier transform and store the coefficients in vTot%pw( ,1)
             CALL cfft( af1, bf1, ivfft, ivfft, ivfft, -1 )
             !            delta = ivfft * delta * 2 / fpi ! * amat(3,3)**2 * ani
             i = 0
             DO  i3 = 0,ivfft - 1
                k3 = i3
                IF (k3 > FLOOR(ivfft/2.) ) k3 = k3 - ivfft
                i = i + 1
                IF ((k3.GE.-stars%mx3).AND.(k3.LE.stars%mx3)) THEN
                   irec3 = stars%ig(stars%kv2(1,irec2),stars%kv2(2,irec2),k3)

                   !                 IF ( (irec2 == 1).AND.(i3 > 0) ) THEN                 ! smooth potential
                   !                   corr = 2.0*mod(abs(k3),2) - 1.0
                   !                   bf1(i) = bf1(i) + delta * corr / k3
                   !                 ENDIF

                   !       ----> only stars within g_max sphere (shz oct.97)
                   IF (irec3.NE.0) THEN
                      !
                      xint = CMPLX(af1(i),bf1(i))*ani
                      nzst1 = stars%nstr(irec3)/stars%nstr2(irec2)
                      vCoul%pw(irec3,1) = vCoul%pw(irec3,1) + xint/nzst1
                   END IF
                ENDIF
             ENDDO
          ENDDO
       ELSEIF (.NOT.input%film) THEN
          vCoul%pw(1,ispin) = CMPLX(0.0,0.0)
          vCoul%pw(2:stars%ng3,ispin)=fpi_const*psq(2:stars%ng3)/(stars%sk3(2:stars%ng3)*stars%sk3(2:stars%ng3))       
       END IF

       CALL timestop("interstitial")

    ENDIF ! mpi%irank == 0
    !     --------------------------------------------
    !     ---> potential in the muffin-tin spheres

    CALL timestart( "MT-spheres" )
#ifdef CPP_MPI
    CALL MPI_BCAST( den%mt, atoms%jmtd*(1+sphhar%nlhd) * atoms%ntype * dimension%jspd, MPI_DOUBLE_PRECISION, 0, mpi%mpi_comm, ierr )
#endif
    CALL vmts( mpi, stars, sphhar, atoms, sym, cell, oneD, vCoul%pw(:,ispin), den%mt(:,0:,:,ispin), &
         yukawa_residual, vCoul%mt(:,0:,:,ispin) )
    CALL timestop( "MT-spheres" )

    IF (mpi%irank == 0) THEN
       !     ---> check continuity of coulomb potential
       IF (input%vchk) THEN
          CALL timestart( "checking" )
          CALL checkDOPAll( input, DIMENSION, sphhar, stars, atoms, sym, vacuum, oneD, &
               cell, vCoul, ispin )
          CALL timestop( "checking" )
       END IF

       !Calculate potential-density integrals if results-variable is present!
       IF (PRESENT(results)) THEN
          IF (input%total) THEN
             CALL timestart( "den-pot integrals" )
             !     CALCULATE THE INTEGRAL OF n*Vcoulomb
             WRITE (6,FMT=8020)
             WRITE (16,FMT=8020)
8020         FORMAT (/,10x,'density-coulomb potential integrals',/)
             !       interstitial first
             !       convolute ufft and pot: F(G) = \sum_(G') U(G - G') V(G')
             CALL convol( stars, vCoul%pw_w(:,ispin), vCoul%pw(:,ispin), stars%ufft )
             results%te_vcoul = 0.0
             CALL int_nv( ispin, stars, vacuum, atoms, sphhar, cell, sym, input, oneD, &
                  vCoul, den, results%te_vcoul )

             WRITE (6,FMT=8030) results%te_vcoul
             WRITE (16,FMT=8030) results%te_vcoul
8030         FORMAT (/,10x,'total density-coulomb potential integral :', t40,f20.10)

             CALL timestop( "den-pot integrals" )
          END IF
       ENDIF
    ENDIF !irank==0

  END SUBROUTINE vgen_coulomb

END MODULE m_vgen_coulomb
