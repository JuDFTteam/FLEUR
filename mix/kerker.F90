--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_kerker
  contains

    SUBROUTINE kerker(field, DIMENSION, mpi, &
                stars, atoms, sphhar, vacuum, input, sym, cell, noco, &
                oneD, inDen, outDen, fsm ,mapmt,mapvac,mapvac2,nmap,nmaph  )
      !Implementation of the Kerker preconditioner by M.Hinzen
      USE m_vgen_coulomb
      USE m_VYukawaFilm
      USE m_juDFT
      USE m_types
      IMPLICIT NONE
      TYPE(t_oneD),      INTENT(in)    :: oneD
      TYPE(t_input),     INTENT(in)    :: input
      TYPE(t_vacuum),    INTENT(in)    :: vacuum
      TYPE(t_noco),      INTENT(in)    :: noco
      TYPE(t_sym),       INTENT(in)    :: sym
      TYPE(t_stars),     INTENT(in)    :: stars
      TYPE(t_cell),      INTENT(in)    :: cell
      TYPE(t_sphhar),    INTENT(in)    :: sphhar
      TYPE(t_field),     INTENT(in)    :: field
      TYPE(t_dimension), INTENT(in)    :: DIMENSION
      TYPE(t_mpi),       INTENT(in)    :: mpi
      TYPE(t_atoms),     INTENT(in)    :: atoms 
      TYPE(t_potden),    INTENT(in)    :: outDen
      TYPE(t_potden),    INTENT(inout) :: inDen
      
      TYPE(t_potden),    INTENT(in) :: outDen,inDen
      REAL, ALLOCATABLE,INTENT(OUT) :: fsm(:)
      INTEGER, INTENT (OUT) :: mapmt,mapvac,mapvac2,nmap,nmaph

      !Locals
      type(t_potden)                :: resDen, vYukawa
    
      
      CALL resDen%init( stars, atoms, sphhar, vacuum, noco, input%jspins, POTDEN_TYPE_DEN )
      CALL vYukawa%init( stars, atoms, sphhar, vacuum, noco, input%jspins, 4 )
      MPI0_b: IF( mpi%irank == 0 ) THEN 
         CALL resDen%subPotDen( outDen, inDen )
         IF( input%jspins == 2 ) CALL resDen%SpinsToChargeAndMagnetisation()
      END IF MPI0_b
#ifdef CPP_MPI
      CALL mpi_bc_potden( mpi, stars, sphhar, atoms, input, vacuum, oneD, noco, resDen )
#endif
      IF ( .NOT. input%film ) THEN
         CALL vgen_coulomb( 1, mpi, DIMENSION, oneD, input, field, vacuum, sym, stars, cell, &
              sphhar, atoms, resDen, vYukawa )
      ELSE
         vYukawa%iter = resDen%iter
         CALL VYukawaFilm( stars, vacuum, cell, sym, input, mpi, atoms, sphhar, DIMENSION, oneD, resDen, &
              vYukawa )
      END IF
      
      MPI0_c: IF( mpi%irank == 0 ) THEN
         IF( input%preconditioning_param /= 0 ) THEN
            resDen%pw(1:stars%ng3,1) = resDen%pw(1:stars%ng3,1) - input%preconditioning_param ** 2 / fpi_const * vYukawa%pw(1:stars%ng3,1)
            DO n = 1, atoms%ntype
               DO lh = 0, sphhar%nlhd
                  resDen%mt(1:atoms%jri(n),lh,n,1) = resDen%mt(1:atoms%jri(n),lh,n,1) &
                       - input%preconditioning_param ** 2 / fpi_const &
                       * vYukawa%mt(1:atoms%jri(n),lh,n,1) * atoms%rmsh(1:atoms%jri(n),n) ** 2
               END DO
            END DO
            resDen%vacz  = resDen%vacz  - input%preconditioning_param ** 2 / fpi_const * vYukawa%vacz
            resDen%vacxy = resDen%vacxy - input%preconditioning_param ** 2 / fpi_const * vYukawa%vacxy
            IF( input%jspins == 2 ) CALL resDen%ChargeAndMagnetisationToSpins()
            ! fix the preconditioned density
            CALL outDen%addPotDen( resDen, inDen )
            CALL qfix(mpi,stars, atoms, sym, vacuum, sphhar, input, cell, oneD, outDen, noco%l_noco, .FALSE., .TRUE., fix )
            CALL resDen%subPotDen( outDen, inDen )
            CALL brysh1( input, stars, atoms, sphhar, noco, vacuum, sym, oneD, &
                 intfac, vacfac, resDen, nmap, nmaph, mapmt, mapvac, mapvac2, fsm )
         END IF
         ! end of preconditioner
      END IF MPI0_c
    END SUBROUTINE kerker
