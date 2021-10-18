!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_doubleCounting

   USE m_constants
   USE m_types

   IMPLICIT NONE

   CONTAINS

   FUNCTION doubleCountingPot(density, ldau, U,J, umatrix, l_spinoffd, l_mix,l_spinAvg, alpha, l_write) RESULT(Vdc)

      !------------------------------------------------------
      ! Calculate the Double Counting Correction in either
      ! the FLL or AMF limit
      ! If l_spinAvg is True the Double counting will be
      ! averaged over the spins
      !------------------------------------------------------

      COMPLEX,             INTENT(IN)  :: density(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_utype),       INTENT(IN)  :: ldau       !LDA+U information
      REAL,                INTENT(IN)  :: umatrix(-lmaxU_const:,-lmaxU_const:,-lmaxU_const:,-lmaxU_const:)
      REAL,                INTENT(IN)  :: U, J
      LOGICAL,             INTENT(IN)  :: l_spinoffd
      LOGICAL,             INTENT(IN)  :: l_mix      !Mix between FLL and AMF
      LOGICAL,             INTENT(IN)  :: l_spinAvg  !Do we want a spin averaged double counting
      REAL,                INTENT(IN)  :: alpha
      LOGICAL, OPTIONAL,   INTENT(IN)  :: l_write

      COMPLEX,ALLOCATABLE :: Vdc(:,:,:)

      REAL :: charge, mag(3)
      COMPLEX :: sigma(2,2,3)
      INTEGER :: spin_dim, ispin,m, spin1,spin2

      sigma(1,2,1)=CMPLX(1.0,0.0)
      sigma(2,1,1)=CMPLX(1.0,0.0)
      sigma(1,2,2)=CMPLX(0.0,-1.0)
      sigma(2,1,2)=CMPLX(0.0,1.0)
      sigma(1,1,3)=CMPLX(1.0,0.0)
      sigma(2,2,3)=CMPLX(-1.0,0.0)

      spin_dim = SIZE(density,3)
      IF(.NOT.l_spinoffd) spin_dim = MIN(2,spin_dim)

      ALLOCATE(Vdc ,mold=density)
      Vdc = cmplx_0

      IF(ldau%l_amf) THEN
         CALL juDFT_error("DADADADA")
      ELSE

         charge = 0.0
         mag = 0.0

         DO ispin = 1, spin_dim
            DO m = -ldau%l, ldau%l
               IF(ispin < 3) THEN
                  charge = charge + REAL(density(m,m,ispin))
                  mag(3) = mag(3) + (-1)**(ispin-1) * REAL(density(m,m,ispin))
               ELSE
                  mag(1) = mag(1) + 2 * REAL(density(m,m,ispin))
                  mag(2) = mag(2) + 2 * AIMAG(density(m,m,ispin))
               ENDIF
            ENDDO
         ENDDO

         IF(spin_dim == 1) THEN
            charge = charge/2.0
            mag = 0.0
         ENDIF

         WRITE(*,*) charge, mag

         DO ispin = 1, spin_dim

            IF(ispin==3) THEN
               spin1 = 2
               spin2 = 1
            ELSE
               spin1 = ispin
               spin2 = ispin
            ENDIF

            DO m = -ldau%l, ldau%l
               IF(spin1 == spin2) Vdc(m,m,ispin) = Vdc(m,m,ispin) + U*(2*charge-1)/2.0 - J*(charge-1)/2.0
               Vdc(m,m,ispin) = Vdc(m,m,ispin) - J/2.0*dot_product(mag,sigma(spin1,spin2,:))
            ENDDO
         ENDDO
      ENDIF

!       IF(PRESENT(l_write)) THEN
!          IF(l_write) THEN
!             WRITE(oUnit,"(/,A)") 'Double counting chemical potential:'
!             IF(l_amf) THEN
!                WRITE(oUnit,9040) 'AMF: ','spin-up','spin-dn','(up+dn)/2','up-dn'
!             ELSE
!                WRITE(oUnit,9040) 'FLL: ','spin-up','spin-dn','(up+dn)/2','up-dn'
!             ENDIF
! 9040        FORMAT(TR3,A4,TR1,A7,TR3,A7,TR3,A9,TR3,A5)
!             WRITE(oUnit,9050) Vdcup,Vdcdn,(Vdcup+Vdcdn)/2.0,Vdcup-Vdcdn
! 9050        FORMAT(TR7,f8.4,TR2,f8.4,TR2,f8.4,TR4,f8.4)
!          ENDIF
!       ENDIF

      ! IF(l_spinAvg) THEN
      !    Vdc(:) = (Vdcup+Vdcdn)/2.0
      ! ELSE
      !    Vdc(1) = Vdcup
      !    IF(SIZE(rho) == 2) THEN
      !       Vdc(2) = Vdcdn
      !    ENDIF
      ! ENDIF

      ! IF(PRESENT(l_write)) THEN
      !    IF(l_write) THEN
      !       WRITE(oUnit,*) "Vdc = ", Vdc(:)
      !    ENDIF
      ! ENDIF

   END FUNCTION doubleCountingPot


   FUNCTION doubleCountingEnergy(U,J,l,l_amf,l_mix,l_spinAvg,rho, alpha) RESULT(Edc)

      !------------------------------------------------------------
      ! Calculate the Double Counting Correction Energy in either
      ! the FLL or AMF limit
      !------------------------------------------------------------

      REAL,                INTENT(IN)  :: U          !Hubbard parameters
      REAL,                INTENT(IN)  :: J
      INTEGER,             INTENT(IN)  :: l
      LOGICAL,             INTENT(IN)  :: l_amf      !Which doubleCounting is used (FLL/AMF)
      LOGICAL,             INTENT(IN)  :: l_mix      !Mix between FLL/AMF
      LOGICAL,             INTENT(IN)  :: l_spinAvg
      REAL,                INTENT(IN)  :: rho(:)     !Trace of the density matrix for each spin
      REAL,                INTENT(IN)  :: alpha

      REAL :: Edc, Edc_AMF, Edc_FLL
      REAL :: nup, ndn


      nup = rho(1)
      IF(SIZE(rho) == 2) THEN
         ndn = rho(2)
      ELSE
         ndn = rho(1)
      ENDIF

      IF(l_spinAvg) THEN
         nup = (nup+ndn)/2.0
         ndn = nup
      ENDIF

      Edc_AMF = U*nup*ndn + (U-J) *l*(nup**2+ndn**2)/(2.0*l+1)

      Edc_FLL =   U/2.0 * (nup+ndn) * (nup+ndn - 1.0) &
                - J/2.0 * nup       * (nup     - 1.0) &
                - J/2.0 * ndn       * (ndn     - 1.0)


      IF(l_amf) THEN
         Edc = Edc_AMF
      ELSE
         Edc = Edc_FLL
      ENDIF

   END FUNCTION doubleCountingEnergy

   FUNCTION doubleCountingMixFactor(mmpmat, l, charge) Result(alpha)

      !---------------------------------------------
      ! Calculate the mixing factor between FLL/AMF
      !---------------------------------------------

      COMPLEX,        INTENT(IN) :: mmpmat(-lmaxU_const:, -lmaxU_const:, :)
      INTEGER,        INTENT(IN) :: l
      REAL,           INTENT(IN) :: charge

      REAL :: alpha

      alpha = 0.0

   END FUNCTION

END MODULE m_doubleCounting