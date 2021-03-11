!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_doubleCounting

   USE m_constants

   IMPLICIT NONE

   CONTAINS

   FUNCTION doubleCountingPot(U,J,l,l_amf,l_spinAvg,rho,l_write) RESULT(Vdc)

      !------------------------------------------------------
      ! Calculate the Double Counting Correction in either
      ! the FLL or AMF limit
      ! If l_spinAvg is True the Double counting will be
      ! averaged over the spins
      !------------------------------------------------------

      REAL,                INTENT(IN)  :: U          !Hubbard parameters
      REAL,                INTENT(IN)  :: J
      INTEGER,             INTENT(IN)  :: l
      LOGICAL,             INTENT(IN)  :: l_amf      !Which doubleCounting is used (FLL/AMF)
      LOGICAL,             INTENT(IN)  :: l_spinAvg  !Do we want a spin averaged double counting
      REAL,                INTENT(IN)  :: rho(:)     !Trace of the density matrix for each spin
      LOGICAL, OPTIONAL,   INTENT(IN)  :: l_write

      REAL :: Vdc(SIZE(rho))
      REAL :: nup, ndn, Vdcup, Vdcdn


      nup = rho(1)
      IF(SIZE(rho) == 2) THEN
         ndn = rho(2)
      ELSE
         ndn = rho(1)
      ENDIF

      IF(l_amf) THEN
         Vdcup = U*ndn+2.0*l/(2.0*l+1)*(U-J)*nup
         Vdcdn = U*nup+2.0*l/(2.0*l+1)*(U-J)*ndn
      ELSE
         Vdcup = U*(nup+ndn - 0.5) - J*(nup - 0.5)
         Vdcdn = U*(nup+ndn - 0.5) - J*(ndn - 0.5)
      ENDIF

      IF(PRESENT(l_write)) THEN
         IF(l_write) THEN
            WRITE(oUnit,"(/,A)") 'Double counting chemical potential:'
            IF(l_amf) THEN
               WRITE(oUnit,9040) 'AMF: ','spin-up','spin-dn','(up+dn)/2','up-dn'
            ELSE
               WRITE(oUnit,9040) 'FLL: ','spin-up','spin-dn','(up+dn)/2','up-dn'
            ENDIF
9040        FORMAT(TR3,A4,TR1,A7,TR3,A7,TR3,A9,TR3,A5)
            WRITE(oUnit,9050) Vdcup,Vdcdn,(Vdcup+Vdcdn)/2.0,Vdcup-Vdcdn
9050        FORMAT(TR7,f8.4,TR2,f8.4,TR2,f8.4,TR4,f8.4)
         ENDIF
      ENDIF

      IF(l_spinAvg) THEN
         Vdc(:) = (Vdcup+Vdcdn)/2.0
      ELSE
         Vdc(1) = Vdcup
         IF(SIZE(rho) == 2) THEN
            Vdc(2) = Vdcdn
         ENDIF
      ENDIF

      IF(PRESENT(l_write)) THEN
         IF(l_write) THEN
            WRITE(oUnit,*) "Vdc = ", Vdc(:)
         ENDIF
      ENDIF

   END FUNCTION doubleCountingPot


   FUNCTION doubleCountingEnergy(U,J,l,l_amf,l_spinAvg,rho) RESULT(Edc)

      !------------------------------------------------------------
      ! Calculate the Double Counting Correction Energy in either
      ! the FLL or AMF limit
      !------------------------------------------------------------

      REAL,                INTENT(IN)  :: U          !Hubbard parameters
      REAL,                INTENT(IN)  :: J
      INTEGER,             INTENT(IN)  :: l
      LOGICAL,             INTENT(IN)  :: l_amf      !Which doubleCounting is used (FLL/AMF)
      LOGICAL,             INTENT(IN)  :: l_spinAvg
      REAL,                INTENT(IN)  :: rho(:)     !Trace of the density matrix for each spin

      REAL :: Edc
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

      IF(l_amf) THEN
         Edc = U*nup*ndn + (U-J) *l*(nup**2+ndn**2)/(2.0*l+1)
      ELSE
         Edc =   U/2.0 * (nup+ndn) * (nup+ndn - 1.0) &
               - J/2.0 * nup       * (nup     - 1.0) &
               - J/2.0 * ndn       * (ndn     - 1.0)
      ENDIF

   END FUNCTION doubleCountingEnergy

END MODULE m_doubleCounting