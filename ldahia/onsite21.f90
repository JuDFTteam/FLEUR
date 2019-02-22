MODULE m_onsite21

   USE m_juDFT
   USE m_constants

   CONTAINS

   SUBROUTINE onsite21(atoms,sym,jspins,noccbd,tetweights,wtkpt,eig,usdus,denCoeffsOffDiag,eigVecCoeffs,gOnsite)

      USE m_types
      USE m_differentiate

      IMPLICIT NONE

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_eigVecCoeffs),      INTENT(IN)     :: eigVecCoeffs
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffDiag
      TYPE(t_usdus),             INTENT(IN)     :: usdus
      TYPE(t_greensf),           INTENT(INOUT)  :: gOnsite

      INTEGER,                   INTENT(IN)     :: jspins
      INTEGER,                   INTENT(IN)     :: noccbd 

      REAL,                      INTENT(IN)     :: wtkpt
      REAL,                      INTENT(IN)     :: tetweights(:,:)
      REAL,                      INTENT(IN)     :: eig(noccbd)



   END SUBROUTINE

END MODULE m_onsite21