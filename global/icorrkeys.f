      MODULE m_icorrkeys
      INTEGER, PARAMETER    ::  icorr_exx    = -3
      INTEGER, PARAMETER    ::  icorr_pbe0   = 12
      INTEGER, PARAMETER    ::  icorr_hf     = -2
!HSE, hybrid exchange functional: J. Chem. Phys. 118, 8207 (2003)
      INTEGER, PARAMETER    ::  icorr_hse    = 13
! only local part of HSE
      INTEGER, PARAMETER    ::  icorr_hseloc = 14
! hybrid functional similar to HSE but with variable screening and mixing parameter
      INTEGER, PARAMETER    ::  icorr_vhse   = 15
      
      CONTAINS
      
      FUNCTION get_exchange_weight(icorr) result(a_ex)
      USE m_hybridmix
      USE m_judft
      IMPLICIT NONE
      INTEGER,INTENT(IN)::icorr
      REAL a_ex
      
      SELECT CASE (icorr)
      CASE (icorr_pbe0)
         a_ex = amix_pbe0
      CASE (icorr_hf   ) 
         a_ex = amix_hf
      CASE ( icorr_hse) 
         a_ex = aMix_HSE
      CASE (icorr_vhse )
         a_ex = aMix_VHSE()
      CASE DEFAULT
         call judft_error('xc functional can not be identified')
      END SELECT
      END
      
      END MODULE m_icorrkeys

