      MODULE icorrkeys

      INTEGER, PARAMETER    ::  icorr_exx    = -3
      INTEGER, PARAMETER    ::  icorr_hf     = -2
      INTEGER, PARAMETER    ::  icorr_pbe0   = 12
c     HSE, hybrid exchange functional: J. Chem. Phys. 118, 8207 (2003)
      INTEGER, PARAMETER    ::  icorr_hse    = 13
      ! only local part of HSE
      INTEGER, PARAMETER    ::  icorr_hseloc = 14
      ! hybrid functional similar to HSE but with variable screening and mixing parameter
      INTEGER, PARAMETER    ::  icorr_vhse   = 15

      END MODULE icorrkeys

