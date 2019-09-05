!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_calcDenCoeffs

CONTAINS

SUBROUTINE calcDenCoeffs(atoms,sphhar,sym,we,noccbd,eigVecCoeffs,ispin,denCoeffs)

   USE m_juDFT
   USE m_types_setup
   USE m_types_cdnval, ONLY: t_eigVecCoeffs,t_denCoeffs
   USE m_rhomt
   USE m_rhonmt
   USE m_rhomtlo
   USE m_rhonmtlo

   IMPLICIT NONE

   TYPE(t_atoms),        INTENT(IN)    :: atoms
   TYPE(t_sphhar),       INTENT(IN)    :: sphhar
   TYPE(t_sym),          INTENT(IN)    :: sym
   TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs
   TYPE(t_denCoeffs),    INTENT(INOUT) :: denCoeffs

   INTEGER,              INTENT(IN)    :: noccbd
   REAL,                 INTENT(IN)    :: we(noccbd)

   INTEGER,              INTENT(IN)    :: ispin

   !--->          set up coefficients for the spherical and
   CALL timestart("cdnval: rhomt")
   CALL rhomt(atoms,we,noccbd,eigVecCoeffs,denCoeffs,ispin)
   CALL timestop("cdnval: rhomt")

   !--->          non-spherical m.t. density
   CALL timestart("cdnval: rhonmt")
   CALL rhonmt(atoms,sphhar,we,noccbd,sym,eigVecCoeffs,denCoeffs,ispin)
   CALL timestop("cdnval: rhonmt")

   !--->          set up coefficients of the local orbitals and the
   !--->          flapw - lo cross terms for the spherical and
   !--->          non-spherical mt density
   CALL timestart("cdnval: rho(n)mtlo")
   CALL rhomtlo(atoms,noccbd,we,eigVecCoeffs,denCoeffs,ispin)

   CALL rhonmtlo(atoms,sphhar,noccbd,we,eigVecCoeffs,denCoeffs,ispin)
   CALL timestop("cdnval: rho(n)mtlo")

END SUBROUTINE calcDenCoeffs

END MODULE m_calcDenCoeffs
