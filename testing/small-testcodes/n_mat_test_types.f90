MODULE m_constants
   IMPLICIT NONE
   INTEGER, PARAMETER :: lmaxU_const = 3
   COMPLEX, PARAMETER :: cmplx_0 = (0.0, 0.0)
END MODULE m_constants

MODULE m_types_radfun
   IMPLICIT NONE
   TYPE t_radfun
      REAL, ALLOCATABLE :: integral(:, :, :, :, :)
   END TYPE t_radfun
END MODULE m_types_radfun

MODULE m_types_abc
   IMPLICIT NONE
   TYPE t_abc
      COMPLEX, ALLOCATABLE :: cof(:, :, :, :)
   END TYPE t_abc
END MODULE m_types_abc

MODULE m_types
   IMPLICIT NONE
   TYPE t_utype
      INTEGER :: atomType = 0
      INTEGER :: l = -1
   END TYPE t_utype
   TYPE t_opctype
      INTEGER :: atomType = 0
      INTEGER :: l = -1
   END TYPE t_opctype
   TYPE t_atoms
      INTEGER :: n_u = 0
      INTEGER :: n_opc = 0
      INTEGER :: n_hia = 0
      INTEGER, ALLOCATABLE :: neq(:), firstAtom(:)
      TYPE(t_utype), ALLOCATABLE :: lda_u(:)
      TYPE(t_opctype), ALLOCATABLE :: lda_opc(:)
   END TYPE t_atoms
   TYPE t_sym
   END TYPE t_sym
END MODULE m_types

MODULE m_symMMPmat
   USE m_constants, ONLY: lmaxU_const
   USE m_types, ONLY: t_sym
   IMPLICIT NONE
CONTAINS
   FUNCTION symMMPmat(mmpmat, sym, natom, l) RESULT(mmpmatSym)
      COMPLEX, INTENT(IN) :: mmpmat(-lmaxU_const:, -lmaxU_const:)
      TYPE(t_sym), INTENT(IN) :: sym
      INTEGER, INTENT(IN) :: natom, l
      COMPLEX :: mmpmatSym(-lmaxU_const:lmaxU_const, -lmaxU_const:lmaxU_const)
      mmpmatSym = mmpmat
   END FUNCTION symMMPmat
END MODULE m_symMMPmat
