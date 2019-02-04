!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_mixvector
  use m_types
  implicit none
  PRIVATE
  !Here we store the pointers used for metric
  TYPE(t_oneD),POINTER   :: oneD
  TYPE(t_input),POINTER  :: input
  TYPE(t_vacuum),POINTER :: vacuum
  TYPE(t_noco),POINTER   :: noco
  TYPE(t_sym),POINTER    :: sym
  TYPE(t_stars),POINTER  :: stars
  TYPE(t_cell),POINTER   :: cell
  TYPE(t_sphhar),POINTER :: sphhar
  TYPE(t_atoms),POINTER  :: atoms  =>null()  
  INTEGER                :: mmap, mmaph, nmaph, nmap, mapmt, mapvac, mapvac2
  real                   :: intfac,vacfac
  
  TYPE,PUBLIC:: t_mixvector
     REAL,ALLOCATABLE       :: vec(:)
     LOGICAL                :: l_pot=.false. !Is this a potential?
   CONTAINS
     PROCEDURE :: init=>mixvector_init
     procedure :: alloc=>mixvector_alloc
     PROCEDURE :: from_density=>mixvector_from_density
     PROCEDURE :: to_density=>mixvector_to_density
     PROCEDURE :: apply_metric=>mixvector_metric
  END TYPE t_mixvector

  INTERFACE assignment(=)
     MODULE PROCEDURE assign_vectors
  END INTERFACE assignment(=)

  INTERFACE OPERATOR (*)
     MODULE PROCEDURE multiply_scalar
  END INTERFACE OPERATOR (*)
  INTERFACE OPERATOR (+)
     MODULE PROCEDURE add_vectors
  END INTERFACE OPERATOR (+)
  INTERFACE OPERATOR (-)
     MODULE PROCEDURE subtract_vectors
  END INTERFACE OPERATOR (-)
  INTERFACE OPERATOR (.dot.)
     MODULE PROCEDURE multiply_dot
  END INTERFACE OPERATOR (.dot.)

  public :: Operator(+),operator(-),operator(*),operator(.dot.)
  public :: assignment(=)
  
  
CONTAINS

  SUBROUTINE mixvector_from_density(vec,den)
    USE m_types
    USE m_brysh1
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(INOUT)    :: vec
    TYPE(t_potden),    INTENT(in)       :: Den

    CALL brysh1( input, stars, atoms, sphhar, noco, vacuum, sym, oneD, &
         intfac, vacfac, Den, nmap, nmaph, mapmt, mapvac, mapvac2, vec%vec)
  END SUBROUTINE mixvector_from_density

  SUBROUTINE mixvector_to_density(vec,den)
    USE m_types
    USE m_brysh2
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(IN)    :: vec
    TYPE(t_potden),    INTENT(OUT)       :: Den

    CALL brysh2( input, stars, atoms, sphhar, noco, vacuum, sym, vec%vec,oneD,den)
  END SUBROUTINE mixvector_to_density


  FUNCTION mixvector_metric(vec)RESULT(mvec)
    USE m_types
    USE m_metric
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(IN)    :: vec
    TYPE(t_mixvector)                :: mvec

    mvec=vec
    CALL metric( cell, atoms, vacuum, sphhar, input, noco, stars, sym, oneD, &
         mmap, nmaph, mapmt, mapvac2, vec%vec, mvec%vec, vec%l_pot )
  END FUNCTION mixvector_metric

  SUBROUTINE mixvector_init(vec,oneD_i,input_i,vacuum_i,noco_i,sym_i,stars_i,cell_i,sphhar_i,atoms_i)
    USE m_types
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(OUT)    :: vec
    TYPE(t_oneD),INTENT(IN),TARGET   :: oneD_i
    TYPE(t_input),INTENT(IN),TARGET  :: input_i
    TYPE(t_vacuum),INTENT(IN),TARGET :: vacuum_i
    TYPE(t_noco),INTENT(IN),TARGET   :: noco_i
    TYPE(t_sym),INTENT(IN),TARGET    :: sym_i
    TYPE(t_stars),INTENT(IN),TARGET  :: stars_i
    TYPE(t_cell),INTENT(IN),TARGET   :: cell_i
    TYPE(t_sphhar),INTENT(IN),TARGET :: sphhar_i
    TYPE(t_atoms),INTENT(IN),TARGET  :: atoms_i

    if(.not.associated(atoms)) then
    !Store pointers to data-types
    oneD=>oneD_i;input=>input_i;vacuum=>vacuum_i;noco=>noco_i
    sym=>sym_i;stars=>stars_i;cell=>cell_i;sphhar=>sphhar_i;atoms=>atoms_i

    !In systems without inversions symmetry the interstitial star-
    !coefficients are complex. Thus twice as many numbers have to be
    !stored.
    intfac=MERGE(1,2,sym%invs)

    !The corresponding is true for the coeff. of the warping vacuum
    !density depending on the two dimensional inversion.
    vacfac=MERGE(1,2,sym%invs2)

    mmaph = intfac * stars%ng3 + atoms%ntype * ( sphhar%nlhd + 1 ) * atoms%jmtd + &
         vacfac * vacuum%nmzxyd * ( oneD%odi%n2d - 1 ) * vacuum%nvac + vacuum%nmzd * vacuum%nvac
    mmap  =mmaph * input%jspins
    !in a non-collinear calculations extra space is needed for the
    !off-diag. part of the density matrix. these coeff. are generally
    !complex independ of invs and invs2.
    IF ( noco%l_noco ) THEN
       mmap = mmap + 2 * stars%ng3 + 2 * vacuum%nmzxyd * ( oneD%odi%n2d - 1 ) * vacuum%nvac + &
            2 * vacuum%nmzd * vacuum%nvac
       IF (noco%l_mtnocopot) mmap= mmap+ 2*atoms%ntype * ( sphhar%nlhd + 1 ) * atoms%jmtd 
    END IF

    ! LDA+U (start)
    PRINT *,"MIXING of LDA+U missing....."
    !n_mmpTemp = inDen%mmpMat
    !n_u_keep = atoms%n_u
    !IF ( atoms%n_u > 0 ) CALL u_mix( input, atoms, inDen%mmpMat, outDen%mmpMat )
    !IF ( ANY( inDen%mmpMat(:,:,:,:) /= 0.0 ) ) THEN
    !    !In an LDA+U caclulation, also the density matrix is included in the
    !    !supervectors (sm,fsm) if no linear mixing is performed on it.
    !    IF ( input%ldauLinMix ) THEN
    !       atoms%n_u = 0
    !    ELSE
    !       mmap = mmap + 7 * 7 * 2 * atoms%n_u * input%jspins ! add 7*7 complex numbers per atoms%n_u and spin
    !    END IF
    ! ELSE
    !    atoms%n_u = 0
    ! END IF
 endif
 call vec%alloc()
 SUBROUTINE mixvector_alloc(vec)
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(OUT)    :: vec
    ALLOCATE( vec%vec(mmap) )
  END SUBROUTINE mixvector_alloc

  !The operators
  SUBROUTINE assign_vectors(vec,vecin)
    TYPE(t_mixvector),INTENT(OUT)::vec
    TYPE(t_mixvector),INTENT(IN) ::vecin
    vec=vecin
  END SUBROUTINE assign_vectors

  FUNCTION multiply_scalar(scalar,vec)RESULT(vecout)
    TYPE(t_mixvector),INTENT(IN)::vec
    REAL,INTENT(IN)             ::scalar
    TYPE(t_mixvector)           ::vecout

    vecout=vec
    vecout%vec=vecout%vec*scalar
  END FUNCTION multiply_scalar

  FUNCTION add_vectors(vec1,vec2)RESULT(vecout)
    TYPE(t_mixvector),INTENT(IN)::vec1,vec2
    TYPE(t_mixvector)           ::vecout

    vecout=vec1
    vecout%vec=vec1%vec+vec2%vec
  END FUNCTION add_vectors

  FUNCTION multiply_dot(vec1,vec2)RESULT(dprod)
    TYPE(t_mixvector),INTENT(IN)::vec1,vec2
    REAL                        ::dprod

    dprod=dot_PRODUCT(vec1%vec,vec2%vec)
  END FUNCTION multiply_dot

  FUNCTION subtract_vectors(vec1,vec2)RESULT(vecout)
    TYPE(t_mixvector),INTENT(IN)::vec1,vec2
    TYPE(t_mixvector)           ::vecout

    vecout=vec1
    vecout%vec=vec1%vec-vec2%vec
  END FUNCTION subtract_vectors
end MODULE m_types_mixvector
