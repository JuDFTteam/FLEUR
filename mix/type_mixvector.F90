!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
  !--------------------------------------------------------------------------------
MODULE m_type_mixvector
  TYPE t_mixvector
     REAL,ALLOCATABLE(:) :: vec(:)
     INTEGER             :: mmap, mmaph, nmaph, nmap, mapmt, mapvac, mapvac2,intfac
     LOGICAL             :: l_pot !Is this a potential?
     !Here we store the pointers used for metric
     TYPE(t_oneD),POINTER   :: oneD
     TYPE(t_input),POINTER  :: input
     TYPE(t_vacuum),POINTER :: vacuum
     TYPE(t_noco),POINTER   :: noco
     TYPE(t_sym),POINTER    :: sym
     TYPE(t_stars),POINTER  :: stars
     TYPE(t_cell),POINTER   :: cell
     TYPE(t_sphhar),POINTER :: sphhar
     TYPE(t_atoms),POINTER  :: atoms    
   CONTAINS
     PROCEDURE init=>mixvector_init
     PROCEDURE from_density=>mixvector_from_density
     PROCEDURE to_density=>mixvector_to_density
     PROCEDURE apply_metric=>mixvector_metric
  END TYPE t_mixvector

  INTERFACE assignement(=)
     MODULE PROCEDURE assign_vectors
  END INTERFACE assignement

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
     
  PRIVATE
CONTAINS

  SUBROUTINE mixvector_from_density(vec,den)
    USE m_types
    USE m_brysh1
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(INOUT)    :: vec
    TYPE(t_potden),    INTENT(in)       :: Den

    CALL brysh1( vec%input, vec%stars, vec%atoms, vec%sphhar, vec%noco, vec%vacuum, vec%sym, vec%oneD, &
         vec%intfac, vec%vacfac, Den, vec%nmap, vec%nmaph, vec%mapmt, vec%mapvac, vec%mapvac2, vec%vec)
  END SUBROUTINE mixvector_from_density

  SUBROUTINE mixvector_to_density(vec,den)
    USE m_types
    USE m_brysh2
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(IN)    :: vec
    TYPE(t_potden),    INTENT(OUT)       :: Den

    CALL brysh2( vec%input, vec%stars, vec%atoms, vec%sphhar, vec%noco, vec%vacuum, vec%sym, vec%vec,vec%oneD,den)
  END SUBROUTINE mixvector_to_density


  FUNCTION mixvector_metric(vec)RESULT(mvec)
    USE m_types
    USE m_metric
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(IN)    :: vec
    TYPE(t_mixvector)                :: mvec

    mvec=vec
    CALL metric( vec%cell, vec%atoms, vec%vacuum, vec%sphhar, vec%input, vec%noco, vec%stars, vec%sym, vec%oneD, &
         vec%mmap, vec%nmaph, vec%mapmt, vec%mapvac2, vec%vec, mvec%vec, vec%l_pot )
  END FUNCTION mixvector_metric
    
  SUBROUTINE mixvector_init(vec,oneD,input,vacuum,noco,sym,stars,cell,sphhar,atoms)
    USE m_types
    IMPLICIT NONE
    CLASS(t_mixvector),INTENT(OUT)    :: vec
    TYPE(t_oneD),INTENT(IN),TARGET   :: oneD
    TYPE(t_input),INTENT(IN),TARGET  :: input
    TYPE(t_vacuum),INTENT(IN),TARGET :: vacuum
    TYPE(t_noco),INTENT(IN),TARGET   :: noco
    TYPE(t_sym),INTENT(IN),TARGET    :: sym
    TYPE(t_stars),INTENT(IN),TARGET  :: stars
    TYPE(t_cell),INTENT(IN),TARGET   :: cell
    TYPE(t_sphhar),INTENT(IN),TARGET :: sphhar
    TYPE(t_atoms),INTENT(IN),TARGET  :: atoms


    !Store pointers to data-types
    vec%oneD=>oneD;vec%input=>input;vec%vacuum=>vacuum;vec%noco=>noco
    vec%sym=>sym;vec%stars=>stars;vec%cell=>cell;vec%sphhar=>sphhar;vec%atoms=>atoms
     
    !In systems without inversions symmetry the interstitial star-
    !coefficients are complex. Thus twice as many numbers have to be
    !stored.
    vec%intfac=MERGE(1.0,2.0,sym%invs)
    
    !The corresponding is true for the coeff. of the warping vacuum
    !density depending on the two dimensional inversion.
    vec%vacfac=MERGE(1.0,2.0,sym%invs2)
     
    vec%mmaph = vec%intfac * stars%ng3 + atoms%ntype * ( sphhar%nlhd + 1 ) * atoms%jmtd + &
         vacfac * vacuum%nmzxyd * ( oneD%odi%n2d - 1 ) * vacuum%nvac + vacuum%nmzd * vacuum%nvac
    vec%mmap  =vec%mmaph * input%jspins
    !in a non-collinear calculations extra space is needed for the
    !off-diag. part of the density matrix. these coeff. are generally
    !complex independ of invs and invs2.
    IF ( noco%l_noco ) THEN
       vec%mmap = vec%mmap + 2 * stars%ng3 + 2 * vacuum%nmzxyd * ( oneD%odi%n2d - 1 ) * vacuum%nvac + &
            2 * vacuum%nmzd * vacuum%nvac
       IF (noco%l_mtnocopot) vec%mmap= vec%mmap+ 2*atoms%ntype * ( sphhar%nlhd + 1 ) * atoms%jmtd 
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
     
    ALLOCATE( vec%vec(mmap) )
     
  END SUBROUTINE mixvector_init

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

      vecout=vecin
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
    
    FUNCTION add_vectors(vec1,vec2)RESULT(vecout)
      TYPE(t_mixvector),INTENT(IN)::vec1,vec2
      TYPE(t_mixvector)           ::vecout

      vecout=vec1
      vecout%vec=vec1%vec-vec2%vec
    END FUNCTION add_vectors
