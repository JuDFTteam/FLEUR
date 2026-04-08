MODULE m_types_BEC_old
    ! type introduced to simplify Born effective charge output
    USE m_types
    IMPLICIT NONE
    PUBLIC
    TYPE t_BEC_old
        complex, allocatable ::  BEC_element(:,:,:)
        complex, allocatable ::  BEC_contributions(:,:,:,:)

    CONTAINS
        !PROCEDURE :: init_BEC_element
        !PROCEDURE :: init_BEC_contributions
        !GENERIC   :: init=>init_BEC_element,init_BEC_contributions
        PROCEDURE :: init
    END TYPE t_BEC_old

CONTAINS
    SUBROUTINE init(this,fi)
        implicit none
        class(t_BEC_old),INTENT(INOUT)        :: this
        TYPE(t_fleurinput),INTENT(IN)     :: fi

        ALLOCATE(this%BEC_element(fi%atoms%nat,3,3))
        this%BEC_element = cmplx(0.0,0.0)
        ALLOCATE(this%BEC_contributions(fi%atoms%nat,3,3,8+fi%atoms%nat))
        this%BEC_contributions = cmplx(0.0,0.0)

    END SUBROUTINE init
        
END MODULE m_types_BEC_old