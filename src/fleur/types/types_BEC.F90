MODULE m_types_BEC

    USE m_types
    IMPLICIT NONE
    PUBLIC
    TYPE t_BEC
        complex, allocatable ::  BEC_element(:,:,:)
        complex, allocatable ::  BEC_contributions(:,:,:)

    CONTAINS
        PROCEDURE :: init_BEC_element
        GENERIC   :: init=>init_BEC_element
    END TYPE t_BEC

CONTAINS
    SUBROUTINE init_BEC_element(this,fi)
        implicit none
        class(t_BEC),INTENT(INOUT)        :: this
        TYPE(t_fleurinput),INTENT(IN)     :: fi

        ALLOCATE(this%BEC_element(fi%atoms%ntype,3,3))
        this%BEC_element = cmplx(0.0,0.0)
        print*,"IM INIT"

    END SUBROUTINE init_BEC_element
        
END MODULE m_types_BEC