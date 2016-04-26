MODULE m_slabdim
  USE m_juDFT
CONTAINS
  SUBROUTINE slab_dim(   atoms,nsld)
    !***********************************************************************
    !     This subroutine calculates  the number of layers in the slab 
    !     
    !                                   Yury Koroteev 2003-09-30
    !***********************************************************************
    !                     ABBREVIATIONS
    !
    ! natd                  : in, the number of atoms in the film
    ! pos(3,natd)         : in, the coordinates of atoms in the film
    ! ntypd,ntype           : in, the number of mt-sphere types
    ! neq(ntypd)            : in, the number of mt-spheres of the same type
    !-----------------------------------------------------------------------
    ! nsld                  : out, the number of layers in the film
    !-----------------------------------------------------------------------
    ! znz(nsl)              : work, the z-ordinate of mt-spheres in 
    !                               the nsl-layer 
    !-----------------------------------------------------------------------
    !
    USE m_types
    IMPLICIT NONE

    TYPE(t_atoms),INTENT(IN)   :: atoms
    !	..
    !       ..Scalar Argument
    INTEGER, INTENT (OUT) :: nsld
    !       ..
    !       ..Array Arguments
    !       ..
    !       ..Local Scalars 
    INTEGER  iz,i,j,na,nz
    REAL    zs
    !       ..
    !       ..Local Arrays 
    REAL    znz(atoms%natd)
    !       ..
    !    ----------------------------------------------
    REAL,PARAMETER:: epsz=1.e-3 
    !    ----------------------------------------------
    !
    ! --->  Calculate the number of the film layers (nsld)
    !
    znz(1) = atoms%pos(3,1)
    nz = 1
    na = 0
    DO i=1,atoms%ntype
       DO j=1,atoms%neq(i)
          na = na + 1
          zs = atoms%pos(3,na)
          DO iz=1,nz
             IF(ABS(zs-znz(iz)).LT.epsz) CYCLE
          ENDDO
          nz = nz+1
          znz(nz) = zs
       ENDDO
    ENDDO
    nsld = nz
    IF(nsld>atoms%natd)   CALL juDFT_error("nsld.GT.atoms%natd ",calledby="slab_dim")
    !
  END SUBROUTINE slab_dim
END MODULE m_slabdim

