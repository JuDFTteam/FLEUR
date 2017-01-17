MODULE m_slabgeom
  USE m_juDFT
CONTAINS
  SUBROUTINE slabgeom(atoms,cell,nsld,&
       nsl,zsl,nmtsl,nslat,volsl,volintsl)
    !***********************************************************************
    !     This subroutine calculates  z-coordinates of film layers, 
    !     a number of mt-pheres in each layer, and they typs.
    !                                   Yury Koroteev 2003-09-30
    !***********************************************************************
    !                     ABBREVIATIONS
    !
    ! natd                  : in, the number of atoms in the film
    ! pos(3,natd)           : in, the coordinates of atoms in the film
    ! ntypd,ntype           : in, the number of mt-sphere types
    ! z1                    : in, half the film thickness (0.5*D_tilde)
    ! neq(ntypd)            : in, the number of mt-spheres of the same type
    ! area                  : in, the area of the surface unit cell
    ! volmts(ntypd)         : in, the volume of mt-spheres
    ! nsld                  : in, the number of layers in the film
    !-----------------------------------------------------------------------
    ! nsl                   : in, the number of layers in the film
    ! zsl(2,nsld)           : out, z-coordinates of the layers
    ! nmtsl(ntypd,nsld)     : out, the number of mt-spheres of the ntypd-
    !                                type in the nsl-layer of the film
    ! nslat(natd,nsld)      : out, 
    ! volsl(nsld)           : out, the volume of film layers  
    ! volintsl(nsld)        : out, the volume of mt-spheres
    !
    !-----------------------------------------------------------------------
    ! znz(nsl)              : work, the z-ordinate of mt-spheres in 
    !                               the nsl-layer 
    !-----------------------------------------------------------------------
    !
    USE m_types
    IMPLICIT NONE
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     ..Scalar Argument
    INTEGER, INTENT  (IN) :: nsld
    INTEGER, INTENT (OUT) :: nsl
    !     ..
    !     ..Array Arguments
    INTEGER, INTENT (OUT) :: nmtsl(atoms%ntype,nsld),nslat(atoms%nat,nsld)
    REAL,    INTENT (OUT) :: zsl(2,nsld),volsl(nsld)  
    REAL,    INTENT (OUT) :: volintsl(nsld)
    !     ..
    !     ..Local Scalars 
    INTEGER  iz,i,j,na,isum,mt,n,nz
    REAL    epsz,half,zs,w,del,vmt
    !     ..
    !     ..Local Arrays 
    REAL    znz(nsld)
    !     ..
    !    ------------------------------------------------------------------
    DATA epsz/1.e-3/ half/0.5/
    !    ----------------------------------------------
    !
    ! --->  Calculate the number of the film layers (nsl)
    !

    znz(1) = atoms%pos(3,1)
    nz = 1
    na = 0
    DO  i=1,atoms%ntype
       DO  j=1,atoms%neq(i)
          na = na + 1
          zs = atoms%pos(3,na)
          
          IF(any(ABS(zs-znz(:nz)).LT.epsz)) CYCLE
          nz = nz+1
          znz(nz) = zs
       ENDDO
    ENDDO

    nsl = nz
    IF (nsl.GT.nsld) THEN
       WRITE(*,*) 'nsl =',nsl,' > nsld =',nsld
       CALL juDFT_error("nsl>nsld ",calledby ="slabgeom")
    ENDIF
    !
    ! ---> Order the film layers
    !
    DO  i=1,nsl
       DO  j=i,nsl
          IF(znz(j).LT.znz(i)) THEN
             w      = znz(i)
             znz(i) = znz(j)
             znz(j) = w
          ENDIF
       ENDDO
    ENDDO
    !
    ! ---> Construct the z-coordinates of the film layers ( zsl(2,nsl) )
    !
    zsl(1,1) = -cell%z1
    DO i=1,nsl-1
       zsl(2,i) = (znz(i) + znz(i+1)) * half
       zsl(1,i+1) = zsl(2,i)
    ENDDO
    zsl(2,nsl) = cell%z1
    ! 
    ! ---> Calculate a number of mt-spheres of the same type
    ! ---> (nmtsl) in each layer of the film
    !
    DO i=1,nsl
       del = ABS( zsl(2,i) - zsl(1,i) )
       volsl(i) = del*cell%area
       n = 0
       vmt = 0.0
       DO j=1,atoms%ntype
          isum = 0
          DO mt=1,atoms%neq(j)
             n = n + 1
             zs = atoms%pos(3,n)
             IF((zsl(1,i).LT.zs).AND.(zs.LT.zsl(2,i)))  THEN
                isum=isum+1
                nslat(n,i)=1
             ELSE
                nslat(n,i)=0
             ENDIF
          ENDDO
          nmtsl(j,i) = isum
          vmt = vmt + isum*atoms%volmts(j)
       ENDDO
       volintsl(i) = volsl(i) - vmt
    ENDDO
    !
  END SUBROUTINE slabgeom
END MODULE m_slabgeom

