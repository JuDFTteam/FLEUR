MODULE m_spratm

!------------------------------------------------------------------
!
!    This is the driver subroutine for a full-relativistic spin-polarized
!    core charge and spin density calculation using the 4 fold coupled
!    dirac equation : collection of references on this story see in
!    H. Ebert, J.Phys.: Condens. Matter 1 (1989) 9111.
!
!    Attention : The algorithm uses Ry-units. Therefore , potential
!                is multiplied by factor of 2 before it is used in the
!                dirac equation
!
!---> input:
!    Vr      =   spherical potential
!    Br      =   spherical magn. field
!    z       =   atomic charge
!    rnot    =   radial mesh starting point
!    dx      =   radial mesh logariphmic increment
!    jtop    =   upper bond for core radial mesh
!---> i/o
!    ectab   =   atomic energy levels for (\kappa,\mu) (in Hr)
!---> output
!    sume    =   sum of atomic eigenvalues (in Hr)
!    rhochr  =   core charge density
!    rhospn  =   core spin density
!
!........................................................ spratm

CONTAINS

   SUBROUTINE spratm(msh,vr,br,z,rnot,dx,jtop,ectab,ntab,ltab,sume,rhochr,rhospn)

      USE m_core

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: msh,jtop
      REAL,    INTENT (IN) :: dx,rnot,z
      REAL,    INTENT (OUT):: sume
      INTEGER, INTENT (IN) :: ntab(100),ltab(100)
      REAL,    INTENT (IN) :: br(msh),vr(msh)
      REAL,    INTENT (OUT):: rhochr(msh),rhospn(msh)
      REAL,    INTENT (INOUT):: ectab(100)

      REAL rr,stval
      INTEGER ic,ir,nshell,n_old,l_old

      REAL bt(msh),vt(msh)
      INTEGER nqntab(15),lqntab(15)

!      INTRINSIC exp,log

      nshell = 0
      ic = 0 ; n_old = -1 ; l_old = -1
      DO WHILE (ntab(ic+1).GT.0) 
         ic = ic + 1
         IF  (ntab(ic).NE.n_old) THEN
            nshell = nshell + 1
            nqntab(nshell) = ntab(ic)
            lqntab(nshell) = ltab(ic)
            n_old = ntab(ic)
            l_old = ltab(ic)
         ELSE IF (ltab(ic).NE.l_old) THEN
            nshell = nshell + 1
            nqntab(nshell) = ntab(ic)
            lqntab(nshell) = ltab(ic)
            n_old = ntab(ic)
            l_old = ltab(ic)
         END IF
      END DO
      
      ! Hr -> Ry
      ic = 0
      DO ic = 1, 100
         ectab(ic) = 2.0 * ectab(ic)
      END DO
      
      ! potential and field redefinition
      rr = rnot
      DO ir = 1,msh
         vt(ir) = 2.*vr(ir)/rr
         bt(ir) = 2.*br(ir)/rr
         rr = rr*exp(dx)
      END DO
      stval = log(rnot)

      CALL core(msh,vt,bt,z,stval,dx,nshell,nqntab,lqntab,jtop,ectab,rhochr,rhospn)

      ! Ry -> Hr
      sume = 0.0
      DO ic = 1,100
         ectab(ic) = ectab(ic)/2.
         sume = sume + ectab(ic)
      END DO

   END SUBROUTINE spratm
END MODULE m_spratm
