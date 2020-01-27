      MODULE m_spratm
c------------------------------------------------------------------
c
c    This is the driver subroutine for a full-relativistic spin-polarized
c    core charge and spin density calculation using the 4 fold coupled
c    dirac equation : collection of references on this story see in
c    H. Ebert, J.Phys.: Condens. Matter 1 (1989) 9111.
c
c    Attention : The algorithm uses Ry-units. Therefore , potential
c                is multiplied by factor of 2 before it is used in the
c                dirac equation
c
c---> input:
c    Vr      =   spherical potential
c    Br      =   spherical magn. field
c    z       =   atomic charge
c    rnot    =   radial mesh starting point
c    dx      =   radial mesh logariphmic increment
c    jtop    =   upper bond for core radial mesh
c---> i/o
c    ectab   =   atomic energy levels for (\kappa,\mu) (in Hr)
c---> output
c    sume    =   sum of atomic eigenvalues (in Hr)
c    rhochr  =   core charge density
c    rhospn  =   core spin density
c
c........................................................ spratm
      CONTAINS
      SUBROUTINE spratm(
     >                  msh,vr,br,z,rnot,dx,jtop,ectab,ntab,ltab,
     <                  sume,rhochr,rhospn)
c
      USE m_core
      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: msh,jtop
      REAL,    INTENT (IN) :: dx,rnot,z
      REAL,    INTENT (OUT):: sume
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ntab(100),ltab(100)
      REAL,    INTENT (IN) :: br(msh),vr(msh)
      REAL,    INTENT (OUT):: rhochr(msh),rhospn(msh)
      REAL,    INTENT (INOUT):: ectab(100)
C     ..
C     .. Local Scalars ..
      REAL rr,stval
      INTEGER ic,ir,nshell,n_old,l_old
C     ..
C     .. Local Arrays ..
      REAL bt(msh),vt(msh)
      INTEGER nqntab(15),lqntab(15)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC exp,log
C     ..
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
        ELSEIF (ltab(ic).NE.l_old) THEN
           nshell = nshell + 1
           nqntab(nshell) = ntab(ic)
           lqntab(nshell) = ltab(ic)
           n_old = ntab(ic)
           l_old = ltab(ic)
        ENDIF
      ENDDO
c Hr -> Ry
      ic = 0
      DO ic = 1,100
         ectab(ic) = 2.*ectab(ic)
      END DO
c potential and field redefinition
      rr = rnot
      DO ir = 1,msh
         vt(ir) = 2.*vr(ir)/rr
         bt(ir) = 2.*br(ir)/rr
         rr = rr*exp(dx)
      END DO
      stval = log(rnot)
c
      CALL core(
     >          msh,vt,bt,z,stval,dx,nshell,nqntab,lqntab,jtop,
     X          ectab,
     <          rhochr,rhospn)

c Ry -> Hr
      sume = 0.0
      DO ic = 1,100
         ectab(ic) = ectab(ic)/2.
         sume = sume + ectab(ic)
      END DO

      END SUBROUTINE spratm
      END MODULE m_spratm
