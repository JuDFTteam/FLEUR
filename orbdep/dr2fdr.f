      MODULE m_dr2fdr
      CONTAINS

      SUBROUTINE dr2fdr(function,rmsh,jri,
     <                  deriv)
c
c     Construct r**2 * df(r)/dr ; input 'function' is on a mesh (rmsh)
c     with 'jri' points and is assumed to be multiplied by r**2.
c     difcub performs analytic derivative of Lagrangian of 3rd order.
c
      USE m_differentiate,ONLY:difcub
      IMPLICIT NONE

! Arguments ...

      INTEGER, INTENT (IN)  :: jri
      REAL,    INTENT (IN)  :: function(jri),rmsh(jri)
      REAL,    INTENT (OUT) :: deriv(jri)

! Locals ...

      INTEGER ir
      REAL faux(jri),xi


c
c take derivative of r**2 f(r): faux = d[r^2 f(r)]/dr
c first point
c
      xi = rmsh(1)
      faux(1) = difcub(rmsh(1),function(1),xi)
c
c 2nd to last-2 
c
      DO ir = 2, jri - 2
         xi = rmsh(ir)
         faux(ir) = difcub(rmsh(ir-1),function(ir-1),xi)
      END DO
c
c last-1
c
      ir = jri - 1
      xi = rmsh(ir)
      faux(ir) = difcub(rmsh(jri-3),function(jri-3),xi)
c
c last point
c
      ir = jri
      xi = rmsh(ir)
      faux(ir) = difcub(rmsh(jri-3),function(jri-3),xi)
c
c calculate r^2 df(r)/dr = d[r^2 f(r)]/dr - 2 r f(r) 
c
      DO ir = 1, jri
         deriv(ir) = faux(ir) - 2.0 * function(ir) / rmsh(ir)
      END DO
 
      RETURN
      END SUBROUTINE dr2fdr
      END MODULE m_dr2fdr
