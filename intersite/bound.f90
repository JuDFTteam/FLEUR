MODULE m_bound
   !------------------------------------------------------------------------------
   !
   ! MODULE:  m_bound
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>   calculates the factor between u_l(r) and R_l(r) for energy E
   !>   where u_l(r) is the normed solution of the radial equation and
   !>   R_l(r) is matched to bessel functions at the mt-radius
   !
   ! REVISION HISTORY:
   ! February 2019 - Initial Version
   !------------------------------------------------------------------------------


   CONTAINS

   SUBROUTINE bound(e,rmt,alpha,u_rmt,w,lmax)

      USE m_SphBessel

      IMPLICIT NONE

      COMPLEX,                   INTENT(IN)  :: e
      REAL,                      INTENT(IN)  :: r_mt
      REAL,                      INTENT(IN)  :: u_rmt !value of u_rmt at the mt-radius
      REAL,                      INTENT(IN)  :: w !logarithmic derivative at mt-radius
      REAL, dimension(0:lmax-2), INTENT(OUT) :: alpha(:)
      INTEGER,                   INTENT(IN)  :: lmax


      CALL SphBesselComplex(jl,nl,hl,SQRT(e)*r_mt,lmax)

      DO l = 0, lmax-2

         tl = 1/sqrt(e) * ((l/r_mt-w)*jl(l) - SQRT(e)*jl(l+1))/((l/r_mt-w)*hl(l) - SQRT(e)*hl(l+1))
         alpha(l) = 1/u_rmt * (jl(l)-ImagUnit*SQRT(e)*hl(l)*tl)

      ENDDO

   END SUBROUTINE bound

END MODULE m_bound