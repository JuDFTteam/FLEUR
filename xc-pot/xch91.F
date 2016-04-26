      MODULE m_xch91
c.....-----------------------------------------------------------------
c     gga91 exchange for a spin-unpolarized electronic system
c.....-----------------------------------------------------------------
c     input d: density
c           s:  abs(grad d)/(2*kf*d)
c           u:  (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)
c           v: (laplacian d)/(d*(2*kf)**2)
c     output:  exchange energy (ex) and potential (vx) in ry.

c       kf=cbrt(3*pai**2*d).
c.....-----------------------------------------------------------------
      CONTAINS
      SUBROUTINE xch91(
     >                 d,s,u,v,
     <                 exl,exg,vxl,vxg)

      IMPLICIT NONE

      REAL, INTENT (IN)  :: d,s,u,v
      REAL, INTENT (OUT) :: exl,exg,vxl,vxg

      REAL a,a1,a2,a3,a4,ax,b1,ex,f,fac,fs,fss,
     +     p0,p1,p10,p11,p2,p3,p4,p5,p6,p7,p8,p9,s2,s3,s4,
     +     thrd,thrd4,vx,a4s2,expsm
c.....-----------------------------------------------------------------
c     ..
      DATA a1,a2,a3,a4/0.19645,0.27430e0,0.15084e0,100.e0/
      DATA ax,a,b1/-0.7385588e0,7.7956e0,0.004e0/
      DATA thrd,thrd4/0.333333333333e0,1.33333333333e0/
c.....-----------------------------------------------------------------
c     expsm: argument of exponential-smallest.
#ifdef CPP_DEC
      expsm=-300.e0
#else
      expsm=-500.e0
#endif

      fac = ax*d**thrd
      s2 = s*s
      s3 = s2*s
      s4 = s2*s2
      p0 = 1.e0/sqrt(1.e0+a*a*s2)
      p1 = log(a*s+1.e0/p0)
      a4s2=max(-a4*s2,expsm)
      p2 = exp(a4s2)
      p3 = 1.e0/ (1.e0+a1*s*p1+b1*s4)
      p4 = 1.e0 + a1*s*p1 + (a2-a3*p2)*s2
      f = p3*p4
      ex = fac*f
c  local exchange exl
      exl = fac
      exg = ex - exl

c  energy done. now the potential:
      p5 = b1*s2 - (a2-a3*p2)
      p6 = a1*s* (p1+a*s*p0)
      p7 = 2.e0* (a2-a3*p2) + 2.e0*a3*a4*s2*p2 - 4.e0*b1*s2*f
      fs = p3* (p3*p5*p6+p7)
      p8 = 2.e0*s* (b1-a3*a4*p2)
      p9 = a1*p1 + a*a1*s*p0* (3.e0-a*a*s2*p0*p0)
      p10 = 4.e0*a3*a4*s*p2* (2.e0-a4*s2) - 8.e0*b1*s*f - 4.e0*b1*s3*fs
      p11 = -p3*p3* (a1*p1+a*a1*s*p0+4.e0*b1*s3)
      fss = p3*p3* (p5*p9+p6*p8) + 2.e0*p3*p5*p6*p11 + p3*p10 + p7*p11
      vx = fac* (thrd4*f- (u-thrd4*s3)*fss-v*fs)
c  local exchange vxl:
      vxl = fac*thrd4
      vxg = vx - vxl

c in ry and energy density.
      exl = exl*2.e0*d
      exg = exg*2.e0*d
      vxl = vxl*2.
      vxg = vxg*2.

      END SUBROUTINE xch91
      END MODULE m_xch91
