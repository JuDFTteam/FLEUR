      MODULE m_corg91
c.....-----------------------------------------------------------------
c     gga91 correlation
c.....-----------------------------------------------------------------
      CONTAINS
      SUBROUTINE corg91(fk,sk,gz,ec,ecrs,eczta,rs,zta,t,uu,vv,ww,h,
     +                  dvcup,dvcdn)
c.....-----------------------------------------------------------------
c     input
c           rs: seitz radius
c         zta: relative spin polarization
c            t: abs(grad d)/(d*2.*ks*gz)
c           uu: (grad d)*grad(abs(grad d))/(d**2 * (2*ks*gz)**3)
c           vv: (laplacian d)/(d * (2*ks*gz)**2)
c           ww: (grad d)*(gradzta)/(d * (2*ks*gz)**2
c      output
c                  h: nonlocal part of correlation energy per electron
c          dvcup,-dn: nonlocal parts of correlation potentials.

c      with ks=sqrt(4*kf/pai), gz=[(1+zta)**(2/3)+(1-zta)**(2/3)]/2, &
c           kf=cbrt(3*pai**2*d).
c.....-----------------------------------------------------------------
c.....-----------------------------------------------------------------
c     .. previously untyped names ..
      IMPLICIT NONE

      REAL dvcdn,dvcup,ec,ecrs,eczta,fk,gz,h,rs,sk,t,uu,vv,ww,zta,a4,
     +     alf,argmx,b,b2,b2fac,bec,bet,bg,c1,c2,c3,c4,c5,c6,cc,cc0,
     +     ccrs,coeff,comm,cx,delt,fact0,fact1,fact2,fact3,fact4,fact5,
     +     gm,gz3,gz4,h0,h0b,h0bt,h0rs,h0rst,h0t,h0tt,h0z,h0zt,h1,h1rs,
     +     h1rst,h1t,h1tt,h1z,h1zt,hrs,hrst,ht,htt,hz,hzt,pon,pref,q4,
     +     q5,q6,q7,q8,q9,r0,r1,r2,r3,r4,rs2,rs3,rsthrd,t2,t4,t6,thrd2,
     +     thrdm,xnu,r1t2,expsm
c     ..
      DATA xnu,cc0,cx,alf/15.75592e0,0.00423500,-0.001667212e0,0.0900/
      DATA c1,c2,c3,c4/0.00256800,0.02326600,7.389e-6,8.723e0/
      DATA c5,c6,a4/0.472e0,7.389e-2,100.00/
      DATA thrdm,thrd2/-0.333333333333e0,0.666666666667e0/
c.....-----------------------------------------------------------------
c     expsm: argument of exponential-smallest.
#ifdef CPP_DEC
      expsm=-300.
#else
      expsm=-500.
#endif
      argmx = 174.0
      bet = xnu*cc0
      delt = 2.e0*alf/bet
      gz3 = gz**3
      gz4 = gz3*gz
      pon = -delt*ec/ (gz3*bet)
      IF (pon.gt.argmx) THEN
          b = 0.
      ELSE
          b = delt/ (exp(pon)-1.00)
      ENDIF
      b2 = b*b
      t2 = t*t
      t4 = t2*t2
      t6 = t4*t2
      rs2 = rs*rs
      rs3 = rs2*rs
      q4 = 1.00 + b*t2
      q5 = 1.00 + b*t2 + b2*t4
      q6 = c1 + c2*rs + c3*rs2
      q7 = 1.00 + c4*rs + c5*rs2 + c6*rs3
      cc = -cx + q6/q7
      r0 = (sk/fk)**2
      r1 = a4*r0*gz4
      coeff = cc - cc0 - 3.e0*cx/7.00
      r2 = xnu*coeff*gz3
c+tagu
      r1t2=max(-r1*t2,expsm)
      r3 = exp(r1t2)
c-tagu
      h0 = gz3* (bet/delt)*log(1.00+delt*q4*t2/q5)
      h1 = r3*r2*t2
      h = h0 + h1
c  local correlation option:
c     h=0.0e0

c  energy done. now the potential:

      ccrs = (c2+2.*c3*rs)/q7 - q6* (c4+2.*c5*rs+3.*c6*rs2)/q7**2
      rsthrd = rs/3.e0
      r4 = rsthrd*ccrs/coeff
      gm = ((1.00+zta)**thrdm- (1.00-zta)**thrdm)/3.00
      IF (pon.gt.argmx) THEN
          b2fac = 0.
      ELSE
          b2fac = b2* (delt/b+1.00)
      ENDIF
      bg = -3.e0*ec*b2fac/ (bet*gz4)
      bec = b2fac/ (bet*gz3)
      q8 = q5*q5 + delt*q4*q5*t2
      q9 = 1.00 + 2.00*b*t2
      h0b = -bet*gz3*b*t6* (2.e0+b*t2)/q8
      h0rs = -rsthrd*h0b*bec*ecrs
      fact0 = 2.e0*delt - 6.00*b
      fact1 = q5*q9 + q4*q9*q9
      h0bt = 2.e0*bet*gz3*t4* ((q4*q5*fact0-delt*fact1)/q8)/q8
      h0rst = rsthrd*t2*h0bt*bec*ecrs
      h0z = 3.e0*gm*h0/gz + h0b* (bg*gm+bec*eczta)
      h0t = 2.*bet*gz3*q9/q8
      h0zt = 3.e0*gm*h0t/gz + h0bt* (bg*gm+bec*eczta)
      fact2 = q4*q5 + b*t2* (q4*q9+q5)
      fact3 = 2.e0*b*q5*q9 + delt*fact2
      h0tt = 4.e0*bet*gz3*t* (2.e0*b/q8- (q9* (fact3/q8))/q8)
      h1rs = r3*r2*t2* (-r4+r1*t2/3.e0)
      fact4 = 2.e0 - r1*t2
      h1rst = r3*r2*t2* (2.e0*r4* (1.00-r1*t2)-thrd2*r1*t2*fact4)
      h1z = gm*r3*r2*t2* (3.e0-4.e0*r1*t2)/gz
      h1t = 2.e0*r3*r2* (1.00-r1*t2)
      h1zt = 2.e0*gm*r3*r2* (3.e0-11.00*r1*t2+4.e0*r1*r1*t4)/gz
      h1tt = 4.e0*r3*r2*r1*t* (-2.e0+r1*t2)
      hrs = h0rs + h1rs
      hrst = h0rst + h1rst
      ht = h0t + h1t
      htt = h0tt + h1tt
      hz = h0z + h1z
      hzt = h0zt + h1zt
      comm = h + hrs + hrst + t2*ht/6.00 + 7.00*t2*t*htt/6.00
      pref = hz - gm*t2*ht/gz
      fact5 = gm* (2.e0*ht+t*htt)/gz
      comm = comm - pref*zta - uu*htt - vv*ht - ww* (hzt-fact5)
      dvcup = comm + pref
      dvcdn = comm - pref

c  local correlation option:
c     dvcup=0.0e0
c     dvcdn=0.0e0

      END SUBROUTINE corg91
      END MODULE m_corg91
