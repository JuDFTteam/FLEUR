      MODULE m_outint
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: outint
      CONTAINS
      SUBROUTINE outint(
     >                  msh,e,fkap,cs,cis,s,vr,z,rn,rnot,h,d,
     <                  a0,b0,a,b,
     <                  ki,nodes)
c
c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
c-----x  perform the outward integration.   this routine is a     x----
c-----x  derivative of start1/diff1.   much unnecessary indexing  x----
c-----x  has been removed among other things.    dale koelling    x----

c     with adams' predictor and corrector. more stable than milne's
c       for GGA.
cc      Handbook of math.func. p.896. T.Asada Feb.20,1998.  
c-----x----x----x----x----x----x----x----x----x----x----x----x----x----

      IMPLICIT NONE

c     .. Arguments ..
      INTEGER,INTENT(IN)  :: msh
      INTEGER,INTENT(OUT) :: ki,nodes
      REAL   ,INTENT(IN)  :: e,fkap,cs,cis,s,z,rn,rnot,h,d
      REAL   ,INTENT(IN)  :: vr(msh)
      REAL   ,INTENT(OUT) :: a0(5),b0(5),a(msh),b(msh)
c     ..
c     .. local scalars ..
      REAL astr,atk,bstr,btk,da1,da2,da3,da4,db1,db2,db3,db4,g,gm,
     +     gp,h3,q11,q22,r,rp12,rp21
      REAL a1s,b1s,r21s,r210,r12s,r120,det1
      INTEGER i,k,kim,n
      REAL :: da(5),db(5)
c     ..
c     .. intrinsic functions ..
      INTRINSIC abs,sign,sqrt
c     ..
c.....------------------------------------------------------------------
c     .. equivalences ..
      EQUIVALENCE (da1,da(1)),(da2,da(2)),(da3,da(3)),(da4,da(4)),
     &  (da5,da(5))
      EQUIVALENCE (db1,db(1)),(db2,db(2)),(db3,db(3)),(db4,db(4)),
     &  (db5,db(5))

      REAL t18,t58,t98,t198,t378,t558,t598, da5,db5
c.....------------------------------------------------------------------

      h3 = h/3.0e0
      n = msh
cta+
      t558=55.e0/8.e0
      t598=59.e0/8.e0
      t378=37.e0/8.e0
      t98=9.e0/8.e0

      t198=19.e0/8.e0
      t58=5.e0/8.e0
      t18=1.e0/8.e0
cta-
      g = sqrt(fkap**2- (cis*z)**2)
      gp = g + s*fkap
      gm = g - s*fkap
      q11 = -gm
      q22 = -gp
c     first-order Frobenius start, a = astr + a1s*r, b = bstr + b1s*r,
c     with r*v = -z + v1*r near the nucleus
      astr = sqrt(abs(gp))
      bstr = sqrt(abs(gm))
      r210 = cis*z
      r21s = (cis* (e*rnot-vr(1))-r210)/rnot
      r12s = -r21s - 2*cs
      r120 = -r210
      det1 = (1.0-q11)* (1.0-q22) - r120*r210
      a1s = ((1.0-q22)*r12s*bstr+r120*r21s*astr)/det1
      b1s = ((1.0-q11)*r21s*astr+r210*r12s*bstr)/det1
      r = rnot
      do 10 i = 1,5
         a(i) = astr + a1s*r
         b(i) = bstr + b1s*r
         rp21 = cis* (e*r-vr(i))
         rp12 = -rp21 - 2*cs*r
         a0(i) = h* (q11*a(i)+rp12*b(i))
         b0(i) = h* (q22*b(i)+rp21*a(i))
         r = r*d
   10 continue
      r = rnot
      do 40 i = 1,5
         rp21 = cis* (e*r-vr(i))
         rp12 = -rp21 - 2*cs*r
         da(i) = h3* (q11*a(i)+rp12*b(i))
         db(i) = h3* (q22*b(i)+rp21*a(i))
         r = r*d
   40 continue
c     end of starting procedure

      nodes = 0
      r = rn

c..... find the classical turning point........

      ki = n
   50 ki = ki - 1
      if(ki.le.10) go to 60
      r = r/d
      if(e*r.lt.vr(ki)) go to 50
   60 ki = ki + 1
      if(ki.ge.n) ki = ki - 1
      kim = ki + 1
      if(kim.gt.n) kim = n
      r = rnot* (d**3)
      k = 4
   70 k = k + 1
      r = r*d
      rp21 = cis* (e*r-vr(k))
      rp12 = -rp21 - 2.0e0*cs*r

c     adams' extrapolation formula for predictor.
      atk=a(k-1)+ t558*da4-t598*da3+t378*da2-t98*da1
      btk=b(k-1)+ t558*db4-t598*db3+t378*db2-t98*db1

      da1 = da2
      da2 = da3
      da3 = da4

      db1 = db2
      db2 = db3
      db3 = db4

      da4 = h3* (q11*atk+rp12*btk)
      db4 = h3* (q22*btk+rp21*atk)

c     adams interpolation formula for corrector.
      a(k)=a(k-1)+ t98*da4+t198*da3-t58*da2+t18*da1
      b(k)=b(k-1)+ t98*db4+t198*db3-t58*db2+t18*db1

      da4=h3*(q11*a(k)+rp12*b(k))
      db4=h3*(q22*b(k)+rp21*a(k))

      if(k.ge.kim) go to 80

      nodes = nodes + 0.501e0*abs(sign(1.0e0,a(k))-sign(1.0e0,a(k-1)))
      go to 70

   80 continue

      RETURN
      END SUBROUTINE
      END
