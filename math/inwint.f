      MODULE m_inwint
      CONTAINS
      SUBROUTINE inwint(
     >                  e,fl,ki,fkap,cs,cis,s,z,h,dd,rn,rnot,msh,vr,
     X                  a,b,
     <                  kj)
c
c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
c-----x  perform the  inward integration.   this routine is a     x----
c-----x  derivative of start2/diff2.   much unnecessary indexing  x----
c-----x  has been removed among other things.    dale koelling    x----
c
c     adams' procedure. more stable for GGA. Feb.20,1998. T.A.
c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
      use m_juDFT
      IMPLICIT NONE

c     .. Arguments ..
      INTEGER,INTENT (IN)  :: msh,ki
      INTEGER,INTENT (OUT) :: kj
      REAL,   INTENT (IN)  :: e,fl,fkap,cs,cis,s,z,h,dd,rn,rnot
      REAL,   INTENT (IN)  :: vr(msh)
c     ..
      REAL,   INTENT (INOUT) :: a(msh),b(msh)
c     ..
c     .. Local Scalars ..
      REAL atk,btk,d,da1,da2,da3,da4,db1,db2,db3,db4,eabsv,eps,f1,f2,f3,
     +     f4,f5,fg1,fg2,fg3,g,g1,g2,g3,g4,g5,h3,q11,q22,r,rp12,rp21
      INTEGER I,k,l,n
      REAL :: da(5),db(5)
c     ..
c     .. intrinsic functions ..
      INTRINSIC abs,mod,sqrt
c     ..
c     .. equivalences ..
cta+
      EQUIVALENCE (da1,da(1)), (da2,da(2)), (da3,da(3)), (da4,da(4)),
     &  (da5,da(5))
      EQUIVALENCE (db1,db(1)), (db2,db(2)), (db3,db(3)), (db4,db(4)),
     &  (db5,db(5))

      REAL c1,c5,c9,c19, t9,t37,t55,t59, da5,db5
cta-
c     ..
c     .. data statements ..
c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
c-----x this version modified to accept positive energies.        x----
c-----x                              d.d.koelling   7/7/77      x----x
c-----x----x----x----x----x----x----x----x----x----x----x----x----x----
c---* eps    value used to determine the practical infinity
      DATA f1/1.045833333e0/,f2/2.691666666e0/,f3/1.1e0/
      DATA f4/.4416666666e0/,f5/.07916666666e0/
      DATA g1/.9666666666e0/,g2/4.133333333e0/,g3/.8e0/
      DATA g4/.1333333333e0/,g5/.03333333333e0/
      DATA fg1/.9333333333e0/,fg2/4.266666666e0/,fg3/1.6e0/
      DATA eps/75.0e0/
c     ..
      h3 = -h/3.0e0
      n = msh
cta+
      t55=55.e0/8.e0
      t59=59.e0/8.e0
      t37=37.e0/8.e0
      t9=9.e0/8.e0

      c19=19.e0/8.e0
      c9=9.e0/8.e0
      c5=5.e0/8.e0
      c1=1.e0/8.e0
cta-

c     e: energy
c     ki: the grid point corresponding to classical terning point
cc        at which potential exceeded the energy.
cc        ki=249 for Fe and l=0.
c     kj: undefined at the beginning. finally the farthest grid 
cc        point which can be considered as infinity for the 
cc        level to be obtained.

      d = 1.0e0/dd
      g = sqrt(fkap**2- (cis*z)**2)
      q11 = s*fkap - g
      q22 = -g - s*fkap
      kj = ki + 10
      r = rnot* (dd**kj)

   10 kj = kj + 1
      IF (kj>msh)  CALL juDFT_error("inward integration failed",calledby
     +     ="inwint",hint
     +     ="This often indicates that the potential is unphysical")

c     eps is a big energy. thus if the next statement is fulfilled,
cc    the grid is considered to be infinitely far.

      IF (r* (vr(kj)-e*r).gt.eps) GOTO 20
      r = r*dd
      IF (kj.lt.n) go to 10
      b(kj) = mod(fl,2.0e0)
      a(kj) = 1.0e0 - b(kj)
      go to 30
   20 a(kj) = 1.0e0
      eabsv = abs(e)
      b(kj) = s*sqrt(eabsv/ (2*cs*cs+e))
   30 continue
c  30 do 40 l = 1,4

c     kj=338 for Fe l=0 and e=-333.
c     a(kj)=1.0

      do 40 l = 1,4
         k = kj - l
         a(k) = a(kj)
         b(k) = b(kj)
   40 continue

c     a(k)=1.0 for k=kj-4,..,kj.

      do 70 i = 1,6
         k = kj + 1
         r = rn*d** (n-k)
         do 60 l = 1,5
            k = k - 1
            r = r*d
            rp21 = cis* (e*r-vr(k))
            rp12 = -rp21 - 2.0e0*cs*r
            da(l) = h3* (q11*a(k)+rp12*b(k))
            db(l) = h3* (q22*b(k)+rp21*a(k))
   60    continue
         l = kj - 1
         a(l) = a(kj) + f1*da1 + f2*da2 + f4*da4 - (f3*da3+f5*da(5))
         b(l) = b(kj) + f1*db1 + f2*db2 + f4*db4 - (f3*db3+f5*db(5))
         l = l - 1
         a(l) = a(kj) + g1*da1 + g2*da2 + g3*da3 + g4*da4 - g5*da(5)
         b(l) = b(kj) + g1*db1 + g2*db2 + g3*db3 + g4*db4 - g5*db(5)
         l = l - 1
         a(l) = a(kj) + 1.0125e0*da1 + 3.825e0*da2 + 2.7e0*da3 +
     +          1.575e0*da4 - .1125e0*da(5)
         b(l) = b(kj) + 1.0125e0*db1 + 3.825e0*db2 + 2.7e0*db3 +
     +          1.575e0*db4 - .1125e0*db(5)
         l = l - 1
         a(l) = a(kj) + fg1*da1 + fg2*da2 + fg3*da3 + fg2*da4 +
     +          fg1*da(5)
         b(l) = b(kj) + fg1*db1 + fg2*db2 + fg3*db3 + fg2*db4 +
     +          fg1*db(5)
   70 continue

c     a(kj)=1.0, a(kj-1),..,a(kj-4) have been obtained.

c     adams' modification here on.

      k = kj - 3
      r = rn*d** (n-k)
   80 k = k - 1
      r = r*d
      rp21 = cis* (e*r-vr(k))
      rp12 = -rp21 - 2.0e0*cs*r

cta+
c     adams' predictor.
      atk=a(k+1) + t55*da4-t59*da3+t37*da2-t9*da1
      btk=b(k+1) + t55*db4-t59*db3+t37*db2-t9*db1

      da1 = da2
      da2 = da3
      da3 = da4

      db1 = db2
      db2 = db3
      db3 = db4

cta-
      da4=h3*(q11*atk+rp12*btk)
      db4=h3*(q22*btk+rp21*atk)

cta+
c     adams interpolation formula 25.5.5 (Handbook of math.functions).
      a(k)=a(k+1) + c9*da4+c19*da3-c5*da2+c1*da1
      b(k)=b(k+1) + c9*db4+c19*db3-c5*db2+c1*db1
cta-

      da4 = h3* (q11*a(k)+rp12*b(k))
      db4 = h3* (q22*b(k)+rp21*a(k))

      if(k.gt.ki) go to 80

      return
      END SUBROUTINE
      end


