      MODULE m_od_phasy
      CONTAINS
      SUBROUTINE od_phasy(
     >                  ntypd,n3d,natd,lmaxd,ntype,neq,lmax,
     >                  taual,bmat,kv3,k,odi,ods,
     <                  pylm)
c*********************************************************************
c calculates 4pi*i**l/nop(3)*sum(R){exp(iRG(taual-taur)*conjg(ylm(RG)) }
c but for chiral groups of symmetries, as in phasy1.F    
c     Y.Mokrousov   august,2003
c ********************************************************************
      USE m_constants
      USE m_ylm
      USE m_od_chirot
      USE m_types, ONLY : od_inp, od_sym

      IMPLICIT NONE
c
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,n3d,natd,lmaxd
      INTEGER, INTENT (IN) :: ntype,k
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),kv3(3,n3d)
      REAL,    INTENT (IN) :: bmat(3,3),taual(3,natd)
      COMPLEX, INTENT (OUT):: pylm( (lmaxd+1)**2, ntypd )
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      COMPLEX sf,ci
      REAL x
      INTEGER j,l,m,n,na,lm
C     ..
C     .. Local Arrays ..
      COMPLEX ciall(0:lmaxd),ylm( (lmaxd+1)**2 )
      REAL rg(3)
      COMPLEX phas(ods%nop)
      REAL kr(3,ods%nop)

C     ..
C     .. Intrinsic Functions ..
      INTRINSIC cmplx,conjg,cos,sin
C     ..

      ci = cmplx(0.0,1.0)
      ciall(0) = fpi_const/ods%nop

      DO 10 l = 1,lmaxd
         ciall(l) = ciall(0)*ci**l
   10 CONTINUE
      na = 1
      DO 70 n = 1,ntype
         DO lm = 1, (lmax(n)+1)**2
               pylm(lm,n) = cmplx(0.,0.)
         ENDDO
         CALL od_chirot(odi,ods,bmat,kv3(1,k),kr,phas)         
         DO 60 j = 1,ods%nop
            rg(1) = kr(1,j)*bmat(1,1) + kr(2,j)*bmat(2,1) +
     +              kr(3,j)*bmat(3,1)
            rg(2) = kr(1,j)*bmat(1,2) + kr(2,j)*bmat(2,2) +
     +              kr(3,j)*bmat(3,2)
            rg(3) = kr(1,j)*bmat(1,3) + kr(2,j)*bmat(2,3) +
     +              kr(3,j)*bmat(3,3)            
            CALL ylm4(
     >                lmax(n),rg,
     <                ylm)
            x = tpi_const* (kr(1,j)*taual(1,na) + kr(2,j)*taual(2,na) +
     +                                        kr(3,j)*taual(3,na))
            sf = cmplx(cos(x),sin(x))*phas(j)
      
            DO l = 0,lmax(n)
               DO m = -l,l
                  lm = l*(l+1) + m + 1 
                  pylm(lm,n) = pylm(lm,n) +
     +                          ciall(l)*sf*conjg(ylm(lm))
               ENDDO
            ENDDO
   60    CONTINUE
         na = na + neq(n)
   70 CONTINUE

      RETURN
      END SUBROUTINE od_phasy
      END MODULE m_od_phasy
