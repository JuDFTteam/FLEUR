      MODULE m_strgndim
      use m_juDFT
      CONTAINS
      SUBROUTINE strgn1_dim(
     >                  gmax,bmat,invs,zrfs,mrot,tau,nop,nop2,
     <                  k1d,k2d,k3d,n3d,n2d,odd)
c
c     *********************************************************
c     generate two- and three-dimensional stars
c     for slab geometry
c     e. wimmer   nov.1984    c.l.fu  1987
c     *********************************************************
      USE m_boxdim
      USE m_spgrot
      USE m_types, ONLY : od_dim
      USE m_ifft, ONLY : i2357
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nop,nop2
      REAL,    INTENT (IN)  :: gmax
      INTEGER, INTENT (OUT) :: k1d,k2d,k3d,n3d,n2d

      INTEGER, INTENT (IN) :: mrot(3,3,nop)
      REAL,    INTENT (IN) :: tau(3,nop),bmat(3,3)
      LOGICAL, INTENT (IN) :: invs,zrfs
c
      INTEGER mx1,mx2,mx3
      INTEGER m
      INTEGER ng2,ng3,n,j,k,k1,k2,k3,m0,mxx1,mxx2
      REAL    arltv1,arltv2,arltv3,s
      INTEGER kv(3),kr(3,nop),inv_du(nop)
      REAL    g(3)
      INTEGER , ALLOCATABLE :: kv2(:,:)
c-odim
      TYPE (od_dim), INTENT (INOUT) :: odd
c+odim
      INTRINSIC int

C---> Determine Gmax box of size mx1, mx2, mx3,
c     for which |G(mx1,mx2,mx3)| < Gmax
c     arltv(i) length of reciprical lattice vector along direction (i)
C
      CALL boxdim(bmat,arltv1,arltv2,arltv3)

      mx1 = int(gmax/arltv1) + 1
      mx2 = int(gmax/arltv2) + 1
      mx3 = int(gmax/arltv3) + 1
      ALLOCATE ( kv2(2,(2*mx1+1)*(2*mx2+1)) )
      DO n=1,nop
         inv_du(n) = n ! dummy array for spgrot
      ENDDO
c
c     two-dimensional stars
c
      mxx1 = 0
      mxx2 = 0
      ng2 = 0
      kv(3) = 0
      DO 40 k1 = mx1,-mx1,-1
         kv(1) = k1
         DO 30 k2 = mx2,-mx2,-1
            kv(2) = k2
            DO j = 1,2
               g(j) = kv(1)*bmat(1,j) + kv(2)*bmat(2,j)
            END DO
            s = sqrt(g(1)**2+g(2)**2)
c
c--->   stars should be within the g_max-sphere !   (Oct.97) sbluegel
            IF (s.LE.gmax) THEN
c
               CALL spgrot(nop,.true.,mrot,tau,inv_du,
     >                     kv,
     <                     kr)
               DO n = 1,nop2
                  IF (mxx1.lt.kr(1,n)) mxx1 = kr(1,n)
                  IF (mxx2.lt.kr(2,n)) mxx2 = kr(2,n)
               ENDDO
               DO k = 1,ng2 
                  DO n = 1,nop2
                     IF (kr(1,n).EQ.kv2(1,k) .AND.
     +                   kr(2,n).EQ.kv2(2,k)) GOTO 30
                  ENDDO
               ENDDO
c--->    new representative found
               ng2 = ng2 + 1
c               IF (ng2.GT.n2dd)  CALL juDFT_error("n2dd",calledby="strgn_dim")
               DO j = 1,2
                  kv2(j,ng2) = kv(j)
               END DO
            END IF
   30    CONTINUE
   40 CONTINUE
c
      n2d = ng2
#ifdef CPP_AIX
      k1d = i2357(mxx1)
      k2d = i2357(mxx2)
      k3d = i2357(mx3)
#else
      k1d = mxx1
      k2d = mxx2
      k3d = mx3
#endif
      IF (odd%d1) THEN
        k1d = max(k1d,mx1)
        k2d = max(k2d,mx2)
      ENDIF
c
c     three dimensional stars
      ng3 = 0
      m0 = -mx3
      IF (zrfs .OR. invs) m0 = 0
      IF (odd%d1 .and. odd%zrfs) m0 = 0
      DO k2 = 1,ng2
         DO k3 = m0,mx3
            DO j = 1,3
               g(j) = kv2(1,k2)*bmat(1,j) + kv2(2,k2)*bmat(2,j) +
     +                k3*bmat(3,j)
            END DO
            s = sqrt(g(1)**2+g(2)**2+g(3)**2)
            IF (s.LE.gmax) THEN
               ng3 = ng3 + 1
            END IF
         END DO
      END DO
      n3d = ng3
c-odim    modifications due to one-dimensionality, YM
      IF (odd%d1) THEN
         j = 0
         m0 = k3d
         IF (odd%chi.NE.1 .AND. (odd%invs .OR. odd%zrfs))
     &        CALL juDFT_error
     +        ("chiral symmetries are not consistent with invs or zrfs"
     +        ,calledby ="strgn_dim")
         IF (odd%zrfs) m0 = 0
         IF (odd%chi.EQ.1) THEN
            odd%nop = odd%rot
            IF (odd%invs .OR. odd%zrfs) odd%nop = odd%nop*2
            DO k3 = -m0,k3d
               DO m = -odd%M,odd%M
                  IF (MOD(m,odd%rot).EQ.0) THEN
                     j = j+1
                  END IF
               END DO
            END DO
         ELSEIF (odd%chi.NE.1) THEN
            odd%nop = odd%chi
            DO k3 = -m0,k3d
               DO m = -odd%M,odd%M
                  IF (MOD(m+(odd%rot)*k3,odd%chi).EQ.0) THEN
                     j = j+1
                  END IF
               END DO
            END DO
         END IF
         odd%n2d = j
         odd%nq2 = j
         odd%kimax2 = odd%n2d - 1
      ELSE
         odd%n2d = n2d
         odd%nq2 = n2d
         odd%nop = nop
      END IF
c+odim
c     listing
      WRITE (6,8000) gmax,n3d,n2d,k1d,k2d,k3d
 8000 FORMAT (' gmax=',f10.6/' n3d=  ',i5/' n2d=  ',i5/' k1d=  ',
     +       i5/' k2d=  ',i5/' k3d=  ',i5/)

      DEALLOCATE (kv2)
      END SUBROUTINE strgn1_dim
!---------------------------------------------------------------------
      SUBROUTINE strgn2_dim(
     >                  gmax,bmat,invs,zrfs,mrot,tau,nop,
     <                  k1d,k2d,k3d,n3d,n2d)
c
c     *********************************************************
c     generate three-dimensional stars for bulk
c     based on  e. wimmer   nov.1984    c.l.fu  1987
c     *********************************************************
      USE m_boxdim
      USE m_spgrot
      USE m_ifft, ONLY : i2357
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nop
      REAL,    INTENT (IN)  :: gmax
      INTEGER, INTENT (OUT) :: k1d,k2d,k3d,n3d,n2d

      INTEGER, INTENT (IN) :: mrot(3,3,nop)
      REAL,    INTENT (IN) :: tau(3,nop),bmat(3,3)
      LOGICAL, INTENT (IN) :: invs,zrfs
c
      INTEGER mx1,mx2,mx3,ng3,n,j,k,k1,k2,k3,mxx1,mxx2,mxx3
      REAL arltv1,arltv2,arltv3,s,gmax2
      INTEGER kr(3,nop),kv(3),inv_du(nop)

      REAL    g(3)
      INTEGER , ALLOCATABLE :: kv3(:,:),ig(:,:,:)

      INTRINSIC int

C---> Determine Gmax box of size mx1, mx2, mx3,
c     for which |G(mx1,mx2,mx3)| < Gmax
c     arltv(i) length of reciprical lattice vector along direction (i)
C
      CALL boxdim(bmat,arltv1,arltv2,arltv3)

      mx1 = int(gmax/arltv1) + 1
      mx2 = int(gmax/arltv2) + 1
      mx3 = int(gmax/arltv3) + 1
      ALLOCATE ( kv3(3,(2*mx1+1)*(2*mx2+1)*(2*mx3+1)) )
      ALLOCATE ( ig(-mx1:mx1,-mx2:mx2,-mx3:mx3) )
      DO n=1,nop
         inv_du(n) = n ! dummy array for spgrot
      ENDDO
c
c     three dimensional stars
c
      mxx1 = 0
      mxx2 = 0
      mxx3 = 0
      ng3 = 0
      gmax2 = gmax * gmax
      ig = 0
      x_dim: DO k1 = mx1,-mx1,-1
        kv(1) = k1
        y_dim: DO k2 = mx2,-mx2,-1
          kv(2) = k2
          z_dim: DO k3 = mx3,-mx3,-1
            IF ( ig(k1,k2,k3) .NE. 0 ) CYCLE z_dim
            kv(3) = k3

            DO j = 1,3
               g(j) = kv(1)*bmat(1,j) + kv(2)*bmat(2,j) +
     +                kv(3)*bmat(3,j)
            END DO
            s = g(1)**2 + g(2)**2 + g(3)**2
c
c--->   stars should be within the g_max-sphere !   (Oct.97) sbluegel
            IF (s.LE.gmax2) THEN
c--->    new representative found
c
               CALL spgrot(nop,.true.,mrot,tau,inv_du,
     >                     kv,
     <                     kr)
               ng3 = ng3 + 1
               DO j = 1,3
                  kv3(j,ng3) = kv(j)
               END DO
               DO n = 1,nop
                  IF (mxx1.lt.kr(1,n)) mxx1 = kr(1,n)
                  IF (mxx2.lt.kr(2,n)) mxx2 = kr(2,n)
                  IF (mxx3.lt.kr(3,n)) mxx3 = kr(3,n)
                  ig(kr(1,n),kr(2,n),kr(3,n)) = ng3
               ENDDO
            END IF
          ENDDO z_dim
        ENDDO y_dim
      ENDDO x_dim
c
      n2d = 2
      n3d = ng3
#ifdef CPP_AIX
      k1d = i2357(mxx1)
      k2d = i2357(mxx2)
      k3d = i2357(mxx3)
      write(*,*) k1d,mxx1
#else
      k1d = mxx1
      k2d = mxx2
      k3d = mxx3
#endif
c
c     listing
      WRITE (6,8000) gmax,n3d,n2d,k1d,k2d,k3d
 8000 FORMAT (' gmax=',f10.6/' n3d=  ',i5/' n2d=  ',i5/' k1d=  ',
     +       i5/' k2d=  ',i5/' k3d=  ',i5/)

      DEALLOCATE (kv3,ig)
      END SUBROUTINE strgn2_dim
      END MODULE m_strgndim
