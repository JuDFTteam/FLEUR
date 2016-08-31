      MODULE m_lhcal
      use m_juDFT

!*********************************************************************
!     determines the lattice harmonics for the given local
!     symmetry defined in terms of the rotation matrices.
!
!     input:  lmax     max l for lattice harmonics
!             nrot     number of operations
!             orth     rotation matrices (cartesian coordinates)
!                          ( 1,1 1,2 1,3 )( x )   ( x' )
!                          ( 2,1 2,2 2,3 )( y ) = ( y' )
!                          ( 3,1 3,2 3,3 )( z )   ( z' )
!             memd     max. members/lattice harmonic dimensioned for
!             nlhd     max. number of lattice harmonics dimensioned for
!
!     output: nlh      number of harmonic from l=0 to l=lmax
!             lnu      l value of harmonic
!             mem      number of members in lattice harmonic
!             lmnu     lm index of each member of lattice harmonic
!             c        coefficients of lattice harmonic
!
!            m. weinert                           1989/98
!*********************************************************************

      CONTAINS 
      SUBROUTINE lhcal(
     >                 memd,nlhd,lmax,nrot,orth,
     <                 nlh,lnu,mem,lmnu,c)
!DEC$ NOOPTIMIZE
      USE m_gaussp
      USE m_gtest
      USE m_ylm
      IMPLICIT NONE

!---> Arguments
      INTEGER, INTENT (IN) :: memd,nlhd
      INTEGER, INTENT (IN) :: lmax           !--->    max. l to consider
      INTEGER, INTENT (IN) :: nrot           !--->    number of
      REAL,    INTENT (IN) :: orth(3,3,nrot) !--->    rotation matrices

      INTEGER, INTENT(OUT) :: nlh
      INTEGER, INTENT(OUT) :: lnu((lmax+1)**2),mem((lmax+1)**2)
      INTEGER, INTENT(OUT) :: lmnu(memd,(lmax+1)**2)
      COMPLEX, INTENT(OUT) ::    c(memd,(lmax+1)**2)

!---> Locals
      INTEGER :: j,l,lh,lm,lm0,lmmax,l2,m,mems,mp,n,nn
      REAL    :: s
      REAL    :: v(3),vr(3),ovlp(0:2*lmax)
      COMPLEX :: a(0:lmax,0:2*lmax,0:lmax)
      COMPLEX, DIMENSION((LMAX+1)**2) :: ylm,ylmr,ylms

      INTEGER :: ngpts                                           ! gaussian
      REAL    :: vgauss(3,(2*lmax+1)*(lmax+1 + mod(lmax+1,2)))   ! integration
      REAL    :: wt((2*lmax+1)*( lmax+1 + mod(lmax+1,2) ))       ! points

      REAL, PARAMETER :: del = 1.e-5 ! parameter for machine roundoff, etc
!
!--->    generate gaussian points and test
!
!$OMP MASTER
      ngpts = (2*lmax+1)*( lmax+1 + mod(lmax+1,2) )
      CALL gaussp(
     >            lmax,
     <            vgauss,wt)
      CALL gtest(
     >           lmax,ngpts,vgauss,wt)

!--->    initialize
      lmmax = (lmax+1)**2
      a = cmplx(0.0,0.0)
!
!===>    loop over gaussian integration points
!
      DO nn = 1, ngpts
         v(:) = vgauss(:,nn)
         CALL ylm4(
     >             lmax,v,
     <             ylm)
         ylms(1:lmmax) = ylm(1:lmmax)
!--->    apply rotations
         DO n =2, nrot
            vr = matmul( orth(:,:,n), v )
            CALL ylm4(
     >                lmax,vr,
     <                ylmr)
            ylms(:) = ylms(:) + ylmr(:)
         ENDDO
         ylms = ylms/nrot
!--->    obtain coefficients
         DO l = 0, lmax
            lm0 = l*(l+1)+1
            DO mp = 0,l
               a(mp,0,l) = a(mp,0,l) +
     +                  wt(nn)*conjg(ylm(lm0+mp))*real(ylms(lm0))
            ENDDO
            DO m = 1, l
               lm = lm0+m
               DO mp = 0, l
                  a(mp,m,l)   = a(mp,m,l) +
     +                  wt(nn)*conjg(ylm(lm0+mp))* real(ylms(lm))
                  a(mp,m+l,l) = a(mp,m+l,l) +
     +                  wt(nn)*conjg(ylm(lm0+mp))*aimag(ylms(lm))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!===>    orthonormalize the projections for each l (many maybe zero)
!
      nlh = 0
      DO l = 0, lmax
         lm0 = l*(l+1)+1
         l2  = 2*l
         DO m = 0, l2
!--->       calculate overlaps with previous harmonics for this l
            DO j = 0, m-1
               s = real( conjg(a(0,j,l)) * a(0,m,l))
               DO mp=1,l
                  s = s + 2*real( conjg(a(mp,j,l)) * a(mp,m,l) )
               ENDDO
               ovlp(j) = s
            ENDDO
!--->       Gram-Schmidt
            DO j = 0,m-1
               a(0:l,m,l) = a(0:l,m,l) - ovlp(j)*a(0:l,j,l)
            ENDDO
!--->       normalize
            s = real( conjg(a(0,m,l)) * a(0,m,l))
            DO mp = 1, l
               s = s + 2*real( conjg(a(mp,m,l)) * a(mp,m,l))
            ENDDO
            IF (s.GT.del) THEN
               s = sqrt(s)
               a(0:l,m,l) = a(0:l,m,l)/s
!--->          store lattice harmonic
               mems = 0
               nlh = nlh + 1
               IF (nlh>nlhd+1)  CALL juDFT_error("nlhd too small",
     +              calledby="lhcal")
               lnu(nlh) = l
               IF (abs(a(0,m,l)).GT.del) THEN
                  mems = mems+1
                  IF (mems>memd)  CALL juDFT_error("memd too small"
     +                 ,calledby ="lhcal")
                  c(mems,nlh) = a(0,m,l)
                  lmnu(mems,nlh) = lm0
               ENDIF
               DO mp=1,l
                  IF( abs(a(mp,m,l)).GT.del) THEN
                     mems = mems + 1
                     IF(mems>memd)  CALL juDFT_error("memd too small"
     +                    ,calledby ="lhcal")
                     c(mems,nlh) = a(mp,m,l)
                     lmnu(mems,nlh) = lm0 + mp
                     mems = mems + 1
                     IF (mems>memd)  CALL juDFT_error("memd too small"
     +                    ,calledby ="lhcal")
                     c(mems,nlh) = ((-1.)**mp)*conjg(c(mems-1,nlh))
                     lmnu(mems,nlh) = lm0 - mp
                  ENDIF
               ENDDO
               mem(nlh) = mems
            ELSE
               a(0:l,m,l) = cmplx(0.0,0.0)
            ENDIF
         ENDDO   ! m = 0, l2
      ENDDO      ! l = 0, lmax
!
!===>    test of lattice harmonics using an arbitary point
!
      v(1) = sqrt(2.0)
      v(2) = sqrt(5.0)
      v(3) = sqrt(17.0)
!---> generate lattice harmonic for this point
      CALL ylm4(
     >          lmax,v,
     <          ylm)
      DO lh = 1, nlh
         ylms(lh) = cmplx(0.0,0.0)
         DO m = 1, mem(lh)
            ylms(lh) = ylms(lh) + c(m,lh)*ylm(lmnu(m,lh))
         ENDDO
      ENDDO
!---> rotate point and generate lattice harmonic
      DO n = 2, nrot
         vr = matmul( orth(:,:,n), v )
         CALL ylm4(
     >             lmax,vr,
     <             ylm)
         DO lh = 1, nlh
            ylmr(lh) = cmplx(0.0,0.0)
            DO m = 1, mem(lh)
               ylmr(lh) = ylmr(lh) + c(m,lh)*ylm(lmnu(m,lh))
            ENDDO
            IF ( abs(ylms(lh)-ylmr(lh)).GT.del ) THEN
               WRITE (6,'(/," error for operation",i3)') n
               WRITE (6,'(  " lattice harmonic   ",i3)') lh
               WRITE (6,'(4f12.6)') ylms(lh),ylmr(lh)
                CALL juDFT_error("k_lv(Rr)",calledby="lhcal")
            ENDIF
         ENDDO
      ENDDO
!$OMP END MASTER
      RETURN
      END SUBROUTINE lhcal
      END MODULE m_lhcal
