!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_abcrot

      CONTAINS

      SUBROUTINE abcrot(
     >                 ntypd,natd,neigd,lmaxd,lmd,llod,nlod,ntype,neq,
     >                 neig,lmax,nlo,llo,nop,ngopr,mrot,invsat,invsatnr,
     >                 bmat,odi,ods,
     X                 acof,bcof,ccof)
C     ***************************************************************
C     * This routine transforms a/b/cof which are given wrt rotated *
C     * MT functions (according to invsat/ngopr) into a/b/cof wrt   *
C     * unrotated MT functions. Needed for GW calculations.         *
C     *                                                             *
C     * Christoph Friedrich Mar/2005                                *
C     ***************************************************************
      USE m_dwigner
      use m_savewigner
      USE m_types
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,natd,neigd,lmd,llod,nlod,ntype,nop
      INTEGER, INTENT (IN) :: lmaxd,neig
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),ngopr(natd),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: invsat(natd),invsatnr(natd)

      REAL,    INTENT (IN) :: bmat(3,3)
      COMPLEX, INTENT (INOUT) :: acof(neigd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: bcof(neigd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: ccof(-llod:llod,neigd,nlod,natd)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
c+odim
C     ..
C     .. Local Scalars ..
      INTEGER itype,ineq,iatom,iop,ilo,i,l,m,lm,lmp,ifac
C     ..
C     .. Local Arrays ..
c***** COMPLEX, ALLOCATABLE :: d_wgn(:,:,:,:) !put into module m_savewigner
C

      IF ( .NOT.ALLOCATED(d_wgn) ) THEN    !calculate d_wgn only once
        PRINT*,"calculate wigner-matrix"
        IF (.NOT.odi%d1) THEN
          ALLOCATE (d_wgn(-lmaxd:lmaxd,-lmaxd:lmaxd,lmaxd,nop))
          d_wgn = CMPLX(0.0,0.0) ! Initialization is done to avoid complaints by Lord Valgrind
          CALL d_wigner(nop,mrot,bmat,lmaxd,d_wgn)
        ELSE
          ALLOCATE (d_wgn(-lmaxd:lmaxd,-lmaxd:lmaxd,lmaxd,ods%nop))
          d_wgn = CMPLX(0.0,0.0) ! Initialization is done to avoid complaints by Lord Valgrind
          CALL d_wigner(ods%nop,ods%mrot,bmat,lmaxd,d_wgn)
        ENDIF
      ENDIF

      iatom=0
      DO itype=1,ntype
        DO ineq=1,neq(itype)
          iatom=iatom+1
          IF (.NOT.odi%d1) THEN
             iop=ngopr(iatom)
          ELSE
             iop=ods%ngopr(iatom)
          ENDIF
C                                    l                        l    l
C inversion of spherical harmonics: Y (pi-theta,pi+phi) = (-1)  * Y (theta,phi)
C                                    m                             m
          ifac = 1
          IF(invsat(iatom).EQ.2) THEN
            IF (.NOT.odi%d1) THEN
               iop=ngopr(invsatnr(iatom))
            ELSE
               iop=ods%ngopr(invsatnr(iatom))
            ENDIF
            ifac = -1 
          ENDIF
          DO l=1,lmax(itype)
c  replaced d_wgn by conjg(d_wgn),FF October 2006
            DO i=1,neig
              acof(i,l**2:l*(l+2),iatom) = ifac**l * matmul(
     &                                 conjg(d_wgn(-l:l,-l:l,l,iop)),
     &                                 acof(i,l**2:l*(l+2),iatom))
              bcof(i,l**2:l*(l+2),iatom) = ifac**l * matmul(
     &                                 conjg(d_wgn(-l:l,-l:l,l,iop)),
     &                                 bcof(i,l**2:l*(l+2),iatom))
            ENDDO
          ENDDO
          DO ilo=1,nlo(itype)
            l=llo(ilo,itype)
            IF(l.gt.0) THEN
              DO i=1,neig
                ccof(-l:l,i,ilo,iatom) = ifac**l * matmul(
     &                               conjg(d_wgn(-l:l,-l:l,l,iop)),
     &                               ccof(-l:l,i,ilo,iatom))
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      END SUBROUTINE abcrot

      END MODULE m_abcrot
