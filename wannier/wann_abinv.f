!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_abinv
      contains
      SUBROUTINE wann_abinv(
     >               ntypd,natd,neigd,lmaxd,lmd,llod,nlod,ntype,neq,
     >               neig,lmax,nlo,llo,invsat,invsatnr,bkpt,taual,
     X               acof,bcof,ccof)
C     ***************************************************************
C     Transform acof,bcof,ccof in case of atoms related by inversion
c     symmetry to obtain the coefficients in the global frame.
c     Based on abcrot.
c     Frank Freimuth
C     ***************************************************************
      use m_constants, only:pimach
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,natd,neigd,lmd,llod,nlod,ntype
      INTEGER, INTENT (IN) :: lmaxd,neig
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: neq(ntypd),lmax(ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd)
      INTEGER, INTENT (IN) :: invsat(natd),invsatnr(natd)
      real,intent(in)      :: bkpt(3)
      REAL,    INTENT (IN) :: taual(3,natd)

      COMPLEX, INTENT (INOUT) :: acof(neigd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: bcof(neigd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: ccof(-llod:llod,neigd,nlod,natd)

C     .. Local Scalars ..
      INTEGER :: itype,ineq,iatom,iop,ilo,i,l,m,lm,lmp,ifac
      integer :: n,nn,jatom,ie,ll1
      real    :: tpi,arg
      complex :: fac
C     ..

      tpi=2.0*pimach()

      iatom=0
      DO itype=1,ntype
        DO ineq=1,neq(itype)
          iatom=iatom+1
          IF(invsat(iatom).ne.2) cycle
          DO l=1,lmax(itype),2
            DO i=1,neig
              acof(i,l**2:l*(l+2),iatom) = (-1)**l *
     &                       acof(i,l**2:l*(l+2),iatom)
              bcof(i,l**2:l*(l+2),iatom) = (-1)**l * 
     &                       bcof(i,l**2:l*(l+2),iatom)
            ENDDO
          ENDDO
          DO ilo=1,nlo(itype)
            l=llo(ilo,itype)
            IF(l.gt.0) THEN
              if(mod(l,2).eq.0)cycle 
              DO i=1,neig
                ccof(-l:l,i,ilo,iatom) = (-1)**l * 
     &                       ccof(-l:l,i,ilo,iatom)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

c$$$      iatom = 0
c$$$      DO n = 1,ntype
c$$$         DO nn = 1,neq(n)
c$$$            iatom = iatom + 1
c$$$            IF (invsat(iatom).EQ.1) THEN
c$$$               jatom = invsatnr(iatom)
c$$$               arg=    (taual(1,jatom)+taual(1,iatom))*bkpt(1)
c$$$               arg=arg+(taual(2,jatom)+taual(2,iatom))*bkpt(2)
c$$$               arg=arg+(taual(3,jatom)+taual(3,iatom))*bkpt(3)
c$$$               arg=arg*tpi
c$$$               fac=cmplx(cos(arg),sin(arg))
c$$$               DO ilo = 1,nlo(n)
c$$$                  l = llo(ilo,n)
c$$$                  DO m = -l,l
c$$$                     DO ie = 1,neig
c$$$                        ccof(m,ie,ilo,jatom) = fac *
c$$$     +                              ccof(m,ie,ilo,jatom)
c$$$                     ENDDO
c$$$                  ENDDO
c$$$               ENDDO
c$$$               DO l = 0,lmax(n)
c$$$                  ll1 = l* (l+1)
c$$$                  DO m =-l,l
c$$$                     lm  = ll1 + m
c$$$                     DO ie = 1,neig
c$$$                        acof(ie,lm,jatom) = fac *
c$$$     *                              acof(ie,lm,jatom)
c$$$                     ENDDO
c$$$                     DO ie = 1,neig
c$$$                        bcof(ie,lm,jatom) = fac *
c$$$     *                              bcof(ie,lm,jatom)
c$$$                     ENDDO
c$$$                  ENDDO
c$$$               ENDDO
c$$$            ENDIF
c$$$         ENDDO
c$$$      ENDDO

      END subroutine
      end module
      

    
