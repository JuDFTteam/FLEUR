!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_writegw
      CONTAINS
      SUBROUTINE write_gw(
     >                    ntype,nop,nwd,jspins,natd,
     >                    ncst,neq,lmax,mrot,amat,bmat,rkmax,
     >                    taual,zatom,vol,scale,neigd,lmaxd,
     >                    nlod,llod,nlo,llo,l_soc)

      USE m_constants, ONLY : tpi_const
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: ntype,nop,nwd,jspins,natd,neigd,lmaxd
      INTEGER, INTENT (IN) :: neq(ntype),lmax(ntype),mrot(3,3,nop)
      INTEGER, INTENT (IN) :: ncst(ntype)
      INTEGER, INTENT (IN) :: nlod,llod,nlo(ntype),llo(nlod,ntype)
      REAL,    INTENT (IN) :: vol,scale
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3),rkmax
      REAL,    INTENT (IN) :: zatom(ntype),taual(3,natd)
      LOGICAL, INTENT (IN) :: l_soc

      INTEGER na,itype,i,iop,lmaxx,isp,j,lc,lsum,ll,ilo,nindx,lc_prev
      INTEGER , ALLOCATABLE :: l(:,:)
      REAL bbmat(3,3)

!
! Generate input file CLASS for subsequent GW calculation
! 19.9.2003 Arno Schindlmayr
!
      WRITE(6,'(A)') 'Info: Write out CLASS for GW.'
      OPEN(15,file='CLASS',status='unknown',action='write')
      na = 1
      DO itype = 1, ntype
        DO i = 1, neq(itype)
          WRITE(15,'(2(1x,i3))') na,itype
          na = na + 1
        ENDDO
      ENDDO
      CLOSE(15)
!
! Generate input file SYMOPS for subsequent GW calculation
! 19.9.2003 Arno Schindlmayr
!
      WRITE(6,'(A)') 'Info: Write out SYMOPS for GW.'
      bbmat = bmat/tpi_const
      OPEN(15,file='SYMOPS',status='unknown',action='write')
      WRITE(15,*) nop
      DO iop = 1, nop
        WRITE(15,*) iop
        DO i = 1, 3
          WRITE(15,*) MATMUL(amat(i,:),MATMUL(mrot(:,:,iop),bbmat))
        ENDDO
      ENDDO
      CLOSE(15)
!
! Generate input file LATTC for subsequent GW calculation
! 22.9.2003 Arno Schindlmayr
!
      WRITE(6,'(A)') 'Info: Write out LATTC for GW.'
      OPEN(15,file='LATTC',status='unknown',action='write')
      WRITE(15,'(e24.16)') scale
      WRITE(15,'(3e24.16)') amat(:,:)/scale
      WRITE(15,'(e24.16)') rkmax
      WRITE(15,*) ' ------------------------------------------- '
      lmaxx = MAXVAL(lmax(1:ntype))
      ALLOCATE (l(0:lmaxx,ntype))
      WRITE(15,'(2i4," ! nbas lmxamx (max l for augmentation)")')
     &  sum(neq(1:ntype)),lmaxx
      WRITE(15,*) ' ------------------------------------------- '
      DO isp = 1, jspins
        WRITE(15,'(" -- ibas lmxa konf(s) konf(p) konf(d)... ",
     &             " isp=",i2)') isp
        na = 1
        DO itype = 1, ntype
          DO i = 0, lmax(itype)
            l(i,itype) = i+1
          ENDDO
          SELECT CASE (ncst(itype))
            CASE (0)
              CONTINUE
            CASE (1)
              l(0,itype) = 2
            CASE (2)
              l(0,itype) = 3
            CASE (4)
              l(0:1,itype) = (/3,3/)
            CASE (5)
              l(0:1,itype) = (/4,3/)
            CASE (7)
              l(0:1,itype) = (/4,4/)
            CASE (9)
              l(0:2,itype) = (/4,4,4/)
            CASE (10)
              l(0:2,itype) = (/5,4,4/)
            CASE (12)
              l(0:2,itype) = (/5,5,4/)
            CASE (14)
              l(0:2,itype) = (/5,5,5/)
            CASE (15)
              l(0:2,itype) = (/6,5,5/)
            CASE (17)
              l(0:2,itype) = (/6,6,5/)
            CASE (19)
              l(0:3,itype) = (/6,6,5,5/)
            CASE (21)
              l(0:3,itype) = (/6,6,6,5/)
            CASE (22)
              l(0:3,itype) = (/7,6,6,5/)
            CASE (24)
              l(0:3,itype) = (/7,7,6,5/)
            CASE (26)
              l(0:3,itype) = (/7,7,6,6/)
            CASE DEFAULT
              l(:,itype) = 0
          END SELECT
          DO i = 1, neq(itype)
            IF (l(0,itype).GT.0) THEN
              WRITE(15,'(3x,99i4)') na,lmax(itype),
     &          l(0:lmax(itype),itype)
            ELSE
              WRITE(15,'(3x,2i4,3x,a)') na,lmax(itype),
     &          'WARNING: Unrecognized number of core levels!'
            ENDIF
            na = na + 1
          ENDDO
        ENDDO
      ENDDO
      CLOSE(15)
c      do i=1,ntype
c        write(*,*) ncst(i),sum(l(0:lmax(i),i))-(lmaxx+1)*(lmaxx+2)/2
c      enddo
!
! Generate input file NLAindx for subsequent GW calculation
! 29.9.2003 Arno Schindlmayr
!
      WRITE(6,'(A)') 'Info: Write out NLAindx for GW.'
      OPEN(15,file='NLAindx',status='unknown',action='write')
      WRITE(15,'(a)') '----NLAindx start---------------'
      lsum = 0
      DO itype = 1, ntype
        lsum = lsum + neq(itype)*2*(lmax(itype)+1)**2
        DO j = 1, nlo(itype)
          lsum = lsum + neq(itype)*(2*llo(j,itype)+1)
        ENDDO
      ENDDO
      WRITE(15,'(i6)') lsum
      lsum = 0
      DO j = 1, 2
        na = 0
        DO itype = 1, ntype
          DO i = 1, neq(itype)
            na = na + 1
            DO lc = 0, lmax(itype)
              ll=l(lc,itype)
              DO ilo = 1, nlo(itype)             ! Here, LOs are assumed to represent lowest valence states. 
                IF(llo(ilo,itype).eq.lc) ll=ll+1 ! states. This only concerns the basis function label
              ENDDO                              ! at the end of each line which is not read in.
              WRITE(15,'(i6,i3,i4,i6,3x,i2,a)') j,lc,na,lsum,
     &        ll,'SPDFGHIJKLMNO'(lc+1:lc+1)//'_'//'pd'(j:j)
              lsum = lsum + 2*lc+1
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! now for the local orbitals
      na = 0
      DO itype = 1, ntype
        DO i = 1, neq(itype)
          lc_prev = -1
          DO j = 1, nlo(itype)
            lc=llo(j,itype)
            IF(lc.eq.lc_prev) THEN
              nindx=nindx+1
            ELSE
              nindx=3
            ENDIF
            WRITE(15,'(i6,i3,i4,i6,3x,i2,a)') nindx,lc,na+i,lsum,
     &        l(lc,itype),'SPDFGHIJKLMNO'(lc+1:lc+1)//'_'//'l'
            lsum    = lsum + 2*lc+1
            lc_prev = lc
          ENDDO
        ENDDO
        na = na + neq(itype)
      ENDDO
      !
      CLOSE(15)
!
! Generate input file gwa for subsequent GW calculation
! 10.10.2003 Arno Schindlmayr
!
      WRITE(6,'(A)') 'Info: Write out gwa for GW.'
      OPEN(15,file='gwa',status='unknown',action='write',
     &     form='unformatted')
      WRITE(15) jspins,                                 ! nsp
     &          na,                                     ! nbas
     &          ntype,                                  ! nclass
     &          lmaxx,                                  ! lmxamx
     &          nlod
      WRITE(15) ((itype,i=1,neq(itype)),itype=1,ntype), ! iclass
     &          lmax(1:ntype),                          ! lmxa
     &          l(0:lmaxx,1:ntype),                     ! konf
     &          zatom(1:ntype),                         ! zz
     &          taual(:,1:na),                          ! bas
     &          scale,                                  ! alat
     &          amat,                                   ! plat
     &          vol,neigd,lmaxd,
     &          nlo(1:ntype),(llo(1:nlo(i),i),i=1,ntype)
      WRITE(15)
#ifdef CPP_INVERSION
     &  .true.,
#else
     &  .false.,
#endif      
     &  l_soc
      CLOSE(15) ! Subroutine eigen will append an additional record to gwa if gw=2.
      DEALLOCATE (l)


      END SUBROUTINE write_gw
      END MODULE m_writegw
