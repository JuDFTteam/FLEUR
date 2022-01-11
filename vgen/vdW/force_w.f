      MODULE m_forcew
c ************************************************************
c Printing force components
c ************************************************************
      CONTAINS
      SUBROUTINE force_w(
     >                   film,ntypd,jspd,natd,nwdd,nop,layerd,nlod,
     >                   ntype,jspins,neq,pos,tau,mrot,ngopr,l_f,
     >                   invtab,l_geo,amat,bmat,zatom,odi,ods,odd,
     X                   force,force_old,tote)
c
      USE m_od_types, ONLY : od_inp, od_sym, od_dim
      USE m_geo
      USE m_relax
      USE m_vdWfleur
      IMPLICIT NONE
C ..
C ..  Scalar Arguments ..
      LOGICAL,INTENT(IN)   :: film
      INTEGER, INTENT (IN) :: ntypd,jspd,natd,nwdd,nop,layerd,nlod
      INTEGER, INTENT (IN) :: ntype,jspins
      REAL,    INTENT (INOUT) :: tote  ! changed in case of vdW
      LOGICAL, INTENT (IN) :: l_f
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: ngopr(natd),neq(ntypd)
      INTEGER, INTENT (IN) :: mrot(3,3,nop),invtab(nop)
      REAL,    INTENT (IN) :: tau(3,nop),pos(3,natd)
      LOGICAL, INTENT (IN) :: l_geo(ntypd)
      REAL,    INTENT (INOUT) :: force(3,ntypd,jspd) ! changed in case of vdW
      REAL,    INTENT (INOUT) ::  force_old(3,ntypd)
      REAL   ,INTENT(IN)      :: amat(3,3),bmat(3,3),zatom(ntype)
c-odim
      TYPE (od_inp), INTENT (IN) :: odi
      TYPE (od_sym), INTENT (IN) :: ods
      TYPE (od_dim), INTENT (IN) :: odd
c+odim
C     ..
C     .. Local Scalars ..
      REAL zero,sum,e_vdW
      INTEGER i,jsp,n,nat1
      REAL eps_force
      LOGICAL :: l_new,l_vdW
C     ..
C     .. Local Arrays ..
      REAL forcetot(3,ntype),f_vdW(3,ntype)
C     ..
C     .. Data statements ..
      DATA zero/0.000/
      l_vdW = .false.
      l_vdW = .true.
C     ..
c
c     write spin-dependent forces
c
      nat1 = 1
      DO n = 1,ntype
         IF (l_geo(n)) THEN
         IF (jspins.EQ.2) THEN
            DO jsp = 1,jspins
               WRITE (6,FMT=8000) jsp,n, (pos(i,nat1),i=1,3),
     +           (force(i,n,jsp),i=1,3)
               WRITE (16,FMT=8000) jsp,n, (pos(i,nat1),i=1,3),
     +           (force(i,n,jsp),i=1,3)
            END DO
         END IF
 8000    FORMAT ('SPIN-',i1,1x,'FORCE FOR ATOM TYPE=',i3,2x,'X=',f7.3,
     +          3x,'Y=',f7.3,3x,'Z=',f7.3,5x,' FX_SP =',f9.6,' FY_SP =',
     +          f9.6,' FZ_SP =',f9.6)
         ENDIF
         nat1 = nat1 + neq(n)
      END DO
c
c     write total forces
c
      WRITE  (6,8005)
      WRITE (16,8005)
 8005 FORMAT (/,' ***** TOTAL FORCES ON ATOMS ***** ',/)
      nat1 = 1
      DO n = 1,ntype
         IF (l_geo(n)) THEN
c
         DO i = 1,3
            forcetot(i,n) = zero
         END DO
         DO jsp = 1,jspins
            DO i = 1,3
               forcetot(i,n) = forcetot(i,n) + force(i,n,jsp)
            END DO
         END DO
c
         WRITE (6,FMT=8010) n, (pos(i,nat1),i=1,3),
     +     (forcetot(i,n),i=1,3)
         WRITE (16,FMT=8010) n, (pos(i,nat1),i=1,3),
     +     (forcetot(i,n),i=1,3)
 8010    FORMAT (' TOTAL FORCE FOR ATOM TYPE=',i3,2x,'X=',f7.3,3x,'Y=',
     +          f7.3,3x,'Z=',f7.3,/,22x,' FX_TOT=',f9.6,
     +          ' FY_TOT=',f9.6,' FZ_TOT=',f9.6)
c
         ENDIF
         nat1 = nat1 + neq(n)
      END DO

      sum=0.0
      DO n = 1,ntype
        IF (l_geo(n)) THEN
          DO i = 1,3
            sum = max(sum,(forcetot(i,n) - force_old(i,n))**2)
          ENDDO
        ENDIF
      ENDDO 
      sum=sqrt(sum)
c-roa
      eps_force=0.00001
      open(88,file='eps_force',form='formatted',status='old',err=188)
      read(88,'(f20.8)') eps_force
      close (88)
  188 continue
c-roa
 
      WRITE (6,8020) eps_force,sum
 8020 FORMAT ('eps_force=',f8.5,'max=',f8.5)

      INQUIRE(file ="relax_inp",exist= l_new)
      IF (l_new) THEN
        CALL relax(film,pos,neq,mrot,tau,amat,bmat,ngopr,invtab
     $        ,forcetot)
      ELSE
        IF ((sum<eps_force).AND.l_f) THEN
!
! vdW
!
          IF (l_vdW) THEN

            CALL vdW_fleur(
     >                     ntype,natd,neq,film,amat,bmat,pos,zatom,
     >                     nop,ngopr,invtab,mrot,
     <                     e_vdW,f_vdW)

            DO n = 1,ntype
              IF (l_geo(n)) THEN
                forcetot(:,n) = forcetot(:,n) + f_vdW(:,n)
              ENDIF
            ENDDO

            tote = tote + e_vdW 
          ENDIF

          CALL geo(
     >             ntypd,natd,nwdd,nop,layerd,nlod,
     >             tau,mrot,ngopr,invtab,odi,ods,odd,tote,
     X             forcetot)
        ENDIF
      ENDIF

      END SUBROUTINE force_w
      END MODULE m_forcew
