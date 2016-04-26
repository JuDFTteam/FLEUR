      MODULE m_forcew
! ************************************************************
! Printing force components
! ************************************************************
      CONTAINS
      SUBROUTINE force_w(&
     &                   input,atoms,sym,results,cell,oneD)
      USE m_geo
      USE m_relax
      USE m_types
      IMPLICIT NONE

      TYPE(t_results),INTENT(IN)   :: results
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_atoms),INTENT(IN)     :: atoms
!     ..
!     .. Local Scalars ..
      REAL,PARAMETER:: zero=0.0
      REAL sum
      INTEGER i,jsp,n,nat1,ierr
      REAL eps_force
      LOGICAL :: l_new
!     ..
!     .. Local Arrays ..
      REAL forcetot(3,atoms%ntypd)
!
!     write spin-dependent forces
!
      nat1 = 1
      DO n = 1,atoms%ntype
         IF (atoms%l_geo(n)) THEN
         IF (input%jspins.EQ.2) THEN
            DO jsp = 1,input%jspins
               WRITE (6,FMT=8000) jsp,n, (atoms%pos(i,nat1),i=1,3),&
     &           (results%force(i,n,jsp),i=1,3)
               WRITE (16,FMT=8000) jsp,n, (atoms%pos(i,nat1),i=1,3),&
     &           (results%force(i,n,jsp),i=1,3)
            END DO
         END IF
 8000    FORMAT ('SPIN-',i1,1x,'FORCE FOR ATOM TYPE=',i3,2x,'X=',f7.3,&
     &          3x,'Y=',f7.3,3x,'Z=',f7.3,5x,' FX_SP =',f9.6,' FY_SP =',&
     &          f9.6,' FZ_SP =',f9.6)
         ENDIF
         nat1 = nat1 + atoms%neq(n)
      END DO
!
!     write total forces
!
      WRITE  (6,8005)
      WRITE (16,8005)
 8005 FORMAT (/,' ***** TOTAL FORCES ON ATOMS ***** ',/)
      nat1 = 1
      DO n = 1,atoms%ntype
         IF (atoms%l_geo(n)) THEN
!
         DO i = 1,3
            forcetot(i,n) = zero
         END DO
         DO jsp = 1,input%jspins
            DO i = 1,3
               forcetot(i,n) = forcetot(i,n) + results%force(i,n,jsp)
            END DO
         END DO
!
         WRITE (6,FMT=8010) n, (atoms%pos(i,nat1),i=1,3),&
     &     (forcetot(i,n),i=1,3)
         WRITE (16,FMT=8010) n, (atoms%pos(i,nat1),i=1,3),&
     &     (forcetot(i,n),i=1,3)
 8010    FORMAT (' TOTAL FORCE FOR ATOM TYPE=',i3,2x,'X=',f7.3,3x,'Y=',&
     &          f7.3,3x,'Z=',f7.3,/,22x,' FX_TOT=',f9.6,&
     &          ' FY_TOT=',f9.6,' FZ_TOT=',f9.6)
!
         ENDIF
         nat1 = nat1 + atoms%neq(n)
      END DO

      sum=0.0
      DO n = 1,atoms%ntype
        IF (atoms%l_geo(n)) THEN
          DO i = 1,3
            sum = max(sum,(forcetot(i,n) - results%force_old(i,n))**2)
          ENDDO
        ENDIF
      ENDDO 
      sum=sqrt(sum)
!-roa
      eps_force=0.00001
      open(88,file='eps_force',form='formatted',status='old',err=188)
      read(88,'(f20.8)') eps_force
      close (88)
  188 continue
!-roa
 
      WRITE (6,8020) eps_force,sum
 8020 FORMAT ('eps_force=',f8.5,'max=',f8.5)

      INQUIRE(file ="relax_inp",exist= l_new)
      IF (l_new) THEN
        CALL relax(input%film,atoms%pos,atoms%neq,sym%mrot,sym%tau,cell%amat,cell%bmat,atoms%ngopr,sym%invtab&
     &        ,forcetot)
      ELSE

         IF ((sum<eps_force).AND.input%l_f) THEN
            CALL geo(&
     &           atoms,sym,cell,&
     &           oneD,input,results%tote,forcetot)
         END IF
      ENDIF

      END SUBROUTINE force_w
      END MODULE m_forcew
