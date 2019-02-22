      MODULE m_forcew
! ************************************************************
! Printing force components
! ************************************************************
      CONTAINS
        SUBROUTINE force_w(mpi,input,atoms,sym,results,cell,oneD,vacuum)
      USE m_geo
      USE m_relax
      USE m_types
      USE m_xmlOutput
      use m_relaxation
      IMPLICIT NONE
      TYPE(t_mpi),INTENT(IN)       :: mpi
      TYPE(t_results),INTENT(IN)   :: results
      TYPE(t_oneD),INTENT(IN)      :: oneD
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_sym),INTENT(IN)       :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_vacuum),INTENT(IN)    :: vacuum
!     ..
!     .. Local Scalars ..
      REAL,PARAMETER:: zero=0.0
      REAL sum
      INTEGER i,jsp,n,nat1,ierr
      REAL eps_force
      LOGICAL :: l_new,l_relax
!     ..
!     .. Local Arrays ..
      REAL forcetot(3,atoms%ntype)
      CHARACTER(LEN=20) :: attributes(7)
!
!     write spin-dependent forces
!
      IF (mpi%irank==0) THEN
      nat1 = 1
      DO n = 1,atoms%ntype
         IF (atoms%l_geo(n)) THEN
         IF (input%jspins.EQ.2) THEN
            DO jsp = 1,input%jspins
               WRITE (6,FMT=8000) jsp,n, (atoms%pos(i,nat1),i=1,3),&
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
 8005 FORMAT (/,' ***** TOTAL FORCES ON ATOMS ***** ',/)
      IF (input%l_f) CALL openXMLElement('totalForcesOnRepresentativeAtoms',(/'units'/),(/'Htr/bohr'/))
      nat1 = 1
      DO n = 1,atoms%ntype
         IF (atoms%l_geo(n)) THEN
            DO i = 1,3
               forcetot(i,n) = zero
            END DO
            DO jsp = 1,input%jspins
               DO i = 1,3
                  forcetot(i,n) = forcetot(i,n) + results%force(i,n,jsp)
               END DO
            END DO

            WRITE (6,FMT=8010) n, (atoms%pos(i,nat1),i=1,3),&
     &        (forcetot(i,n),i=1,3)
  8010       FORMAT (' TOTAL FORCE FOR ATOM TYPE=',i3,2x,'X=',f7.3,3x,'Y=',&
     &              f7.3,3x,'Z=',f7.3,/,22x,' FX_TOT=',f9.6,&
     &              ' FY_TOT=',f9.6,' FZ_TOT=',f9.6)

            WRITE(attributes(1),'(i0)') n
            WRITE(attributes(2),'(f12.6)') atoms%pos(1,nat1)
            WRITE(attributes(3),'(f12.6)') atoms%pos(2,nat1)
            WRITE(attributes(4),'(f12.6)') atoms%pos(3,nat1)
            WRITE(attributes(5),'(f12.8)') forcetot(1,n)
            WRITE(attributes(6),'(f12.8)') forcetot(2,n)
            WRITE(attributes(7),'(f12.8)') forcetot(3,n)
            IF (input%l_f) THEN
               CALL writeXMLElementFormPoly('forceTotal',(/'atomType','x       ','y       ','z       ',&
                                                           'F_x     ','F_y     ','F_z     '/),attributes,&
                                            reshape((/8,1,1,1,3,3,3,6,12,12,12,12,12,12/),(/7,2/)))
            END IF
         END IF
         nat1 = nat1 + atoms%neq(n)
      END DO
      IF (input%l_f) CALL closeXMLElement('totalForcesOnRepresentativeAtoms')

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
   ENDIF
   l_relax=sum<eps_force
#ifdef CPP_MPI
   CALL MPI_BCAST(l_relax,1,MPI_LOGICAL,0,ierr)
#endif
   IF (l_relax.and.input%l_f) CALL relaxation(mpi,input,atoms,cell,sym,forcetot,results%tote)
      
 END SUBROUTINE force_w
END MODULE m_forcew
