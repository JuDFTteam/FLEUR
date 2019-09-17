!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_forcew
  ! ************************************************************
  ! Printing force components
  ! ************************************************************
CONTAINS
  SUBROUTINE force_w(mpi,input,atoms,sym,results,cell,oneD,vacuum)
    USE m_types
    USE m_xmlOutput
    USE m_relaxation
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    !     ..
    !     .. Local Scalars ..
    REAL sum
    INTEGER i,jsp,n,nat1,ierr
    REAL eps_force
    LOGICAL :: l_new,l_relax
    !     ..
    !     .. Local Arrays ..
    REAL forcetot(3,atoms%ntype)
    CHARACTER(LEN=20) :: attributes(7)
#ifdef CPP_MPI
    include 'mpif.h'
#endif
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
8000         FORMAT ('SPIN-',i1,1x,'FORCE FOR ATOM TYPE=',i3,2x,'X=',f7.3,&
                  &          3x,'Y=',f7.3,3x,'Z=',f7.3,5x,' FX_SP =',f9.6,' FY_SP =',&
                  &          f9.6,' FZ_SP =',f9.6)
          ENDIF
          nat1 = nat1 + atoms%neq(n)
       END DO
       !
       !     write total forces
       !
       WRITE  (6,8005)
8005   FORMAT (/,' ***** TOTAL FORCES ON ATOMS ***** ',/)
       IF (input%l_f) CALL openXMLElement('totalForcesOnRepresentativeAtoms',(/'units'/),(/'Htr/bohr'/))
       nat1 = 1
       forcetot = 0.0
       DO n = 1,atoms%ntype
          IF (atoms%l_geo(n)) THEN
             DO jsp = 1,input%jspins
                DO i = 1,3
                   forcetot(i,n) = forcetot(i,n) + results%force(i,n,jsp)
                END DO
             END DO

             WRITE (6,FMT=8010) n, (atoms%pos(i,nat1),i=1,3),&
                  &        (forcetot(i,n),i=1,3)
8010         FORMAT (' TOTAL FORCE FOR ATOM TYPE=',i3,2x,'X=',f7.3,3x,'Y=',&
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
                     RESHAPE((/8,1,1,1,3,3,3,6,12,12,12,12,12,12/),(/7,2/)))
             END IF
          END IF
          nat1 = nat1 + atoms%neq(n)
       END DO
       IF (input%l_f) CALL closeXMLElement('totalForcesOnRepresentativeAtoms')


       !Check convergence of force by comparing force with old_force
       sum=MAXVAL(ABS(forcetot - results%force_old))
       results%force_old(:,:)=forcetot !Store for next iteration
       results%force=0.0
       l_relax=sum<input%force_converged
       IF (.NOT.l_relax) THEN
          WRITE (6,8020) input%force_converged,sum
8020      FORMAT ('No new postions, force convergence required=',f8.5,'; max force distance=',f8.5)
       END IF
    ENDIF
#ifdef CPP_MPI
    CALL MPI_BCAST(l_relax,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif
    IF (l_relax.AND.input%l_f) CALL relaxation(mpi,input,atoms,cell,sym,forcetot,results%tote)

  END SUBROUTINE force_w
END MODULE m_forcew
