!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_forcew
   !-----------------------------------------------------------------------------
   ! Printing force components
   !-----------------------------------------------------------------------------

#ifdef CPP_MPI
   USE mpi 
#endif

CONTAINS
   SUBROUTINE force_w(fmpi,input,atoms,sym,results,cell,oneD,vacuum)
      USE m_types
      USE m_constants
      USE m_xmlOutput
      USE m_relaxation
      USE m_rotate_forces
      USE m_vdWfleur_grimme
      IMPLICIT NONE

      TYPE(t_mpi),     INTENT(IN)    :: fmpi
      TYPE(t_results), INTENT(INOUT) :: results
      TYPE(t_oneD),    INTENT(IN)    :: oneD
      TYPE(t_input),   INTENT(IN)    :: input
      TYPE(t_sym),     INTENT(IN)    :: sym
      TYPE(t_cell),    INTENT(IN)    :: cell
      TYPE(t_atoms),   INTENT(IN)    :: atoms
      TYPE(t_vacuum),  INTENT(IN)    :: vacuum

      REAL maxAbsForceDist
      INTEGER i, jsp, n, nat1, ierr
      REAL eps_force,e_vdW
      LOGICAL :: l_new, l_forceConverged
      REAL,ALLOCATABLE:: f_vdW(:,:)
      REAL forcetot(3,atoms%ntype)
      CHARACTER(LEN=20) :: attributes(7)

      ! Write spin-dependent forces

      IF (fmpi%irank==0) THEN
         nat1 = 1
         DO n = 1, atoms%ntype
            IF (atoms%l_geo(n)) THEN
               IF (input%jspins.EQ.2) THEN
                  DO jsp = 1,input%jspins
                     WRITE (oUnit,FMT=8000) jsp, n, (atoms%pos(i,nat1),i=1,3), &
                                                    (results%force(i,n,jsp),i=1,3)
                  END DO
               END IF

8000           FORMAT ('SPIN-',i1,1x,'FORCE FOR ATOM TYPE=',i3,2x,&
                       'X=',f7.3,3x,'Y=',f7.3,3x,'Z=',f7.3,5x,&
                       ' FX_SP =',f9.6,' FY_SP =',f9.6,' FZ_SP =',f9.6)
            END IF
            nat1 = nat1 + atoms%neq(n)
         END DO

         ! Write total forces
         WRITE  (oUnit,8005)

8005     FORMAT (/,' ***** TOTAL FORCES ON ATOMS ***** ',/)
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

               WRITE (oUnit,FMT=8010) n, (atoms%pos(i,nat1),i=1,3), &
                                         (forcetot(i,n),i=1,3)

8010           FORMAT (' TOTAL FORCE FOR ATOM TYPE=',i3,2x,&
                       'X=',f7.3,3x,'Y=',f7.3,3x,'Z=',f7.3,/,22x,&
                       ' FX_TOT=',f9.6,' FY_TOT=',f9.6,' FZ_TOT=',f9.6)

               WRITE(attributes(1),'(i0)') n
               WRITE(attributes(2),'(f12.6)') atoms%pos(1,nat1)
               WRITE(attributes(3),'(f12.6)') atoms%pos(2,nat1)
               WRITE(attributes(4),'(f12.6)') atoms%pos(3,nat1)
               WRITE(attributes(5),'(f12.8)') forcetot(1,n)
               WRITE(attributes(6),'(f12.8)') forcetot(2,n)
               WRITE(attributes(7),'(f12.8)') forcetot(3,n)

               IF (input%l_f) THEN
                  CALL writeXMLElementFormPoly('forceTotal',(/'atomType',&
                       'x       ','y       ','z       ',&
                       'F_x     ','F_y     ','F_z     '/),attributes,&
                        RESHAPE((/8,1,1,1,3,3,3,6,12,12,12,12,12,12/),(/7,2/)))
               END IF
            END IF
            nat1 = nat1 + atoms%neq(n)
         END DO
         IF (input%l_f) CALL closeXMLElement('totalForcesOnRepresentativeAtoms')

         ! Check convergence of force by comparing force with old_force
         maxAbsForceDist=MAXVAL(ABS(forcetot - results%force_old))
         results%force_old(:,:)=forcetot !Store for next iteration
         results%force=0.0

         l_forceConverged=maxAbsForceDist<input%force_converged
         l_forceConverged=l_forceConverged.and.(results%last_distance.LE.input%mindistance)
         l_forceConverged=l_forceConverged.and.(results%last_distance.GE.0.0)
         ! In the first iteration results%last_distance is initialized to -1.0.

         IF (.NOT.l_forceConverged) THEN
            WRITE (oUnit,8020) input%force_converged,maxAbsForceDist
8020        FORMAT ('No new postions, force convergence required=',f8.5,'; max force distance=',f8.5)
         END IF
      END IF

#ifdef CPP_MPI
      CALL MPI_BCAST(l_forceConverged,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
#endif
      IF (l_forceConverged.AND.input%l_f.AND.(input%f_level.GE.0).AND.(fmpi%irank==0)) THEN
         CALL rotate_forces(atoms%ntype,atoms%ntype,atoms%nat,sym%nop,results%tote,&
                            cell%omtil,atoms%neq,sym%mrot,cell%amat,cell%bmat,&
                            atoms%taual,sym%tau,forcetot)
      END IF

      IF (l_forceConverged.and.btest(input%vdW,0)) THEN
         ALLOCATE(f_vdW,mold=forcetot)
         call vdW_fleur_grimme(input,atoms,sym,cell,e_vdW,f_vdW)
         forcetot=forcetot+f_vdW
         results%tote=results%tote+e_vdW
      ENDIF

      IF (l_forceConverged.AND.input%l_f) CALL relaxation(fmpi,input,atoms,cell,sym,oneD,vacuum,forcetot,results%tote)

   END SUBROUTINE force_w

END MODULE m_forcew
