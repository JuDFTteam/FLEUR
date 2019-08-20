!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_relaxation
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  PUBLIC relaxation !This is the interface. Below there are internal subroutines for bfgs, simple mixing, CG ...

CONTAINS
  SUBROUTINE relaxation(mpi,input,atoms,cell,sym,force_new,energies_new)
    !This routine uses the current force,energies and atomic positions to 
    !generate a displacement in a relaxation step. 
    !The history is taken into account by read_relax from m_relaxio
    !After generating new positions the code stops
    USE m_types
    USE m_relaxio
    USE m_mixing_history
#ifdef CPP_MPI
    INCLUDE 'mpif.h'
#endif
    TYPE(t_mpi),INTENT(IN)   :: mpi
    TYPE(t_input),INTENT(IN) :: input
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_cell),INTENT(IN)  :: cell
    REAL,INTENT(in)          :: force_new(:,:),energies_new !data for this iteration

    REAL,ALLOCATABLE :: pos(:,:,:),force(:,:,:),energies(:)
    REAL,ALLOCATABLE :: displace(:,:),old_displace(:,:)
    INTEGER          :: n,ierr
    LOGICAL          :: l_conv

    IF (mpi%irank==0) THEN
       ALLOCATE(pos(3,atoms%ntype,1)); 
       DO n=1,atoms%ntype
          pos(:,n,1)=atoms%pos(:,SUM(atoms%neq(:n-1))+1)
       END DO
       ALLOCATE(force(3,atoms%ntype,1)); force(:,:,1)=force_new
       ALLOCATE(energies(1));energies(1)=energies_new
       ALLOCATE(displace(3,atoms%ntype),old_displace(3,atoms%ntype))

       !Remove force components that are not selected for relaxation
       DO n=1,atoms%ntype
          IF (atoms%l_geo(n)) THEN
             force(:,n,1)=force(:,n,1)*REAL(atoms%relax(:,n))
          ELSE
             force(:,n,1)=0.0
          ENDIF
       ENDDO
       
       ! add history 
       CALL read_relax(pos,force,energies)

       !determine new positions
       IF (SIZE(energies)==1.OR.input%forcemix==0) THEN
          !no history present simple step
          ! choose a reasonable first guess for scaling
          ! this choice is based on a Debye temperature of 330K;
          ! modify as needed
          !alpha = (250.0/(MAXVAL(atoms%zatom)*input%xa))*((330./input%thetad)**2)
          CALL simple_step(input%forcealpha,0.25,force,displace)
       ELSE IF (input%forcemix==1) THEN
          CALL simple_cg(pos,force,displace)
       ELSE IF (input%forcemix==2) THEN
          CALL simple_bfgs(pos,force,displace)
       ELSE
          CALL juDFT_error('unkown mixing scheme for forces', calledby='relaxation')
       END IF

       !Check for convergence of forces/displacements
       l_conv=.TRUE.
       DO n=1,atoms%ntype
          IF (DOT_PRODUCT(force(:,n,SIZE(force,3)),force(:,n,SIZE(force,3)))>input%epsforce**2) l_conv=.FALSE.
          IF (DOT_PRODUCT(displace(:,n),displace(:,n))>input%epsdisp**2) l_conv=.FALSE.
       ENDDO

       !New displacements relative to positions in inp.xml
       CALL read_displacements(atoms,old_displace)
       displace=displace+old_displace

       !Write file
       CALL write_relax(pos,force,energies,displace)


    ENDIF
#ifdef CPP_MPI
    CALL MPI_BCAST(l_conv,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif
    IF (l_conv) THEN
       CALL judft_end("Structual relaxation: Done",mpi%irank)
    ELSE
       CALL mixing_history_reset(mpi)
       CALL judft_end("Structual relaxation: new displacements generated",mpi%irank)
    END IF
  END SUBROUTINE relaxation



  SUBROUTINE simple_step(alpha,maxdisp,force,displace)
    !-----------------------------------------------
    IMPLICIT NONE
    REAL,INTENT(in)  :: alpha,maxdisp
    REAL,INTENT(in)  :: force(:,:,:)
    REAL,INTENT(OUT) :: displace(:,:)

    real :: corr
    
    displace = alpha*force(:,:,SIZE(force,3))
    corr=maxdisp/maxval(abs(displace))
    if (corr<1.0) displace = corr*alpha*force(:,:,size(force,3))
    
  END SUBROUTINE simple_step

  SUBROUTINE simple_bfgs(pos,force,shift)
    !-----------------------------------------------
    !  Simple BFGS method to calculate shift out of old positions and forces
    !-----------------------------------------------
    IMPLICIT NONE
    REAL,INTENT(in)  :: pos(:,:,:),force(:,:,:)
    REAL,INTENT(OUT) :: shift(:,:)

    INTEGER         :: n,i,j,hist_length,n_force
    REAL,ALLOCATABLE:: h(:,:)
    REAL,ALLOCATABLE:: p(:),y(:),v(:)
    REAL            :: py,yy,gamma

    n_force=3*SIZE(pos,2)
    ALLOCATE(h(n_force,n_force))
    ALLOCATE(p(n_force),y(n_force),v(n_force))

    !calculate approx. Hessian
    !initialize H
    h = 0.0
    DO n = 1,n_force
       h(n,n) = 1.0
    ENDDO
    !loop over all iterations (including current)
    hist_length=SIZE(pos,3)
    DO n = 2,hist_length
       ! differences
       p(:) = RESHAPE(pos(:,:,n)-pos(:,:,n-1),(/SIZE(p)/))
       y(:) = RESHAPE(force(:,:,n)-force(:,:,n-1),(/SIZE(p)/))
       ! get necessary inner products and H|y>
       py = DOT_PRODUCT(p,y)
       v = MATMUL(y,h)
       yy = DOT_PRODUCT(y,v)
       !check that update will leave h positive definite;
       IF (py <= 0.0) THEN
          WRITE (6,*) '  bfgs: <p|y> < 0'
          WRITE (6,*) '  check convergence of forces'
          !Starting over with initial hessian
          h = 0.0
          DO j = 1,n_force
             h(j,j) = 1.0
          ENDDO
          CYCLE 
       ELSE
          !update h
          IF (n == 2) THEN
             gamma = py/yy
          ELSE
             gamma = 1.0
          ENDIF
          DO j = 1,n_force
             DO i = 1,n_force
                h(i,j) = (h(i,j) - (v(i)*p(j)+p(i)*v(j))/py)*gamma + (1.+gamma*yy/py)*p(i)*p(j)/py
             ENDDO
          ENDDO
       ENDIF
    ENDDO
    y(:) = RESHAPE(force(:,:,hist_length),(/SIZE(p)/))
    shift = RESHAPE(MATMUL(y,h),SHAPE(shift))
  END SUBROUTINE simple_bfgs

  SUBROUTINE simple_cg(pos,force,shift)
    !-----------------------------------------------
    IMPLICIT NONE
    REAL,INTENT(in)  :: pos(:,:,:),force(:,:,:)
    REAL,INTENT(OUT) :: shift(:,:)

    REAL                :: f1(3,SIZE(pos,2)),f2(3,SIZE(pos,2))
    INTEGER             :: n_old

    n_old = SIZE(pos,3)-1

    f1 = (force(:,:,n_old+1)-force(:,:,n_old))/(pos(:,:,n_old+1)-pos(:,:,n_old))
    f2 = force(:,:,n_old+1)-f1*pos(:,:,n_old+1)
    shift = -1.*f2/f1-force(:,:,n_old+1)
  END SUBROUTINE simple_cg
END MODULE m_relaxation
