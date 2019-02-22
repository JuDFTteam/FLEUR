MODULE m_relaxation
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  integer:: input_force_relax=3
  public relaxation
  
CONTAINS
  
  SUBROUTINE relaxation(mpi,input,atoms,cell,sym,force_new,energies_new)
    USE m_types
    use m_relaxio
    use m_broyd_io
    TYPE(t_mpi),INTENT(IN)   :: mpi
    TYPE(t_input),INTENT(IN) :: input
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_sym),INTENT(IN)   :: sym
    TYPE(t_cell),INTENT(IN)  :: cell
    REAL,INTENT(in)  :: force_new(:,:),energies_new
    
    REAL,ALLOCATABLE :: pos(:,:,:),force(:,:,:),energies(:)
    REAL,ALLOCATABLE :: displace(:,:),old_displace(:,:)
    REAL             :: alpha
    INTEGER          :: n
  
    IF (mpi%irank==0) THEN
       ALLOCATE(pos(3,atoms%ntype,1)); 
       DO n=1,atoms%ntype
          pos(:,n,1)=atoms%taual(:,SUM(atoms%neq(:n-1))+1)
       END DO
       ALLOCATE(force(3,atoms%ntype,1)); force(:,:,1)=force_new
       ALLOCATE(energies(1));energies(1)=energies_new
       ALLOCATE(displace(3,atoms%ntype),old_displace(3,atoms%ntype))
    
    ! add history 
    CALL read_relax(pos,force,energies)
    
    !determine new positions
    IF (SIZE(energies)==1.OR.input_force_relax==0) THEN
       !no history present simple step
       ! choose a reasonable first guess for scaling
       ! this choice is based on a Debye temperature of 330K;
       ! modify as needed
       alpha = (250.0/(MAXVAL(atoms%zatom)*input%xa))*((330./input%thetad)**2)
       CALL simple_step(alpha,force,displace)
    ELSEIF (input_force_relax==1) THEN
       CALL simple_cg(pos,force,displace)
    ELSE
       CALL simple_bfgs(pos,force,displace)
    ENDIF
    
    CALL read_displacements(atoms,old_displace)
    DO n=1,atoms%ntype
       PRINT *,"OD:",old_displace(:,n)
       PRINT *,"ND:",displace(:,n)
    END DO

    displace=displace+old_displace
    
    !Write file
    CALL write_relax(pos,force,energies,displace)
 ENDIF
 CALL resetBroydenHistory()
 CALL judft_end("Structual relaxation done",0)
  END SUBROUTINE relaxation
  


  SUBROUTINE simple_step(alpha,force,displace)
    !-----------------------------------------------
    IMPLICIT NONE
    REAL,INTENT(in)  :: alpha
    REAL,INTENT(in)  :: force(:,:,:)
    REAL,INTENT(OUT) :: displace(:,:)
    
    
    displace = alpha*force(:,:,SIZE(force,3))
  END SUBROUTINE simple_step
  
  SUBROUTINE simple_bfgs(pos,force,shift)
    !-----------------------------------------------
    !  Simple BFGS method to calculate shift out of old positions and forces
    !-----------------------------------------------
    IMPLICIT NONE
    REAL,INTENT(in)  :: pos(:,:,:),force(:,:,:)
    real,INTENT(OUT) :: shift(:,:)

    INTEGER         :: n,i,j,hist_length,n_force
    REAL,ALLOCATABLE:: h(:,:)
    REAL,ALLOCATABLE:: p(:),y(:),v(:)
    REAL            :: py,yy,gamma

    n_force=3*size(pos,2)
    allocate(h(n_force,n_force))
    allocate(p(n_force),y(n_force),v(n_force))

    !calculate approx. Hessian
    !initialize H
    h = 0.0
    DO n = 1,n_force
       h(n,n) = 1.0
    ENDDO
    !loop over all iterations (including current)
    hist_length=size(pos,3)
    DO n = 2,hist_length
       ! differences
       p(:) = RESHAPE(pos(:,:,n)-pos(:,:,n-1),(/SIZE(p)/))
       y(:) = RESHAPE(force(:,:,n)-force(:,:,n-1),(/SIZE(p)/))
       ! get necessary inner products and H|y>
       py = dot_PRODUCT(p,y)
       v = MATMUL(y,h)
       yy = dot_PRODUCT(y,v)
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
    shift = reshape(MATMUL(y,h),shape(shift))
  END SUBROUTINE simple_bfgs
 
  SUBROUTINE simple_cg(pos,force,shift)
    !-----------------------------------------------
    IMPLICIT NONE
    REAL,intent(in)  :: pos(:,:,:),force(:,:,:)
    real,INTENT(OUT) :: shift(:,:)

    REAL                :: f1(3,SIZE(pos,2)),f2(3,SIZE(pos,2))
    INTEGER             :: n_old
    
    n_old = SIZE(pos,3)-1
    
    f1 = (force(:,:,n_old+1)-force(:,:,n_old))/(pos(:,:,n_old+1)-pos(:,:,n_old))
    f2 = force(:,:,n_old+1)-f1*pos(:,:,n_old+1)
    shift = -1.*f2/f1-force(:,:,n_old+1)
  END SUBROUTINE simple_cg
END MODULE m_relaxation
