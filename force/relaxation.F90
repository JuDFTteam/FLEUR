!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_relaxation
   USE m_judft

#ifdef CPP_MPI
   USE mpi
#endif

   IMPLICIT NONE

   PRIVATE
   PUBLIC relaxation !This is the interface. Below there are internal subroutines for bfgs, simple mixing, CG ...

CONTAINS
   SUBROUTINE relaxation(fmpi,input,atoms,cell,sym ,vacuum,force_new,energies_new)
      ! This routine uses the current force,energies and atomic positions to
      ! generate a displacement in a relaxation step.
      ! The history is taken into account by read_relax from m_relaxio
      ! After generating new positions the code stops

      USE m_types
      USE m_constants
      USE m_relaxio
      USE m_mixing_history
      USE m_chkmt
      USE m_types_xml
      USE m_xsf_io

      TYPE(t_mpi),    INTENT(IN) :: fmpi
      TYPE(t_input),  INTENT(IN) :: input
      TYPE(t_atoms),  INTENT(IN) :: atoms
      TYPE(t_sym),    INTENT(IN) :: sym

      TYPE(t_vacuum), INTENT(IN) :: vacuum
      TYPE(t_cell),   INTENT(IN) :: cell
      REAL,           INTENT(IN) :: force_new(:,:), energies_new !data for this iteration

      REAL, ALLOCATABLE :: pos(:,:,:), force(:,:,:), energies(:)
      REAL, ALLOCATABLE :: displace(:,:), old_displace(:,:), tempDisplace(:,:)
      REAL, ALLOCATABLE :: totalDisplace(:,:)
      REAL              :: dispAll(3,atoms%nat), overlap(0:atoms%ntype,atoms%ntype)
      REAL              :: dispLength, maxDisp, limitDisp
      INTEGER           :: iType, ierr, numDispReduce
      LOGICAL           :: l_conv

      CHARACTER(len=100):: filename_add

      ! To calculate the current displacement
      TYPE(t_xml)   :: xml
      TYPE(t_atoms) :: atoms_non_displaced
      TYPE(t_atoms) :: tempAtoms

      IF (fmpi%irank==0) THEN
         filename_add = ""
         IF (judft_was_argument("-add_name")) filename_add = TRIM(judft_string_for_argument("-add_name"))//"_"
         CALL xml%init(filename_add)
         ALLOCATE(pos(3,atoms%ntype,1));

         DO iType = 1, atoms%ntype
            pos(:,iType,1)=atoms%pos(:,atoms%firstAtom(iType))
         END DO

         ALLOCATE(force(3,atoms%ntype,1)); force(:,:,1)=force_new
         ALLOCATE(energies(1));energies(1)=energies_new
         ALLOCATE(displace(3,atoms%ntype),old_displace(3,atoms%ntype))
         ALLOCATE(tempDisplace(3,atoms%ntype),totalDisplace(3,atoms%ntype))
         totalDisplace=0.0
         ! Remove force components that are not selected for relaxation
         DO iType = 1, atoms%ntype
            IF (atoms%l_geo(iType)) THEN
               force(:,iType,1)=force(:,iType,1)*REAL(atoms%relax(:,iType))
            ELSE
               force(:,iType,1)=0.0
            END IF
         END DO

         ! Add history
         CALL read_relax(pos,force,energies)

         ! Determine new positions
         IF (SIZE(energies)==1.OR.input%forcemix==0) THEN
            ! No history present simple step
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
            CALL juDFT_error('Unknown mixing scheme for forces', calledby='relaxation')
         END IF

         ! Check for convergence of forces/displacements
         maxDisp = 0.0
         l_conv = .TRUE.
         DO iType = 1, atoms%ntype
            IF (DOT_PRODUCT(force(:,iType,SIZE(force,3)),force(:,iType,SIZE(force,3)))>input%epsforce**2) l_conv=.FALSE.
            dispLength = SQRT(DOT_PRODUCT(displace(:,iType),displace(:,iType)))
            maxDisp = MAX(maxDisp,dispLength)
            IF (dispLength>input%epsdisp) l_conv=.FALSE.
         END DO

         ! Limit the maximal displacement in a single force relaxation step to limitDisp = 0.2 a_0.
         limitDisp = 0.2
         IF(maxDisp.GT.limitDisp) THEN
            displace(:,:) = limitDisp*displace(:,:) / maxDisp
         END IF

         ! New displacements relative to positions in inp.xml
         !CALL read_displacements(atoms,old_displace)
         CALL atoms_non_displaced%read_xml(xml)
         CALL xml%freeResources()
         CALL atoms_non_displaced%init(cell)

         DO iType = 1, atoms%ntype
            old_displace(:,iType) = atoms%pos(:,atoms%firstAtom(iType)) - &
                                    atoms_non_displaced%pos(:,atoms%firstAtom(iType))
         END DO

         tempAtoms = atoms_non_displaced
         tempDisplace = MATMUL(cell%bmat,totalDisplace)/tpi_const

         CALL rotate_to_all_sites(tempDisplace,atoms_non_displaced,cell,sym,dispAll)
         tempAtoms%taual(:,:)=atoms_non_displaced%taual(:,:)+dispAll(:,:)
         tempAtoms%pos=MATMUL(cell%amat,tempAtoms%taual)

         numDispReduce = 0
         overlap=1.0
         DO WHILE(ANY(overlap.GT.0.0))
            overlap = 0.0
            totalDisplace=displace+old_displace

            tempAtoms = atoms_non_displaced
            tempDisplace = MATMUL(cell%bmat,totalDisplace)/tpi_const

            CALL rotate_to_all_sites(tempDisplace,atoms_non_displaced,cell,sym,dispAll)
            tempAtoms%taual(:,:)=atoms_non_displaced%taual(:,:)+dispAll(:,:)
            tempAtoms%pos=MATMUL(cell%amat,tempAtoms%taual)
            CALL chkmt(tempAtoms,input,vacuum,cell ,.TRUE.,overlap=overlap)

            IF (ANY(overlap.GT.0.0)) THEN
               numDispReduce = numDispReduce + 1
               IF (numDispReduce.GE.3) THEN
                  CALL juDFT_warn("Strong MT spheres crash in structural relaxation")
               END IF
               displace(:,:) = 0.5 * displace(:,:)
               WRITE(oUnit,*) 'Automatically reducing atom displacements because MT spheres crash into each other!'
               WRITE(*,*) 'Automatically reducing atom displacements because MT spheres crash into each other!'
               ! TODO: Why two calls?
            END IF
         END DO

         ! Write relax file
         CALL write_relax(pos,force,energies,totalDisplace)

         ! Structure in xsf-format
         OPEN (55,file="struct-relax.xsf",status='replace')
         CALL xsf_WRITE_atoms(55,tempAtoms,input%film,cell%amat)
         CLOSE (55)
      END IF

#ifdef CPP_MPI
      CALL MPI_BCAST(l_conv,1,MPI_LOGICAL,0,fmpi%mpi_comm,ierr)
#endif

      IF (l_conv) THEN
         CALL judft_end("Structural relaxation: Done",fmpi%irank)
      ELSE
         CALL mixing_history_reset(fmpi)
         CALL judft_end("Structural relaxation: new displacements generated",fmpi%irank)
      END IF

   END SUBROUTINE relaxation

   SUBROUTINE simple_step(alpha,maxdisp,force,displace)

      IMPLICIT NONE

      REAL, INTENT(IN)  :: alpha,maxdisp
      REAL, INTENT(IN)  :: force(:,:,:)
      REAL, INTENT(OUT) :: displace(:,:)

      REAL :: corr

      displace = alpha*force(:,:,SIZE(force,3))
      corr=maxdisp/maxval(abs(displace))
      IF (corr<1.0) displace = corr*alpha*force(:,:,size(force,3))

   END SUBROUTINE simple_step

   SUBROUTINE simple_bfgs(pos,force,shift)
      !--------------------------------------------------------------------------
      !  Simple BFGS method to calculate shift out of old positions and forces
      !--------------------------------------------------------------------------
      USE m_constants

      REAL,INTENT(IN)  :: pos(:,:,:),force(:,:,:)
      REAL,INTENT(OUT) :: shift(:,:)

      INTEGER           :: n, i, j, hist_length, n_force
      REAL, ALLOCATABLE :: h(:,:)
      REAL, ALLOCATABLE :: p(:), y(:), v(:)
      REAL              :: py, yy, gamma

      n_force=3*SIZE(pos,2)
      ALLOCATE(h(n_force,n_force))
      ALLOCATE(p(n_force),y(n_force),v(n_force))

      ! Calculate approx. Hessian

      ! Initialize H
      h = 0.0

      DO n = 1,n_force
         h(n,n) = 1.0
      END DO

      ! Loop over all iterations (including current)
      hist_length=SIZE(pos,3)

      DO n = 2,hist_length
         ! Differences
         p(:) = RESHAPE(pos(:,:,n)-pos(:,:,n-1),(/SIZE(p)/))
         y(:) = RESHAPE(force(:,:,n-1)-force(:,:,n),(/SIZE(p)/))

         ! Get necessary inner products and H|y>
         py = DOT_PRODUCT(p,y)
         v = MATMUL(y,h)
         yy = DOT_PRODUCT(y,v)

         ! Check that update will leave h positive definite:
         IF (py <= 0.0) THEN
            WRITE (oUnit,*) '  bfgs: <p|y> < 0'
            WRITE (oUnit,*) '  check convergence of forces'
            ! Starting over with initial hessian
            h = 0.0
            DO j = 1,n_force
               h(j,j) = 1.0
            END   DO
            CYCLE
         ELSE
            ! Update h
            IF (n == 2) THEN
               gamma = py/yy
            ELSE
               gamma = 1.0
            END IF

            DO j = 1,n_force
               DO i = 1,n_force
                  h(i,j) = (h(i,j) - (v(i)*p(j)+p(i)*v(j))/py)*gamma + (1.+gamma*yy/py)*p(i)*p(j)/py
               END DO
            END DO
         END IF
      END DO

      y(:) = RESHAPE(force(:,:,hist_length),(/SIZE(p)/))
      shift = RESHAPE(MATMUL(y,h),SHAPE(shift))
   END SUBROUTINE simple_bfgs

   SUBROUTINE simple_cg(pos,force,shift)

      REAL, INTENT(IN)  :: pos(:,:,:),force(:,:,:)
      REAL, INTENT(OUT) :: shift(:,:)

      REAL               :: f1(3,SIZE(pos,2)),f2(3,SIZE(pos,2))
      REAL               :: dist(3,SIZE(pos,2))
      REAL               :: eps
      INTEGER            :: n_old, i, j

      eps = 1.0e-9

      n_old = SIZE(pos,3)-1

      dist(:,:) = pos(:,:,n_old+1)-pos(:,:,n_old)

      DO i = 1, SIZE(pos,2)
         DO j = 1, 3
            IF(ABS(dist(j,i)).LT.eps) dist(j,i) = eps ! To avoid calculation of 0.0/0.0 below.
         END DO
      END DO

      f1 = (force(:,:,n_old+1)-force(:,:,n_old))/dist
      f2 = force(:,:,n_old+1)-f1*pos(:,:,n_old+1)
      shift = -1.*f2/f1-force(:,:,n_old+1)
   END SUBROUTINE simple_cg

END MODULE m_relaxation
