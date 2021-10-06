!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_unfold_band_kpts

CONTAINS

  SUBROUTINE build_primitive_cell(banddos,p_cell,cell)
    USE m_types
    USE m_inv3
    USE m_constants, ONLY : tpi_const
    implicit none
    TYPE(t_banddos),INTENT(IN)  :: banddos
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_cell),INTENT(INOUT)  :: p_cell
	REAL    :: unfold(3,3)  !this variable should be given in the input xml
	REAL    :: inv_unfold(3,3)
    REAL    :: inv_unfold_det
	INTEGER :: i
	unfold=banddos%unfoldTransMat
	unfold(1,1)=banddos%unfoldTransMat(1,1)*banddos%s_cell_x
	unfold(2,2)=banddos%unfoldTransMat(2,2)*banddos%s_cell_y
	unfold(3,3)=banddos%unfoldTransMat(3,3)*banddos%s_cell_z
	
	CALL inv3(unfold,inv_unfold,inv_unfold_det)

    DO i =1,3
	p_cell%amat(:,i)=matmul(inv_unfold,cell%amat(:,i))
    END DO
    CALL inv3(p_cell%amat,p_cell%bmat,p_cell%omtil)
    p_cell%bmat=p_cell%bmat*tpi_const
  END SUBROUTINE  build_primitive_cell
!---------- the following routines are not used anymore (but instructive)-----
  SUBROUTINE unfold_band_kpts(banddos,p_cell,cell,p_kpts,kpts)
    USE m_types
    USE m_inv3
    USE m_constants, ONLY : tpi_const

    implicit none

    TYPE(t_banddos),INTENT(IN)  :: banddos
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_cell),INTENT(INOUT)  :: p_cell
    TYPE(t_kpts),INTENT(INOUT)  :: p_kpts
    TYPE(t_kpts),INTENT(INOUT)  :: kpts

    CALL build_primitive_cell(banddos,p_cell,cell)

    p_kpts=kpts
    !write(1088,*) 'banddos%unfoldband: ', banddos%unfoldband
    !write(1088,*) 'brav. matrix: '
    !write(1088,'(f15.8,f15.8,f15.8)') cell%amat(1,1), cell%amat(1,2), cell%amat(1,3)
    !write(1088,'(f15.8,f15.8,f15.8)') cell%amat(2,1), cell%amat(2,2), cell%amat(2,3)
    !write(1088,'(f15.8,f15.8,f15.8)') cell%amat(3,1), cell%amat(3,2), cell%amat(3,3)
  END SUBROUTINE unfold_band_kpts

  SUBROUTINE find_supercell_kpts(banddos,p_cell,cell,p_kpts,kpts)
    USE m_types
    USE m_juDFT
    USE m_inv3
    implicit none

    TYPE(t_banddos),INTENT(IN)  :: banddos
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_cell),INTENT(IN)     :: p_cell
    TYPE(t_kpts),INTENT(IN)     :: p_kpts
    TYPE(t_kpts),INTENT(INOUT)  :: kpts

    INTEGER :: i,m1,m2,m3
    REAL    :: rez_inv_to_internal(3,3)
    REAL    :: rez_inv_det
    REAL    :: list(13,p_kpts%nkpt)  !cartesion coordinates for k,K,m
    REAL    :: pc_kpoint_i(3)    !primitive cell kpoint internal
    REAL    :: sc_kpoint_i(3)    !super cell kpoint internal
    REAL    :: pc_kpoint_c(3)    !primitive cell kpoint cartesian
    REAL    :: sc_kpoint_c(3)    !super cell kpoint cartesian
    REAL    :: eps(3)
    REAL    :: eps_r, eps_kpt
    LOGICAL :: representation_found
    REAL    ::kpt_dist

    eps = 1.0e-10
    eps_r = 0.000000001

    CALL inv3(cell%bmat,rez_inv_to_internal,rez_inv_det)
    !write(1088,*) p_kpts%specialPoints
    !write(333,'(3f15.8)')p_kpts%bk
    kpt_dist=0
    DO i= 1,size(list,2)
		pc_kpoint_c(1)=p_kpts%bk(1,i)*p_cell%bmat(1,1)+p_kpts%bk(2,i)*p_cell%bmat(2,1)+p_kpts%bk(3,i)*p_cell%bmat(3,1)
		pc_kpoint_c(2)=p_kpts%bk(1,i)*p_cell%bmat(1,2)+p_kpts%bk(2,i)*p_cell%bmat(2,2)+p_kpts%bk(3,i)*p_cell%bmat(3,2)
		pc_kpoint_c(3)=p_kpts%bk(1,i)*p_cell%bmat(1,3)+p_kpts%bk(2,i)*p_cell%bmat(2,3)+p_kpts%bk(3,i)*p_cell%bmat(3,3)
		list(1,i)=pc_kpoint_c(1)
		list(2,i)=pc_kpoint_c(2)
		list(3,i)=pc_kpoint_c(3)
	   !----------------------- method internal coordintes --------------------
	    sc_kpoint_i(:)=matmul(pc_kpoint_c,rez_inv_to_internal)
	    pc_kpoint_i(:)=p_kpts%bk(1:3,i)
	    !sc_kpoint_i(:) = sc_kpoint_i(:) + 0.5
	    m1 = FLOOR(sc_kpoint_i(1))
	    m2 = FLOOR(sc_kpoint_i(2))
	    m3 = FLOOR(sc_kpoint_i(3))
	    m1=0
	    m2=0
	    m3=0
	    sc_kpoint_i(1) = sc_kpoint_i(1) - m1
	    sc_kpoint_i(2) = sc_kpoint_i(2) - m2
	    sc_kpoint_i(3) = sc_kpoint_i(3) - m3
	    !sc_kpoint_i(:) = sc_kpoint_i(:) - 0.5
	    list(4,i)=sc_kpoint_i(1)
	    list(5,i)=sc_kpoint_i(2)
	    list(6,i)=sc_kpoint_i(3)
	    list(7,i)=m1
	    list(8,i)=m2
	    list(9,i)=m3 !this whole block is to move kpoints into first BZ within -0.5 to 0.5

	!  	kpts%bk(:,i)=matmul(rez_inv_to_internal,pc_kpoint_c)
	    !-------------saving old kpts----------
	    list(11:13,i)=kpts%bk(:,i)
  	    !------finished---------
	    kpts%bk(:,i)=list(4:6,i)

	IF (i>1) THEN
	kpt_dist=kpt_dist+sqrt(dot_product(list(1:3,i)-list(1:3,i-1),list(1:3,i)-list(1:3,i-1)))
	END IF
	list(10,i)=kpt_dist
    END DO
    !write(91,'(3f15.8)') kpts%bk
    !write(92,*) kpts%wtkpt
    ALLOCATE (kpts%sc_list(13,p_kpts%nkpt))
    kpts%specialPointIndices(:) = p_kpts%specialPointIndices(:)
    kpts%sc_list=list
    write(90,'(10f15.8)') kpts%sc_list
  END SUBROUTINE find_supercell_kpts
!----------------------------------------------------------------
 SUBROUTINE calculate_plot_w_n(banddos,cell,kpts,zMat,lapw,i_kpt,jsp,eig,results,input,atoms,unfoldingBuffer,fmpi,l_soc,smat_unfold,zso)
	USE m_types
	USE m_juDFT
	USE m_inv3
	USE m_types_mpimat
    USE m_constants
	implicit none

    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_atoms),INTENT(IN)     :: atoms
	TYPE(t_banddos),INTENT(IN)   :: banddos
	TYPE(t_results),INTENT(INOUT)  :: results
	TYPE(t_cell),INTENT(IN)     :: cell
	TYPE(t_kpts),INTENT(IN)     :: kpts
	CLASS(t_mat),OPTIONAL,INTENT(INOUT)  :: smat_unfold
	CLASS(t_mat),INTENT(IN)     :: zMat
	TYPE(t_lapw),INTENT(IN)     :: lapw
    TYPE(t_mpi),INTENT(IN)      :: fmpi
	TYPE(t_cell)            :: p_cell
	LOGICAL, INTENT(IN)     :: l_soc
	INTEGER, INTENT(IN)	    :: i_kpt,jsp
	REAL, INTENT(IN)	    :: eig(:)
	COMPLEX, INTENT(IN), OPTIONAL     ::zso(:,:,:)
    COMPLEX, INTENT(INOUT)  :: unfoldingBuffer(:,:,:)
	INTEGER :: i,j,k,l,n
	INTEGER :: na,n_i,nn,nk,nki,gi,lo
	REAL, ALLOCATABLE	  ::w_n(:)
	COMPLEX, ALLOCATABLE  ::w_n_c(:)
	REAL, ALLOCATABLE	  ::w_n_sum(:)
	COMPLEX, ALLOCATABLE  ::w_n_c_sum(:)
    LOGICAL :: method_rubel = .FALSE. 
    LOGICAL :: write_to_file = .false.
    CLASS(t_mat), ALLOCATABLE :: zMat_s
	REAL    :: unfold(3,3)  !this variable should be given in the input xml
	REAL    :: multiple(3)
	REAL    :: inv_unfold(3,3)
    REAL    :: inv_unfold_det
    REAL    :: eps_r=0.0000000001
!---------combining matrix input and unfolding factor input-----------	
	unfold=banddos%unfoldTransMat
	unfold(1,1)=banddos%unfoldTransMat(1,1)*banddos%s_cell_x
	unfold(2,2)=banddos%unfoldTransMat(2,2)*banddos%s_cell_y
	unfold(3,3)=banddos%unfoldTransMat(3,3)*banddos%s_cell_z
    CALL inv3(unfold,inv_unfold,inv_unfold_det)
        method_rubel = .NOT.banddos%unfoldUseOlap

	CALL build_primitive_cell(banddos,p_cell,cell)

	IF (.not. method_rubel) THEN
		DO j = 1, lapw%nv(jsp)
		  DO i = 1, j-1
	      		IF(smat_unfold%l_real) THEN
				smat_unfold%data_r(j,i) = smat_unfold%data_r(i,j)
	      		ELSE
				smat_unfold%data_c(j,i) = CONJG(smat_unfold%data_c(i,j))
	      		END IF
		   END DO
		END DO
	END IF
!   	write_to_file=.true.
	IF (write_to_file) THEN
		IF (i_kpt==1) THEN
			IF (jsp==1) OPEN (679,file='bands_sc_old.1',status='unknown') !This is kind of my birthday 6 july 1992 (S.R.)
			IF (jsp==2) OPEN (680,file='bands_sc_old.2',status='unknown')
		END IF
	END IF

	IF (zmat%l_real) THEN
		ALLOCATE(w_n(zMat%matsize2))
	        w_n = 0
!	    IF (method_rubel) THEN
		ALLOCATE(w_n_sum(zMat%matsize2))
	        w_n_sum = 0
!	    END IF
	ELSE
		ALLOCATE(w_n_c(zMat%matsize2))
		w_n_c=0
!	    IF (method_rubel) THEN
		ALLOCATE(w_n_c_sum(zMat%matsize2))
		w_n_c_sum=0
!	    END IF
	END IF
!---------create zmat_s--- smat*zmat---------------------
	select type(zMat)
		type is (t_mat)
		allocate(t_mat::zMat_s)
		select type(zMat_s)
             		type is (t_mat)
	     		zMat_s=zMat
		end select
		type is (t_mpimat)
		allocate(t_mpimat::zMat_s)
		select type(zMat_s)
             		type is (t_mpimat)
	     		zMat_s=zMat
		end select
	end select
!---------------------------------------------------------
!		write(345,'(3I6)') lapw%gvec(:,:,jsp)
	write (*,*)results%ef
    write (*,*) i_kpt
	IF (.not. method_rubel) THEN
!          IF (fmpi%n_size==1) THEN
!             call smat_unfold%multiply(zMat,zMat_s)
!          ELSE
!             call smat_unfold%mpimat_multiply(zMat,zMat_s)
!          ENDIF
	   call smat_unfold%multiply(zMat,zMat_s)
    END IF
   !$omp parallel default(none) private(j,n_i,nn,na,lo,nk,nki,gi,multiple) shared(zmat,method_rubel,jsp,lapw,w_n_sum, w_n_c_sum,inv_unfold,eps_r,w_n,w_n_c,atoms,zmat_s,write_to_file,i_kpt,kpts,eig,results,unfoldingBuffer,zso,l_soc)
   !$omp do
	DO i=1,zMat%matsize2
!	        write (*,*) 'here i work 1 -', i
		IF (method_rubel) THEN
			DO j=1,lapw%nv(jsp)
				IF (zmat%l_real) THEN
					w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat%data_r(j,i)
				!	write(*,*) 'zMat is real'
				ELSE
					IF (l_soc) THEN
						w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zso(j,i,jsp))*zso(j,i,jsp)
					ELSE 
						w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
					END IF
				!	write(*,*) 'zMat is complex', j,i
				END IF
				multiple=matmul(inv_unfold,lapw%gvec(:,j,jsp))
				IF ((abs(modulo(multiple(1),1.0))<eps_r).AND.&
					&(abs(modulo(multiple(2),1.0))<eps_r).AND.&
					&(abs(modulo(multiple(3),1.0))<eps_r)) THEN    
					IF (zmat%l_real) THEN
						w_n(i)=w_n(i)+zMat%data_r(j,i)*zMat%data_r(j,i)
						!write(*,*) 'zMat is real'
					ELSE
						IF (l_soc) THEN
							w_n_c(i)=w_n_c(i)+CONJG(zso(j,i,jsp))*zso(j,i,jsp)
						ELSE
							w_n_c(i)=w_n_c(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
						END IF
						!write(*,*) 'zMat is complex - restricted sum'
					END IF
				END IF
			END DO
!------------------LO's------------------------
			na=0
			!write(*,*) 'start lo', i
			DO n_i=1,atoms%ntype
				DO nn=1,atoms%neq(n_i)
					na=na+1
					DO lo=1,atoms%nlo(n_i)
						nk=lapw%nkvec(lo,na)
						DO nki=1,nk
							gi=lapw%kvec(nki,lo,na)
							j=lapw%nv(jsp)+lapw%index_lo(lo,na)+nki
							IF (zmat%l_real) THEN
								w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat%data_r(j,i)
							ELSE
								w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
							END IF
							multiple=matmul(inv_unfold,lapw%gvec(:,gi,jsp))
							IF ((abs(modulo(multiple(1),1.0))<eps_r).AND.&
								&(abs(modulo(multiple(2),1.0))<eps_r).AND.&
								&(abs(modulo(multiple(3),1.0))<eps_r)) THEN 
								IF (zmat%l_real) THEN
									w_n(i)=w_n(i)+zMat%data_r(j,i)*zMat%data_r(j,i)
									!write(*,*) zMat%data_r(j,i)*zMat%data_r(j,i)
								ELSE
									w_n_c(i)=w_n_c(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
									!write (*,*) CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
								END IF
							END IF
						END DO
					END DO
				END DO
			END DO
			!write(*,*) 'finished lo', i
!--------------------------LO's finished----------------
		ELSE
			!write (*,*) 'start else'
			!write (*,*) 'lapw%nv',lapw%nv(jsp),'j',j
			!DO j=1,lapw%nv(jsp)
			!        write (*,*) 'test loop', j
			!END DO
			DO j=1,lapw%nv(jsp)
				!write (*,*) 'start do',j
				!DO k=1,zMat%matsize1
				IF (zmat%l_real) THEN
					!w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat%data_r(k,i)*smat_unfold%data_r(j,k)
					w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat_s%data_r(j,i)
					!write (*,*) 'weight sum real'
				ELSE
					!w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(k,i)*smat_unfold%data_c(j,k)
					w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat_s%data_c(j,i)
				END IF
!				END DO
!				write (*,*) lapw%gvec(:,j,jsp)
!				write (*,*) kpts%sc_list(:,i_kpt)
!				write (*,*) banddos%s_cell_x,banddos%s_cell_y,banddos%s_cell_z
				!CALL juDFT_error('debugging stop, unfolding')
				multiple=matmul(inv_unfold,lapw%gvec(:,j,jsp))
				IF ((abs(modulo(multiple(1),1.0))<eps_r).AND.&
					&(abs(modulo(multiple(2),1.0))<eps_r).AND.&
					&(abs(modulo(multiple(3),1.0))<eps_r)) THEN  
					IF (zmat%l_real) THEN
						w_n(i)=w_n(i)+zMat%data_r(j,i)*zMat_s%data_r(j,i)
					ELSE
						w_n_c(i)=w_n_c(i)+CONJG(zMat%data_c(j,i))*zMat_s%data_c(j,i)
						write (*,*) CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
					END IF
				END IF
			END DO
!------------------LO's------------------------
				na=0
				DO n_i=1,atoms%ntype
					DO nn=1,atoms%neq(n_i)
						na=na+1
						DO lo=1,atoms%nlo(n_i)
						nk=lapw%nkvec(lo,na)
						DO nki=1,nk
							gi=lapw%kvec(nki,lo,na)
							j=lapw%nv(jsp)+lapw%index_lo(lo,na)+nki
								IF (zmat%l_real) THEN
									w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat_s%data_r(j,i)
								ELSE
									w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat_s%data_c(j,i)
								END IF
							multiple=matmul(inv_unfold,lapw%gvec(:,gi,jsp))
							IF ((abs(modulo(multiple(1),1.0))<eps_r).AND.&
								&(abs(modulo(multiple(2),1.0))<eps_r).AND.&
								&(abs(modulo(multiple(3),1.0))<eps_r)) THEN
								IF (zmat%l_real) THEN
									w_n(i)=w_n(i)+zMat%data_r(j,i)*zMat_s%data_r(j,i)
								ELSE
									w_n_c(i)=w_n_c(i)+CONJG(zMat%data_c(j,i))*zMat_s%data_c(j,i)
								END IF
							END IF
						END DO
					END DO
				END DO
			END DO
!--------------------------LO's finished----------------
		END IF
		!IF (method_rubel) THEN
		IF (write_to_file) THEN
			IF (zmat%l_real) THEN
				IF (w_n(i)/w_n_sum(i)<0) w_n(i)=0   ! delete negative entries
				IF (jsp==1) write(679,'(3f15.8)') kpts%sc_list(10,i_kpt), ((eig(i)-results%ef)*hartree_to_ev_const),w_n(i)/w_n_sum(i)
				IF (jsp==2) write(680,'(3f15.8)') kpts%sc_list(10,i_kpt), ((eig(i)-results%ef)*hartree_to_ev_const),w_n(i)/w_n_sum(i)
				IF ((w_n(i)/w_n_sum(i)>1).or.(w_n(i)/w_n_sum(i)<0)) write(*,*) 'w_n/sum larger 1 or smaller 0', w_n(i)/w_n_sum(i), 'eigenvalue',eig(i)
			ELSE
				IF (real(w_n_c(i))<0) w_n_c(i)=0    ! delete negative entries
				IF (jsp==1) write(679,'(4f15.8)') kpts%sc_list(10,i_kpt), ((eig(i)-results%ef)*hartree_to_ev_const),w_n_c(i)/w_n_c_sum(i)
				IF (jsp==2) write(680,'(4f15.8)') kpts%sc_list(10,i_kpt), ((eig(i)-results%ef)*hartree_to_ev_const),w_n_c(i)/w_n_c_sum(i)
				IF ((abs(w_n_c(i)/w_n_c_sum(i))>1).or.(real(w_n_c(i))<0)) write(*,*) 'w_n_c/sum larger 1 or smaller 0', w_n_c(i)/w_n_c_sum(i), 'eigenvalue',eig(i)
			END IF
		END IF
		IF (zmat%l_real) THEN
			IF (w_n(i)/w_n_sum(i)<0) w_n(i)=0   ! delete negative entries
			unfoldingBuffer(i,i_kpt,jsp)=w_n(i)/w_n_sum(i)
			IF ((w_n(i)/w_n_sum(i)>1).or.(w_n(i)/w_n_sum(i)<0)) write(*,*) 'w_n/sum larger 1 or smaller 0', w_n(i)/w_n_sum(i), 'eigenvalue',eig(i)
		ELSE
			IF (real(w_n_c(i))<0) w_n_c(i)=0    ! delete negative entries
			unfoldingBuffer(i,i_kpt,jsp)=w_n_c(i)/w_n_c_sum(i)
			IF ((abs(w_n_c(i)/w_n_c_sum(i))>1).or.(real(w_n_c(i))<0)) write(*,*) 'w_n_c/sum larger 1 or smaller 0', w_n_c(i)/w_n_c_sum(i), 'eigenvalue',eig(i)
		END IF
	END DO
   !$omp end do
   !$omp end parallel
	write (*,*) 'finished',i_kpt
	IF (i_kpt==kpts%nkpt) THEN
		IF (write_to_file .AND. jsp==1) CLOSE (679)
		IF (jsp==input%jspins) THEN
			IF (write_to_file .AND. jsp==2) CLOSE (680)
			!kpts%bk(:,:)=kpts%sc_list(11:13,:)
			write(*,*) 'Unfolded Bandstructure calculated succesfully, calledby=calculate_plot_w_n'
			!CALL juDFT_error('Unfolded Bandstructure created succesfully - use band_sc.gnu to plot', calledby='calculate_plot_w_n')
		END IF
	END IF
 END SUBROUTINE

SUBROUTINE write_band_sc(banddos,cell,kpts,results,eFermiPrev)
     USE m_types
     USE m_juDFT
	 USE m_constants
	 USE m_inv3
     IMPLICIT NONE
	TYPE(t_results),INTENT(IN)  :: results
	TYPE(t_banddos),INTENT(IN)  :: banddos
	TYPE(t_kpts),INTENT(IN)     :: kpts
    REAL, INTENT(IN) :: eFermiPrev
	INTEGER :: i,i_kpt,jsp
	TYPE(t_cell),INTENT(IN)     :: cell
	TYPE(t_cell) :: p_cell


	REAL    :: kpt_dist
	REAL    :: list(4,kpts%nkpt)
!-------------build primitive cell ----------
	p_cell=cell
	DO i =1,3
		p_cell%amat(1,i)=cell%amat(1,i)/banddos%s_cell_x
		p_cell%amat(2,i)=cell%amat(2,i)/banddos%s_cell_y
		p_cell%amat(3,i)=cell%amat(3,i)/banddos%s_cell_z
	END DO
		CALL inv3(p_cell%amat,p_cell%bmat,p_cell%omtil)
		p_cell%bmat=p_cell%bmat*tpi_const

!-------------- calculate distance ------------
	kpt_dist=0
	DO i=1,size(list,2)
		list(1,i)=kpts%bk(1,i)*p_cell%bmat(1,1)+kpts%bk(2,i)*p_cell%bmat(2,1)+kpts%bk(3,i)*p_cell%bmat(3,1)
		list(2,i)=kpts%bk(1,i)*p_cell%bmat(1,2)+kpts%bk(2,i)*p_cell%bmat(2,2)+kpts%bk(3,i)*p_cell%bmat(3,2)
		list(3,i)=kpts%bk(1,i)*p_cell%bmat(1,3)+kpts%bk(2,i)*p_cell%bmat(2,3)+kpts%bk(3,i)*p_cell%bmat(3,3)
		IF (i>1) THEN
			kpt_dist=kpt_dist+sqrt(dot_product(list(1:3,i)-list(1:3,i-1),list(1:3,i)-list(1:3,i-1)))
		END IF
		list(4,i)=kpt_dist
	END DO
!--------------------------------------
	OPEN (679,file='bands_sc.1',status='unknown') !This is kind of my birthday 6 july 1992 (S.R.)
	IF (SIZE(results%unfolding_weights,3)==2) OPEN (680,file='bands_sc.2',status='unknown')
        DO jsp=1,SIZE(results%unfolding_weights,3)
		DO i_kpt=1,SIZE(results%unfolding_weights,2)
			DO i=1,results%neig(i_kpt,jsp)
				IF (jsp==1) write(679,'(4f15.8)') list(4,i_kpt), ((results%eig(i,i_kpt,1)-eFermiPrev)*hartree_to_ev_const),results%unfolding_weights(i,i_kpt,1)
				IF (jsp==2) write(680,'(4f15.8)') list(4,i_kpt), ((results%eig(i,i_kpt,2)-eFermiPrev)*hartree_to_ev_const),results%unfolding_weights(i,i_kpt,2)
			END DO
		END DO
	END DO
	CLOSE (679)
	IF (SIZE(results%unfolding_weights,3)==2) CLOSE (680)
	write(*,*) 'Unfolded Bandstructure written succesfully - use band_sc.gnu to plot, calledby=write_band_sc',eFermiPrev
END SUBROUTINE

!---- new subroutine in gnuplot.F90
!    SUBROUTINE write_gnu_sc_old(nosyp,d,ssy,input)
!    USE m_types
!	USE m_juDFT
!      IMPLICIT NONE
!
!      TYPE(t_input),INTENT(IN) :: input
!      INTEGER, INTENT (IN) :: nosyp
!      REAL,    INTENT (IN) :: d(nosyp)
!      CHARACTER(len=1), INTENT (IN) :: ssy(nosyp)
!
!      INTEGER n,aoff,adel
!      CHARACTER(LEN=200) tempTitle
!      aoff = iachar('a')-1
!      adel = iachar('a')-iachar('A')
!      !write(*,*) aoff,adel
!
!      OPEN (27,file='band_sc.gnu',status='unknown')
!      WRITE (27,*) 'reset'
!      WRITE (27,900)
!      WRITE (27,901)
!      WRITE (27,902)
!      WRITE (27,903)
!      WRITE(tempTitle,'(10a)') input%comment
!      IF(TRIM(ADJUSTL(tempTitle)).EQ.'') THEN
!         tempTitle = "Fleur Bandstructure"
!      END IF
!      WRITE (27,904) TRIM(ADJUSTL(tempTitle))
!      DO n = 1, nosyp
!        WRITE (27,905) d(n),d(n)
!      ENDDO
!      WRITE (27,906) d(1),d(nosyp)
!!
!! nomal labels
!!
!      IF (iachar(ssy(1)) < aoff ) THEN
!        WRITE (27,907) ssy(1),d(1),achar(92)
!      ELSE
!        WRITE (27,907) " ",d(1),achar(92)
!      ENDIF
!      DO n = 2, nosyp-1
!        IF (iachar(ssy(n)) < aoff ) THEN
!          WRITE (27,908) ssy(n),d(n),achar(92)
!        ELSE
!          WRITE (27,908) " ",d(n),achar(92)
!        ENDIF
!      ENDDO
!      IF (iachar(ssy(nosyp)) < aoff ) THEN
!        WRITE (27,909) ssy(nosyp),d(nosyp)
!      ELSE
!        WRITE (27,909) " ",d(nosyp)
!      ENDIF
!!
!! greek labels
!!
!      DO n = 1, nosyp
!        IF (iachar(ssy(n)) > aoff ) THEN
!          WRITE (27,914) achar(iachar(ssy(n))-adel),d(n)
!        ENDIF
!      ENDDO
!!
!! now write the rest
!!
!      WRITE (27,910)
!      WRITE (27,*) 'set palette model RGB'
!      WRITE (27,*) 'set palette defined (-2 "black", -1 "white" ,0 "white",',achar(92)
!      WRITE (27,*) '0.67 "light-blue",1 "blue")'
!      WRITE (27,*) 'set cbrange [-2:1]'
!      WRITE (27,*) 'unset colorbox'
!      WRITE (27,*) 'size1(x)=0.9*x**(0.4)'
!      WRITE (27,*) 'color1(x)=0.3+x/2.4'
!      WRITE (27,*) 'size2(x)=0.35*(1-x**(0.01))'
!      WRITE (27,*) 'color2(x)=1.15*(x-1)'
!      WRITE (27,*) 'e_f=0.000000 #fermi energy is already corrected when using hdf5'
!      WRITE (27,911) d(nosyp)+0.00001,achar(92)
!      IF (input%jspins == 2) THEN
!	WRITE (27,912) achar(92)
!	WRITE (27,916) achar(92)
!      END IF
!      WRITE (27,913) achar(92)
!      WRITE (27,915)
!      CLOSE (27)
!
! 900  FORMAT ('set terminal postscript enhanced color "Times-Roman" 20')
! 901  FORMAT ('set xlabel ""')
! 902  FORMAT ('set ylabel "E - E_F (eV)"')
! 903  FORMAT ('set nokey')
! 904  FORMAT ('set title "',a,'"')
! 905  FORMAT ('set arrow from',f9.5,', -9.0 to',f9.5,',  5.0 nohead')
! 906  FORMAT ('set arrow from',f9.5,', 0.0 to',f9.5,', 0.0 nohead lt 3')
! 907  FORMAT ('set xtics ("',a1,'"',f9.5,', ',a)
! 908  FORMAT ('           "',a1,'"',f9.5,', ',a)
! 909  FORMAT ('           "',a1,'"',f9.5,'  )')
! 910  FORMAT ('set ytics -8,2,4')
! 911  FORMAT ('plot [0:',f9.5,'] [-9:5] ',a)
! 912  FORMAT ('"bands_sc.2" using 1:($2-e_f):(size1($3)):(color1($3))  w p pt 7 ps variable lc palette, ',a)
! 916  FORMAT ('"bands_sc.2" using 1:($2-e_f):(size2($3)):(color2($3)) w p pt 7 ps variable lc palette,',a)
! 913  FORMAT ('"bands_sc.1" using 1:($2-e_f):(size1($3)):(color1($3))  w p pt 7 ps variable lc palette, ',a)
! 915  FORMAT ('"bands_sc.1" using 1:($2-e_f):(size2($3)):(color2($3)) w p pt 7 ps variable lc palette')
! 914  FORMAT ('set label "',a1,'" at ',f9.5,', -9.65 center font "Symbol,20"')
!      END SUBROUTINE write_gnu_sc_old
END MODULE m_unfold_band_kpts
