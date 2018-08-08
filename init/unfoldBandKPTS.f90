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

    INTEGER :: i
    DO i =1,3
	p_cell%amat(1,i)=cell%amat(1,i)/banddos%s_cell_x
	p_cell%amat(2,i)=cell%amat(2,i)/banddos%s_cell_y 
	p_cell%amat(3,i)=cell%amat(3,i)/banddos%s_cell_z 
!	p_cell%amat(i,1)=cell%amat(i,1)/banddos%s_cell_x
!	p_cell%amat(i,2)=cell%amat(i,2)/banddos%s_cell_y 
!	p_cell%amat(i,3)=cell%amat(i,3)/banddos%s_cell_z 
    END DO
    CALL inv3(p_cell%amat,p_cell%bmat,p_cell%omtil)
    p_cell%bmat=p_cell%bmat*tpi_const
    p_cell%latnam=cell%latnam
  END SUBROUTINE  build_primitive_cell

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
    write(1088,*) 'banddos%unfoldband: ', banddos%unfoldband
    write(1088,*) 'brav. matrix: '
    write(1088,'(f15.8,f15.8,f15.8)') cell%amat(1,1), cell%amat(1,2), cell%amat(1,3)
    write(1088,'(f15.8,f15.8,f15.8)') cell%amat(2,1), cell%amat(2,2), cell%amat(2,3)
    write(1088,'(f15.8,f15.8,f15.8)') cell%amat(3,1), cell%amat(3,2), cell%amat(3,3)
    write(1088,*) 'brav. rez. matrix: '
    write(1088,'(f15.8,f15.8,f15.8)') cell%bmat(1,1), cell%bmat(1,2), cell%bmat(1,3)
    write(1088,'(f15.8,f15.8,f15.8)') cell%bmat(2,1), cell%bmat(2,2), cell%bmat(2,3)
    write(1088,'(f15.8,f15.8,f15.8)') cell%bmat(3,1), cell%bmat(3,2), cell%bmat(3,3)
    write(1088,*) ' primitive brav. matrix: '
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%amat(1,1), p_cell%amat(1,2), p_cell%amat(1,3)
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%amat(2,1), p_cell%amat(2,2), p_cell%amat(2,3)
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%amat(3,1), p_cell%amat(3,2), p_cell%amat(3,3)
    write(1088,*) 'primitive brav. rez. matrix: '
    write(89,'(3f15.8)') p_cell%bmat
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%bmat(1,1), p_cell%bmat(1,2), p_cell%bmat(1,3)
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%bmat(2,1), p_cell%bmat(2,2), p_cell%bmat(2,3)
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%bmat(3,1), p_cell%bmat(3,2), p_cell%bmat(3,3)
    write(1088,'(a,i7,a,i7)') 'kpts%nkpt',kpts%nkpt,'   p_kpts%nkpt',p_kpts%nkpt
    write(1088,*) kpts%specialPoints
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
    REAL    :: list(10,p_kpts%nkpt)  !cartesion coordinates for k,K,m
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
    eps_kpt = 0.00000001

    CALL inv3(cell%bmat,rez_inv_to_internal,rez_inv_det)
    write(1088,*) p_kpts%specialPoints
    write(333,'(3f15.8)')p_kpts%bk
    kpt_dist=0
    DO i= 1,size(list,2)
	!        pc_kpoint_c(1)=p_kpts%bk(1,i)*p_cell%bmat(1,1)+p_kpts%bk(2,i)*p_cell%bmat(1,2)+p_kpts%bk(3,i)*p_cell%bmat(1,3)
	!        pc_kpoint_c(2)=p_kpts%bk(1,i)*p_cell%bmat(2,1)+p_kpts%bk(2,i)*p_cell%bmat(2,2)+p_kpts%bk(3,i)*p_cell%bmat(2,3)
	!        pc_kpoint_c(3)=p_kpts%bk(1,i)*p_cell%bmat(3,1)+p_kpts%bk(2,i)*p_cell%bmat(3,2)+p_kpts%bk(3,i)*p_cell%bmat(3,3)
		pc_kpoint_c(1)=p_kpts%bk(1,i)*p_cell%bmat(1,1)+p_kpts%bk(2,i)*p_cell%bmat(2,1)+p_kpts%bk(3,i)*p_cell%bmat(3,1)
		pc_kpoint_c(2)=p_kpts%bk(1,i)*p_cell%bmat(1,2)+p_kpts%bk(2,i)*p_cell%bmat(2,2)+p_kpts%bk(3,i)*p_cell%bmat(3,2)
		pc_kpoint_c(3)=p_kpts%bk(1,i)*p_cell%bmat(1,3)+p_kpts%bk(2,i)*p_cell%bmat(2,3)+p_kpts%bk(3,i)*p_cell%bmat(3,3)
		list(1,i)=pc_kpoint_c(1)
		list(2,i)=pc_kpoint_c(2)
		list(3,i)=pc_kpoint_c(3)
	!!!!------- finding kpts in primitive rez. unit cell ----- 
	!	representation_found=.false.
	!m_loop:	DO m1= -banddos%s_cell_x,banddos%s_cell_x
	!		DO m2= -banddos%s_cell_y,banddos%s_cell_y
	!			DO m3= -banddos%s_cell_z,banddos%s_cell_z
	!				pc_kpoint_c(1)=list(1,i)-m1*cell%bmat(1,1)-m2*cell%bmat(1,2)-m3*cell%bmat(1,3)
	!				pc_kpoint_c(2)=list(2,i)-m1*cell%bmat(2,1)-m2*cell%bmat(2,2)-m3*cell%bmat(2,3)
	!				pc_kpoint_c(3)=list(3,i)-m1*cell%bmat(3,1)-m2*cell%bmat(3,2)-m3*cell%bmat(3,3)
	!!				IF (         (dot_product(pc_kpoint_c(:)+eps(:), cell%bmat(:,1)) >= 0).AND.((dot_product(pc_kpoint_c(:)+eps(:), cell%bmat(:,1)) < dot_product(cell%bmat(:,1), cell%bmat(:,1)))) &
	!!				     & .AND. (dot_product(pc_kpoint_c(:)+eps(:), cell%bmat(:,2)) >= 0).AND.((dot_product(pc_kpoint_c(:)+eps(:), cell%bmat(:,2)) < dot_product(cell%bmat(:,2), cell%bmat(:,2)))) &
	!!				     & .AND. (dot_product(pc_kpoint_c(:)+eps(:), cell%bmat(:,3)) >= 0).AND.((dot_product(pc_kpoint_c(:)+eps(:), cell%bmat(:,3)) < dot_product(cell%bmat(:,3), cell%bmat(:,3))))) THEN
	!				IF (all((matmul(rez_inv_to_internal,pc_kpoint_c)+eps(:))>=0).and.all((matmul(rez_inv_to_internal,pc_kpoint_c)+eps(:))<1)) THEN
	!					list(4,i)=pc_kpoint_c(1)
	!					list(5,i)=pc_kpoint_c(2)
	!					list(6,i)=pc_kpoint_c(3)
	!					list(7,i)=-m1
	!					list(8,i)=-m2
	!					list(9,i)=-m3
	!					representation_found=.true.
	!				END IF
	 !      			        IF (representation_found) EXIT m_loop
	!			END DO
	!		END DO
	!	END DO m_loop
	 !       IF (.not.representation_found) THEN
	  !      write(*,'(a,f15.8,f15.8,f15.8)') 'No representation found for the following kpoint:',list(1,i),list(2,i),list(3,i)
	   !     END IF
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
eps_kpt=0
	    sc_kpoint_i(1) = sc_kpoint_i(1) - m1 + eps_kpt
	    sc_kpoint_i(2) = sc_kpoint_i(2) - m2 + eps_kpt
	    sc_kpoint_i(3) = sc_kpoint_i(3) - m3 + eps_kpt
	    !sc_kpoint_i(:) = sc_kpoint_i(:) - 0.5
	    list(4,i)=sc_kpoint_i(1)
	    list(5,i)=sc_kpoint_i(2)
	    list(6,i)=sc_kpoint_i(3)
	    list(7,i)=m1
	    list(8,i)=m2
	    list(9,i)=m3 !this whole block is to move kpoints into first BZ within -0.5 to 0.5

	!  	kpts%bk(:,i)=matmul(rez_inv_to_internal,pc_kpoint_c)
	    kpts%bk(:,i)=list(4:6,i)
	
	IF (i>1) THEN
	kpt_dist=kpt_dist+sqrt(dot_product(list(1:3,i)-list(1:3,i-1),list(1:3,i)-list(1:3,i-1)))
	END IF
	list(10,i)=kpt_dist
    END DO
    write(91,'(3f15.8)') kpts%bk
    write(92,*) kpts%wtkpt
    ALLOCATE (kpts%sc_list(10,p_kpts%nkpt))
    kpts%sc_list=list
    write(90,'(10f15.8)') kpts%sc_list
  END SUBROUTINE find_supercell_kpts

 SUBROUTINE calculate_plot_w_n(banddos,cell,kpts,smat_unfold,zMat,lapw,i_kpt,jsp,eig,results,input,atoms)
	USE m_types
	USE m_juDFT
	USE m_inv3
	USE m_types_mpimat
        USE m_constants
	implicit none

        TYPE(t_input),INTENT(IN) :: input
        TYPE(t_atoms),INTENT(IN)     :: atoms
	TYPE(t_banddos),INTENT(IN)  :: banddos
	TYPE(t_results),INTENT(IN)  :: results
	TYPE(t_cell),INTENT(IN)     :: cell
	TYPE(t_kpts),INTENT(IN)     :: kpts
	CLASS(t_mat),INTENT(INOUT)  :: smat_unfold
	CLASS(t_mat),INTENT(IN)     :: zMat
	TYPE(t_lapw),INTENT(IN)     :: lapw
	TYPE(t_cell)      :: p_cell
	INTEGER, INTENT(IN)	    :: i_kpt,jsp
	REAL, INTENT(IN)	    :: eig(:)
	INTEGER :: i,j,k,l,n
	INTEGER :: na,n_i,nn,nk,nki,gi,lo
	REAL, ALLOCATABLE	::w_n(:)
	COMPLEX, ALLOCATABLE    ::w_n_c(:)
	REAL, ALLOCATABLE	::w_n_sum(:)
	COMPLEX, ALLOCATABLE    ::w_n_c_sum(:)
        LOGICAL :: method_rubel=.false.

	CALL build_primitive_cell(banddos,p_cell,cell)

        DO j = 1, lapw%nv(jsp)
           DO i = 1, j-1
              IF(smat_unfold%l_real) THEN
                 smat_unfold%data_r(j,i) = smat_unfold%data_r(i,j)
              ELSE
                 smat_unfold%data_c(j,i) = CONJG(smat_unfold%data_c(i,j))
              END IF
           END DO
        END DO
	IF (i_kpt==1) THEN
		IF (jsp==1) OPEN (679,file='bands_sc.1',status='unknown') !This is kind of my birthday 6 july 1992 (S.R.)
		IF (jsp==2) OPEN (680,file='bands_sc.2',status='unknown')
	END IF

!		write(*,*) 'real zmat size dim 1:', size(zMat%data_r,1), 'dim2:', size(zMat%data_r,2)
!		write(*,*) 'smat dim1', size(smat_unfold%data_r,1), 'dim2', size(smat_unfold%data_r,2),'data',smat_unfold%data_r(2,2)
!		write(222,'(234f15.8)') zMat%data_r
!		write(223,'(234f15.8)') smat_unfold%data_r

!	method_rubel=.true.    !this switch is to switch between overlap matrix and rubel method (without overlap matrix)

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
!		write(345,'(3I6)') lapw%gvec(:,:,jsp)
	write (*,*)results%ef
        write (*,*) i_kpt
	DO i=1,zMat%matsize2
		IF (method_rubel) THEN
			DO j=1,lapw%nv(jsp)
				IF (zmat%l_real) THEN
					w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat%data_r(j,i)
!						write(*,*) 'zMat is real'
				ELSE
					w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
!						write(*,*) 'zMat is complex'
				END IF
				IF ((modulo(lapw%gvec(1,j,jsp)+NINT(kpts%sc_list(7,i_kpt)),banddos%s_cell_x)==0).AND.&
				     &(modulo(lapw%gvec(2,j,jsp)+NINT(kpts%sc_list(8,i_kpt)),banddos%s_cell_y)==0).AND.&
				     &(modulo(lapw%gvec(3,j,jsp)+NINT(kpts%sc_list(9,i_kpt)),banddos%s_cell_z)==0)) THEN
					IF (zmat%l_real) THEN
						w_n(i)=w_n(i)+zMat%data_r(j,i)*zMat%data_r(j,i)
!							write(*,*) 'zMat is real'
					ELSE
						w_n_c(i)=w_n_c(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
!							write(*,*) 'zMat is complex'
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
								w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat%data_r(j,i)
							ELSE
								w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
							END IF
							IF ((modulo(lapw%gvec(1,gi,jsp)+NINT(kpts%sc_list(7,i_kpt)),banddos%s_cell_x)==0).AND.&
							   &(modulo(lapw%gvec(2,gi,jsp)+NINT(kpts%sc_list(8,i_kpt)),banddos%s_cell_y)==0).AND.&
							   &(modulo(lapw%gvec(3,gi,jsp)+NINT(kpts%sc_list(9,i_kpt)),banddos%s_cell_z)==0)) THEN
								IF (zmat%l_real) THEN
									w_n(i)=w_n(i)+zMat%data_r(j,i)*zMat%data_r(j,i)
								ELSE
									w_n_c(i)=w_n_c(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(j,i)
								END IF
							END IF
						END DO
					END DO
				END DO
			END DO
!--------------------------LO's finished----------------
		ELSE
			DO j=1,lapw%nv(jsp)
				DO k=1,zMat%matsize1
					IF (zmat%l_real) THEN
						w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat%data_r(k,i)*smat_unfold%data_r(j,k)
					ELSE
						w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(k,i)*smat_unfold%data_c(j,k)
					END IF
				END DO
				IF ((modulo(lapw%gvec(1,j,jsp)+NINT(kpts%sc_list(7,i_kpt)),banddos%s_cell_x)==0).AND.&
				   &(modulo(lapw%gvec(2,j,jsp)+NINT(kpts%sc_list(8,i_kpt)),banddos%s_cell_y)==0).AND.&
				   &(modulo(lapw%gvec(3,j,jsp)+NINT(kpts%sc_list(9,i_kpt)),banddos%s_cell_z)==0)) THEN
					DO k=1,zMat%matsize1
						IF (zmat%l_real) THEN
							w_n(i)=w_n(i)+zMat%data_r(j,i)*zMat%data_r(k,i)*smat_unfold%data_r(j,k)
						ELSE
							w_n_c(i)=w_n_c(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(k,i)*smat_unfold%data_c(j,k)
						END IF
					END DO
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
							DO k=1,zMat%matsize1
								IF (zmat%l_real) THEN
									w_n_sum(i)=w_n_sum(i)+zMat%data_r(j,i)*zMat%data_r(k,i)*smat_unfold%data_r(j,k)
								ELSE
									w_n_c_sum(i)=w_n_c_sum(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(k,i)*smat_unfold%data_c(j,k)
								END IF
							END DO
							IF ((modulo(lapw%gvec(1,gi,jsp)+NINT(kpts%sc_list(7,i_kpt)),banddos%s_cell_x)==0).AND.&
							   &(modulo(lapw%gvec(2,gi,jsp)+NINT(kpts%sc_list(8,i_kpt)),banddos%s_cell_y)==0).AND.&
							   &(modulo(lapw%gvec(3,gi,jsp)+NINT(kpts%sc_list(9,i_kpt)),banddos%s_cell_z)==0)) THEN
								DO k=1,zMat%matsize1
									IF (zmat%l_real) THEN
										w_n(i)=w_n(i)+zMat%data_r(j,i)*zMat%data_r(k,i)*smat_unfold%data_r(j,k)
									ELSE
										w_n_c(i)=w_n_c(i)+CONJG(zMat%data_c(j,i))*zMat%data_c(k,i)*smat_unfold%data_c(j,k)
									END IF
		    						END DO
							END IF
						END DO
					END DO
				END DO
			END DO
!--------------------------LO's finished----------------
		END IF
!		IF (method_rubel) THEN
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
!		ELSE
!			IF (zmat%l_real) THEN
!				IF (jsp==1) write(679,'(3f15.8)') kpt_dist, ((eig(i)-results%ef)*hartree_to_ev_const),w_n(i)
!				IF (jsp==2) write(680,'(3f15.8)') kpt_dist, ((eig(i)-results%ef)*hartree_to_ev_const),w_n(i)
!				IF ((w_n(i)>1).or.(w_n(i)<0)) write(*,*) 'w_n larger 1 or smaller 0', w_n(i), 'eigenvalue',eig(i)
!			ELSE
!				IF (jsp==1) write(679,'(4f15.8)') kpt_dist, ((eig(i)-results%ef)*hartree_to_ev_const),w_n_c(i)
!				IF (jsp==2) write(680,'(4f15.8)') kpt_dist, ((eig(i)-results%ef)*hartree_to_ev_const),w_n_c(i)
!				IF ((abs(w_n_c(i))>1).or.(real(w_n_c(i))<0)) write(*,*) 'w_n_c larger 1 or smaller 0', w_n_c(i), 'eigenvalue',eig(i)
!	        	END IF
!		END IF			
	END DO
	IF (i_kpt==kpts%nkpt) THEN
		IF (jsp==1) CLOSE (679)
		IF (jsp==input%jspins) THEN
			IF (jsp==2) CLOSE (680)
			CALL juDFT_error('Unfolded Bandstructure created succesfully - use band_sc.gnu to plot', calledby='calculate_plot_w_n')
		END IF
	END IF
 END SUBROUTINE
      	
      SUBROUTINE write_gnu_sc(nosyp,d,ssy,input)
      	USE m_types
	USE m_juDFT
      IMPLICIT NONE

      TYPE(t_input),INTENT(IN) :: input
      INTEGER, INTENT (IN) :: nosyp
      REAL,    INTENT (IN) :: d(nosyp)
      CHARACTER(len=1), INTENT (IN) :: ssy(nosyp)
      
      INTEGER n,aoff,adel
      CHARACTER(LEN=200) tempTitle
      aoff = iachar('a')-1
      adel = iachar('a')-iachar('A')
      !write(*,*) aoff,adel 

      OPEN (27,file='band_sc.gnu',status='unknown')
      WRITE (27,*) 'reset'
      WRITE (27,900)
      WRITE (27,901)
      WRITE (27,902)
      WRITE (27,903)
      WRITE(tempTitle,'(10a)') input%comment
      IF(TRIM(ADJUSTL(tempTitle)).EQ.'') THEN
         tempTitle = "Fleur Bandstructure"
      END IF
      WRITE (27,904) TRIM(ADJUSTL(tempTitle))
      DO n = 1, nosyp
        WRITE (27,905) d(n),d(n)
      ENDDO
      WRITE (27,906) d(1),d(nosyp)
!
! nomal labels
!
      IF (iachar(ssy(1)) < aoff ) THEN
        WRITE (27,907) ssy(1),d(1),achar(92)
      ELSE
        WRITE (27,907) " ",d(1),achar(92)
      ENDIF
      DO n = 2, nosyp-1
        IF (iachar(ssy(n)) < aoff ) THEN 
          WRITE (27,908) ssy(n),d(n),achar(92)
        ELSE
          WRITE (27,908) " ",d(n),achar(92)
        ENDIF
      ENDDO
      IF (iachar(ssy(nosyp)) < aoff ) THEN
        WRITE (27,909) ssy(nosyp),d(nosyp)
      ELSE
        WRITE (27,909) " ",d(nosyp)
      ENDIF
!
! greek labels
!
      DO n = 1, nosyp
        IF (iachar(ssy(n)) > aoff ) THEN
          WRITE (27,914) achar(iachar(ssy(n))-adel),d(n)
        ENDIF
      ENDDO
!
! now write the rest
!
      WRITE (27,910)
      WRITE (27,*) 'set palette model RGB'
      WRITE (27,*) 'set palette defined (-2 "black", -1 "white" ,0 "white",',achar(92)
      WRITE (27,*) '0.67 "light-blue",1 "blue")'
      WRITE (27,*) 'set cbrange [-2:1]'
      WRITE (27,*) 'unset colorbox'
      WRITE (27,*) 'size1(x)=0.9*x**(0.4)'
      WRITE (27,*) 'color1(x)=0.3+x/2.4'
      WRITE (27,*) 'size2(x)=0.35*(1-x**(0.01))'
      WRITE (27,*) 'color2(x)=1.15*(x-1)'
      WRITE (27,911) d(nosyp)+0.00001,achar(92)
      IF (input%jspins == 2) THEN
	WRITE (27,912) achar(92)
	WRITE (27,916) achar(92)
      END IF
      WRITE (27,913) achar(92)
      WRITE (27,915)
      CLOSE (27)

 900  FORMAT ('set terminal postscript enhanced color "Times-Roman" 20')
 901  FORMAT ('set xlabel ""')
 902  FORMAT ('set ylabel "E - E_F (eV)"')
 903  FORMAT ('set nokey')
 904  FORMAT ('set title "',a,'"')
 905  FORMAT ('set arrow from',f9.5,', -9.0 to',f9.5,',  5.0 nohead')
 906  FORMAT ('set arrow from',f9.5,', 0.0 to',f9.5,', 0.0 nohead lt 3')
 907  FORMAT ('set xtics ("',a1,'"',f9.5,', ',a)
 908  FORMAT ('           "',a1,'"',f9.5,', ',a)
 909  FORMAT ('           "',a1,'"',f9.5,'  )')
 910  FORMAT ('set ytics -8,2,4')
 911  FORMAT ('plot [0:',f9.5,'] [-9:5] ',a)
 912  FORMAT ('"bands_sc.2" using 1:($2-6.00):(size1($3)):(color1($3))  w p pt 7 ps variable lc palette, ',a)
 916  FORMAT ('"bands_sc.2" using 1:($2-6.00):(size2($3)):(color2($3)) w p pt 7 ps variable lc palette,',a)
 913  FORMAT ('"bands_sc.1" using 1:($2-6.00):(size1($3)):(color1($3))  w p pt 7 ps variable lc palette, ',a)
 915  FORMAT ('"bands_sc.1" using 1:($2-6.00):(size2($3)):(color2($3)) w p pt 7 ps variable lc palette')
 914  FORMAT ('set label "',a1,'" at ',f9.5,', -9.65 center font "Symbol,20"')
      END SUBROUTINE write_gnu_sc
END MODULE m_unfold_band_kpts
