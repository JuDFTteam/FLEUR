MODULE m_unfold_band_kpts

CONTAINS

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
   
    INTEGER :: i
    DO i =1,3
	p_cell%amat(1,i)=cell%amat(1,i)/banddos%s_cell_x
	p_cell%amat(2,i)=cell%amat(2,i)/banddos%s_cell_y 
	p_cell%amat(3,i)=cell%amat(3,i)/banddos%s_cell_z 
    END DO
    CALL inv3(p_cell%amat,p_cell%bmat,p_cell%omtil)
    p_cell%bmat=p_cell%bmat*tpi_const
    p_cell%latnam=cell%latnam

    p_kpts=kpts
    write(1088,*) 'banddos%unfoldband: ', banddos%unfoldband
    write(1088,*) 'brav. matrix: '
    write(1088,'(3f15.8)') cell%amat
    write(1088,*) 'brav. rez. matrix: '
    write(1088,'(3f15.8)') cell%bmat
    write(1088,*) ' primitive brav. matrix: '
    write(1088,'(3f15.8)') p_cell%amat
    write(1088,*)'written as normal matrix (zeile, spalte)'
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%amat(1,1), p_cell%amat(1,2), p_cell%amat(1,3)
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%amat(2,1), p_cell%amat(2,2), p_cell%amat(2,3)
    write(1088,'(f15.8,f15.8,f15.8)') p_cell%amat(3,1), p_cell%amat(3,2), p_cell%amat(3,3)
    write(1088,*) 'primitive brav. rez. matrix: '
    write(89,'(3f15.8)') p_cell%bmat
    write(1088,'(3f15.8)') p_cell%bmat
    write(1088,'(a,i7,a,i7)') 'kpts%nkpt',kpts%nkpt,'   p_kpts%nkpt',p_kpts%nkpt
    write(1088,*) kpts%specialPoints


  END SUBROUTINE unfold_band_kpts
  
  SUBROUTINE find_supercell_kpts(banddos,p_cell,cell,p_kpts,kpts)
    USE m_types
    USE m_juDFT
    implicit none

    TYPE(t_banddos),INTENT(IN)  :: banddos
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_cell),INTENT(IN)     :: p_cell
    TYPE(t_kpts),INTENT(IN)     :: p_kpts
    TYPE(t_kpts),INTENT(INOUT)  :: kpts
    
    INTEGER :: i,m1,m2,m3
    REAL    :: list(9,p_kpts%nkpt)  !cartesion coordinates for k,K,m
    REAL    :: pc_kpoint_i(3)    !primitive cell kpoint internal
    REAL    :: sc_kpoint_i(3)    !super cell kpoint internal
    REAL    :: pc_kpoint_c(3)    !primitive cell kpoint cartesian
    REAL    :: sc_kpoint_c(3)    !super cell kpoint cartesian
    LOGICAL :: representation_found
    write(1088,'(a,i7,a,i7)') 'kpts%nkpt',kpts%nkpt,'   p_kpts%nkpt',p_kpts%nkpt
    write(1088,'(3f15.8)') p_kpts%specialPoints
    write(88,'(3f15.8)') p_kpts%bk

    DO i= 1,size(list,2)
        pc_kpoint_c(1)=p_kpts%bk(1,i)*p_cell%bmat(1,1)+p_kpts%bk(2,i)*p_cell%bmat(1,2)+p_kpts%bk(3,i)*p_cell%bmat(1,3)
        pc_kpoint_c(2)=p_kpts%bk(2,i)*p_cell%bmat(2,1)+p_kpts%bk(2,i)*p_cell%bmat(2,2)+p_kpts%bk(3,i)*p_cell%bmat(2,3)
        pc_kpoint_c(3)=p_kpts%bk(3,i)*p_cell%bmat(3,1)+p_kpts%bk(2,i)*p_cell%bmat(3,2)+p_kpts%bk(3,i)*p_cell%bmat(3,3)
	list(1,i)=pc_kpoint_c(1)
	list(2,i)=pc_kpoint_c(2)
	list(3,i)=pc_kpoint_c(3)
	representation_found=.false.
m_loop:	DO m1= -banddos%s_cell_x,banddos%s_cell_x
		DO m2= -banddos%s_cell_y,banddos%s_cell_y
			DO m3= -banddos%s_cell_z,banddos%s_cell_z
				pc_kpoint_c(1)=list(1,i)-m1*cell%bmat(1,1)-m2*cell%bmat(1,2)-m3*cell%bmat(1,3)
				pc_kpoint_c(2)=list(2,i)-m1*cell%bmat(2,1)-m2*cell%bmat(2,2)-m3*cell%bmat(2,3)
				pc_kpoint_c(3)=list(3,i)-m1*cell%bmat(3,1)-m2*cell%bmat(3,2)-m3*cell%bmat(3,3)
				IF ((dot_product(pc_kpoint_c(:), cell%bmat(:,1)) >= 0).AND.(dot_product(pc_kpoint_c(:), cell%bmat(:,1)) < dot_product(cell%bmat(:,1), cell%bmat(:,1))) &
				     & .AND. (dot_product(pc_kpoint_c(:), cell%bmat(:,2)) >= 0).AND.(dot_product(pc_kpoint_c(:), cell%bmat(:,2)) < dot_product(cell%bmat(:,2), cell%bmat(:,2))) &
				     & .AND. (dot_product(pc_kpoint_c(:), cell%bmat(:,3)) >= 0).AND.(dot_product(pc_kpoint_c(:), cell%bmat(:,3)) < dot_product(cell%bmat(:,3), cell%bmat(:,3)))) THEN
					list(4,i)=pc_kpoint_c(1)
					list(5,i)=pc_kpoint_c(2)
					list(6,i)=pc_kpoint_c(3)
					list(7,i)=-m1
					list(8,i)=-m2
					list(9,i)=-m3
     				     write(*,'(a,f15.8,f15.8,f15.8,6l,3i)') 'representation found for the following kpoint:',list(1,i),list(2,i),list(3,i),(dot_product(pc_kpoint_c(:),&
      				       & cell%bmat(:,1)) >= 0),(dot_product(pc_kpoint_c(:), cell%bmat(:,1)) <= dot_product(cell%bmat(:,1), cell%bmat(:,1)) ) &
				     & , (dot_product(pc_kpoint_c(:), cell%bmat(:,2)) >= 0),(dot_product(pc_kpoint_c(:), cell%bmat(:,2)) <= dot_product(cell%bmat(:,2), cell%bmat(:,2))) &
				     & ,(dot_product(pc_kpoint_c(:), cell%bmat(:,3)) >= 0),(dot_product(pc_kpoint_c(:), cell%bmat(:,3)) <= dot_product(cell%bmat(:,3), cell%bmat(:,3))),-m1,-m2,-m3
					representation_found=.true.
				END IF
       			        IF (representation_found) EXIT m_loop
			END DO
		END DO
	END DO m_loop
        IF (.not.representation_found) THEN
        write(*,'(a,f15.8,f15.8,f15.8)') 'No representation found for the following kpoint:',list(1,i),list(2,i),list(3,i)
        END IF 
    END DO
    write(90,'(9f15.8)') list
    CALL juDFT_error('Not yet implemented', calledby='find_supercell_kpts')
  END SUBROUTINE find_supercell_kpts

END MODULE m_unfold_band_kpts
