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
    write(1088,*) 'primitive brav. rez. matrix: '
    write(1088,'(3f15.8)') p_cell%bmat
    write(1088,'(a,i7,a,i7)') 'kpts%nkpt',kpts%nkpt,'   p_kpts%nkpt',p_kpts%nkpt


  END SUBROUTINE unfold_band_kpts
  
  SUBROUTINE find_supercell_kpts(banddos,p_cell,cell,p_kpts,kpts)
    USE m_types
    USE m_juDFT
    implicit none

    TYPE(t_banddos),INTENT(IN)  :: banddos
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_cell),INTENT(IN)  :: p_cell
    TYPE(t_kpts),INTENT(IN)  :: p_kpts
    TYPE(t_kpts),INTENT(INOUT)  :: kpts

    CALL juDFT_error('Not yet implemented', calledby='find_supercell_kpts')

  END SUBROUTINE find_supercell_kpts

END MODULE m_unfold_band_kpts
