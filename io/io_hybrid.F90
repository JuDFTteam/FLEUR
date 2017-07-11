module m_io_hybrid
  use m_io_matrix
  use m_judft
  use m_types
  implicit none
  !private
  integer,save :: id_olap,id_z,id_v_x
  !public:: open_hybrid_io,read_cmt,write_cmt
contains

  subroutine open_hybrid_io(hybrid,dimension,atoms,l_real)
    implicit none
    TYPE(t_hybrid),INTENT(IN)   :: hybrid
    TYPE(t_dimension),INTENT(IN):: dimension
    TYPE(t_atoms),INTENT(IN)    :: atoms
    LOGICAL,INTENT(IN)          :: l_real

    
    OPEN(unit=777,file='cmt',form='unformatted',access='direct',&
         &     recl=dimension%neigd*hybrid%maxlmindx*atoms%nat*16)
    print *,"Open olap.mat"
    id_olap=OPEN_MATRIX(l_real,dimension%nbasfcn,1,"olap.mat")
    print *,"Open z.mat"
    id_z=OPEN_MATRIX(l_real,dimension%nbasfcn,1,"z.mat")
    print *,"Open v_x.mat"
    id_v_x=OPEN_MATRIX(l_real,dimension%nbasfcn,1,"v_x.mat")
    
  end subroutine open_hybrid_io
  
  subroutine write_cmt(cmt,nk)
    implicit none
    complex,INTENT(IN):: cmt(:,:,:)
    integer,INTENT(IN):: nk

    write(777,rec=nk) cmt
  end subroutine write_cmt

  subroutine read_cmt(cmt,nk)
    implicit none
    complex,INTENT(OUT):: cmt(:,:,:)
    integer,INTENT(IN):: nk

    read(777,rec=nk) cmt
  end subroutine read_cmt


  subroutine read_olap(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(INOUT):: mat
    INTEGER,INTENT(IN)           :: rec
    
    CALL read_matrix(mat,rec,id_olap)
  END subroutine read_olap

    subroutine write_olap(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(IN)   :: mat
    INTEGER,INTENT(IN)           :: rec
    
    CALL write_matrix(mat,rec,id_olap)
  END subroutine write_olap

  subroutine read_z(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(INOUT):: mat
    INTEGER,INTENT(IN)           :: rec
    
    CALL read_matrix(mat,rec,id_z)
  END subroutine read_z

    subroutine write_z(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(IN)   :: mat
    INTEGER,INTENT(IN)           :: rec
    
    CALL write_matrix(mat,rec,id_z)
  END subroutine write_z

  subroutine read_v_x(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(INOUT):: mat
    INTEGER,INTENT(IN)           :: rec
    
    CALL read_matrix(mat,rec,id_v_x)
  END subroutine read_v_x

  subroutine write_v_x(mat,rec)
    implicit none
    TYPE(t_mat),INTENT(IN)   :: mat
    INTEGER,INTENT(IN)           :: rec
    
    CALL write_matrix(mat,rec,id_v_x)
  END subroutine write_v_x

  
 

end module m_io_hybrid
