MODULE m_eig66_mpi
#include "juDFT_env.h"
  USE m_eig66_data
  USE m_types
  USE m_judft
#ifdef CPP_MPI
  use mpi
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC open_eig,read_eig,write_eig,close_eig,reset_eig
CONTAINS

  SUBROUTINE priv_find_data(id,d)
    INTEGER,INTENT(IN)::id
    TYPE(t_data_mpi),POINTER,ASYNCHRONOUS:: d

    CLASS(t_data),POINTER   ::dp
    CALL eig66_find_data(dp,id)
    SELECT TYPE(dp)
    TYPE is (t_data_mpi)
       d=>dp
       CLASS default
       CALL judft_error("BUG: wrong datatype in eig66_mpi")
    END SELECT
  END SUBROUTINE priv_find_data


  SUBROUTINE open_eig(id,mpi_comm,nmat,neig,nkpts,jspins,create,l_real,l_soc,l_noco,n_size_opt,filename)
    USE,INTRINSIC::iso_c_binding
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id,mpi_comm,nmat,neig,nkpts,jspins
    LOGICAL, INTENT(IN) :: l_noco,create,l_real,l_soc
    INTEGER,INTENT(IN),OPTIONAL:: n_size_opt
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
#ifdef CPP_MPI
    INTEGER:: isize,e,slot_size,local_slots
    INTEGER,PARAMETER::mcored=27 !there should not be more that 27 core states
    TYPE(t_data_MPI),POINTER,ASYNCHRONOUS :: d

    CALL priv_find_data(id,d)
    CALL eig66_data_storedefault(d,jspins,nkpts,nmat,neig,l_real.and..not.l_soc,l_soc)

    IF (PRESENT(n_size_opt)) d%n_size=n_size_opt
    IF (ALLOCATED(d%pe_ev)) THEN
       IF (create) CALL reset_eig(id,l_soc)
       IF (PRESENT(filename)) CALL judft_error("Storing of data not implemented for MPI case",calledby="eig66_mpi.F")
       RETURN !everything already done!
    ENDIF

    CALL timestart("create data spaces in ei66_mpi")
    CALL MPI_COMM_RANK(MPI_COMM,d%irank,e)
    CALL MPI_COMM_SIZE(MPI_COMM,isize,e)

    CALL create_maps(d,isize,nkpts,jspins,neig,d%n_size)
    local_slots=COUNT(d%pe_basis==d%irank)
    !Now create the windows

    !Window for neig
    slot_size=1
    CALL priv_create_memory(1,local_slots,d%neig_handle,d%neig_data)
    d%neig_data=0

    !The eigenvalues
    d%size_eig=neig
    CALL priv_create_memory(d%size_eig,local_slots,d%eig_handle,real_data_ptr=d%eig_data)
    d%eig_data=1E99
    !The w_iks
    CALL priv_create_memory(d%size_eig,local_slots,d%w_iks_handle,real_data_ptr=d%w_iks_data)
    d%w_iks_data=1E99

    !The eigenvectors
    local_slots=COUNT(d%pe_ev==d%irank)
    slot_size=nmat
    IF (l_real.AND..NOT.l_soc) THEN
       CALL priv_create_memory(slot_size,local_slots,d%zr_handle,real_data_ptr=d%zr_data)
    else
       CALL priv_create_memory(slot_size,local_slots,d%zc_handle,cmplx_data_ptr=d%zc_data)
    ENDIF
    IF (PRESENT(filename).AND..NOT.create) CALL judft_error("Storing of data not implemented for MPI case",calledby="eig66_mpi.F")
    CALL MPI_BARRIER(MPI_COMM,e)
    CALL timestop("create data spaces in ei66_mpi")
  CONTAINS
    SUBROUTINE priv_create_memory(slot_size,local_slots,handle,int_data_ptr,real_data_ptr,cmplx_data_ptr)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: slot_size,local_slots
      INTEGER,POINTER,OPTIONAL,ASYNCHRONOUS  :: int_data_ptr(:)
      REAL   ,POINTER,OPTIONAL,ASYNCHRONOUS  :: real_data_ptr(:)
      COMPLEX,POINTER,OPTIONAL,ASYNCHRONOUS  :: cmplx_data_ptr(:)
      INTEGER,INTENT(OUT)          :: handle
#ifdef CPP_MPI
      TYPE(c_ptr)::ptr
      INTEGER:: e
      INTEGER(MPI_ADDRESS_KIND) :: length
      INTEGER                   :: type_size

      length=0   
      IF (present(real_data_ptr)) THEN
          length=length+1
          CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,type_size,e)
      ENDIF
      IF (present(cmplx_data_ptr)) THEN
          length=length+1
          CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,type_size,e)
      ENDIF
      IF (present(int_data_ptr)) THEN 
          length=length+1
          CALL MPI_TYPE_SIZE(MPI_INTEGER,type_size,e)
      ENDIF
      if (length.ne.1) call judft_error("Bug in eig66_mpi:create_memory") 
      length=MAX(1,slot_size*local_slots)
 
#ifdef CPP_MPI_ALLOC      
      length=length*type_size
      CALL MPI_ALLOC_MEM(length,MPI_INFO_NULL,ptr,e)
      IF (e.NE.0) CPP_error("Could not allocated MPI-Data in eig66_mpi")
#endif	
      IF (PRESENT(real_data_ptr)) THEN
#ifdef CPP_MPI_ALLOC         
         CALL C_F_POINTER(ptr,real_data_ptr,(/length/type_size/))
#else
         ALLOCATE(real_data_ptr(length))
#endif         
      	CALL MPI_WIN_CREATE(real_data_ptr, length,slot_size*type_size,Mpi_INFO_NULL, MPI_COMM,handle, e)
    ELSEIF(PRESENT(int_data_ptr)) THEN
#ifdef CPP_MPI_ALLOC
       CALL C_F_POINTER(ptr,int_data_ptr,(/length/type_size/))
#else
       ALLOCATE(int_data_ptr(length))
#endif         
      	CALL MPI_WIN_CREATE(int_data_ptr, length,slot_size*type_size,Mpi_INFO_NULL, MPI_COMM,handle, e)
    ELSE
#ifdef CPP_MPI_ALLOC       
       CALL C_F_POINTER(ptr,cmplx_data_ptr,(/length/type_size/))
#else
       ALLOCATE(cmplx_data_ptr(length))
#endif   
       CALL MPI_WIN_CREATE(cmplx_data_ptr, length,slot_size*type_size,Mpi_INFO_NULL, MPI_COMM,handle, e)
    ENDIF
#endif
    END SUBROUTINE priv_create_memory


#endif
  END SUBROUTINE open_eig
  SUBROUTINE close_eig(id,delete,filename)
    INTEGER,INTENT(IN)         :: id
    LOGICAL,INTENT(IN),OPTIONAL:: delete
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL::filename
    TYPE(t_data_MPI),POINTER,ASYNCHRONOUS :: d
    CALL priv_find_data(id,d)

    IF (PRESENT(delete)) THEN
       IF (delete) WRITE(*,*) "No deallocation of memory implemented in eig66_mpi"
    ENDIF
    IF (PRESENT(filename)) CALL judft_error("Storing of data not implemented for MPI case",calledby="eig66_mpi.F")
  END SUBROUTINE close_eig

  SUBROUTINE read_eig(id,nk,jspin,neig,eig,w_iks,list,zmat)
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: id,nk,jspin
    INTEGER, INTENT(OUT),OPTIONAL  :: neig
    REAL,    INTENT(OUT),OPTIONAL  :: eig(:),w_iks(:)
    INTEGER, INTENT(IN),OPTIONAL   :: list(:)
    TYPE(t_mat),OPTIONAL  :: zmat

#ifdef CPP_MPI
    INTEGER                   :: pe,tmp_size,e
    INTEGER(MPI_ADDRESS_KIND) :: slot
    INTEGER                   :: n1,n2,n3,n
    INTEGER,ALLOCATABLE,ASYNCHRONOUS       :: tmp_int(:)
    REAL,ALLOCATABLE,ASYNCHRONOUS          :: tmp_real(:)
    COMPLEX,ALLOCATABLE,ASYNCHRONOUS       :: tmp_cmplx(:)
    TYPE(t_data_MPI),POINTER,ASYNCHRONOUS :: d
    CALL priv_find_data(id,d)
    pe=d%pe_basis(nk,jspin)
    slot=d%slot_basis(nk,jspin)
    IF (PRESENT(neig))THEN
       CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%neig_handle,e)
       ! Get current values
       CALL  MPI_GET(neig,1,MPI_INTEGER,pe,slot,1,MPI_INTEGER,d%neig_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%neig_handle,e)

    ENDIF
    IF (PRESENT(eig).or.PRESENT(w_iks)) THEN
       ALLOCATE(tmp_real(MIN(SIZE(eig),d%size_eig)))
       IF (PRESENT(eig)) THEN
          CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%eig_handle,e)
          CALL MPI_GET(tmp_real,SIZE(tmp_real),MPI_DOUBLE_PRECISION,pe,slot,size(tmp_real),MPI_DOUBLE_PRECISION,d%eig_handle,e)
          CALL MPI_WIN_UNLOCK(pe,d%eig_handle,e)
          eig(:size(tmp_real))=tmp_real
       END IF
       IF (PRESENT(w_iks)) THEN
          CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%w_iks_handle,e)
          CALL MPI_GET(tmp_real,size(tmp_real),MPI_DOUBLE_PRECISION,pe,slot,size(tmp_real),MPI_DOUBLE_PRECISION,d%w_iks_handle,e)
          CALL MPI_WIN_UNLOCK(pe,d%w_iks_handle,e)
          w_iks(:SIZE(tmp_real))=tmp_real
       END IF
       DEALLOCATE(tmp_real)
    ENDIF

    IF (PRESENT(zmat)) THEN
       tmp_size=zmat%matsize1
       ALLOCATE(tmp_real(tmp_size))
       ALLOCATE(tmp_cmplx(tmp_size))
       DO n=1,zmat%matsize2
          IF (PRESENT(list)) THEN
             IF (n>SIZE(list)) CYCLE
             n1=list(n)
          END IF
          slot=d%slot_ev(nk,jspin,n1)
          pe=d%pe_ev(nk,jspin,n1)
          
          if (zmat%l_real) THEN
             if (.not.d%l_real) THEN
                CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%zc_handle,e)
                CALL MPI_GET(tmp_cmplx,tmp_size,MPI_DOUBLE_COMPLEX,pe,slot,tmp_size,MPI_DOUBLE_COMPLEX,d%zc_handle,e)
                CALL MPI_WIN_UNLOCK(pe,d%zc_handle,e)
                !print *, nk,jspin,n1,"r PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_cmplx(1)
                zmat%data_r(:,n)=REAL(tmp_cmplx)
             else
                CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%zr_handle,e)
                CALL MPI_GET(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,d%zr_handle,e)
                CALL MPI_WIN_UNLOCK(pe,d%zr_handle,e)
                !print *, nk,jspin,n1,"r PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_real(1)
                zmat%data_r(:,n)=tmp_real
             endif
          ELSE
             if (d%l_real) call judft_error("Could not read complex data, only real data is stored",calledby="eig66_mpi%read_eig")
             CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%zc_handle,e)
             CALL MPI_GET(tmp_cmplx,tmp_size,MPI_DOUBLE_COMPLEX,pe,slot,tmp_size,MPI_DOUBLE_COMPLEX,d%zc_handle,e)
             CALL MPI_WIN_UNLOCK(pe,d%zc_handle,e)
             !print *, nk,jspin,n1,"r PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_cmplx(1)
             zmat%data_c(:,n)=tmp_cmplx
          ENDIF
       ENDDO
    ENDIF

#endif
  END SUBROUTINE read_eig

  SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,eig,w_iks,n_size,n_rank,zmat)
    INTEGER, INTENT(IN)          :: id,nk,jspin
    INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
    INTEGER, INTENT(IN),OPTIONAL :: neig,neig_total
    REAL,    INTENT(IN),OPTIONAL :: eig(:),w_iks(:)
    TYPE(t_mat),INTENT(IN),OPTIONAL :: zmat

#ifdef CPP_MPI
    INTEGER                   :: pe,tmp_size,e
    INTEGER(MPI_ADDRESS_KIND) :: slot
    INTEGER                   :: n1,n2,n3,n,nn
    INTEGER,ALLOCATABLE,ASYNCHRONOUS       :: tmp_int(:)
    REAL,ALLOCATABLE,ASYNCHRONOUS          :: tmp_real(:)
    COMPLEX,ALLOCATABLE,ASYNCHRONOUS       :: tmp_cmplx(:)
    LOGICAL                   :: acc
    TYPE(t_data_MPI),POINTER,ASYNCHRONOUS :: d

    INTEGER:: irank,ierr

    CALL priv_find_data(id,d)

    CALL MPI_COMM_RANK(MPI_COMM_WORLD,irank,ierr)

    pe=d%pe_basis(nk,jspin)
    slot=d%slot_basis(nk,jspin)
    !write the number of eigenvalues 
    !only one process needs to do it
    IF (PRESENT(neig_total)) THEN
       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%neig_handle,e)
       ALLOCATE(tmp_int(1))
       tmp_int(1)=neig_total
       CALL MPI_PUT(tmp_int,1,MPI_INTEGER,pe,slot,1,MPI_INTEGER,d%neig_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%neig_handle,e)
       DEALLOCATE(tmp_int)
    ENDIF

    !write the eigenvalues 
    !only one process needs to do it
    IF (PRESENT(eig).OR.PRESENT(w_iks)) THEN
       ALLOCATE(tmp_real(d%size_eig))
       tmp_real=1E99
       if (PRESENT(EIG)) THEN
          tmp_real(:SIZE(eig)) = eig(:SIZE(eig))
          CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%eig_handle,e)
          CALL MPI_PUT(tmp_real,d%size_eig,MPI_DOUBLE_PRECISION,pe,slot,d%size_eig,MPI_DOUBLE_PRECISION,d%eig_handle,e)
          CALL MPI_WIN_UNLOCK(pe,d%eig_handle,e)
       END if
       IF (PRESENT(w_iks)) THEN
          tmp_real(:size(w_iks))=w_iks
          CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%w_iks_handle,e)
          CALL MPI_PUT(tmp_real,d%size_eig,MPI_DOUBLE_PRECISION,pe,slot,d%size_eig,MPI_DOUBLE_PRECISION,d%w_iks_handle,e)
          CALL MPI_WIN_UNLOCK(pe,d%w_iks_handle,e)
       END IF
       DEALLOCATE(tmp_real)
    ENDIF

    !write the eigenvectors
    !all procceses participate 
    IF (PRESENT(zmat)) THEN
       tmp_size=zmat%matsize1
       ALLOCATE(tmp_real(tmp_size))
       ALLOCATE(tmp_cmplx(tmp_size))
       DO n=1,zmat%matsize2
          n1=n-1
          IF (PRESENT(n_size)) n1=n_size*n1
          IF (PRESENT(n_rank)) n1=n1+n_rank
          IF (n1+1>size(d%slot_ev,3)) EXIT
          slot=d%slot_ev(nk,jspin,n1+1)
          pe=d%pe_ev(nk,jspin,n1+1)
          IF (zmat%l_real) THEN
             if (.not.d%l_real) THEN
                tmp_cmplx=zmat%data_r(:,n)
                CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%zc_handle,e)
                CALL MPI_PUT(tmp_cmplx,tmp_size,MPI_DOUBLE_COMPLEX,pe,slot,tmp_size,MPI_DOUBLE_COMPLEX,d%zc_handle,e)
                CALL MPI_WIN_UNLOCK(pe,d%zc_handle,e)
               else
                tmp_real=zmat%data_r(:,n)
                CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%zr_handle,e)
                CALL MPI_PUT(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,d%zr_handle,e)
                CALL MPI_WIN_UNLOCK(pe,d%zr_handle,e)
             endif
          ELSE
             if (d%l_real) CALL juDFT_error("Could not write complex data to file prepared for real data",calledby="eig66_mpi%write_eig")
             tmp_cmplx=zmat%data_c(:,n)
             CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%zc_handle,e)
             CALL MPI_PUT(tmp_cmplx,tmp_size,MPI_DOUBLE_COMPLEX,pe,slot,tmp_size,MPI_DOUBLE_COMPLEX,d%zc_handle,e)
             CALL MPI_WIN_UNLOCK(pe,d%zc_handle,e)
          ENDIF
       ENDDO
    ENDIF

#endif
  END SUBROUTINE write_eig

  SUBROUTINE reset_eig(id,l_soc)
    INTEGER, INTENT(IN)        :: id
    LOGICAL, INTENT(IN)        :: l_soc
#ifdef CPP_MPI
    TYPE(t_data_MPI),POINTER,ASYNCHRONOUS :: d
    CALL priv_find_data(id,d)

    d%neig_data=0
    d%eig_data=1E99
    d%w_iks_data=1E99
    if (d%l_real.and..not.l_soc) THEN
       d%zr_data=0.0
    else
       d%zc_data=0.0
    endif
#endif
  END SUBROUTINE reset_eig


#ifdef CPP_MPI
  SUBROUTINE create_maps(d,isize,nkpts,jspins,neig,n_size)
    IMPLICIT NONE
    TYPE(t_data_MPI),INTENT(INOUT),ASYNCHRONOUS:: d
    INTEGER,INTENT(IN):: isize,nkpts,jspins,neig,n_size

    INTEGER:: nk,j,n1,n2,n,pe,n_members
    INTEGER::used(0:isize)

    ALLOCATE(d%pe_basis(nkpts,jspins),d%slot_basis(nkpts,jspins))
    ALLOCATE(d%pe_ev(nkpts,jspins,neig),d%slot_ev(nkpts,jspins,neig))

    !basis contains a total of nkpts*jspins entries
    d%pe_basis=-1
    d%pe_ev=-1
    used=0
    n_members=isize/n_size !no of k-points in parallel
    DO j=1,jspins
       DO nk=1,nkpts
          n1=nk+(j-1)*nkpts-1
          pe=MOD(n1,n_members)*n_size
          d%pe_basis(nk,j)=pe
          d%slot_basis(nk,j)=used(pe)
          used(pe)=used(pe)+1
       ENDDO
    ENDDO
    used=0
    DO n=1,neig
       DO j=1,jspins
          DO nk=1,nkpts
             n1=nk+(j-1)*nkpts-1
             !eigenvectors have more entries
             !pe=MOD(n1,n_members)*n_size+MOD(n,n_size)
             pe=MOD(n1,n_members)*n_size+MOD(n-1,n_size)
             d%pe_ev(nk,j,n)=pe
             d%slot_ev(nk,j,n)=used(pe)
             used(pe)=used(pe)+1
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE create_maps
#endif

END MODULE m_eig66_mpi
