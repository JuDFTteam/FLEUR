MODULE m_eig66_mpi
#include "juDFT_env.h"
  USE m_eig66_data
#ifdef CPP_MPI
  use mpi
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC open_eig,read_eig,write_eig,close_eig,write_dos,read_dos
CONTAINS

  SUBROUTINE priv_find_data(id,d)
    INTEGER,INTENT(IN)::id
    TYPE(t_data_mpi),POINTER:: d

    CLASS(t_data),POINTER   ::dp
    CALL eig66_find_data(dp,id)
    SELECT TYPE(dp)
    TYPE is (t_data_mpi)
       d=>dp
       CLASS default
       CALL judft_error("BUG: wrong datatype in eig66_mpi")
    END SELECT
  END SUBROUTINE priv_find_data


  SUBROUTINE open_eig(id,mpi_comm,nmat,neig,nkpts,jspins,lmax,nlo,ntype,create,nlotot,l_noco,n_size_opt,l_dos,l_mcd,l_orb,filename,layers,nstars,ncored,nsld,nat)
    USE,INTRINSIC::iso_c_binding
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: id,mpi_comm,nmat,neig,nkpts,jspins,nlo,ntype,lmax,nlotot
    LOGICAL, INTENT(IN) :: l_noco,create
    INTEGER,INTENT(IN),OPTIONAL:: n_size_opt
    LOGICAL,INTENT(IN),OPTIONAL ::l_dos,l_mcd,l_orb
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
    INTEGER,INTENT(IN),OPTIONAL :: layers,nstars,ncored,nsld,nat
#ifdef CPP_MPI
    INTEGER:: isize,e,slot_size,local_slots
    INTEGER,PARAMETER::mcored=27 !there should not be more that 27 core states
    TYPE(t_data_MPI),POINTER :: d

    CALL priv_find_data(id,d)
    CALL eig66_data_storedefault(d,jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype,l_dos,l_mcd,l_orb)

    IF (PRESENT(n_size_opt)) d%n_size=n_size_opt
    IF (ALLOCATED(d%pe_ev)) THEN
       IF (create) THEN
          d%neig_data=0
          d%eig_data=1E99
          d%int_data=9999999
          d%real_data=1E99
#if defined(CPP_INVERSION)&&!defined(CPP_SOC)
     d%zr_data=0.0
#else
          d%zc_data=0.0
#endif
          d%qal_data=0.0
          d%qvac_data=0.0
          d%qvlay_data=0.0
          d%qstars_data=0.0
          d%ksym_data=0.0
          d%jsym_data=0.0
          d%mcd_data=0.0
          d%qintsl_data=0.0
          d%qmtsl_data=0.0
          d%qmtp_data=0.0
          d%orbcomp_data=0.0
       ENDIF
       IF (PRESENT(filename)) CALL priv_readfromfile()
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

    !The integer values
    d%size_k=nmat
    slot_size=(5+3*d%size_k+1+nlotot)
    CALL priv_create_memory(slot_size,local_slots,d%int_handle,d%int_data)
    d%int_data=9999999

    !The real values
    d%size_el=(1+lmax)*ntype
    d%size_ello=nlo*ntype
    slot_size=(6+d%size_el+d%size_ello)
    CALL priv_create_memory(slot_size,local_slots,d%real_handle,real_data_ptr=d%real_data)
    d%real_data=1E99

    !The eigenvalues
    d%size_eig=neig
    CALL priv_create_memory(d%size_eig,local_slots,d%eig_handle,real_data_ptr=d%eig_data)
    d%eig_data=1E99
    !The eigenvectors
    local_slots=COUNT(d%pe_ev==d%irank)
    slot_size=nmat

#if defined(CPP_INVERSION)&&!defined(CPP_SOC)
    CALL priv_create_memory(slot_size,local_slots,d%zr_handle,real_data_ptr=d%zr_data)
#else
    CALL priv_create_memory(slot_size,local_slots,d%zc_handle,cmplx_data_ptr=d%zc_data)
#endif
    !Data for DOS etc
    IF (d%l_dos) THEN
       local_slots=COUNT(d%pe_basis==d%irank)
       CALL priv_create_memory(4*ntype*neig,local_slots,d%qal_handle,real_data_ptr=d%qal_data)
       CALL priv_create_memory(neig*2,local_slots,d%qvac_handle,real_data_ptr=d%qvac_data)
       CALL priv_create_memory(neig,local_slots,d%qis_handle,real_data_ptr=d%qis_data)
       CALL priv_create_memory(neig*max(1,layers)*2,local_slots,d%qvlay_handle,real_data_ptr=d%qvlay_data)
       CALL priv_create_memory(max(1,nstars)*neig*max(1,layers)*2,local_slots,d%qstars_handle,cmplx_data_ptr=d%qstars_data)
       CALL priv_create_memory(neig,local_slots,d%jsym_handle,d%jsym_data)
       CALL priv_create_memory(neig,local_slots,d%ksym_handle,d%ksym_data)
       IF (l_mcd) CALL priv_create_memory(3*ntype*mcored*neig,local_slots,d%mcd_handle,real_data_ptr=d%mcd_data)
       IF (l_orb) THEN
          CALL priv_create_memory(nsld*neig,local_slots,d%qintsl_handle,real_data_ptr=d%qintsl_data)
          CALL priv_create_memory(nsld*neig,local_slots,d%qmtsl_handle,real_data_ptr=d%qmtsl_data)
          CALL priv_create_memory(nat*neig,local_slots,d%qmtp_handle,real_data_ptr=d%qmtp_data)
          CALL priv_create_memory(23*nat*neig,local_slots,d%orbcomp_handle,real_data_ptr=d%orbcomp_data)
       ENDIF
    ELSE
       ALLOCATE(d%qal_data(1),d%qvac_data(1),d%qis_data(1),d%qvlay_data(1),d%qstars_data(1),&
            d%jsym_data(1),d%ksym_data(1),d%mcd_data(1),d%qintsl_data(1),d%qmtsl_data(1),&
            d%qmtp_data(1),d%orbcomp_data(1))
    ENDIF
    IF (PRESENT(filename).AND..NOT.create) CALL priv_readfromfile()
    CALL timestop("create data spaces in ei66_mpi")
  CONTAINS
    SUBROUTINE priv_create_memory(slot_size,local_slots,handle,int_data_ptr,real_data_ptr,cmplx_data_ptr)
      IMPLICIT NONE
      INTEGER,INTENT(IN)           :: slot_size,local_slots
      INTEGER,POINTER,INTENT(OUT),OPTIONAL  :: int_data_ptr(:)
      REAL   ,POINTER,INTENT(OUT),OPTIONAL  :: real_data_ptr(:)
      COMPLEX,POINTER,INTENT(OUT),OPTIONAL  :: cmplx_data_ptr(:)
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
      length=slot_size*local_slots
 
      length=length*type_size

      CALL MPI_ALLOC_MEM(length,MPI_INFO_NULL,ptr,e)
      IF (e.NE.0) CPP_error("Could not allocated MPI-Data in eig66_mpi")
	
      IF (present(real_data_ptr)) THEN	
      	CALL C_F_POINTER(ptr,real_data_ptr,(/length/type_size/))
      	CALL MPI_WIN_CREATE(real_data_ptr, length,slot_size*type_size,Mpi_INFO_NULL, MPI_COMM,handle, e)
      ELSEIF(present(int_data_ptr)) THEN
      	CALL C_F_POINTER(ptr,int_data_ptr,(/length/type_size/))
      	CALL MPI_WIN_CREATE(int_data_ptr, length,slot_size*type_size,Mpi_INFO_NULL, MPI_COMM,handle, e)
      ELSE
      	CALL C_F_POINTER(ptr,cmplx_data_ptr,(/length/type_size/))
      	CALL MPI_WIN_CREATE(cmplx_data_ptr, length,slot_size*type_size,Mpi_INFO_NULL, MPI_COMM,handle, e)
      ENDIF
#endif
    END SUBROUTINE priv_create_memory


    SUBROUTINE priv_readfromfile()
      USE m_eig66_DA,ONLY:open_eig_DA=>open_eig,read_eig_DA=>read_eig,close_eig_da=>close_eig
      INTEGER:: jspin,nk,i,ii,iii,nv,tmp_id
      REAL   :: wk,bk3(3),evac(2)
      INTEGER :: k1(nmat),k2(nmat),k3(nmat),kveclo(nlotot)
      REAL    :: eig(neig),ello(nlo,ntype),el(lmax,ntype)
#ifdef CPP_INVERSION
      REAL    :: z(nmat,neig)
#else
      COMPLEX :: z(nmat,neig)
#endif
      !only do this with PE=0
      IF (d%irank==0) THEN
         tmp_id=eig66_data_newid(DA_mode)
         IF (d%l_dos) CPP_error("Could not read DOS data")
         CALL open_eig_DA(tmp_id,nmat,neig,nkpts,jspins,lmax,nlo,ntype,nlotot,.FALSE.,.FALSE.,.FALSE.,.FALSE.,filename)
         DO jspin=1,jspins
            DO nk=1,nkpts
               CALL read_eig_DA(tmp_id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z)
               CALL write_eig(id,nk,jspin,ii,ii,nv,nmat,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z)
            ENDDO
         ENDDO
         CALL close_eig_DA(tmp_id)
      ENDIF
    END SUBROUTINE priv_readfromfile
#endif
  END SUBROUTINE open_eig
  SUBROUTINE close_eig(id,delete,filename)
    INTEGER,INTENT(IN)         :: id
    LOGICAL,INTENT(IN),OPTIONAL:: delete
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL::filename
    TYPE(t_data_MPI),POINTER :: d
    CALL priv_find_data(id,d)

    IF (PRESENT(delete)) THEN
       IF (delete) WRITE(*,*) "No deallocation of memory implemented in eig66_mpi"
    ENDIF
    IF (PRESENT(filename)) CALL priv_writetofile()
  CONTAINS
    SUBROUTINE priv_writetofile()
      USE m_eig66_DA,ONLY:open_eig_DA=>open_eig,write_eig_DA=>write_eig,close_eig_DA=>close_eig
      IMPLICIT NONE

      INTEGER:: nlotot,nk,jspin,nv,i,ii,tmp_id
      REAL   :: wk,bk3(3),evac(2)
      INTEGER :: k1(d%nmat),k2(d%nmat),k3(d%nmat),kveclo(d%nlotot)
      REAL    :: eig(d%neig),ello(d%nlo,d%ntype),el(d%lmax,d%ntype)
#ifdef CPP_INVERSION
      REAL    :: z(d%nmat,d%neig)
#else
      COMPLEX :: z(d%nmat,d%neig)
#endif
      nlotot=d%nlotot

      IF (d%irank==0) THEN
         tmp_id=eig66_data_newid(DA_mode)
         IF (d%l_dos) CPP_error("Could not write DOS data")
         CALL open_eig_DA(tmp_id,d%nmat,d%neig,d%nkpts,d%jspins,d%lmax,d%nlo,d%ntype,d%nlotot,.FALSE.,.FALSE.,.FALSE.,.FALSE.,filename)
         DO jspin=1,d%jspins
            DO nk=1,d%nkpts
               CALL read_eig(id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z)
               CALL write_eig_DA(tmp_id,nk,jspin,ii,ii,nv,i,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z)

            ENDDO
         ENDDO
         CALL close_eig_DA(tmp_id)
      ENDIF
      CALL eig66_remove_data(id)
    END SUBROUTINE priv_writetofile

  END SUBROUTINE close_eig

  SUBROUTINE read_eig(id,nk,jspin,nv,nmat,k1,k2,k3,bk3,wk,neig,eig,el,&
       ello,evac,kveclo,n_start,n_end,z)
    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: id,nk,jspin
    INTEGER, INTENT(OUT),OPTIONAL  :: nv,nmat
    INTEGER, INTENT(OUT),OPTIONAL  :: neig
    REAL,    INTENT(OUT),OPTIONAL  :: eig(:)
    INTEGER, INTENT(OUT),OPTIONAL  :: k1(:),k2(:),k3(:),kveclo(:)
    REAL,    INTENT(OUT),OPTIONAL  :: evac(:),ello(:,:),el(:,:)
    REAL,    INTENT(OUT),OPTIONAL  :: bk3(:),wk
    INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
    CLASS(*),TARGET,INTENT(OUT),OPTIONAL  :: z(:,:)

#ifdef CPP_MPI
    INTEGER                   :: pe,tmp_size,e
    INTEGER(MPI_ADDRESS_KIND) :: slot
    INTEGER                   :: n1,n2,n3,n
    INTEGER,ALLOCATABLE       :: tmp_int(:)
    REAL,ALLOCATABLE          :: tmp_real(:)
    COMPLEX,ALLOCATABLE       :: tmp_cmplx(:)
    TYPE(t_data_MPI),POINTER :: d
    CALL priv_find_data(id,d)
    pe=d%pe_basis(nk,jspin)
    slot=d%slot_basis(nk,jspin)
    IF (PRESENT(neig))THEN
       CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%neig_handle,e)
       ! Get current values
       CALL  MPI_GET(neig,1,MPI_INTEGER,pe,slot,1,MPI_INTEGER,d%neig_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%neig_handle,e)

    ENDIF
    !read the integer values
    IF (ANY((/PRESENT(nv),PRESENT(nmat),PRESENT(k1),PRESENT(k2),PRESENT(k3),PRESENT(kveclo)/))) THEN
       tmp_size=4+3*d%size_k
       IF (PRESENT(kveclo)) tmp_size=tmp_size+SIZE(kveclo)
       ALLOCATE(tmp_int(tmp_size))
       CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%int_handle,e)
       ! Get current values
       CALL  MPI_GET(tmp_int,tmp_size,MPI_INTEGER,pe,slot,tmp_size,MPI_INTEGER,d%int_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%int_handle,e)
       !IF (present(neig)) neig=tmp_int(1)
       IF (PRESENT(nv))   nv=tmp_int(2)
       IF (PRESENT(nmat)) nmat=tmp_int(3)
       IF (PRESENT(k1))   k1=tmp_int(4+1:4+SIZE(k1))
       IF (PRESENT(k2))   k2=tmp_int(4+d%size_k+1:4+d%size_k+SIZE(k2))
       IF (PRESENT(k3))   k3=tmp_int(4+2*d%size_k+1:4+2*d%size_k+SIZE(k3))
       IF (PRESENT(kveclo)) kveclo=tmp_int(4+3*d%size_k+1:4+3*d%size_k+SIZE(kveclo))

    ENDIF
    !read the real-values
    IF (ANY((/PRESENT(wk),PRESENT(bk3),PRESENT(el),PRESENT(ello),PRESENT(evac)/))) THEN
       tmp_size=6+d%size_el+d%size_ello
       ALLOCATE(tmp_real(tmp_size))
       CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%real_handle,e)
       ! Get current values
       CALL  MPI_GET(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,d%real_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%real_handle,e)
       IF (PRESENT(wk))   wk=tmp_real(1)
       IF (PRESENT(bk3))  bk3=tmp_real(2:4)
       IF (PRESENT(evac)) evac=tmp_real(5:6)
       IF (PRESENT(el))   el=RESHAPE(tmp_real(6+1:6+SIZE(el)),SHAPE(el))
       IF (PRESENT(ello)) ello=RESHAPE(tmp_real(6+d%size_el+1:6+d%size_el+SIZE(ello)),SHAPE(ello))
       DEALLOCATE(tmp_real)
    ENDIF
    IF (PRESENT(eig)) THEN
       CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%eig_handle,e)
       ALLOCATE(tmp_real(d%size_eig))
       CALL MPI_GET(tmp_real,d%size_eig,MPI_DOUBLE_PRECISION,pe,slot,d%size_eig,MPI_DOUBLE_PRECISION,d%eig_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%eig_handle,e)
       n1=1;n3=1;n2=SIZE(eig)
       IF (PRESENT(n_start)) n1=n_start
       IF (PRESENT(n_end)) n2=n_end
       eig(:n2-n1+1)=tmp_real(n1:n2)
       DEALLOCATE(tmp_real)
    ENDIF

    IF (PRESENT(z)) THEN
       tmp_size=SIZE(z,1)
       ALLOCATE(tmp_real(tmp_size))
       ALLOCATE(tmp_cmplx(tmp_size))
       DO n=1,SIZE(z,2)
          n1=n
          IF (PRESENT(n_start)) n1=n_start+n-1
          IF (PRESENT(n_end)) THEN
             IF (n1>n_end) CYCLE
          ENDIF
          slot=d%slot_ev(nk,jspin,n1)
          pe=d%pe_ev(nk,jspin,n1)
          SELECT TYPE(z)
          TYPE IS(REAL)
#ifdef CPP_SOC
             CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%zc_handle,e)
             CALL MPI_GET(tmp_cmplx,tmp_size,MPI_DOUBLE_COMPLEX,pe,slot,tmp_size,MPI_DOUBLE_COMPLEX,d%zc_handle,e)
             CALL MPI_WIN_UNLOCK(pe,d%zc_handle,e)
             !print *, nk,jspin,n1,"r PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_cmplx(1)
             z(:,n)=REAL(tmp_cmplx)
#else
             CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%zr_handle,e)
             CALL MPI_GET(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,d%zr_handle,e)
             CALL MPI_WIN_UNLOCK(pe,d%zr_handle,e)
             !print *, nk,jspin,n1,"r PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_real(1)
             z(:,n)=tmp_real
#endif
          TYPE IS (COMPLEX)
             CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%zc_handle,e)
             CALL MPI_GET(tmp_cmplx,tmp_size,MPI_DOUBLE_COMPLEX,pe,slot,tmp_size,MPI_DOUBLE_COMPLEX,d%zc_handle,e)
             CALL MPI_WIN_UNLOCK(pe,d%zc_handle,e)
             !print *, nk,jspin,n1,"r PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_cmplx(1)
             z(:,n)=tmp_cmplx
          END SELECT
       ENDDO
    ENDIF

#endif
  END SUBROUTINE read_eig

  SUBROUTINE write_eig(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk3,wk, &
       eig,el,ello,evac,                     &
       nlotot,kveclo,n_size,n_rank,z)
    INTEGER, INTENT(IN)          :: id,nk,jspin
    INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
    REAL,    INTENT(IN),OPTIONAL :: wk
    INTEGER, INTENT(IN),OPTIONAL :: neig,nv,nmat,nlotot,neig_total
    INTEGER, INTENT(IN),OPTIONAL :: k1(:),k2(:),k3(:),kveclo(:)
    REAL,    INTENT(IN),OPTIONAL :: bk3(3),eig(:),el(:,:)
    REAL,    INTENT(IN),OPTIONAL :: evac(:),ello(:,:)
    CLASS(*),INTENT(IN),OPTIONAL :: z(:,:)

#ifdef CPP_MPI
    INTEGER                   :: pe,tmp_size,e
    INTEGER(MPI_ADDRESS_KIND) :: slot
    INTEGER                   :: n1,n2,n3,n
    INTEGER,ALLOCATABLE       :: tmp_int(:)
    REAL,ALLOCATABLE          :: tmp_real(:)
    COMPLEX,ALLOCATABLE       :: tmp_cmplx(:)
    LOGICAL                   :: acc
    TYPE(t_data_MPI),POINTER :: d

    CALL priv_find_data(id,d)

    pe=d%pe_basis(nk,jspin)
    slot=d%slot_basis(nk,jspin)
    !write the number of eigenvalues values
    IF (PRESENT(neig_total)) THEN
       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%neig_handle,e)
       ALLOCATE(tmp_int(1))
       tmp_int(1)=neig_total
       CALL MPI_PUT(tmp_int,1,MPI_INTEGER,pe,slot,1,MPI_INTEGER,d%neig_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%neig_handle,e)
       DEALLOCATE(tmp_int)
    ENDIF

    IF (ANY((/PRESENT(nv),PRESENT(nmat),PRESENT(nlotot),PRESENT(k1),PRESENT(k2),PRESENT(k3),PRESENT(kveclo)/))) THEN
       tmp_size=5+3*d%size_k
       IF (PRESENT(kveclo)) tmp_size=tmp_size+SIZE(kveclo)
       ALLOCATE(tmp_int(tmp_size))
       tmp_int=9999999
       tmp_int(1)=0
       IF (PRESENT(nv))   tmp_int(2)=nv
       IF (PRESENT(nmat)) tmp_int(3)=nmat
       IF (PRESENT(nlotot)) tmp_int(4)=nlotot
       IF (PRESENT(k1))   tmp_int(4+1:4+SIZE(k1))=k1
       IF (PRESENT(k2))   tmp_int(4+d%size_k+1:4+d%size_k+SIZE(k2))=k2
       IF (PRESENT(k3))   tmp_int(4+2*d%size_k+1:4+2*d%size_k+SIZE(k3))=k3
       IF (PRESENT(kveclo)) tmp_int(4+3*d%size_k+1:4+3*d%size_k+SIZE(kveclo))=kveclo
       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%int_handle,e)
       CALL MPI_ACCUMULATE(tmp_int,tmp_size,MPI_INTEGER,pe,slot,tmp_size,MPI_INTEGER,MPI_MIN,d%int_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%int_handle,e)
    ENDIF
    !write the real-values
    IF (ANY((/PRESENT(wk),PRESENT(bk3),PRESENT(el),PRESENT(ello),PRESENT(evac)/))) THEN
       tmp_size=6+d%size_el+d%size_ello
       ALLOCATE(tmp_real(tmp_size))
       tmp_real=1E99
       IF (PRESENT(wk))   tmp_real(1)=wk
       IF (PRESENT(bk3))  tmp_real(2:4)=bk3
       IF (PRESENT(evac)) tmp_real(5:6)=evac
       IF (PRESENT(el))   tmp_real(6+1:6+SIZE(el))=RESHAPE(el,(/SIZE(el)/))
       IF (PRESENT(ello)) tmp_real(6+d%size_el+1:6+d%size_el+SIZE(ello))=RESHAPE(ello,(/SIZE(ello)/))

       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%real_handle,e)
       CALL MPI_ACCUMULATE(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,MPI_MIN,d%real_handle,e)
       CALL MPI_WIN_UNLOCK(pe,d%real_handle,e)
       DEALLOCATE(tmp_real)
    ENDIF
    IF (PRESENT(eig)) THEN
       ALLOCATE(tmp_real(d%size_eig))
       tmp_real=1E99
       n1=1;n3=1
       IF (PRESENT(n_rank)) n1=n_rank+1
       IF (PRESENT(n_size)) n3=n_size
       n2=SIZE(eig)*n3+n1-1
       tmp_real(n1:n2:n3)=eig
       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%eig_handle,e)
       IF (n3.ne.1) THEN
          CALL MPI_ACCUMULATE(tmp_real,d%size_eig,MPI_DOUBLE_PRECISION,pe,slot,d%size_eig,MPI_DOUBLE_PRECISION,MPI_MIN,d%eig_handle,e)
       ELSE
          CALL MPI_PUT(tmp_real,d%size_eig,MPI_DOUBLE_PRECISION,pe,slot,d%size_eig,MPI_DOUBLE_PRECISION,d%eig_handle,e)
       ENDIF
       CALL MPI_WIN_UNLOCK(pe,d%eig_handle,e)
       DEALLOCATE(tmp_real)
    ENDIF
    IF (PRESENT(z)) THEN
       tmp_size=SIZE(z,1)
       ALLOCATE(tmp_real(tmp_size))
       ALLOCATE(tmp_cmplx(tmp_size))
       DO n=1,SIZE(z,2)
          n1=n-1
          IF (PRESENT(n_size)) n1=n_size*n1
          IF (PRESENT(n_rank)) n1=n1+n_rank
          slot=d%slot_ev(nk,jspin,n1+1)
          pe=d%pe_ev(nk,jspin,n1+1)
          !print *, "PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_real(1)
          SELECT TYPE(z)
          TYPE IS(REAL)
#ifdef CPP_SOC
             tmp_cmplx=z(:,n)
             CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%zc_handle,e)
             CALL MPI_PUT(tmp_cmplx,tmp_size,MPI_DOUBLE_COMPLEX,pe,slot,tmp_size,MPI_DOUBLE_COMPLEX,d%zc_handle,e)
             CALL MPI_WIN_UNLOCK(pe,d%zc_handle,e)
#else
             tmp_real=z(:,n)
             CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%zr_handle,e)
             CALL MPI_PUT(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,d%zr_handle,e)
             CALL MPI_WIN_UNLOCK(pe,d%zr_handle,e)
#endif
          TYPE IS(COMPLEX)
             tmp_cmplx=z(:,n)
             CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%zc_handle,e)
             CALL MPI_PUT(tmp_cmplx,tmp_size,MPI_DOUBLE_COMPLEX,pe,slot,tmp_size,MPI_DOUBLE_COMPLEX,d%zc_handle,e)
             CALL MPI_WIN_UNLOCK(pe,d%zc_handle,e)
          END SELECT
       ENDDO
    ENDIF

#endif
  END SUBROUTINE write_eig

#ifdef CPP_MPI
  SUBROUTINE priv_put_data(pe,slot,DATA,handle)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: pe,slot
    CLASS(*),INTENT(IN) :: DATA(:)
    INTEGER,INTENT(IN)  :: handle

    INTEGER             :: len,e
    INTEGER,ALLOCATABLE :: int_tmp(:)
    REAL,ALLOCATABLE    :: real_tmp(:)
    COMPLEX,ALLOCATABLE:: cmplx_tmp(:)
    INCLUDE 'mpif.h'
    len=SIZE(DATA)
    SELECT TYPE(DATA)
    TYPE IS (INTEGER)
       ALLOCATE(int_tmp(len))
       int_tmp=DATA
       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,handle,e)
       CALL MPI_PUT(int_tmp,len,MPI_INTEGER,pe,int(slot,MPI_ADDRESS_KIND),len,MPI_INTEGER,handle,e)
       CALL MPI_WIN_UNLOCK(pe,handle,e)
    TYPE is (REAL)
       ALLOCATE(real_tmp(len))
       real_tmp=DATA
       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,handle,e)
       CALL MPI_PUT(real_tmp,len,MPI_DOUBLE_PRECISION,pe,int(slot,MPI_ADDRESS_KIND),len,MPI_DOUBLE_PRECISION,handle,e)
       CALL MPI_WIN_UNLOCK(pe,handle,e)
    TYPE is (COMPLEX)
       ALLOCATE(cmplx_tmp(len))
       cmplx_tmp=DATA
       CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,handle,e)
       CALL MPI_PUT(cmplx_tmp,len,MPI_DOUBLE_COMPLEX,pe,int(slot,MPI_ADDRESS_KIND),len,MPI_DOUBLE_COMPLEX,handle,e)
       CALL MPI_WIN_UNLOCK(pe,handle,e)
    END SELECT
  END SUBROUTINE priv_put_data

  SUBROUTINE priv_get_data(pe,slot,len,handle,idata,rdata,cdata)
    IMPLICIT NONE
    INTEGER,INTENT(IN)  :: pe,slot,len
    INTEGER,INTENT(OUT),optional :: iDATA(len)
    REAL,INTENT(OUT),optional    :: rDATA(len)
    COMPLEX,INTENT(OUT),optional :: cDATA(len)
    INTEGER,INTENT(IN)  :: handle

    INTEGER             :: e
    INTEGER,ALLOCATABLE :: int_tmp(:)
    REAL,ALLOCATABLE    :: real_tmp(:)
    COMPLEX,ALLOCATABLE:: cmplx_tmp(:)
    INCLUDE 'mpif.h'

    IF (present(idata)) THEN
       ALLOCATE(int_tmp(len))
       CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,handle,e)
       CALL MPI_GET(int_tmp,len,MPI_INTEGER,pe,int(slot,MPI_ADDRESS_KIND),len,MPI_INTEGER,handle,e)
       CALL MPI_WIN_UNLOCK(pe,handle,e)
       iDATA=int_tmp
    ELSE IF (PRESENT(rdata)) THEN
       ALLOCATE(real_tmp(len))
       CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,handle,e)
       CALL MPI_GET(real_tmp,len,MPI_DOUBLE_PRECISION,pe,int(slot,MPI_ADDRESS_KIND),len,MPI_DOUBLE_PRECISION,handle,e)
       CALL MPI_WIN_UNLOCK(pe,handle,e)
       rDATA=real_tmp
    ELSE IF (PRESENT(cdata)) THEN
       ALLOCATE(cmplx_tmp(len))
       CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,handle,e)
       CALL MPI_GET(cmplx_tmp,len,MPI_DOUBLE_COMPLEX,pe,int(slot,MPI_ADDRESS_KIND),len,MPI_DOUBLE_COMPLEX,handle,e)
       CALL MPI_WIN_UNLOCK(pe,handle,e)
       cDATA=cmplx_tmp
    ELSE
       call judft_error("BUG in priv_get_data")
    ENDIF

  END SUBROUTINE priv_get_data
#endif

  SUBROUTINE write_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(IN)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
    COMPLEX,INTENT(IN)           :: qstars(:,:,:,:)
    INTEGER,INTENT(IN)           :: ksym(:),jsym(:)
    REAL,INTENT(IN),OPTIONAL     :: mcd(:,:,:)
    REAL,INTENT(IN),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)
#ifdef CPP_MPI
    TYPE(t_data_MPI),POINTER :: d
    INTEGER:: pe,slot

    CALL priv_find_data(id,d)
    pe=d%pe_basis(nk,jspin)
    slot=d%slot_basis(nk,jspin)

    CALL priv_put_data(pe,slot,RESHAPE(qal,(/SIZE(qal)/)),d%qal_handle)
    CALL priv_put_data(pe,slot,RESHAPE(qvac,(/SIZE(qvac)/)),d%qvac_handle)
    CALL priv_put_data(pe,slot,RESHAPE(qis,(/SIZE(qis)/)),d%qis_handle)
    CALL priv_put_data(pe,slot,RESHAPE(qvlay,(/SIZE(qvlay)/)),d%qvlay_handle)
    CALL priv_put_data(pe,slot,RESHAPE(qstars,(/SIZE(qstars)/)),d%qstars_handle)
    CALL priv_put_data(pe,slot,RESHAPE(ksym,(/SIZE(ksym)/)),d%ksym_handle)
    CALL priv_put_data(pe,slot,RESHAPE(jsym,(/SIZE(jsym)/)),d%jsym_handle)
    IF (d%l_mcd.AND.PRESENT(mcd))  CALL priv_put_data(pe,slot,RESHAPE(mcd,(/SIZE(mcd)/)),d%mcd_handle)
    IF (d%l_orb.AND.PRESENT(qintsl)) THEN
       CALL priv_put_data(pe,slot,RESHAPE(qintsl,(/SIZE(qintsl)/)),d%qintsl_handle)
       CALL priv_put_data(pe,slot,RESHAPE(qmtsl,(/SIZE(qmtsl)/)),d%qmtsl_handle)
       CALL priv_put_data(pe,slot,RESHAPE(qmtp,(/SIZE(qmtp)/)),d%qmtp_handle)
       CALL priv_put_data(pe,slot,RESHAPE(orbcomp,(/SIZE(orbcomp)/)),d%orbcomp_handle)
    ENDIF
#endif
  END SUBROUTINE write_dos

  SUBROUTINE read_dos(id,nk,jspin,qal,qvac,qis,qvlay,qstars,ksym,jsym,mcd,qintsl,qmtsl,qmtp,orbcomp)
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: id,nk,jspin
    REAL,INTENT(out)              :: qal(:,:,:),qvac(:,:),qis(:),qvlay(:,:,:)
    COMPLEX,INTENT(out)           :: qstars(:,:,:,:)
    INTEGER,INTENT(out)           :: ksym(:),jsym(:)
    REAL,INTENT(out),OPTIONAL     :: mcd(:,:,:)
    REAL,INTENT(out),OPTIONAL     :: qintsl(:,:),qmtsl(:,:),qmtp(:,:),orbcomp(:,:,:)
#ifdef CPP_MPI
    TYPE(t_data_MPI),POINTER :: d
    INTEGER:: pe,slot

    CALL priv_find_data(id,d)
    pe=d%pe_basis(nk,jspin)
    slot=d%slot_basis(nk,jspin)

    CALL priv_get_data(pe,slot,SIZE(qal),d%qal_handle,rdata=qal)
    CALL priv_get_data(pe,slot,SIZE(qvac),d%qvac_handle,rdata=qvac)
    CALL priv_get_data(pe,slot,SIZE(qis),d%qis_handle,rdata=qis)
    CALL priv_get_data(pe,slot,SIZE(qvlay),d%qvlay_handle,rdata=qvlay)
    CALL priv_get_data(pe,slot,SIZE(qstars),d%qstars_handle,cdata=qstars)
    CALL priv_get_data(pe,slot,SIZE(ksym),d%ksym_handle,idata=ksym)
    CALL priv_get_data(pe,slot,SIZE(jsym),d%jsym_handle,idata=jsym)
    IF (d%l_mcd.AND.PRESENT(mcd))  CALL priv_get_data(pe,slot,SIZE(mcd),d%mcd_handle,rdata=mcd)
    IF (d%l_orb.AND.PRESENT(qintsl)) THEN
       CALL priv_get_data(pe,slot,SIZE(qintsl),d%qintsl_handle,rdata=qintsl)
       CALL priv_get_data(pe,slot,SIZE(qmtsl),d%qmtsl_handle,rdata=qmtsl)
       CALL priv_get_data(pe,slot,SIZE(qmtp),d%qmtp_handle,rdata=qmtp)
       CALL priv_get_data(pe,slot,SIZE(orbcomp),d%orbcomp_handle,rdata=orbcomp)
    ENDIF
#endif
  END SUBROUTINE read_dos


#ifdef CPP_MPI
  SUBROUTINE create_maps(d,isize,nkpts,jspins,neig,n_size)
    IMPLICIT NONE
    TYPE(t_data_MPI),INTENT(INOUT):: d
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
             pe=MOD(n1,n_members)*n_size+MOD(n,n_size)
             d%pe_ev(nk,j,n)=pe
             d%slot_ev(nk,j,n)=used(pe)
             used(pe)=used(pe)+1
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE create_maps


#endif
END MODULE m_eig66_mpi
