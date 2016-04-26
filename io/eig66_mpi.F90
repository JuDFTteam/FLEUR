MODULE m_eig66_mpi
#include "juDFT_env.h"
    use m_eig66_data
    IMPLICIT NONE
#ifdef CPP_MPI
    include 'mpif.h'
    PRIVATE

#endif



    PUBLIC open_eig,read_eig,write_eig,close_eig
CONTAINS

  subroutine priv_find_data(id,d)
        INTEGER,INTENT(IN)::id
        TYPE(t_data_mpi),POINTER:: d

        class(t_data),pointer   ::dp
        call eig66_find_data(dp,id)
        select type(dp)
            type is (t_data_mpi)
            d=>dp
            class default
            call judft_error("BUG: wrong datatype in eig66_mpi")
        END SELECT
    END subroutine
    

    SUBROUTINE open_eig(id,mpi_comm,nmat,neig,nkpts,jspins,lmax,nlo,ntype,create,nlotot,l_noco,n_size_opt,filename)
        USE,INTRINSIC::iso_c_binding
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: id,mpi_comm,nmat,neig,nkpts,jspins,nlo,ntype,lmax,nlotot
        LOGICAL, INTENT(IN) :: l_noco,create
        INTEGER,INTENT(IN),OPTIONAL:: n_size_opt
        CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename
#ifdef CPP_MPI
        INTEGER:: isize,e,slot_size,type_size,local_slots
        INTEGER(MPI_ADDRESS_KIND):: length
        TYPE(c_ptr)::ptr
        TYPE(t_data_MPI),pointer :: d

        call priv_find_data(id,d)

        IF (present(n_size_opt)) d%n_size=n_size_opt
        IF (allocated(d%pe_ev)) THEN
            IF (create) THEN
                d%neig_data=0
                d%eig_data=1E99
                d%int_data=9999999
                d%real_data=1E99
#ifdef CPP_INVERSION
                d%zr_data=0.0
#else
                d%zc_data=0.0
#endif
            ENDIF
            if (present(filename)) call priv_readfromfile()
            RETURN !everything already done!
        ENDIF

        CALL timestart("preperation")
        CALL MPI_COMM_RANK(MPI_COMM,d%irank,e)
        CALL MPI_COMM_SIZE(MPI_COMM,isize,e)
        d%jspins=jspins
        d%nkpts=nkpts
        d%nmat=nmat
        d%neig=neig
        d%nlotot=nlotot
        d%lmax=lmax
        d%nlo=nlo
        d%ntype=ntype


        CALL timestart("map creation")
        CALL create_maps(d,isize,nkpts,jspins,neig,d%n_size)
        CALL timestop("map creation")
        local_slots=count(d%pe_basis==d%irank)
        !Now create the windows

        !Window for neig
        CALL MPI_TYPE_SIZE(MPI_INTEGER,type_size,e)
        slot_size=1*type_size
        length=slot_size*local_slots
        CALL timestop("preperation")
        !Allocate memory
        CALL timestart("Allocate")
        CALL MPI_ALLOC_MEM(length,MPI_INFO_NULL,ptr,e)
        CALL c_f_pointer(ptr,d%neig_data,(/length/type_size/))
        IF (e.NE.0) CPP_error("Could not allocated MPI-Data (neig) in eig66_mpi")
        !create window
        d%neig_data=0
        CALL timestop("Allocate")
        CALL timestart("win create")
        CALL MPI_WIN_CREATE(d%neig_data,length,slot_size,Mpi_INFO_NULL,MPI_COMM, d%neig_handle, e)
        CALL timestop("win create")


        !The integer values
        d%size_k=nmat
        slot_size=(5+3*d%size_k+1+nlotot)*type_size
        length=slot_size*local_slots
        !Allocate memory
        CALL timestart("Allocate")
        CALL MPI_ALLOC_MEM(length,MPI_INFO_NULL,ptr,e)
        CALL c_f_pointer(ptr,d%int_data,(/length/type_size/))
        IF (e.NE.0) CPP_error("Could not allocated MPI-Data (int) in eig66_mpi")
        !create window
        d%int_data=9999999
        CALL timestop("Allocate")
        CALL timestart("win create")
        CALL MPI_WIN_CREATE(d%int_data, length,slot_size,Mpi_INFO_NULL, MPI_COMM, d%int_handle, e)
        CALL timestop("win create")
        CALL timestart("Allocate")

        !The real values
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,type_size,e)
        d%size_el=(1+lmax)*ntype
        d%size_ello=nlo*ntype
        slot_size=(6+d%size_el+d%size_ello)*type_size
        length=slot_size*local_slots
        !Allocate memory
        CALL MPI_ALLOC_MEM(length,MPI_INFO_NULL,ptr,e)
        CALL c_f_pointer(ptr,d%real_data,(/length/type_size/))
        IF (e.NE.0) CPP_error("Could not allocated MPI-Data (real) in eig66_mpi")
        !create window
        d%real_data=1E99
        CALL timestop("Allocate")
        CALL timestart("win create")
        CALL MPI_WIN_CREATE(d%real_data, length,slot_size,Mpi_INFO_NULL, MPI_COMM, d%real_handle, e)
        CALL timestop("win create")
        CALL timestart("Allocate")

        !The eigenvalues
        d%size_eig=neig
        slot_size=d%size_eig*type_size
        length=slot_size*local_slots
        !Allocate memory
        CALL MPI_ALLOC_MEM(length,MPI_INFO_NULL,ptr,e)
        CALL c_f_pointer(ptr,d%eig_data,(/length/type_size/))
        IF (e.NE.0) CPP_error("Could not allocated MPI-Data (real) in eig66_mpi")
        !create window
        d%eig_data=1E99
        CALL timestop("Allocate")
        CALL timestart("win create")
        CALL MPI_WIN_CREATE(d%eiG_data, length,slot_size,Mpi_INFO_NULL, MPI_COMM, d%eig_handle, e)
        CALL timestop("win create")
        !The eigenvectors
#ifndef CPP_INVERSION
        CALL MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,type_size,e)
#endif
        local_slots=count(d%pe_ev==d%irank)
        slot_size=nmat*type_size
        length=slot_size*local_slots
        !Allocate memory
        CALL timestart("Allocate")
        CALL MPI_ALLOC_MEM(length,MPI_INFO_NULL,ptr,e)
        IF (e.NE.0) CPP_error("Could not allocated MPI-Data (real) in eig66_mpi")

#ifdef CPP_INVERSION
        !create window

        CALL c_f_pointer(ptr,d%zr_data,(/length/type_size/))
        CALL timestop("Allocate")
        CALL timestart("win create")
        CALL MPI_WIN_CREATE(d%zr_data, length,slot_size,Mpi_INFO_NULL, MPI_COMM,d%zr_handle, e)
        CALL timestop("win create")

#ifdef CPP_SOC
        CALL judft_error("SOC+INVERSION can not be used with eigenvalues stored in memory")
#endif
#else
        !create window

        CALL c_f_pointer(ptr,d%zc_data,(/length/type_size/))
        CALL timestop("Allocate")
        CALL timestart("win create")

        CALL MPI_WIN_CREATE(d%zc_data, length,slot_size,Mpi_INFO_NULL, MPI_COMM, d%zc_handle, e)
        CALL timestop("win create")
#endif
        if (present(filename).and..not.create) call priv_readfromfile()
        contains
        subroutine priv_readfromfile()
        use m_eig66_DA,ONLY:open_eig_DA=>open_eig,read_eig_DA=>read_eig,close_eig_da=>close_eig
        integer:: jspin,nk,i,ii,iii,nv,tmp_id
        real   :: wk,bk3(3),evac(2)
        integer :: k1(nmat),k2(nmat),k3(nmat),kveclo(nlotot)
        real    :: eig(neig),ello(nlo,ntype),el(lmax,ntype)
#ifdef CPP_INVERSION
        real    :: z(nmat,neig)
#else
        complex :: z(nmat,neig)
#endif
        !only do this with PE=0
        if (d%irank==0) THEN
            tmp_id=eig66_data_newid(DA_mode)
            call open_eig_DA(tmp_id,nmat,neig,nkpts,jspins,lmax,nlo,ntype,nlotot,.false.,filename)
            DO jspin=1,jspins
                DO nk=1,nkpts
                    call read_eig_DA(tmp_id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z)
                    call write_eig(id,nk,jspin,ii,ii,nv,nmat,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z)
                ENDDO
            ENDDO
            call close_eig_DA(tmp_id)
        ENDIF
        end subroutine
#endif
    END SUBROUTINE open_eig
    SUBROUTINE close_eig(id,delete,filename)
    INTEGER,INTENT(IN)         :: id
    LOGICAL,INTENT(IN),OPTIONAL:: delete
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL::filename
    TYPE(t_data_MPI),pointer :: d
    call priv_find_data(id,d)

    IF (present(delete)) THEN
         IF (delete) WRITE(*,*) "No deallocation of memory implemented in eig66_mpi"
    ENDIF
    if (present(filename)) call priv_writetofile()
    contains
    subroutine priv_writetofile()
    use m_eig66_DA,ONLY:open_eig_DA=>open_eig,write_eig_DA=>write_eig,close_eig_DA=>close_eig
    implicit none

    integer:: nlotot,nk,jspin,nv,i,ii,tmp_id
    real   :: wk,bk3(3),evac(2)
    integer :: k1(d%nmat),k2(d%nmat),k3(d%nmat),kveclo(d%nlotot)
    real    :: eig(d%neig),ello(d%nlo,d%ntype),el(d%lmax,d%ntype)
#ifdef CPP_INVERSION
    real    :: z(d%nmat,d%neig)
#else
    complex :: z(d%nmat,d%neig)
#endif
    nlotot=d%nlotot

    if (d%irank==0) THEN
    tmp_id=eig66_data_newid(DA_mode)
    call open_eig_DA(tmp_id,d%nmat,d%neig,d%nkpts,d%jspins,d%lmax,d%nlo,d%ntype,d%nlotot,.false.,filename)
    DO jspin=1,d%jspins
         DO nk=1,d%nkpts
             call read_eig(id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z)
             call write_eig_DA(tmp_id,nk,jspin,ii,ii,nv,i,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z)

         ENDDO
    ENDDO
    CALL close_eig_DA(tmp_id)
    ENDIF
    call eig66_remove_data(id)
    end subroutine

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
        TYPE(t_data_MPI),pointer :: d
        call priv_find_data(id,d)
        pe=d%pe_basis(nk,jspin)
        slot=d%slot_basis(nk,jspin)
        IF (present(neig))THEN
            CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%neig_handle,e)
            ! Get current values
            CALL  MPI_GET(neig,1,MPI_INTEGER,pe,slot,1,MPI_INTEGER,d%neig_handle,e)
            CALL MPI_WIN_UNLOCK(pe,d%neig_handle,e)

        ENDIF
        !read the integer values
        IF (any((/present(nv),present(nmat),present(k1),present(k2),present(k3),present(kveclo)/))) THEN
            tmp_size=4+3*d%size_k
            IF (present(kveclo)) tmp_size=tmp_size+size(kveclo)
            ALLOCATE(tmp_int(tmp_size))
            CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%int_handle,e)
            ! Get current values
            CALL  MPI_GET(tmp_int,tmp_size,MPI_INTEGER,pe,slot,tmp_size,MPI_INTEGER,d%int_handle,e)
            CALL MPI_WIN_UNLOCK(pe,d%int_handle,e)
            !IF (present(neig)) neig=tmp_int(1)
            IF (present(nv))   nv=tmp_int(2)
            IF (present(nmat)) nmat=tmp_int(3)
            IF (present(k1))   k1=tmp_int(4+1:4+size(k1))
            IF (present(k2))   k2=tmp_int(4+d%size_k+1:4+d%size_k+size(k2))
            IF (present(k3))   k3=tmp_int(4+2*d%size_k+1:4+2*d%size_k+size(k3))
            IF (present(kveclo)) kveclo=tmp_int(4+3*d%size_k+1:4+3*d%size_k+size(kveclo))

        ENDIF
        !read the real-values
        IF (any((/present(wk),present(bk3),present(el),present(ello),present(evac)/))) THEN
            tmp_size=6+d%size_el+d%size_ello
            ALLOCATE(tmp_real(tmp_size))
            CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%real_handle,e)
            ! Get current values
            CALL  MPI_GET(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,d%real_handle,e)
            CALL MPI_WIN_UNLOCK(pe,d%real_handle,e)
            IF (present(wk))   wk=tmp_real(1)
            IF (present(bk3))  bk3=tmp_real(2:4)
            IF (present(evac)) evac=tmp_real(5:6)
            IF (present(el))   el=reshape(tmp_real(6+1:6+size(el)),shape(el))
            IF (present(ello)) ello=reshape(tmp_real(6+d%size_el+1:6+d%size_el+size(ello)),shape(ello))
            DEALLOCATE(tmp_real)
        ENDIF
        IF (present(eig)) THEN
            CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%eig_handle,e)
            ALLOCATE(tmp_real(d%size_eig))
            CALL MPI_GET(tmp_real,d%size_eig,MPI_DOUBLE_PRECISION,pe,slot,d%size_eig,MPI_DOUBLE_PRECISION,d%eig_handle,e)
            CALL MPI_WIN_UNLOCK(pe,d%eig_handle,e)
            n1=1;n3=1;n2=size(eig)
            IF (present(n_start)) n1=n_start
            IF (present(n_end)) n2=n_end
            eig(:n2-n1+1)=tmp_real(n1:n2)
            DEALLOCATE(tmp_real)
        ENDIF

        IF (present(z)) THEN
            tmp_size=size(z,1)
            ALLOCATE(tmp_real(tmp_size))
            ALLOCATE(tmp_cmplx(tmp_size))
            DO n=1,size(z,2)
                n1=n
                IF (present(n_start)) n1=n_start+n-1
                IF (present(n_end)) THEN
                    IF (n1>n_end) CYCLE
                ENDIF
                slot=d%slot_ev(nk,jspin,n1)
                pe=d%pe_ev(nk,jspin,n1)
                SELECT TYPE(z)
                  TYPE IS(real)
                     CALL MPI_WIN_LOCK(MPI_LOCK_SHARED,pe,0,d%zr_handle,e)
                     CALL MPI_GET(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,d%zr_handle,e)
                     CALL MPI_WIN_UNLOCK(pe,d%zr_handle,e)
                     !print *, nk,jspin,n1,"r PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_real(1)
                     z(:,n)=tmp_real
                  TYPE IS (complex)
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
        TYPE(t_data_MPI),pointer :: d

        call priv_find_data(id,d)

        pe=d%pe_basis(nk,jspin)
        slot=d%slot_basis(nk,jspin)
        !write the number of eigenvalues values
        IF (present(neig_total)) THEN
            CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%neig_handle,e)
            ALLOCATE(tmp_int(1))
            tmp_int(1)=neig_total
            CALL MPI_PUT(tmp_int,1,MPI_INTEGER,pe,slot,1,MPI_INTEGER,d%neig_handle,e)
            CALL MPI_WIN_UNLOCK(pe,d%neig_handle,e)
            DEALLOCATE(tmp_int)
        ENDIF

        IF (any((/present(nv),present(nmat),present(nlotot),present(k1),present(k2),present(k3),present(kveclo)/))) THEN
            tmp_size=5+3*d%size_k
            IF (present(kveclo)) tmp_size=tmp_size+size(kveclo)
            ALLOCATE(tmp_int(tmp_size))
            tmp_int=9999999
            tmp_int(1)=0
            IF (present(nv))   tmp_int(2)=nv
            IF (present(nmat)) tmp_int(3)=nmat
            IF (present(nlotot)) tmp_int(4)=nlotot
            IF (present(k1))   tmp_int(4+1:4+size(k1))=k1
            IF (present(k2))   tmp_int(4+d%size_k+1:4+d%size_k+size(k2))=k2
            IF (present(k3))   tmp_int(4+2*d%size_k+1:4+2*d%size_k+size(k3))=k3
            IF (present(kveclo)) tmp_int(4+3*d%size_k+1:4+3*d%size_k+size(kveclo))=kveclo
            CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%int_handle,e)
            CALL MPI_ACCUMULATE(tmp_int,tmp_size,MPI_INTEGER,pe,slot,tmp_size,MPI_INTEGER,MPI_MIN,d%int_handle,e)
            CALL MPI_WIN_UNLOCK(pe,d%int_handle,e)
        ENDIF
        !write the real-values
        IF (any((/present(wk),present(bk3),present(el),present(ello),present(evac)/))) THEN
            tmp_size=6+d%size_el+d%size_ello
            ALLOCATE(tmp_real(tmp_size))
            tmp_real=1E99
            IF (present(wk))   tmp_real(1)=wk
            IF (present(bk3))  tmp_real(2:4)=bk3
            IF (present(evac)) tmp_real(5:6)=evac
            IF (present(el))   tmp_real(6+1:6+size(el))=reshape(el,(/size(el)/))
            IF (present(ello)) tmp_real(6+d%size_el+1:6+d%size_el+size(ello))=reshape(ello,(/size(ello)/))

            CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%real_handle,e)
            CALL MPI_ACCUMULATE(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,MPI_MIN,d%real_handle,e)
            CALL MPI_WIN_UNLOCK(pe,d%real_handle,e)
            DEALLOCATE(tmp_real)
        ENDIF
        IF (present(eig)) THEN
            ALLOCATE(tmp_real(d%size_eig))
            tmp_real=1E99
            n1=1;n3=1
            IF (present(n_rank)) n1=n_rank+1
            IF (present(n_size)) n3=n_size
            n2=size(eig)*n3+n1-1
            tmp_real(n1:n2:n3)=eig
            CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%eig_handle,e)
            CALL MPI_ACCUMULATE(tmp_real,d%size_eig,MPI_DOUBLE_PRECISION,pe,slot,d%size_eig,MPI_DOUBLE_PRECISION,MPI_MIN,d%eig_handle,e)
            CALL MPI_WIN_UNLOCK(pe,d%eig_handle,e)
            DEALLOCATE(tmp_real)
        ENDIF
        IF (present(z)) THEN
            tmp_size=size(z,1)
            ALLOCATE(tmp_real(tmp_size))
            ALLOCATE(tmp_cmplx(tmp_size))
            DO n=1,size(z,2)
                n1=n-1
                IF (present(n_size)) n1=n_size*n1
                IF (present(n_rank)) n1=n1+n_rank
                slot=d%slot_ev(nk,jspin,n1+1)
                pe=d%pe_ev(nk,jspin,n1+1)
                !print *, "PE:",pe," Slot: ",slot," Size:",tmp_size,tmp_real(1)
                SELECT TYPE(z)
                TYPE IS(REAL)
                    tmp_real=z(:,n)
                    CALL MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,pe,0,d%zr_handle,e)
                    CALL MPI_PUT(tmp_real,tmp_size,MPI_DOUBLE_PRECISION,pe,slot,tmp_size,MPI_DOUBLE_PRECISION,d%zr_handle,e)
                    CALL MPI_WIN_UNLOCK(pe,d%zr_handle,e)
               TYPE IS(complex)
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
                    pe=MOD(n1,n_members)*n_size+mod(n,n_size)
                    d%pe_ev(nk,j,n)=pe
                    d%slot_ev(nk,j,n)=used(pe)
                    used(pe)=used(pe)+1
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE


#endif
END MODULE m_eig66_mpi
