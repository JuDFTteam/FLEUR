module m_eig66_mem
#include "juDFT_env.h"
! Do the IO of the eig-file into memory
! The eig-file is split into four arrays:
! eig_int contains the basis-set information/integers (nv,nmat,ne,k1,k2,k3,kveclo)
! eig_real contains the basis-set information/real (el,evac,ello,bkpt,wtkpt)
! eig_eig contains the eigenvalues
! eig_vec contains the eigenvectors
! The record number is given by nrec=nk+(jspin-1)*nkpts
    use m_eig66_data
    implicit none
    contains

    subroutine priv_find_data(id,d)
        INTEGER,INTENT(IN)::id
        TYPE(t_data_mem),pointer,intent(out):: d

        class(t_data),pointer   ::dp
        call eig66_find_data(dp,id)
        select type(dp)
            type is (t_data_mem)
            d=>dp
            class default
            call judft_error("BUG: wrong datatype in eig66_mem")
        END SELECT
    END subroutine

    subroutine open_eig(id,nmat,neig,nkpts,jspins,lmax,nlo,ntype,l_create,nlotot,l_noco,filename)
    INTEGER, INTENT(IN) :: id,nmat,neig,nkpts,jspins,nlo,ntype,lmax,nlotot
    LOGICAL, INTENT(IN) :: l_noco,l_create
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: filename

    integer:: length
    TYPE(t_data_mem),pointer:: d
    call priv_find_data(id,d)

    if (allocated(d%eig_int)) then
          if (.not.l_create) THEN
             if (present(filename)) call priv_readfromfile()
             return
          endif
          call close_eig(id,.true.)

    endif

    call eig66_data_storedefault(d,jspins,nkpts,nmat,neig,lmax,nlotot,nlo,ntype)

    !d%eig_int
    length=3 !nv+nmat+ne
    length=length+nmat*3 !k1,k2,k3
    length=length+nlotot !kveclo

    allocate(d%eig_int(length,jspins*nkpts))

    !d%eig_real
    length=3+1+2 !bk,wk,evac
    length=length+lmax*ntype !el
    length=length+nlo*ntype  !ello
    ALLOCATE(d%eig_real(length,jspins*nkpts))
    !d%eig_eig
    length=jspins
    if (l_noco) length=1
    allocate(d%eig_eig(neig,jspins*nkpts))
    !d%eig_vec
#ifdef CPP_INVERSION
    allocate(d%eig_vecr(nmat*neig,length*nkpts))
#ifdef CPP_SOC
    call judft_error("SOC+INVERSION can not be used with eigenvalues stored in memory")
#endif
#else
    allocate(d%eig_vecc(nmat*neig,length*nkpts))
#endif
    if (present(filename)) call priv_readfromfile()
    contains
    subroutine priv_readfromfile()
    use m_eig66_da,ONLY:open_eig_IO=>open_eig,read_eig_IO=>read_eig,close_eig_IO=>close_eig
    integer:: jspin,nk,i,ii,iii,nv,tmp_id
    real   :: wk,bk3(3),evac(2)
    integer :: k1(nmat),k2(nmat),k3(nmat),kveclo(nlotot)
    real    :: eig(neig),ello(d%nlo,d%ntype),el(d%lmax,d%ntype)
#ifdef CPP_INVERSION
     real    :: z(nmat,neig)
#else
     complex :: z(nmat,neig)
#endif
     tmp_id=eig66_data_newid(DA_mode)
     call open_eig_IO(tmp_id,nmat,neig,nkpts,jspins,d%lmax,d%nlo,d%ntype,nlotot,.false.,filename)
        DO jspin=1,jspins
           DO nk=1,nkpts
              call read_eig_IO(tmp_id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z)
              call write_eig(id,nk,jspin,ii,ii,nv,i,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z)
           ENDDO
       ENDDO
   call close_eig_IO(tmp_id)
   end subroutine

    end subroutine open_eig

    subroutine close_eig(id,delete,filename)
    integer,intent(in)         :: id
    logical,intent(in),optional::delete
    character(len=*),optional,intent(in)::filename
    TYPE(t_data_mem),pointer:: d
    call priv_find_data(id,d)

    if (present(filename)) call priv_writetofile()

    if (present(delete)) THEN
      if (delete) THEN
       if (allocated(d%eig_int)) deallocate(d%eig_int)
       if (allocated(d%eig_real)) deallocate(d%eig_real)
       if (allocated(d%eig_eig)) deallocate(d%eig_eig)
       if (allocated(d%eig_vecr)) deallocate(d%eig_vecr)
       if (allocated(d%eig_vecc)) deallocate(d%eig_vecc)
      endif
    endif
    contains
    subroutine priv_writetofile()
    use m_eig66_DA,ONLY:open_eig_DA=>open_eig,write_eig_DA=>write_eig,close_eig_DA=>close_eig
    implicit none

    integer:: nlotot,nk,jspin,nv,i,ii,tmp_id
    real   :: wk,bk3(3),evac(2)
    integer :: k1(d%nmat),k2(d%nmat),k3(d%nmat),kveclo(size(d%eig_int,1)-3-3*d%nmat)
    real    :: eig(size(d%eig_eig,1)),ello(d%nlo,d%ntype),el(d%lmax,d%ntype)
#ifdef CPP_INVERSION
    real    :: z(d%nmat,size(d%eig_eig,1))
#else
    complex :: z(d%nmat,size(d%eig_eig,1))
#endif
    tmp_id=eig66_data_newid(DA_mode)
     call open_eig_DA(tmp_id,d%nmat,d%neig,d%nkpts,d%jspins,d%lmax,d%nlo,d%ntype,d%nlotot,.false.,filename)
        DO jspin=1,d%jspins
           DO nk=1,d%nkpts
              call read_eig(id,nk,jspin,nv,i,k1,k2,k3,bk3,wk,ii,eig,el,ello,evac,kveclo,z=z)
              call write_eig_DA(tmp_id,nk,jspin,ii,ii,nv,i,k1,k2,k3,bk3,wk,eig,el,ello,evac,nlotot,kveclo,z=z)
           ENDDO
       ENDDO
   call close_eig_DA(tmp_id)
   call eig66_remove_data(id)
    end subroutine
    end subroutine close_eig

    subroutine read_eig(id,nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,&
                        ello,evac,kveclo,n_start,n_end,z)
      implicit none
      INTEGER, INTENT(IN)            :: id,nk,jspin
      INTEGER, INTENT(OUT),OPTIONAL  :: nv,nmat
      INTEGER, INTENT(OUT),OPTIONAL  :: neig
      REAL,    INTENT(OUT),OPTIONAL  :: eig(:)
      INTEGER, INTENT(OUT),OPTIONAL  :: k1(:),k2(:),k3(:),kveclo(:)
      REAL,    INTENT(OUT),OPTIONAL  :: evac(:),ello(:,:),el(:,:)
      REAL,    INTENT(OUT),OPTIONAL  :: bk(:),wk
      INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
      CLASS(*),INTENT(OUT),OPTIONAL  :: z(:,:)

      INTEGER::nrec
      TYPE(t_data_mem),pointer:: d
      call priv_find_data(id,d)

      nrec=nk+(jspin-1)*d%nkpts
      ! data from d%eig_int
      if (present(nv)) nv=d%eig_int(1,nrec)
      if (present(nmat)) nmat=d%eig_int(2,nrec)
      if (present(neig)) then
          neig=d%eig_int(3,nrec)
      endif
      if (present(k1)) then
         if (.not.present(k2).or..not.present(k3)) call juDFT_error("BUG: always read k1,k2,k3")
         k1=d%eig_int(3+1:3+d%nmat,nrec)
         k2=d%eig_int(3+d%nmat+1:3+2*d%nmat,nrec)
         k3=d%eig_int(3+2*d%nmat+1:3+3*d%nmat,nrec)
      endif
      if (present(kveclo)) kveclo=d%eig_int(4+3*d%nmat:3+3*d%nmat+size(kveclo),nrec)

      !data from d%eig_real
      if (present(bk)) bk=d%eig_real(1:3,nrec)
      if (present(wk)) wk=d%eig_real(4,nrec)
      if (present(evac)) evac=d%eig_real(5:6,nrec)
      if (present(el)) el=reshape(d%eig_real(7:7+size(el)-1,nrec),shape(el))
      if (present(ello)) ello=reshape(d%eig_real(size(d%eig_real,1)-size(ello)+1:,nrec),shape(ello))

      !data from d%eig_eig
      if (present(eig)) THEN
           eig=0.0
           eig=d%eig_eig(:size(eig),nrec)
           !print *,"R-eig:",nrec,shape(eig)
           !print*,"R-eig(data):",shape(d%eig_eig)
           !print*,"R:",eig
      ENDIF
      !data from d%eig_vec

      if (present(z)) then
          write(*,*) "R-Z:",nrec,shape(z)
          SELECT TYPE(z)
             TYPE is (real)
               if (.not.allocated(d%eig_vecr)) call juDFT_error("BUG: can not read real vectors from memory")
               z=reshape(d%eig_vecr(:size(z),nrec),shape(z))
             TYPE is (complex)
                if (.not.allocated(d%eig_vecc)) call juDFT_error("BUG: can not read complex vectors from memory")
                z=reshape(d%eig_vecc(:size(z),nrec),shape(z))
          END SELECT
      endif
    end subroutine read_eig

    
    subroutine write_eig(id,nk,jspin,neig,neig_total,nv,nmat,k1,k2,k3,bk,wk, &
                         eig,el,ello,evac,                     &
                         nlotot,kveclo,n_size,n_rank,z)
     INTEGER, INTENT(IN)          :: id,nk,jspin
     INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
     REAL,    INTENT(IN),OPTIONAL :: wk
     INTEGER, INTENT(IN),OPTIONAL :: neig,neig_total,nv,nmat,nlotot
     INTEGER, INTENT(IN),OPTIONAL :: k1(:),k2(:),k3(:),kveclo(:)
     REAL,    INTENT(IN),OPTIONAL :: bk(3),eig(:),el(:,:)
     REAL,    INTENT(IN),OPTIONAL :: evac(:),ello(:,:)
     CLASS(*),INTENT(IN),OPTIONAL :: z(:,:)
     INTEGER::nrec
     TYPE(t_data_mem),pointer:: d
     call priv_find_data(id,d)

      nrec=nk+(jspin-1)*d%nkpts
      ! data from d%eig_int
      if (present(nv)) d%eig_int(1,nrec)=nv
      if (present(nmat)) d%eig_int(2,nrec)=nmat
      if (present(neig)) THEN
          if (present(neig_total)) THEN
                 if (neig.ne.neig_total) STOP "BUG in eig_mem"
                 d%eig_int(3,nrec)=neig_total
          else
              STOP "BUG2 in eig_mem"
          endif
      endif

      if (present(k1)) then
         if (.not.present(k2).or..not.present(k3)) call juDFT_error("BUG: always write k1,k2,k3")
         d%eig_int(3+1:3+d%nmat,nrec)=k1
         d%eig_int(3+d%nmat+1:3+2*d%nmat,nrec)=k2
         d%eig_int(3+2*d%nmat+1:3+3*d%nmat,nrec)=k3
      endif
      if (present(kveclo)) d%eig_int(4+3*d%nmat:3+3*d%nmat+size(kveclo),nrec)=kveclo

      !data from d%eig_real
      if (present(bk)) d%eig_real(1:3,nrec)=bk
      if (present(wk)) d%eig_real(4,nrec)=wk
      if (present(evac)) d%eig_real(5:6,nrec)=evac
      if (present(el)) d%eig_real(7:7+size(el)-1,nrec)=reshape(el,(/size(el)/))
      if (present(ello)) d%eig_real(size(d%eig_real,1)-size(ello)+1:,nrec)=reshape(ello,(/size(ello)/))
      !data from d%eig_eig
      if (present(eig)) then
               d%eig_eig(:size(eig),nrec)=eig
               !print*,"W-eig:",nrec,shape(eig)
               !print*,"W:",eig
      endif
      !data from d%eig_vec
      if (present(z)) then
         write(*,*) "W-Z:",nrec,shape(z)

         SELECT TYPE(z)
           TYPE IS (real)
              if (.not.allocated(d%eig_vecr)) call juDFT_error("BUG: can not write real vectors to memory")
              d%eig_vecr(:size(z),nrec)=reshape(real(z),(/size(z)/))
           TYPE IS(complex)
              if (.not.allocated(d%eig_vecc)) call juDFT_error("BUG: can not write complex vectors to memory")
              d%eig_vecc(:size(z),nrec)=reshape(cmplx(z),(/size(z)/))
         END SELECT
      endif


    end subroutine write_eig


end module m_eig66_mem
