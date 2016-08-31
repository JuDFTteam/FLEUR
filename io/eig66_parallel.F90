module m_eig66_parallel
#include "juDFT_env.h"
! Do the IO of the eig-file into memory, different from simple approach by
! using k-point parallelism
! The eig-file is split into four arrays:
! eig_int contains the basis-set information/integers (nv,nmat,ne,k1,k2,k3,kveclo)
! eig_real contains the basis-set information/real (el,evac,ello,bkpt,wtkpt)
! eig_eig contains the eigenvalues
! eig_vec contains the eigenvectors
! The record number is given by nrec=nk+(jspin-1)*nkpts
    implicit none
    INTEGER,ALLOCATABLE :: eig_int(:,:)
    REAL,ALLOCATABLE    :: eig_real(:,:)
    REAL,ALLOCATABLE    :: eig_eig(:,:)
    REAL,ALLOCATABLE    :: eig_vecr(:,:)
    COMPLEX,ALLOCATABLE :: eig_vecc(:,:)
    integer             :: size_k,isize

    INTEGER:: ntypes,lmax,nlo,jspins,nkpts

    contains
    subroutine open_eig(nmat,neig,nkpts_in,jspins_in,lmax_in,nlo_in,ntype_in,l_create,nlotot,l_noco)
    INTEGER, INTENT(IN) :: nmat,neig,nkpts_in,jspins_in,nlo_in,ntype_in,lmax_in,nlotot
    LOGICAL, INTENT(IN) :: l_noco,l_create

    integer:: length,e

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,isize,e)

    if (mod(nkpts_in/isize).ne.0) call judft_error("Could not distribute k-points in eig66_parallel")

    if (allocated(eig_int)) then
          if (.not.l_create) return
          call close_eig(.true.)
    endif


    ntypes=ntype_in
    lmax=lmax_in
    nlo=nlo_in
    jspins=jspins_in
    nkpts=nkpts_in/isize

    !eig_int
    length=3 !nv+nmat+ne
    length=length+nmat*3 !k1,k2,k3
    length=length+nlotot !kveclo
    size_k=nmat
    allocate(eig_int(length,jspins*nkpts))

    !eig_real
    length=3+1+2 !bk,wk,evac
    length=length+lmax*ntypes !el
    length=length+nlo*ntypes  !ello
    ALLOCATE(eig_real(length,jspins*nkpts))
    !eig_eig
    length=jspins
    if (l_noco) length=1
    allocate(eig_eig(neig,length*nkpts))
    !eig_vec
#ifdef CPP_INVERSION
    allocate(eig_vecr(nmat*neig,length*nkpts))
#ifdef CPP_SOC
    call judft_error("SOC+INVERSION can not be used with eigenvalues stored in memory")
#endif
#else
    allocate(eig_vecc(nmat*neig,length*nkpts))
#endif
    end subroutine open_eig

    subroutine close_eig(delete)
    logical,intent(in),optional::delete
    if (present(delete)) THEN
      if (delete) THEN
       if (allocated(eig_int)) deallocate(eig_int)
       if (allocated(eig_real)) deallocate(eig_real)
       if (allocated(eig_eig)) deallocate(eig_eig)
       if (allocated(eig_vecr)) deallocate(eig_vecr)
       if (allocated(eig_vecc)) deallocate(eig_vecc)
      endif
    endif

    end subroutine close_eig

    subroutine read_eig(nk,jspin,nv,nmat,k1,k2,k3,bk,wk,neig,eig,el,&
                        ello,evac,kveclo,n_start,n_end,zr,zc)
      implicit none
      INTEGER, INTENT(IN)            :: nk,jspin
      INTEGER, INTENT(OUT),OPTIONAL  :: nv,nmat
      INTEGER, INTENT(OUT),OPTIONAL  :: neig
      REAL,    INTENT(OUT),OPTIONAL  :: eig(:)
      INTEGER, INTENT(OUT),OPTIONAL  :: k1(:),k2(:),k3(:),kveclo(:)
      REAL,    INTENT(OUT),OPTIONAL  :: evac(:),ello(:,:),el(:,:)
      REAL,    INTENT(OUT),OPTIONAL  :: bk(:),wk
      INTEGER, INTENT(IN),OPTIONAL   :: n_start,n_end
      REAL,    INTENT(OUT),OPTIONAL  :: zr(:,:)
      COMPLEX, INTENT(OUT),OPTIONAL  :: zc(:,:)

      INTEGER::nrec

      nrec=nk/isize+(jspin-1)*nkpts/isize+1
      ! data from eig_int
      if (present(nv)) nv=eig_int(1,nrec)
      if (present(nmat)) nmat=eig_int(2,nrec)
      if (present(neig)) then
          neig=eig_int(3,nrec)
      endif
      if (present(k1)) then
         if (.not.present(k2).or..not.present(k3)) call juDFT_error("BUG: always read k1,k2,k3")
         k1=eig_int(3+1:3+size_k,nrec)
         k2=eig_int(3+size_k+1:3+2*size_k,nrec)
         k3=eig_int(3+2*size_k+1:3+3*size_k,nrec)
      endif
      if (present(kveclo)) kveclo=eig_int(4+3*size_k:3+3*size_k+size(kveclo),nrec)

      !data from eig_real
      if (present(bk)) bk=eig_real(1:3,nrec)
      if (present(wk)) wk=eig_real(4,nrec)
      if (present(evac)) evac=eig_real(5:6,nrec)
      if (present(el)) el=reshape(eig_real(7:7+size(el)-1,nrec),shape(el))
      if (present(ello)) ello=reshape(eig_real(size(eig_real,1)-size(ello)+1:,nrec),shape(ello))

      !data from eig_eig
      if (present(eig)) THEN
           eig=0.0
           eig=eig_eig(:size(eig),nrec)
      ENDIF
      !data from eig_vec
      if (present(zr)) then
         write(*,*) "R-Z:",nrec,shape(zr)

         if (.not.allocated(eig_vecr)) call juDFT_error("BUG: can not read real vectors from memory")
         zr=reshape(eig_vecr(:size(zr),nrec),shape(zr))

      endif
      if (present(zc)) then
          write(*,*) "R-ZC:",nrec,shape(zc)

         if (.not.allocated(eig_vecc)) call juDFT_error("BUG: can not read complex vectors from memory")
         zc=reshape(eig_vecc(:size(zc),nrec),shape(zc))
      endif
    end subroutine read_eig

    
    subroutine write_eig(nk,jspin,neig,nv,nmat,k1,k2,k3,bk,wk, &
                         eig,el,ello,evac,                     &
                         nlotot,kveclo,n_size,n_rank,zr,zc)
     INTEGER, INTENT(IN)          :: nk,jspin
     INTEGER, INTENT(IN),OPTIONAL :: n_size,n_rank
     REAL,    INTENT(IN),OPTIONAL :: wk
     INTEGER, INTENT(IN),OPTIONAL :: neig,nv,nmat,nlotot
     INTEGER, INTENT(IN),OPTIONAL :: k1(:),k2(:),k3(:),kveclo(:)
     REAL,    INTENT(IN),OPTIONAL :: bk(3),eig(:),el(:,:)
     REAL,    INTENT(IN),OPTIONAL :: evac(:),ello(:,:)
     REAL,    INTENT(IN),OPTIONAL :: zr(:,:)
     COMPLEX,INTENT(IN),OPTIONAL  :: zc(:,:)
     INTEGER::nrec

      nrec=(nk+(jspin-1)*nkpts)/isize+1
      ! data from eig_int
      if (present(nv)) eig_int(1,nrec)=nv
      if (present(nmat)) eig_int(2,nrec)=nmat
      if (present(neig)) THEN
                 eig_int(3,nrec)=neig
      endif
      if (present(k1)) then
         if (.not.present(k2).or..not.present(k3)) call juDFT_error("BUG: always write k1,k2,k3")
         eig_int(3+1:3+size_k,nrec)=k1
         eig_int(3+size_k+1:3+2*size_k,nrec)=k2
         eig_int(3+2*size_k+1:3+3*size_k,nrec)=k3
      endif
      if (present(kveclo)) eig_int(4+3*size_k:3+3*size_k+size(kveclo),nrec)=kveclo

      !data from eig_real
      if (present(bk)) eig_real(1:3,nrec)=bk
      if (present(wk)) eig_real(4,nrec)=wk
      if (present(evac)) eig_real(5:6,nrec)=evac
      if (present(el)) eig_real(7:7+size(el)-1,nrec)=reshape(el,(/size(el)/))
      if (present(ello)) eig_real(size(eig_real,1)-size(ello)+1:,nrec)=reshape(ello,(/size(ello)/))
      !data from eig_eig
      if (present(eig)) eig_eig(:size(eig),nrec)=eig
      !data from eig_vec
      if (present(zr)) then
         write(*,*) "W-ZR:",nrec,shape(zr)

         if (.not.allocated(eig_vecr)) call juDFT_error("BUG: can not write real vectors to memory")
         eig_vecr(:size(zr),nrec)=reshape(zr,(/size(zr)/))
      endif
      if (present(zc)) then
         write(*,*) "W-ZC:",nrec,shape(zc)

         if (.not.allocated(eig_vecc)) call juDFT_error("BUG: can not write complex vectors to memory")
         eig_vecc(:size(zc),nrec)=reshape(zc,(/size(zc)/))
      endif


    end subroutine write_eig


end module m_eig66_mem
