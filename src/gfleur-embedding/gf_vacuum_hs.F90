!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! Module that contains code to diagonalize the vacuum Hamiltonian
Module m_gf_vacuum_hs
      use m_juDFT
	use m_juDFT
    implicit none
    private
	!defaults for vacuum grid
	integer,parameter:: nmzxy=100
	integer,parameter:: nmz=250
	real,parameter   :: delz=0.1

	complex,allocatable,save :: h_mat(:,:,:)
	real,allocatable,save    :: uz(:,:),udz(:,:),duz(:,:),dudz(:,:),ddnv(:,:),eigval(:,:)

    public gf_vacuum_diagonalize,gf_generate_embpot, gf_vacuum_getHS
	contains
    subroutine gf_vacuum_diagonalize(fmpi,vacuum,input,nococonv,lapw,lapw_gf,stars,cell,sym,kpts,enpara,jspins,nk,l_noco)
    ! Set up the vacuum Hamiltonian and diagonalize it, store resulting eigenvectors and eigenvalues
    ! for later use
	use m_gf_types
	use m_gf_iodop,only:gf_iodop_readvacuum,GF_POTFILE
    implicit none
    type(t_lapw),intent(in) ::  lapw
    type(t_lapw_gf),intent(in) :: lapw_gf
    type(t_mpi),intent(in)  ::  fmpi
    type(t_vacuum),intent(in) :: vacuum
    type(t_input),intent(in)  :: input
    type(t_nococonv),intent(in) :: nococonv
    type(t_stars),intent(in):: stars
    type(t_cell),intent(in) ::  cell
    type(t_sym),intent(in)  ::   sym
    type(t_kpts),intent(in) ::  kpts
    type(t_enpara),intent(in):: enpara
    integer,intent(in)      :: jspins,nk
    logical,intent(in)      :: l_noco

	integer :: n,nv2,jspin,jsps
	real    :: wronk,evac_array(2,jspins)
	complex,dimension(lapw_gf%nv2(jspins),lapw_gf%nv2(jspins)) :: tuuv,tudv,tduv,tddv
	complex,allocatable::s_mat(:,:)

    !workspace for LAPACK
    integer:: lwork ,liwork,lrwork,i,ii
    complex,allocatable :: work(:)
    real,allocatable ::    rwork(:)
    integer,allocatable :: iwork(:)
    real:: rwork_d(1)
    complex::work_d(1)
    integer::iwork_d(1)

    if (.true.) then !All PE calculate vacuum!
		nv2=lapw_gf%nv2_tot

	    if (allocated(h_mat)) then
    		deallocate(h_mat,eigval,uz,duz,udz,dudz,ddnv)
    	endif
    	allocate(h_mat(2*nv2,2*nv2,jspins),eigval(2*nv2,jspins))
    	allocate(uz(nv2,jspins),udz(nv2,jspins),duz(nv2,jspins),dudz(nv2,jspins),ddnv(nv2,jspins))
    	allocate(s_mat(2*nv2,2*nv2))
    	IF (l_noco) THEN
       	  jsps=1
        ELSE
          jsps=jspins
        ENDIF
	    do jspin=1,jsps
		!enpara%evac -0.1!
		   call gf_vacuum_getHS(fmpi,vacuum,input,nococonv,stars,lapw,lapw_gf,jspin,jspins,nk,l_noco,sym,kpts,enpara%evac(1,jspin),cell,h_mat(:,:,jspin),s_mat(:,:))
		   !Call Lapack to determine workspace
	 	   if (jspin==1) then
	 		  CALL zhegvd(1,'V','U',2*nv2,h_mat(:,:,jspin),2*nv2   &
                 ,s_mat,2*nv2,eigval,work_d,-1,rwork_d,-1,iwork_d             &
                 ,-1,n)
     		  lwork = work_d(1)
     		  ALLOCATE(work(lwork))
     		  lrwork=rwork_d(1)
     		  ALLOCATE(rwork(lrwork))
     		  liwork=iwork_d(1)
     		  ALLOCATE(iwork(liwork))
     	    endif
     	  !Diagonalize

     	  CALL zhegvd(1,'V','U',2*nv2,h_mat(:,:,jspin),2*nv2     &
               ,s_mat,2*nv2,eigval(:,jspin),work,lwork,rwork,lrwork,iwork          &
               ,liwork,n)
          if (n/=0) call juDFT_error("Could not diagonalize Vacuum Hamiltonian")

      enddo
	endif
	!No broadcasting needed as all PE do the vacuum
	!call priv_mpi_bcast_HS(mpi,jspins)

	end subroutine

	subroutine priv_mpi_bcast_HS(mpi,jspins)
	! subroutine to distribute diagonalized vacuum Hamiltonian on all PE
	use m_gf_types
	implicit none
	type(t_mpi),intent(in)::mpi
	integer,intent(in)    ::jspins
#ifdef CPP_MPI
	integer:: nv2,ierr(3)

	if (mpi%pe0) then
	   nv2=size(uz,1)
	endif
	call MPI_BCAST(nv2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    if (.not.mpi%pe0) then
       if (allocated(h_mat)) then
    	deallocate(h_mat,eigval,uz,duz,udz,dudz,ddnv)
       endif
       allocate(h_mat(2*nv2,2*nv2,jspins),eigval(2*nv2,jspins))
       allocate(uz(nv2,jspins),udz(nv2,jspins),duz(nv2,jspins),dudz(nv2,jspins),ddnv(nv2,jspins))
    endif
    CALL MPI_BCAST(h_mat,size(h_mat),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(eigval,size(eigval),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(uz,size(uz),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(udz,size(udz),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(duz,size(duz),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(dudz,size(dudz),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    CALL MPI_BCAST(ddnv,size(ddnv),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
	end subroutine

	subroutine 	gf_vacuum_getHS(fmpi,vacuum,input,nococonv,stars,lapw,lapw_gf,jspin,jspins,nk,l_noco,sym,kpts,enpara_vac,cell,h_mat,s_mat,l_sym)
	use m_gf_types
	use m_vacfun !modern interface
	use m_gf_iodop
    implicit none
    type(t_mpi),intent(in)::     fmpi
    type(t_vacuum),intent(in)::  vacuum
    type(t_input),intent(in)::   input
    type(t_nococonv),intent(in)::nococonv
    type(t_lapw),intent(in)::    lapw
    type(t_lapw_gf),intent(in):: lapw_gf
    type(t_stars),intent(in)::   stars
    type(t_cell),intent(in)::    cell
    type(t_sym),intent(in)::     sym
    type(t_kpts),intent(in)::    kpts
    real,intent(in)         ::  enpara_vac
    integer,intent(in)      ::   jspins,jspin,nk
    logical,intent(in)      ::   l_noco
    complex,intent(out)     ::   s_mat(:,:),h_mat(:,:)
    logical,optional,intent(in):: l_sym

	!locals
	logical :: symmetrize
	integer :: nv2,n,nn,jsp
	integer,allocatable :: kvac(:,:,:)
	complex,allocatable :: vxy4(:,:,:,:)
	real    :: wronk
	complex,dimension(lapw_gf%nv2_tot,lapw_gf%nv2_tot) :: tuuv,tudv,tduv,tddv
    complex,allocatable     ::   vz(:,:,:)
    real,allocatable        ::   vz_in(:,:)
    complex,allocatable     ::   vxy(:,:,:)
    real::evac_array(2,jspins)
    symmetrize=.true.
    if (present(l_sym)) symmetrize=l_sym
    allocate(vz(nmz,2,4))
	if (l_noco) then
	    allocate(vz_in(nmz,4),vxy(nmzxy,stars%ng2-1,3))
	else
	    allocate(vz_in(nmz,jspins),vxy(nmzxy,stars%ng2-1,jspins))
	endif
	!read the vacuum potential
	call gf_iodop_readvacuum(GF_POTFILE,vxy,vz_in)
	vz=0.0
	vz(:,1,:size(vz_in,2))=vz_in(:,:)

	nv2=lapw_gf%nv2(jspin)
    evac_array=enpara_vac

    !modern vacfun interface
    ALLOCATE(kvac(2,SIZE(lapw_gf%k1p,1),SIZE(lapw_gf%k1p,2)))
    kvac(1,:,:)=lapw_gf%k1p
    kvac(2,:,:)=lapw_gf%k2p
    ALLOCATE(vxy4(SIZE(vxy,1),SIZE(vxy,2),1,SIZE(vxy,3)))
    vxy4(:,:,1,:)=vxy
    call vacfun(fmpi,vacuum,stars,input,nococonv,jspin,jspin,           &
        	cell,1,evac_array,kpts%bk(:,nk),                                &
        	vxy4,vz,kvac,lapw_gf%nv2,                                       &
        	tuuv(:nv2,:nv2),tddv(:nv2,:nv2),tudv(:nv2,:nv2),tduv(:nv2,:nv2),&
        	uz(:nv2,:),duz(:nv2,:),udz(:nv2,:),dudz(:nv2,:),ddnv(:nv2,:),wronk)
    DEALLOCATE(kvac,vxy4)
     if (l_noco) then
            CALL juDFT_error("noco not yet supported in the gfleur port",&
                             calledby="gf_vacuum_hs")
     endif
     !Setup of Hamiltonian
     h_mat(:nv2,:nv2)=tuuv
     h_mat(nv2+1:,:nv2)=tduv
     h_mat(:nv2,nv2+1:)=tudv
     h_mat(nv2+1:,nv2+1:)=tddv
     !Add contribution of surface normal derivative
     jsp=jspin
     DO n=1,nv2
        nn=n
        if (n>lapw_gf%nv2(jspin)) then
             nn=n-lapw_gf%nv2(jspin)
             jsp=1
        endif
     	h_mat(n,n)=h_mat(n,n)+0.5*uz(nn,jsp)*duz(nn,jsp)
     	h_mat(nv2+n,nv2+n)=h_mat(nv2+n,nv2+n)+0.5*udz(nn,jsp)*dudz(nn,jsp)
     	if (symmetrize) THEN
     	    h_mat(nv2+n,n)=h_mat(nv2+n,n)+0.25*(uz(nn,jsp)*dudz(nn,jsp)+duz(nn,jsp)*udz(nn,jsp))
     	    h_mat(n,nv2+n)=h_mat(n,nv2+n)+0.25*(uz(nn,jsp)*dudz(nn,jsp)+duz(nn,jsp)*udz(nn,jsp))
     	else
     	    h_mat(nv2+n,n)=h_mat(n,nv2+n)+0.5*duz(nn,jsp)*udz(nn,jsp)
            h_mat(n,nv2+n)=h_mat(nv2+n,n)+0.5*uz(nn,jsp)*dudz(nn,jsp)
        endif
     enddo
     !Setup of Overlap
     s_mat=0.0
     DO n=1,nv2
     	s_mat(n,n)=1
     enddo
     DO n=1,lapw_gf%nv2(jspin)
        s_mat(nv2+n,nv2+n)=ddnv(n,jspin)
     enddo
     if (l_noco) then
       DO n=1,lapw_gf%nv2(2)
          s_mat(nv2+lapw_gf%nv2(1)+n,nv2+lapw_gf%nv2(1)+n)=ddnv(n,2)
       enddo
     endif



	end subroutine



	subroutine gf_generate_embpot(en,jspin,sigma)
	use m_gf_energies
	use m_gf_math
	implicit none
	integer,intent(in)::en,jspin
	complex,intent(out)::sigma(:,:)
	!Locals
	integer::j,nv2,i
	complex::z,enedeno
	complex,dimension(size(sigma,1),2*size(sigma,2)):: tmp,tmp1


	if (.not.allocated(h_mat)) call juDFT_error("Hamiltonian not diagonalized in Vacuum?")

    nv2=size(h_mat,1)/2
	if (size(sigma,1).ne.nv2) THEN
	     print *,size(sigma,1),size(h_mat,1)
	     call juDFT_error("BUG:Wrong size of Sigma?",calledby="gf_vacuum_hs%gf_generate_sigma")
	endif

	!Construct the projected Green's function
	z=gf_z(en,0)

	tmp1=0.0;tmp=0.0
	DO j=1,nv2
		enedeno = 1.0/(z-eigval(j,jspin))
    	tmp(:nv2,j)= h_mat(:nv2,j,jspin)*uz(:,jspin)+h_mat(nv2+1:,j,jspin)*udz(:,jspin)
        tmp1(:nv2,j) = tmp(:nv2,j)*enedeno
        enedeno = 1.0/(z-eigval(j+nv2,jspin))
    	tmp(:nv2,j+nv2)= h_mat(:nv2,j+nv2,jspin)*uz(:,jspin)+h_mat(nv2+1:,j+nv2,jspin)*udz(:,jspin)
        tmp1(:nv2,j+nv2) = tmp(:nv2,j+nv2)*enedeno
    ENDDO





	sigma=matmul(tmp1,transpose(conjg(tmp)))
    !Invert surface Green function for embedding potential

	sigma=mat_inverse(sigma)

	end subroutine



end module m_gf_vacuum_hs
