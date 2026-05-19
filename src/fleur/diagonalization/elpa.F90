!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_elpa
   use m_judft
   use m_constants
   use m_types_solver
   use m_types_mpimat
   use m_types_mat
#ifdef CPP_ELPA
   use elpa
#endif     
#ifdef CPP_MPI 
   use mpi 
#endif
   implicit none
   private
   type, extends(t_solver):: t_solver_elpa
   contains
      procedure        :: solve_gev => elpa_gev  !solver for generalized eigenvalue problem
      procedure        :: solve_std_dp => elpa_diag
      procedure        :: solve_std_sp => elpa_diag_sp
      procedure        :: to_std => elpa_to_std
      procedure        :: backtrans => elpa_recover  !transform the Eigenvalue back to the generalized problem
   end type
#ifdef CPP_ELPA   
   class(elpa_t), pointer :: elpa_obj
#endif   
   logical,save:: firstcall=.true.   
   public get_solver_elpa

contains

#ifdef CPP_ELPA
   subroutine check_elpa_err(err, where)
      implicit none
      integer, intent(in) :: err
      character(len=*), intent(in) :: where
      character(len=256) :: errmsg

      if (err /= ELPA_OK) then
         write(errmsg, '(a,a,a,i0)') 'ELPA call failed in ', trim(where), ', err=', err
         call juDFT_error(trim(errmsg), calledby='elpa')
      end if
   end subroutine check_elpa_err

   subroutine dump_mat_layout(where, mat)
      implicit none
      character(len=*), intent(in) :: where
      class(t_mat), intent(in) :: mat
      integer :: rank, ierr

      select type (mat)
      type is (t_mpimat)
         rank = -1
#ifdef CPP_MPI
         call MPI_COMM_RANK(mat%blacsdata%mpi_com, rank, ierr)
#endif
         write(oUnit, '(a,1x,a,1x,a,1x,i0)') 'ELPA layout dump', trim(where), 'rank', rank
         write(oUnit, '(a,2(1x,i0),a,2(1x,i0),a,2(1x,i0))') '  global', mat%global_size1, mat%global_size2, &
                                                            'local', mat%matsize1, mat%matsize2, &
                                                            'grid', mat%blacsdata%nprow, mat%blacsdata%npcol
         write(oUnit, '(a,2(1x,i0),a,4(1x,i0))') '  prow/pcol', mat%blacsdata%myrow, mat%blacsdata%mycol, &
                                                 'desc', mat%blacsdata%blacs_desc(5), mat%blacsdata%blacs_desc(6), &
                                                 mat%blacsdata%blacs_desc(7), mat%blacsdata%blacs_desc(8)
      class default
         write(oUnit, '(a,1x,a)') 'ELPA layout dump', trim(where)
         write(oUnit, '(a,l1,a,2(1x,i0))') '  real', mat%l_real, 'size', mat%matsize1, mat%matsize2
      end select
      call flush(oUnit)
   end subroutine dump_mat_layout

   subroutine dump_elpa_obj_settings(where)
      implicit none
      character(len=*), intent(in) :: where
      integer :: rank, ierr, err
      integer :: na, nev, nrows, ncols, nblk, prow, pcol

      rank = -1
#ifdef CPP_MPI
      if (associated(elpa_obj)) then
         call elpa_obj%get("mpi_comm_parent", ierr, err)
         if (err == ELPA_OK) call MPI_COMM_RANK(ierr, rank, err)
      end if
#endif

      if (.not.associated(elpa_obj)) then
         write(oUnit, '(a,1x,a,1x,a,1x,i0)') 'ELPA obj dump', trim(where), 'rank', rank
         write(oUnit, '(a)') '  elpa_obj not associated'
         return
      end if

      call elpa_obj%get("na", na, err)
      call elpa_obj%get("nev", nev, err)
      call elpa_obj%get("local_nrows", nrows, err)
      call elpa_obj%get("local_ncols", ncols, err)
      call elpa_obj%get("nblk", nblk, err)
      call elpa_obj%get("process_row", prow, err)
      call elpa_obj%get("process_col", pcol, err)

      write(oUnit, '(a,1x,a,1x,a,1x,i0)') 'ELPA obj dump', trim(where), 'rank', rank
      write(oUnit, '(a,2(1x,i0),a,2(1x,i0),a,3(1x,i0))') '  na/nev', na, nev, 'local', nrows, ncols, 'nblk/prow/pcol', nblk, prow, pcol
      call flush(oUnit)
   end subroutine dump_elpa_obj_settings
#endif
   function get_solver_elpa() result(solver)
      class(t_solver), allocatable :: solver
      allocate(t_solver_elpa :: solver)
      solver%name = "elpa"
#ifdef CPP_ELPA
      solver%available = .true.
#else
      solver%available = .false.
#endif
      solver%parallel = .true.
      solver%serial = .false.
      solver%generalized = .true.
      solver%standard = .true.      
#ifdef CPP_ELPA_PATCH
      solver%transform = .true.
#else
      solver%transform = .false.
#endif   
      solver%GPU = .true.
#ifdef CPP_ELPA_SP
solver%single_precision = .true.
#else      
   solver%single_precision = .false.
#endif      
      solver%use_sp = .false.
   end function



   subroutine create_elpa_obj(hmat, ne)
      !$ use omp_lib
      implicit none
      class(t_mat), intent(IN)              :: hmat
      integer, intent(IN)                  :: ne
      
#ifdef CPP_ELPA
      integer           :: np, myid
      integer           :: err
      TYPE(t_mpimat) :: tmp
      
      if (firstcall) then
         err = elpa_init(20211125)
         call check_elpa_err(err, 'elpa_init')
         firstcall = .false.
         elpa_obj=>null()
      end if
      if (associated(elpa_obj)) return

      call timestart("ELPA SETUP")
      elpa_obj => elpa_allocate(err)
      call check_elpa_err(err, 'elpa_allocate')
      if (.not.associated(elpa_obj)) call juDFT_error('ELPA allocation returned null object', calledby='elpa')

      !Some settings are set for all matrices
      call elpa_obj%set("local_nrows", hmat%matsize1, err)
      call check_elpa_err(err, 'set(local_nrows)')
      call elpa_obj%set("local_ncols", hmat%matsize2, err)
      call check_elpa_err(err, 'set(local_ncols)')

      !$ call elpa_obj%set("omp_threads", omp_get_num_threads(),err)
      call elpa_obj%set("timings", 1, err)
      call check_elpa_err(err, 'set(timings)')
      !Some other settings depend on matrix type
      select type (hmat)
      type is (t_mpimat)
         call MPI_BARRIER(hmat%blacsdata%mpi_com, err)
         call MPI_COMM_SIZE(hmat%blacsdata%mpi_com, np, err)
         call MPI_COMM_RANK(hmat%blacsdata%mpi_com, myid, err)
         ! Blocking factor
         if (hmat%blacsdata%blacs_desc(5) .ne. hmat%blacsdata%blacs_desc(6)) &
            call judft_error("Different block sizes for rows/columns not supported")
         call elpa_obj%set("na", hmat%global_size1, err)
         call check_elpa_err(err, 'set(na)')
         call elpa_obj%set("nblk", hmat%blacsdata%blacs_desc(5), err)
         call check_elpa_err(err, 'set(nblk)')
         call elpa_obj%set("mpi_comm_parent", hmat%blacsdata%mpi_com, err)
         call check_elpa_err(err, 'set(mpi_comm_parent)')
         call elpa_obj%set("process_row", hmat%blacsdata%myrow, err)
         call check_elpa_err(err, 'set(process_row)')
         call elpa_obj%set("process_col", hmat%blacsdata%mycol, err)
         call check_elpa_err(err, 'set(process_col)')
         call elpa_obj%set("blacs_context", hmat%blacsdata%blacs_desc(2), err)
         call check_elpa_err(err, 'set(blacs_context)')
      type is (t_mat)
         call judft_bug("Elpa solver not available for non-distributed matrices")
         call elpa_obj%set("na", hmat%matsize1, err)
         call check_elpa_err(err, 'set(na)')
         call elpa_obj%set("nblk", hmat%matsize1, err)
         call check_elpa_err(err, 'set(nblk)')
         call elpa_obj%set("mpi_comm_parent", MPI_COMM_SELF, err)
         call check_elpa_err(err, 'set(mpi_comm_parent)')
         call elpa_obj%set("process_row", 0, err)
         call check_elpa_err(err, 'set(process_row)')
         call elpa_obj%set("process_col", 0, err)
         call check_elpa_err(err, 'set(process_col)')
         !Generate a blacs context for this PE only
         call tmp%init(.true.,1,1,MPI_COMM_SELF,MPIMAT_2D_BLOCK_CYCLIC)
         call elpa_obj%set("blacs_context", tmp%blacsdata%blacs_desc(2), err)
         call check_elpa_err(err, 'set(blacs_context)')
      end select
        ! Set the number of eigenvalues
      call elpa_obj%set("nev", ne, err)
      call check_elpa_err(err, 'set(nev)')

      err = elpa_obj%setup()
      call check_elpa_err(err, 'setup')

#if defined(CPP_GPU)||defined(_OPENACC)
      call elpa_obj%set("nvidia-gpu", 1, err)
      call check_elpa_err(err, 'set(nvidia-gpu)')
      !call elpa_obj%set("verbose",1,err)
      err=elpa_obj%setup_gpu()
      call check_elpa_err(err, 'setup_gpu')
      !print *,"ELPA-GPU-err:",err
      call elpa_obj%set("solver", ELPA_SOLVER_1STAGE, err)
      call check_elpa_err(err, 'set(solver)')
#else
      call elpa_obj%set("solver", ELPA_SOLVER_2STAGE, err)
      call check_elpa_err(err, 'set(solver)')
#endif
      
    
   call timestop("ELPA SETUP")
#endif
   end subroutine

   subroutine elpa_gev(self, hmat, smat, ne, eig, zmat, ikpt)
      use m_types_mat
      use m_judft

      implicit none
      class(t_solver_elpa)                  :: self
      class(t_mat), intent(INOUT) :: hmat, smat
      integer, intent(INOUT) :: ne
      class(t_mat), allocatable, intent(OUT)   :: zmat
      real, intent(OUT)   :: eig(:)
      integer, intent(IN) :: ikpt
#ifdef CPP_ELPA
      integer           :: num, np, myid,kernel
      integer           :: err
      integer           :: i
      real, allocatable      :: eig2(:)
      class(t_mat),allocatable        :: ev_dist
      !Update elpa object
      call create_elpa_obj(hmat, ne)
      allocate(ev_dist,mold=hmat)

      call timestart("ELPA GEV")
      select type(hmat)
         type is (t_mpimat)  !we need some data from t_mpimat
            allocate (eig2(hmat%global_size1), stat=err) ! The eigenvalue array
         type is (t_mat)
            allocate (eig2(hmat%matsize1), stat=err) ! The eigenvalue array
      end select         
      if (err .ne. 0) call juDFT_error('Failed to allocated "eig2"', calledby='elpa')

      call ev_dist%init(hmat)! Eigenvectors
      if (err .ne. 0) call juDFT_error('Failed to allocated "ev_dist"', calledby='elpa')

      call hmat%u2l()
      call smat%u2l()
      call elpa_obj%timer_start("ELPA")
      if (hmat%l_real) then
         call elpa_obj%generalized_eigenvectors(hmat%data_r, smat%data_r, eig2, ev_dist%data_r, .false., err)
         call check_elpa_err(err, 'generalized_eigenvectors(real)')
      else
         call elpa_obj%generalized_eigenvectors(hmat%data_c, smat%data_c, eig2, ev_dist%data_c, .false., err)
         call check_elpa_err(err, 'generalized_eigenvectors(complex)')
      end if
      call elpa_obj%timer_stop("ELPA")
      !Useful for debugging
      !call mpi_comm_rank(MPI_COMM_WORLD,myid,err)
      !if (myid == 0) then
      !   call elpa_obj%get("complex_kernel", kernel)
      !   print *, "elpa uses "//elpa_int_value_to_string("complex_kernel", kernel)//" kernel"
      !   call elpa_obj%print_times("ELPA")
      !endif

      
      eig(:ne) = eig2(:ne)
      deallocate (eig2)
      
      select type(hmat)
         type is (t_mpimat)
            !
            !     Redistribute eigenvectors  from ScaLAPACK distribution to each process, i.e. for
            !     process i these are eigenvectors i+1, np+i+1, 2*np+i+1...
            !     Only num=num2/np eigenvectors per process
            !
         call MPI_COMM_RANK(hmat%blacsdata%mpi_com, myid, err)
         call MPI_COMM_SIZE(hmat%blacsdata%mpi_com, np, err)
            
            num = ne
            ne = 0
            do i = myid + 1, num, np
               ne = ne + 1
            end do
            !
            !
            allocate (t_mpimat::zmat)
            call zmat%init(hmat%l_real, hmat%global_size1,hmat%global_size1 , hmat%blacsdata%mpi_com, MPIMAT_ROWCYCLIC)
            call zmat%copy(ev_dist, 1, 1)
         type is (t_mat)
            allocate(t_mat::zmat)
            call zmat%init(hmat%l_real,hmat%matsize1,ne)
            if (zmat%l_real) THEN
               zmat%data_r(:,:)=ev_dist%data_r(:,:ne)
            else      
               zmat%data_c(:,:)=ev_dist%data_c(:,:ne)
            endif
      end select   
      call timestop("ELPA GEV")
      call elpa_deallocate(elpa_obj, err)
      call check_elpa_err(err, 'elpa_deallocate')
      if (associated(elpa_obj)) elpa_obj=>null()
            
#endif

   end subroutine elpa_gev


   subroutine elpa_diag(self, hmat, ne, eig, zmat)
      !Simple driver to solve Generalized Eigenvalue Problem using LAPACK routine
      implicit none
      class(t_solver_elpa)            :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)
#ifdef CPP_ELPA
      real,allocatable:: eig2(:)
      class(t_mat),allocatable :: ev_dist
      integer :: err,myid,num,np,i

      !Update elpa object
      call create_elpa_obj(hmat, ne)
      allocate(ev_dist,mold=hmat)
      
      call timestart("ELPA STD")
      select type(hmat)
         type is (t_mpimat)  !we need some data from t_mpimat
            allocate (eig2(hmat%global_size1), stat=err) ! The eigenvalue array
         type is (t_mat)   
            allocate (eig2(hmat%matsize1), stat=err) ! The eigenvalue array
      end select
      if (err .ne. 0) call juDFT_error('Failed to allocated "eig2"', calledby='elpa')
      call ev_dist%init(hmat)! Eigenvectors
      if (err .ne. 0) call juDFT_error('Failed to allocated "ev_dist"', calledby='elpa')

      call hmat%u2l()
      call elpa_obj%timer_start("ELPA")
      if (hmat%l_real) then
         call elpa_obj%eigenvectors(hmat%data_r, eig2, ev_dist%data_r, err)
         call check_elpa_err(err, 'eigenvectors(real)')
      else
         call elpa_obj%eigenvectors(hmat%data_c, eig2, ev_dist%data_c, err)
         call check_elpa_err(err, 'eigenvectors(complex)')
      end if
      call elpa_obj%timer_stop("ELPA")
      ! END of ELPA stuff
      !
      !     Each process has all eigenvalues in output
      eig(:ne) = eig2(:ne)
      deallocate (eig2)
      select type (hmat)
         type is (t_mpimat)
         !
         !
         !     Redistribute eigenvectors  from ScaLAPACK distribution to each process, i.e. for
         !     process i these are eigenvectors i+1, np+i+1, 2*np+i+1...
         !     Only num=num2/np eigenvectors per process
         !
         call MPI_COMM_RANK(hmat%blacsdata%mpi_com, myid, err)
         call MPI_COMM_SIZE(hmat%blacsdata%mpi_com, np, err)
         
         num = ne
         ne = 0
         do i = myid + 1, num, np
            ne = ne + 1
         end do
         allocate(t_mpimat::zmat)
         call zmat%init(hmat%l_real, hmat%global_size1, hmat%global_size1, hmat%blacsdata%mpi_com, MPIMAT_ROWCYCLIC)
         call zmat%copy(ev_dist, 1, 1)
      type is (t_mat)
         allocate(t_mat::zmat)
         call zmat%init(hmat%l_real,hmat%matsize1,ne)
         if (zmat%l_real) then
            zmat%data_r(:,:)=ev_dist%data_r(:,:ne)
         else
            zmat%data_c(:,:)=ev_dist%data_c(:,:ne)
         end if
      end select
      call timestop("ELPA STD")
      call elpa_deallocate(elpa_obj, err)
      call check_elpa_err(err, 'elpa_deallocate')
      if (associated(elpa_obj)) elpa_obj=>null()
#endif
   end subroutine
   subroutine elpa_diag_sp(self, hmat, ne, eig, zmat)
      !Simple driver to solve Generalized Eigenvalue Problem using LAPACK routine
      implicit none
      class(t_solver_elpa)            :: self
      class(t_mat), intent(INOUT)  :: hmat
      integer, intent(INOUT)      :: ne
      class(t_mat), allocatable, intent(OUT)    :: zmat
      real, intent(OUT)           :: eig(:)


      integer, parameter:: sp = selected_real_kind(6)
      real(kind=sp),allocatable:: eig2(:)
      class(t_mat),allocatable :: ev_dist
      integer :: err,myid,num,np,i
      
#ifdef CPP_ELPA_SP
      !Update elpa object
      call create_elpa_obj(hmat, ne)
      allocate(ev_dist,mold=hmat)
      

      call timestart("ELPA STD-SP")
      select type(hmat)
            type is (t_mpimat)  !we need some data from t_mpimat
            allocate (eig2(hmat%global_size1), stat=err) ! The eigenvalue array
            type is (t_mat)
            allocate (eig2(hmat%matsize1), stat=err) ! The eigenvalue array
      end select
      if (err .ne. 0) call juDFT_error('Failed to allocated "eig2"', calledby='elpa')
      call ev_dist%init(hmat)! Eigenvectors
      if (err .ne. 0) call juDFT_error('Failed to allocated "ev_dist"', calledby='elpa')

      call hmat%u2l()
      call elpa_obj%timer_start("ELPA")
      if (hmat%l_real) then
         block
            real(kind=sp),allocatable:: mat(:,:),z(:,:)
            mat=hmat%data_r
            allocate(z(size(ev_dist%data_r,1),size(ev_dist%data_r,2)))
            call elpa_obj%eigenvectors(mat, eig2, z, err)
            call check_elpa_err(err, 'eigenvectors_sp(real)')
            ev_dist%data_r=z
         end block
      else
         block
            complex(kind=sp),allocatable:: mat(:,:),z(:,:)
            mat=hmat%data_c
            allocate(z(size(ev_dist%data_c,1),size(ev_dist%data_c,2)))
            call elpa_obj%eigenvectors(mat, eig2, z, err)
            call check_elpa_err(err, 'eigenvectors_sp(complex)')
            ev_dist%data_c=z
         end block
      end if
      call elpa_obj%timer_stop("ELPA")
      !     Each process has all eigenvalues in output
      eig(:ne) = eig2(:ne)
      deallocate (eig2)
      
      select type (hmat)
      type is (t_mpimat)
            !
            !     Redistribute eigenvectors  from ScaLAPACK distribution to each process, i.e. for
            !     process i these are eigenvectors i+1, np+i+1, 2*np+i+1...
            !     Only num=num2/np eigenvectors per process
            !
            call MPI_COMM_RANK(hmat%blacsdata%mpi_com, myid, err)
            call MPI_COMM_SIZE(hmat%blacsdata%mpi_com, np, err)
            num = ne
            ne = 0
            do i = myid + 1, num, np
               ne = ne + 1
            end do
            allocate(t_mpimat::zmat)
            call zmat%init(hmat%l_real, hmat%global_size1, hmat%global_size1, hmat%blacsdata%mpi_com, MPIMAT_ROWCYCLIC)
            call zmat%copy(ev_dist, 1, 1)
      type is (t_mat)
            allocate(t_mat::zmat)
            call zmat%init(hmat%l_real,hmat%matsize1,ne)
            if (zmat%l_real) then
               zmat%data_r(:,:)=ev_dist%data_r(:,:ne)
            else
               zmat%data_c(:,:)=ev_dist%data_c(:,:ne)
            end if
      end select

      call timestop("ELPA STD-SP")
      call elpa_deallocate(elpa_obj, err)
      call check_elpa_err(err, 'elpa_deallocate')
      if (associated(elpa_obj)) elpa_obj=>null()
#endif
   end subroutine

   subroutine elpa_to_std(self, hmat, smat, ne)
      !Simple driver to transform Generalized Eigenvalue Problem to Standard problem using LAPACK routine

      class(t_solver_elpa) :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer, intent(IN)  :: ne
      integer            :: err
      real, allocatable :: tmp_real(:,:)
      complex, allocatable :: tmp_cmplx(:,:)

      call timestart("ELPA REDUCTION")
      call create_elpa_obj(hmat, ne)
#ifdef CPP_ELPA
      
#ifndef CPP_ELPA_PATCH
      ! Convert from upper to lower triangular storage for ELPA
      call timestart("ELPA REDUCTION U2L")
      call hmat%u2l()
      call smat%u2l()
      call timestop("ELPA REDUCTION U2L")
      ! Step 1: Cholesky decomposition of S: S = U^T*U, store U in smat
      call timestart("ELPA REDUCTION CHOLESKY")
      if (hmat%l_real) then
         call elpa_obj%cholesky(smat%data_r, err)
         call check_elpa_err(err, 'cholesky(real)')
         ! Step 2: Invert U in-place: smat <- inv(U)
         call elpa_obj%invert_triangular(smat%data_r, err)
         call check_elpa_err(err, 'invert_triangular(real)')
      else
         call elpa_obj%cholesky(smat%data_c, err)
         call check_elpa_err(err, 'cholesky(complex)')
         call elpa_obj%invert_triangular(smat%data_c, err)
         call check_elpa_err(err, 'invert_triangular(complex)')
      end if
      call timestop("ELPA REDUCTION CHOLESKY")
      ! Steps 3+4: A <- inv(U)^T * A * inv(U)
      call timestart("ELPA REDUCTION HERMITIAN_MULTIPLY")
      select type(hmat)
      type is (t_mpimat)
         select type(smat)
         type is (t_mpimat)
            if (hmat%l_real) then
               block
                  real, allocatable :: tmp(:,:)
                  allocate(tmp(hmat%matsize1, hmat%matsize2))
                  ! tmp <- inv(U)^T * A
                  call elpa_obj%hermitian_multiply('U', 'F', hmat%global_size1, &
                       smat%data_r, hmat%data_r, hmat%matsize1, hmat%matsize2, &
                       tmp, hmat%matsize1, hmat%matsize2, err)
                  call check_elpa_err(err, 'hermitian_multiply(real)')
                  hmat%data_r = tmp
               end block
            else
               block
                  complex, allocatable :: tmp(:,:)
                  allocate(tmp(hmat%matsize1, hmat%matsize2))
                  ! tmp <- inv(U)^H * A
                  call elpa_obj%hermitian_multiply('U', 'F', hmat%global_size1, &
                       smat%data_c, hmat%data_c, hmat%matsize1, hmat%matsize2, &
                       tmp, hmat%matsize1, hmat%matsize2, err)
                  call check_elpa_err(err, 'hermitian_multiply(complex)')
                  hmat%data_c = tmp
               end block
            end if
         end select
      end select
      call timestop("ELPA REDUCTION HERMITIAN_MULTIPLY")
      call timestart("ELPA REDUCTION TRMM")
      select type(hmat)
      type is (t_mpimat)
         select type(smat)
         type is (t_mpimat)
            if (hmat%l_real) then
               ! A <- A * inv(U)
               call pdtrmm('R', 'U', 'N', 'N', hmat%global_size1, hmat%global_size1, &
                            1.0d0, smat%data_r, 1, 1, smat%blacsdata%blacs_desc, &
                            hmat%data_r, 1, 1, hmat%blacsdata%blacs_desc)
            else
               ! A <- A * inv(U)
               call pztrmm('R', 'U', 'N', 'N', hmat%global_size1, hmat%global_size1, &
                            (1.0d0,0.0d0), smat%data_c, 1, 1, smat%blacsdata%blacs_desc, &
                            hmat%data_c, 1, 1, hmat%blacsdata%blacs_desc)
            end if
         end select
      end select
      call timestop("ELPA REDUCTION TRMM")
#else
      ! Fallback to old private API when CPP_ELPA_PATCH is not defined
      call timestart("ELPA REDUCTION U2L")
      call hmat%u2l()
      call smat%u2l()
      call timestop("ELPA REDUCTION U2L")
      call timestart("ELPA REDUCTION TRANSFORM_GENERALIZED")
      IF (hmat%l_real) THEN
         allocate(tmp_real(size(hmat%data_r,1),size(hmat%data_r,2)), stat=err)
         call elpa_obj%elpa_transform_generalized_double(hmat%data_r, smat%data_r, tmp_real, .false., err)
         call check_elpa_err(err, 'elpa_transform_generalized_double')
      else
         allocate(tmp_cmplx(size(hmat%data_c,1),size(hmat%data_c,2)), stat=err)
         call elpa_obj%elpa_transform_generalized_double_complex(hmat%data_c, smat%data_c, tmp_cmplx, .false., err)
         call check_elpa_err(err, 'elpa_transform_generalized_double_complex')
      endif
   call timestop("ELPA REDUCTION TRANSFORM_GENERALIZED")
#endif
      call timestop("ELPA REDUCTION")
#endif      
   end subroutine   

   subroutine elpa_recover(self, smat, zmat)

      class(t_solver_elpa) :: self
      class(t_mat), intent(INOUT)  :: zmat, smat
      integer :: error, err, ne

      real, allocatable :: tmp_real(:,:)
      complex, allocatable :: tmp_cmplx(:,:)
      type(t_mpimat):: tmp_mat, dist_zmat

      call timestart("ELPA BACKTRANSFORM")

#ifdef CPP_ELPA

      select type(zmat)
      type is (t_mpimat)
         ne = zmat%global_size2
      type is (t_mat)
         call judft_error("elpa.F90:Back transformation not implemented for non-distributed matrices")
      end select

      ! Back-transform eigenvectors: Q <- inv(U) * Q
#ifndef CPP_ELPA_PATCH
      call timestart("ELPA BACKTRANSFORM TRMM")
      select type(zmat)
      type is (t_mpimat)
         select type(smat)
         type is (t_mpimat)
            call dist_zmat%init(smat)
            call dist_zmat%copy(zmat, 1, 1)
            if (smat%l_real) then
               call pdtrmm('L', 'U', 'N', 'N', smat%global_size1, dist_zmat%global_size2, &
                            1.0d0, smat%data_r, 1, 1, smat%blacsdata%blacs_desc, &
                            dist_zmat%data_r, 1, 1, dist_zmat%blacsdata%blacs_desc)
            else
               call pztrmm('L', 'U', 'N', 'N', smat%global_size1, dist_zmat%global_size2, &
                            (1.0d0,0.0d0), smat%data_c, 1, 1, smat%blacsdata%blacs_desc, &
                            dist_zmat%data_c, 1, 1, dist_zmat%blacsdata%blacs_desc)
            end if
            call zmat%copy(dist_zmat, 1, 1)
         end select
      end select
      call timestop("ELPA BACKTRANSFORM TRMM")
#else
      ! Fallback to old private API
   if (associated(elpa_obj)) then
      call elpa_deallocate(elpa_obj, err)
      call check_elpa_err(err, 'elpa_deallocate (recreate for backtransform)')
      elpa_obj => null()
   end if
   call create_elpa_obj(smat, ne)
   call tmp_mat%init(zmat)
   call tmp_mat%copy(zmat, 1, 1)
      if (smat%l_real) then
         allocate(tmp_real(size(zmat%data_r,1),size(zmat%data_r,2)), stat=err)
         if (err /= 0) call juDFT_error('alloc tmp_real failed in elpa_recover', calledby='elpa')
         call timestart("ELPA BACKTRANSFORM API REAL")
         call elpa_obj%elpa_transform_back_generalized_double(smat%data_r, tmp_mat%data_r, tmp_real, error)
         call timestop("ELPA BACKTRANSFORM API REAL")
         call check_elpa_err(error, 'elpa_transform_back_generalized_double')
      else
         allocate(tmp_cmplx(size(zmat%data_c,1),size(zmat%data_c,2)), stat=err)
         if (err /= 0) call juDFT_error('alloc tmp_cmplx failed in elpa_recover', calledby='elpa')
         call timestart("ELPA BACKTRANSFORM API CMPLX")
         call elpa_obj%elpa_transform_back_generalized_double_complex(smat%data_c, tmp_mat%data_c, tmp_cmplx, error)
         call timestop("ELPA BACKTRANSFORM API CMPLX")
         call check_elpa_err(error, 'elpa_transform_back_generalized_double_complex')
      endif
      call elpa_deallocate(elpa_obj, err)
      call check_elpa_err(err, 'elpa_deallocate')
      if (associated(elpa_obj)) elpa_obj => null()

      ! Copy only the nev columns of the full temporary matrix back to zmat.
      select type(zmat)
      type is (t_mpimat)
      select type(smat)
         type is (t_mpimat)      
         call zmat%free()
         call zmat%init(smat%l_real, smat%global_size1, smat%global_size1, smat%blacsdata%mpi_com, MPIMAT_ROWCYCLIC)
         call zmat%copy(tmp_mat, 1, 1)
      end select
      end select
#endif

      call timestop("ELPA BACKTRANSFORM")

#endif

   end subroutine


end module m_elpa
