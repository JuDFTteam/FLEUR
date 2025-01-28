! Copyright (c) 2024 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
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

   function get_solver_elpa() result(solver)
      type(t_solver_elpa), pointer::solver
      allocate (solver)
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



   subroutine create_elpa_obj(hmat)
      !$ use omp_lib
      implicit none
      class(t_mat), intent(IN)              :: hmat
      
#ifdef CPP_ELPA
      integer           :: np, myid
      integer           :: err
      TYPE(t_mpimat) :: tmp
      
      if (firstcall) then
         err = elpa_init(20180525)
         firstcall = .false.
         elpa_obj=>null()
      end if
      if (associated(elpa_obj)) return
      elpa_obj => elpa_allocate()

      !Some settings are set for all matrices
      call elpa_obj%set("local_nrows", hmat%matsize1, err)
      call elpa_obj%set("local_ncols", hmat%matsize2, err)
#if defined(CPP_GPU)||defined(_OPENACC)
      call elpa_obj%set("gpu_hermitian_multiply", 1, err)
      !call elpa_obj%set("cannon_for_generalized",0,err)
      call elpa_obj%set("nvidia-gpu", 1, err)
      call elpa_obj%setup_gpu()
#else
      call elpa_obj%set("solver", ELPA_SOLVER_2STAGE)
#endif
      !$ call elpa_obj%set("omp_threads", omp_get_num_threads(),err)
      call elpa_obj%set("timings", 1, err)
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
         call elpa_obj%set("nblk", hmat%blacsdata%blacs_desc(5), err)
         call elpa_obj%set("mpi_comm_parent", hmat%blacsdata%mpi_com, err)
         call elpa_obj%set("process_row", hmat%blacsdata%myrow, err)
         call elpa_obj%set("process_col", hmat%blacsdata%mycol, err)
         call elpa_obj%set("blacs_context", hmat%blacsdata%blacs_desc(2), err)
      type is (t_mat)
         call judft_bug("Elpa solver not available for non-distributed matrices")
         call elpa_obj%set("na", hmat%matsize1, err)
         call elpa_obj%set("nblk", hmat%matsize1, err)
         call elpa_obj%set("mpi_comm_parent", MPI_COMM_SELF, err)
         call elpa_obj%set("process_row", 1, err)
         call elpa_obj%set("process_col", 1, err)
         !Generate a blacs context for this PE only
         call tmp%init(.true.,1,1,MPI_COMM_SELF,MPIMAT_2D_BLOCK_CYCLIC)
         call elpa_obj%set("blacs_context", tmp%blacsdata%blacs_desc(2), err)
      end select
      err = elpa_obj%setup()
#endif
   end subroutine

   subroutine elpa_gev(self, hmat, smat, ne, eig, zmat)
      use m_types_mat
      use m_judft

      implicit none
      class(t_solver_elpa)                  :: self
      class(t_mat), intent(INOUT) :: hmat, smat
      integer, intent(INOUT) :: ne
      class(t_mat), allocatable, intent(OUT)   :: zmat
      real, intent(OUT)   :: eig(:)
#ifdef CPP_ELPA
      integer           :: num, np, myid
      integer           :: err
      integer           :: i
      real, allocatable      :: eig2(:)
      class(t_mat),allocatable        :: ev_dist
      !Update elpa object
      call create_elpa_obj(hmat)
      call elpa_obj%set("nev", ne, err)
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
      else
         call elpa_obj%generalized_eigenvectors(hmat%data_c, smat%data_c, eig2, ev_dist%data_c, .false., err)
      end if
      call elpa_obj%timer_stop("ELPA")
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
      call elpa_deallocate(elpa_obj)
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
      integer :: err,myid,num,np,i

      !Update elpa object
      call create_elpa_obj(hmat)
      call elpa_obj%set("nev", ne, err)
      
      call timestart("ELPA STD")
      select type(hmat)
         type is (t_mpimat)  !we need some data from t_mpimat
            allocate (eig2(hmat%global_size1), stat=err) ! The eigenvalue array
         type is (t_mat)   
            allocate (eig2(hmat%matsize1), stat=err) ! The eigenvalue array
      end select
      if (err .ne. 0) call juDFT_error('Failed to allocated "eig2"', calledby='elpa')

      call zmat%init(hmat)! Eigenvectors
      if (err .ne. 0) call juDFT_error('Failed to allocated "zmat"', calledby='elpa')

      call hmat%u2l()
      call elpa_obj%timer_start("ELPA")
      if (hmat%l_real) then
         call elpa_obj%eigenvectors(hmat%data_r, eig2, zmat%data_r, err)
      else
         call elpa_obj%eigenvectors(hmat%data_c, eig2, zmat%data_c, err)
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
         !
      end select
      call timestop("ELPA STD")
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
      integer :: err,myid,num,np,i
      
#ifdef CPP_ELPA_SP
      !Update elpa object
      call create_elpa_obj(hmat)
      call elpa_obj%set("nev", ne, err)
            

      call timestart("ELPA STD-SP")
      select type(hmat)
            type is (t_mpimat)  !we need some data from t_mpimat
            allocate (eig2(hmat%global_size1), stat=err) ! The eigenvalue array
            type is (t_mat)
            allocate (eig2(hmat%matsize1), stat=err) ! The eigenvalue array
      end select
      if (err .ne. 0) call juDFT_error('Failed to allocated "eig2"', calledby='elpa')

      call zmat%init(hmat)! Eigenvectors
      if (err .ne. 0) call juDFT_error('Failed to allocated "zmat"', calledby='elpa')

      call hmat%u2l()
      call elpa_obj%timer_start("ELPA")
      if (hmat%l_real) then
         block
            real(kind=sp),allocatable:: mat(:,:),z(:,:)
            mat=hmat%data_r
            allocate(z(size(zmat%data_r,1),size(zmat%data_r,2)))
            call elpa_obj%eigenvectors(mat, eig2, z, err)
            zmat%data_r=z
         end block
      else
         block
            complex(kind=sp),allocatable:: mat(:,:),z(:,:)
            mat=hmat%data_c
            allocate(z(size(zmat%data_c,1),size(zmat%data_c,2)))
            call elpa_obj%eigenvectors(mat, eig2, z, err)
            zmat%data_c=z
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
            num = ne
            ne = 0
            do i = myid + 1, num, np
               ne = ne + 1
            end do
            !
      end select

      call timestop("ELPA STD-SP")
#endif
   end subroutine

   subroutine elpa_to_std(self, hmat, smat)
      !Simple driver to transform Generalized Eigenvalue Problem to Standard problem using LAPACK routine

      class(t_solver_elpa) :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer            :: err,n
      logical :: decomposed

      call create_elpa_obj(hmat)
      
      decomposed=.false.
#ifdef CPP_ELPA_PATCH 
      IF (hmat%l_real) THEN
         call elpa_obj%elpa_transform_generalized_d(hmat%data_r,smat%data_r,decomposed,err)
      else   
         call elpa_obj%elpa_transform_generalized_dc(hmat%data_c,smat%data_c,decomposed,err)
      endif
#endif      
      
   end subroutine   

   subroutine elpa_recover(self, smat, zmat)

      class(t_solver_elpa) :: self
      class(t_mat), intent(INOUT)  :: zmat, smat
      integer :: error

      type(t_mat):: tmp_mat
      type(t_mpimat):: tmp_mpimat
#ifdef CPP_ELPA_PATCH
      if (smat%l_real) THEN
         call elpa_obj%elpa_transform_back_generalized_d(smat%data_r, zmat%data_r, error)
      else
         call elpa_obj%elpa_transform_back_generalized_dc(smat%data_c, zmat%data_c, error)
      endif   
      call elpa_deallocate(elpa_obj)
      if (associated(elpa_obj)) elpa_obj=>null()
#endif

      select type(zmat)
      type is (t_mpimat)
         call tmp_mpimat%init(zmat%l_real,zmat%global_size1,zmat%global_size2,zmat%blacsdata%mpi_com,MPIMAT_ROWCYCLIC)
         call tmp_mpimat%copy(zmat,1,1)
         zmat=tmp_mpimat
      type is (t_mat)
         call tmp_mat%init(zmat)
         call tmp_mat%copy(zmat,1,1)
         zmat=tmp_mat
      end select   
      
   end subroutine


end module m_elpa
