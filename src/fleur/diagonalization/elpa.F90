! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!
! @authors: Miriam Hinzen, Gregor Michalicek
! Added MPI implementation, DW 2018
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
   end type
   class(elpa_t), pointer :: elpa_obj
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
      solver%standard = .false.
      solver%single_precision = .false.
      solver%transform = .false.
      solver%GPU = .true.
      solver%use_sp = .false.
   end function

   subroutine init()
      integer:: err
      if (firstcall) then
         err = elpa_init(20180525)
         firstcall = .false.
         elpa_obj=>null()
      end if
   end subroutine

   subroutine create_elpa_obj(hmat)
      implicit none
      class(t_mat), intent(IN)              :: hmat
      
#ifdef CPP_ELPA
      integer           :: np, myid
      integer           :: err
      
      call init()
      if (associated(elpa_obj)) return
      select type (hmat)
      type IS (t_mpimat)
            call MPI_BARRIER(hmat%blacsdata%mpi_com, err)
            call MPI_COMM_SIZE(hmat%blacsdata%mpi_com, np, err)
            call MPI_COMM_RANK(hmat%blacsdata%mpi_com, myid, err)
            elpa_obj => elpa_allocate()

            ! Blocking factor
            if (hmat%blacsdata%blacs_desc(5) .ne. hmat%blacsdata%blacs_desc(6)) &
               call judft_error("Different block sizes for rows/columns not supported")
            call elpa_obj%set("na", hmat%global_size1, err)
            call elpa_obj%set("local_nrows", hmat%matsize1, err)
            call elpa_obj%set("local_ncols", hmat%matsize2, err)
            call elpa_obj%set("nblk", hmat%blacsdata%blacs_desc(5), err)
            call elpa_obj%set("mpi_comm_parent", hmat%blacsdata%mpi_com, err)
            call elpa_obj%set("process_row", hmat%blacsdata%myrow, err)
            call elpa_obj%set("process_col", hmat%blacsdata%mycol, err)
            call elpa_obj%set("blacs_context", hmat%blacsdata%blacs_desc(2), err)
            call elpa_obj%set("timings", 1, err)
            err = elpa_obj%setup()

#if defined(CPP_GPU)||defined(_OPENACC)
            call elpa_obj%set("gpu_hermitian_multiply", 1, err)
            !call elpa_obj%set("cannon_for_generalized",0,err)
            call elpa_obj%set("nvidia-gpu", 1, err)
            call elpa_obj%setup_gpu()
            !print *, "ELPA for GPU"
            !if (myid .eq. 0) call elpa_obj%store_settings("save_to_disk.txt", err)
#else
            call elpa_obj%set("solver", ELPA_SOLVER_2STAGE)
#endif
      class DEFAULT
         call judft_error("Wrong type (2) in elpa")
      end select
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
      type(t_mpimat)        :: ev_dist
      call init()
      !Update elpa object
      call create_elpa_obj(hmat)
      call elpa_obj%set("nev", ne, err)
      
      call timestart("ELPA GEV")
      select type(hmat)
            type is (t_mpimat)  !we need some data from t_mpimat
            
            allocate (eig2(hmat%global_size1), stat=err) ! The eigenvalue array
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
            !if (myid .eq. 0) call elpa_obj%print_times("ELPA")
            call MPI_BARRIER(hmat%blacsdata%mpi_com, err)
            !CALL elpa_uninit()
            ! END of ELPA stuff
            !
            !     Each process has all eigenvalues in output
            eig(:ne) = eig2(:ne)
            deallocate (eig2)
            !
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
            !
            allocate (t_mpimat::zmat)
            call zmat%init(hmat%l_real, hmat%global_size1,hmat%global_size1 , hmat%blacsdata%mpi_com, MPIMAT_ROWCYCLIC)
            call zmat%copy(ev_dist, 1, 1)
            end select

      call timestop("ELPA GEV")
      call elpa_deallocate(elpa_obj)
      if (associated(elpa_obj)) elpa_obj=>null()
            
#endif

   end subroutine elpa_gev

   subroutine elpa_to_std(self, hmat, smat)
      !Simple driver to transform Generalized Eigenvalue Problem to Standard problem using LAPACK routine

      class(t_solver_elpa) :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer            :: err,n

      class(t_mat),allocatable :: tmp

      allocate(tmp,mold=hmat)
      call tmp%init(hmat)

      select type(hmat)
      type is (t_mpimat)
         n=hmat%global_size1
      end select   

      call init()
      call create_elpa_obj(hmat)
      ! 1. Calculate Cholesky factorization of Matrix S = U**T * U
      !    and invert triangular matrix U.
      ! 2. Calculate U**-T * H * U**-1
            ! 2a. tmp = U**-T * H
       ! 2b. tmp2 = ev_dist**T
      IF (hmat%l_real) THEN
         CALL elpa_obj%cholesky(smat%data_r, err)
         CALL elpa_obj%invert_triangular(smat%data_r, err)
         CALL elpa_obj%hermitian_multiply('U', 'L', n, smat%data_r, hmat%data_r, &
                                          smat%matsize1, smat%matsize2, tmp%data_r, hmat%matsize1, hmat%matsize2, err)
         CALL elpa_obj%hermitian_multiply('U', 'U', n, smat%data_r, tmp%data_r, &
                                          smat%matsize1, smat%matsize2, hmat%data_r, hmat%matsize1, hmat%matsize2, err)
      ELSE
         CALL elpa_obj%cholesky(smat%data_c, err)
         CALL elpa_obj%invert_triangular(smat%data_c, err)
      
         CALL elpa_obj%hermitian_multiply('U', 'L', n, smat%data_c, hmat%data_c, &
                                          smat%matsize1, smat%matsize2, ev_dist%data_c, hmat%matsize1, hmat%matsize2, err)
         CALL elpa_obj%hermitian_multiply('U', 'U', n, smat%data_c, tmp2_c, &
                                          smat%matsize1, smat%matsize2, hmat%data_c, hmat%matsize1, hmat%matsize2, err)
      ENDIF
   end subroutine   
end module m_elpa
