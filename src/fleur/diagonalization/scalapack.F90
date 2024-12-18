!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
#ifdef FLEUR_USE_SCOREP
#include 'scorep/SCOREP_User.inc'
#endif
module m_scalapack
   use m_juDFT
   use m_constants
   use m_types_mpimat
   use m_types_mat
   use m_types_solver
#ifdef CPP_MPI
   use mpi
#endif
   implicit none

   type, extends(t_solver)::t_solver_scalapack
   contains
      procedure        :: solve_gev => scalapack_gev  !solver for generalized eigenvalue problem
      procedure        :: to_std => scalapack_reduction     !transform the H of the generalized problem to a std problem
      procedure        :: backtrans => scalapack_recover  !transform the Eigenvalue back to the generalized problem
   end type
   public :: get_solver_scalapack

contains

   function get_solver_scalapack() result(solver)
      type(t_solver_scalapack), pointer::solver
      allocate (solver)
      solver%name = "scalapack"
#ifdef CPP_SCALAPACK
      solver%available = .true.
#else
      solver%available = .false.
#endif
      solver%parallel = .true.
      solver%serial = .false.
      solver%generalized = .true.
      solver%standard = .false.
      solver%single_precision = .false.
      solver%transform = .true.
      solver%GPU = .false.
   end function

   subroutine scalapack_gev(self, hmat, smat, ne, eig, zmat)
      !
      !----------------------------------------------------
      !- Parallel eigensystem solver - driver routine; gb99
      !  Uses the SCALAPACK for the actual diagonalization
      !
      ! hmat ..... Hamiltonian matrix
      ! smat ..... overlap matrix
      ! ne ....... number of ev's searched (and found) on this node
      !            On input, overall number of ev's searched,
      !            On output, local number of ev's found
      ! eig ...... all eigenvalues, output
      ! ev ....... local eigenvectors, output
      !
      !----------------------------------------------------
      !

      implicit none
      class(t_solver_scalapack) :: self
      class(t_mat), intent(INOUT)    :: hmat, smat
      class(t_mat), allocatable, intent(OUT)::zmat
      real, intent(out)              :: eig(:)
      integer, intent(INOUT)         :: ne

#ifdef CPP_SCALAPACK
      !...  Local variables
      !
      integer i, ierr, err
      integer, allocatable :: iwork(:)
      real, allocatable :: rwork(:)
      integer              :: lrwork

      !
      !  ScaLAPACK things
      character(len=1)    :: uplo
      integer              :: num, num1, num2, liwork, lwork2, np0, mq0, np, myid
      integer              :: iceil, numroc, nn, nb
      integer, allocatable :: ifail(:), iclustr(:)
      real                 :: abstol, orfac = 1.e-4, dlamch
      real, allocatable     :: eig2(:), gap(:)
      real, allocatable :: work2_r(:)
      complex, allocatable :: work2_c(:)

      type(t_mpimat):: ev_dist

      external iceil, numroc
      external dlamch

#ifdef FLEUR_USE_SCOREP

      SCOREP_RECORDING_OFF()
#endif
      select type (hmat)
      type IS (t_mpimat)
         select type (smat)
         type IS (t_mpimat)

            allocate (eig2(hmat%global_size1))

            call MPI_COMM_RANK(hmat%blacsdata%mpi_com, myid, ierr)
            call MPI_COMM_SIZE(hmat%blacsdata%mpi_com, np, ierr)

            num = ne !no of states solved for

            abstol = 2.0*dlamch('S') ! PDLAMCH gave an error on ZAMpano

            call ev_dist%init(hmat)

            !smat%blacs_desc(2)    = hmat%blacs_desc(2)
            !ev_dist%blacs_desc(2) = hmat%blacs_desc(2)
            !smat%blacs_desc=hmat%blacs_desc
            !ev_dist%blacs_desc=hmat%blacs_desc

            nb = hmat%blacsdata%blacs_desc(5)! Blocking factor
            if (nb .ne. hmat%blacsdata%blacs_desc(6)) call judft_error("Different block sizes for rows/columns not supported")

            !
            nn = max(max(hmat%global_size1, nb), 2)
            np0 = numroc(nn, nb, 0, 0, hmat%blacsdata%nprow)
            mq0 = numroc(max(max(ne, nb), 2), nb, 0, 0, hmat%blacsdata%npcol)
            if (hmat%l_real) then
               lwork2 = 5*hmat%global_size1 + max(5*nn, np0*mq0 + 2*nb*nb) + &
                        iceil(ne, hmat%blacsdata%nprow*hmat%blacsdata%npcol)*nn
               allocate (work2_r(lwork2 + 10*hmat%global_size1), stat=err) ! Allocate more in case of clusters
            else
               lwork2 = hmat%global_size1 + max(nb*(np0 + 1), 3)
               allocate (work2_c(lwork2), stat=err)
            end if
            if (err .ne. 0) then
               write (*, *) 'work2  :', err, lwork2
               call juDFT_error('Failed to allocated "work2"', calledby='chani')
            end if

            liwork = 6*max(max(hmat%global_size1, hmat%blacsdata%nprow*hmat%blacsdata%npcol + 1), 4)
            allocate (iwork(liwork), stat=err)
            if (err .ne. 0) then
               write (*, *) 'iwork  :', err, liwork
               call juDFT_error('Failed to allocated "iwork"', calledby='chani')
            end if
            allocate (ifail(hmat%global_size1), stat=err)
            if (err .ne. 0) then
               write (*, *) 'ifail  :', err, hmat%global_size1
               call juDFT_error('Failed to allocated "ifail"', calledby='chani')
            end if
            allocate (iclustr(2*hmat%blacsdata%nprow*hmat%blacsdata%npcol), stat=err)
            if (err .ne. 0) then
               write (*, *) 'iclustr:', err, 2*hmat%blacsdata%nprow*hmat%blacsdata%npcol
               call juDFT_error('Failed to allocated "iclustr"', calledby='chani')
            end if
            allocate (gap(hmat%blacsdata%nprow*hmat%blacsdata%npcol), stat=err)
            if (err .ne. 0) then
               write (*, *) 'gap    :', err, hmat%blacsdata%nprow*hmat%blacsdata%npcol
               call juDFT_error('Failed to allocated "gap"', calledby='chani')
            end if
            !
            !     Compute size of workspace
            !
            if (hmat%l_real) then
               uplo = 'U'
               call pdsygvx(1, 'V', 'I', 'U', hmat%global_size1, hmat%data_r, 1, 1, &
                            hmat%blacsdata%blacs_desc, smat%data_r, 1, 1, smat%blacsdata%blacs_desc, &
                            0.0, 1.0, 1, num, abstol, num1, num2, eig2, orfac, ev_dist%data_r, 1, 1, &
                            ev_dist%blacsdata%blacs_desc, work2_r, -1, iwork, -1, ifail, iclustr, gap, ierr)
               if (work2_r(1) .gt. lwork2) then
                  lwork2 = work2_r(1)
                  deallocate (work2_r)
                  allocate (work2_r(lwork2 + 20*hmat%global_size1), stat=err) ! Allocate even more in case of clusters
                  if (err .ne. 0) then
                     write (*, *) 'work2  :', err, lwork2
                     call juDFT_error('Failed to allocated "work2"', calledby='chani')
                  end if
               end if
            else
               lrwork = 4*hmat%global_size1 + max(5*nn, np0*mq0) + &
                        iceil(ne, hmat%blacsdata%nprow*hmat%blacsdata%npcol)*nn
               ! Allocate more in case of clusters
               allocate (rwork(lrwork + 10*hmat%global_size1), stat=ierr)
               if (err /= 0) then
                  write (*, *) 'ERROR: chani.F: Allocating rwork failed'
                  call juDFT_error('Failed to allocated "rwork"', calledby='chani')
               end if

               call pzhegvx(1, 'V', 'I', 'U', hmat%global_size1, hmat%data_c, 1, 1, &
                            hmat%blacsdata%blacs_desc, smat%data_c, 1, 1, smat%blacsdata%blacs_desc, &
                            0.0, 1.0, 1, num, abstol, num1, num2, eig2, orfac, ev_dist%data_c, 1, 1, &
                            ev_dist%blacsdata%blacs_desc, work2_c, -1, rwork, -1, iwork, -1, ifail, iclustr, &
                            gap, ierr)
               if (abs(work2_c(1)) .gt. lwork2) then
                  lwork2 = work2_c(1)
                  deallocate (work2_c)
                  allocate (work2_c(lwork2), stat=err)
                  if (err /= 0) then
                     write (*, *) 'ERROR: chani.F: Allocating rwork failed:', lwork2
                     call juDFT_error('Failed to allocated "work2"', calledby='chani')
                  end if
               end if
               if (rwork(1) .gt. lrwork) then
                  lrwork = rwork(1)
                  deallocate (rwork)
                  ! Allocate even more in case of clusters
                  allocate (rwork(lrwork + 20*hmat%global_size1), stat=err)
                  if (err /= 0) then
                     write (*, *) 'ERROR: chani.F: Allocating rwork failed: ', lrwork + 20*hmat%global_size1
                     call juDFT_error('Failed to allocated "rwork"', calledby='chani')
                  end if
               end if
            end if
            if (iwork(1) .gt. liwork) then
               liwork = iwork(1)
               deallocate (iwork)
               allocate (iwork(liwork), stat=err)
               if (err /= 0) then
                  write (*, *) 'ERROR: chani.F: Allocating iwork failed: ', liwork
                  call juDFT_error('Failed to allocated "iwork"', calledby='chani')
               end if
            end if
            !
            !     Now solve generalized eigenvalue problem
            !
            call timestart("SCALAPACK call")
            if (hmat%l_real) then
               call pdsygvx(1, 'V', 'I', 'U', hmat%global_size1, hmat%data_r, 1, 1, &
                            hmat%blacsdata%blacs_desc, smat%data_r, 1, 1, smat%blacsdata%blacs_desc, &
                            1.0, 1.0, 1, num, abstol, num1, num2, eig2, orfac, ev_dist%data_r, 1, 1, &
                            ev_dist%blacsdata%blacs_desc, work2_r, lwork2, iwork, liwork, ifail, iclustr, &
                            gap, ierr)
            else
               call pzhegvx(1, 'V', 'I', 'U', hmat%global_size1, hmat%data_c, 1, 1, &
                            hmat%blacsdata%blacs_desc, smat%data_c, 1, 1, smat%blacsdata%blacs_desc, &
                            1.0, 1.0, 1, num, abstol, num1, num2, eig2, orfac, ev_dist%data_c, 1, 1, &
                            ev_dist%blacsdata%blacs_desc, work2_c, lwork2, rwork, lrwork, iwork, liwork, &
                            ifail, iclustr, gap, ierr)
               deallocate (rwork)
            end if
            call timestop("SCALAPACK call")
            if (ierr .ne. 0) then
               !IF (ierr /= 2) WRITE (oUnit,*) myid,' error in pzhegvx/pdsygvx, ierr=',ierr
               !IF (ierr <= 0) WRITE (oUnit,*) myid,' illegal input argument'
               if (mod(ierr, 2) /= 0) then
                  !WRITE (oUnit,*) myid,'some eigenvectors failed to converge'
                  eigs: do i = 1, ne
                     if (ifail(i) /= 0) then
                        !WRITE (oUnit,*) myid,' eigenvector',ifail(i), 'failed to converge'
                     else
                        exit eigs
                     end if
                  end do eigs
                  !CALL CPP_flush(oUnit)
               end if
               if (mod(ierr/4, 2) .ne. 0) then
                  !WRITE(oUnit,*) myid,' only',num2,' eigenvectors converged'
                  !CALL CPP_flush(oUnit)
               end if
               if (mod(ierr/8, 2) .ne. 0) then
                  !WRITE(oUnit,*) myid,' PDSTEBZ failed to compute eigenvalues'
                  call judft_warn("SCALAPACK failed to solve eigenvalue problem", calledby="scalapack.f90")
               end if
               if (mod(ierr/16, 2) .ne. 0) then
                  !WRITE(oUnit,*) myid,' B was not positive definite, Cholesky failed at',ifail(1)
                  call judft_warn("SCALAPACK failed: B was not positive definite. "//new_line("a")// &
                                  "Order of the smallest minor which is not positive definite:"//int2str(ifail(1)) &
                                  , calledby="scalapack.f90")
               end if
            end if
            if (num2 < num1) then
               !IF (myid ==0) THEN
               write (oUnit, *) 'Not all eigenvalues wanted are found'
               write (oUnit, *) 'number of eigenvalues/vectors wanted', num1
               write (oUnit, *) 'number of eigenvalues/vectors found', num2
               !CALL CPP_flush(oUnit)
               !ENDIF
            end if
            !
            !     Each process has all eigenvalues in output
            eig(:num2) = eig2(:num2)
            deallocate (eig2)
            !
            !
            !     Redistribute eigenvectors  from ScaLAPACK distribution to each process, i.e. for
            !     process i these are eigenvectors i+1, np+i+1, 2*np+i+1...
            !     Only num=num2/np eigenvectors per process
            !
            num = floor(real(num2)/np)
            if (myid .lt. num2 - (num2/np)*np) num = num + 1
            ne = 0
            do i = myid + 1, num2, np
               ne = ne + 1
               !eig(ne)=eig2(i)
            end do
            allocate (t_mpimat::zmat)
            call zmat%init(ev_dist%l_real, ev_dist%global_size1, ev_dist%global_size1, &
                           ev_dist%blacsdata%mpi_com, MPIMAT_ROWCYCLIC)
            call zmat%copy(ev_dist, 1, 1)
         class DEFAULT
            call judft_error("Wrong type (1) in scalapack")
         end select
      class DEFAULT
         call judft_error("Wrong type (2) in scalapack")
      end select
      call ev_dist%free()
#ifdef FLEUR_USE_SCOREP
      SCOREP_RECORDING_ON()
#endif
#endif
   end subroutine scalapack_gev

   subroutine scalapack_reduction(self, hmat, smat)
      !Simple driver to transform Generalized Eigenvalue Problem to Standard problem using LAPACK routine

      class(t_solver_scalapack) :: self
      class(t_mat), intent(INOUT)  :: hmat, smat
      integer            :: info, n

#ifdef CPP_SCALAPACK
      real :: scale
      select type (hmat)
      type is (t_mpimat)
         select type (smat)
         type is (t_mpimat)
            !Transform to standard problem using SCALAPACK
            if (hmat%l_real) then
               call pdpotrf('U', smat%global_size1, smat%data_r, 1, 1, smat%blacsdata%blacs_desc, info)
               if (info .ne. 0) then
                  write (*, *) 'Error in pdpotrf: info =', info
                  call juDFT_error("1 Reduction failed", calledby="scalapack_reduction")
               end if
               call pdsygst(1, 'U', smat%global_size1, hmat%data_r, 1, 1, hmat%blacsdata%blacs_desc, &
                            smat%data_r, 1, 1, smat%blacsdata%blacs_desc, scale, info)
               if (abs(scale - 1) > 1e-10) call judft_error("Scale parameter not implemented in scalapack_reduction")
               if (info .ne. 0) then
                  write (oUnit, *) 'Error in pdsygst: info =', info
                  call juDFT_error("2 Reduction failed", calledby="scalapack_reduction")
               end if
            else
               call pzpotrf('U', smat%global_size1, smat%data_c, 1, 1, smat%blacsdata%blacs_desc, info)
               if (info .ne. 0) then
                  write (*, *) 'Error in pzpotrf: info =', info
                  call juDFT_error("3 Reduction failed", calledby="scalapack_reduction")
               end if
               call pzhegst(1, 'U', smat%global_size1, hmat%data_c, 1, 1, smat%blacsdata%blacs_desc, &
                            smat%data_c, 1, 1, smat%blacsdata%blacs_desc, scale, info)
               if (abs(scale - 1) > 1e-10) call judft_error("Scale parameter not implemented in scalapack_reduction")
               if (info .ne. 0) then
                  write (oUnit, *) 'Error in pzhegst: info =', info
                  call juDFT_error("4 Reduction failed", calledby="scalapack_reduction")
               end if
            end if
         class default
            call judft_bug("Wrong matrix type in scalapack")
         end select
      class default
         call judft_bug("Wrong matrix type in scalapack")
      end select
#endif
   end subroutine

   subroutine scalapack_recover(self, smat, zmat)

      class(t_solver_scalapack) :: self
      class(t_mat), intent(INOUT)  :: zmat, smat
      integer :: m, n, info
#ifdef CPP_SCALAPACK

      select type (smat)
      type is (t_mpimat)
         select type (zmat)
         type is (t_mpimat)
            n = smat%global_size1
            m = zmat%global_size2
            ! --> recover the generalized eigenvectors z by solving z' = l^t * z
            if (smat%l_real) then
               call pdtrtrs('U', 'N', 'N', n, m, smat%data_r, 1, 1, smat%blacsdata%blacs_desc, &
                            zmat%data_r, 1, 1, zmat%blacsdata%blacs_desc, info)
            else
               call pztrtrs('U', 'N', 'N', n, m, smat%data_c, 1, 1, smat%blacsdata%blacs_desc, &
                            zmat%data_c, 1, 1, zmat%blacsdata%blacs_desc, info)
            end if
            if (info .ne. 0) then
               write (oUnit, *) 'Error in p?trtrs: info =', info
               call juDFT_error("Recovery failed", calledby="scalapack_recover")
            end if
         class default
            call judft_bug("Wrong matrix type in scalapack")
         end select
      class default
         call judft_bug("Wrong matrix type in scalapack")
      end select
#endif
   end subroutine

end module m_scalapack
