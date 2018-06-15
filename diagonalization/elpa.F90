!-------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_elpa
CONTAINS
  SUBROUTINE elpa_diag(hmat,smat,ne,eig,ev)
    !
    !----------------------------------------------------
    !- Parallel eigensystem solver - driver routine based on chani; dw'12
    !  Uses the ELPA for the actual diagonalization
    !
    !
    ! hmat ..... Hamiltonian matrix
    ! smat ..... overlap matrix
    ! ne ....... number of ev's searched (and found) on this node
    !            On input, overall number of ev's searched,
    !            On output, local number of ev's found
    ! eig ...... eigenvalues, output
    ! ev ....... eigenvectors, output
    !
    !----------------------------------------------------

#include"cpp_double.h"
    USE m_juDFT
    USE m_types_mpimat
    USE m_types
#ifdef CPP_ELPA  
    USE elpa1
#ifdef CPP_ELPA2
    USE elpa2
#endif
#ifdef CPP_ELPA_201705003    
    USE elpa
#endif
#endif
    IMPLICIT NONE

    CLASS(t_mat),INTENT(INOUT)    :: hmat,smat
    CLASS(t_mat),ALLOCATABLE,INTENT(OUT)::ev
    REAL,INTENT(out)              :: eig(:)
    INTEGER,INTENT(INOUT)         :: ne
    
#ifdef CPP_ELPA  
    INCLUDE 'mpif.h'
    !...  Local variables
    !
    INTEGER           :: num,num2,myrow,mycol
    INTEGER           :: nb, myid, np
    INTEGER           :: n_col, n_row
    LOGICAL           :: ok
    INTEGER           :: err
    INTEGER           :: mpi_comm_rows,mpi_comm_cols
    INTEGER           :: i,k

    ! large distributed Matrices
    REAL,ALLOCATABLE     :: eig2(:)
    TYPE(t_mpimat)       :: ev_dist
    REAL,ALLOCATABLE     :: tmp2_r(:,:)
    COMPLEX,ALLOCATABLE  :: tmp2_c(:,:)
    INTEGER, EXTERNAL    :: numroc, indxl2g  !SCALAPACK functions
#ifdef CPP_ELPA_201705003    
    INTEGER :: kernel
    CLASS(elpa_t),pointer :: elpa_obj

    err = elpa_init(20170403)
    elpa_obj => elpa_allocate()
#endif

    SELECT TYPE(hmat)
    TYPE IS (t_mpimat)
    SELECT TYPE(smat)
    TYPE IS (t_mpimat)

    CALL MPI_BARRIER(hmat%mpi_com,err)    
    CALL MPI_COMM_RANK(hmat%mpi_com,myid,err)
    CALL MPI_COMM_SIZE(hmat%mpi_com,np,err)
    myrow = myid/hmat%npcol
    mycol = myid -(myid/hmat%npcol)*hmat%npcol  


    !Create communicators for ELPA
#if defined (CPP_ELPA_201705003)
    mpi_comm_rows = -1
    mpi_comm_cols = -1
#elif defined (CPP_ELPA_201605004) || defined (CPP_ELPA_201605003)||defined(CPP_ELPA_NEW)
    err=get_elpa_row_col_comms(hmat%mpi_com, myrow, mycol,mpi_comm_rows, mpi_comm_cols)
#else
    CALL get_elpa_row_col_comms(hmat%mpi_com, myrow, mycol,mpi_comm_rows, mpi_comm_cols)
#endif
    !print *,"creating ELPA comms  --  done"

    num2=ne !no of states solved for
    
    ALLOCATE ( eig2(hmat%global_size1), stat=err ) ! The eigenvalue array for ScaLAPACK
    IF (err.NE.0) CALL juDFT_error('Failed to allocated "eig2"', calledby ='elpa')

    CALL ev_dist%init(hmat%l_real,hmat%global_size1,hmat%global_size2,hmat%mpi_com,.TRUE.)! Eigenvectors for ScaLAPACK
    IF (err.NE.0) CALL juDFT_error('Failed to allocated "ev_dist"',calledby ='elpa')

    IF (hmat%l_real) THEN
       ALLOCATE ( tmp2_r(hmat%matsize1,hmat%matsize2), stat=err ) ! tmp_array
    ELSE
       ALLOCATE ( tmp2_c(hmat%matsize1,hmat%matsize2), stat=err ) ! tmp_array
    ENDIF
    IF (err.NE.0) CALL juDFT_error('Failed to allocated "tmp2"', calledby ='elpa')

    smat%blacs_desc=hmat%blacs_desc
    ev_dist%blacs_desc=hmat%blacs_desc
    
    
    nb=hmat%blacs_desc(5)! Blocking factor
    IF (nb.NE.hmat%blacs_desc(6)) CALL judft_error("Different block sizes for rows/columns not supported")

#ifdef CPP_ELPA_201705003    
    CALL elpa_obj%set("na", hmat%global_size1, err)
    CALL elpa_obj%set("nev", num2, err)
    CALL elpa_obj%set("local_nrows", hmat%matsize1, err)
    CALL elpa_obj%set("local_ncols", hmat%matsize2, err)
    CALL elpa_obj%set("nblk", nb, err)
    CALL elpa_obj%set("mpi_comm_parent", hmat%mpi_com, err)
    CALL elpa_obj%set("process_row", myrow, err)
    CALL elpa_obj%set("process_col", mycol, err)
#ifdef CPP_ELPA2
    CALL elpa_obj%set("solver", ELPA_SOLVER_2STAGE)
#else
    CALL elpa_obj%set("solver", ELPA_SOLVER_1STAGE)
#endif
    err = elpa_obj%setup()
    if (myid == 0) then
       call elpa_obj%get("complex_kernel", kernel)
        print *, "elpa uses " // elpa_int_value_to_string("complex_kernel", kernel) // " kernel"
    endif
#endif
    !print *,"Before elpa"

    !ELPA -start here
    ! Solive generalized preblem
    !
    ! 1. Calculate Cholesky factorization of Matrix S = U**T * U
    !    and invert triangular matrix U
    !
    ! Please note: cholesky_complex/invert_trm_complex are not trimmed for speed.
    ! The only reason having them is that the Scalapack counterpart
    ! PDPOTRF very often fails on higher processor numbers for unknown reasons!

#if defined(CPP_ELPA_201705003)
    IF (hmat%l_real) THEN
       CALL elpa_obj%cholesky(smat%data_r, err)
       CALL elpa_obj%invert_triangular(smat%data_r, err)
    ELSE
       CALL elpa_obj%cholesky(smat%data_c, err)
       CALL elpa_obj%invert_triangular(smat%data_c, err)
    ENDIF
#elif defined(CPP_ELPA_201605003) || defined(CPP_ELPA_201605004)
    IF (hmat%l_real) THEN
       ok=cholesky_real(smat%global_size1,smat%data_r,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false.)
       ok=invert_trm_real(smat%global_size1,smat%data_r,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false.)
    ELSE
       ok=cholesky_complex(smat%global_size1,smat%data_c,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false.)
       ok=invert_trm_complex(smat%global_size1,smat%data_c,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false.)
    ENDIF
#elif defined CPP_ELPA_NEW
    IF (hmat%l_real) THEN
       CALL cholesky_real(smat%global_size1,smat%data_r,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false., ok)
       CALL invert_trm_real(smat%global_size1,smat%data_r,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false., ok)
    ELSE
       CALL cholesky_complex(smat%global_size1,smat%data_c,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false.,ok)
       CALL invert_trm_complex(smat%global_size1,smat%data_c,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false.,ok)
    ENDIF
#else
    IF (hmat%l_real) THEN
       CALL cholesky_real(smat%global_size1,smat%data_r,smat%matsize1,nb,mpi_comm_rows,mpi_comm_cols,.false.)
       CALL invert_trm_real(smat%global_size1,smat%data_r,smat%matsize1,nb,mpi_comm_rows,mpi_comm_cols,.false.)
    ELSE
       CALL cholesky_complex(smat%global_size1,smat%data_c,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false.)
       CALL invert_trm_complex(smat%global_size1,smat%data_c,smat%matsize1,nb,smat%matsize2,mpi_comm_rows,mpi_comm_cols,.false.)
    ENDIF
#endif

    ! 2. Calculate U**-T * H * U**-1

    ! 2a. ev_dist = U**-T * H

    ! H is only set in the upper half, solve_evp_real needs a full matrix
    ! Set lower half from upper half

    DO i=1,hmat%matsize2
       ! Get global column corresponding to i and number of local rows up to
       ! and including the diagonal, these are unchanged in H
       n_col = indxl2g(i,     nb, mycol, 0, hmat%npcol)
       n_row = numroc (n_col, nb, myrow, 0, hmat%nprow)
       IF (hmat%l_real) THEN
          hmat%data_r(n_row+1:hmat%matsize1,i) = 0.d0 
       ELSE
          hmat%data_c(n_row+1:hmat%matsize1,i) = 0.d0
       ENDIF
    ENDDO
 

    IF (hmat%l_real) THEN
       CALL pdtran(hmat%global_size1,hmat%global_size1,1.d0,hmat%data_r,1,1,&
                         hmat%blacs_desc,0.d0,ev_dist%data_r,1,1,ev_dist%blacs_desc)
    ELSE
       CALL pztranc(hmat%global_size1,hmat%global_size2,cmplx(1.d0,0.d0),hmat%data_c,1,1,&
                         hmat%blacs_desc,cmplx(0.d0,0.d0),ev_dist%data_c,1,1,ev_dist%blacs_desc)
    ENDIF

    
    DO i=1,hmat%matsize2
       ! Get global column corresponding to i and number of local rows up to
       ! and including the diagonal, these are unchanged in H
       n_col = indxl2g(i,     nb, mycol, 0, hmat%npcol)
       n_row = numroc (n_col, nb, myrow, 0, hmat%nprow)
       IF (hmat%l_real) THEN
          hmat%data_r(n_row+1:hmat%matsize1,i) = ev_dist%data_r(n_row+1:ev_dist%matsize1,i)
       ELSE
          hmat%data_c(n_row+1:hmat%matsize1,i) = ev_dist%data_c(n_row+1:ev_dist%matsize1,i)
       ENDIF
    ENDDO
 

#if defined (CPP_ELPA_201705003)
    IF (hmat%l_real) THEN
       CALL elpa_obj%hermitian_multiply('U','L',hmat%global_size1,smat%data_r,hmat%data_r,&
                     smat%matsize1,smat%matsize2,ev_dist%data_r,hmat%matsize1,hmat%matsize2,err)
    ELSE
       CALL elpa_obj%hermitian_multiply('U','L',hmat%global_size1,smat%data_c,hmat%data_c,&
                     smat%matsize1,smat%matsize2,ev_dist%data_c,hmat%matsize1,hmat%matsize2,err)
    ENDIF
#elif defined (CPP_ELPA_201605004)
    IF (hmat%l_real) THEN
       ok=elpa_mult_at_b_real('U', 'L',hmat%global_size1,hmat%global_size1,smat%data_r,smat%matsize1,smat%matsize2,&
                     hmat%data_r,hmat%matsize1,hmat%matsize2,nb,mpi_comm_rows, mpi_comm_cols,&
                     ev_dist%data_r,ev_dist%matsize1,ev_dist%matsize2)
    ELSE
       ok=mult_ah_b_complex('U', 'L',hmat%global_size1,hmat%global_size1,smat%data_c,smat%matsize1,smat%matsize2,&
                     hmat%data_c,hmat%matsize1,hmat%matsize2,nb,mpi_comm_rows, mpi_comm_cols,&
                     ev_dist%data_c,ev_dist%matsize1,ev_dist%matsize2)
    ENDIF
#elif defined (CPP_ELPA_201605003)
    IF (hmat%l_real) THEN
       ok=mult_at_b_real('U', 'L',hmat%global_size1,hmat%global_size1,smat%data_r,smat%matsize1,&
                     hmat%data_r,hmat%matsize1,nb,mpi_comm_rows, mpi_comm_cols,ev_dist%data_r,ev_dist%matsize1)
    ELSE
       ok=mult_ah_b_complex('U', 'L',hmat%global_size1,hmat%global_size1,smat%data_c,smat%matsize1,&
                     hmat%data_c,hmat%matsize1,nb,mpi_comm_rows, mpi_comm_cols,ev_dist%data_c,ev_dist%matsize1)
    ENDIF
#else
    IF (hmat%l_real) THEN
       CALL mult_at_b_real('U', 'L',hmat%global_size1,hmat%global_size1,smat%data_r,smat%matsize1,&
                     hmat%data_r,hmat%matsize1,nb,mpi_comm_rows, mpi_comm_cols,ev_dist%data_r,ev_dist%matsize1)
    ELSE
       CALL mult_ah_b_complex('U', 'L',hmat%global_size1,hmat%global_size1,smat%data_c,smat%matsize1,&
                     hmat%data_c,hmat%matsize1,nb,mpi_comm_rows, mpi_comm_cols,ev_dist%data_c,ev_dist%matsize1)
    ENDIF
#endif

    ! 2b. tmp2 = eigvec**T
    IF (hmat%l_real) THEN
       CALL pdtran(ev_dist%global_size1,ev_dist%global_size1,1.d0,ev_dist%data_r,1,1,&
                          ev_dist%blacs_desc,0.d0,tmp2_r,1,1,ev_dist%blacs_desc)
    ELSE
       CALL pztranc(ev_dist%global_size1,ev_dist%global_size1,cmplx(1.0,0.0),ev_dist%data_c,1,1,&
                          ev_dist%blacs_desc,cmplx(0.d0,0.d0),tmp2_c,1,1,ev_dist%blacs_desc)
    ENDIF

    ! 2c. A =  U**-T * tmp2 ( = U**-T * Aorig * U**-1 )
#if defined (CPP_ELPA_201705003)
    IF (hmat%l_real) THEN
       CALL elpa_obj%hermitian_multiply('U','U',smat%global_size1,smat%data_r,tmp2_r,&
                     smat%matsize1,smat%matsize2,hmat%data_r,hmat%matsize1,hmat%matsize2,err)
    ELSE
       CALL elpa_obj%hermitian_multiply('U','U',smat%global_size1,smat%data_c,tmp2_c,&
                     smat%matsize1,smat%matsize2,hmat%data_c,hmat%matsize1,hmat%matsize2,err)
    ENDIF
#elif defined (CPP_ELPA_201605004)
    IF (hmat%l_real) THEN
       ok=elpa_mult_at_b_real('U', 'U',smat%global_size1,smat%global_size1,smat%data_r,smat%matsize1,smat%matsize2,&
                     tmp2_r,SIZE(tmp2_r,1),SIZE(tmp2_r,2),nb,mpi_comm_rows, mpi_comm_cols,&
                     hmat%data_r,hmat%matsize1,hmat%matsize2)
    ELSE
       ok=mult_ah_b_complex('U', 'U',smat%global_size1,smat%global_size1,smat%data_c,smat%matsize1,smat%matsize2,&
                     tmp2_c,SIZE(tmp2_c,1),SIZE(tmp2_c,2),nb,mpi_comm_rows, mpi_comm_cols,&
                     hmat%data_c,hmat%matsize1,hmat%matsize2)
    ENDIF
#elif defined (CPP_ELPA_201605003)
    IF (hmat%l_real) THEN
       ok=mult_at_b_real('U', 'U',smat%global_size1,smat%global_size1,smat%data_r,smat%matsize1,&
                     tmp2_r,SIZE(tmp2_r,1),nb,mpi_comm_rows, mpi_comm_cols,hmat%data_r,hmat%matsize1)
    ELSE
       ok=mult_ah_b_complex('U', 'U',smat%global_size1,smat%global_size1,smat%data_c,smat%matsize1,&
                     tmp2_c,SIZE(tmp2_c,1),nb,mpi_comm_rows, mpi_comm_cols,hmat%data_c,hmat%matsize1)
    ENDIF
#else
    IF (hmat%l_real) THEN
       CALL mult_at_b_real('U', 'U',smat%global_size1,smat%global_size1,smat%data_r,smat%matsize1,&
                     tmp2_r,SIZE(tmp2_r,1),nb,mpi_comm_rows, mpi_comm_cols,hmat%data_r,hmat%matsize1)
    ELSE
       CALL mult_ah_b_complex('U', 'U',smat%global_size1,smat%global_size1,smat%data_c,smat%matsize1,&
                     tmp2_c,SIZE(tmp2_c,1),nb,mpi_comm_rows, mpi_comm_cols,hmat%data_c,hmat%matsize1)
    ENDIF
#endif

    ! A is only set in the upper half, solve_evp_real needs a full matrix
    ! Set lower half from upper half

    IF (hmat%l_real) THEN
       CALL pdtran(hmat%global_size1,hmat%global_size1,1.d0,hmat%data_r,1,1,&
                          hmat%blacs_desc,0.d0,ev_dist%data_r,1,1,ev_dist%blacs_desc)
    ELSE
       CALL pztranc(hmat%global_size1,hmat%global_size1,cmplx(1.0,0.0),hmat%data_c,1,1,&
                          hmat%blacs_desc,cmplx(0.d0,0.d0),ev_dist%data_c,1,1,ev_dist%blacs_desc)
    ENDIF


    DO i=1,hmat%matsize2
       ! Get global column corresponding to i and number of local rows up to
       ! and including the diagonal, these are unchanged in A
       n_col = indxl2g(i,     nb, mycol, 0, hmat%npcol)
       n_row = numroc (n_col, nb, myrow, 0, hmat%nprow)
       IF (hmat%l_real) THEN
          hmat%data_r(n_row+1:hmat%matsize1,i) = ev_dist%data_r(n_row+1:ev_dist%matsize1,i)
       ELSE
          hmat%data_c(n_row+1:hmat%matsize1,i) = ev_dist%data_c(n_row+1:ev_dist%matsize1,i)
       ENDIF
    ENDDO

!print*, "matrix bfore diag" 
!open(unit = 99+myid,file="devel_matr_"//achar(myid+48)) 
!    write(99+myid,*) hmat%data_c
!close(unit = 99+myid) 

    ! 3. Calculate eigenvalues/eigenvectors of U**-T * A * U**-1
    !    Eigenvectors go to eigvec
#if defined (CPP_ELPA_201705003)
    IF (hmat%l_real) THEN
       CALL elpa_obj%eigenvectors(hmat%data_r, eig2, ev_dist%data_r, err)
    ELSE
       CALL elpa_obj%eigenvectors(hmat%data_c, eig2, ev_dist%data_c, err)
    ENDIF
#elif defined(CPP_ELPA_201605003) || defined(CPP_ELPA_201605004)
#ifdef CPP_ELPA2
    IF (hmat%l_real) THEN
       ok=solve_evp_real_2stage(hmat%global_size1,num2,hmat%data_r,hmat%matsize1,&
             eig2,ev_dist%data_r,ev_dist%matsize1, nb,ev_dist%matsize2, mpi_comm_rows, mpi_comm_cols,hmat%mpi_com)
    ELSE
       ok=solve_evp_complex_2stage(hmat%global_size1,num2,hmat%data_c,hmat%matsize1,&
             eig2,ev_dist%data_c,ev_dist%matsize1, nb,ev_dist%matsize2, mpi_comm_rows, mpi_comm_cols,hmat%mpi_com)
    ENDIF
#else
    IF (hmat%l_real) THEN
       ok=solve_evp_real_1stage(hmat%global_size1,num2,hmat%data_r,hmat%matsize1,&
             eig2,ev_dist%data_r,ev_dist%matsize1, nb,ev_dist%matsize2, mpi_comm_rows, mpi_comm_cols)
    ELSE
       ok=solve_evp_complex_1stage(hmat%global_size1,num2,hmat%data_c,hmat%matsize1,&
             eig2,ev_dist%data_c,ev_dist%matsize1, nb,ev_dist%matsize2, mpi_comm_rows, mpi_comm_cols)
    ENDIF
#endif
#elif defined CPP_ELPA_NEW
#ifdef CPP_ELPA2
    IF (hmat%l_real) THEN
       err=solve_evp_real_2stage(hmat%global_size1,num2,hmat%data_r,hmat%matsize1,&
             eig2,ev_dist%data_r,ev_dist%matsize1, nb,ev_dist%matsize2, mpi_comm_rows, mpi_comm_cols,hmat%mpi_com)
    ELSE
       err=solve_evp_complex_2stage(hmat%global_size1,num2,hmat%data_c,hmat%matsize1,&
             eig2,ev_dist%data_c,ev_dist%matsize1, nb,ev_dist%matsize2, mpi_comm_rows, mpi_comm_cols,hmat%mpi_com)
    ENDIF
#else
    IF (hmat%l_real) THEN
       err=solve_evp_real_1stage(hmat%global_size1,num2,hmat%data_r,hmat%matsize1,&
             eig2,ev_dist%data_r,ev_dist%matsize1, nb,ev_dist%matsize2, mpi_comm_rows, mpi_comm_cols)
    ELSE
       err=solve_evp_complex_1stage(hmat%global_size1,num2,hmat%data_c,hmat%matsize1,&
             eig2,ev_dist%data_c,ev_dist%matsize1, nb,ev_dist%matsize2, mpi_comm_rows, mpi_comm_cols)
    ENDIF
#endif
#else
#ifdef CPP_ELPA2
    IF (hmat%l_real) THEN
       CALL solve_evp_real_2stage(hmat%global_size1,num2,hmat%data_r,hmat%matsize1,&
             eig2,ev_dist%data_r,ev_dist%matsize1, nb, mpi_comm_rows, mpi_comm_cols,hmat%mpi_com)
    ELSE
       CALL solve_evp_complex_2stage(hmat%global_size1,num2,hmat%data_c,hmat%matsize1,&
             eig2,ev_dist%data_c,ev_dist%matsize1, nb, mpi_comm_rows, mpi_comm_cols,hmat%mpi_com)
    ENDIF
#else
    IF (hmat%l_real) THEN
       CALL solve_evp_real(hmat%global_size1,num2,hmat%data_r,hmat%matsize1,&
             eig2,ev_dist%data_r,ev_dist%matsize1, nb,mpi_comm_rows, mpi_comm_cols)
    ELSE
       CALL solve_evp_complex(hmat%global_size1,num2,hmat%data_c,hmat%matsize1,&
             eig2,ev_dist%data_c,ev_dist%matsize1, nb,mpi_comm_rows, mpi_comm_cols)
    ENDIF
#endif
#endif

!print*, "eigenvalues" 
!open(unit = 99+myid,file="devel_eigVal_"//achar(myid+48)) 
!do i = 1, size(eig2)
!    write(99+myid,*)i, eig2(i)
!enddo
!close(unit = 99+myid) 


    ! 4. Backtransform eigenvectors: Z = U**-1 * eigvec

    ! mult_ah_b_complex needs the transpose of U**-1, thus tmp2 = (U**-1)**T
    IF (hmat%l_real) THEN
       CALL pdtran(smat%global_size1,smat%global_size1,1.d0,smat%data_r,1,1,&
                         smat%blacs_desc,0.d0,tmp2_r,1,1,smat%blacs_desc)
    ELSE
       CALL pztranc(smat%global_size1,smat%global_size1,cmplx(1.d0,0.d0),smat%data_c,1,1,&
                         smat%blacs_desc,cmplx(0.d0,0.d0),tmp2_c,1,1,smat%blacs_desc)
    ENDIF

#if defined (CPP_ELPA_201705003)
    IF (hmat%l_real) THEN
       CALL elpa_obj%hermitian_multiply('L','N',num2,tmp2_r,ev_dist%data_r,&
                     ev_dist%matsize1,ev_dist%matsize2,hmat%data_r,hmat%matsize1,hmat%matsize2,err)
    ELSE
       CALL elpa_obj%hermitian_multiply('L','N',num2,tmp2_c,ev_dist%data_c,&
                     ev_dist%matsize1,ev_dist%matsize2,hmat%data_c,hmat%matsize1,hmat%matsize2,err)
    ENDIF
#elif defined (CPP_ELPA_201605004)
    IF (hmat%l_real) THEN
       ok=elpa_mult_at_b_real('L', 'N',hmat%global_size1,num2,tmp2_r,hmat%matsize1,hmat%matsize2,&
                     ev_dist%data_r,ev_dist%matsize1,ev_dist%matsize2,nb,mpi_comm_rows, mpi_comm_cols,&
                     hmat%data_r,hmat%matsize1,hmat%matsize2)
    ELSE
       ok=mult_ah_b_complex('L', 'N',hmat%global_size1,num2,tmp2_c,hmat%matsize1,hmat%matsize2,&
                     ev_dist%data_c,ev_dist%matsize1,ev_dist%matsize2,nb,mpi_comm_rows, mpi_comm_cols,&
                     hmat%data_c,hmat%matsize1,hmat%matsize2)
    ENDIF
#elif defined (CPP_ELPA_201605003)
    IF (hmat%l_real) THEN
       ok=elpa_mult_at_b_real('L', 'N',hmat%global_size1,num2,tmp2_r,hmat%matsize1,&
                     ev_dist%data_r,ev_dist%matsize1,nb,mpi_comm_rows, mpi_comm_cols,&
                     hmat%data_r,hmat%matsize1)
    ELSE
       ok=mult_ah_b_complex('L', 'N',hmat%global_size1,num2,tmp2_c,hmat%matsize1,&
                     ev_dist%data_c,ev_dist%matsize1,nb,mpi_comm_rows, mpi_comm_cols,&
                     hmat%data_c,hmat%matsize1)
    ENDIF
#else
    IF (hmat%l_real) THEN
       CALL mult_at_b_real('L', 'N',hmat%global_size1,num2,tmp2_r,hmat%matsize1,&
                     ev_dist%data_r,ev_dist%matsize1,nb,mpi_comm_rows, mpi_comm_cols,&
                     hmat%data_r,hmat%matsize1)
    ELSE
       CALL mult_ah_b_complex('L', 'N',hmat%global_size1,num2,tmp2_c,hmat%matsize1,&
                     ev_dist%data_c,ev_dist%matsize1,nb,mpi_comm_rows, mpi_comm_cols,&
                     hmat%data_c,hmat%matsize1)
    ENDIF
#endif

#if defined (CPP_ELPA_201705003)
    CALL elpa_deallocate(elpa_obj)
    CALL elpa_uninit()
#endif
    ! END of ELPA stuff
#if ( !defined (CPP_ELPA_201705003))
    CALL MPI_COMM_FREE(mpi_comm_rows,err)
    CALL MPI_COMM_FREE(mpi_comm_cols,err)
#endif
    !
    !     Put those eigenvalues expected by chani to eig, i.e. for
    !     process i these are eigenvalues i+1, np+i+1, 2*np+i+1...
    !     Only num=num2/np eigenvalues per process
    !
    num=FLOOR(REAL(num2)/np)
    IF (myid.LT.num2-(num2/np)*np) num=num+1
    ne=0
    DO i=myid+1,num2,np
       ne=ne+1
       eig(ne)=eig2(i)
    ENDDO
    DEALLOCATE(eig2)
    !
    !     Redistribute eigvec from ScaLAPACK distribution to each process
    !     having all eigenvectors corresponding to his eigenvalues as above
    !
    ALLOCATE(t_mpimat::ev)
    CALL ev%init(hmat%l_real,hmat%global_size1,hmat%global_size1,hmat%mpi_com,.FALSE.)
    CALL ev%copy(hmat,1,1)
    CLASS DEFAULT
      call judft_error("Wrong type (1) in scalapack")
    END SELECT
    CLASS DEFAULT
      call judft_error("Wrong type (2) in scalapack")
    END SELECT

    IF (hmat%l_real) THEN
        DEALLOCATE (tmp2_r)
    ELSE
        DEALLOCATE (tmp2_c)
    ENDIF

#endif
  END SUBROUTINE elpa_diag 
END MODULE m_elpa
