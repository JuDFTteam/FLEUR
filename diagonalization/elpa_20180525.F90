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
    USE m_juDFT
    USE m_types_mpimat
    USE m_types
    USE elpa
    IMPLICIT NONE

    CLASS(t_mat),INTENT(INOUT)    :: hmat,smat
    CLASS(t_mat),ALLOCATABLE,INTENT(OUT)::ev
    REAL,INTENT(out)              :: eig(:)
    INTEGER,INTENT(INOUT)         :: ne
    
    !...  Local variables
    !
    INTEGER           :: num, np,myid
    INTEGER           :: err
    INTEGER           :: i
    REAL,ALLOCATABLE      :: eig2(:)
    TYPE(t_mpimat)        :: ev_dist
    INTEGER               :: kernel
    CLASS(elpa_t),pointer :: elpa_obj


    SELECT TYPE(hmat)
    TYPE IS (t_mpimat)
    SELECT TYPE(smat)
    TYPE IS (t_mpimat)
       CALL MPI_BARRIER(hmat%blacsdata%mpi_com,err)    
       CALL MPI_COMM_SIZE(hmat%blacsdata%mpi_com,np,err)
       CALL MPI_COMM_RANK(hmat%blacsdata%mpi_com,myid,err)
       err = elpa_init(20180525)
       elpa_obj => elpa_allocate()
       
       ALLOCATE ( eig2(hmat%global_size1), stat=err ) ! The eigenvalue array
       IF (err.NE.0) CALL juDFT_error('Failed to allocated "eig2"', calledby ='elpa')

       CALL ev_dist%init(hmat)! Eigenvectors
       IF (err.NE.0) CALL juDFT_error('Failed to allocated "ev_dist"',calledby ='elpa')
       
       ! Blocking factor
       IF (hmat%blacsdata%blacs_desc(5).NE.hmat%blacsdata%blacs_desc(6)) CALL judft_error("Different block sizes for rows/columns not supported")
       CALL elpa_obj%set("na", hmat%global_size1, err)
       CALL elpa_obj%set("nev", ne, err)
       CALL elpa_obj%set("local_nrows", hmat%matsize1, err)
       CALL elpa_obj%set("local_ncols", hmat%matsize2, err)
       CALL elpa_obj%set("nblk",hmat%blacsdata%blacs_desc(5), err)
       CALL elpa_obj%set("mpi_comm_parent", hmat%blacsdata%mpi_com, err)
       CALL elpa_obj%set("process_row", hmat%blacsdata%myrow, err)
       CALL elpa_obj%set("process_col", hmat%blacsdata%mycol, err)
       CALL elpa_obj%set("blacs_context", hmat%blacsdata%blacs_desc(2), err)
       CALL elpa_obj%set("solver", ELPA_SOLVER_2STAGE)
       err = elpa_obj%setup()

       CALL hmat%generate_full_matrix()
       CALL smat%generate_full_matrix()
       
       IF (hmat%l_real) THEN
          CALL elpa_obj%generalized_eigenvectors(hmat%data_r,smat%data_r,eig2, ev_dist%data_r, .FALSE.,err)
       ELSE
          CALL elpa_obj%generalized_eigenvectors(hmat%data_c,smat%data_c,eig2, ev_dist%data_c, .FALSE., err)
       ENDIF
       
       CALL elpa_deallocate(elpa_obj)
       CALL elpa_uninit()
       ! END of ELPA stuff
       !
       !     Put those eigenvalues expected by chani to eig, i.e. for
       !     process i these are eigenvalues i+1, np+i+1, 2*np+i+1...
       !
       num=ne
       ne=0
       DO i=myid+1,num,np
          ne=ne+1
          eig(ne)=eig2(i)
       ENDDO
       DEALLOCATE(eig2)
       !
       !     Redistribute eigvec from ScaLAPACK distribution to each process
       !     having all eigenvectors corresponding to his eigenvalues as above
       !
       ALLOCATE(t_mpimat::ev)
       CALL ev%init(hmat%l_real,hmat%global_size1,hmat%global_size1,hmat%blacsdata%mpi_com,.FALSE.)
       CALL ev%copy(ev_dist,1,1)
    CLASS DEFAULT
       CALL judft_error("Wrong type (1) in scalapack")
    END SELECT
 CLASS DEFAULT
    CALL judft_error("Wrong type (2) in scalapack")
 END SELECT
 
END SUBROUTINE elpa_diag
END MODULE m_elpa
