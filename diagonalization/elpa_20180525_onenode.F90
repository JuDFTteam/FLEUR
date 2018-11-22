!-------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_elpa_onenode
CONTAINS
  SUBROUTINE elpa_diag_onenode(hmat,smat,ne,eig,ev)
    !
    !----------------------------------------------------
    ! Sigensystem solver 
    !  Uses the ELPA (1 node version)
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
    ! U.Alekseeva     Nov. 2018
    !----------------------------------------------------
    USE m_juDFT
    USE m_types_mat
    USE m_types
#ifdef CPP_ELPA_ONENODE
    USE elpa
#endif
    IMPLICIT NONE

    CLASS(t_mat),INTENT(INOUT)    :: hmat,smat
    CLASS(t_mat),ALLOCATABLE,INTENT(OUT)::ev
    REAL,INTENT(OUT)              :: eig(:)
    INTEGER,INTENT(INOUT)         :: ne
    
#ifdef CPP_ELPA_ONENODE
    !...  Local variables
    !
    INTEGER           :: err
    REAL,ALLOCATABLE  :: eig2(:)
    TYPE(t_mat)       :: ev_dist
    INTEGER           :: kernel
    CLASS(elpa_t),pointer :: elpa_obj

    err = elpa_init(20180525)
    elpa_obj => elpa_allocate()
       
    ALLOCATE ( eig2(hmat%matsize1), stat=err ) ! The eigenvalue array
    IF (err.NE.0) CALL juDFT_error('Failed to allocated "eig2"', calledby ='elpa')

    CALL ev_dist%init(hmat)! Eigenvectors
    IF (err.NE.0) CALL juDFT_error('Failed to allocated "ev_dist"',calledby ='elpa')
       
    CALL elpa_obj%set("na", hmat%matsize1, err)
    CALL elpa_obj%set("nev", ne, err)
    CALL elpa_obj%set("local_nrows", hmat%matsize1, err)
    CALL elpa_obj%set("local_ncols", hmat%matsize2, err)
    CALL elpa_obj%set("nblk",hmat%matsize1, err)
    CALL elpa_obj%set("blacs_context", -1, err)
    err = elpa_obj%setup()

    CALL hmat%add_transpose(hmat)       
    CALL smat%add_transpose(smat)       

    IF (hmat%l_real) THEN
        CALL elpa_obj%generalized_eigenvectors(hmat%data_r,smat%data_r,eig2, ev_dist%data_r, .FALSE.,err)
    ELSE
        CALL elpa_obj%generalized_eigenvectors(hmat%data_c,smat%data_c,eig2, ev_dist%data_c, .FALSE., err)
    ENDIF
       
    CALL elpa_deallocate(elpa_obj)
    CALL elpa_uninit()
    ! END of ELPA stuff

    eig(1:ne) = eig2(1:ne)
    DEALLOCATE(eig2)
       
    ALLOCATE(t_mat::ev)
    CALL ev%alloc(hmat%l_real,hmat%matsize1,ne)
    CALL ev%copy(ev_dist,1,1)

#endif
 
END SUBROUTINE elpa_diag_onenode
END MODULE m_elpa_onenode
