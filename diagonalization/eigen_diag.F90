!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_diag
  USE m_juDFT
! the parameters are set to negative values to indicate that a particular solver is not compiled
#ifdef CPP_ELEMENTAL
    USE m_elemental
#endif
#ifdef CPP_SCALAPACK
    USE m_chani
#endif
#ifdef CPP_ELPA
    USE m_elpa
#endif
    IMPLICIT NONE
#ifdef CPP_ELPA
  INTEGER,PARAMETER:: diag_elpa=1
#else
  INTEGER,PARAMETER:: diag_elpa=-1
#endif
#ifdef CPP_ELEMENTAL
  INTEGER,PARAMETER:: diag_elemental=2
#else
  INTEGER,PARAMETER:: diag_elemental=-2
#endif
#ifdef CPP_SCALAPACK
  INTEGER,PARAMETER:: diag_scalapack=3
#else
  INTEGER,PARAMETER:: diag_scalapack=-3
#endif
#ifdef CPP_MAGMA
  INTEGER,PARAMETER:: diag_magma=6
#else
  INTEGER,PARAMETER:: diag_magma=-6
#endif
  INTEGER,PARAMETER:: diag_lapack=4
  INTEGER,PARAMETER:: diag_lapack2=5
CONTAINS
  SUBROUTINE eigen_diag(jsp,eig_id,it,atoms,dimension,matsize,mpi, n_rank,n_size,ne,nk,lapw,input,nred,sub_comm,&
       sym,matind,kveclo, noco,cell,bkpt,el,jij,l_wu,oneD,td,ud, eig,a,b,z,ne_found)
    USE m_zsymsecloc
    USE m_aline
    USE m_alinemuff
    USE m_types
    USE m_franza
#ifdef CPP_MAGMA
    USE m_magma
#endif
#ifdef CPP_ELPA
    USE m_elpa
#endif
#ifdef CPP_SCALAPACK
    USE m_chani
#endif
#ifdef CPP_elemental
    USE m_elemental
#endif
    IMPLICIT NONE
#ifdef CPP_MPI    
    include 'mpif.h'
#endif
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: dimension
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_jij),INTENT(IN)       :: jij
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_lapw),INTENT(INOUT)   :: lapw !might be modified in aline
    INTEGER, INTENT(IN) :: jsp,eig_id,it,matsize 
    INTEGER, INTENT(IN) :: n_rank,n_size  ,nk   ,nred,sub_comm
    INTEGER, INTENT(IN) :: matind(dimension%nbasfcn,2),kveclo(atoms%nlotot)
    INTEGER,INTENT(IN)  :: ne
    INTEGER,INTENT(OUT) :: ne_found
    REAL,INTENT(IN)     :: el(:,:,:)
    LOGICAL, INTENT(IN) :: l_wu   
    REAL,INTENT(INOUT)  :: bkpt(3)
    TYPE(t_tlmplm),INTENT(IN) :: td
    TYPE(t_usdus),INTENT(IN)  :: ud

    REAL,INTENT(OUT) :: eig(:)

#ifndef CPP_INVERSION
    COMPLEX, ALLOCATABLE,INTENT(INOUT) :: a(:)
    COMPLEX, ALLOCATABLE,INTENT(INOUT) :: b(:)
    COMPLEX, ALLOCATABLE,INTENT(OUT) :: z(:,:)
#else
    REAL, ALLOCATABLE,INTENT(INOUT) :: a(:)
    REAL, ALLOCATABLE,INTENT(INOUT) :: b(:)
    REAL, ALLOCATABLE,INTENT(OUT) :: z(:,:)
#endif

    !Locals
    REAL :: time1
    REAL :: time2
    INTEGER :: ndim,err,n,nn,i
    LOGICAL :: parallel
    CHARACTER(len=20)::f

#if 1==2
    !This is only needed for debugging
    print *,n_rank,lapw%nmat
    write(f,'(a,i0)') "a.",n_rank
    open(99,file=f)
    write(f,'(a,i0)') "b.",n_rank
    open(98,file=f)
    i=0
    DO n=n_rank+1,lapw%nmat,n_size
       DO nn=1,n
          i=i+1
          write(99,'(2(i10,1x),4(f15.4,1x))') n,nn,a(i)
          write(98,'(2(i10,1x),4(f15.4,1x))') n,nn,b(i)
       ENDDO
    ENDDO
    CALL MPI_BARRIER(MPI_COMM_WORLD,err)
    close(99)
     close(98)
    STOP 'DEBUG'
#endif
    
    !
    !  The array z will contain the eigenvectors and is allocated now
    !
    IF (n_size.NE.1) THEN
       ndim = CEILING(real(dimension%neigd)/n_size)
       ALLOCATE ( z(lapw%nmat,ndim), STAT = err )
    ELSE
       ndim = dimension%neigd
       ALLOCATE ( z(dimension%nbasfcn,ndim), STAT = err )
    ENDIF
    IF (err.NE.0) THEN
       WRITE (*,*) 'eigen: error during allocation of the'
       WRITE (*,*) 'eigenvecs',err,'  size: ',dimension%nbasfcn*ndim
       CALL juDFT_error("eigen: Error during allocation of the eigenvecs",calledby ="eigen")
    ENDIF

    !l_wu selects a full diagonalization step or a direct call of aline with a subspace diagonalization only
    IF (.NOT.l_wu) THEN
       CALL timestart("Diagonalization")
       !Select the solver
       parallel=(n_size>1)
       ne_found=ne
       SELECT CASE (priv_select_solver(parallel))
#ifdef CPP_ELPA
       CASE (diag_elpa)
          CALL elpa(lapw%nmat,n,SUB_COMM,a,b,z,eig,ne_found)
#endif
#ifdef CPP_ELEMENTAL
       CASE (diag_elemental)
          IF (it==1) THEN !switch between direct solver and iterative solver
             CALL elemental(lapw%nmat,dimension%nbasfcn/n_size,SUB_COMM,a,b,z,eig,ne_found,1)
          ELSE
             CALL elemental(lapw%nmat,,dimension%nbasfcn/n_size,SUB_COMM,a,b,z,eig,ne_found,0)
          ENDIF

#endif
#ifdef CPP_SCALAPACK
       CASE (diag_scalapack)
          CALL chani(lapw%nmat,dimension%nbasfcn/n_size,ndim, n_rank,n_size,SUB_COMM,mpi%mpi_comm, a,b,z,eig,ne_found)
#endif
#ifdef CPP_MAGMA
       CASE (diag_magma)
          CALL magma_diag(lapw%nmat,a,b,z,eig,ne_found)
#endif
       CASE (diag_lapack2)
          if (noco%l_ss) call juDFT_error("zsymsecloc not tested with noco%l_ss")
          if (input%gw>1) call juDFT_error("zsymsecloc not tested with input%gw>1")
          CALL zsymsecloc(jsp,input,lapw,bkpt,atoms,kveclo, sym,cell, dimension,matsize,ndim,&
                jij,matind,nred, a,b, z,eig,ne_found)
       CASE (diag_lapack)
          CALL franza(dimension%nbasfcn,ndim, lapw%nmat,&
               (sym%l_zref.AND.(atoms%nlotot.EQ.0)), jij%l_j,matind,nred, a,b,input%gw, z,eig,ne_found)
       CASE DEFAULT
          !This should only happen if you select a solver by hand which was not compiled against
          print*, "You selected a diagonalization scheme without compiling for it"
          CALL priv_solver_error(priv_select_solver(parallel),parallel)
       END SELECT
       CALL timestop("Diagonalization")
       !
    ELSE
       call timestart("aline")
       CALL aline(eig_id,nk,atoms,dimension,sym,cell,input,jsp,el,&
            ud,a,b,lapw,td,noco,oneD,bkpt,z,eig,ne_found)
       call timestop("aline")
    ENDIF
    !--->         SECOND VARIATION STEP
    IF (input%secvar) THEN
       !--->           compute and diagonalize the second variation
       !--->           hamiltonian
       call timestart("second variation diagonalization")

       CALL aline_muff(atoms,dimension,sym, cell, jsp,ne_found, ud,td, bkpt,lapw, z,eig)
       call timestop("second variation diagonalization")
    END IF
  END SUBROUTINE eigen_diag

  FUNCTION priv_select_solver(parallel) result(diag_solver)
    LOGICAL,INTENT(IN):: parallel
    INTEGER           :: diag_solver
    diag_solver=-99

    !Determine the default solver
    IF (parallel) THEN
#ifdef CPP_ELPA
       diag_solver=diag_elpa
#else
       diag_solver=diag_scalapack
#endif
    ELSE
#ifdef CPP_MAGMA
       diag_solver=diag_magma
#else
       diag_solver=diag_lapack
#endif
    ENDIF

    !check if a special solver was requested
    IF (juDFT_was_argument("-elpa")) diag_solver=diag_elpa
    IF (juDFT_was_argument("-scalapack")) diag_solver=diag_scalapack
    IF (juDFT_was_argument("-elemental")) diag_solver=diag_elemental
    IF (juDFT_was_argument("-lapack")) diag_solver=diag_lapack
    IF (juDFT_was_argument("-lapack2")) diag_solver=diag_lapack2
    IF (juDFT_was_argument("-magma")) diag_solver=diag_magma
    !Check if solver is possible
    if (diag_solver<0) call priv_solver_error(diag_solver,parallel)
    if (any((/diag_lapack,diag_lapack2/)==diag_solver).and.parallel) call priv_solver_error(diag_solver,parallel)
    if (any((/diag_elpa,diag_elemental,diag_scalapack/)==diag_solver).and..not.parallel) call priv_solver_error(diag_solver,parallel)

  END FUNCTION priv_select_solver

  SUBROUTINE priv_solver_error(diag_solver,parallel)
    IMPLICIT NONE
    INTEGER,INTENT(IN):: diag_solver
    LOGICAL,INTENT(IN)::parallel
    print *,diag_solver,parallel
    SELECT CASE(diag_solver)
    CASE (diag_elpa)
       IF (parallel) THEN
          CALL juDFT_error("You did not compile with the ELPA solver and hence can not use it")
       ELSE
          CALL juDFT_error("The ELPA solver can not be used in serial")
       ENDIF
    CASE (diag_elemental)
       IF (parallel) THEN
          CALL juDFT_error("You did not compile with the ELEMENTAL solver and hence can not use it")
       ELSE
          CALL juDFT_error("The ELEMENTAL solver can not be used in serial")
       ENDIF
    CASE (diag_scalapack)
       IF (parallel) THEN
          CALL juDFT_error("You did not compile with the SCALAPACK solver and hence can not use it")
       ELSE
          CALL juDFT_error("The SCALAPACK solver can not be used in serial")
       ENDIF
    CASE (diag_lapack)
       IF (parallel) THEN
          CALL juDFT_error("The LAPACK solver can not be used in parallel")
       ENDIF
    CASE (diag_lapack2)
       IF (parallel) THEN
          CALL juDFT_error("The LAPACK2 solver can not be used in parallel")
       ENDIF  
    CASE (diag_magma)
       CALL juDFT_error("You have not compiled with MAGMA support")
    CASE DEFAULT
       CALL judft_error("You have selected an unkown eigensolver")
    END SELECT
  END SUBROUTINE priv_solver_error


END MODULE m_eigen_diag
