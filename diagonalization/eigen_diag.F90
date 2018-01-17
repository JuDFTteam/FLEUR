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
    PRIVATE
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
  INTEGER,PARAMETER:: diag_debugout=99
  PUBLIC eigen_diag,parallel_solver_available
CONTAINS

  LOGICAL FUNCTION parallel_solver_available()
    parallel_solver_available=any((/diag_elpa,diag_elemental,diag_scalapack/)>0)
  END FUNCTION parallel_solver_available

  SUBROUTINE eigen_diag(jsp,eig_id,it,atoms,DIMENSION,mpi, n_rank,n_size,ne,nk,lapw,input,nred,sub_comm,&
       sym,l_zref, noco,cell,bkpt,el,jij,l_wu,oneD,td,ud, eig,ne_found,smat,hmat,zMat)
    USE m_lapack_diag
    USE m_aline
    USE m_alinemuff
    USE m_types
#ifdef CPP_MAGMA
    USE m_magma
#endif
#ifdef CPP_ELPA
    USE m_elpa
#endif
#ifdef CPP_SCALAPACK
    USE m_chani
#endif
#ifdef CPP_ELEMENTAL
    USE m_elemental
#endif
    IMPLICIT NONE
#ifdef CPP_MPI    
    include 'mpif.h'
#endif
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_dimension),INTENT(IN)  :: dimension
    TYPE(t_oneD),INTENT(IN)       :: oneD
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_jij),INTENT(IN)        :: jij
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(INOUT)    :: lapw !might be modified in aline
    TYPE(t_lapwmat),INTENT(INOUT) :: hmat,smat
    TYPE(t_zMat),INTENT(INOUT)    :: zMat
    INTEGER, INTENT(IN) :: jsp,eig_id,it 
    INTEGER, INTENT(IN) :: n_rank,n_size  ,nk   ,nred,sub_comm
    !INTEGER, INTENT(IN) :: matind(dimension%nbasfcn,2),kveclo(atoms%nlotot)
    INTEGER,INTENT(IN)  :: ne
    INTEGER,INTENT(OUT) :: ne_found
    REAL,INTENT(IN)     :: el(:,:,:)
    LOGICAL, INTENT(IN) :: l_wu,l_zref
    REAL,INTENT(INOUT)  :: bkpt(3)
    TYPE(t_tlmplm),INTENT(IN) :: td
    TYPE(t_usdus),INTENT(IN)  :: ud

    REAL,INTENT(OUT) :: eig(:)

    !Locals
    INTEGER :: ndim,err,n,nn,i,ndim1
    LOGICAL :: parallel
    CHARACTER(len=20)::f
    LOGICAL :: l_real

    l_real=smat%l_real


    !
    !  The array z will contain the eigenvectors and is allocated now
    !
    IF (n_size.NE.1) THEN
       ndim = CEILING(real(dimension%neigd)/n_size)
       ndim1=lapw%nmat
     ELSE
       ndim = dimension%neigd
       ndim1=dimension%nbasfcn   
    ENDIF
    if (l_real) THEN
       ALLOCATE ( zMat%z_r(ndim1,ndim), STAT = err )
    ELSE
       ALLOCATE ( zMat%z_c(ndim1,ndim), STAT = err )
    END if

    zMat%nbasfcn = ndim1
    zMat%nbands = ndim
    zMat%l_real = l_real
    
    IF (err.NE.0) THEN
       WRITE (*,*) 'eigen: error during allocation of the'
       WRITE (*,*) 'eigenvecs',err,'  size: ',dimension%nbasfcn*ndim
       CALL juDFT_error("eigen: Error during allocation of the eigenvecs",calledby ="eigen")
    ENDIF
    if (l_real) THEN
       zMat%z_r=0.0
    else
       zMat%z_c=0.0
    endif
   
    !l_wu selects a full diagonalization step or a direct call of aline with a subspace diagonalization only
    IF (.NOT.l_wu) THEN
       CALL timestart("Diagonalization")
       !Select the solver
       parallel=(n_size>1)
       ne_found=ne
       SELECT CASE (priv_select_solver(parallel))
#ifdef CPP_ELPA
       CASE (diag_elpa)
          CALL MPI_COMM_RANK(sub_comm,n,err)
          IF (smat%l_real) THEN
              CALL elpa_diag(lapw%nmat,SUB_COMM,hmat%data_r,smat%data_r,zMat%z_r,eig,ne_found)
          ELSE
              CALL elpa_diag(lapw%nmat,SUB_COMM,hmat%data_c,smat%data_c,zMat%z_c,eig,ne_found)
          ENDIF
#endif
#ifdef CPP_ELEMENTAL
       CASE (diag_elemental)
          IF (it==1) THEN !switch between direct solver and iterative solver
             CALL ELEMENTAL(lapw%nmat,DIMENSION%nbasfcn/n_size,SUB_COMM,eig,ne_found,1,smat,hmat,zMat)
          ELSE
             CALL ELEMENTAL(lapw%nmat,DIMENSION%nbasfcn/n_size,SUB_COMM,eig,ne_found,0,smat,hmat,zMat)
          ENDIF

#endif
#ifdef CPP_SCALAPACK
       CASE (diag_scalapack)
          CALL chani(lapw%nmat,ndim, n_rank,n_size,SUB_COMM,mpi%mpi_comm,eig,ne_found,smat,hmat,zMat)
#endif
#ifdef CPP_MAGMA
       CASE (diag_magma)
          if (l_real) THEN
             call juDFT_error("REAL diagonalization not implemented in magma.F90")
          else
             print *,"Start magma_diag"
             CALL magma_diag(lapw%nmat,eig,ne_found,a_c=hmat%data_c,b_c=smat%data_c,z_c=zMat%z_c)
          endif
#endif
       CASE (diag_lapack)
          CALL lapack_diag(hmat,smat,eig,zmat)
       CASE (diag_debugout)
          CALL priv_debug_out(smat,hmat)
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
            ud,lapw,td,noco,oneD,bkpt,eig,ne_found,zMat,hmat,smat)
       call timestop("aline")
    ENDIF
    !--->         SECOND VARIATION STEP
    IF (input%secvar) THEN
       !--->           compute and diagonalize the second variation
       !--->           hamiltonian
       call timestart("second variation diagonalization")

       CALL aline_muff(atoms,dimension,sym, cell, jsp,ne_found, ud,td, bkpt,lapw,eig,zMat%z_r,zMat%z_c,l_real)
       call timestop("second variation diagonalization")
    END IF
  END SUBROUTINE eigen_diag
  SUBROUTINE priv_debug_out(smat,hmat)
    USE m_types
    use m_judft
    IMPLICIT NONE
    TYPE(t_lapwmat),INTENT(INOUT) :: hmat,smat
    !small subroutine that does only wite the matrix to a file
    INTEGER:: i,ii
    OPEN(999,file="hmat")
    OPEN(998,file="smat")
    DO i=1,hmat%matsize1
       DO ii=1,i
          IF (hmat%l_real) THEN
             WRITE(999,"(2i6,f15.6)") ii,i,hmat%data_r(ii,i)
             WRITE(998,"(2i6,f15.6)") ii,i,smat%data_r(ii,i)
          ELSE
             WRITE(999,"(2i6,2f15.6)") ii,i,hmat%data_c(ii,i)
             WRITE(998,"(2i6,2f15.6)") ii,i,smat%data_c(ii,i)
          ENDIF
       END DO
    ENDDO
    CLOSE(999)
    CLOSE(998)
    CALL judft_error("STOP in eigen_diag:debug_diag")
  END SUBROUTINE priv_debug_out
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
    IF (juDFT_was_argument("-magma")) diag_solver=diag_magma
    IF (juDFT_was_argument("-debugdiag")) diag_solver=diag_debugout
    
    !Check if solver is possible
    if (diag_solver<0) call priv_solver_error(diag_solver,parallel)
    IF (ANY((/diag_lapack,diag_magma/)==diag_solver).AND.parallel) CALL priv_solver_error(diag_solver,parallel)
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
    CASE (diag_magma)
       CALL juDFT_error("You have not compiled with MAGMA support")
    CASE DEFAULT
       CALL judft_error("You have selected an unkown eigensolver")
    END SELECT
  END SUBROUTINE priv_solver_error


END MODULE m_eigen_diag
