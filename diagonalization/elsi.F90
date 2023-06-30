!-------------------------------------------------------------------------------
! Copyright (c) 2023 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_elsi
CONTAINS
   SUBROUTINE elsi_diag(solver,hmat, smat, ne, eig, ev)
      
      USE m_juDFT
      USE m_types_mpimat
      USE m_types_mat
#ifdef CPP_ELSI
      USE elsi
      USE mpi
#endif
      IMPLICIT NONE

      CLASS(t_mat), INTENT(INOUT)    :: hmat, smat
      CLASS(t_mat), ALLOCATABLE, INTENT(OUT)::ev
      REAL, INTENT(out)              :: eig(:)
      INTEGER, INTENT(INOUT)         :: ne
      INTEGER,INTENT(IN)             :: solver

#ifdef CPP_ELSI
      !...  Local variables
      !
      INTEGER           :: blk,nev,myid,np,ierr,i
      INTEGER,PARAMETER :: BLACS_DENSE=0
      TYPE(elsi_handle) :: eh
      CLASS(t_mat),ALLOCATABLE :: evec
      real,allocatable         :: eig_tmp(:)


      call timestart("ELSI")

      call hmat%u2l()
      call smat%u2l()

      SELECT TYPE (hmat)
      TYPE IS (t_mpimat)
         SELECT TYPE (smat)
         TYPE IS (t_mpimat)
            if (hmat%blacsdata%blacs_desc(5)==hmat%blacsdata%blacs_desc(6)) then
                blk=hmat%blacsdata%blacs_desc(5)
            else
                call judft_error("BUG: in ELSI the row/column blocksize must be equal")
            endif    
            call elsi_init(eh,solver,1,BLACS_DENSE,hmat%global_size1,1.0*ne,ne)
            call elsi_set_mpi(eh,hmat%blacsdata%mpi_com)
            call elsi_set_blacs(eh,hmat%blacsdata%blacs_desc(2),blk)
            ALLOCATE(t_mpimat::evec)
            CALL evec%init(hmat)
            allocate(eig_tmp(hmat%global_size1))
         CLASS DEFAULT
            call judft_error("BUG: Inconsistent matrixes in ELSI call")
         END SELECT
      TYPE IS (t_mat)
          SELECT TYPE(smat)
          TYPE IS (t_mat)
              call elsi_init(eh,solver,0,BLACS_DENSE,hmat%matsize1,1.*ne,ne)
              allocate(t_mat::evec)
              call evec%init(hmat)
              allocate(eig_tmp(hmat%matsize1))

          CLASS DEFAULT
            call judft_error("BUG: Inconsistent matrixes in ELSI call")
         END SELECT  
      END SELECT

      call elsi_set_elpa_gpu(eh,1)
      !Now perform diagonalization
      call elsi_set_output(eh,3)
      call elsi_set_output_unit(eh,7)
      call elsi_set_illcond_check(eh, 0)
      call elsi_reinit(eh)
      IF (hmat%l_real) THEN
        call elsi_ev_real(eh,hmat%data_r,smat%data_r,eig_tmp,evec%data_r)
      ELSE  
        call elsi_ev_complex(eh,hmat%data_c,smat%data_c,eig_tmp,evec%data_c)
      ENDIF  
      
      !Copy data into correct data structures
      eig=eig_tmp(:size(eig))
      SELECT TYPE(evec)
      TYPE IS (t_mat)
         allocate(t_mat::ev)
         call ev%init(evec%l_real,evec%matsize1,ne)
         if (evec%l_real) THEN
            ev%data_r=evec%data_r(:,:ne)
         else   
            ev%data_c=evec%data_c(:,:ne)
         endif
      TYPE IS (t_mpimat)
         ALLOCATE(t_mpimat::ev)
         CALL ev%init(evec%l_real,evec%global_size1,evec%global_size1,evec%blacsdata%mpi_com,.FALSE.)
         CALL ev%copy(evec,1,1)
         !determine ev assigned to this rank
         nev=ne
         ne=0
         CALL MPI_COMM_RANK(evec%blacsdata%mpi_com,myid,ierr)
         CALL MPI_COMM_SIZE(evec%blacsdata%mpi_com,np,ierr)
         DO i=myid+1,nev,np
            ne=ne+1
         !   eig(ne)=eigenvalues(i)
         ENDDO


      END select   
      call elsi_finalize(eh)
      call timestop("ELSI")
     
#endif
   END SUBROUTINE elsi_diag
END MODULE m_elsi
