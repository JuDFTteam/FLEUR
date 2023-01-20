!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_redist_matrix
CONTAINS
  !> Collect Hamiltonian or overlap matrix to final form
  !!
  !! In the collinear case, this routine just copies mat(1,1) into the final matrix.
  !! If the matrices are distributed, the copy includes a redistribution into the block-cylic form needed by
  !! the diagonalization.
  !! In the non-collinear case, the 2x2 array of matrices is combined into the final matrix. Again a redistribution will happen in the parallel case


  SUBROUTINE eigen_redist_matrix(fmpi,lapw,atoms,mat,mat_final,mat_final_templ)
   USE m_types
   USE m_types_mpimat
   IMPLICIT NONE
   TYPE(t_mpi),INTENT(IN)    :: fmpi
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_atoms),INTENT(IN)  :: atoms
    CLASS(t_mat),INTENT(INOUT):: mat(:,:)
    CLASS(t_mat),INTENT(INOUT):: mat_final
    CLASS(t_mat),INTENT(IN),OPTIONAL :: mat_final_templ

    INTEGER:: m

    !determine final matrix size and allocate the final matrix
    m=lapw%nv(1)+atoms%nlotot
    IF (SIZE(mat)>1) m=m+lapw%nv(2)+atoms%nlotot
    IF (.NOT.PRESENT(mat_final_templ)) THEN
       CALL mat_final%init(mat(1,1)%l_real,m,m,fmpi%diag_sub_comm,.TRUE.) !here the .true. creates a block-cyclic scalapack distribution
    ELSE
       CALL mat_final%init(mat_final_templ)
    ENDIF
    !up-up component (or only component in collinear case)
    IF (SIZE(mat)==1) THEN
       CALL mat_final%move(mat(1,1))
       CALL mat(1,1)%free()
       RETURN
    ENDIF

    CALL mat_final%copy(mat(1,1),1,1)
    CALL mat(1,1)%free()

    !down-down component
    CALL mat_final%copy(mat(2,2),lapw%nv(1)+atoms%nlotot+1,lapw%nv(1)+atoms%nlotot+1)
    CALL mat(2,2)%free()

    if (lapw%nv(1).ne.lapw%nv(2).and.atoms%nlotot>0) call priv_copy_lapwLO_part(mat(2,1),mat(1,2),lapw%nv,atoms%nlotot,fmpi)


    !Now collect off-diagonal parts
    IF (fmpi%n_size == 1 ) THEN
       CALL mat(1,2)%add_transpose(mat(2,1))
    ELSE
       CALL mingeselle(mat(2,1),mat(1,2))
    ENDIF
    CALL mat_final%copy(mat(1,2),1,lapw%nv(1)+atoms%nlotot+1)
    CALL mat(1,2)%free()
    CALL mat(2,1)%free()

  END SUBROUTINE eigen_redist_matrix

  subroutine priv_copy_lapwLO_Part(m1,m2,nv,nlotot,fmpi)
   USE m_types
   USE m_types_mpimat
#ifdef CPP_MPI
   use mpi 
#endif   
   implicit none
   CLASS(t_mat),target,INTENT(INOUT):: m1,m2
   integer,intent(in)               :: nv(2),nlotot
   TYPE(t_mpi),INTENT(IN)           :: fmpi
   
   integer                          :: blocksize,nstart,nstop,noff,i,ii,ierr
   class(t_mat),pointer             :: m_to,m_from
   COMPLEX,ALLOCATABLE              :: tmp(:,:)
   
   if (m1%matsize1>m1%matsize2) THEN
      m_to=>m2
      m_from=>m1
   else
      m_to=>m1
      m_from=>m2
   endif   

  
   select type(m_from)
      type is(t_mat)   
      blocksize=abs(nv(1)-nv(2))
      nstart=m_from%matsize1-nlotot+1-blocksize
      nstop=m_from%matsize1-nlotot+1
      noff=m_from%matsize2-nlotot+1
      ! Do a simple copy
         m_to%data_c(noff:,nstart:nstop)=transpose(conjg(m_from%data_c(nstart:nstop,noff:)))
      type is(t_mpimat)
#ifdef CPP_MPI       
         blocksize=abs(nv(1)-nv(2))
         nstart=m_from%global_size1-nlotot+1-blocksize
         nstop=m_from%global_size1-nlotot+1
         noff=m_from%global_size2-nlotot+1
  
         !In parallel case create a matrix containing the block of the matrix
         ALLOCATE(tmp(blocksize+1,nlotot))
         tmp=0
         !Fill it with all data locally available
         DO i=0,nlotot-1
            IF (mod(i+noff-1,fmpi%n_size)==fmpi%irank) THEN
               ii=(i+noff-1)/fmpi%n_size+1
               tmp(:,i+1)=m_from%data_c(nstart:nstop,ii)
            ENDIF   
         enddo
         !send around (+conjgTranspose)
         tmp=conjg(transpose(tmp))
         CALL mpi_allreduce(MPI_IN_PLACE,tmp,size(tmp),MPI_DOUBLE_COMPLEX,mpi_sum,fmpi%mpi_comm,ierr)
         !Select data relevant for local matrix
         DO i=nstart,nstop
            IF (mod(i-1,fmpi%n_size)==fmpi%irank) THEN
               ii=(i-1)/fmpi%n_size+1
               m_to%data_c(noff:,ii)=tmp(:,i-nstart+1)
            ENDIF
         ENDDO
#endif         
   END SELECT      

  end subroutine
END MODULE m_eigen_redist_matrix
