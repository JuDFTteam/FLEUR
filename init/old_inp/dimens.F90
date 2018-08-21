!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dimens
  USE m_juDFT
  USE m_utility
  private
  public :: dimens
CONTAINS
  SUBROUTINE dimens(&
       &                  mpi,input,sym,stars,&
       &                  atoms,sphhar,dimension,vacuum,&
       &                  obsolete,kpts,oneD,hybrid)

    USE m_types
    USE m_dimen7
    USE m_firstglance
    USE m_writeOutHeader
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(INOUT) :: mpi
    TYPE(t_input),INTENT(INOUT) :: input
    TYPE(t_sym),INTENT(INOUT) :: sym
    TYPE(t_stars),INTENT(INOUT) :: stars 
    TYPE(t_atoms),INTENT(INOUT) :: atoms
    TYPE(t_sphhar),INTENT(INOUT) :: sphhar
    TYPE(t_dimension),INTENT(INOUT) :: dimension
    TYPE(t_vacuum),INTENT(INOUT) :: vacuum
    TYPE(t_obsolete),INTENT(INOUT) :: obsolete
    TYPE(t_kpts),INTENT(INOUT) :: kpts
    TYPE(t_oneD),INTENT(INOUT) :: oneD
    TYPE(t_hybrid),INTENT(INOUT) :: hybrid
 
    TYPE(t_cell)     :: cell

    LOGICAL l_kpts,l_qpts,l_inpexist,ldum
    INTEGER n1,n2,n3,n4,n5,n6,n7,n8(3),n9,n10(3),i,j
    INTEGER i_vec(33)


#ifdef CPP_MPI
    INCLUDE 'mpif.h'
    INTEGER ierr(3)
#endif
    oneD%odd%d1=.TRUE.
    l_kpts=.TRUE.

    IF (mpi%irank.EQ.0) WRITE (6,*) 'Your parameters: '

#ifdef CPP_MPI
    CALL MPI_BARRIER(mpi%Mpi_comm,ierr)
#endif

201 IF (mpi%irank == 0) THEN
       IF (l_kpts) WRITE (6,*) ' No fl7para-file found, '
       WRITE (6,*) ' invoking dimen7... '
       !call first_glance to generate k-points
       CALL first_glance(n1,n2,n3,n5,n6,input%itmax,l_kpts,l_qpts,ldum,n7,n8,n10)

       CALL dimen7(input,sym,stars,atoms,sphhar,dimension,vacuum,obsolete,kpts,&
                   oneD,hybrid,cell)
    ENDIF
    !     in case of a parallel calculation we have to broadcast
#ifdef CPP_MPI
    i_vec = (/sym%nop,stars%mx1,stars%mx2,stars%mx3,stars%ng3,stars%ng2,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft&
         &     ,atoms%ntype,atoms%nat,atoms%jmtd,sphhar%ntypsd,sphhar%nlhd,sphhar%memd,atoms%lmaxd,dimension%jspd,vacuum%nvacd,dimension%nvd,dimension%nv2d&
         &     ,1,kpts%nkpt,dimension%nstd,dimension%neigd,dimension%msh,dimension%ncvd,vacuum%layerd,atoms%nlod,atoms%llod,input%itmax/)
    CALL MPI_BCAST(i_vec,33,MPI_INTEGER,0,mpi%Mpi_comm,ierr)
    sym%nop=i_vec(1);stars%mx1=i_vec(2);stars%mx2=i_vec(3);stars%mx3=i_vec(4);stars%ng3=i_vec(5)
    stars%ng2 = i_vec(6);stars%kq1_fft=i_vec(7);stars%kq2_fft=i_vec(8);stars%kq3_fft=i_vec(9)
    stars%kxc1_fft = i_vec(10);stars%kxc2_fft = i_vec(11);stars%kxc3_fft = i_vec(12)
    atoms%ntype = i_vec(13);atoms%nat =i_vec(14);atoms%jmtd=i_vec(15);sphhar%ntypsd=i_vec(16)
    sphhar%nlhd = i_vec(17);sphhar%memd=i_vec(18);atoms%lmaxd=i_vec(19);dimension%jspd=i_vec(20)
    vacuum%nvacd=i_vec(21);dimension%nvd=i_vec(22);dimension%nv2d=i_vec(23)
    kpts%nkpt = i_vec(25); dimension%nstd=i_vec(26);dimension%neigd=i_vec(27);dimension%msh=i_vec(28)
    dimension%ncvd=i_vec(29);vacuum%layerd=i_vec(30);atoms%nlod=i_vec(31);atoms%llod=i_vec(32)
    input%itmax=i_vec(33)
    CALL MPI_BCAST(oneD%odd%d1,1,MPI_LOGICAL,0,mpi%Mpi_comm,ierr)
    !      IF (odd%d1) THEN
    i_vec(:7) = (/oneD%odd%mb,oneD%odd%M,oneD%odd%m_cyl,oneD%odd%chi,oneD%odd%rot,oneD%odd%nop&
         &        ,oneD%odd%n2d/)
    CALL MPI_BCAST(i_vec,7,MPI_INTEGER,0,mpi%Mpi_comm,ierr)
    oneD%odd%mb = i_vec(1);oneD%odd%M = i_vec(2);oneD%odd%m_cyl=i_vec(3)
    oneD%odd%chi = i_vec(4);oneD%odd%rot = i_vec(5);oneD%odd%nop=i_vec(6)
    oneD%odd%n2d= i_vec(7)
    !      ELSE
    !         odd%nop = nop
    !      ENDIF
#endif
    dimension%nspd=(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
    vacuum%nmzd = 250
    vacuum%nmzxyd = 100


  END SUBROUTINE dimens

END MODULE m_dimens
