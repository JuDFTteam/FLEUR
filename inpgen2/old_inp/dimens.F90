!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_dimens
  USE m_juDFT
  private
  public :: dimens
CONTAINS
  SUBROUTINE dimens(&
       &                  input,sym,stars,&
       &                  atoms,sphhar,dimension,vacuum,&
       &                  kpts,oneD,hybrid)

    USE m_types_input
    USE m_types_sym
    USE m_types_stars
    USE m_types_atoms
    USE m_types_sphhar
    USE m_types_dimension
    USE m_types_vacuum
    USE m_types_kpts
    USE m_types_oned
    USE m_types_hybrid
    USE m_types_cell
    USE m_dimen7
    USE m_firstglance
    IMPLICIT NONE
    TYPE(t_input),INTENT(INOUT) :: input
    TYPE(t_sym),INTENT(INOUT) :: sym
    TYPE(t_stars),INTENT(INOUT) :: stars 
    TYPE(t_atoms),INTENT(INOUT) :: atoms
    TYPE(t_sphhar),INTENT(INOUT) :: sphhar
    TYPE(t_dimension),INTENT(INOUT) :: dimension
    TYPE(t_vacuum),INTENT(INOUT) :: vacuum
    TYPE(t_kpts),INTENT(INOUT) :: kpts
    TYPE(t_oneD),INTENT(INOUT) :: oneD
    TYPE(t_hybrid),INTENT(INOUT) :: hybrid
 
    TYPE(t_cell)     :: cell

    LOGICAL l_kpts,l_qpts,l_inpexist,ldum
    INTEGER n1,n2,n3,n4,n5,n6,n7,n8(3),n9,n10(3),i,j
    INTEGER i_vec(33)


    oneD%odd%d1=.TRUE.
    l_kpts=.TRUE.

       IF (l_kpts) WRITE (6,*) ' No fl7para-file found, '
       WRITE (6,*) ' invoking dimen7... '
       !call first_glance to generate k-points
       CALL first_glance(n1,n2,n3,n5,n6,input%itmax,l_kpts,l_qpts,ldum,n7,n8,n10)

       CALL dimen7(input,sym,stars,atoms,sphhar,dimension,vacuum,kpts,&
                   oneD,hybrid,cell)
    !     in case of a parallel calculation we have to broadcast
    dimension%nspd=(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1)
    vacuum%nmzd = 250
    vacuum%nmzxyd = 100


  END SUBROUTINE dimens

END MODULE m_dimens
