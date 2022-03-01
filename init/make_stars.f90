!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_make_stars
  USE m_juDFT
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: make_stars
CONTAINS
  SUBROUTINE make_stars(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot,oneD,noco,fmpi)
    USE m_od_strgn1
    USE m_strgn
    USE m_stepf
    USE m_strgndim
    USE m_lapwdim
    USE m_types_sym
    USE m_types_atoms
    USE m_types_vacuum
    USE m_types_sphhar
    USE m_types_input
    USE m_types_cell
    USE m_types_xcpot
    USE m_types_oned
    USE m_types_mpi
    USE m_types_noco
    USE m_mpi_bc_tool
    USE m_types_stars

    CLASS(t_stars),INTENT(INOUT) :: stars
    TYPE(t_sym),INTENT(in)::sym
    TYPE(t_atoms),INTENT(in)::atoms
    TYPE(t_vacuum),INTENT(in)::vacuum
    TYPE(t_sphhar),INTENT(in)::sphhar
    TYPE(t_input),INTENT(inout)::input
    TYPE(t_cell),INTENT(in)::cell
    CLASS(t_xcpot),INTENT(in)::xcpot
    TYPE(t_oneD),INTENT(inout)::oneD
    TYPE(t_noco),INTENT(in)::noco
    TYPE(t_mpi),INTENT(in)::fmpi
    ! Generate stars
    INTEGER :: kimax,kimax2


    ! Dimensioning of stars
    IF (fmpi%irank==0) THEN
       IF (input%film) THEN
          CALL strgn1_dim(fmpi%irank==0,input%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
               sym%tau,sym%nop,sym%nop2,stars%mx1,stars%mx2,stars%mx3,&
               stars%ng3,stars%ng2,oneD%odd)

       ELSE
          CALL strgn2_dim(fmpi%irank==0,input%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
               sym%tau,sym%nop,stars%mx1,stars%mx2,stars%mx3,&
               stars%ng3,stars%ng2)
          oneD%odd%n2d = stars%ng2
          oneD%odd%nq2 = stars%ng2
          oneD%odd%nop = sym%nop
       END IF
       stars%gmax=input%gmax
       kimax2= (2*stars%mx1+1)* (2*stars%mx2+1)-1
       kimax = (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1)-1
       IF (oneD%odd%d1) THEN
          oneD%odd%k3 = stars%mx3
          oneD%odd%nn2d = (2*(oneD%odd%k3)+1)*(2*(oneD%odd%M)+1)
       ELSE
          oneD%odd%k3 = 0
          oneD%odd%M = 0
          oneD%odd%nn2d = 1
          oneD%odd%mb = 0
       END IF
       ALLOCATE (stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
       ALLOCATE (stars%ig2(stars%ng3))
       ALLOCATE (stars%kv2(2,stars%ng2),stars%kv3(3,stars%ng3))
       ALLOCATE (stars%nstr2(stars%ng2),stars%nstr(stars%ng3))
       ALLOCATE (stars%sk2(stars%ng2),stars%sk3(stars%ng3),stars%phi2(stars%ng2))
       ALLOCATE (stars%igfft(0:kimax,2),stars%igfft2(0:kimax2,2))
       ALLOCATE (stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
       ALLOCATE (stars%pgfft(0:kimax),stars%pgfft2(0:kimax2))
       ALLOCATE (stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1),stars%ustep(stars%ng3))

       stars%kv3(:,:) = 0
       stars%sk2(:) = 0.0
       stars%sk3(:) = 0.0
       stars%phi2(:) = 0.0

       
       ! Missing xc functionals initializations
       IF (xcpot%needs_grad()) THEN
          ALLOCATE (stars%ft2_gfx(0:kimax2),stars%ft2_gfy(0:kimax2))
          ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
               oneD%pgft1xy(0:oneD%odd%nn2d-1),&
               oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
       ELSE
          ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
          ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
               oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
       END IF
       oneD%odd%nq2 = oneD%odd%n2d
       oneD%odi%nq2 = oneD%odd%nq2


       CALL timestart("strgn")
       IF (input%film) THEN
          IF (oneD%odd%d1) THEN
             CALL od_strgn1(xcpot,cell,sym,oneD)
          END IF
          CALL strgn1(fmpi%irank==0,stars,oneD,sym,atoms,vacuum,sphhar,input,cell,xcpot)
       ELSE
          CALL strgn2(fmpi%irank==0,stars,oneD,sym,atoms,vacuum,sphhar,input,cell,xcpot)
       END IF

       CALL lapw_fft_dim(cell,input,noco,stars)

       !count number of stars in 2*rkmax (stars are ordered)
       associate(i=>stars%ng3_fft)
       DO i=stars%ng3,1,-1
         IF ( stars%sk3(i).LE.2.0*input%rkmax ) EXIT
       END DO
       end associate 

       CALL timestop("strgn")
    ENDIF

    CALL stars%mpi_bc(fmpi%mpi_comm)

    CALL timestart("stepf")
    CALL stepf(sym,stars,atoms,oneD,input,cell,vacuum,fmpi)
    CALL mpi_bc(stars%ustep,0,fmpi%mpi_comm)
    CALL mpi_bc(stars%ufft,0,fmpi%mpi_comm)
    CALL timestop("stepf")



  END SUBROUTINE make_stars
END MODULE m_make_stars
