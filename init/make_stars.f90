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
   SUBROUTINE make_stars(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot,oneD,noco,mpi)
     USE m_od_strgn1
     USE m_strgn
     USE m_stepf
     USE m_prpqfftmap
     USE m_prpqfft
     USE m_strgndim
     use m_lapwdim
     use m_types_sym
     use m_types_atoms
     use m_types_vacuum
     use m_types_sphhar
     use m_types_input
     use m_types_cell
     use m_types_xcpot
     use m_types_oned
     USE m_types_mpi
     use m_types_noco


     class(t_stars),intent(INOUT) :: stars
     type(t_sym),intent(in)::sym
     type(t_atoms),intent(in)::atoms
     type(t_vacuum),intent(in)::vacuum
     type(t_sphhar),intent(in)::sphhar
     type(t_input),intent(inout)::input
     type(t_cell),intent(in)::cell
     class(t_xcpot),intent(in)::xcpot
     TYPE(t_oneD),INTENT(inout)::oneD
     type(t_noco),intent(in)::noco
     type(t_mpi),intent(in)::mpi
     ! Generate stars

    ! Dimensioning of stars
  
  IF (input%film) THEN
     CALL strgn1_dim(input%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
          sym%tau,sym%nop,sym%nop2,stars%mx1,stars%mx2,stars%mx3,&
          stars%ng3,stars%ng2,oneD%odd)
     
  ELSE
     CALL strgn2_dim(input%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
          sym%tau,sym%nop,stars%mx1,stars%mx2,stars%mx3,&
          stars%ng3,stars%ng2)
     oneD%odd%n2d = stars%ng2
     oneD%odd%nq2 = stars%ng2
     oneD%odd%nop = sym%nop
  END IF
  stars%gmax=input%gmax
  stars%kimax2= (2*stars%mx1+1)* (2*stars%mx2+1)-1
  stars%kimax = (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1)-1
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
  ALLOCATE (stars%igfft(0:stars%kimax,2),stars%igfft2(0:stars%kimax2,2))
  ALLOCATE (stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3))
  ALLOCATE (stars%pgfft(0:stars%kimax),stars%pgfft2(0:stars%kimax2))
  ALLOCATE (stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1),stars%ustep(stars%ng3))
  
  stars%sk2(:) = 0.0
  stars%phi2(:) = 0.0

  ! Initialize xc fft box
  
  CALL prp_xcfft_box(xcpot%gmaxxc,cell%bmat,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft)
   
  ! Missing xc functionals initializations
  IF (xcpot%needs_grad()) THEN
     ALLOCATE (stars%ft2_gfx(0:stars%kimax2),stars%ft2_gfy(0:stars%kimax2))
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
     CALL strgn1(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
     IF (oneD%odd%d1) THEN
        CALL od_strgn1(xcpot,cell,sym,oneD)
     END IF
  ELSE
     CALL strgn2(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
  END IF

  CALL lapw_fft_dim(cell,input,noco,stars)

  
  ALLOCATE (stars%igq_fft(0:stars%kq1_fft*stars%kq2_fft*stars%kq3_fft-1))
  ALLOCATE (stars%igq2_fft(0:stars%kq1_fft*stars%kq2_fft-1))
  
  ! Set up pointer for backtransformation from g-vector in positive 
  ! domain of carge density fftibox into stars
  CALL prp_qfft_map(stars,sym,input,stars%igq2_fft,stars%igq_fft)
  CALL prp_qfft(stars,cell,noco,input)
  
  CALL timestop("strgn") 

  CALL timestart("stepf") 
  CALL stepf(sym,stars,atoms,oneD,input,cell,vacuum,mpi)
  CALL timestop("stepf") 
  


END SUBROUTINE make_stars
END MODULE m_make_stars
