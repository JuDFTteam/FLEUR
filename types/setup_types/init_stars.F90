!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_initstars
CONTAINS
  SUBROUTINE init_stars(stars,cell,input,noco,sym,xcpot,atoms,vacuum,oneD,mpi,sphhar)
    USE m_types
    USE m_lapwdim
    USE m_prpqfft
    USE m_prpxcfftmap
    USE m_spg2set
    USE m_stepf
    USE m_strgn
    USE m_judft
    USE m_strgndim
    USE m_prpxcfft
    IMPLICIT NONE
    CLASS(t_stars),INTENT(INOUT):: stars
    TYPE(t_cell),INTENT(INOUT)  :: cell
    TYPE(t_input),INTENT(IN) :: input
    TYPE(t_noco),INTENT(IN)  :: noco
    TYPE(t_sym),INTENT(IN)   :: sym
    CLASS(t_xcpot),INTENT(INOUT) :: xcpot
    TYPE(t_atoms),INTENT(INOUT)  :: atoms
    TYPE(t_vacuum),INTENT(IN)  :: vacuum
    TYPE(t_oneD),INTENT(INOUT)  :: oneD
    TYPE(t_mpi),INTENT(IN)  :: mpi
    TYPE(t_sphhar),INTENT(IN)::sphhar
    
    
     
    CALL lapw_fft_dim(cell,input,noco,stars)
    ! Dimensioning of stars

    !IF (input%film.OR.(sym%namgrp.NE.'any ')) THEN
    IF (input%film) THEN
       CALL strgn1_dim(stars%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
            sym%tau,sym%nop,sym%nop2,stars%mx1,stars%mx2,stars%mx3,&
            stars%ng3,stars%ng2,oneD%odd)
    ELSE
       CALL strgn2_dim(stars%gmax,cell%bmat,sym%invs,sym%zrfs,sym%mrot,&
            sym%tau,sym%nop,stars%mx1,stars%mx2,stars%mx3,&
            stars%ng3,stars%ng2)
    END IF
    stars%kimax2= (2*stars%mx1+1)* (2*stars%mx2+1)-1
    stars%kimax = (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1)-1

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
    IF (xcpot%is_gga()) THEN
       ALLOCATE (stars%ft2_gfx(0:stars%kimax2),stars%ft2_gfy(0:stars%kimax2))
    ELSE
       ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
    END IF

    ! Generate stars
    !IF (input%film.OR.(sym%namgrp.NE.'any ')) THEN
    IF (input%film) THEN
       CALL strgn1(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
       !IF (oneD%odd%d1) THEN
       !   CALL od_strgn1(xcpot,cell,sym,oneD)
       !END IF
    ELSE
       CALL strgn2(stars,sym,atoms,vacuum,sphhar,input,cell,xcpot)
    END IF

    CALL prp_qfft(stars,cell,noco,input)

    CALL prp_xcfft(stars,input,cell,xcpot)

    CALL stepf(sym,stars,atoms,oneD,input,cell,vacuum,mpi)

  END SUBROUTINE init_stars

     
END MODULE m_initstars
