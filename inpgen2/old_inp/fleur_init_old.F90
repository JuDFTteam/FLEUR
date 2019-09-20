!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_fleur_init_old
  IMPLICIT NONE
CONTAINS
  !> Collection of code for old-style inp-file treatment
  SUBROUTINE fleur_init_old(&
       input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,&
       sliceplot,banddos,enpara,xcpot,kpts,hybrid,&
       oneD,grid)
    USE m_types_input
    USE m_types_dimension
    USE m_types_atoms
    USE m_types_sphhar
    USE m_types_cell
    USE m_types_stars
    USE m_types_sym
    USE m_types_noco
    USE m_types_vacuum
    USE m_types_sliceplot
    USE m_types_banddos
    USE m_types_enpara
    USE m_types_xcpot_inbuild_nofunction
    USE m_types_kpts
    USE m_types_hybrid
    USE m_types_oned


    USE m_judft
    USE m_dimens
    USE m_inped
    USE m_setup
    USE m_constants

    IMPLICIT NONE
    !     Types, these variables contain a lot of data!
    TYPE(t_input)    ,INTENT(INOUT):: input
    TYPE(t_dimension),INTENT(OUT)  :: DIMENSION
    TYPE(t_atoms)    ,INTENT(OUT)  :: atoms
    TYPE(t_sphhar)   ,INTENT(OUT)  :: sphhar
    TYPE(t_cell)     ,INTENT(OUT)  :: cell
    TYPE(t_stars)    ,INTENT(OUT)  :: stars
    TYPE(t_sym)      ,INTENT(OUT)  :: sym
    TYPE(t_noco)     ,INTENT(OUT)  :: noco
    TYPE(t_vacuum)   ,INTENT(OUT)  :: vacuum
    TYPE(t_sliceplot),INTENT(INOUT):: sliceplot
    TYPE(t_banddos)  ,INTENT(OUT)  :: banddos
    TYPE(t_enpara)   ,INTENT(OUT)  :: enpara
    TYPE(t_xcpot_inbuild_nf),INTENT(OUT)  :: xcpot
    TYPE(t_kpts)     ,INTENT(INOUT):: kpts
    TYPE(t_hybrid)   ,INTENT(OUT)  :: hybrid
    TYPE(t_oneD)     ,INTENT(OUT)  :: oneD
    INTEGER,INTENT(OUT)::grid(3)


    !     .. Local Scalars ..
    INTEGER    :: i,n,l,m1,m2,isym,iisym,pc,iAtom,iType
    COMPLEX    :: cdum
    CHARACTER(len=4)              :: namex
    CHARACTER(len=12)             :: relcor, tempNumberString
    CHARACTER(LEN=20)             :: filename
    REAL                          :: a1(3),a2(3),a3(3)
    REAL                          :: dtild, phi_add
    LOGICAL                       :: l_found, l_kpts, l_exist, l_krla
    character(len=4)              :: latnam,namgrp

    namex = '    '
    relcor = '            '

    CALL dimens(input,sym,stars,atoms,sphhar,DIMENSION,vacuum,&
         kpts,oneD,hybrid)
    stars%kimax2= (2*stars%mx1+1)* (2*stars%mx2+1)-1
    stars%kimax = (2*stars%mx1+1)* (2*stars%mx2+1)* (2*stars%mx3+1)-1
    !-odim
    IF (oneD%odd%d1) THEN
       oneD%odd%k3 = stars%mx3
       oneD%odd%nn2d = (2*(oneD%odd%k3) + 1)*(2*(oneD%odd%M) + 1)
    ELSE
       oneD%odd%k3 = 0 ; oneD%odd%M =0 ; oneD%odd%nn2d = 1
       oneD%odd%mb = 0
    ENDIF
    !-odim
    ALLOCATE ( atoms%nz(atoms%ntype),atoms%relax(3,atoms%ntype),atoms%nlhtyp(atoms%ntype))
    ALLOCATE ( sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd),stars%ustep(stars%ng3) )
    ALLOCATE ( stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3),stars%ig2(stars%ng3) )
    ALLOCATE ( atoms%jri(atoms%ntype),stars%kv2(2,stars%ng2),stars%kv3(3,stars%ng3),sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd) )
    ALLOCATE (sym%mrot(3,3,sym%nop),sym%tau(3,sym%nop))
    ALLOCATE ( atoms%lmax(atoms%ntype),sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd))!,sym%mrot(3,3,sym%nop) )
    ALLOCATE ( atoms%ncv(atoms%ntype),atoms%neq(atoms%ntype),sym%ngopr(atoms%nat) )
    ALLOCATE ( sphhar%nlh(sphhar%ntypsd),sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd) )
    ALLOCATE ( stars%nstr2(stars%ng2),sym%ntypsy(atoms%nat),stars%nstr(stars%ng3) )
    ALLOCATE ( stars%igfft(0:stars%kimax,2),stars%igfft2(0:stars%kimax2,2),atoms%nflip(atoms%ntype) )
    ALLOCATE ( atoms%econf(atoms%ntype) )
    ALLOCATE ( vacuum%izlay(vacuum%layerd,2) )
    ALLOCATE ( sym%invarop(atoms%nat,sym%nop),sym%invarind(atoms%nat) )
    ALLOCATE ( sym%multab(sym%nop,sym%nop),sym%invtab(sym%nop) )
    ALLOCATE ( sym%invsat(atoms%nat),sym%invsatnr(atoms%nat) )
    ALLOCATE ( atoms%lnonsph(atoms%ntype) )
    ALLOCATE ( atoms%dx(atoms%ntype),atoms%pos(3,atoms%nat))!,sym%tau(3,sym%nop) )
    ALLOCATE ( atoms%rmsh(atoms%jmtd,atoms%ntype),atoms%rmt(atoms%ntype),stars%sk2(stars%ng2),stars%sk3(stars%ng3) )
    ALLOCATE ( stars%phi2(stars%ng2) )
    ALLOCATE ( atoms%taual(3,atoms%nat),atoms%volmts(atoms%ntype),atoms%zatom(atoms%ntype) )
    ALLOCATE ( stars%rgphs(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,-stars%mx3:stars%mx3)  )
    ALLOCATE ( kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt) )
    ALLOCATE ( stars%pgfft(0:stars%kimax),stars%pgfft2(0:stars%kimax2) )
    ALLOCATE ( stars%ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1) )
    ALLOCATE ( atoms%bmu(atoms%ntype) )
    ALLOCATE ( atoms%l_geo(atoms%ntype) )
    ALLOCATE ( atoms%nlo(atoms%ntype),atoms%llo(atoms%nlod,atoms%ntype) )
    ALLOCATE ( atoms%lo1l(0:atoms%llod,atoms%ntype),atoms%nlol(0:atoms%llod,atoms%ntype),atoms%lapw_l(atoms%ntype) )
    ALLOCATE ( noco%alphInit(atoms%ntype),noco%alph(atoms%ntype),noco%beta(atoms%ntype),noco%l_relax(atoms%ntype) )
    ALLOCATE ( noco%b_con(2,atoms%ntype),atoms%lda_u(atoms%ntype),atoms%l_dulo(atoms%nlod,atoms%ntype) )
    ALLOCATE ( sym%d_wgn(-3:3,-3:3,3,sym%nop) )
    ALLOCATE ( atoms%ulo_der(atoms%nlod,atoms%ntype) )
    ALLOCATE ( kpts%ntetra(4,kpts%ntet), kpts%voltet(kpts%ntet))
    !+odim
    ALLOCATE ( oneD%ig1(-oneD%odd%k3:oneD%odd%k3,-oneD%odd%M:oneD%odd%M) )
    ALLOCATE ( oneD%kv1(2,oneD%odd%n2d),oneD%nstr1(oneD%odd%n2d) )
    ALLOCATE ( oneD%ngopr1(atoms%nat),oneD%mrot1(3,3,oneD%odd%nop),oneD%tau1(3,oneD%odd%nop) )
    ALLOCATE ( oneD%invtab1(oneD%odd%nop),oneD%multab1(oneD%odd%nop,oneD%odd%nop) )
    ALLOCATE ( oneD%igfft1(0:oneD%odd%nn2d-1,2),oneD%pgfft1(0:oneD%odd%nn2d-1) )
    stars%sk2(:) = 0.0 ; stars%phi2(:) = 0.0
    !-odim

    ! HF/hybrid functionals/EXX
    ALLOCATE ( hybrid%nindx(0:atoms%lmaxd,atoms%ntype) )

    input%l_coreSpec = .FALSE.





    CALL inped(atoms,vacuum,input,banddos,xcpot,sym,&
         cell,sliceplot,noco,&
         stars,oneD,hybrid,kpts,a1,a2,a3,namex,relcor,latnam,namgrp,grid)
    !
    IF (xcpot%needs_grad()) THEN
       ALLOCATE (stars%ft2_gfx(0:stars%kimax2),stars%ft2_gfy(0:stars%kimax2))
       ALLOCATE (oneD%pgft1x(0:oneD%odd%nn2d-1),oneD%pgft1xx(0:oneD%odd%nn2d-1),&
            oneD%pgft1xy(0:oneD%odd%nn2d-1),&
            oneD%pgft1y(0:oneD%odd%nn2d-1),oneD%pgft1yy(0:oneD%odd%nn2d-1))
    ELSE
       ALLOCATE (stars%ft2_gfx(0:1),stars%ft2_gfy(0:1))
       ALLOCATE (oneD%pgft1x(0:1),oneD%pgft1xx(0:1),oneD%pgft1xy(0:1),&
            oneD%pgft1y(0:1),oneD%pgft1yy(0:1))
    ENDIF
    oneD%odd%nq2 = stars%ng2!oneD%odd%n2d
    oneD%odi%nq2 = oneD%odd%nq2
    !-odim
    namex=xcpot%get_name()
    l_krla = xcpot%data%krla.EQ.1

    !Call xcpot%init(namex,l_krla,atoms%ntype)

    CALL setup(atoms,kpts,&
         sym,oneD,input,cell,&
         enpara,latnam,namgrp)

    banddos%l_orb = .FALSE.
    banddos%orbCompAtom = 0

    ALLOCATE(noco%socscale(atoms%ntype))
    xcpot%lda_atom(:) = .FALSE.
    noco%socscale(:) = 1.0

  END SUBROUTINE fleur_init_old
END MODULE m_fleur_init_old
