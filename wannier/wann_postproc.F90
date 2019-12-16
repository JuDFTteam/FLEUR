!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_wann_postproc
CONTAINS 
  SUBROUTINE wann_postproc(&
       stars,vacuum,atoms,sphhar,input,kpts,sym,mpi,&
       lapw,oneD,noco,cell,vTot,enpara,sliceplot,eig_id,l_real,&
       wann, fullnkpts, l_proj,ef,l_sgwf,fullnqpts)
    !     <          fermi_weights)

    !***********************************************
    !     Collection of those subroutines which may
    !     be called after the Wannier functions have
    !     been computed.
    !     Frank Freimuth
    !***********************************************
    USE m_types
    USE m_wann_dipole
    USE m_wann_dipole2
    USE m_wann_wannierize
    USE m_wann_hopping
    USE m_wann_plot_um_dat
    USE m_wannier_to_lapw
    USE m_wann_plot_from_lapw
    USE m_wann_nabla_rs
    USE m_wann_pauli_rs
    USE m_wann_nabla_pauli_rs
    USE m_wann_socmat_rs
    USE m_wann_perpmag_rs
    USE m_types
    USE m_wann_wigner_seitz
    USE m_wann_get_mp
    USE m_wann_get_kpts
    USE m_wann_fft4
    USE m_wann_fft5
    USE m_wann_rmat
    USE m_wann_nocoplot
#ifdef CPP_TOPO
    USE m_wann_torque_rs
    USE m_wann_offdiposop_rs
    USE m_wann_fft6
    USE m_wann_fft3
    USE m_wann_nedrho
#endif

#ifdef DCPP_WANN_EXT
    USE m_wannier_to_lapw_kpts
    USE m_wannier_lapw_gfleur
#endif
    IMPLICIT NONE

    
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_lapw),INTENT(IN)      :: lapw
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_potden),INTENT(IN)    :: vTot
    TYPE(t_enpara),INTENT(IN)    :: enpara
    type(t_sliceplot),INTENT(IN) :: sliceplot

    TYPE(t_wann), INTENT(in) :: wann
    LOGICAL,      INTENT(in) :: l_real
    INTEGER,      INTENT(in) :: eig_id
    INTEGER,      INTENT(in) :: fullnkpts,fullnqpts
    LOGICAL,      INTENT(in) :: l_proj
    REAL,    INTENT (in) :: ef
    LOGICAL, INTENT (in) :: l_sgwf
    !      real,intent(inout)   :: fermi_weights(:,:,:) !(input%neig,nkptd,jspd)

    CHARACTER(len=12) :: fending
    INTEGER :: i,nkpts,ikpt,nkqpts,iqpt
    REAL    :: delta3,time_lapw_expand,delta2,time_lapw_plot
    LOGICAL :: l_need_fft,l_file
    INTEGER :: hopmin_z,hopmax_z
    INTEGER :: hopmin_y,hopmax_y
    INTEGER :: hopmin_x,hopmax_x
    INTEGER             :: ii
    INTEGER             :: rvecnum,rvecind,num(3),int_dummy
    INTEGER,ALLOCATABLE :: rvec(:,:),ndegen(:)
    REAL                :: scale
    REAL,ALLOCATABLE    :: kpoints(:,:)!(3,fullnkpts)
    INTEGER             :: r1,r2,r3,spinspin,num_angl
    LOGICAL             :: l_nocosoc

    ALLOCATE(kpoints(3,fullnkpts))
    l_nocosoc=.FALSE.
    IF(noco%l_soc)l_nocosoc=.TRUE.
    IF(noco%l_noco)l_nocosoc=.TRUE.

    !***************************************************
    !     Read in the kpoints from w90kpts or kpts.
    !***************************************************      
    CALL wann_get_kpts(input,kpts,&
         wann%l_bzsym,input%film,oneD%odi%d1,.FALSE.,&
         nkpts,kpoints)
    IF (nkpts /= fullnkpts)&
         CALL juDFT_error ('mismatch in number of kpoints',&
         calledby="wann_postproc")

    CALL wann_get_kpts(input,kpts,&
         wann%l_bzsym,input%film,oneD%odi%d1,.TRUE.,&
         nkpts,kpoints)
    !*********************************************************
    !     Find out the structure of k-point set.
    !*********************************************************
    CALL wann_get_mp(&
         fullnkpts,kpoints,&
         num)

    IF(wann%l_wannierize.AND.mpi%irank==0)THEN
#ifndef CPP_WANN
       WRITE(*,*) 'At this point a wannierization has to be performed'
       WRITE(*,*) 'but the Wannier90 library is not linked!'
       CALL juDFT_error ('Wannierization without Wannier90 library',&
            calledby="wann_postproc")
#else
       CALL wann_wannierize(&
            input%film,wann%l_bzsym,input%jspins,&
            atoms%nat,atoms%pos,cell%amat,cell%bmat,atoms%ntype,atoms%neq,atoms%zatom)
#endif
    ENDIF

    IF(wann%l_dipole2.AND.mpi%irank==0)THEN
       CALL wann_dipole2(&
            input%jspins,atoms%pos,cell%omtil,atoms%nat,&
            (noco%l_soc.OR.noco%l_noco))
    ENDIF

    IF(wann%l_dipole.AND.mpi%irank==0)THEN
       CALL wann_dipole(&
            input%jspins,cell%omtil,atoms%nat,atoms%pos,cell%amat,&
            atoms%ntype,atoms%neq,atoms%zatom)
    ENDIF

#ifdef CPP_TOPO
    IF(wann%l_nedrho.AND.mpi%irank==0)THEN
       CALL wann_nedrho(&
            fullnkpts,cell%area,ef,&
            fermi_weights)
    ENDIF
#endif

    IF(wann%l_byenergy) i=3
    IF(wann%l_bynumber) i=2
    IF(wann%l_byindex)  i=1

    IF(wann%l_plotw90)THEN
       CALL juDFT_error("not in this version",&
            calledby ="wann_postproc")

       !         call wann_plotw90(i,wann%band_min,wann%band_max,numbands,nwfs,
       !     >   atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,atoms%ntype,
       !     >      input%neig,atoms%nat,sym%nop,lapw%dim_nvd(),jspd,lapw%dim_nbasfcn(),atoms%llod,atoms%nlod,atoms%ntype,
       !     >      nwdd,cell%omtil,atoms%nlo,atoms%llo,atoms%lapw_l,sym%invtab,sym%mrot,sym%ngopr,atoms%neq,atoms%lmax,
       !     >      sym%invsat,sym%invsatnr,kpts%nkpt,atoms%taual,atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,
       !     >      noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,mpi%irank,mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,
       !     >      size(atoms%rmsh,1),sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,
       !     >      sym%ntypsy,input%jspins,kpts%nkpt,atoms%dx,stars%ng2,atoms%rmsh,sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,
       !     >      stars%ustep,stars%ig,stars%mx1,stars%mx2,stars%mx3,stars%rgphs,sliceplot%slice,sliceplot%kk,sliceplot%nnne,
       !     >      cell%z1,lapw%dim_nv2d(),vacuum%nmzxy,vacuum%nmz,vacuum%delz,stars%ig2,cell%area,sym%tau,atoms%zatom,stars%ng2,sym%nop2,
       !     >      cell%volint,sym%symor,atoms%pos,ef,wann%l_bzsym,irecl)

    ENDIF

    l_need_fft=.FALSE.
    IF(wann%l_hopping.OR.wann%l_nablars.OR.&
         wann%l_nablapaulirs.OR.wann%l_pauli.OR.&
         wann%l_perpmagrs.OR.wann%l_socmatrs.OR.&
         wann%l_torquers.OR.wann%l_offdiposoprs.OR.&
         wann%l_socspicomrs.OR.wann%l_spindisprs.OR.&
         wann%l_anglmomrs .OR.wann%l_perturbrs.OR.&
         wann%l_orbcomprs.OR.wann%l_rmat)l_need_fft=.TRUE.

    IF(l_need_fft.AND.mpi%irank==0)THEN

       IF(.FALSE.)THEN !specify r-mesh by its boundaries
          hopmin_z=-5;hopmax_z=5
          hopmin_x=0;hopmax_x=0
          hopmin_y=0;hopmax_y=0
          rvecnum=(hopmax_z-hopmin_z+1)
          IF(.NOT.oneD%odi%d1.AND.input%film)THEN
             hopmin_x=-5;hopmax_x=5
             hopmin_y=-5;hopmax_y=5
             hopmin_z=0;     hopmax_z=0
          ELSE
             hopmin_x=-5;hopmax_x=5
             hopmin_y=-5;hopmax_y=5
          ENDIF
          rvecnum=        (hopmax_z-hopmin_z+1)
          rvecnum=rvecnum*(hopmax_y-hopmin_y+1)
          rvecnum=rvecnum*(hopmax_x-hopmin_x+1)

          ALLOCATE(rvec(3,rvecnum))
          rvecind=0
          DO r3=hopmin_z,hopmax_z
             DO r2=hopmin_y,hopmax_y
                DO r1=hopmin_x,hopmax_x
                   rvecind=rvecind+1
                   IF (rvecind > rvecnum)&
                        CALL juDFT_error ('mismatch in number of kpoints',&
                        calledby="wann_postproc")

                   rvec(1,rvecind)=r1
                   rvec(2,rvecind)=r2
                   rvec(3,rvecind)=r3
                ENDDO !r1
             ENDDO !r2
          ENDDO !r3
       ELSE !determine optimal r-mesh
          CALL wann_wigner_seitz(&
               .TRUE.,num,cell%amat,0,&
               rvecnum,rvec,ndegen)
          ALLOCATE(rvec(3,rvecnum))
          ALLOCATE(ndegen(rvecnum))
          CALL wann_wigner_seitz(&
               .FALSE.,num,cell%amat,rvecnum,&
               int_dummy,rvec,ndegen)

          OPEN(333,file='wig_vectors',recl=1000)
          DO ii=1,rvecnum
             WRITE(333,*)ii,rvec(1,ii),rvec(2,ii),rvec(3,ii),&
                  ndegen(ii)
          ENDDO
          CLOSE(333)

       ENDIF

    ENDIF

    IF(wann%l_hopping.AND.mpi%irank==0) THEN
       CALL wann_hopping(&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            l_nocosoc,wann%band_min,wann%band_max,&
            input%neig,wann%l_socmmn0,wann%l_ndegen,ndegen,&
            wann%wan90version,wann%l_unformatted)
    ENDIF

    IF(wann%l_nablars.AND.mpi%irank==0) THEN
       CALL wann_nabla_rs(&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,&
            input%neig,.FALSE.,&
            wann%wan90version)
    ENDIF

    IF(wann%l_orbcomprs.AND.mpi%irank==0)THEN
       num_angl=9
       IF(wann%l_oc_f)num_angl=16
       CALL wann_fft5(&
            .FALSE.,&
            num_angl,wann%oc_num_orbs,&
            'WF1.orbcomp','orbcomp.1',&
            rvecnum,rvec,kpoints,&
            1,fullnkpts,input%film,&
            noco%l_soc,wann%band_min,wann%band_max,input%neig,&
            wann%wan90version)
       IF( input%jspins.EQ.2 )THEN
          spinspin=2
          IF(noco%l_soc.OR.noco%l_noco)spinspin=1
          CALL wann_fft5(&
               .FALSE.,&
               num_angl,wann%oc_num_orbs,&
               'WF2.orbcomp','orbcomp.2',&
               rvecnum,rvec,kpoints,&
               spinspin,fullnkpts,input%film,&
               noco%l_soc,wann%band_min,wann%band_max,input%neig,&
               wann%wan90version)
       ENDIF
    ENDIF

    IF(wann%l_rmat)THEN
       CALL wann_rmat(&
            cell%bmat,cell%amat,&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            l_nocosoc,wann%band_min,wann%band_max,&
            input%neig,.FALSE.,wann%wan90version)
    ENDIF


    IF(wann%l_anglmomrs.AND.mpi%irank==0)THEN
       CALL wann_fft4(&
            'WF1.anglmom',&
            'anglmomrs.1',.FALSE.,&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,input%neig,&
            .FALSE.,wann%wan90version)
    ENDIF

#ifdef CPP_TOPO
    IF(wann%l_offdiposoprs.AND.mpi%irank==0)THEN
       CALL wann_offdiposop_rs(&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,input%neig,&
            .FALSE.)
    ENDIF

    IF(wann%l_spindisprs.AND.mpi%irank==0)THEN
       CALL wann_fft5(&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,input%neig,&
            .FALSE.)
    ENDIF

    IF(wann%l_perturbrs.AND.mpi%irank==0)THEN
       CALL wann_fft3(&
            'WF1.perturb'  ,'perturbrs.1'  ,.FALSE.,&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,input%neig,&
            .FALSE.)
    ENDIF

    IF(wann%l_socspicomrs.AND.mpi%irank==0)THEN
       IF(noco%l_soc)THEN
          CALL wann_fft4(&
               'WF1.socspicom','socspicomrs.1',.TRUE.,&
               rvecnum,rvec,kpoints,&
               input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
               noco%l_soc,wann%band_min,wann%band_max,input%neig,&
               .FALSE.)
       ELSE
          CALL wann_fft6(&
               rvecnum,rvec,kpoints,&
               input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
               noco%l_soc,wann%band_min,wann%band_max,input%neig,&
               .FALSE.)
       ENDIF
    ENDIF



    IF(wann%l_torquers.AND.mpi%irank==0) THEN
       CALL wann_torque_rs(&
            atoms%ntype,atoms%neq,rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,&
            input%neig,.FALSE.)
    ENDIF
#endif
    IF (wann%l_nablapaulirs.AND.mpi%irank==0)THEN
       CALL wann_nabla_pauli_rs(&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,&
            input%neig,.FALSE.,wann%wan90version)
    ENDIF

    IF (wann%l_pauli.AND.mpi%irank==0)THEN
       CALL wann_pauli_rs(&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            l_nocosoc,wann%band_min,wann%band_max,&
            input%neig,.FALSE.,wann%l_ndegen,ndegen,&
            wann%wan90version,wann%l_unformatted)
    ENDIF

    IF (wann%l_perpmagrs.AND.mpi%irank==0)THEN
       CALL wann_perpmag_rs(&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,&
            input%neig,.FALSE.,wann%l_ndegen,ndegen,wann%wan90version,&
            wann%l_unformatted)
    ENDIF

    IF (wann%l_socmatrs.AND.mpi%irank==0)THEN
       CALL wann_socmat_rs(&
            rvecnum,rvec,kpoints,&
            input%jspins,fullnkpts,wann%l_bzsym,input%film,oneD%odi%d1,&
            noco%l_soc,wann%band_min,wann%band_max,&
            input%neig,.FALSE.,wann%wan90version)
    ENDIF

    IF(wann%l_plot_umdat)THEN
       CALL wann_plot_um_dat(&
            stars,vacuum,atoms,sphhar,input,sym,mpi,&
            lapw,oneD,noco,cell,vTot,enpara,eig_id,l_real,&
            mpi%mpi_comm,i,wann%band_min,wann%band_max,noco%l_soc,&
            atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,atoms%ntype,&
            input%neig,atoms%nat,sym%nop,lapw%dim_nvd(),input%jspins,lapw%dim_nbasfcn(),atoms%llod,atoms%nlod,atoms%ntype,&
            cell%omtil,atoms%nlo,atoms%llo,atoms%lapw_l,sym%invtab,sym%mrot,sym%ngopr,atoms%neq,atoms%lmax,&
            sym%invsat,sym%invsatnr,kpts%nkpt,atoms%taual,atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,&
            noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,mpi%irank,mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,&
            SIZE(atoms%rmsh,1),sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,&
            sym%ntypsy,input%jspins,kpts%nkpt,atoms%dx,stars%ng2,atoms%rmsh,sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,&
            stars%ustep,stars%ig,stars%mx1,stars%mx2,stars%mx3,stars%rgphs,sliceplot%slice,sliceplot%kk,sliceplot%nnne,&
            cell%z1,lapw%dim_nv2d(),vacuum%nmzxy,vacuum%nmz,vacuum%delz,stars%ig2,cell%area,sym%tau,atoms%zatom,stars%ng2,sym%nop2,&
            cell%volint,sym%symor,atoms%pos,ef,wann%l_bzsym,wann%l_proj_plot,&
            wann%wan90version)
    ENDIF


    CALL CPU_TIME(delta2)
    IF(wann%l_lapw.AND.mpi%irank==0)THEN
       CALL wannier_to_lapw(&
            mpi%mpi_comm,eig_id,l_real,&
            input,lapw,oneD,noco,sym,cell,atoms,stars,vacuum,sphhar,&
            vTot,&
            noco%l_soc,wann%unigrid,i,wann%band_min,wann%band_max,&
            atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,atoms%ntype,&
            input%neig,atoms%nat,sym%nop,lapw%dim_nvd(),input%jspins,lapw%dim_nbasfcn(),atoms%llod,atoms%nlod,atoms%ntype,&
            cell%omtil,atoms%nlo,atoms%llo,atoms%lapw_l,sym%invtab,sym%mrot,sym%ngopr,atoms%neq,atoms%lmax,&
            sym%invsat,sym%invsatnr,kpts%nkpt,atoms%taual,atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,&
            noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,mpi%irank,mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,&
            SIZE(atoms%rmsh,1),sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,&
            sym%ntypsy,input%jspins,kpts%nkpt,atoms%dx,stars%ng2,atoms%rmsh,sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,&
            stars%ustep,stars%ig,stars%mx1,stars%mx2,stars%mx3,stars%rgphs,sliceplot%slice,sliceplot%kk,sliceplot%nnne,&
            cell%z1,lapw%dim_nv2d(),vacuum%nmzxy,vacuum%nmz,vacuum%delz,stars%ig2,cell%area,sym%tau,atoms%zatom,stars%ng2,sym%nop2,&
            cell%volint,sym%symor,atoms%pos,ef,wann%l_bzsym,wann%l_proj_plot,&
            wann%wan90version)
    ENDIF

    IF(wann%l_lapw_kpts)THEN
#ifdef DCPP_WANN_EXT
       CALL wannier_to_lapw_kpts(&
            unigrid,i,wann%band_min,wann%band_max,&
            atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,atoms%ntype,&
            input%neig,atoms%nat,sym%nop,lapw%dim_nvd(),input%jspins,lapw%dim_nbasfcn(),atoms%llod,atoms%nlod,atoms%ntype,&
            nwdd,cell%omtil,atoms%nlo,atoms%llo,atoms%lapw_l,sym%invtab,sym%mrot,sym%ngopr,atoms%neq,atoms%lmax,&
            sym%invsat,sym%invsatnr,kpts%nkpt,atoms%taual,atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,&
            noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,mpi%irank,mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,&
            SIZE(atoms%rmsh,1),sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,&
            sym%ntypsy,input%jspins,kpts%nkpt,atoms%dx,stars%ng2,atoms%rmsh,sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,&
            stars%ustep,stars%ig,stars%mx1,stars%mx2,stars%mx3,stars%rgphs,sliceplot%slice,sliceplot%kk,sliceplot%nnne,&
            cell%z1,lapw%dim_nv2d(),vacuum%nmzxy,vacuum%nmz,vacuum%delz,stars%ig2,cell%area,sym%tau,atoms%zatom,stars%ng2,sym%nop2,&
            cell%volint,sym%symor,atoms%pos,ef,wann%l_bzsym,wann%l_proj_plot,irecl)
#else
       CALL juDFT_error("not yet tested in this release",calledby&
            ="wann_postproc")
#endif
    ENDIF
    IF(wann%l_lapw_gfleur)THEN
#ifdef DCPP_WANN_EXT
       CALL wannier_lapw_gfleur(&
            gfthick,gfcut,i,wann%band_min,wann%band_max,&
            atoms%l_dulo,noco%l_noco,noco%l_ss,atoms%lmaxd,atoms%ntype,&
            input%neig,atoms%nat,sym%nop,lapw%dim_nvd(),input%jspins,lapw%dim_nbasfcn(),atoms%llod,atoms%nlod,atoms%ntype,&
            nwdd,cell%omtil,atoms%nlo,atoms%llo,atoms%lapw_l,sym%invtab,sym%mrot,sym%ngopr,atoms%neq,atoms%lmax,&
            sym%invsat,sym%invsatnr,kpts%nkpt,atoms%taual,atoms%rmt,cell%amat,cell%bmat,cell%bbmat,noco%alph,&
            noco%beta,noco%qss,stars%sk2,stars%phi2,oneD%odi,oneD%ods,mpi%irank,mpi%isize,stars%ng3,vacuum%nmzxyd,vacuum%nmzd,&
            SIZE(atoms%rmsh,1),sphhar%nlhd,stars%ng3,vacuum%nvac,sym%invs,sym%invs2,input%film,sphhar%nlh,atoms%jri,sphhar%ntypsd,&
            sym%ntypsy,input%jspins,kpts%nkpt,atoms%dx,stars%ng2,atoms%rmsh,sliceplot%e1s,sliceplot%e2s,atoms%ulo_der,&
            stars%ustep,stars%ig,stars%mx1,stars%mx2,stars%mx3,stars%rgphs,sliceplot%slice,sliceplot%kk,sliceplot%nnne,&
            cell%z1,lapw%dim_nv2d(),vacuum%nmzxy,vacuum%nmz,vacuum%delz,stars%ig2,cell%area,sym%tau,atoms%zatom,stars%ng2,sym%nop2,&
            cell%volint,sym%symor,atoms%pos,ef,wann%l_bzsym,wann%l_proj_plot,irecl)
#else
       CALL juDFT_error("not yet tested in this release",calledby&
            ="wann_postproc")
#endif
    ENDIF

    CALL CPU_TIME(delta3)
    time_lapw_expand=delta3-delta2

    CALL CPU_TIME(delta2)
    IF(wann%l_plot_lapw.AND.mpi%irank==0)THEN
       CALL juDFT_error("not yet tested in this release",calledby&
            ="wann_postproc")
       !         call wann_plot_from_lapw(
       !     >     lapw%dim_nv2d(),input%jspins,oneD%odi,oneD%ods,stars%ng3,vacuum%nmzxyd,stars%ng2,
       !     >     sphhar%ntypsd,
       !     >     atoms%ntype,atoms%lmaxd,size(atoms%rmsh,1),atoms%nat,vacuum%nmzd,atoms%neq,stars%ng3,vacuum%nvac,
       !     >     vacuum%nmz,vacuum%nmzxy,stars%ng2,sym%nop,sym%nop2,cell%volint,input%film,sliceplot%slice,sym%symor,
       !     >     sym%invs,sym%invs2,cell%z1,vacuum%delz,sym%ngopr,sym%ntypsy,atoms%jri,atoms%pos,atoms%zatom,
       !     >     atoms%lmax,sym%mrot,sym%tau,atoms%rmsh,sym%invtab,cell%amat,cell%bmat,cell%bbmat,sliceplot%nnne,sliceplot%kk,
       !     >     atoms%nlod,atoms%llod,lmd,cell%omtil,atoms%nlo,atoms%llo)
    ENDIF
    CALL CPU_TIME(delta3)
    time_lapw_plot=delta3-delta2

    WRITE(6,*)"time_lapw_expand=",time_lapw_expand
    WRITE(6,*)"time_lapw_plot=",time_lapw_plot

    IF(wann%l_finishnocoplot.AND.mpi%irank==0) THEN
       !         write(*,*)'doing the UNK mixing'
       !         write(*,*)'noco%alph',noco%alph
       !         write(*,*)'noco%beta',noco%beta
       !         nkqpts=fullnkpts
       !         if(l_sgwf) nkqpts=fullnkpts*fullnqpts 
       CALL wann_nocoplot(atoms,sliceplot%slice,sliceplot%nnne,&!nslibd&
            cell%amat,cell%bmat,fullnkpts,oneD%odi,input%film,&
            atoms%nat,atoms%ntype,SIZE(atoms%rmsh,1),atoms%ntype,atoms%neq,atoms%pos,&
            atoms%jri,atoms%rmsh,noco%alph,noco%beta,fullnqpts,noco%qss,&
            cell%z1,atoms%zatom)
    ENDIF

  END SUBROUTINE wann_postproc
END MODULE m_wann_postproc
