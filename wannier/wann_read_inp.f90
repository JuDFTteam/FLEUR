!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_wann_read_inp

contains
subroutine wann_read_inp(input,l_p0,wann)
!********************************************
!     Read the Wannier input file 'wann_inp'.
!     Frank Freimuth
!********************************************
   use m_judft
   use m_types_setup

   implicit none

   TYPE(t_input), intent(in)   :: input
   TYPE(t_wann), intent(inout) :: wann
   logical,intent(in)          :: l_p0

   logical           :: l_file,l_orbcompinp
   integer           :: i,ios,n
   character(len=30) :: task
   real              :: version_real

!-----some defaults
   wann%l_atomlist=.false.
   wann%l_ndegen=.false.
   wann%l_orbitalmom=.false.
   wann%l_orbcomp=.false.
   wann%l_orbcomprs=.false.
   wann%l_perturbrs=.false.
   wann%l_denmat=.false.                                             
   wann%l_perturb=.false.
   wann%l_nedrho=.false.
   wann%l_anglmomrs=.false.
   wann%l_anglmom=.false.
   wann%l_socspicom=.false.
   wann%l_socspicomrs=.false.
   wann%l_offdiposoprs=.false.
   wann%l_offdiposop=.false.
   wann%l_spindisp=.false.
   wann%l_spindisprs=.false.
   wann%l_torque=.false.
   wann%l_torquers=.false.
   wann%l_berry=.false.
   wann%l_perpmag=.false.
   wann%l_perpmagrs=.false.
   wann%l_perpmagat=.false.
   wann%l_perpmagatrs=.false.
   wann%l_socmatrs=.false.
   wann%l_socmat=.false.
   wann%l_socmatvec=.false.
   wann%l_soctomom=.false.
   wann%l_kptsreduc2=.false.
   wann%l_nablars=.false.
   wann%l_nablapaulirs=.false.
   wann%l_surfcurr=.false.
   wann%l_updown=.false.
   wann%l_ahe=.false.
   wann%l_she=.false.
   wann%l_rmat=.false.
   wann%l_nabla=.false.
   wann%l_socodi=.false.
   wann%l_pauli=.false.
   wann%l_pauliat=.false.
   wann%l_potmat=.false.
   wann%l_projgen=.false.
   wann%l_stopopt=.false.
   wann%l_kpointgen=.false.
   wann%l_w90kpointgen=.false.
   wann%l_plot_symm=.false.
   wann%l_socmmn0=.false.
   wann%l_bzsym=.false.
   wann%l_hopping=.false.
   wann%l_kptsreduc=.false.
   wann%l_prepwan90=.false.
   wann%l_plot_umdat=.false.
   wann%l_wann_plot=.false.
   wann%l_bynumber=.false.
   wann%l_matrixmmn=.false.
   wann%l_matrixamn=.false.
   wann%l_projmethod=.false.
   wann%l_wannierize=.false.
   wann%l_plotw90=.false.
   wann%l_byenergy=.false.
   wann%l_proj_plot=.false.
   wann%l_bestproj=.false.
   wann%l_ikptstart=.false.
   wann%l_lapw=.false.
   wann%l_plot_lapw=.false.
   wann%l_fermi=.false.
   wann%l_dipole=.false.
   wann%l_dipole2=.false.
   wann%l_dipole3=.false.
   wann%l_mmn0=.false.
   wann%l_mmn0at=.false.
   wann%l_manyfiles=.false.
   wann%l_collectmanyfiles=.false.
   wann%l_ldauwan=.false.
   wann%l_lapw_kpts=.false.
   wann%l_lapw_gfleur=.false.
   wann%l_unformatted=.false.
   wann%l_finishnocoplot=.false.
   wann%l_finishgwf=.false.
   wann%l_skipkov=.false.
   wann%l_matrixuHu=.false.
   wann%l_matrixuHu_dmi=.false.
   wann%ikptstart=1
   wann%wan90version=2 ! Set the standard to Wannier90-1.2
   wann%l_mmn0_unf_to_spn_unf=.false.
   wann%l_mmn0_to_spn_unf=.false.
   wann%l_mmn0_to_spn=.false.
   wann%l_mmn0_to_spn2=.false.
   wann%l_mmn0_unf_to_spn=.false.
   wann%l_perpmag_unf_to_tor_unf=.false.
   wann%l_perpmag_to_tor_unf=.false.
   wann%l_perpmag_to_tor=.false.
   wann%l_perpmag_unf_to_tor=.false.
   wann%l_hsomtxvec_unf_to_lmpzsoc_unf=.false.
   wann%l_hsomtxvec_to_lmpzsoc_unf=.false.
   wann%l_hsomtxvec_to_lmpzsoc=.false.
   wann%l_hsomtxvec_unf_to_lmpzsoc=.false.
   wann%l_hsomtx_unf_to_hsoc_unf=.false.
   wann%l_hsomtx_to_hsoc_unf=.false.
   wann%l_hsomtx_to_hsoc=.false.
   wann%l_hsomtx_unf_to_hsoc=.false.

!-----read the input file 'wann_inp'
   l_file=.false.
   inquire(file='wann_inp',exist=l_file)
   IF ((.not.l_file).AND.(.NOT.input%l_inpXML)) THEN
      CALL juDFT_error ("wann_inp not found", calledby="wann_read_inp")
   END IF
   IF (l_file) THEN
      wann%band_min(1:2)=-1
      wann%band_max(1:2)=-1
      wann%l_byindex=.false.
      open(916,file='wann_inp',form='formatted')
      i=0
      do 
         i=i+1
         read(916,'(a)',iostat=ios)task
         if(ios.ne.0)exit
         if(l_p0) write(6,*)"line ",i,":",task
         if(task(1:1).eq.'!')cycle
         if(index(task,'=F').ne.0.or.index(task,'=f').ne.0)cycle
         read(task,*,iostat=ios)task
         if(index(task,'=').ne.0)then
            task=task(1:(index(task,'=')-1))
         endif
         if(l_p0) write(6,*)"==>key: ",task
         if(trim(task).eq.'endjobs')then
            exit
         elseif(trim(task).eq.'nabla')then
            wann%l_nabla=.true.
         elseif(trim(task).eq.'kptsreduc2')then
            wann%l_kptsreduc2=.true.
            backspace(916)
            read(916,*,iostat=ios)task,wann%mhp(1),wann%mhp(2),wann%mhp(3)
            if (ios /= 0) CALL juDFT_error ("error reading mhp", calledby="wann_read_inp")
            if(l_p0)write(6,*)"mhp=",wann%mhp(1),wann%mhp(2),wann%mhp(3)
         elseif(trim(task).eq.'nablars')then
            wann%l_nablars=.true.
         elseif(trim(task).eq.'nablapaulirs')then
            wann%l_nablapaulirs=.true.
         elseif(trim(task).eq.'socspicom')then
            wann%l_socspicom=.true.
         elseif(trim(task).eq.'orbitalmom')then
            wann%l_orbitalmom=.true.
         elseif(trim(task).eq.'anglmom')then
            wann%l_anglmom=.true.
         elseif(trim(task).eq.'ndegen')then
            wann%l_ndegen=.true.
         elseif(trim(task).eq.'unformatted')then
            wann%l_unformatted=.true.
         elseif(trim(task).eq.'orbcomp')then
            wann%l_orbcomp=.true.
         elseif(trim(task).eq.'orbcomprs')then
            wann%l_orbcomprs=.true.
         elseif(trim(task).eq.'denmat')then
            wann%l_denmat=.true.
         elseif(trim(task).eq.'nedrho')then
            wann%l_nedrho=.true.
         elseif(trim(task).eq.'perturb')then
            wann%l_perturb=.true.
         elseif(trim(task).eq.'perturbrs')then
            wann%l_perturbrs=.true.
         elseif(trim(task).eq.'rmat')then
            wann%l_rmat=.true.
         elseif(trim(task).eq.'anglmomrs')then
            wann%l_anglmomrs=.true.
         elseif(trim(task).eq.'socspicomrs')then
            wann%l_socspicomrs=.true.
         elseif(trim(task).eq.'she')then
            wann%l_she=.true.
         elseif(trim(task).eq.'berry')then
            wann%l_berry=.true.
         elseif(trim(task).eq.'offdiposop')then
            wann%l_offdiposop=.true.
         elseif(trim(task).eq.'offdiposoprs')then
            wann%l_offdiposoprs=.true.
         elseif(trim(task).eq.'spindisp')then
            wann%l_spindisp=.true.
         elseif(trim(task).eq.'spindisprs')then
            wann%l_spindisprs=.true.
         elseif(trim(task).eq.'torque')then
            wann%l_torque=.true.
         elseif(trim(task).eq.'torquers')then
            wann%l_torquers=.true.
         elseif(trim(task).eq.'perpmag')then
            wann%l_perpmag=.true.
         elseif(trim(task).eq.'perpmagrs')then
            wann%l_perpmagrs=.true.
         elseif(trim(task).eq.'perpmagat')then
            wann%l_perpmagat=.true.
         elseif(trim(task).eq.'perpmagatrs')then
            wann%l_perpmagatrs=.true.
         elseif(trim(task).eq.'mmn0_unf_to_spn_unf')then
            wann%l_mmn0_unf_to_spn_unf=.true.
         elseif(trim(task).eq.'mmn0_to_spn_unf')then
            wann%l_mmn0_to_spn_unf=.true.
         elseif(trim(task).eq.'mmn0_to_spn')then
            wann%l_mmn0_to_spn=.true.
         elseif(trim(task).eq.'mmn0_to_spn2')then
            wann%l_mmn0_to_spn2=.true.
         elseif(trim(task).eq.'mmn0_unf_to_spn')then
            wann%l_mmn0_unf_to_spn=.true.
         elseif(trim(task).eq.'perpmag_unf_to_tor_unf')then
            wann%l_perpmag_unf_to_tor_unf=.true.
         elseif(trim(task).eq.'perpmag_to_tor_unf')then
            wann%l_perpmag_to_tor_unf=.true.
         elseif(trim(task).eq.'perpmag_to_tor')then
            wann%l_perpmag_to_tor=.true.
         elseif(trim(task).eq.'perpmag_unf_to_tor')then
            wann%l_perpmag_unf_to_tor=.true.
         elseif(trim(task).eq.'hsomtxvec_unf_to_lmpzsoc_unf')then
            wann%l_hsomtxvec_unf_to_lmpzsoc_unf=.true.
         elseif(trim(task).eq.'hsomtxvec_to_lmpzsoc_unf')then
            wann%l_hsomtxvec_to_lmpzsoc_unf=.true.
         elseif(trim(task).eq.'hsomtxvec_to_lmpzsoc')then
            wann%l_hsomtxvec_to_lmpzsoc=.true.
         elseif(trim(task).eq.'hsomtxvec_unf_to_lmpzsoc')then
            wann%l_hsomtxvec_unf_to_lmpzsoc=.true.  
         elseif(trim(task).eq.'hsomtx_unf_to_hsoc_unf')then
            wann%l_hsomtx_unf_to_hsoc_unf=.true.
         elseif(trim(task).eq.'hsomtx_to_hsoc_unf')then
            wann%l_hsomtx_to_hsoc_unf=.true.
         elseif(trim(task).eq.'hsomtx_to_hsoc')then
            wann%l_hsomtx_to_hsoc=.true.
         elseif(trim(task).eq.'hsomtx_unf_to_hsoc')then
            wann%l_hsomtx_unf_to_hsoc=.true.
         elseif(trim(task).eq.'socmat')then
            wann%l_socmat=.true.
         elseif(trim(task).eq.'socmatvec')then
            wann%l_socmatvec=.true.  
         elseif(trim(task).eq.'socmatrs')then
            wann%l_socmatrs=.true.
         elseif(trim(task).eq.'soctomom')then
            wann%l_soctomom=.true.
         elseif(trim(task).eq.'surfcurr')then
            wann%l_surfcurr=.true.
         elseif(trim(task).eq.'lapw_kpts')then
            wann%l_lapw_kpts=.true.
         elseif(trim(task).eq.'updown')then
            wann%l_updown=.true.
         elseif(trim(task).eq.'stopopt')then
            wann%l_stopopt=.true.
         elseif(trim(task).eq.'projgen')then
            wann%l_projgen=.true.
         elseif(trim(task).eq.'kpointgen')then
            wann%l_kpointgen=.true.
         elseif(trim(task).eq.'potmat')then
            wann%l_potmat=.true.
         elseif(trim(task).eq.'w90kpointgen')then
            wann%l_w90kpointgen=.true.
         elseif(trim(task).eq.'lapw_gfleur')then
            wann%l_lapw_gfleur=.true.
            backspace(916)
            read(916,*,iostat=ios)task,wann%gfthick,wann%gfcut
            if (ios /= 0) CALL juDFT_error ("error reading gfcut", calledby="wann_read_inp")
            if(l_p0)write(6,*)"gfcut=",wann%gfthick,wann%gfcut
         elseif(trim(task).eq.'lapw')then
            wann%l_lapw=.true.
            backspace(916)
            read(916,*,iostat=ios)task,wann%unigrid(:)
            if (ios /= 0) CALL juDFT_error ("error reading unigrid", calledby="wann_read_inp")
            if(l_p0)write(6,*)"unigrid=",wann%unigrid(:)
         elseif(trim(task).eq.'plot_lapw')then
            wann%l_plot_lapw=.true.
         elseif(trim(task).eq.'bzsym')then
            wann%l_bzsym=.true.
         elseif(trim(task).eq.'mmn0')then
            wann%l_mmn0=.true.
         elseif(trim(task).eq.'mmn0at')then
            wann%l_mmn0at=.true.
         elseif(trim(task).eq.'manyfiles')then
            wann%l_manyfiles=.true.
         elseif(trim(task).eq.'collectmanyfiles')then
            wann%l_collectmanyfiles=.true.
         elseif(trim(task).eq.'bestproj')then
            wann%l_bestproj=.true.
         elseif(trim(task).eq.'pauli')then
            wann%l_pauli=.true.
         elseif(trim(task).eq.'pauliat')then
            wann%l_pauliat=.true.
         elseif(trim(task).eq.'proj_plot')then
            wann%l_proj_plot=.true.
         elseif(trim(task).eq.'hopping')then
            wann%l_hopping=.true.
         elseif(trim(task).eq.'plot_symm')then
            wann%l_plot_symm=.true.
         elseif(trim(task).eq.'kptsreduc')then
            wann%l_kptsreduc=.true.
         elseif(trim(task).eq.'fermi')then
            wann%l_fermi=.true.
         elseif(trim(task).eq.'prepwan90')then
            wann%l_prepwan90=.true.
         elseif(trim(task).eq.'plot_umdat')then
            wann%l_plot_umdat=.true.
         elseif(trim(task).eq.'wann_plot')then
            wann%l_wann_plot=.true.
         elseif(trim(task).eq.'bynumber')then
            wann%l_bynumber=.true.
         elseif(trim(task).eq.'matrixmmn')then
            wann%l_matrixmmn=.true.
         elseif(trim(task).eq.'projmethod')then
            wann%l_projmethod=.true.
         elseif(trim(task).eq.'matrixamn')then
            wann%l_matrixamn=.true.
         elseif(trim(task).eq.'wannierize')then
            wann%l_wannierize=.true.
         elseif(trim(task).eq.'plotw90')then
            wann%l_plotw90=.true.
         elseif(trim(task).eq.'dipole')then
            wann%l_dipole=.true.
         elseif(trim(task).eq.'dipole3')then
            wann%l_dipole3=.true.
         elseif(trim(task).eq.'ldauwan')then
            wann%l_ldauwan=.true.
         elseif(trim(task).eq.'byenergy')then
            wann%l_byenergy=.true.
         elseif(trim(task).eq.'finishnocoplot') then
            wann%l_finishnocoplot=.true.
         elseif(trim(task).eq.'finishgwf') then
            wann%l_finishgwf=.true.
         elseif(trim(task).eq.'skipkov') then
            wann%l_skipkov=.true.
         elseif(trim(task).eq.'matrixuhu') then
            wann%l_matrixuHu=.true.
         elseif(trim(task).eq.'matrixuhu-dmi') then
            wann%l_matrixuHu_dmi=.true.
         elseif(trim(task).eq.'wan90version')then
            backspace(916)
            read(916,*,iostat=ios)task,version_real
            if (ios /= 0) CALL judft_error("error reading wan90version", calledby="wann_read_inp")
            if(abs(version_real-1.1).lt.1.e-9)then
               wann%wan90version=1
            elseif(abs(version_real-1.2).lt.1.e-9)then
               wann%wan90version=2
            elseif(abs(version_real-2.0).lt.1.e-9)then
               wann%wan90version=3
            else
              CALL judft_error ("chosen w90 version unknown", calledby="wann_read_inp")
            endif
         elseif(trim(task).eq.'atomlist')then
            wann%l_atomlist=.true.
            backspace(916)
            read(916,*,iostat=ios)task,wann%atomlist_num
            if (ios /= 0) CALL judft_error ("error reading atomlist_num", calledby="wann_read_inp")
            if(allocated(wann%atomlist))deallocate(wann%atomlist)
            allocate(wann%atomlist(wann%atomlist_num))
            backspace(916)
            read(916,*,iostat=ios)task,wann%atomlist_num,wann%atomlist
            if (ios /= 0) CALL judft_error ("error reading atomlist", calledby="wann_read_inp")
            if(l_p0)write(6,*)"atomlist_num=",wann%atomlist_num
            if(l_p0)write(6,*)"atomlist=",wann%atomlist
         elseif(trim(task).eq.'byindex')then
            wann%l_byindex=.true.
            backspace(916)
            read(916,*,iostat=ios)task,wann%band_min(1),wann%band_max(1)
            if (ios /= 0) CALL juDFT_error("error reading byindex,band_min,band_max", calledby="wann_read_inp")
            if(l_p0)write(6,*)"band_min=",wann%band_min(1)
            if(l_p0)write(6,*)"band_max=",wann%band_max(1)
            if(wann%band_min(2).eq.-1)then
               wann%band_min(2)=wann%band_min(1)
               wann%band_max(2)=wann%band_max(1)
            endif
         elseif(trim(task).eq.'byindex2')then
            wann%l_byindex=.true.
            backspace(916)
            read(916,*,iostat=ios)task,wann%band_min(2),wann%band_max(2)
            if (ios /= 0) CALL juDFT_error ("error reading byindex2,band_min2,band_max", calledby="wann_read_inp")
            if(l_p0)write(6,*)"band_min2=",wann%band_min(2)
            if(l_p0)write(6,*)"band_max2=",wann%band_max(2)
            if(wann%band_min(1).eq.-1)then
               wann%band_min(1)=wann%band_min(2)
               wann%band_max(1)=wann%band_max(2)
            endif
         elseif(trim(task).eq.'ikptstart')then
            wann%l_ikptstart=.true.
            backspace(916)
            read(916,*,iostat=ios)task,wann%ikptstart
            if (ios /= 0) CALL juDFT_error ("error reading ikptstart", calledby="wann_read_inp")
            if(l_p0)write(6,*)"ikptstart=",wann%ikptstart
         else
            write(6,*)"unrecognized key: ",task
            CALL juDFT_error ("unrecognized key in wann_inp", calledby="wann_read_inp")
         endif
      enddo

      if (ios /= 0) CALL juDFT_error ("error reading wann_inp", calledby="wann_read_inp")
      if(l_p0.and.ios.lt.0)write(6,*)"end of wann_inp reached"
      close(916)

   ELSE IF (input%l_inpXML) THEN

      DO i = 1, SIZE(wann%jobList)
         task = TRIM(ADJUSTL(wann%jobList(i)))
         if(l_p0) write(6,*)"task ",i,":",task
         if(task(1:1).eq.'!')cycle
         if(index(task,'=F').ne.0.or.index(task,'=f').ne.0)cycle
         read(task,*,iostat=ios)task
         if(index(task,'=').ne.0)then
            task=task(1:(index(task,'=')-1))
         endif
         if(l_p0) write(6,*)"==>key: ",task
         if(trim(task).eq.'endjobs') THEN
            EXIT
         elseif(trim(task).eq.'nabla')then
            wann%l_nabla=.true.
!         elseif(trim(task).eq.'kptsreduc2')then
!            wann%l_kptsreduc2=.true.
!            backspace(916)
!            read(916,*,iostat=ios)task,wann%mhp(1),wann%mhp(2),wann%mhp(3)
!            if (ios /= 0) CALL juDFT_error ("error reading mhp", calledby="wann_read_inp")
!            if(l_p0)write(6,*)"mhp=",wann%mhp(1),wann%mhp(2),wann%mhp(3)
         elseif(trim(task).eq.'nablars')then
            wann%l_nablars=.true.
         elseif(trim(task).eq.'nablapaulirs')then
            wann%l_nablapaulirs=.true.
         elseif(trim(task).eq.'socspicom')then
            wann%l_socspicom=.true.
         elseif(trim(task).eq.'orbitalmom')then
            wann%l_orbitalmom=.true.
         elseif(trim(task).eq.'anglmom')then
            wann%l_anglmom=.true.
         elseif(trim(task).eq.'ndegen')then
            wann%l_ndegen=.true.
         elseif(trim(task).eq.'unformatted')then
            wann%l_unformatted=.true.
         elseif(trim(task).eq.'orbcomp')then
            wann%l_orbcomp=.true.
         elseif(trim(task).eq.'orbcomprs')then
            wann%l_orbcomprs=.true.
         elseif(trim(task).eq.'denmat')then
            wann%l_denmat=.true.
         elseif(trim(task).eq.'nedrho')then
            wann%l_nedrho=.true.
         elseif(trim(task).eq.'perturb')then
            wann%l_perturb=.true.
         elseif(trim(task).eq.'perturbrs')then
            wann%l_perturbrs=.true.
         elseif(trim(task).eq.'rmat')then
            wann%l_rmat=.true.
         elseif(trim(task).eq.'anglmomrs')then
            wann%l_anglmomrs=.true.
         elseif(trim(task).eq.'socspicomrs')then
            wann%l_socspicomrs=.true.
         elseif(trim(task).eq.'she')then
            wann%l_she=.true.
         elseif(trim(task).eq.'berry')then
            wann%l_berry=.true.
         elseif(trim(task).eq.'offdiposop')then
            wann%l_offdiposop=.true.
         elseif(trim(task).eq.'offdiposoprs')then
            wann%l_offdiposoprs=.true.
         elseif(trim(task).eq.'spindisp')then
            wann%l_spindisp=.true.
         elseif(trim(task).eq.'spindisprs')then
            wann%l_spindisprs=.true.
         elseif(trim(task).eq.'torque')then
            wann%l_torque=.true.
         elseif(trim(task).eq.'torquers')then
            wann%l_torquers=.true.
         elseif(trim(task).eq.'perpmag')then
            wann%l_perpmag=.true.
         elseif(trim(task).eq.'perpmagrs')then
            wann%l_perpmagrs=.true.
         elseif(trim(task).eq.'perpmagat')then
            wann%l_perpmagat=.true.
         elseif(trim(task).eq.'perpmagatrs')then
            wann%l_perpmagatrs=.true.
         elseif(trim(task).eq.'mmn0_unf_to_spn_unf')then
            wann%l_mmn0_unf_to_spn_unf=.true.
         elseif(trim(task).eq.'mmn0_to_spn_unf')then
            wann%l_mmn0_to_spn_unf=.true.
         elseif(trim(task).eq.'mmn0_to_spn')then
            wann%l_mmn0_to_spn=.true.
         elseif(trim(task).eq.'mmn0_to_spn2')then
            wann%l_mmn0_to_spn2=.true.
         elseif(trim(task).eq.'mmn0_unf_to_spn')then
            wann%l_mmn0_unf_to_spn=.true.
         elseif(trim(task).eq.'perpmag_unf_to_tor_unf')then
            wann%l_perpmag_unf_to_tor_unf=.true.
         elseif(trim(task).eq.'perpmag_to_tor_unf')then
            wann%l_perpmag_to_tor_unf=.true.
         elseif(trim(task).eq.'perpmag_to_tor')then
            wann%l_perpmag_to_tor=.true.
         elseif(trim(task).eq.'perpmag_unf_to_tor')then
            wann%l_perpmag_unf_to_tor=.true.
         elseif(trim(task).eq.'hsomtxvec_unf_to_lmpzsoc_unf')then
            wann%l_hsomtxvec_unf_to_lmpzsoc_unf=.true.
         elseif(trim(task).eq.'hsomtxvec_to_lmpzsoc_unf')then
            wann%l_hsomtxvec_to_lmpzsoc_unf=.true.
         elseif(trim(task).eq.'hsomtxvec_to_lmpzsoc')then
            wann%l_hsomtxvec_to_lmpzsoc=.true.
         elseif(trim(task).eq.'hsomtxvec_unf_to_lmpzsoc')then
            wann%l_hsomtxvec_unf_to_lmpzsoc=.true.  
         elseif(trim(task).eq.'hsomtx_unf_to_hsoc_unf')then
            wann%l_hsomtx_unf_to_hsoc_unf=.true.
         elseif(trim(task).eq.'hsomtx_to_hsoc_unf')then
            wann%l_hsomtx_to_hsoc_unf=.true.
         elseif(trim(task).eq.'hsomtx_to_hsoc')then
            wann%l_hsomtx_to_hsoc=.true.
         elseif(trim(task).eq.'hsomtx_unf_to_hsoc')then
            wann%l_hsomtx_unf_to_hsoc=.true.
         elseif(trim(task).eq.'socmat')then
            wann%l_socmat=.true.
         elseif(trim(task).eq.'socmatvec')then
            wann%l_socmatvec=.true.
         elseif(trim(task).eq.'socmatrs')then
            wann%l_socmatrs=.true.
         elseif(trim(task).eq.'soctomom')then
            wann%l_soctomom=.true.
         elseif(trim(task).eq.'surfcurr')then
            wann%l_surfcurr=.true.
         elseif(trim(task).eq.'lapw_kpts')then
            wann%l_lapw_kpts=.true.
         elseif(trim(task).eq.'updown')then
            wann%l_updown=.true.
         elseif(trim(task).eq.'stopopt')then
            wann%l_stopopt=.true.
         elseif(trim(task).eq.'projgen')then
            wann%l_projgen=.true.
         elseif(trim(task).eq.'kpointgen')then
            wann%l_kpointgen=.true.
         elseif(trim(task).eq.'potmat')then
            wann%l_potmat=.true.
         elseif(trim(task).eq.'w90kpointgen')then
            wann%l_w90kpointgen=.true.
!         elseif(trim(task).eq.'lapw_gfleur')then
!            wann%l_lapw_gfleur=.true.
!            backspace(916)
!            read(916,*,iostat=ios)task,wann%gfthick,wann%gfcut
!            if (ios /= 0) CALL juDFT_error ("error reading gfcut", calledby="wann_read_inp")
!            if(l_p0)write(6,*)"gfcut=",wann%gfthick,wann%gfcut
!         elseif(trim(task).eq.'lapw')then
!            wann%l_lapw=.true.
!            backspace(916)
!            read(916,*,iostat=ios)task,wann%unigrid(:)
!            if (ios /= 0) CALL juDFT_error ("error reading unigrid", calledby="wann_read_inp")
!            if(l_p0)write(6,*)"unigrid=",wann%unigrid(:)
         elseif(trim(task).eq.'plot_lapw')then
            wann%l_plot_lapw=.true.
         elseif(trim(task).eq.'bzsym')then
            wann%l_bzsym=.true.
         elseif(trim(task).eq.'mmn0')then
            wann%l_mmn0=.true.
         elseif(trim(task).eq.'mmn0at')then
            wann%l_mmn0at=.true.
         elseif(trim(task).eq.'manyfiles')then
            wann%l_manyfiles=.true.
         elseif(trim(task).eq.'collectmanyfiles')then
            wann%l_collectmanyfiles=.true.
         elseif(trim(task).eq.'bestproj')then
            wann%l_bestproj=.true.
         elseif(trim(task).eq.'pauli')then
            wann%l_pauli=.true.
         elseif(trim(task).eq.'pauliat')then
            wann%l_pauliat=.true.
         elseif(trim(task).eq.'proj_plot')then
            wann%l_proj_plot=.true.
         elseif(trim(task).eq.'hopping')then
            wann%l_hopping=.true.
         elseif(trim(task).eq.'plot_symm')then
            wann%l_plot_symm=.true.
         elseif(trim(task).eq.'kptsreduc')then
            wann%l_kptsreduc=.true.
         elseif(trim(task).eq.'fermi')then
            wann%l_fermi=.true.
         elseif(trim(task).eq.'prepwan90')then
            wann%l_prepwan90=.true.
         elseif(trim(task).eq.'plot_umdat')then
            wann%l_plot_umdat=.true.
         elseif(trim(task).eq.'wann_plot')then
            wann%l_wann_plot=.true.
         elseif(trim(task).eq.'bynumber')then
            wann%l_bynumber=.true.
         elseif(trim(task).eq.'matrixmmn')then
            wann%l_matrixmmn=.true.
         elseif(trim(task).eq.'projmethod')then
            wann%l_projmethod=.true.
         elseif(trim(task).eq.'matrixamn')then
            wann%l_matrixamn=.true.
         elseif(trim(task).eq.'wannierize')then
            wann%l_wannierize=.true.
         elseif(trim(task).eq.'plotw90')then
            wann%l_plotw90=.true.
         elseif(trim(task).eq.'dipole')then
            wann%l_dipole=.true.
         elseif(trim(task).eq.'dipole3')then
            wann%l_dipole3=.true.
         elseif(trim(task).eq.'ldauwan')then
            wann%l_ldauwan=.true.
         elseif(trim(task).eq.'byenergy')then
            wann%l_byenergy=.true.
         elseif(trim(task).eq.'finishnocoplot') then
            wann%l_finishnocoplot=.true.
         elseif(trim(task).eq.'finishgwf') then
            wann%l_finishgwf=.true.
         elseif(trim(task).eq.'skipkov') then
            wann%l_skipkov=.true.
         elseif(trim(task).eq.'matrixuhu') then
            wann%l_matrixuHu=.true.
         elseif(trim(task).eq.'matrixuhu-dmi') then
            wann%l_matrixuHu_dmi=.true.
!         elseif(trim(task).eq.'wan90version')then
!            backspace(916)
!            read(916,*,iostat=ios)task,version_real
!            if (ios /= 0) CALL judft_error("error reading wan90version", calledby="wann_read_inp")
!            if(abs(version_real-1.1).lt.1.e-9)then
!               wann%wan90version=1
!            elseif(abs(version_real-1.2).lt.1.e-9)then
!               wann%wan90version=2
!            elseif(abs(version_real-2.0).lt.1.e-9)then
!               wann%wan90version=3
!            else
!              CALL judft_error ("chosen w90 version unknown", calledby="wann_read_inp")
!            endif
!         elseif(trim(task).eq.'atomlist')then
!            wann%l_atomlist=.true.
!            backspace(916)
!            read(916,*,iostat=ios)task,wann%atomlist_num
!            if (ios /= 0) CALL judft_error ("error reading atomlist_num", calledby="wann_read_inp")
!            if(allocated(wann%atomlist))deallocate(wann%atomlist)
!            allocate(wann%atomlist(wann%atomlist_num))
!            backspace(916)
!            read(916,*,iostat=ios)task,wann%atomlist_num,wann%atomlist
!            if (ios /= 0) CALL judft_error ("error reading atomlist", calledby="wann_read_inp")
!            if(l_p0)write(6,*)"atomlist_num=",wann%atomlist_num
!            if(l_p0)write(6,*)"atomlist=",wann%atomlist
!         elseif(trim(task).eq.'ikptstart')then
!            wann%l_ikptstart=.true.
!            backspace(916)
!            read(916,*,iostat=ios)task,wann%ikptstart
!            if (ios /= 0) CALL juDFT_error ("error reading ikptstart", calledby="wann_read_inp")
!            if(l_p0)write(6,*)"ikptstart=",wann%ikptstart
         else
            write(6,*)"unrecognized key: ",task
            CALL juDFT_error ("unrecognized key in wannier jobList", calledby="wann_read_inp")
         endif
      enddo

      IF (wann%l_byindex) THEN
         if(l_p0)write(6,*)"band_min1=",wann%band_min(1)
         if(l_p0)write(6,*)"band_max1=",wann%band_max(1)
         if(l_p0)write(6,*)"band_min2=",wann%band_min(2)
         if(l_p0)write(6,*)"band_max2=",wann%band_max(2)
      END IF

   END IF ! l_file


!-----input file for orbital decomposition
   if(wann%l_orbcomp.or.wann%l_orbcomprs)then
      inquire(file='orbcomp_inp',exist=l_orbcompinp)
      if(l_orbcompinp)then
         open(159,file='orbcomp_inp')
         read(159,*)wann%oc_num_orbs,wann%l_oc_f
         if(allocated(wann%oc_orbs))deallocate(wann%oc_orbs)
         allocate(wann%oc_orbs(wann%oc_num_orbs))
         do n=1,wann%oc_num_orbs
            read(159,*)wann%oc_orbs(n)
         enddo
         close(159)
      else !default is all atoms including f
!         wann%oc_num_orbs=natd
         wann%l_oc_f=.true.
         if(allocated(wann%oc_orbs))deallocate(wann%oc_orbs)
         allocate(wann%oc_orbs(wann%oc_num_orbs))
         do n=1,wann%oc_num_orbs
            wann%oc_orbs(n)=n
         enddo
      endif
   endif

!-----default atom list: all atoms
   if(.not.wann%l_atomlist)then
     if(allocated(wann%atomlist))deallocate(wann%atomlist)
     allocate(wann%atomlist(wann%atomlist_num))
     do n=1,wann%atomlist_num
       wann%atomlist(n)=n
     enddo
   endif      

end subroutine wann_read_inp

end module m_wann_read_inp
