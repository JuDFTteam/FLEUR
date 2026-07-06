!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_iodop 
      USE hdf5 
      USE m_hdf_tools
      use m_juDFT 
      IMPLICIT NONE
!***********************************************************************
!     This subroutine contains subroutines for the IO of potentials for 
!      gfleur                                                           
!     This version uses the hdf data format for the gf_pot-file!        
!                                  Daniel Wortmann, Juelich, 2002       
!***********************************************************************
      PRIVATE 
      INTEGER,PARAMETER,PUBLIC      :: GF_POTFILE = 1 
      INTEGER,PARAMETER,PUBLIC      :: GF_CDNFILE = 2 
      INTEGER,PARAMETER,PUBLIC      :: GF_CDNSTARTFILE = 3 
      INTEGER,PARAMETER,PUBLIC      :: GF_CDNDIFFFILE = 4 
                                                                        
      INTEGER(hid_t),TARGET        :: cdn_id = 0 
      INTEGER(hid_t),TARGET        :: pot_id = 0 
      integer,save                 :: iodop_subcom=-1
                                                                        
                                                                        
      PUBLIC :: gf_loddop,gf_wrtdop,gf_wrtcoul,gf_loddop_name                          &
     &     ,gf_lodcoul,gf_renamepot,gf_iodop_readpseudo                 &
     &     ,gf_iodop_writepseudo,gf_iodop_writevacuum,gf_iodop_readvacuum

      CONTAINS 


      subroutine gf_iodop_readvacuum(file_type,vxy,vz,evac,old)
      ! Read the vacuum potential from potential file
      USE m_hdf_tools
      implicit none
      integer,intent(in)  :: file_type
      real,intent(out)::vz(:,:)
      complex,intent(out)::vxy(:,:,:)
      real,optional,intent(out)::evac
      logical,optional,intent(in):: old
      !locals
      INTEGER(HID_T) :: fid,gid,did
      INTEGER        :: dims(4),hdferr
      logical        :: oldfile

      oldfile=.false.
      if (present(old)) oldfile=old

      fid = priv_openfile(file_type,old=oldfile)
      CALL io_gopen(fid,"vacuum",gid,hdferr)
      if (hdferr/=0) then
             call juDFT_warn("Not a potential file for a surface calculation")
             vxy=0.0
             vz=0.0
             evac=0.0
             return
      endif
      !read vz
      CALL io_dopen(gid,"vz",did,hdferr)
      CALL io_datadim(did,dims(1:2))
      if (dims(1)/=size(vz,1)) call juDFT_error("Wrong dimension in vz")
      call io_read(did,(/1,1/),(/dims(1),1/),vz(:,1))
      if (size(vz,2)>1) THEN
      	if (dims(2)>1) THEN
      		call io_read(did,(/1,2/),(/dims(1),1/),vz(:,2))!read second spin
      		if (size(vz,2)>2) THEN
      		    !now check noco-part
      		    if (dims(2)>2) then
      		        call io_read(did,(/1,3/),(/dims(1),2/),vz(:,3:4))
      		    else
      		        vz(:,3:4)=0.0
      		    endif
      		endif
        else
            vz(:,2)=vz(:,1)
            if (size(vz,2)>2) vz(:,3:4)=0.0
            if (file_type==GF_CDNFILE) vz=vz/2.0
        endif
      endif
      CALL io_dclose(did,hdferr)
      !read vxy
      CALL io_dopen(gid,"vxy",did,hdferr)
      CALL io_datadim(did,dims)
      if (any(dims(2:3)/=(/size(vxy,1),size(vxy,2)/))) then
      		write(*,*) dims(2:3)
      		write(*,*)  size(vxy,1),size(vxy,2)
            call juDFT_error("Wrong dimension in vxy")
      endif
      dims(1)=1
      CALL io_READ(did,(/-1,1,1,1/),(/1,dims(2),dims(3),1/),vxy(:,:,1))
      if (size(vxy,3)>1) THEN
         if (dims(4)>1) then
         	CALL io_READ(did,(/-1,1,1,2/),(/1,dims(2),dims(3),1/),vxy(:,:,2))
         	if (size(vxy,3)>2) then !read noco-part
         	    if (dims(4)>2) THEN
         	      CALL io_READ(did,(/-1,1,1,3/),(/1,dims(2),dims(3),1/),vxy(:,:,3))
         	    else
         	      vxy(:,:,3)=0.0
         	    endif
         	endif
         else
            vxy(:,:,2)=vxy(:,:,1)
            if (size(vxy,3)>2) vxy(:,:,3)=0.0
            if (file_type==GF_CDNFILE) vxy=vxy/2.
         endif
      endif
      CALL io_dclose(did,hdferr)

	  if (present(evac)) call io_read_att(gid,"evac",evac)

      CALL io_gclose(gid,hdferr)
      CALL  priv_closefile(file_type)
      end subroutine


      subroutine gf_iodop_writevacuum(file_type,vxy,vz,mpi_subcom,evac)
      ! Write the vacuum potential to the potential file
      USE m_hdf_tools
      use hdf5
      implicit none
      integer,intent(in)::file_type
      real,intent(in)     ::vz(:,:)
      complex,intent(in)  ::vxy(:,:,:)
      real,intent(in),optional::evac
	  INTEGER             :: mpi_subcom
      !locals
      INTEGER(HID_T) :: fid,gid,dvz,dvxy
      INTEGER        :: hdferr


      fid = priv_openfile(file_type,.FALSE.,mpi_subcom)
      IF (io_groupexists(fid,"vacuum")) THEN
         CALL io_gopen(fid,"vacuum",gid,hdferr)
         call io_dopen(gid,"vz",dvz,hdferr)
         call io_dopen(gid,"vz",dvxy,hdferr)
      ELSE
      	CALL io_gcreate(fid,"vacuum",gid,hdferr)
      	call io_createvar(gid,"vz", H5T_NATIVE_DOUBLE,(/size(vz,1),size(vz,2)/),dvz)
      	call io_createvar(gid,"vxy", H5T_NATIVE_DOUBLE,(/2,size(vxy,1),size(vxy,2),size(vxy,3)/),dvxy)
      endif
      call io_write(dvz,(/1,1/),(/size(vz,1),size(vz,2)/),vz)
      call io_write(dvxy,(/-1,1,1,1/),(/1,size(vxy,1),size(vxy,2),size(vxy,3)/),vxy)
      if (present(evac)) call io_write_att(gid,"evac",evac)
      CALL io_dclose(dvz,hdferr)
      CALL io_dclose(dvxy,hdferr)

      CALL io_gclose(gid,hdferr)
	  CALL  priv_closefile(file_type,mpi_subcom)
      end subroutine
                                                                        
      !<-- S:gf_iodop_readpseudo(layer,psq)                             
                                                                        
      SUBROUTINE gf_iodop_readpseudo(layer,psq) 
!-----------------------------------------------                        
!                                                                       
!           (last modified:08-07-10) D. Wortmann                        
!-----------------------------------------------                        
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)              :: layer 
      COMPLEX,ALLOCATABLE,INTENT(OUT) :: psq(:,:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: fid,gid,did 
      INTEGER        :: dims(3),hdferr 
      !>                                                                
      fid = priv_openfile(GF_CDNFILE,.FALSE.) 
      CALL io_gopen(fid,io_layername(layer),gid,hdferr) 
      CALL io_dopen(gid,"pseudocharge",did,hdferr) 
      CALL io_datadim(did,dims) 
      ALLOCATE(psq(dims(2),dims(3))) 
      dims(1) = 1 
      CALL io_READ(did,(/-1,1,1/),dims,psq) 
      CALL io_dclose(did,hdferr) 
      CALL io_gclose(gid,hdferr) 
      CALL  priv_closefile(GF_CDNFILE)
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:gf_iodop_writepseudo(layer,psq,mpi_subcom)                 
                                                                        
      SUBROUTINE gf_iodop_writepseudo(layer,psq,mpi_subcom) 
!-----------------------------------------------                        
!                                                                       
!           (last modified:08-07-10) D. Wortmann                        
!-----------------------------------------------                        
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)             :: layer 
      COMPLEX,INTENT(IN)             :: psq(:,:) 
      INTEGER                        :: mpi_subcom 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: fid,gid,did 
      INTEGER        :: hdferr 
      !>                                                                
      fid = priv_openfile(GF_CDNFILE,.FALSE.,mpi_subcom) 
      CALL io_gopen(fid,io_layername(layer),gid,hdferr) 
      IF (io_dataexists(gid,"pseudocharge")) THEN 
         CALL io_dopen(gid,"pseudocharge",did,hdferr) 
      ELSE 
         CALL io_createvar(gid,"pseudocharge", H5T_NATIVE_DOUBLE,(/2    &
     &        ,SIZE(psq,1),SIZE(psq,2)/),did)                           
      ENDIF 
      CALL io_WRITE(did,(/-1,1,1/),(/1,SIZE(psq,1),SIZE(psq,2)/),psq) 
      CALL io_dclose(did,hdferr) 
      CALL io_gclose(gid,hdferr) 
      CALL  priv_closefile(GF_CDNFILE,mpi_subcom) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:gf_loddop(file_type,layer,jspins,atoms,stars,sphhar,vr,vpw,
                                                                        
      SUBROUTINE gf_loddop(file_type,layer,jspins,                      &
     &     atoms,stars,sphhar,                                          &
     &     vr,vpw,noco,old)
!***********************************************************************
!     This subroutine reads in the potential or the charge density      
!     needed by the green fct part of FLEUR.                            
!                                                                       
!                                                                       
!                                  Daniel Wortmann, Juelich, 2001       
!***********************************************************************
      USE m_gf_types 
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)       :: file_TYPE 
      INTEGER,       INTENT(IN) :: jspins,layer 
      TYPE(t_atoms), INTENT(IN)::atoms 
      TYPE(t_stars), INTENT(IN)::stars 
      TYPE(t_sphhar),INTENT(IN)::sphhar 
      type(t_noco),intent(in),optional  :: noco
      REAL,          INTENT(INOUT):: vr(:,0:,:,:) 
      COMPLEX,       INTENT(INOUT):: vpw(:,:) 
      LOGICAL ,      INTENT(IN),optional   :: old
                                                                        
! locals                                                                
!     scalars                                                           
!     the hdf-IDs                                                       
      INTEGER(HID_T)::fid,pwid,pwodid,vrid,gid,ggid 
      INTEGER::hdferr,dims(3)
!     needed for input checking                                         
      INTEGER::ntype_in,n3d_in,jrid_in,nlhd_in,jspins_in
                                                                        
      !for old charge density                                           
      INTEGER :: iter 
      CHARACTER(LEN = 7) ::gname 

      !we might want to rotate the spin
      real :: theta,phi
      complex :: u(2,2),ut(2,2),tmp(2,2)
      integer :: n

                                                                        
!     temporary array for interstitial pot                              
      REAL, ALLOCATABLE:: pwd_in(:,:,:) 
      REAL,ALLOCATABLE :: vr_tmp(:,:,:,:) 
!                                                                       
!  lmax might have been changed!                                        
!                                                                       
      vr(:,0:,:,:)=0.0 
!                                                                       
!                                                                       
!  OPEN FILE and get ID's                                               
      fid = priv_openfile(file_type,old) 
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr) 
      CALL io_gopen (ggid,"pot",gid,hdferr) 
      CALL io_dopen(gid, 'vr', vrid, hdferr) 
      CALL io_dopen(gid, 'pw', pwid, hdferr) 
      CALL io_check('gf_iodop$gf_loddop:open ',hdferr)

!     check no of spins in file
      call io_datadim(pwid,dims)
      jspins_in=dims(2)
!                                                                       
!                                                                       
!      read MT-part from file                                           
!                                                                       
      CALL io_read_att(vrid,"ntype",ntype_in) 
      IF (ntype_in/=atoms%ntype) THEN 
         WRITE (*,*) 'Wrong no of atoms' 
         CALL juDFT_error('Wrong MT-Potential') 
      ENDIF 
      CALL io_read_att(vrid,"jrid",jrid_in) 
      CALL io_read_att(vrid,"nlhd",nlhd_in) 
      vr=0.0 
      !read spins seperately                                            
      nlhd_in=min(nlhd_in,sphhar%nlhd) 
      jrid_in=min(jrid_in,size(vr,1)) 
      ALLOCATE(vr_tmp(jrid_in,nlhd_in+1,ntype_in,1)) 
      CALL io_read(vrid,(/1,1,1,1/),                                    &
     &     (/jrid_in,nlhd_in+1,ntype_in,1/)                             &
     &     ,vr_tmp)                                                     
      vr(:jrid_in,0:nlhd_in,:ntype_in,1:1)=vr_tmp 
      IF (jspins>1) THEN
         if (jspins_in<2) then
              vr(:jrid_in,0:nlhd_in,:ntype_in,2:2)=vr(:jrid_in,0:nlhd_in,:ntype_in,1:1)
         else
            CALL io_read(vrid,(/1,1,1,2/),                                 &
     &      (/jrid_in,nlhd_in+1,ntype_in,1/)                             &
     &      ,vr_tmp)
            vr(:jrid_in,0:nlhd_in,:ntype_in,2:2)=vr_tmp
         endif
      ENDIF 
      DEALLOCATE(vr_tmp) 
!                                                                       
!     read interstitial part                                            
!                                                                       
      CALL io_read_att(pwid, "n3d", n3d_in) 

      vpw=0.0
                                                                        
      ALLOCATE(pwd_in(n3d_in,1,2))
      CALL io_read(pwid,(/1,1,1/),                                      &
     &        (/n3d_in,1,2/),pwd_in)
      vpw(1:min(n3d_in,stars%nq3),1)=                             &
     &     cmplx(pwd_in(1:min(n3d_in,stars%nq3),1,1),pwd_in(1:min(n3d_in&
     &     ,stars%nq3),1,2))
     if (jspins>1) Then
         if (jspins_in<2) then
              vpw(1:min(n3d_in,stars%nq3),2)=vpw(1:min(n3d_in,stars%nq3),1)
         else
            CALL io_read(pwid,(/1,2,1/),(/n3d_in,1,2/),pwd_in)
            vpw(1:min(n3d_in,stars%nq3),2)=cmplx(pwd_in(1:min(n3d_in,stars%nq3),1,1),pwd_in(1:min(n3d_in,stars%nq3),1,2))
         endif
      endif
      DEALLOCATE(pwd_in)
      !if second spin was missing renormalize data
      if (jspins.ne.jspins_in) then
            call juDFT_warn("Generate spin-polarized calculation from nonpolarized data")
            if (file_type==GF_CDNFILE) THEN
                vpw=vpw/2.
                vr=vr/2.
            endif
      endif
      if(present(noco)) then
      IF (noco%l_noco ) THEN
         IF (io_dataexists(gid,'pwod')) THEN 
            CALL io_dopen(gid, 'pwod', pwodid, hdferr) 
            ALLOCATE(pwd_in(n3d_in,1,2)) 
            CALL io_read(pwodid,(/1,1,1/),                              &
     &           (/n3d_in,1,2/),pwd_in)                                 
            vpw(1:min(n3d_in,stars%nq3),3)=                             &
     &         cmplx(pwd_in(1:min(n3d_in,stars%nq3),1,1),               &
     &          pwd_in(1:min(n3d_in,stars%nq3),1,2))                    
            DEALLOCATE(pwd_in) 
            CALL io_dclose(pwodid,hdferr) 
         ELSE 
            WRITE(6,*) 'Off-diagonal potential matrix could not be      &
     &           read and are now set to zero'                          
         ENDIF
         !now rotate the interstitial charge if needed!!!
         IF (noco%alpha_int.ne.0.0.or.noco%beta_int.ne.0.0) THEN
            call juDFT_warn("Interstitial spin is rotated, only use once in self-consistency")
            theta=noco%beta_int
            phi=noco%alpha_int
            !<-- Now rotate the spin-quantisation axis if needed
            u(1,1) =  EXP(-cmplx(0,1)*phi/2)*COS(theta/2)
            u(1,2) = -EXP(-cmplx(0,1)*phi/2)*SIN(theta/2)
            u(2,1) =  EXP(cmplx(0,1)*phi/2)*SIN(theta/2)
            u(2,2) =  EXP(cmplx(0,1)*phi/2)*COS(theta/2)
            ut     = TRANSPOSE(CONJG(u))
            DO n = 1,size(vpw,1)
               tmp(1,1)=vpw(n,1)
               tmp(2,2)=vpw(n,2)
               tmp(2,1)=vpw(n,3)
               tmp(1,2)=conjg(vpw(n,3))
               tmp = MATMUL(MATMUL(u,tmp),ut)
               vpw(n,1)=tmp(1,1)
               vpw(n,2)=tmp(2,2)
               vpw(n,3)=tmp(2,1)
             ENDDO
         ENDIF 
      ENDIF 
      endif
                                                                        
                                                                        
      CALL priv_checkstarinfo(gid,stars) 
!  close all id's                                                       
      CALL io_dclose(vrid,hdferr) 
      CALL io_dclose(pwid,hdferr) 
      CALL io_gclose(gid,hdferr) 
      CALL io_gclose(ggid,hdferr) 
      CALL io_check('gf_iodop$gf_loddop:close ',hdferr) 
      CALL priv_closefile(file_type) 
      END SUBROUTINE gf_loddop 
                                                                        
      !>                                                                
                                                                        
                                                                        
                                                                        
      !<-- S:gf_wrtdop(file_type,jspins,gfinp,atoms,stars,sphhar,vr,vpw)
                                                                        
      SUBROUTINE gf_wrtdop(file_type,layer,jspins,                      &
     &     gfinp,atoms,stars,sphhar,                                    &
     &     vr,vpw,l_noco,mpi_subcom,l_valence_only)
!***********************************************************************
!     This subroutine writes in the potential or the charge density     
!     needed by the green fct part of FLEUR.                            
!                                                                       
!                                                                       
!     added writing of planar averaged interstitial charges             
!                                                                       
!                                  Daniel Wortmann, Juelich, 2001       
!***********************************************************************
      USE hdf5 
      USE m_gf_types 
      USE m_hdf_tools 
      USE m_gf_stars 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)        :: file_type 
      INTEGER,INTENT(IN)        :: layer 
      INTEGER,       INTENT(IN) :: jspins 
      TYPE(t_gfinp), INTENT(IN) :: gfinp 
      TYPE(t_atoms), INTENT(IN) :: atoms 
      TYPE(t_stars), INTENT(IN) :: stars 
      TYPE(t_sphhar),INTENT(IN) :: sphhar 
      INTEGER,INTENT(IN)        :: mpi_subcom 
      REAL,          INTENT(IN) :: vr(:,0:,:,:) 
      COMPLEX,       INTENT(IN) :: vpw(:,:) 
      LOGICAL, INTENT(IN)       :: l_noco
      logical,intent(in),optional:: l_valence_only
                                                                        
! locals                                                                
                                                                        
!     the hdf-IDs                                                       
      INTEGER(HID_T)::fid,vrid,pwid,pwodid,gid,avpwid,ffid 
      INTEGER       ::vrdim(4),pwdim(3),avdim(2) 
      INTEGER       ::pwstart(3),pwcount(3) 
      INTEGER       ::hdferr 
      LOGICAL       ::lexist 
      REAL          :: av(2*stars%mx3+1,jspins+1)
      logical       :: l_valence

      l_valence=.false.
      if (present(l_valence_only)) l_valence=l_valence_only
      write(*,*) file_type,layer,l_valence
!  OPEN FILE and get ID's                                               
      ffid=priv_openfile(file_type,mpi_subcom=mpi_subcom) 
      IF (io_groupexists(ffid,io_layername(layer))) THEN 
         CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
      ELSE 
         CALL io_gcreate(ffid,io_layername(layer),fid,hdferr) 
      ENDIF 
! create Dataspaces+Datasets
      if (io_groupexists(fid,"pot")) THEN
         CALL io_gopen(fid,"pot",gid,hdferr)
      ELSE
         CALL io_gcreate(fid,"pot",gid,hdferr)
      ENDIF

      vrdim = (/maxval(atoms%jri),maxval(sphhar%nlh)+1,atoms%ntype      &
     &     ,jspins/)
      IF (l_valence) THEN
        CALL io_createvar(gid,"vr_valence", H5T_NATIVE_DOUBLE,vrdim,vrid)
      else
        CALL io_createvar(gid,"vr", H5T_NATIVE_DOUBLE,vrdim,vrid)
                                                                        
        pwdim=(/ stars%nq3,jspins,2 /)
        CALL io_createvar(gid,"pw", H5T_NATIVE_DOUBLE,pwdim,pwid)
                                                                        
        avdim = (/ 2*stars%mx3+3,jspins+1 /)
        CALL io_createvar(gid,"pwav", H5T_NATIVE_DOUBLE,avdim,avpwid)
                                                                        
        IF (l_noco) THEN
! off-diagonalpart                                                      
           pwdim=(/ stars%nq3,1,2 /)
           CALL io_createvar(gid,"pwod", H5T_NATIVE_DOUBLE,pwdim,pwodid)
        ENDIF
      endif
                                                                        
!    save all attributes                                                
!                                                                       
      CALL io_write_att(vrid, "ntype",atoms%ntype) 
      CALL io_write_att(vrid, "jrid",maxval(atoms%jri)) 
      CALL io_write_att(vrid, "nlhd",maxval(sphhar%nlh)) 
      if (.not.l_valence) CALL io_write_att(pwid, "n3d",stars%nq3)
                                                                        
!                                                                       
!      write MT-part to file                                            
!                                                                       
                                                                        
      ! Write the spins seperately                                      
      vrdim(4)=1 
      CALL io_write(vrid,(/1,1,1,1/),vrdim,                             &
     &     vr(:vrdim(1),0:vrdim(2)-1,:vrdim(3),1:1))                    
      IF (jspins>1) THEN 
         CALL io_write(vrid,(/1,1,1,2/),vrdim,                          &
     &     vr(:vrdim(1),0:vrdim(2)-1,:vrdim(3),2:2))                    
      ENDIF 


      if (.not.l_valence) then
!                                                                       
!      interstitial part                                                
!                                                                       
      pwstart=(/1,1,1/) 
      pwcount=(/ stars%nq3,jspins,1 /) 
                                                                        
      CALL io_write(pwid,pwstart,pwcount,real(vpw(:stars%nq3,:jspins))) 
      pwstart=(/1,1,2/) 
      CALL io_write(pwid,pwstart,pwcount,aimag(vpw(:stars%nq3           &
     &     ,:jspins)))                                                  
      IF (l_noco) THEN 
         pwstart=(/1,1,1/) 
         pwcount=(/ stars%nq3,1,1 /) 
         CALL io_write(pwodid,pwstart,pwcount,                          &
     &        real(vpw(:stars%nq3,3)))                                  
         pwstart=(/1,1,2/) 
         CALL io_write(pwodid,pwstart,pwcount,                          &
     &      aimag(vpw(:stars%nq3,3)))                                   
      ENDIF 
      CALL priv_writestarinfo(gid,stars) 
                                                                        
      !<-- Write planar av. of charge                                   
      av = priv_getPlanar(stars,vpw) 
      CALL io_write(avpwid,(/1,1/),shape(av),av)

      !>                                                                
!                                                                       
!  close all id's                                                       
      CALL io_dclose(pwid,hdferr) 
      CALL io_dclose(avpwid,hdferr) 
      IF (l_noco)  CALL io_dclose(pwodid,hdferr) 
      endif
      CALL io_dclose(vrid,hdferr)
      CALL io_gclose(gid,hdferr) 
                                                                        
      !<-- Save the stars if not in the file                            
      IF (.NOT.io_groupexists(fid,"stars")) THEN 
         CALL gf_WRITE_stars(fid,stars) 
      ENDIF 
      !>                                                                
      CALL io_gclose(fid,hdferr) 
      CALL io_check('gf_iodop%gf_wrtdop ',hdferr) 
      CALL priv_closefile(file_type,mpi_subcom) 
      END SUBROUTINE gf_wrtdop 
                                                                        
      !>                                                                
                                                                        
      !<-- S:priv_writestarinfo(gid,vpw)                                
                                                                        
      SUBROUTINE priv_writestarinfo(gid,stars) 
!******************************************                             
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_gf_types 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN) :: gid 
      TYPE(t_stars),INTENT(IN)  :: stars 
      !>                                                                
      !<--Locals                                                        
      INTEGER(HID_t)      :: vid 
      INTEGER             :: nx,ny,nz,hdferr 
      !>                                                                
      !<--open dataset+write attributes                                 
      nx=stars%mx1*2+1 
      ny=stars%mx2*2+1 
      nz=stars%mx3*2+1 
      CALL io_createvar(gid,"starinfo", H5T_NATIVE_INTEGER,(/nx,ny,nz/) &
     &     ,vid)                                                        
      CALL io_write_att(vid, "mx1",stars%mx1) 
      CALL io_write_att(vid, "mx2",stars%mx2) 
      CALL io_write_att(vid, "mx3",stars%mx3) 
      !>                                                                
      CALL io_WRITE(vid,(/1,1,1/),(/nx,ny,nz/),                         &
     &     stars%ig(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,          &
     &     -stars%mx3:stars%mx3))                                       
      CALL io_dclose(vid,hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:priv_checkstarinfo(gid,stars)                              
      SUBROUTINE priv_checkstarinfo(gid,stars) 
!******************************************                             
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_gf_types 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER(HID_T),INTENT(IN) :: gid 
      TYPE(t_stars),INTENT(IN)  :: stars 
      !>                                                                
      !<--Locals                                                        
      INTEGER(HID_t)      :: vid 
      INTEGER             :: nx,ny,nz,x,y,z,hdferr 
      INTEGER,ALLOCATABLE :: ig(:,:,:) 
      !>                                                                
      !<--open dataset+write attributes                                 
      CALL io_dopen(gid,"starinfo",vid,hdferr) 
      IF (hdferr/=0) THEN 
         WRITE(6,*) 'No starinfo found in gf_pot' 
         WRITE(6,*) 'Potential expandion into stars not checked!' 
         WRITE(*,*) 'Warning: starcheck in gf_iodop failed' 
         RETURN 
      ENDIF 
                                                                        
      CALL io_read_att(vid, "mx1",nx) 
      CALL io_read_att(vid, "mx2",ny) 
      CALL io_read_att(vid, "mx3",nz) 
      !>                                                                
      ALLOCATE(ig(-nx:nx,-ny:ny,-nz:nz)) 
      CALL io_read(vid,(/1,1,1/),(/2*nx+1,2*ny+1,2*nz+1/),              &
     &     ig(-nx:nx,-ny:ny,-nz:nz))                                    
      CALL io_dclose(vid,hdferr) 
      !Check if mapping is ok!                                          
                                                                        
      DO x=-min(stars%mx1,nx),min(stars%mx1,nx) 
         DO y=-min(stars%mx2,ny),min(stars%mx2,ny) 
            DO z=-MIN(stars%mx3,nz),MIN(stars%mx3,nz) 
             IF ((ig(x,y,z)>0).AND.(stars%ig(x,y,z)>0)) THEN 
               IF (ig(x,y,z)/=stars%ig(x,y,z)) THEN
                 write(*,*) "x,y,z",x,y,z
                 write(*,*) ig(x,y,z),stars%ig(x,y,z)
                 CALL juDFT_error(          &
     &              'Mapping error! Check the stars used in gf_pot')
                 endif
             ENDIF 
            ENDDO 
         ENDDO 
      ENDDO 
      DEALLOCATE(ig) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
                                                                        
      !<-- F: priv_getPlanar(stars,pw)results(planar)                   
                                                                        
      FUNCTION priv_getPlanar(stars,pw)RESULT(planar) 
!-----------------------------------------------                        
!     calculate the planar average of the charge/potential in pw        
!             (last modified: 09-11-19) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_fft_singleton 
      USE m_gf_types 
      USE m_constants 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(IN) :: stars 
      COMPLEX,INTENT(IN)       :: pw(:,:) 
      REAL                     :: planar(-stars%mx3:stars%mx3           &
     &     ,SIZE(pw,2)+1)                                               
      !>                                                                
      !<-- Locals                                                       
      COMPLEX :: vz(0:2*stars%mx3) 
      INTEGER :: n,jspin 
      REAL    :: dz 
      !>                                                                
      DO jspin = 1,SIZE(pw,2) 
         vz    = 0.0 
         DO n  =-stars%mx3,stars%mx3 
            IF (stars%ig(0,0,n) == 0) CYCLE 
            IF (n<0) THEN 
               vz(n+2*stars%mx3+1) = pw(stars%ig(0,0,n),jspin) 
            ELSE 
               vz(n)               = pw(stars%ig(0,0,n),jspin) 
            ENDIF 
         ENDDO 
         vz = fft(vz,inv = .TRUE.) 
         DO n  = -stars%mx3,stars%mx3 
            IF (n<0) THEN 
               planar(n,jspin+1) = REAL(vz(n+2*stars%mx3+1)) 
            ELSE 
               planar(n,jspin+1) = REAL(vz(n)) 
            ENDIF 
         ENDDO 
      ENDDO 
                                                                        
      !Add z-coordinate as first row                                    
      dz = stars%sk3(stars%ig(0,0,1)) 
      dz = 2*pimach()/dz/REAL(2*stars%mx3+1) 
      DO n  =-stars%mx3,stars%mx3 
         planar(n,1) = n*dz 
      ENDDO 
                                                                        
      END FUNCTION 
                                                                        
      !>                                                                
                                                                        
      !<-- F: priv_openFile(file.mpi)                                   
                                                                        
      FUNCTION priv_openFile(file,old,mpi_subcom,veryold)
!-----------------------------------------------                        
!  opens/creates the gf_pot.hdf and gf_cdn.hdf files                    
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      USE hdf5 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)              :: file 
      LOGICAL,INTENT(IN),OPTIONAL     :: old,veryold
      INTEGER,INTENT(IN),OPTIONAL     :: mpi_subcom 
      INTEGER(hid_t)                  :: priv_openFile 
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T),POINTER   :: ID 
      CHARACTER(LEN = 16)      :: filename 
      LOGICAL                  :: lexist,o 
      INTEGER                  :: hdferr,mpierr,irank,n
      INTEGER(hid_t)           :: access_prp,rw_prp 
      !>                                                                
#ifdef CPP_MPI                                                          
      INCLUDE 'mpif.h' 
#endif                                                                  
      access_prp = H5P_DEFAULT_f 
      IF (PRESENT(mpi_subcom)) THEN 
        if(iodop_subcom==-1) iodop_subcom=mpi_subcom
#ifdef CPP_HDFMPI
        if(iodop_subcom.ne.mpi_subcom) call juDFT_warn("Parallelisation error in gf_iodop")
#endif
#ifdef CPP_MPI                                                          
         CALL MPI_COMM_RANK(mpi_subcom,irank,mpierr) 
         IF (irank>0) THEN 
            !wait for other PE to finish                                
            CALL MPI_RECV(hdferr,0,MPI_INTEGER,irank-1,1                &
     &           ,mpi_subcom,MPI_STATUS_IGNORE,mpierr)                  
         ENDIF 
#endif                                                                  
         rw_prp = H5F_ACC_RDWR_f 
      ELSE 
         rw_prp = H5F_ACC_RDONLY_F 
      ENDIF 
                                                                        
      o = .FALSE. 
      IF (PRESENT(old)) o = old
      IF (PRESENT(veryold)) o = veryold
                                                                        
      IF (o) THEN 
         IF (present(mpi_subcom)) CALL juDFT_error                        &
     &        ("Can not open old-file for writing")
         n= priv_lastfileNO(file)
         if (present(veryold)) then
           if (veryold) n=n-1
         endif
         IF (file == gf_potfile) THEN 
            ID    => pot_id 
            WRITE(filename,"(a,i0,a)") "gf_pot",n   &
     &           ,".hdf"                                                
         ELSEIF (file == gf_cdnfile) THEN 
            ID        => cdn_id 
            WRITE(filename,"(a,i0,a)") "gf_cdn",n   &
     &           ,".hdf"                                                
         ENDIF 
      ELSE 
         IF (file == gf_potfile) THEN 
            ID       => pot_id 
            filename ="gf_pot.hdf" 
         ELSEIF (file   == gf_cdnfile) THEN 
            ID   => cdn_id 
            filename  ="gf_cdn.hdf" 
         ELSEIF (file   == gf_cdnstartfile) THEN 
            ID   => cdn_id 
            filename  ="gf_cdn.start.hdf" 
         ELSEIF (file   == gf_cdndifffile) THEN 
            ID   => cdn_id 
            filename  ="gf_cdn.diff.hdf" 
         ENDIF 
      ENDIF 
      IF (ID /= 0) CALL juDFT_error("File was opened before??") 
                                                                        
      IF (.NOT.PRESENT(mpi_subcom)) THEN 
         CALL io_hdfopen(filename,rw_prp,id                              &
     &        ,hdferr,access_prp)                                       
      ELSE 
         INQUIRE(FILE = filename,EXIST=lexist) 
         IF (lexist) THEN 
            CALL io_hdfopen(filename,rw_prp,id                           &
     &           ,hdferr,access_prp)                                    
         ELSE 
            WRITE(*,*) "Create file:",filename 
            CALL h5fcreate_f(filename,H5F_ACC_TRUNC_F,id                &
     &           ,hdferr,h5p_default_f,access_prp)                      
         ENDIF 
      ENDIF 
                                                                        
      IF (access_prp /= H5P_DEFAULT_f)                                  &
     &     CALL h5pclose_f(access_prp,hdferr)                           
      priv_openfile=ID 
                                                                        
      END FUNCTION 
                                                                        
      !>                                                                
                                                                        
      !<-- S: priv_closefile(file,mpi_subcom)                           
                                                                        
      SUBROUTINE priv_CLOSEfile(file,mpi_subcom) 
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE hdf5 
      USE m_gf_types 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)               :: file 
      INTEGER,INTENT(IN),OPTIONAL      :: mpi_subcom 
      INTEGER                          :: hdferr,mpierr,irank,isize 
      INTEGER(hid_t),POINTER           :: ID 
#ifdef CPP_MPI                                                          
      include "mpif.h" 
#endif                                                                  
      IF (file == gf_potfile) THEN 
         ID    => pot_id 
      ELSE 
         ID  => cdn_id 
      ENDIF 
      CALL io_hdfclose(id,hdferr) 
      id = 0 
#ifdef CPP_MPI                                                          
      IF (PRESENT(mpi_subcom)) THEN 
        if(iodop_subcom==-1) iodop_subcom=mpi_subcom
        if(iodop_subcom.ne.mpi_subcom) call juDFT_warn("Parallelisation error in gf_iodop")
         CALL MPI_COMM_RANK(mpi_subcom,irank,mpierr) 
         CALL MPI_COMM_SIZE(mpi_subcom,isize,mpierr) 
         IF (irank<isize-1) THEN 
            !send a message to next PE in Communicator                  
            CALL MPI_SEND(1,0,MPI_INTEGER,irank+1,1                     &
     &           ,mpi_subcom,mpierr)                                    
         ENDIF 
      ENDIF 
#endif                                                                  
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
      !Used only fro the construction of the pot-boundary condition     
      !<-- S:gf_wrtcoul(file_type,layer,vpw,qpw)                        
                                                                        
      SUBROUTINE  gf_wrtcoul(file_type,layer,stars,vpw,qpw,vr)
!******************************************                             
!     writes the coulomb potential to the file                          
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_hdf_tools 
      use m_juDFT 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      INTEGER,INTENT(IN)          :: file_type 
      TYPE(t_stars),INTENT(IN)    :: stars 
      INTEGER,INTENT(IN)          :: layer 
      COMPLEX,OPTIONAL,INTENT(IN) :: vpw(:) 
      COMPLEX,OPTIONAL,INTENT(IN) :: qpw(:,:) 
      REAL,INTENT(IN),OPTIONAL    :: vr(:,:,:,:)

      !>                                                                
                                                                        
      !<-- Locals                                                       
      INTEGER(HID_T) :: fid,avid,pwid,gid,ggid,vrid
      INTEGER        :: hdferr,n,dim(3) 
      LOGICAL        :: old 
      REAL           :: av(2*stars%mx3+1,2) 
      !>                                                                
      CALL io_hdfopen("gf_pot.hdf",H5F_ACC_RDWR_F,fid,hdferr) 
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr) 
      IF (hdferr /= 0)  CALL juDFT_error("Could not write Coulpot") 

      IF(io_groupexists(ggid,"coulpot")) THEN
         CALL io_gopen(ggid,"coulpot",gid,hdferr)
      ELSE
         CALL io_gcreate(ggid,"coulpot",gid,hdferr) 
      ENDIF

      !Check if old data is of correct size was disabled!!
      IF (PRESENT(vpw)) THEN 
         IF (io_dataexists(gid,"pwav")) THEN
             CALL io_dopen(gid,"pwav",avid,hdferr)
         ELSE
            CALL io_createvar(gid,"pwav",H5T_NATIVE_DOUBLE,(/2       &
     &              *stars%mx3+1,2/),avid)                              
         ENDIF
         IF (io_dataexists(gid,"coul")) THEN
             CALL io_dopen(gid,"coul",pwid,hdferr)
         ELSE
            CALL io_createvar(gid,"coul",H5T_NATIVE_DOUBLE,(/size(vpw,1),2/),pwid)
         ENDIF 
         CALL io_WRITE(pwid,(/1,1/),(/SIZE(vpw),1/),REAL(vpw)) 
         CALL io_WRITE(pwid,(/1,2/),(/SIZE(vpw),1/),AIMAG(vpw)) 
         av = priv_getPlanar(stars,RESHAPE(vpw,(/SIZE(vpw),1/))) 
                                                                        
         CALL io_WRITE(avid,(/1,1/),SHAPE(av),av) 
         CALL io_dclose(avid,hdferr) 
         CALL io_dclose(pwid,hdferr) 
      ENDIF 
      IF (PRESENT(qpw)) THEN 
         IF (io_dataexists(gid,"charge")) THEN
            CALL io_dopen(gid,"charge",pwid,hdferr) 
         ELSE 
            CALL io_createvar(gid,"charge", H5T_NATIVE_DOUBLE,(/SIZE(qpw&
     &           ,1),SIZE(qpw,2),2/),pwid)                              
         ENDIF 
         CALL io_WRITE(pwid,(/1,1,1/),(/SIZE(qpw,1),SIZE(qpw,2),1/)     &
     &        ,REAL(qpw))                                               
         CALL io_WRITE(pwid,(/1,1,2/),(/SIZE(qpw,1),SIZE(qpw,2),1/)     &
     &        ,AIMAG(qpw))                                              
         CALL io_dclose(pwid,hdferr) 
      ENDIF 
      IF (PRESENT(vr)) THEN
         IF (io_dataexists(gid,"vr")) THEN
            CALL io_dopen(gid,"vr",vrid,hdferr)
         ELSE
            CALL io_createvar(gid,"vr", H5T_NATIVE_DOUBLE,(/SIZE(vr,1),  &
     &       SIZE(vr,2),size(vr,3),size(vr,4)/),vrid)
         ENDIF
         CALL io_WRITE(vrid,(/1,1,1,1/),(/SIZE(vr,1),SIZE(vr,2),         &
     &    SIZE(vr,3),SIZE(vr,4)/),vr)
         CALL io_dclose(vrid,hdferr)
      ENDIF
      CALL io_gclose(gid,hdferr) 
      CALL io_gclose(ggid,hdferr) 
      CALL io_hdfclose(fid,hdferr) 
      CALL io_check('gf_iodop$gf_wrtcoul ',hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_lodcoul(file_type,layer,vpw,qpw)                        
                                                                        
      SUBROUTINE  gf_lodcoul(file_type,layer,vpw,qpw) 
!******************************************                             
!     reads the coulomb potential and the pseudo-charge from a file     
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)           :: file_TYPE 
      INTEGER,INTENT(IN)           :: layer 
      COMPLEX,OPTIONAL,INTENT(OUT) :: vpw(:) 
      COMPLEX,OPTIONAL,INTENT(OUT) :: qpw(:,:) 
      !>                                                                
                                                                        
      !<-- Locals                                                       
      INTEGER(HID_T) :: fid,pwid,gid,ggid 
      INTEGER        :: hdferr,DIM(3),n 
      !>                                                                
      CALL io_hdfopen("gf_pot.hdf",H5F_ACC_RDONLY_F,fid,hdferr) 
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr) 
      CALL io_gopen(ggid,"coulpot",gid,hdferr) 
      IF (PRESENT(vpw)) THEN 
         CALL io_dopen(gid, 'coul', pwid, hdferr) 
         CALL io_datadim(pwid,DIM) 
         n = DIM(1) 
         IF (n /= SIZE(vpw)) THEN 
            WRITE(*,*) "Warning: Wrong dimension of coul potential" 
            vpw = 0.0 
            n   = MIN(n,SIZE(vpw)) 
         ENDIF 
         CALL io_READ(pwid,(/1,-1/),(/n,1/),vpw(:n)) 
         CALL io_dclose(pwid,hdferr) 
      ENDIF 
      IF (PRESENT(qpw)) THEN 
         CALL io_dopen(gid, 'charge', pwid, hdferr) 
         CALL io_datadim(pwid,DIM) 
         IF (DIM(1) /= SIZE(qpw,1).OR.DIM(2) /= SIZE(qpw,2)) THEN 
            WRITE(*,*) "Warning: Wrong dimension of pseudo-charge" 
            qpw = 0.0 
            DIM(1) = MIN(DIM(1),SIZE(qpw,1)) 
            DIM(2) = MIN(DIM(2),SIZE(qpw,2)) 
         ENDIF 
                                                                        
         CALL io_READ(pwid,(/1,1,-1/),(/DIM(1),DIM(2),1/)               &
     &        ,qpw)                                                     
         CALL io_dclose(pwid,hdferr) 
      ENDIF 
      CALL io_gclose(gid,hdferr) 
      CALL io_gclose(ggid,hdferr) 
      CALL io_hdfclose(fid,hdferr) 
                                                                        
      CALL io_check('gf_iodop%lodcoul ',hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:gf_renamepot(file_type)                                    
                                                                        
      SUBROUTINE gf_renamepot(file_type,mpi_subcom,layer)
!******************************************                             
!     rename the gf_potfile to a new file                               
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER , INTENT(IN) :: file_TYPE,mpi_subcom,layer
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n,irank
      CHARACTER(LEN = 20) :: filename 
                                                                        
      CHARACTER(LEN = 11)  :: file
      CHARACTER(LEN = 200) :: command 
      LOGICAL             :: exists

      !>                                                                
#ifdef CPP_MPI
      include "mpif.h"
      integer:: mpierr
      CALL MPI_COMM_RANK(mpi_subcom,irank,mpierr)
#else
      irank=0
#endif
      if (irank==0.and.layer==1) then
      !find a filename                                                  
      	IF (file_type == gf_cdnfile) THEN
        	 file      ="gf_cdn"
      	ELSE
         	file ="gf_pot"
      	ENDIF
                                                                        
                                                                        
      	DO n=1,999
         	WRITE(filename,"(a6,i0,a4)") trim(file),n,".hdf"
         	INQUIRE(FILE =filename,EXIST= exists)
         	exit
      	ENDDO
      	WRITE(filename,"(a6,i0,a4)") trim(file),priv_lastfileNo(file_type)+1,".hdf"
#if defined(CPP_AIX)&&defined(CPP_MPI)
	    filename=trim(filename)//"\0"
	    file=trim(file)//".hdf\0"
        call rename(file,filename)
#else
      	WRITE(command,"(a,a,a,a)") "mv ",trim(file),".hdf ",filename
      	WRITE(*,*) "execute:",command
      	CALL system(command)
#endif
      endif


#ifdef CPP_MPI
     call MPI_BARRIER(mpi_subcom,mpierr)
#endif
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- F: priv_lastfileNo(file_type)                                
      FUNCTION priv_lastfileNo(file_type) 
!-----------------------------------------------                        
!  searches last file-no of an old file in the current directory        
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: file_TYPE 
      INTEGER                :: priv_lastfileNO 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      CHARACTER(LEN = 6)  :: file 
      CHARACTER(LEN = 20) :: filename 
      LOGICAL             :: exists 
      !>                                                                
      IF (file_type == gf_cdnfile) THEN 
         file      ="gf_cdn" 
      ELSE 
         file ="gf_pot" 
      ENDIF 
      DO n = 999,1,-1 
         WRITE(filename,"(a6,i0,a4)") file,n,".hdf" 
         INQUIRE(FILE =filename,EXIST= exists) 
         IF (exists) THEN 
            priv_lastfileNO=n 
            RETURN 
         ENDIF 
      ENDDO 
      priv_lastfileno = 0 
      END FUNCTION 
      !>


      SUBROUTINE gf_loddop_name(file_name,layer,jspins,atoms,stars,sphhar,vr,vpw)
!***********************************************************************
!     This subroutine reads in the potential or the charge density
!     needed by the green fct part of FLEUR.
!
!
!                                  Daniel Wortmann, Juelich, 2001
!***********************************************************************
      USE m_gf_types
      USE hdf5
      USE m_hdf_tools
      IMPLICIT NONE
      character(len=*),INTENT(IN)    :: file_name
      INTEGER,       INTENT(IN)      :: jspins,layer
      TYPE(t_atoms), INTENT(IN)      ::atoms
      TYPE(t_stars), INTENT(IN)      ::stars
      TYPE(t_sphhar),INTENT(IN)      ::sphhar
      REAL,optional,          INTENT(INOUT):: vr(:,0:,:,:)
      COMPLEX,optional,       INTENT(INOUT):: vpw(:,:)


! locals
!     scalars
!     the hdf-IDs
      INTEGER(HID_T)::fid,pwid,pwodid,vrid,gid,ggid
      INTEGER::hdferr,dims(3)
!     needed for input checking
      INTEGER::ntype_in,n3d_in,jrid_in,nlhd_in,jspins_in


      integer :: n


!     temporary array for interstitial pot
      REAL, ALLOCATABLE:: pwd_in(:,:,:)
      REAL,ALLOCATABLE :: vr_tmp(:,:,:,:)
!
!  lmax might have been changed!
!
      vr(:,0:,:,:)=0.0
!
!
!  OPEN FILE and get ID's
      call io_hdfopen(file_name,H5F_ACC_RDONLY_F,fid,hdferr)
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr)
      CALL io_gopen (ggid,"pot",gid,hdferr)
      CALL io_dopen(gid, 'vr', vrid, hdferr)
      CALL io_dopen(gid, 'pw', pwid, hdferr)
      CALL io_check('gf_iodop$gf_loddop:open ',hdferr)

!     check no of spins in file
      call io_datadim(pwid,dims)
      jspins_in=dims(2)
      if (present(vr)) then
!
!
!      read MT-part from file
!
      CALL io_read_att(vrid,"ntype",ntype_in)
      IF (ntype_in/=atoms%ntype) THEN
         WRITE (*,*) 'Wrong no of atoms'
         CALL juDFT_error('Wrong MT-Potential')
      ENDIF
      CALL io_read_att(vrid,"jrid",jrid_in)
      CALL io_read_att(vrid,"nlhd",nlhd_in)
      vr=0.0
      !read spins seperately
      nlhd_in=min(nlhd_in,sphhar%nlhd)
      jrid_in=min(jrid_in,size(vr,1))
      ALLOCATE(vr_tmp(jrid_in,nlhd_in+1,ntype_in,1))
      CALL io_read(vrid,(/1,1,1,1/),                                    &
     &     (/jrid_in,nlhd_in+1,ntype_in,1/)                             &
     &     ,vr_tmp)
      vr(:jrid_in,0:nlhd_in,:ntype_in,1:1)=vr_tmp
      IF (jspins>1) THEN
         if (jspins_in<2) then
              vr(:jrid_in,0:nlhd_in,:ntype_in,2:2)=vr(:jrid_in,0:nlhd_in,:ntype_in,1:1)
         else
            CALL io_read(vrid,(/1,1,1,2/),                                 &
     &      (/jrid_in,nlhd_in+1,ntype_in,1/)                             &
     &      ,vr_tmp)
            vr(:jrid_in,0:nlhd_in,:ntype_in,2:2)=vr_tmp
         endif
      ENDIF
      DEALLOCATE(vr_tmp)
      endif
!
!     read interstitial part
!
      if(present(vpw)) then
      CALL io_read_att(pwid, "n3d", n3d_in)

      vpw=0.0
                                                                        
      ALLOCATE(pwd_in(n3d_in,1,2))
      CALL io_read(pwid,(/1,1,1/),                                      &
     &        (/n3d_in,1,2/),pwd_in)
      vpw(1:min(n3d_in,stars%nq3),1)=                             &
     &     cmplx(pwd_in(1:min(n3d_in,stars%nq3),1,1),pwd_in(1:min(n3d_in&
     &     ,stars%nq3),1,2))
     if (jspins>1) Then
         if (jspins_in<2) then
              vpw(1:min(n3d_in,stars%nq3),2)=vpw(1:min(n3d_in,stars%nq3),1)
         else
            CALL io_read(pwid,(/1,2,1/),(/n3d_in,1,2/),pwd_in)
            vpw(1:min(n3d_in,stars%nq3),2)=cmplx(pwd_in(1:min(n3d_in,stars%nq3),1,1),pwd_in(1:min(n3d_in,stars%nq3),1,2))
         endif
      endif
      DEALLOCATE(pwd_in)
      endif


!  close all id's
      CALL io_dclose(vrid,hdferr)
      CALL io_dclose(pwid,hdferr)
      CALL io_gclose(gid,hdferr)
      CALL io_gclose(ggid,hdferr)
      CALL io_check('gf_iodop$gf_loddop:close ',hdferr)
      CALL io_hdfclose(fid,hdferr)
      END SUBROUTINE
                                                                        
#ifdef OLDSTUFF_NOT_USED_ANYMORE                                        
                                                                        
      !<-- S:gf_removepot(file_type,layer)                              
                                                                        
      SUBROUTINE gf_removepot(file_type,layer) 
!******************************************                             
!     rename the current pot-group in a sc-cycle                        
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     ::  file_TYPE 
      INTEGER,INTENT(IN)     :: layer 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T) :: fid,ggid 
      INTEGER        :: hdferr 
      !>                                                                
      fid = priv_openfile(file_type,mpi) 
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr) 
      CALL io_gdelete(ggid,'pot',hdferr) 
      CALL io_gclose(ggid,hdferr) 
      CALL io_check("remove-pot from:",hdferr) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
                                                                        
      !<-- S:gf_wrtcoul(file_type,layer,vpw,qpw)                        
                                                                        
      SUBROUTINE  gf_wrtcoul(file_TYPE,layer,stars,vpw,qpw,vr)
!******************************************                             
!     writes the coulomb potential to the file                          
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_hdf_tools 
      use m_juDFT 
      USE m_gf_types
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)          :: file_type 
      INTEGER,INTENT(IN)          :: layer 
      TYPE(t_stars),INTENT(IN)    :: stars 
      COMPLEX,OPTIONAL,INTENT(IN) :: vpw(:) 
      COMPLEX,OPTIONAL,INTENT(IN) :: qpw(:,:) 
      REAL,INTENT(IN),OPTIONAL    :: vr(:,:,:,:)
      !>                                                                
                                                                        
      !<-- Locals                                                       
      INTEGER(HID_T) :: fid,pwid,avid,gid,ggid,vrid
      INTEGER        :: hdferr,n,dim(3) 
      LOGICAL        :: old 
      REAL           :: av(2*stars%mx3+3) 
      !>                                                                
      fid = priv_openfile(file_type) 
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr) 
      IF (hdferr /= 0)  CALL juDFT_error("Could not write Coulpot") 
      old = io_groupexists(ggid,"coulpot") 
      IF (old) THEN 
         CALL io_gopen(ggid,"coulpot",gid,hdferr) 
         IF(PRESENT(vpw)) THEN 
            IF (io_dataexists(gid,"coul")) THEN 
               CALL io_dopen(gid,"coul",pwid,hdferr) 
               CALL io_datadim(pwid,dim) 
               CALL io_dclose(pwid,hdferr) 
               IF (DIM(1) / = SIZE(vpw)) THEN 
                  old   = .FALSE. 
                  CALL io_gclose(gid,hdferr) 
                  CALL io_gdelete(ggid,'coulpot',hdferr) 
                  CALL io_gcreate(ggid,"coulpot",gid,hdferr) 
               ENDIF 
            ELSE 
               old = .FALSE. 
            ENDIF 
         ENDIF 
         IF(PRESENT(qpw)) THEN 
            IF (io_dataexists(gid,"charge")) THEN 
               CALL io_dopen(gid,"charge",pwid,hdferr) 
               CALL io_datadim(pwid,dim) 
               CALL io_dclose(pwid,hdferr) 
               IF (DIM(1)/ = SIZE(qpw,1).OR.DIM(2) /=SIZE(qpw,2)) THEN 
                  old = .FALSE. 
                  CALL io_gclose(gid,hdferr) 
                  CALL io_gdelete(ggid,'coulpot',hdferr) 
                  CALL io_gcreate(ggid,"coulpot",gid,hdferr) 
               ENDIF 
            ELSE 
               old=.FALSE. 
            ENDIF 
                                                                        
         ENDIF 
         IF (.NOT.old) THEN 
            CALL io_gclose(gid,hdferr) 
            CALL io_gdelete(ggid,'coulpot',hdferr) 
            CALL io_gcreate(ggid,"coulpot",gid,hdferr) 
         ENDIF 
      ELSE 
         CALL io_gcreate(ggid,"coulpot",gid,hdferr) 
      ENDIF 
      IF (PRESENT(vpw)) THEN 
         IF (old) THEN 
            CALL io_dopen(gid,"coul",pwid,hdferr) 
            CALL io_dopen(gid,"pwav",avid,hdferr) 
         ELSE 
            CALL io_createvar(gid,"coul", H5T_NATIVE_DOUBLE,(/SIZE(vpw) &
     &           ,2/),pwid)                                             
            CALL io_createvar(gid,"pwav",H5T_NATIVE_DOUBLE,(/2*stars%mx3&
     &           +3/),avid)                                             
         ENDIF 
         CALL io_WRITE(pwid,(/1,1/),(/SIZE(vpw),1/),REAL(vpw)) 
         CALL io_WRITE(pwid,(/1,2/),(/SIZE(vpw),1/),AIMAG(vpw)) 
         av = priv_getPlanar(stars,vpw) 
         CALL io_WRITE(avid,(/1/),SHAPE(av),av) 
         CALL io_dclose(avid,hdferr) 
         CALL io_dclose(pwid,hdferr) 
      ENDIF 
      IF (PRESENT(qpw)) THEN 
         IF (old) THEN 
            CALL io_dopen(gid,"charge",pwid,hdferr) 
         ELSE 
            CALL io_createvar(gid,"charge", H5T_NATIVE_DOUBLE,(/SIZE(qpw&
     &           ,1),SIZE(qpw,2),2/),pwid)                              
         ENDIF 
         CALL io_WRITE(pwid,(/1,1,1/),(/SIZE(qpw,1),SIZE(qpw,2),1/)     &
     &        ,REAL(qpw))                                               
         CALL io_WRITE(pwid,(/1,1,2/),(/SIZE(qpw,1),SIZE(qpw,2),1/)     &
     &        ,AIMAG(qpw))                                              
         CALL io_dclose(pwid,hdferr) 
      ENDIF 
      IF (PRESENT(vr)) THEN
         IF (old) THEN
            CALL io_dopen(gid,"vr",vrid,hdferr)
         ELSE
            CALL io_createvar(gid,"vr", H5T_NATIVE_DOUBLE,(/SIZE(vr,1),  &
     &       SIZE(vr,2),size(vr,3),size(vr,4)/),vrid)
         ENDIF
         CALL io_WRITE(vrid,(/1,1,1,1/),(/SIZE(vr,1),SIZE(vr,2),         &
     &    SIZE(vr,3),SIZE(vr,4)/),vr)
         CALL io_dclose(vrid,hdferr)
      ENDIF
      CALL io_gclose(gid,hdferr) 
      CALL io_gclose(ggid,hdferr) 
      CALL io_check('gf_iodop$gf_wrtcoul ',hdferr) 
                                                                        
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_lodcoul(file_type,layer,vpw,qpw)                        
                                                                        
      SUBROUTINE  gf_lodcoul(file_type,layer,vpw,qpw) 
!******************************************                             
!     reads the coulomb potential and the pseudo-charge from a file     
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)           :: file_TYPE 
      INTEGER,INTENT(IN)           :: layer 
      COMPLEX,OPTIONAL,INTENT(OUT) :: vpw(:) 
      COMPLEX,OPTIONAL,INTENT(OUT) :: qpw(:,:) 
      !>                                                                
                                                                        
      !<-- Locals                                                       
      INTEGER(HID_T) :: fid,pwid,gid,ggid 
      INTEGER        :: hdferr,DIM(3),n 
      !>                                                                
      fid=priv_openfile(file_type) 
      CALL io_gopen(fid,io_layername(layer),ggid,hdferr) 
      CALL io_gopen(ggid,"coulpot",gid,hdferr) 
      IF (PRESENT(vpw)) THEN 
         CALL io_dopen(gid, 'coul', pwid, hdferr) 
         CALL io_datadim(pwid,DIM) 
         n = DIM(1) 
         IF (n /= SIZE(vpw)) THEN 
            WRITE(*,*) "Warning: Wrong dimension of coul potential" 
            vpw = 0.0 
            n   = MIN(n,SIZE(vpw)) 
         ENDIF 
         CALL io_READ(pwid,(/1,-1/),(/n,1/),vpw(:n)) 
         CALL io_dclose(pwid,hdferr) 
      ENDIF 
      IF (PRESENT(qpw)) THEN 
         CALL io_dopen(gid, 'charge', pwid, hdferr) 
         CALL io_datadim(pwid,DIM) 
         IF (DIM(1) /= SIZE(qpw,1).OR.DIM(2) /= SIZE(qpw,2)) THEN 
            WRITE(*,*) "Warning: Wrong dimension of pseudo-charge" 
            qpw = 0.0 
            DIM(1) = MIN(DIM(1),SIZE(qpw,1)) 
            DIM(2) = MIN(DIM(2),SIZE(qpw,2)) 
         ENDIF 
                                                                        
         CALL io_READ(pwid,(/1,1,-1/),(/DIM(1),DIM(2),1/)               &
     &        ,qpw)                                                     
         CALL io_dclose(pwid,hdferr) 
      ENDIF 
      CALL io_gclose(gid,hdferr) 
      CALL io_gclose(ggid,hdferr) 
      CALL io_check('gf_iodop%lodcoul ',hdferr) 
      CALL priv_closefile(file_type) 
      END SUBROUTINE 
                                                                        
      !>                                                                
#endif                                                                  
      END                                           
