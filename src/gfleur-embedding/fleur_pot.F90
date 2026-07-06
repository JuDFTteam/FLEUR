!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur_pot 
      use m_juDFT
          IMPLICIT NONE
!-------------------------------------------------------------          
!Module provides an interface to the potential subroutines of FLEUR     
!as used in the GF-code                                                 
!------------------------------------------------------------           
      CONTAINS 
      !<-- S: fleur_xcpot(jspins,atoms,stars,sphhar,vr,vpw)             
                                                                        
      SUBROUTINE fleur_xcpot(layer,jspins,atoms,stars,sphhar,xcpot,     &
     &      sym,cell,vr                                                 &
     &     ,vpw)                                                        
!-----------------------------------------------                        
!Subroutine that encapsulates the XC-potential generation of FLEUR      
!           (last modified: 05-03-23) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_iodop 
      USE m_visxc 
      USE m_visxcg 
      USE m_vmtxcg 
      USE m_vmtxc
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)        :: layer 
      INTEGER,INTENT(IN)        :: jspins 
      TYPE(t_atoms),INTENT(IN)  :: atoms 
      TYPE(t_stars),INTENT(IN)  :: stars 
      TYPE(t_sphhar),INTENT(IN) :: sphhar 
      TYPE(t_xcpot),INTENT(INOUT) :: xcpot 
      TYPE(t_sym),INTENT(IN)    :: sym 
      TYPE(t_cell),INTENT(IN)   :: cell 
      REAL   ,INTENT(INOUT)     :: vr(:,0:,:,:) 
      COMPLEX,INTENT(INOUT)     :: vpw(:,:) 
      !>                                                                
      !<-- Locals                                                       
                                                                        
      REAL,    ALLOCATABLE      :: rho(:,:,:,:) 
      COMPLEX, ALLOCATABLE      :: qpw(:,:) 
      REAL,    ALLOCATABLE      :: excr(:,:,:) 
      COMPLEX,ALLOCATABLE       :: vpw_w(:,:),cdom(:),excpw(:) 
      COMPLEX,ALLOCATABLE       :: vxpw_w(:,:),vxpw(:,:) 
      REAL,ALLOCATABLE          :: vxr(:,:,:,:)
      REAL,POINTER              :: ufft(:) 
                                                                        
                                                                        
      INTEGER           :: ifftd,ifftxc3d,n 
      INTEGER           :: ichsmrg = 0 
                                        !rhmn: rho-min.                 
      REAL              :: rhmn = 1.E+10 
      !dimensions                                                       
      INTEGER:: ntypsd,nlhd,jmtd,lmaxd,memd 
                     !!calculate from lmaxd                             
      INTEGER:: nspd 
                                                                        
      !>                                                                
                                                                        
      !<--Generate the step-function!                                   
                                                                        
                                                                        
      ALLOCATE(ufft(0:27*stars%mx1*stars%mx2*stars%mx3-1)) 
      ufft = 1.0 
                                                                        
      !>                                                                
                                                                        
      !<--Allocate                                                      
      ntypsd=size(sphhar%clnu,3) 
      nlhd=size(sphhar%clnu,2)-1 
      jmtd=size(atoms%rmsh,1) 
      lmaxd=maxval(atoms%lmax) 
      memd=size(sphhar%clnu,1) 
      nspd=(lmaxd+1+mod(lmaxd+1,2))*(2*lmaxd+1) 
      ALLOCATE(excr(jmtd,0:nlhd,atoms%ntype)) 
      ALLOCATE (rho(SIZE(atoms%rmsh,1),0:SIZE(sphhar%clnu,2)-1          &
     &     ,atoms%ntype,jspins))               
      ALLOCATE(vxr(SIZE(atoms%rmsh,1),0:SIZE(sphhar%clnu,2)-1          &
     &     ,atoms%ntype,jspins))

      ALLOCATE (qpw(stars%nq3,jspins)) 
      ALLOCATE( excpw(stars%nq3),vpw_w(stars%nq3,jspins),cdom(stars%nq3)&
     &     )                                                            
      cdom=0.0
      ALLOCATE(vxpw(stars%nq3,jspins),vxpw_w( stars%nq3,jspins))                                         
      ifftxc3d = xcpot%kxc1d*xcpot%kxc2d*xcpot%kxc3d 
      ifftd=27*stars%mx1*stars%mx2*stars%mx3 
      !>                                                                
      !<-- load charge density                                          
      CALL gf_loddop(GF_CDNFILE,layer,jspins,                           &
     &        atoms,stars,sphhar,                                       &
     &        rho,qpw)
      !>                                                                
      !<-- Coulomb potential is equal for both spins                    
      vr(:,:,:,jspins) = vr(:,:,:,1) 
      vpw(:,jspins) = vpw(:,1) 
      !>                                                                
      !<-- interstitial region                                          
                !not used in this code!                                 
      vpw_w=0.0 
      IF ((xcpot%igrd==0).AND.(xcpot%icorr/=-1)) THEN 
         CALL visxc(                                                    &
     &        ifftd,stars%mx1,stars%mx2,stars%mx3,stars%nq3,jspins,     &
     &        qpw,cdom,.FALSE.,stars%kimax,stars%igfft,stars%pgfft,ufft,&
     &        xcpot%icorr,.FALSE.,xcpot%krla,jspins,stars%nq3,stars%nstr&
     &        ,vpw,vpw_w,vxpw,vxpw_w,excpw)                                         
      ELSEIF ((xcpot%igrd>0).OR.(xcpot%icorr==-1)) THEN 
        IF (xcpot%lwb) THEN 
          WRITE(6,'('' W+B trick cancelled out. visxcwb uses at'',      &
     &              '' present common block cpgft3. visxcwb needs'',/,  &
     &              '' to be reprogrammed according to visxcg.f'')')    
           CALL juDFT_error("visxcwb",calledby="fleur_pot.F90")
       ELSE 
! ff                                                                    
! ff for GGA: calculate gradient of whole matrix, then project on       
! magnetization-axis in vmatgen.F                                       
          ! NOTE: noco was removed here!!                               
          CALL visxcg(                                                  &
     &         ifftd,stars%mx1,stars%mx2,stars%mx3,stars%nq3,jspins     &
     &         ,sym%nop,xcpot%kxc1d,xcpot%kxc2d,xcpot%kxc3d,ifftxc3d    &
     &         ,stars%kv3,sym%nop2,sym%mrot,cell%bmat                   &
     &         ,xcpot%kxc1_fft,xcpot%kxc2_fft,xcpot%kxc3_fft            &
     &         ,xcpot%nxc3_fft,xcpot%kmxxc_fft,qpw,cdom,stars%kimax     &
     &         ,stars%igfft,stars%pgfft,ufft,xcpot%icorr,.FALSE.,jspins &
     &         ,stars%nq3,stars%nstr,xcpot%igrd,xcpot%ndvgrd            &
     &         ,.FALSE.,.FALSE.,(/0.0,0.0,0.0/)&
     &         ,xcpot%chng,xcpot%lwb,rhmn    &
     &         ,ichsmrg,vpw,vpw_w,vxpw,vxpw_w,excpw)                                
       ENDIF 
      ELSE 
         CALL juDFT_error("something wrong with igrd before visxc",calledby="fleur_pot.F90")
      ENDIF 
                                                                        
      !>                                                                
      !<-- muffin tin spheres region                                    
      IF ((xcpot%igrd==0).AND.(xcpot%icorr/=-1)) THEN 
         CALL vmtxc(                                                    &
     &              jspins,memd,nlhd,ntypsd,jmtd,atoms%ntype,lmaxd,nspd,&
     &        sphhar%clnu,sphhar%mlh,sphhar%nmem,sphhar%llh,sphhar%nlh  &
     &        ,atoms%rmsh,atoms%ntypsy,atoms%jri,atoms%nat,atoms%neq,rho&
     &        ,xcpot%icorr,.FALSE.,xcpot%krla,atoms%ntype,jspins,ntypsd &
     &        ,vr,excr)                                                 
      ELSEIF ((xcpot%igrd>0).OR.(xcpot%icorr==-1)) THEN 
         CALL vmtxcg(                                                   &
     &        jspins,memd,nlhd,ntypsd,jmtd,atoms%ntype,lmaxd,nspd,      &
     &        sphhar%clnu,sphhar%mlh,sphhar%nmem,sphhar%llh,sphhar%nlh  &
     &        ,atoms%rmsh,atoms%ntypsy,atoms%jri,atoms%dx,atoms%nat     &
     &        ,atoms%neq,rho,xcpot%icorr,.FALSE.,xcpot%krla,atoms%ntype &
     &        ,jspins,ntypsd,xcpot%igrd,xcpot%ndvgrd  &
     &        ,xcpot%chng,vxr,vr    &
     &        ,rhmn,ichsmrg,excr)
      ELSE 
         CALL juDFT_error("something wrong with igrd before vmtxc",calledby="fleur_pot.F90")
      ENDIF 
      !>                                                                
      DEALLOCATE(ufft) 
      DEALLOCATE(excr,rho,qpw,excpw,vpw_w,cdom) 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S: fleur_vmts(jspins,stars,atoms,vpw,rho,vr)                 
      SUBROUTINE fleur_vmts(mpi_comm,jspins,stars,atoms,sphhar,sym,cell &
     &     ,vpw,rho,vr)                                                 
!-----------------------------------------------                        
!   interface to vmts of FLEUR                                          
!           (last modified: 05-03-23) D. Wortmann                       
!-----------------------------------------------                        
      USE m_vmts 
      USE m_gf_types 
      USE m_od_types, ONLY : od_inp, od_sym 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)        :: jspins,mpi_comm 
      TYPE(t_stars),INTENT(IN)  :: stars 
      TYPE(t_atoms),INTENT(IN)  :: atoms 
      TYPE(t_sphhar),INTENT(IN) :: sphhar 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_sym),INTENT(IN)    :: sym 
      COMPLEX, INTENT (IN)      :: vpw(:,:) 
      REAL,    INTENT (IN)      :: rho(:,0:,:,:) 
      REAL,    INTENT (OUT)     :: vr(:,0:,:,:) 
      !>                                                                
      !<-- Locals                                                       
      !dimensions                                                       
      INTEGER:: ntypsd,nlhd,jmtd,lmaxd,memd 
      INTEGER :: isize,irank 
      TYPE (od_inp) :: odi 
      TYPE (od_sym) :: ods 
#ifdef CPP_MPI                                                          
      INCLUDE 'mpif.h' 
      INTEGER :: ierr(3) 
      CALL MPI_COMM_RANK (MPI_COMM,irank,ierr) 
      CALL MPI_COMM_SIZE (MPI_COMM,isize,ierr) 
#else                                                                   
      irank = 0;isize = 1 
#endif                                                                  
      odi%d1 = .FALSE. 
      ntypsd = SIZE(sphhar%clnu,3) 
      nlhd=size(sphhar%clnu,2)-1 
      jmtd=size(atoms%rmsh,1) 
      lmaxd=maxval(atoms%lmax) 
      memd=size(sphhar%clnu,1) 
      !>                                                                
      CALL vmts(irank,isize,mpi_comm,                                   &
     &      stars%nq3,jspins,memd,nlhd,ntypsd,jmtd,atoms%ntype,lmaxd,   &
     &      atoms%ntype,stars%nq3,atoms%zatom,atoms%dx,atoms%rmsh       &
     &      ,atoms%jri,atoms%rmt,stars%sk3,atoms%lmax,sphhar%clnu       &
     &      ,sphhar%mlh,sphhar%nmem,sphhar%llh,stars%nstr,atoms%ntypsy  &
     &      ,sphhar%nlh,sym%invtab,sym%nop,atoms%nat,atoms%neq,stars%kv3&
     &     ,sym%mrot,cell%bmat,sym%tau,atoms%taual,sym%symor,vpw,rho,odi&
     &     ,ods,vr)                                                     
      END SUBROUTINE 
      !>                                                                
      END                                           
