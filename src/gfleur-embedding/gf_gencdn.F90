!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_gencdn 
      USE m_constants, ONLY: oUnit
      use m_juDFT
      use m_juDFT 
      IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_gencdn(layer,jspin,l_chargefromstates,              &
     &     jspins,gfinp,input,atoms,cell,sym,kpts,stars,sphhar,mpi,     &
     &     enpara,vr,qpw,rho,qmtl,neigd,l_noco,qtot_el,qtot_nuc)
!*****************************************************************      
! DESC:This subroutine constructes the charge-density from the Green    
! function. It can be considered as the GF-version of cdngen            
!                                                                       
!                          Daniel Wortmann, Sat Feb 23 12:52:14 2002    
!*****************************************************************      
      USE m_gf_cdntot 
      USE m_gf_types 
      USE m_gf_plot 
      USE m_gf_iodop 
      USE m_cored
      IMPLICIT NONE 
!     .. Scalar Arguments ..                                            
      INTEGER, INTENT(IN)          :: jspin,jspins,layer
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_atoms),INTENT(INOUT)  :: atoms 
      TYPE(t_embinp),INTENT(IN)     :: gfinp 
      TYPE(t_sphhar),INTENT(IN)    :: sphhar 
      TYPE(t_sym),INTENT(INOUT)    :: sym 
      TYPE(t_kpts),INTENT(IN)      :: kpts 
      TYPE(t_cell),INTENT(INOUT)   :: cell 
      TYPE(t_stars),INTENT(IN)     :: stars 
      TYPE(t_gfmpi),INTENT(IN)       :: mpi 
      TYPE(t_enpara),INTENT(INOUT) :: enpara 
!Arguments for cdnstates
      LOGICAL,INTENT(IN)           :: l_noco
      LOGICAL,INTENT(IN)           :: l_chargefromstates 
      REAL   ,INTENT(INOUT)        :: qtot_el(:),qtot_nuc(:) 
      INTEGER                      :: neigd 
                                                                        
!     .. Array Arguments ..                                             
      !potential (spherical)                                            
      REAL,    INTENT(IN)       :: vr(:,0:,:,:) 
      ! charge density                                                  
      COMPLEX,INTENT(INOUT)     :: qpw(:,:) 
      REAL,   INTENT(INOUT)     :: rho(:,0:,:,:) 
      ! l-like charge                                                   
      REAL,INTENT(INOUT)        :: qmtl(0:,:) 
!     .. Locals ..                                                      
      !dims                                                             
                                                                        
      INTEGER,PARAMETER         :: nspd = 300 
      INTEGER                   :: msh 
      !loops                                                            
      INTEGER                   :: nk,l,itype 
      LOGICAL                   :: lexist 
      REAL,ALLOCATABLE          :: q_INT(:,:),rh(:,:) 
      REAL                      :: seig 
      real, allocatable         :: vr_fixed(:,:,:,:)
                                                                        
                                                                        
                                                                        
      !determine dims                                                   
      msh = MAXVAL(atoms%jri)+100 
      ALLOCATE(rh(msh,atoms%ntype),q_INT(atoms%ntype,jspins)) 
                                                                        
                                                                        
      !<-- Collect charge from different PEs                            
                                                                        
#ifdef CPP_MPI                                                          
      IF (mpi%isize/=1) THEN 
         CALL gf_MPIcollect(rho(:,0:,:,jspin),qpw(:,jspin),qmtl,mpi     &
     &        ,kpts%nkpts)                                              
      ENDIF 
#endif                                                                  
                                      ! this PE does not have the       
      IF (mpi%k_kpts(1) /= 1) RETURN 
                                      ! correct layer-charge            
                                                                        
      !>                                                                
                                                                        
      ! Spin doubling...                                                
      IF (jspins==1) THEN 
         rho(:,0:,:,1)=2*rho(:,0:,:,1) 
         qpw(:,1)=2* qpw(:,1) 
         qmtl(:,:)=2*qmtl(:,:) 
      ENDIF 
      rho(:,0:,:,jspin)=-1*rho(:,0:,:,jspin) 
      if (l_noco) rho(:,0:,:,jspins)=-1*rho(:,0:,:,jspins)
                                                                        
!add the density from the eig-file                                      
      IF (l_chargefromstates) THEN 
          CALL juDFT_error("chargefromstates not implemented",calledby="gf_gencdn.F90")
      ENDIF 
                                                                        
!Plot l-like charges!                                                   
      WRITE (oUnit,FMT=8000) layer 
 8000 FORMAT (/,i3,5x,'l-like charge',/,t6,'atom',t15,'s',t24,'p',      &
     &     t33,'d',t42,'f',t51,'total')                                 
      DO itype=1,atoms%ntype 
         WRITE (oUnit,FMT = 8100) layer,itype, (qmtl(l,itype),l=0,3),       &
     &        sum(qmtl(:,itype))                                        
      ENDDO 
      WRITE (oUnit,FMT=8200) layer,sum(qmtl(:,:)) 
 8200 FORMAT (/,i3,5x,'total valence charge in MT-spheres:',f12.5) 
 8100 FORMAT (i3,' -->',i2,2x,4f9.5,2x,f9.5) 
      CALL priv_plotPlanar(stars,qpw) 
      !call gf_plot in case a charge plot is needed                     
      CALL gf_plot(layer,stars,cell,atoms,sym,jspins,                   &
     &     qpw(:,:),GF_PLOT_CHARGE,sphhar,rho(:,0:,:,:))                
                                                                        
      INQUIRE(FILE ="gf_cdn.diff.hdf",EXIST = lexist) 
      IF (.NOT.gfinp%l_totalmix.OR.lexist) THEN 
                                                                        
!                                                                       
!     ---> add in core density                                          
!

         if (jspin==1) CALL gf_wrtdop(GF_cdnfile,layer,jspins,                 &
     &           gfinp,atoms,stars,sphhar,                              &
     &           rho,qpw,.FALSE.,mpi%self_subcom,.true.)  !write only valence charge to file


         inquire(file="gf_pot_fixedcore.hdf",exist=lexist)
         if (lexist) then
            write(*,*) "Fixed core potential calculation"
            allocate(vr_fixed(size(vr,1),size(vr,2),size(vr,3),size(vr,4)))
            vr_fixed=0.0
            call gf_loddop_name("gf_pot_fixedcore.hdf",layer,jspins,atoms,stars,sphhar,vr_fixed)
            call priv_cored(input,sphhar,jspin,atoms,vr_fixed,q_int,seig,rho(:,0:,:,:))
            if (l_noco) CALL priv_cored(input,sphhar,2,atoms,vr_fixed,q_int,seig,rho(:,0:,:,:))
            deallocate(vr_fixed)
         else
               CALL priv_cored(input,sphhar,jspin,atoms,vr,q_int,     &
     &              seig,rho(:,0:,:,:))
              if (l_noco) CALL priv_cored(input,sphhar,2,atoms,vr,q_int,seig,rho(:,0:,:,:))
         endif
      ENDIF 
!                                                                       
!     print information on total charges                                
!                                                                       
      CALL priv_plotPlanar(stars,qpw) 
      IF (jspin==jspins.or.l_noco) THEN
         CALL gf_cdntot(layer,mpi%fmpi,jspins,stars,cell,atoms,rho(:,0:,:,:) &
     &        ,qpw(:,:),qtot_el(layer),qtot_nuc(layer))                 
         !<-- generate gf_cdn.diff file for totalcharge mixing          
         IF (gfinp%l_totalmix.AND..NOT.lexist) THEN 
            qpw(:,:) = qpw(:,:)/qtot_el(layer) 
            rho = rho/qtot_el(layer) 
            CALL gf_wrtdop(GF_cdndifffile,layer,jspins,                 &
     &           gfinp,atoms,stars,sphhar,                              &
     &           rho,qpw,.FALSE.,mpi%self_subcom)                       
                                                  !no noco, mpi does not
         ENDIF 
      ENDIF 
      END SUBROUTINE 
      !<-- S:priv_plotPlanar(stars,qpw)                                 
      SUBROUTINE priv_plotPlanar(stars,qpw) 
!-----------------------------------------------                        
!     calculate the planar average of the charge/potential in pw        
!             (last modified: 04-08-13) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_fft_singleton 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(IN) :: stars 
      COMPLEX,INTENT(IN)       :: qpw(:,:) 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      COMPLEX :: vz(0:2*stars%mx3-1,size(qpw,2)) 
      INTEGER :: n 
      !>                                                                
      vz = 0.0 
      DO n =-stars%mx3+1,stars%mx3-1 
         IF (stars%ig(0,0,n) == 0) CYCLE 
         IF (n<0) THEN 
            vz(n+2*stars%mx3,:) = qpw(stars%ig(0,0,n),:) 
         ELSE 
            vz(n,:)             = qpw(stars%ig(0,0,n),:) 
         ENDIF 
      ENDDO 
      DO n = 1,SIZE(qpw,2) 
         vz(:,n) = fft(vz(:,n),inv = .TRUE.) 
      ENDDO 
      WRITE(oUnit,*) "Planar Charge" 
      WRITE(oUnit,*) "z      qpw        " 
      DO n = -stars%mx3+1,stars%mx3-1 
         IF (n<0) THEN 
            WRITE(oUnit,"(i5,1x,3(f0.7,1x))") n,REAL(vz(n+2*stars%mx3,:)) 
         ELSE 
            WRITE(oUnit,"(i5,1x,3(f0.7,1x))") n,REAL(vz(n,:)) 
         ENDIF 
      ENDDO 
                                                                        
      END SUBROUTINE 
                                                                        
      SUBROUTINE priv_cored(input,sphhar,jspin,atoms,vr,q_int,seig,rho)
      !per-type loop over the modern cored (replaces the old EXTERNAL
      !cored call; deliberately NOT cdncore, which would additionally
      !spread core tails into the interstitial)
      USE m_cored
      USE m_gf_types
      IMPLICIT NONE
      TYPE(t_input),INTENT(IN)  :: input
      TYPE(t_sphhar),INTENT(IN) :: sphhar
      INTEGER,INTENT(IN)        :: jspin
      TYPE(t_atoms),INTENT(IN)  :: atoms
      REAL,INTENT(IN)           :: vr(:,0:,1:,:)
      REAL,INTENT(OUT)          :: q_int(:,:)
      REAL,INTENT(OUT)          :: seig
      REAL,INTENT(INOUT)        :: rho(:,0:,:,:)

      REAL :: rhc(atoms%msh,atoms%ntype,input%jspins)
      REAL :: tec(atoms%ntype,input%jspins)
      REAL :: qint2(atoms%ntype,input%jspins)
      REAL :: seig_t
      INTEGER :: iType

      seig=0.0
      qint2=0.0
      rhc=0.0; tec=0.0
      DO iType=1,atoms%ntype
         CALL cored(input,jspin,iType,atoms,rho,sphhar,.FALSE.,          &
     &              vr(:,0,:,jspin),qint2,rhc,tec,seig_t)
         seig=seig+seig_t
      ENDDO
      q_int(:,jspin)=qint2(:,jspin)
      END SUBROUTINE priv_cored

      END                                           
