!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur_interface 
#include "juDFT_env.h"

!----------------------------------------------------------------       
!Module provides interface to FLEUR                                     
!USE m_grp_k                                                            
!USE m_ylm                                                              
!USE radfun                                                             
!----------------------------------------------------------------       
      USE m_ylm,ONLY:fleur_ylm=>ylm4 
      USE m_dsphbs,ONLY:fleur_dsphbs=>dsphbs 
      USE m_sphbes,ONLY:fleur_sphbes=>sphbes 
      USE m_dwigner,ONLY:fleur_d_wigner=>d_wigner 
      USE m_intgr,ONLY:fleur_intgr3=>intgr3 
      USE m_starf,ONLY:fleur_starf3=>starf3 
      USE m_gaunt,ONLY:fleur_gaunt => gaunt1 
      USE m_ifft,ONLY:fleur_ifft235 => ifft235 
      USE m_rwnoco,ONLY:fleur_rw_noco => rw_noco 
      USE m_spgrot 
      IMPLICIT NONE
      CONTAINS 
      !<--S:fleur_grp_k                                                 
      SUBROUTINE fleur_grp_k(mrot,mrot_k,amat,bk,nclass,nirr,char_table,&
     &     grpname,irrname,su)                                          
      USE m_grp_k 
      IMPLICIT NONE 
                                                                        
      INTEGER, INTENT(IN)                :: mrot(:,:,:) 
      REAL, INTENT(IN)                   :: bk(3),amat(3,3) 
      COMPLEX , INTENT(OUT)              :: char_table(:,:) 
      INTEGER, INTENT(OUT)               :: mrot_k(:,:,:),nclass,nirr 
      CHARACTER(LEN = 7), INTENT(OUT)    :: grpname
      CHARACTER(LEN = 5), INTENT(OUT)    :: irrname(:)
      COMPLEX, OPTIONAL, INTENT(OUT)     :: su(:,:,:) 
      CALL grp_k(mrot,mrot_k,amat,bk,nclass,nirr,char_table,            &
     &    grpname,irrname,su)                                           
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:fleur_radfun(l,el,vr,ntyp,atoms,f,g,us,dus,uds,duds,ddn,nod
      SUBROUTINE fleur_radfun(l,el,vr,ntyp,atoms,                       &
     &     f,g,us,dus,uds,duds,ddn,nodeu,noded,wronk)                   
!-----------------------------------------------                        
!   Interface to radfun                                                 
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_radfun 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER, INTENT (IN)      :: l,ntyp 
      REAL,    INTENT (IN)      :: el,vr(:) 
      TYPE(t_atoms),INTENT(IN)  :: atoms 
      INTEGER, INTENT (OUT)    :: noded,nodeu 
      REAL,    INTENT (OUT):: ddn,duds,dus,uds,us,wronk 
      REAL,    INTENT (OUT):: f(:,:),g(:,:) 
      !>                                                                
                                                                        
      CALL radfun(                                                      &
     &           l,el,vr,atoms%jri(ntyp),atoms%rmsh(1                   &
     &           ,ntyp) ,atoms%dx(ntyp),maxval(atoms%jri),f,g,us        &
     &           ,dus,uds,duds,ddn,nodeu,noded,wronk)                   
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:fleur_spgrot(                                              
      SUBROUTINE fleur_spgrot(                                          &
     &                  sym,                                            &
     &                  k,                                              &
     &                  kr,phas)                                        
!-----------------------------------------------                        
!  interface to spgrot from FLEUR                                       
!           (last modified: 05-03-22) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_constants 
      USE m_spgrot 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_sym),INTENT(IN) :: sym 
      INTEGER   ,INTENT(IN)  :: k(3) 
      INTEGER   ,INTENT(OUT) :: kr(:,:) 
      REAL   ,INTENT(OUT)    :: phas(:) 
      !>                                                                
                                                                        
      CALL spgrot(sym%nop,sym%symor,2.*pimach(),sym%mrot,sym%tau        &
     &     ,sym%invtab,k,kr,phas)                                       
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:fleur_mapatom                                              
                                                                        
      SUBROUTINE fleur_mapatom(sym,atoms,                               &
     &     cell,soc,                                                    &
     &     invsatnr,multab,invarop,invarind)                            
!-----------------------------------------------                        
!  interface to mapatom form FLEUR                                      
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_mapatom 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_sym),INTENT(INOUT)   :: sym 
      TYPE(t_cell),INTENT(INOUT)  :: cell 
      TYPE(t_atoms),INTENT(INOUT) :: atoms 
      TYPE(t_soc),INTENT(IN)      :: soc 
      INTEGER, INTENT (OUT)       :: invarind(:) 
      INTEGER, INTENT (OUT)       :: invarop(:,:) 
      INTEGER, INTENT (OUT)       :: invsatnr(:) 
      INTEGER, INTENT (OUT)       :: multab(:,:) 
      !>                                                                
                                                                        
      CALL mapatom(sym%nop,atoms%ntype,atoms%nat,atoms%ntype            &
     &     ,atoms%ntypsy,cell%amat,cell%bmat,sym%mrot,atoms%neq,sym%tau &
     &     ,atoms%taual,.FALSE.,atoms%n_u,soc%soc,soc%theta,soc%phi,    &
     &     sym%invs,atoms%ngopr                                         &
     &     ,atoms%invsat,invsatnr,cell%bbmat,multab,sym%invtab,invarop  &
     &     ,invarind)                                                   
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:fleur_prp_qfft_map                                         
      SUBROUTINE fleur_prp_qfft_map(stars,lapw,sym,igq2_fft,igq_fft) 
!-----------------------------------------------                        
!  Interface to prp_qfft_map                                            
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_prpqfftmap 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(IN) :: stars 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_sym),INTENT(IN)   :: sym 
      INTEGER,INTENT(OUT)      :: igq2_fft(0:),igq_fft(0:) 
      !>                                                                
                                                                        
      CALL prp_qfft_map(                                                &
     &     stars%nq3,lapw%kq1_fft,lapw%kq2_fft,lapw%kq3_fft,sym%nop,    &
     &     stars%kv3,sym%nop2,sym%mrot,.TRUE.,                          &
     &     lapw%kq1_fft,lapw%kq2_fft,lapw%kq3_fft,                      &
     &     lapw%nq3_fft,lapw%kmxq2_fft,lapw%kmxq_fft,                   &
     &     igq2_fft,igq_fft)                                            
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:fleur_prp_qfft(stars,noco,lapw,cell)                       
      SUBROUTINE fleur_prp_qfft(stars,noco,lapw,cell) 
!-----------------------------------------------                        
!    interface to fleur prp_qfft                                        
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_prpqfft 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(IN)   :: stars 
      TYPE(t_noco),INTENT(IN)    :: noco 
      TYPE(t_lapw),INTENT(INOUT) :: lapw 
      TYPE(t_cell),INTENT(IN)    :: cell 
      !>                                                                
      CALL prp_qfft(                                                    &
     &     stars%nq2,stars%nq3,lapw%kq1_fft,lapw%kq2_fft,lapw%kq3_fft,  &
     &     stars%nq2,stars%nq3,stars%nstr2,stars%nstr,stars%gmax        &
     &     ,cell%bmat,stars%sk2,stars%sk3,stars%kv3,                    &
     &     noco%l_ss,                                                   &
     &     lapw%rkmax,                                                  &
     &     lapw%kq1_fft,lapw%kq2_fft,lapw%kq3_fft,                      &
     &     lapw%nq2_fft,lapw%nq3_fft,lapw%kmxq2_fft,lapw%kmxq_fft)      
                      !no l_ss                                          
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: fleur_mksphhar(atoms,cell,sym,sphhar)                     

      SUBROUTINE fleur_mksphhar(atoms,cell,sym,sphhar,gfprep) 
!-----------------------------------------------                        
! Generates the spherical harmonics                                     
! interface to localsym from FLEUR                                      
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_localsym 
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_atoms),INTENT(INOUT) :: atoms 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_sym),INTENT(IN)     ::sym 
      TYPE(t_sphhar),INTENT(OUT) ::sphhar 
      LOGICAL,INTENT(IN),OPTIONAL :: gfprep 
                                                                        
                                                                        
      INTEGER ::ntypsy(atoms%nat) 
      INTEGER::idum(1,1,1) 
      COMPLEX::cdum(1,1,1) 
      ntypsy = atoms%ntypsy 
      !Call for determining the dimensions                              
      sphhar%ntypsd = 0 
      CALL local_sym(                                                   &
     &     maxval(atoms%lmax),atoms%lmax,sym%nop,sym%mrot,sym%tau,      &
     &     atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat          &
     &     ,atoms%taual,sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.TRUE.    &
     &     ,atoms%nlhtyp,atoms%ntypsy,idum(:,1,1),idum(:,:,1),idum(:,:,1&
     &     ),idum(:,:,:),cdum)                                          
                                                                     !no
                                                                        
      IF (PRESENT(gfprep)) THEN 
                                 !in prep-mode simply assign the symetri
         atoms%ntypsy=ntypsy 
      ELSE 
         IF (ANY(ntypsy /= atoms%ntypsy)) CPP_error('Symmetry of atom changed!')
             
      ENDIF 
      !ALLOCATE                                                         
      WRITE(6,*) 'Symmetry types:',sphhar%ntypsd 
      WRITE(6,*) 'Max no of lm:',sphhar%nlhd 
      WRITE(6,*) 'Max no of members of spherical harmonics:',sphhar%memd 
                                                                        
      ALLOCATE(sphhar%nlh(sphhar%ntypsd)) 
      ALLOCATE(sphhar%llh(0:sphhar%nlhd,sphhar%ntypsd)) 
      ALLOCATE(sphhar%nmem(0:sphhar%nlhd,sphhar%ntypsd)) 
      ALLOCATE(sphhar%mlh(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd)) 
      ALLOCATE(sphhar%clnu(sphhar%memd,0:sphhar%nlhd,sphhar%ntypsd)) 
      !Construct spherical harmonics                                    
      CALL local_sym(                                                   &
     &     maxval(atoms%lmax),atoms%lmax,sym%nop,sym%mrot,sym%tau,      &
     &     atoms%nat,atoms%ntype,atoms%neq,cell%amat,cell%bmat          &
     &     ,atoms%taual,sphhar%nlhd,sphhar%memd,sphhar%ntypsd,.FALSE.   &
     &     ,atoms%nlhtyp,atoms%ntypsy,sphhar%nlh,sphhar%llh,sphhar%nmem &
     &     ,sphhar%mlh,sphhar%clnu)                                     
                                                               !allocate
                                                                        
                                                                        
      END SUBROUTINE 

      !>                                                                
                                                                        
      !<-- S:fleur_prp_xcfft(stars,lapw,cell,xcpot)                     
                                                                        
      SUBROUTINE fleur_prp_xcfft(stars,lapw,cell,xcpot) 
!-----------------------------------------------                        
!  interface to prp_xcfft                                               
!  also includes a call to prp_xcfft_box for dimensioning               
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_prpxcfft 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(IN)    :: stars 
      TYPE(T_lapw),INTENT(IN)     :: lapw 
      TYPE(t_cell),INTENT(IN)     :: cell 
      TYPE(t_xcpot),INTENT(INOUT) :: xcpot 
      !>                                                                
                                                                        
      CALL prp_xcfft_box(                                               &
     &     xcpot%gmaxxc,cell%bmat,                                      &
     &     xcpot%kxc1d,xcpot%kxc2d,xcpot%kxc3d)                         
                                                                        
      CALL prp_xcfft(                                                   &
     &     stars%nq3,xcpot%kxc1d,xcpot%kxc2d,xcpot%kxc3d,               &
     &     stars%nq3,stars%nstr,xcpot%gmaxxc,lapw%rkmax_inp,cell%bmat   &
     &     ,stars%sk3,stars%kv3,xcpot%gmaxxc,xcpot%kxc1_fft             &
     &     ,xcpot%kxc2_fft,xcpot%kxc3_fft,xcpot%nxc3_fft,xcpot%kmxxc_fft&
     &     )                                                            
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: fleur_renpara(atoms,jspins,enpara)                        
      SUBROUTINE fleur_renpara(layers,atoms,jspins,enpara) 
!*****************************************************************      
!     DESC:loads enpara file                                            
!     Daniel Wortmann, Wed Aug 28 10:40:06 2002                         
!*****************************************************************      
      USE m_gf_types 
      USE m_enpara,ONLY:r_enpara 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_layers),INTENT(IN)  ::layers 
      INTEGER,INTENT(IN)         ::jspins 
      TYPE(t_atoms),INTENT(IN)   ::atoms(:) 
      TYPE(t_enpara),INTENT(INOUT) ::enpara(:) 
                                                                        
      !>                                                                
                                                                        
      LOGICAL,PARAMETER::film=.TRUE. 
      INTEGER,PARAMETER::nw=1 
      INTEGER::jsp,n,l,layer 
                                                                        
      REAL   ,ALLOCATABLE :: el(:,:),ello(:,:) 
      INTEGER,ALLOCATABLE :: skiplo(:) 
      LOGICAL,ALLOCATABLE :: lchange(:,:),llochg(:,:) 
                                                                        
      OPEN(40,FILE="enpara",STATUS="old") 
      DO layer=1,layers%num_layers 
        !<--Allocate memory                                             
         ALLOCATE(lchange( 0:MAXVAL(atoms(layer)%lmax0),atoms(layer     &
     &        )%ntype))                                                 
        ALLOCATE( enpara(layer)%lchange( 0:MAXVAL(atoms(layer)%lmax0),  &
     &                                   atoms(layer)%ntype,jspins)  )  
        ALLOCATE(llochg( atoms(layer)%nlod,atoms(layer)%ntype)) 
        ALLOCATE( enpara(layer)%llochg ( atoms(layer)%nlod,             &
     &                                   atoms(layer)%ntype,jspins)  )  
        ALLOCATE(el(0:MAXVAL(atoms(layer)%lmax0),atoms(layer)%ntype)) 
        ALLOCATE( enpara(layer)%el( 0:maxval(atoms(layer)%lmax0),       &
     &                                   atoms(layer)%ntype,jspins)  )  
        ALLOCATE(ello(atoms(layer)%nlod, atoms(layer)%ntype)) 
        ALLOCATE( enpara(layer)%ello( atoms(layer)%nlod,                &
     &                                   atoms(layer)%ntype,jspins)  )  
        ALLOCATE(skiplo(atoms(layer)%ntype)) 
        ALLOCATE( enpara(layer)%skiplo( atoms(layer)%ntype,jspins)   ) 
        !>                                                              
        el=0.0;ello=0.0
        !read the enpara(layer)-file                                    
        DO jsp=1,jspins 
         CALL r_enpara(                                                 &
     &        MAXVAL(atoms(layer)%lmax0),atoms(layer)%nlod,             &
     &        atoms(layer)%ntype,                                       &
     &        film,jsp,nw,atoms(layer)%nlo,                             &
     &        atoms(layer)%lmax0,atoms(layer)%neq,                      &
     &        atoms(layer)%zatom,atoms(layer)%llo                       &
     &        ,skiplo,                                                  &
     &         ello,el,enpara(layer)%evac(:,jsp),                       &
     &         lchange                                                  &
     &        ,llochg,enpara(layer)%lchg_v(1,jsp)                       &
     &        ,enpara(layer)%enmix(jsp))                                
         enpara(layer)%el(:,:,jsp)=el 
         enpara(layer)%ello(:,:,jsp)=ello 
         enpara(layer)%skiplo(:,jsp)=skiplo 
         enpara(layer)%lchange(:,:,jsp)=lchange 
         enpara(layer)%llochg(:,:,jsp)=llochg 
         !<-- check if no LO has the same energy parameter as           

         !corresponding LAPW                                            
         DO n = 1,atoms(layer)%ntype 
            DO l = 1,atoms(layer)%nlo(n) 
               IF (ABS(enpara(layer)%el(atoms(layer)%llo(l,n),n,jsp)    &
     &             -enpara(layer)%ello(l,n,jsp))<1E-4) THEN             
                  WRITE(6,*) "layer:",layer," Atom:",n," Spin:",jsp 
                  WRITE(6,*) MAXVAL(atoms(layer)%lmax0) 
                  WRITE(6,*) "LO:",l 
                  WRITE(6,*) "L(LO):",atoms(layer)%llo(l,n) 
                  WRITE(6,*) "Enpara(LO):",enpara(layer)%ello(l,n,jsp) 
                  WRITE(6,*) "Enpara(APW):",enpara(layer)%el(atoms(layer&
     &                 )%llo(l,n),n,jsp)                                
                  CPP_error("Same energy parameter for LO and (L)APW")
               ENDIF 
            ENDDO 
         ENDDO 

         !>                                                             
        ENDDO 
         DEALLOCATE(el,ello,skiplo,lchange,llochg) 
                                                                        
            !layer                                                      
      ENDDO 
      CLOSE(40) 
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: fleur_boxdim(bmat,arltv1,arltv2,arltv3)                   
      SUBROUTINE fleur_boxdim(bmat,arltv1,arltv2,arltv3) 
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL,    INTENT (IN)  :: bmat(3,3) 
      REAL,    INTENT (OUT) :: arltv1,arltv2,arltv3 
      !>                                                                
      EXTERNAL boxdim 
      CALL boxdim(bmat,arltv1,arltv2,arltv3) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: fleur_setlomap(n,atoms)                                   
      SUBROUTINE fleur_setlomap(n,atoms) 
!-----------------------------------------------                        
!interface to setlomap from FLEUR                                       
!           (last modified: 05-03-22) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)           :: n 
      TYPE(t_atoms),INTENT(INOUT)  :: atoms 
      !>                                                                
      EXTERNAL setlomap 
       CALL setlomap(                                                   &
     &           atoms%ntype,atoms%nlod,atoms%llod,                     &
     &           n,atoms%nlo,atoms%llo,                                 &
     &           atoms%lo1l,atoms%nlol,atoms%l_dulo,atoms%ulo_der)      
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: fleur_cored(spins,jspin,atoms,,,,)                        
      SUBROUTINE fleur_cored(jspins,jspin,atoms,vr,                     &
     &     q_int,rh,seig,rho)                                           
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)       :: jspins,jspin 
      TYPE(t_atoms),INTENT(IN) :: atoms 
      REAL   ,INTENT(IN)       :: vr(:,0:,1:,:) 
      REAL,INTENT(OUT)         :: q_INT(:,:) 
      REAL,INTENT(OUT)         :: rh(:,:) 
      REAL,INTENT(OUT)         :: seig 
      REAL,INTENT(INOUT)       :: rho(:,0:,:,:) 
      !>                                                                
      EXTERNAL cored 
      CALL cored(jspins,jspin,atoms%ntype,atoms%neq,atoms%zatom         &
     &     ,atoms%rmsh,atoms%dx,atoms%jri,.FALSE.,.FALSE.,atoms%ncst    &
     &     ,atoms%rmt,rho(:,0:,:,:),SIZE(vr,1),SIZE(vr,4)               &
     &     ,MAXVAL(atoms%jri)+100,SIZE(rho,2)-1                         &
     &     ,MAX(MAXVAL(atoms%ncst),30),atoms%ntype,vr(1:,0,1:,jspin),0  &
     &     ,q_int,rh,seig)                                              
                                                                      !G
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: fleur_cdnval(,,,,)                                        
#ifdef CPP_NEVER
      SUBROUTINE fleur_cdnval(weights,invsatnr,multab,invarop           &
     &     ,invarind,igq_fft,cell,atoms,sphhar,sym,kpts,stars,enpara    &
     &     ,mpi,neigd,msh,nstd,nspd,lapw,jspins,jspin,lmd,llpd,vr,rho   &
     &     ,qpw, odi,ods,sk2,phi2)                                      
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_od_types 
      USE m_gf_types 
      USE m_cdnval 
      USE m_constants 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_cell),INTENT(IN)     :: cell 
      TYPE(t_atoms),INTENT(IN)    :: atoms 
      TYPE(t_sphhar),INTENT(IN)   :: sphhar 
      TYPE(t_sym),INTENT(INOUT)   :: sym 
      TYPE(t_kpts),INTENT(IN)     :: kpts 
      TYPE(t_stars),INTENT(IN)    :: stars 
      TYPE(t_enpara),INTENT(INOUT) :: enpara 
      TYPE(t_mpi),INTENT(IN)      :: mpi 
      TYPE(t_lapw),INTENT(IN)     :: lapw 
      INTEGER,INTENT(IN)          :: neigd,msh,nstd,nspd,jspins,jspin 
      INTEGER,INTENT(IN)          :: llpd,lmd 
      !The potential                                                    
      REAL, INTENT(IN)             :: vr(:,0:,:,:) 
      !The charge                                                       
      REAL,INTENT(INOUT)           :: rho(:,:,:,:) 
      COMPLEX,INTENT(INOUT)        :: qpw(:,:) 
      !kpts-weights                                                     
      REAL,INTENT(IN)              :: weights(:,:,:) 
      INTEGER,INTENT(IN)           :: invarind(:) 
      INTEGER,INTENT(IN)           :: invsatnr(:) 
      !for fft-box                                                      
      INTEGER,INTENT(IN) :: igq_fft(:) 
                                                                        
      INTEGER, INTENT (IN)       :: multab(:,:) 
      INTEGER, INTENT (IN)       :: invarop(:,:) 
                                                                        
      !odim                                                             
      REAL, INTENT(IN),OPTIONAL           :: sk2(:),phi2(:) 
      TYPE (od_inp), INTENT (IN),OPTIONAL :: odi 
      TYPE (od_sym), INTENT (IN),OPTIONAL :: ods 
      !>                                                                
      !<-- Locals                                                       
                                                                        
      !dummy variables                                                  
      REAL,ALLOCATABLE    :: rdum1(:),rdum3(:,:,:),rdum4(:,:,:,:) 
      REAL,ALLOCATABLE    :: rdum5(:,:,:,:,:) 
      REAL                :: timedum(9) 
      INTEGER             :: idum0 
      INTEGER,ALLOCATABLE :: idum2(:,:) 
      COMPLEX,ALLOCATABLE :: cdum2(:,:),cdum3(:,:,:,:),cdum4(:,:,:,:) 
                                                                        
      !output of cdnval                                                 
      REAL      :: chmom(atoms%ntype,jspins),clmom(3,atoms%ntype,jspins) 
      REAL      :: qis(neigd,kpts%nkpts,jspins) 
      COMPLEX   :: cdom(stars%nq3) 
      COMPLEX   :: qa21(atoms%ntype) 
      !LDA+U                                                            
      COMPLEX   :: d_wgn(-3:3,-3:3,3,sym%nop) 
                                      !in fact the last dim             
      COMPLEX   :: n_mmp(-3:3,-3:3,1) 
                                                                        
      INTEGER             :: i 
                                             !should be n_u             
      timedum = 0 
      d_wgn = 0.0 
      !>                                                                
                                                                        
      IF (atoms%n_u /= 0)  CPP_error("LDA+U does not work in gf_cdnval")
                                                                        
      CALL cdnval(0,0,kpts%nkpts,jspin,1,.FALSE.,.FALSE.,               &
     &     .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,0,cell%vol           &
     &     ,cell%volint,atoms%volmts,reshape((/(0.0,i=1,jspins)/),(/    &
     &      jspins,1/)),stars%sk3,.FALSE.,.FALSE.,0,idum2,0,0.0,0.0,0,  &
     &     .FALSE.,1,neigd,kpts%nkpts,1,atoms%ntype,maxval(atoms%ntypsy)&
     &     ,maxval(sphhar%nlh),stars%nq3,maxval(atoms%jri)              &
     &     ,maxval(atoms%lmax),jspins,maxval(sphhar%nmem),lapw%nvd,1    &
     &     ,stars%mx1,stars%mx2,stars%mx3,lapw%kq1_fft,lapw%kq2_fft     &
     &     ,lapw%kq3_fft,mpi%irank,mpi%isize,lapw%nbasfcn,sym%nop       &
     &     ,stars%nq2,atoms%nat,1,1,nspd,lmd,llpd,lapw%kq1_fft          &
     &     ,lapw%kq2_fft,lapw%kq3_fft,lapw%nq3_fft,lapw%kmxq_fft,igq_fft&
     &     ,stars%igfft,stars%pgfft,(/0.0/),weights,vr(:,0,:,:),vr,(/(/ &
     &    0.0/),(/0.0/)/),0,0.0,0.0,0.0,0.0,enpara%lchange,enpara%lchg_v&
     &     ,.FALSE.,(/(.FALSE.,i = 1,atoms%ntype)/),cell%bbmat,2        &
     &     * SQRT(pimach()),invarop,invarind,multab,sym%invtab          &
     &     ,atoms%invsat,invsatnr,atoms%ntype,.FALSE.,sym%zrfs,sym%invs &
     &     ,sym%symor,1,1,1,.FALSE.,jspins,0.01,cell%area,cell%omtil    &
     &     ,cell%z1,cell%amat,cell%bmat,atoms%taual,atoms%pos,sym%nop2  &
     &     ,sphhar%llh,sphhar%nmem,sphhar%mlh,sphhar%clnu               &
     &     ,MAXVAL(atoms%ntypsy),atoms%ngopr,sphhar%nlh,atoms%ntypsy    &
     &     ,sym%mrot,sym%tau,stars%nstr,stars%nstr2,stars%ig,stars%ig2  &
     &     ,stars%rgphs,stars%mx1,stars%mx2,stars%mx3,stars%kv2         &
     &     ,stars%kv3,stars%nq3,stars%nq2,atoms%rmt,atoms%dx,atoms%rmsh &
     &     ,atoms%jri,atoms%lmax,atoms%neq,atoms%nlod,atoms%llod        &
     &     ,atoms%nlo,atoms%llo,enpara%llochg,atoms%lda_u(:)%l,atoms%n_u&
     &     ,enpara%skiplo,atoms%nlol,atoms%lo1l,atoms%l_dulo            &
     &     ,atoms%ulo_der,atoms%lapw_l,0.0,1,1,.FALSE.,(/0.0,0.0/),(/   &
     &     .0,0.0/),msh,nstd,                                           &
     &     odi,ods,sk2,phi2,                                            &
     &     atoms%ncst,atoms%zatom,d_wgn,rdum1,rdum1,(/0.,0.,0./),idum0  &
     &     ,n_mmp,timedum,enpara%ello,enpara%el,enpara%evac,rdum3,qpw   &
     &     ,cdum4,rho,rdum3,rdum4,rdum5,cdom,cdum2,cdum3,qa21,chmom     &
     &     ,clmom,qis)                                                  
                                                                        
      END SUBROUTINE 

      !>                                                                
#endif
      !<-- S: fleur_convndim(gmaxr,ncvd)                                
      SUBROUTINE fleur_convndim(                                        &
     &                     gmaxr,                                       &
     &                     ncvd)                                        
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_convndim 
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)     ::  gmaxr 
      INTEGER,INTENT(OUT)    ::  ncvd 
      !>                                                                
                                                                        
      CALL convn_DIM(                                                   &
     &                     gmaxr,                                       &
     &                     ncvd)                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: fleur_convn(ncvd,ntype,gmax,rmt,ncv)                      
      SUBROUTINE fleur_convn(ncvd,ntype,gmax,rmt,ncv) 
!-----------------------------------------------                        
!                                                                       
!           (last modified: 06-04-27) D. Wortmann                       
!-----------------------------------------------                        
      USE m_convn 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER, INTENT (IN)  :: ntype,ncvd 
      REAL,    INTENT (IN)  :: gmax 
      INTEGER, INTENT (OUT) :: ncv(ntype) 
      REAL,    INTENT (IN)  :: rmt(ntype) 
      !>                                                                
      CALL  convn(                                                      &
     &     ncvd,ntype,gmax,rmt,                                         &
     &     ncv)                                                         
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:fleur_calcenpara(jspins,enpara,atoms,vr)                   
      SUBROUTINE fleur_calcenpara(jspins,mpi,enpara,atoms,vr) 
!-----------------------------------------------                        
!  call lodpot from FLEUR to determine energy parameters                
!           (last modified:08-12-11) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_lodpot 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)             :: jspins 
      TYPE(t_enpara),INTENT(INOUT)   :: enpara 
      TYPE(t_mpi),INTENT(IN)         :: mpi 
      TYPE(t_atoms),INTENT(IN)       :: atoms 
      REAL   ,INTENT(IN)             :: vr(:,:,:,:) 
      !>                                                                
      !<-- Locals                                                       
      !temps for vacuum                                                 
      REAL :: vz(1,1,1) 
      REAL :: evac0(1,1,1) 
      REAL :: evac(1,1,1) 
      !bounds                                                           
      REAL :: bound_lo(1),bound_up(1) 
      REAL :: elup(1),ellow(1) 
      !enparas                                                          
      REAL :: ello(size(enpara%ello,1),size(enpara%ello,2)              &
     &     ,SIZE(enpara%ello,3))                                        
      REAL :: el0(size(enpara%el,1),size(enpara%el,2),size(enpara%el,3) &
     &     ,1)                                                          
      REAL :: el(size(enpara%el,1),size(enpara%el,2),size(enpara%el,3)  &
     &     ,1)                                                          
      !>                                                                
      el0(:,:,:,1) = enpara%el 
      el = 0.0
      ello = 0.0
      elup=0.0;ellow=0.0
      CALL lodpot(                                                      &
     &     mpi%irank,                                                   &
     &     jspins,MAXVAL(atoms%lmax0),MAXVAL(atoms%jri),size(vr,2)-1     &
     &     ,atoms%ntype,1,0,jspins,atoms%ntype,1,.FALSE.,1,0,0,atoms%jri&
     &     ,atoms%dx,atoms%rmt,atoms%rmsh,atoms%lmax,vr,vz,atoms%llo    &
     &     ,atoms%zatom,el0,evac0,enpara%ello,atoms%nlo,MAXVAL(atoms%nlo&
     &     ),atoms%l_dulo,ellow,elup,el,evac,ello,bound_lo,bound_up)    
                                                                        
                                                                        
      enpara%ello = ello 
      enpara%el   = el(:,:,:,1) 
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
