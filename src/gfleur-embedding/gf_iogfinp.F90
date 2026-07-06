!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_iogfinp 
#include "juDFT_env.h"
      IMPLICIT NONE
      CONTAINS 
      !<--S: GF_RGFINP                                                  
                                                                        
      SUBROUTINE GF_RGFINP(                                             &
     &     layers,gfinp,atoms,xcpot,mix,lapw,tsoc,fermi,stars,tnoco)    
!***********************************************************************
!     This subroutine:                                                  
!     - reads the gf_inp file                                           
!     - set some variables to defaults                                  
!     modified to make gf_inp-file more flexible                        
!     Daniel Wortmann, Juelich 2003                                     
!                                                                       
!***********************************************************************
      USE m_gf_types 
      USE m_fleur_INTERFACE,ONLY:fleur_setlomap,fleur_convndim          &
     &     ,fleur_convn                                                 
      IMPLICIT NONE 
      !<-- Scalar Arguments ..                                          
                                                                        
      TYPE(t_layers),INTENT(IN)    :: layers 
      TYPE(t_gfinp), INTENT(INOUT)  :: gfinp 
      TYPE(t_atoms),INTENT(INOUT) :: atoms(:) 
      TYPE(t_xcpot),INTENT(INOUT)   :: xcpot(:) 
      TYPE(t_mix),INTENT(INOUT)     :: mix 
      TYPE(t_lapw),INTENT(INOUT)    :: lapw(:) 
      TYPE(t_soc),INTENT(INOUT)     :: tsoc 
      TYPE(t_fermi),INTENT(INOUT)   :: fermi 
      TYPE(t_stars),INTENT(INOUT)   :: stars(:) 
      TYPE(t_noco), INTENT(INOUT)   :: tnoco(:) 
                                                                        
      !>                                                                
      !<-- locals!                                                      
                                                                        
      !mode switches                                                    
      LOGICAL :: eigen,tmat,gmat,charge,dos,totalmix,writeT,potmix      &
     &     ,inp2plot,savemem,writehs,IEC,gproj,band,spectral,fullgreen  &
     &     ,nohelpregion,overlap,intdos,hdfio
      INTEGER :: curr,nz_dos 
      REAL    :: z_min,z_max 
      CHARACTER(LEN=4)::kpts 
      LOGICAL :: doslayer(1000)
      NAMELIST /MODE/ intdos,totalmix,eigen,tmat,gmat,curr,charge,dos   &
     &     ,writeT,potmix,inp2plot,kpts,savemem,writehs,IEC,gproj,band  &
     &     ,spectral,fullgreen,nz_dos,z_min,z_max,nohelpregion,overlap,hdfio,doslayer
      !CBS                                                              
      REAL   ::eps_current,eps_non_bloch,CBS_bz,d1,d2,kappa_max 
      NAMELIST /CBSMODE/ CBS,eps_current,eps_non_bloch,CBS_bz,d1,d2     &
     &     ,kappa_max                                                   
      LOGICAL:: cbs 
      !BASIS                                                            
      LOGICAL::solwil 
      INTEGER::npw 
      INTEGER :: napw(1000) 
      REAL   ::rkmax 
      NAMELIST /BASIS/solwil,npw,napw,rkmax 
      !EMB                                                              
      LOGICAL::addemb,advanced,surface,embmt,embspectral,simplevacuum
                                                                   !bias
      REAL  ::charge_limit,imag_broad,EField,vacuum_energy
      NAMELIST /EMB/addemb,advanced,surface,charge_limit,embmt,         &
     &       imag_broad,Efield,embspectral,simplevacuum,vacuum_energy
                                                               !,bias   
      !SC                                                               
      REAL::trans,c_sc 
      INTEGER::imix,maxiter,iter,precond
      REAL::alpha,spinf,beta
      LOGICAL::gauss,tria 
      REAL::tkb,delgau 
      REAL::gmax_pot,gmax_decouple,g0max,k_kerker,g0scale
      NAMELIST /SC/iter,trans,c_sc,imix,maxiter,alpha,spinf,gauss,tria  &
     &     ,tkb,delgau,gmax_pot,gmax_decouple,k_kerker,g0max,precond,g0scale
      !SOC                                                              
      LOGICAL::soc 
      REAL::theta,phi 
      NAMELIST /SPINORBIT/soc,theta,phi 
! noco                                                                  
      LOGICAL :: noco,ss,mperp,constr 
      REAL    :: qss(3),mix_b 
      NAMELIST /NONCOLLINEAR/noco,ss,mperp,constr,qss,mix_b,alpha,beta
! noco-per atom                                                         
      LOGICAL :: relax
      REAL    :: b_con(2) 
      NAMELIST /NOCOL/alpha,beta,b_con,relax 
      !for XC-pot                                                       
      CHARACTER*4 ::namex 
      LOGICAL     ::all 
      !for ldaU                                                         
      INTEGER l,layer 
      REAL    u,j 
      NAMELIST /ldaU/ l,u,j 
      !for atom-loop                                                    
      INTEGER::n,i 
      ! for LO's                                                        
      INTEGER::llo_tmp(20,maxval(atoms(:)%ntype)) 
      ! for ncv                                                         
      INTEGER::ncvd 
                                                                        
      !for free input                                                   
      CHARACTER(LEN = 10) :: teststring 
      LOGICAL             :: free_input 
      LOGICAL             :: required(3) 
                                                                        
      !>                                                                
      !<--set defaults                                                  
      CPP_juDFT_timestart_info("Reading the gf_inp")
      !CBS                                                              
      gfinp%l_CBS = .FALSE.;gfinp%eps_current = 1E-4 
      gfinp%eps_non_bloch = 1E-4;gfinp%CBS_bz =-1 
      gfinp%dp1     = 0.0;gfinp%dp2 = 0.0 
      gfinp%kappa_max=1E20 
      !EMB                                                              
      gfinp%l_embmt = .FALSE. 
      gfinp%charge_limit = 0.0 
      gfinp%l_addemb = .FALSE. 
      gfinp%l_adv = .FALSE. 
      !gfinp%bias = 0.0                                                 
      gfinp%l_surface  = .FALSE. 
      gfinp%l_embspectral = .FALSE.
                             !                                          
      gfinp%imag_broad = 0.0
      gfinp%l_simplevacuum=.FALSE.                  !
                           !                                            

      gfinp%Efield   = 0.0 
      !SC                                                               
      gfinp%trans = 0.0;gfinp%c_sc = 0.0 
      mix%imix = 0;mix%maxiter = 99;mix%iter =0; 
      mix%alpha = 0.01;mix%spinf = 1.0
      mix%precond=0;mix%k_kerker=0.0;mix%g0max=0.0
      fermi%gauss = .FALSE.;fermi%tria = .FALSE. 
      fermi%tkb = 0.005;fermi%delgau = 0.0 
      stars%gmax_pot=50.0 
      stars%gmax_decouple=10.0 
      !SPINORBIT                                                        
      tsoc%soc = .FALSE.;tsoc%theta = 0.0;tsoc%phi=0.0 
      !Noco                                                             
      tnoco%l_noco = .FALSE. 
      tnoco%l_ss   = .FALSE. 
      required = .FALSE. 
                                                                        
      !>                                                                
                                                                        
      !<--READ the gf_inp-file                                          
                                                                        
      OPEN(92,FILE='gf_inp',FORM='formatted',STATUS='OLD') 
      free_input = .TRUE. 
      DO WHILE(free_input) 
         !<--read a teststring                                          
         READ(92,'(a10)') teststring 
                                          !remove heading blanks        
         teststring = ADJUSTL(teststring) 
         !>                                                             
                                            !skip comment lines         
         IF (teststring(1:1) =='!') CYCLE 
         IF (teststring(1:7) =='*****XC') THEN 
                                          !****XC marks the end of free 
            free_input       = .FALSE. 
            CYCLE 
         ENDIF 
                       !The current line has to be read again           
         BACKSPACE(92) 
         !<--Now test for the different possible imput-lines            
                                                                        
         IF (teststring(1:4)=='gmax') THEN 
            !READ gmax                                                  
            READ(92,'(5x,f10.5,8x,f10.5)') stars(1)%gmax_inp,           &
     &                                       xcpot(1)%gmaxxc            
            IF (xcpot(1)%gmaxxc == 0)                                   &
     &        xcpot(1:layers%num_layers)%gmaxxc=stars(1)%gmax_inp       
            stars(1:layers%num_layers)%gmax_inp=stars(1)%gmax_inp 
            required(1) = .TRUE. 
         ELSEIF (teststring(1:4)=='&MOD') THEN 
            !READ the mode-switches                                     
            spectral = .TRUE.;fullgreen=.FALSE. 
            inp2plot = .FALSE.;band = .FALSE. 
            eigen = .FALSE.;tmat=.TRUE.;curr=0; 
            charge = .FALSE.;dos=.FALSE.;IEC=.FALSE.;intdos=.FALSE. 
            writeT = .FALSE.;gmat=.TRUE.;potmix=.FALSE.;kpts="none" 
            writeHS = .FALSE.;saveMEM=.FALSE.;gproj=.FALSE.             &
     &           ;nohelpregion = .TRUE.;overlap=.FALSE.                 
            nz_dos = 1;z_min=0.0;z_max=1.0;totalmix= .FALSE.;hdfio=.false.
            doslayer=.true.
            READ(92,NML = MODE) 
            gfinp%l_hdfio=hdfio
            gfinp%l_intdos=intdos 
            gfinp%l_totalmix=totalmix 
!            gfinp%l_overlap = overlap                                  
            gfinp%l_spectral=spectral;gfinp%l_fullgreen=fullgreen 
            gfinp%l_eigen = eigen;gfinp%l_tmat=tmat;gfinp%curr=curr 
            gfinp%l_charge = charge;gfinp%l_dos=dos;gfinp%l_band=band 
            gfinp%l_writeT =writeT;gfinp%l_gmat=gmat; 
            gfinp%l_potmix = potmix;gfinp%l_inp2plot=inp2plot 
            gfinp%l_addselfen=.FALSE. 
            required(2)    = .TRUE.;gfinp%kpts =                        &
     &           kpts;gfinp%l_nohelpregion = nohelpregion               
            gfinp%l_writeHS = writeHS;gfinp%l_savemem=savemem 
            gfinp%l_IEC=IEC;gfinp%l_gproj=gproj 
            gfinp%nz_dos=nz_dos;gfinp%z_min=z_min;gfinp%z_max=z_max 
            ALLOCATE(gfinp%l_doslayer(layers%num_layers))
            gfinp%l_doslayer=doslayer(1:layers%num_layers)
         ELSEIF (teststring(1:4) =='&CBS') THEN 
            !READ the CBS                                               
            kappa_max = 1E20; 
            CBS = .FALSE.;eps_current = 1E-4;eps_non_bloch=1E-4 
            CBS_bz =-1.0;d1 = 0.0;d2 = 0.0 
            READ(92,NML = CBSMODE) 
            gfinp%kappa_max=kappa_max 
            gfinp%l_CBS = CBS;gfinp%eps_current=eps_current 
            gfinp%eps_non_bloch = eps_non_bloch;gfinp%CBS_bz=CBS_bz 
            gfinp%dp1 = d1;gfinp%dp2=d2 
         ELSEIF (teststring(1:4) =='&BAS') THEN 
            !READ the BASIS                                             
            rkmax = 3.5;solwil=.FALSE.;npw=0;napw=15 
            READ(92,NML = BASIS) 
            gfinp%l_solwil = solwil;gfinp%npw =                         &
     &           npw;gfinp%napw(:layers%num_layers) =                   &
     &           napw(:layers%num_layers)                               
            lapw(1:layers%num_layers)%rkmax_inp = rkmax 
                                                                        
            required(3) = .TRUE. 
         ELSEIF (teststring(1:4) =='&EMB') THEN 
            !READ EMB                                                   
            embmt=.FALSE. 
            addemb = .FALSE. 
            simplevacuum=.false.
            embspectral=.FALSE.                  !;bias = 0.0
            advanced = .FALSE. 
            surface  = .FALSE. 
            imag_broad = 0.0
            charge_limit=0.0 
            Efield = 0.0 
            vacuum_energy=1.0
            READ(92,NML = EMB) 
            gfinp%vacuum_energy=vacuum_energy
            gfinp%l_simplevacuum=simplevacuum
            gfinp%l_embspectral=embspectral
            gfinp%charge_limit = charge_limit 
            gfinp%l_addemb = addemb 
            gfinp%l_adv    = advanced 
            gfinp%l_embmt  = embmt 
            !gfinp%bias     = bias;                                     
            gfinp%l_surface = surface 
            gfinp%imag_broad = imag_broad
            gfinp%Efield   = Efield 
         ELSEIF (teststring(1:3) =='&SC') THEN 
            !READ SC                                                    
            gmax_pot = 50.0 
            gmax_decouple = 10.0 
            trans = 0.0;c_sc=0.0;imix=1;maxiter=99;alpha=0.05;spinf=0.2 
            gauss = .FALSE.;tria=.FALSE.;precond=0;k_kerker=0
            tkb = 0.005;delgau=0.0;g0scale=1.0;g0max=0.0
            READ(92,NML = SC) 
            gfinp%trans = trans;gfinp%c_sc=c_sc 
            mix%imix    = imix;mix%maxiter = maxiter;mix%iter=iter 
            mix%alpha   = alpha;mix%spinf = spinf 
            mix%precond=precond;mix%k_kerker=k_kerker;mix%g0max=g0max;mix%g0scale=g0scale
            fermi%gauss = gauss;fermi%tria=tria 
            fermi%tkb   = tkb;fermi%delgau = delgau 
            stars(1:layers%num_layers)%gmax_pot=gmax_pot 
            stars(1:layers%num_layers)%gmax_decouple=gmax_decouple 
         ELSEIF (teststring(1:4) =='&SPI') THEN 
            !READ SPINORBIT                                             
            soc = .FALSE.;theta=0.0;phi=0.0 
            READ(92,NML = SPINORBIT) 
            tsoc%soc = soc;tsoc%theta=theta;tsoc%phi=phi 
         ELSEIF (teststring(1:4) =='&NON') THEN 
            !READ NOCO                                                  
            noco = .FALSE.;relax = .FALSE. 
            ss   = .FALSE.; constr = .FALSE.;mperp = .FALSE. 
            qss = 0.0;mix_b=0.0;alpha=0.0;beta=0.0
            READ(92,NML = NONCOLLINEAR) 
            tnoco%l_noco = noco 
            tnoco%l_ss = ss 
            tnoco%l_mperp = mperp 
            tnoco%l_constr = constr 
            tnoco%qss(1) = qss(1) 
            tnoco%qss(2) = qss(2) 
            tnoco%qss(3) = qss(3) 
            tnoco%mix_b  = mix_b 
            tnoco%alpha_int	= alpha
            tnoco%beta_int = beta
         ELSE 
            WRITE(*,*) "Unknown entry in gf_inp-file" 
            WRITE(*,*) "Line starting with:",teststring 
            CPP_error("ERROR in gf_inp")
         ENDIF 
                                                                        
         !>                                                             
      ENDDO 
      !End of free-format input                                         
      IF (ANY(.NOT.required)) CPP_error("Required input missing in gf_inp")
      !<--read fixed part of gf_inp                                     
                                                                        
      !READ XCPOT                                                       
      READ (92,'(a4,5x,l1)') namex,all 
      IF (namex =='l91 ') xcpot(1)%icorr =-1 
      IF (namex =='wign') xcpot(1)%icorr = 1 
      IF (namex=='mjw ') xcpot(1)%icorr = 2 
      IF (namex=='hl  ') xcpot(1)%icorr = 3 
      IF (namex=='bh  ') xcpot(1)%icorr = 3 
      IF (namex=='vwn ') xcpot(1)%icorr = 4 
      IF (namex=='pz  ') xcpot(1)%icorr = 5 
      IF (namex=='pw91') xcpot(1)%icorr = 6 
!     pbe: easy_pbe (phys.rev.lett.77,3865(1996).                       
!     rpbe: rev_pbe (phys.rev.lett.80,890 (1998).                       
!     Rpbe: Rev_pbe (phys.rev.b.59 7413 (1999).                         
      IF (namex=='pbe ') xcpot(1)%icorr = 7 
      IF (namex=='rpbe') xcpot(1)%icorr = 8 
      IF (namex=='Rpbe'.OR.namex=='    ') xcpot(1)%icorr = 9 
      !defaults for gga                                                 
      xcpot(1)%igrd = 0;xcpot(1)%lwb= .FALSE.;xcpot(1)%ndvgrd= 6 
      xcpot(1)%idsprs=0 
      xcpot(1)%chng=-.100D-11;xcpot(1)%iggachk=0 
      xcpot(1)%idsprs0=0 
      xcpot(1)%idsprsl=0;xcpot(1)%idsprsi=0;xcpot(1)%idsprsv=0 
                      !always non-relativistic!!                        
      xcpot(1)%krla=0 
      IF ((namex=='pw91').OR.(namex=='l91 ').OR.                    &
     &    (namex=='pbe ').OR.(namex=='rpbe').OR.                    &
     &    (namex=='Rpbe')) THEN                                       
                       !read more parameters                            
         IF (all) THEN 
            READ (92,FMT = 7121)                                        &
     &           xcpot(1)%igrd,xcpot(1)%lwb,xcpot(1)%ndvgrd,            &
     &           xcpot(1)%idsprs                                        &
     &           ,xcpot(1)%chng                                         
 7121       FORMAT (5x,i1,5x,l1,8x,i1,8x,i1,6x,d10.3) 
            READ (92,FMT = 7122)                                        &
     &           xcpot(1)%iggachk,xcpot(1)%idsprs0,xcpot(1)%idsprsl,    &
     &           xcpot(1)%idsprsi                                       &
     &           ,xcpot(1)%idsprsv                                      
 7122       FORMAT (8x,i1,9x,i1,9x,i1,9x,i1,9x,i1) 
         ELSE 
            xcpot(1)%igrd = 1 
         ENDIF 
      ENDIF 
      xcpot(1:layers%num_layers)=xcpot(1) 
      llo_tmp=-1 
      !<--READ Information for atoms!                                   
                                                                        
      DO layer=1,layers%num_layers 
       ALLOCATE(atoms(layer)%nlo(atoms(layer)%ntype)) 
       ALLOCATE(atoms(layer)%lmax(atoms(layer)%ntype)) 
       ALLOCATE(atoms(layer)%lmax0(atoms(layer)%ntype)) 
       ALLOCATE(atoms(layer)%lapw_l(atoms(layer)%ntype)) 
       ALLOCATE(atoms(layer)%lda_u(atoms(layer)%ntype)) 
       ALLOCATE(atoms(layer)%lnonsph(atoms(layer)%ntype)) 
       ALLOCATE(tnoco(layer)%alph(atoms(layer)%ntype)) 
       ALLOCATE(tnoco(layer)%beta(atoms(layer)%ntype)) 
       ALLOCATE(tnoco(layer)%l_relax(atoms(layer)%ntype)) 
       ALLOCATE(tnoco(layer)%b_con(2,atoms(layer)%ntype)) 
!       tnoco(layer)%l_noco = .FALSE.                                   
!       tnoco(layer)%l_ss   = .FALSE.                                   
       DO n=1,atoms(layer)%ntype 
                     !READ empty line !Layer and atom type              
         READ(92,*) 
         READ(92,'(5x,i2,9x,i2,8x,i2,7x,i2)') atoms(layer)%lmax(n),     &
     &        atoms(layer)%lnonsph(n),atoms(layer)%lapw_l(n),           &
     &        atoms(layer)%lmax0(n)                                     
         IF(atoms(layer)%lmax0(n)==0)atoms(layer)%lmax0(n)=             &
     &          atoms(layer)%lmax(n)                                    
         IF (atoms(layer)%lmax(n)>atoms(layer)%lmax0(n))                &
     &        CPP_error("lmax can not be larger than lmax0")
         READ(92,'(4x,i2,5x,20i2)')atoms(layer)%nlo(n),(llo_tmp(i,n),i=1&
     &        ,atoms(layer)%nlo(n))                                     
         l=-1 
         READ(92,NML=ldaU) 
         atoms(layer)%lda_u(n)%l=l 
         atoms(layer)%lda_u(n)%u=u 
         atoms(layer)%lda_u(n)%j=j 
         !read the noco-line                                            
         alpha = 0.0;beta = 0.0;b_con = 0.0;relax=.FALSE. 
         IF (tnoco(layer)%l_noco) READ(92,NML = NOCOL) 
         tnoco(layer)%alph(n)    = alpha 
         tnoco(layer)%beta(n)    = beta 
         tnoco(layer)%b_con(:,n) = b_con 
         tNoco(layer)%l_relax(n) = relax 
       ENDDO 
       atoms(layer)%nlod=max(maxval(atoms(layer)%nlo),1) 
       atoms(layer)%llod=maxval(llo_tmp) 
       ALLOCATE(atoms(layer)%llo(atoms(layer)%nlod,atoms(layer)%ntype)  &
     &     ,atoms(layer)%lo1l(0:atoms(layer)%llod,atoms(layer)%ntype)   &
     &     ,atoms(layer)%l_dulo(atoms(layer)%nlod,atoms(layer)%ntype)   &
     &     ,atoms(layer)%ulo_der(atoms(layer)%nlod,atoms(layer)%ntype)  &
     &     ,atoms(layer)%nlol(0:atoms(layer)%llod,atoms(layer)%ntype))  
       atoms(layer)%llo=llo_tmp(:atoms(layer)%nlod,:atoms(layer)%ntype) 
                                                                        
       IF (maxval(atoms(layer)%nlo)>0) THEN 
         DO n=1,atoms(layer)%ntype 
            CALL fleur_setlomap(n,atoms(layer)) 
         ENDDO 
         !determine the total number of lo's : nlotot                   
         atoms(layer)%nlotot = 0 
         DO n = 1, atoms(layer)%ntype 
            DO l = 1,atoms(layer)%nlo(n) 
               atoms(layer)%nlotot =                                    &
     &           atoms(layer)%nlotot+                                   &
     &           atoms(layer)%neq(n)*(2*atoms(layer)%llo(l,n)+1)        
            ENDDO 
         ENDDO 
       ELSE 
         atoms(layer)%llod=0 
         atoms(layer)%nlotot = 0 
       ENDIF 
       ALLOCATE(lapw(layer)%kveclo(atoms(layer)%nlotot)) 
                                                                        
      !>                                                                
      !<--determine convergence parameter for pseudo-charge             
                                                                        
       CALL fleur_convndim(                                             &
     &     stars(layer)%gmax_inp,                                       &
     &     ncvd)                                                        
       ALLOCATE(atoms(layer)%ncv(atoms(layer)%ntype)) 
       CALL fleur_convn(                                                &
     &     ncvd,atoms(layer)%ntype,stars(layer)%gmax_inp,               &
     &     atoms(layer)%rmt,                                            &
     &     atoms(layer)%ncv)                                            
                                                                        
      !>                                                                
                                                                        
      !<-- Perform some checks on input Data&Compiler options           
                                                                        
!     Check for LDA+U:                                                  
                                                                        
      atoms(layer)%n_u = 0 
      DO  n = 1,atoms(layer)%ntype 
         IF (atoms(layer)%lda_u(n)%l>=0)  THEN 
            atoms(layer)%n_u = atoms(layer)%n_u + 1 
            IF (atoms(layer)%nlo(n)>=1) THEN 
               DO i = 1, atoms(layer)%nlo(n) 
                  IF ((abs(atoms(layer)%llo(i,n))                       &
     &                 ==atoms(layer)%lda_u(n)%l) .AND.               &
     &                 (.NOT.atoms(layer)%l_dulo(i,n)) ) WRITE (*,*)    &
     &                 'LO and LDA+U for same l not implemented'        
               ENDDO 
            ENDIF 
         ENDIF 
      ENDDO 
                                                                        
            !layer                                                      
      ENDDO 
                                                                        
      !>                                                                
                                                                        
      !>                                                                
      CLOSE(92) 
                                                                        
      !>end of reading                                                  
                                                                        
      !<--Set up LO's                                                   
                                                                        
      IF ( gfinp%curr/=0 .AND. gfinp%l_cbs ) CPP_error('embtmat: l_CBS&curr')

      CPP_juDFT_timestop_info("Reading the gf_inp")
      !>                                                                
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S: GF_WGFINP                                                  
                                                                        
      SUBROUTINE GF_WGFINP_short(filename,                              &
     &     gfinp,atoms,xcpot,mix,tsoc,fermi,stars,tnoco,rkmax_inp)      
                                                                        
!***********************************************************************
!     This subroutine:                                                  
!       - reads the gf_inp file                                         
!       - set some variables to defaults                                
!                                                                       
!                                          Daniel Wortmann, Juelich 2003
!***********************************************************************
      USE m_gf_types 
      IMPLICIT NONE 
      !<--  Scalar Arguments ..                                         
                                                                        
      CHARACTER, INTENT(IN)      ::filename*(*) 
      TYPE(t_gfinp), INTENT(IN)  ::gfinp 
      TYPE(t_atoms),INTENT(IN)   ::atoms(:) 
      TYPE(t_xcpot),INTENT(IN)   ::xcpot 
      TYPE(t_mix),INTENT(IN)     ::mix 
      TYPE(t_soc),INTENT(IN)     ::tsoc 
      TYPE(t_fermi),INTENT(IN)   ::fermi 
      TYPE(t_stars),INTENT(IN)   ::stars 
      TYPE(t_noco), INTENT(IN)   :: tnoco 
      REAL   ,INTENT(IN)         :: rkmax_inp 
                                                                        
      !>                                                                
      !<-- locals!                                                      
                                                                        
      !mode switches                                                    
      LOGICAL::tmat,gmat,charge 
      NAMELIST /MODE/ tmat,gmat,charge 
      !CBS                                                              
      LOGICAL::CBS 
      REAL   ::d1,d2 
      NAMELIST /CBSMODE/ CBS,d1,d2 
      !EMB                                                              
      LOGICAL::addemb 
      NAMELIST /EMB/addemb 
      !SC                                                               
      INTEGER::imix,maxiter,iter 
      REAL::alpha,spinf 
      LOGICAL::gauss,tria 
      REAL::tkb,delgau 
      NAMELIST /SC/iter,imix,maxiter,alpha,spinf,gauss,tria             &
     &     ,tkb,delgau                                                  
                                                                        
      !for XC-pot                                                       
      CHARACTER*4 ::namex 
      !for atom-loop                                                    
      INTEGER::n,i,layer 
                                                                        
      !>                                                                
      CPP_juDFT_timestop("Writing the gf_inp")
      !<--OPEN the gf_inp-file                                          
                                                                        
      OPEN(92,FILE=filename,FORM='formatted',STATUS='REPLACE') 
      !WRITE the mode-switches                                          
      tmat=gfinp%l_tmat 
      charge=gfinp%l_charge 
      gmat=gfinp%l_gmat 
      WRITE(92,NML=MODE) 
      !WRITE the CBS                                                    
      CBS=gfinp%l_CBS;d1=gfinp%dp1;d2                                   &
     &     =gfinp%dp2;                                                  
      WRITE(92,NML=CBSMODE) 
      !WRITE the BASIS                                                  
      WRITE(92,*) "&BASIS" 
      DO n = 1,size(gfinp%napw) 
         WRITE(92,"(a6,i0,a4,i0,a1)") " napw(",n,") = ",gfinp%napw(n)   &
     &        ,","                                                      
      ENDDO 
      WRITE(92,*) "rkmax=",rkmax_inp,"/" 
      !WRITE gmax                                                       
      WRITE(92,'(a5,f10.5,a8,f10.5)') 'gmax=',stars%gmax_inp,'gmaxxc='  &
     &     ,xcpot%gmaxxc                                                
      !WRITE EMB                                                        
      addemb=gfinp%l_addemb 
      WRITE(92,NML=EMB) 
      !WRITE SC                                                         
                                                                        
      imix = mix%imix;iter=mix%iter;maxiter=mix%maxiter 
      alpha = mix%alpha;spinf = mix%spinf 
      gauss=fermi%gauss;tria=fermi%tria 
      tkb=fermi%tkb;delgau=fermi%delgau 
      WRITE(92,NML=SC) 
      !WRITE empty line                                                 
      WRITE(92,*) '*****XC-POT******' 
      !WRITE XCPOT                                                      
      IF (xcpot%icorr==-1) namex='l91 ' 
      IF (xcpot%icorr==1) namex='wign' 
      IF (xcpot%icorr==2) namex='mjw ' 
      IF (xcpot%icorr==3) namex='bh  ' 
      IF (xcpot%icorr==4) namex='vwn ' 
      IF (xcpot%icorr==5) namex='pz  ' 
      IF (xcpot%icorr==6) namex='pw91' 
      IF (xcpot%icorr==7) namex='pbe ' 
      IF (xcpot%icorr==8) namex='rpbe' 
      IF (xcpot%icorr==9) namex='Rpbe' 
      WRITE (92,'(a4,a5,l1)') namex,',all=',.FALSE. 
      !WRITE Information for atoms!                                     
      DO layer = 1,SIZE(atoms) 
         DO n  = 1,atoms(layer)%ntype 
            !WRITE empty line                                           
            WRITE(92,*) '****Layer:',layer,' Atom:',n,'***' 
            WRITE(92,'(a5,i2,a9,i2,a8,i2)') 'lmax =',atoms(layer)%lmax(n&
     &           ),',lnonsph =',atoms(layer)%lnonsph(n),',lapw_l ='     &
     &           ,atoms(layer)%lapw_l(n)                                
            WRITE(92,'(a4,i2,a5,20i2)') 'nlo =',atoms(layer)%nlo(n)     &
     &           ,',llo =',(atoms(layer)%llo(i,n),i = 1,atoms(layer     &
     &           )%nlo(n))                                              
 8180       FORMAT ('&ldaU l =',i3,',u=',f6.3,',j=',f6.3,' /') 
            IF (atoms(layer)%lda_u(n)%l >= 0) THEN 
               WRITE(92,8180) atoms(layer)%lda_u(n)%l,atoms(layer       &
     &              )%lda_u(n)%u,atoms(layer)%lda_u(n)%j                
            ELSE 
               WRITE(92,8180) -1,0.,0. 
            ENDIF 
         ENDDO 
      ENDDO 
      CLOSE(92) 
      CPP_juDFT_timestop("Writing the gf_inp")
      !>                                                                
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
      END                                           
