!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_types 
!*************************************************************          
!     This module contains definitions for all kind of types            
!*************************************************************          
                                                                        
      USE m_types,ONLY:t_utype,t_tlmplm,t_raddata=>t_usdus
      IMPLICIT NONE
                                                                        
      !!Add two constants for left and right                            
      INTEGER,PARAMETER :: GF_L=1 
      INTEGER,PARAMETER :: GF_R=2 
                                                                        
                                  !this type stores all information abou
      TYPE t_kpts 
                                  !no                                   
       INTEGER ::nkpts 
                                  !(3,nkpts) k-vectors internal units   
       REAL,ALLOCATABLE ::bk(:,:) 
                                  !(nkpts) weights                      
       REAL,POINTER ::weight(:) 
      ENDTYPE 
                                                                        
                                        !enpara information             
      TYPE t_enpara 
                                        !mixing-parameters for both spin
       REAL ::enmix(2) 
                                        !change for l,type,spin         
       LOGICAL,POINTER ::lchange(:,:,:) 
                                        !change for lo,type,spin        
       LOGICAL,POINTER ::llochg(:,:,:) 
                                        !energy-param for l,type,spin   
       REAL,POINTER    ::el(:,:,:) 
                                        !energy-param for lo,type,spin  
       REAL,POINTER    ::ello(:,:,:) 
                                        !for type,spin                  
       INTEGER,POINTER ::skiplo(:,:) 
                                        !change vac,spin                
       LOGICAL         ::lchg_v(2,2) 
                                        !vacuum energy param            
       REAL            ::evac(2,2) 
      ENDTYPE 
                                                                        
                    ! this type stores all the information about atoms  
      TYPE t_atoms 
                                       !no of types                     
       INTEGER ::ntype 
                                       !total-no of atoms               
       INTEGER ::nat 
                                       !dimensions of LO's              
       INTEGER ::nlod,llod,nlotot 
       INTEGER ::n_u 
                                       !No of element                   
       INTEGER,POINTER ::nz(:) 
                                       !atoms per type                  
       INTEGER,POINTER::neq(:) 
                                       !radial grid points              
       INTEGER,POINTER::jri(:) 
                                       !core states                     
       INTEGER,POINTER::ncst(:) 
                                       !a larger lmax for kin.energy    
       INTEGER,POINTER::lmax0(:) 
                                       !lmax                            
       INTEGER,POINTER::lmax(:) 
                                       !lmax non-spherical              
       INTEGER,POINTER::lnonsph(:) 
                                       !expansion of pseudo-charge      
       INTEGER,POINTER::ncv(:) 
                                       !no of LO                        
       INTEGER,POINTER::nlo(:) 
                                       !l of LO (nlo,ntype)             
       INTEGER,POINTER::llo(:,:) 
                                       !lmax for lapw (ntype)           
       INTEGER,POINTER::lapw_l(:) 
                                       !first LO with a given l (max(nlo
       INTEGER,POINTER::lo1l(:,:) 
                                       !??                              
       INTEGER,POINTER::ulo_der(:,:) 
                                       !no of LOs per l (max(nlo1),ntype
       INTEGER,POINTER::nlol(:,:) 
                                       !true if LO is formed by \dot u (
       LOGICAL,POINTER::l_dulo(:,:) 
                                       !no of op that maps atom into    
       INTEGER,POINTER::ngopr(:) 
                                       !represent. (nat)                
                                       !symetry of atom (nat???)        
       INTEGER,POINTER::ntypsy(:) 
                                       !no of sphhar for atom type(ntype
       INTEGER,POINTER ::nlhtyp(:) 
                                       !atom mapped to by inversion (nat
       INTEGER,POINTER ::invsat(:) 
                                       !MT-Radius (ntype)               
       REAL,POINTER::rmt(:) 
                                       !log increment(ntype)            
       REAL,POINTER::dx(:) 
                                       !vol of MT(ntype)                
       REAL,POINTER::volmts(:) 
                                       !radial grid points(max(jri),ntyp
       REAL,POINTER::rmsh(:,:) 
                                       !charge of nucleus(ntype)        
       REAL,POINTER::zatom(:) 
                                       !initial mag moment(ntype)       
       REAL,POINTER::bmu(:) 
                                       !pos of atom (absol) (3,nat)     
       REAL,POINTER::pos(:,:) 
                                       !pos of atom (relat)(3,nat)      
       REAL,POINTER::taual(:,:) 
                                       !lda_u information(ntype)        
       TYPE(t_utype),POINTER::lda_u(:) 
      END TYPE 
                                                                        
                          !Unit-cell information                        
      TYPE t_cell 
                          !name of 2D-lattice type                      
       CHARACTER*3::latnam 
                          !vol of dtilde box                            
       REAL::omtil 
                          !2D area                                      
       REAL::area 
                          !bravais matrix                               
       REAL::amat(3,3) 
                          !rez. bravais matrx                           
       REAL::bmat(3,3) 
                          !square of bbmat                              
       REAL::bbmat(3,3) 
                          !d-value                                      
       REAL::z1 
                          !volume of cell                               
       REAL::vol 
                          !volume of interstitial                       
       REAL::volint 
       REAL:: c
      END TYPE 
                                                                        
                                  !The stars                            
      TYPE t_stars 
                                       !max-length of star              
        REAL :: gmax 
                                       !cutoff of gmax                  
        REAL :: gmax_inp 
                                       !cutoff for correction of        
        REAL :: gmax_decouple 
                                       !coul-pot                        
                                       !cutoff for correction of        
        REAL :: gmax_pot 
                                       !coul-pot                        
                                       !no of 3d-stars                  
        INTEGER ::nq3 
                                       !no of 2d-stars                  
        INTEGER ::nq2 
                                       !No of elements in FFT           
        INTEGER ::kimax 
                                       !No of elements in 2D-FFT        
        INTEGER ::kimax2 
                                       !dim of box                      
        INTEGER ::mx1,mx2,mx3 
                                       !No of elements in z-direction   
        INTEGER ::ngz,izmin,izmax 
                                       !rep. g-vector of star           
        INTEGER,POINTER ::kv3(:,:) 
                                       !length of star                  
        REAL,POINTER    ::sk3(:) 
                                       !mapping of g-vectors to stars   
        INTEGER,POINTER ::ig(:,:,:) 
                                       !No of g-vectors in star         
        INTEGER,POINTER ::nstr(:) 
                                       !rep. g-vector of 2D-star        
        INTEGER,POINTER ::kv2(:,:) 
                                       !length of 2D-star               
        REAL,POINTER    ::sk2(:) 
                                       !No of g-vecs in 2D-star         
        INTEGER,POINTER ::nstr2(:) 
                                       !mapping of                      
        INTEGER,POINTER ::ig2(:) 
                                       !                                
        INTEGER,POINTER ::igz(:) 
                                       !phase phactor of g-vector       
        REAL,POINTER    ::rgphs(:,:,:) 
                                       !mapping of stars to FFT-box     
        INTEGER, POINTER :: igfft(:,:) 
                                       !same for 2D                     
        INTEGER, POINTER :: igfft2(:,:) 
                                       !phasefactors for mapping        
        REAL,   POINTER  :: pgfft(:) 
                                       !same of 2D                      
        REAL,   POINTER  :: pgfft2(:) 
                                       !                                
        REAL,   POINTER  :: pgft2xy(:) 
        REAL,   POINTER  :: pgft2x(:),pgft2y(:) 
        REAL,   POINTER  :: pgft2xx(:),pgft2yy(:) 
      END TYPE 
                                                                        
                           !Data for EX-Potential                       
      TYPE t_xcpot 
        INTEGER::icorr 
        INTEGER::igrd 
        INTEGER::krla 
        INTEGER::kxc1_fft,kxc2_fft,kxc3_fft,nxc3_fft,kmxxc_fft 
        INTEGER::kxc1d,kxc2d,kxc3d 
        INTEGER::idsprs0,idsprsl,idsprsi,idsprsv,iggachk 
        INTEGER:: idsprs,isprsv,ndvgrd 
        LOGICAL::lwb 
        REAL::sprsv,chng 
        REAL::gmaxxc 
      END TYPE 
                                                                        
      TYPE t_potential 
        REAL,ALLOCATABLE    :: vr(:,:,:,:) 
        COMPLEX,ALLOCATABLE :: vpw(:,:) 
      END TYPE 
                                                                        
      TYPE t_charge 
        COMPLEX,ALLOCATABLE :: pwd_new(:,:) 
        REAL,ALLOCATABLE    :: rho_new(:,:,:,:) 
        REAL,ALLOCATABLE    :: qmtl_new(:,:) 
      END TYPE 
                                                                        
                           !Mixing and self-conistency parameters       
      TYPE t_mix 
                           !no of iterations                            
        INTEGER :: iter 
                           !Mixing method                               
        INTEGER :: imix 
                           !Max no of iterations                        
        INTEGER :: maxiter 
                           !Mixing parameter                            
        REAL    :: alpha 
                           !Spin-mixing Parameter                       
        REAL    :: spinf 
        real    :: k_kerker,g0max,g0scale
        integer :: precond
      END TYPE 
                                                                        
                          !Data for the spherical harmonics             
      TYPE t_sphhar 
                                         !No of symmetry types (must    
        INTEGER ::ntypsd 
                                         !equal maxval(atoms%ntypsy)    
                                         !Max no of members of sphhar   
        INTEGER ::memd 
                                         !max of nlh                    
        INTEGER ::nlhd 
                                         !No of sphhar (ntypsd)         
        INTEGER,POINTER ::nlh(:) 
                                         !l's of sphhar (0:nlhd,ntypsd) 
        INTEGER,POINTER ::llh(:,:) 
                                         !No of members in sphhar (0:nlh
        INTEGER,POINTER ::nmem(:,:) 
                                         !lm's of of members (max(nmem),
        INTEGER,POINTER ::mlh(:,:,:) 
                                         !phasefactors (max(nmem),0:nlhd
        COMPLEX,POINTER ::clnu(:,:,:) 
      END TYPE 
                                                                        
                    !symmetry information                               
      TYPE t_sym 
                                        !Symophic group                 
       LOGICAL ::symor 
                                        !2D-inv-sym                     
       LOGICAL ::invs2 
                                        !Inversion-sym                  
       LOGICAL ::invs 
                                        !Z-refls. sym                   
       LOGICAL ::zrfs 
                                        !No of sym ops                  
       INTEGER ::nop 
                                        !No of 2D-sym ops               
       INTEGER ::nop2 
!       INTEGER,POINTER ::mrot(:,:,:)    !Rot-matrices (3,3,nop)        
!       INTEGER,POINTER ::invtab(:)      !inverse operation (nop)       
!       REAL,POINTER   ::tau(:,:)        !translation vectors (3,nop)   
       INTEGER,POINTER::mrot(:,:,:) 
       INTEGER,POINTER::invtab(:) 
       REAL,POINTER::tau(:,:) 
                                        !Name of lattice type           
       CHARACTER*3   :: latnam 
                                        !Name of sym                    
       CHARACTER*4   :: namgrp 
      END TYPE 
                                                                        
      TYPE t_gfinp 
       LOGICAL :: l_gf,l_eigen,l_surface,l_inp2plot,l_savemem,l_writeHS 
       LOGICAL :: l_tmat,l_gmat,l_CBS,l_solwil,l_gproj,l_band,l_embmt 
       LOGICAL :: l_charge,l_writeT,l_addemb,l_dos,l_potmix,l_adv,l_IEC 
       LOGICAL :: l_totalmix,l_hdfio
                             !use spectral representation to get Green f
       LOGICAL :: l_spectral ,l_embspectral,l_simplevacuum
                              ! calculate full-green function instead of
       LOGICAL :: l_fullgreen 
                              ! only projections                        
       LOGICAL :: l_nogno 
       LOGICAL :: l_nohelpregion 
       LOGICAL :: l_addselfen 
       LOGICAL :: l_intdos
       CHARACTER(LEN = 4) ::kpts 
                                ! c1,c2, info about embedding plane     
       REAL   :: dp1,dp2,CBS_bz 
       REAL   :: eps_current,eps_non_bloch,charge_limit,kappa_max 

       REAL   :: imag_broad   !broading by imaginary energy
                            ! should a Efield be applied.               
       REAL   :: Efield ,vacuum_energy
                             !min and max of z-positions of dos         
       REAL   :: z_min,z_max 
                             !no of dos-planes                          
       INTEGER :: nz_dos 
       INTEGER::npw,curr 
       INTEGER,ALLOCATABLE             :: napw(:) 
       LOGICAL,ALLOCATABLE             :: l_doslayer(:)
       REAL   ::trans 
       REAL   ::c_sc 
  !     REAL   ::bias                                                   
      END TYPE t_gfinp 
                                                                        
                                                                        
      TYPE t_layers 
      INTEGER :: num_layers 
      REAL,ALLOCATABLE    :: c1(:),c2(:),d(:),dt(:),c(:) 
      END TYPE t_layers 
                                                                        
      TYPE t_k 
       SEQUENCE 
       INTEGER,ALLOCATABLE:: k1(:,:),k2(:,:),k3(:,:),kp(:,:) 
       REAL,ALLOCATABLE::rk(:,:) 
      END TYPE t_k 
                                                                        
      TYPE t_kp 
       SEQUENCE 
       INTEGER,POINTER     :: k1p(:,:),k2p(:,:) 
       REAL,ALLOCATABLE    :: rkp(:,:) 
      END TYPE t_kp 
                                                                        
      TYPE t_lapw 
                                      !max-length of lapw for this k    
       REAL    :: rkmax 
                                      !input of rkmax                   
       REAL    :: rkmax_inp 
                                      !matrix dimension                 
       INTEGER :: nmat,nmat_sphere
                                      !max no of gvecs                  
       INTEGER :: nvd 
                                      !max no of 2D-gvecs               
       INTEGER :: nv2d 
                                      !total basis-size (2*nv2 for noco)
       INTEGER :: nv2_tot 
                                      !total basis-size (2*nv for noco) 
       INTEGER :: nv_tot,nv_tot_sphere
                                      !max no of basis functions        
       INTEGER :: nbasfcn 
                                      !actual no of g's for this k(jspin
       INTEGER :: nv(2),nv_sphere(2)
                                      !actual no of 2D-g's for this k   
       INTEGER :: nv2(2) 
                                          !dims of fft for charge genera
       INTEGER :: kq1_fft,kq2_fft,kq3_fft 
       INTEGER :: nq2_fft,nq3_fft,kmxq2_fft,kmxq_fft 
       LOGICAL :: l_cylinder 
       !The 2D g-vectors for this particular k-point can be mapped to a 
       !global list defined below                                       
       INTEGER,POINTER :: global2Dlist(:,:) 
       INTEGER,POINTER :: global2Dmap(:,:) 
       integer,allocatable :: g2map(:)
       INTEGER,ALLOCATABLE::kveclo(:) 
       INTEGER         :: g_MAX(2) 
                          !The k-values of the lapw                     
       TYPE(t_k)  :: k 
                           !the kp-values                               
       TYPE(t_kp) :: kp 
      END TYPE t_lapw 
                                                                        
                                                                        
      TYPE t_mpi 
       SEQUENCE 
       INTEGER :: irank,isize 
                                                         !old ...       
       INTEGER :: n_start,n_stride,n_rank,n_size,sub_comm
       LOGICAL :: pe0 
                                                                        
       !for k-parallelization                                           
       INTEGER             :: k_kperPE,k_PEperK,k_rank 
       INTEGER,ALLOCATABLE :: k_kpts(:) 
       INTEGER             :: com_world,self_subcom 
                                           !MPI-Com for IO              
       INTEGER             :: iodop_subcom 
                                               !MPI-COM of all PE       
       INTEGER             :: samelayer_subcom,samek_subcom
                                               !working on same layer   
       !for k+layer parallelization                                     
       INTEGER             :: kl_LayerPerPE 
       INTEGER,ALLOCATABLE :: kl_layers(:) 
       ! for k+en parallelization                                       
       INTEGER             :: ke_ENPerPE 
       INTEGER,ALLOCATABLE :: ke_energies(:) 
      END TYPE t_mpi 
                                                                        
      TYPE t_dims 
        INTEGER::irecl0 
        INTEGER::irecl1 
        INTEGER::jspins 
        INTEGER::n_u 
                        !maxval(abs(llo))                               
        INTEGER::llod 
      END TYPE t_dims 
                                                                        
      TYPE t_soc 
         LOGICAL::soc 
         REAL::theta,phi 
      END TYPE t_soc 
                                                                        
      TYPE t_noco 
      LOGICAL :: l_noco 
      LOGICAL, POINTER :: l_relax(:) 
      REAL, POINTER :: alph(:) 
      REAL, POINTER :: beta(:) 
      REAL, POINTER :: b_con(:,:) 
      LOGICAL :: l_ss,l_mperp,l_constr 
      REAL :: qss(3) 
      REAL :: mix_b 
      real :: alpha_int,beta_int
      END TYPE t_noco 
                                                                        
      TYPE t_fermi 
        LOGICAL ::gauss,tria 
        REAL    ::tkb,delgau 
      END TYPE t_fermi 
                                                                        
                                                                        
      !type for radial wavefunction data
      !Now type t_usdus in m_types.f from FLEUR
      !TYPE t_raddata
      !   REAL,ALLOCATABLE :: ulos(:,:,:),dulos(:,:,:),uulon(:,:,:)      &
     !&     ,dulon(:,:,:),uloulopn(:,:,:,:)
     !    REAL,ALLOCATABLE :: us(:,:,:),uds(:,:,:),dus(:,:,:),           &
     !&        duds(:,: ,:),ddn(:,:,:)
      !END TYPE
                                                                        
      !type for tlmplm
      ! Now in m_types.f from FLEUR
      !TYPE t_tlmplm
      !  COMPLEX, ALLOCATABLE :: tuu(:,:,:),tdd(:,:,:),tud(:,:,:)
      !  COMPLEX, ALLOCATABLE :: tuulo(:,:,:,:),tdulo(:,:,:,:)
      !  COMPLEX, ALLOCATABLE :: tdu(:,:,:),tuloulo(:,:,:,:)
      !  INTEGER, ALLOCATABLE :: ind(:,:,:,:)
      !END TYPE
                                                                        
      END                                           
