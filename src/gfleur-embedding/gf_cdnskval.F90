!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_cdnskval 
      USE m_gf_ab_coef
      USE m_gf_pwden
      USE m_gf_rhocdnmt
      USE m_gf_types
      USE m_radfun
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: gf_cdnskval
#define cdebug                                                          
      CONTAINS 
      recursive SUBROUTINE gf_cdnskval(l_noco,fmpi,                     &
     &     lapw,lapw_gf,jspin,nk,nsymt,enpara,vr,                       &
     &     Gfkt,atoms,stars,sphhar,cell,sym,kpts,                       &
     &     qpw_s,rho_s,qmtl_s,nococonv,l_nomt,l_nopw,ityp)              
                                                                        
!*********************************************************************  
!     Generates the valence charge from the Green function for a single 
!     spin and k-point                                                  
!                                                                       
!     Not implemented:                                                  
!     Noco+LO's                                                         
!                                                                       
!                                      Daniel Wortmann                  
!*********************************************************************  
!     Implementation of local orbitals (in gf_rhocdnmt).                
!                                                                       
!                                      Frank Freimuth, April 2007       
!*********************************************************************  
                                                                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      LOGICAL, INTENT (IN) :: l_noco
      INTEGER, INTENT (IN) :: jspin 
      INTEGER, INTENT (IN) :: nsymt,nk 
                                                                        
      TYPE(t_sym),INTENT(IN)     :: sym 
      TYPE(t_atoms),INTENT(IN)   :: atoms 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_sphhar),INTENT(IN)  :: sphhar 
      TYPE(t_stars),INTENT(IN)   :: stars 
      TYPE(t_mpi),INTENT(IN)     :: fmpi
      TYPE(t_lapw),INTENT(IN)    :: lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      TYPE(t_enpara),INTENT(IN)  :: enpara 
      !noco: spin frames of the atom types; skip the MT / the interstitial part;
      !accumulate the MT density of one atom type only
      TYPE(t_nococonv),INTENT(IN),OPTIONAL :: nococonv
      LOGICAL,INTENT(IN),OPTIONAL :: l_nomt,l_nopw
      INTEGER,INTENT(IN),OPTIONAL :: ityp
      TYPE(t_kpts),INTENT(IN)    :: kpts 
      !>                                                                
                                                                        
      !<-- Array Arguments ..                                           
                                                                        
                                                                        
                                       !jmtd,ntypd,jspd                 
      REAL, INTENT    (IN) :: vr(:,:,:) 
      !green-function                                                   
      COMPLEX,INTENT  (IN) :: GFKT(:,:) 
      !charge-density(Output!)                                          
      COMPLEX,INTENT(INOUT):: qpw_s(:,:) 
                                             !jmtd,0:nlhd,ntypd,jspd)   
      REAL, INTENT (INOUT) :: rho_s(:,0:,:,:) 
                                           !lmaxd,ntype l-like charge   
      REAL, INTENT (INOUT) :: qmtl_s(0:,:) 
                                                                        
      !>                                                                
      !<--Locals                                                        
      REAL :: epar(0:size(enpara%el0,1)-1,                               &
     &         size(enpara%el0,2),size(enpara%el0,3))                     
                                                   !0:lmaxd,ntypd,jspd  
                                                !radial basis           
      REAL,ALLOCATABLE :: ddn(:,:),fg(:,:,:,:,:) 
                                                        !functions      
                                                        !loop-vars      
      INTEGER          :: itype,l 
      REAL,ALLOCATABLE :: duds(:,:),dus(:,:),uds(:,:),us(:,:) 
      REAL  wronk 
                                                !for radfun             
      INTEGER          :: noded,nodeu 
      TYPE(t_usdus) :: ud_loc
      COMPLEX,ALLOCATABLE ::  qpw(:)
      REAL   ,ALLOCATABLE ::  rho(:,:,:,:) 
      REAL   ,ALLOCATABLE :: qmtl(:,:)
      COMPLEX,ALLOCATABLE :: gs(:,:)
      integer             :: nv,nlo,it,ic,t1,t2
      COMPLEX             :: u(2,2),c
      COMPLEX,ALLOCATABLE :: qpw3(:)
      LOGICAL             :: l_pw,l_mt
      !>                                                                

      !noco: the interstitial density is calculated in the global spin
      !frame (diagonal blocks and the off-diagonal part), the MT density
      !of every atom type in its local spin frame
      if (l_noco) then
            nv=lapw_gf%nv_sphere(1)
            nlo=(lapw_gf%nmat_sphere-lapw_gf%nv_tot_sphere)/2
            allocate(gs(nv+nlo,nv+nlo))
            !interstitial, diagonal blocks
            do t1=1,2
               gs=0.0
               gs(:nv,:nv)=gfkt((t1-1)*nv+1:t1*nv,(t1-1)*nv+1:t1*nv)
               call gf_cdnskval(.false.,fmpi,                                   &
                 lapw,lapw_gf,t1,nk,nsymt,enpara,vr,                            &
                 gs,atoms,stars,sphhar,cell,sym,kpts,                           &
                 qpw_s,rho_s,qmtl_s,l_nomt=.true.)
            enddo
            !interstitial, off-diagonal part
            if (size(qpw_s,2)>=3) then
               allocate(qpw3(size(qpw_s,1)))
               qpw3=0.0
               call gf_pwden(lapw,lapw_gf,stars,1,cell%omtil,                   &
                 gfkt(1:nv,nv+1:2*nv),qpw3,gfkt(nv+1:2*nv,1:nv))
               qpw_s(:,3)=qpw_s(:,3)+kpts%wtkpt(nk)*qpw3
               deallocate(qpw3)
            endif
            !MT: local spin channels of every atom type
            do it=1,atoms%ntype
               u=0.0;u(1,1)=1.0;u(2,2)=1.0
               if (present(nococonv)) u=nococonv%umat(it)
               do ic=1,2
                  gs=0.0
                  do t1=1,2
                     do t2=1,2
                        c=CONJG(u(t1,ic))*u(t2,ic)
                        if (ABS(c)<1.E-12) cycle
                        gs(:nv,:nv)=gs(:nv,:nv)+c*gfkt((t1-1)*nv+1:t1*nv,(t2-1)*nv+1:t2*nv)
                        if (nlo>0) then
                           gs(:nv,nv+1:)=gs(:nv,nv+1:)+c*gfkt((t1-1)*nv+1:t1*nv,2*nv+(t2-1)*nlo+1:2*nv+t2*nlo)
                           gs(nv+1:,:nv)=gs(nv+1:,:nv)+c*gfkt(2*nv+(t1-1)*nlo+1:2*nv+t1*nlo,(t2-1)*nv+1:t2*nv)
                           gs(nv+1:,nv+1:)=gs(nv+1:,nv+1:)+c*gfkt(2*nv+(t1-1)*nlo+1:2*nv+t1*nlo,2*nv+(t2-1)*nlo+1:2*nv+t2*nlo)
                        endif
                     enddo
                  enddo
                  call gf_cdnskval(.false.,fmpi,                                &
                    lapw,lapw_gf,ic,nk,nsymt,enpara,vr,                         &
                    gs,atoms,stars,sphhar,cell,sym,kpts,                        &
                    qpw_s,rho_s,qmtl_s,l_nopw=.true.,ityp=it)
               enddo
            enddo
            return
      endif
      l_pw=.TRUE.
      IF (PRESENT(l_nopw)) l_pw=.NOT.l_nopw
      l_mt=.TRUE.
      IF (PRESENT(l_nomt)) l_mt=.NOT.l_nomt

                                                                        
      !<-- debug write out Greenfunction                                
!debug       INTEGER             :: n,nn                                
!debug       character(len = 20) :: filename                            
!debug       write(filename,"(a,i3,i1)") "GFKT.",nk,jspin               
!debug       open(99,file = filename)                                   
!debug       DO n = 1,size(gfkt,1)                                      
!debug          DO nn = 1,size(gfkt,2)                                  
!debug             write(99,"(2(i6,1x),3(f0.8,1x))") n,nn,gfkt(n,nn)    
!debug          enddo                                                   
!debug       enddo                                                      
!debug       close(99)                                                  
      !>                                                                
      ALLOCATE(qpw(size(qpw_s,1)))
      ALLOCATE(rho(size(rho_s,1),0:size(rho_s,2)-1,size(rho_s,3)        &
     &     ,size(rho_s,4)))                                             
      ALLOCATE(qmtl(0:size(qmtl_s,1)-1,size(qmtl_s,2))) 
      qpw = 0.0 
      rho = 0.0 
      qmtl = 0.0 
      epar=enpara%el0 
                                                                        
      !<-- Calculate interstitial charge                                

      IF (l_pw) THEN
         CALL timestart("gf_pwden")
         CALL gf_pwden(                                                 &
     &        lapw,lapw_gf,stars,jspin,cell%omtil,                      &
     &        GFKT,                                                     &
     &        qpw)
         CALL timestop("gf_pwden")
         qpw_s(:,jspin) = qpw_s(:,jspin)+kpts%wtkpt(nk)*qpw
      ENDIF
      deallocate(qpw)
      if (jspin>2 .OR. .NOT.l_mt) return ! No calculation of 3rd component for MTs
                                                                        
      !<--Calculate MT charge                                           

      CALL timestart("gf_rhocdnmt") 
                                                                        
      !<-- Set up radial functions                                      
                                                                        
      ALLOCATE(fg(MAXVAL(atoms%jri),2,0:MAXVAL(atoms%lmax),atoms%ntype,2&
     &     ))                                                           
      ALLOCATE(ddn(0:MAXVAL(atoms%lmax),atoms%ntype))
      ALLOCATE(us(0:MAXVAL(atoms%lmax),atoms%ntype))
      ALLOCATE(dus(0:MAXVAL(atoms%lmax),atoms%ntype))
      ALLOCATE(uds(0:MAXVAL(atoms%lmax),atoms%ntype))
      ALLOCATE(duds(0:MAXVAL(atoms%lmax),atoms%ntype))
      DO itype = 1,atoms%ntype 
         IF (.NOT.ALLOCATED(ud_loc%us)) CALL ud_loc%init(atoms,jspin)
         DO l = 0,atoms%lmax(itype)
            CALL radfun(l,itype,jspin,epar(l,itype,jspin),              &
     &           vr(:,itype,jspin),atoms,                               &
     &           fg(:,:,l,itype,1),fg(:,:,l,itype,2),                   &
     &           ud_loc,nodeu,noded,wronk)
            us(l,itype) = ud_loc%us(l,itype,jspin)
            dus(l,itype) = ud_loc%dus(l,itype,jspin)
            uds(l,itype) = ud_loc%uds(l,itype,jspin)
            duds(l,itype) = ud_loc%duds(l,itype,jspin)
            ddn(l,itype) = ud_loc%ddn(l,itype,jspin)
         ENDDO 
      ENDDO 
                                                                        
                                                                        
                                                                        
      !>                                                                
      CALL gf_rhocdnmt(atoms,fmpi,lapw,lapw_gf,enpara,                 &
     &     gfkt,sphhar,jspin,cell,sym,kpts%bk(:,nk),                    &
     &     fg(:,:,:,:,1),fg(:,:,:,:,2),ddn,vr,                  &
     &     us,dus,uds,duds,                                             &
     &     rho(:,0:,:,:),qmtl)                                          
      CALL timestop("gf_rhocdnmt") 

      !>                                                                
                                                                        
      DEALLOCATE(fg,ddn,us,dus,uds,duds) 
                                                                        

      IF (PRESENT(ityp)) THEN
         rho_s(:,:,ityp,:) = rho_s(:,:,ityp,:)+kpts%wtkpt(nk)*rho(:,:,ityp,:)
         qmtl_s(:,ityp) = qmtl_s(:,ityp)+kpts%wtkpt(nk)*qmtl(:,ityp)
      ELSE
         rho_s = rho_s+kpts%wtkpt(nk)*rho 
         qmtl_s = qmtl_s+kpts%wtkpt(nk)*qmtl 
      ENDIF
      DEALLOCATE(rho,qmtl)
      END SUBROUTINE gf_cdnskval 
      END                                           
