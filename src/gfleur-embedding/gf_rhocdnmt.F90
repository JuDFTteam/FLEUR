!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_rhocdnmt 
      use m_juDFT
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_rhocdnmt(atoms,lapw,enpara,                         &
     &     gmat,sphhar,jspin,cell,sym,bk,                               &
     &     f,g,ddn,vr,                                          &
     &     us,dus,uds,duds,                                             &
     &     rho,qmtl)                                                    
!********************************************************************** 
! This subroutine constructs the spherical and non-spherical            
! charge density inside the muffin-tin spheres. It can be considered as 
! the GF version of rhomt, rhonmt and cdnmt.                            
! Updated version, uses different loop structure+lapack for speed       
! increase of >800%                                                     
!                          Daniel Wortmann, Sun Feb 15 10:51:02 2004    
!********************************************************************** 
! Change of loop structure and implementation of level 3 BLAS (cgemm).  
! Speed-up: typically factor 30.                                        
!                               Frank Freimuth, March 2007              
!********************************************************************** 
! Implementation of local orbitals.                                     
!                                  Frank Freimuth, April 2007           
!********************************************************************** 
! Comment on the lnonsph-array: Its name is misleading! For consistency 
! with the formula used in standard-fleur, it should be set to lmax, as 
! in the code below. This comment refers only to the lnonsph-array      
! defined locally in this subroutine and not to the lnonsph-array within
! the atom-type.                                                        
!               Frank Freimuth, May 2008                                
!********************************************************************** 
! Usually the switch DCPP_LMPSYM should NOT be set, because this switch 
! requires the Green's function to satisfy G(x,x')=G(x',x), which is tru
! only for the full Green's function but not for the contribution to the
! Green's function of a certain k-point.                                
!                                       Frank Freimuth, June 2008       
!********************************************************************** 
      !<--Use                                                           
                                                                        
#include "cpp_double.h"
      USE m_fleur_interface,ONLY: fleur_gaunt 
      USE m_constants,ONLY:Pimach 
      USE m_gf_types 
      USE m_gf_rholocorbs 
      use m_gf_ab_coef
      IMPLICIT NONE 
                                                                        
      !>                                                                
      !<--Arguments                                                     
                                                                        
      INTEGER,INTENT(IN)       :: jspin 
                                                     ! The Green functio
      COMPLEX, INTENT(IN)      :: gmat(:,:) 
      TYPE(T_atoms),INTENT(IN) :: atoms 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_sphhar),INTENT(IN):: sphhar 
      TYPE(t_sym),INTENT(IN)   :: sym 
      TYPE(t_enpara),INTENT(IN):: enpara 
                                                     ! lapw info        
      !radial basis functions:                                          
      REAL,INTENT(IN)             :: ddn(0:,:),f(:,:,0:,:),g(:,:,0:,:) 
                                                    ! Kpt               
      REAL, INTENT(IN)            :: Bk(3) 
      TYPE(t_cell),INTENT(IN)     :: cell 
      REAL,INTENT(IN)     :: us(0:,:),dus(0:,:),uds(0:,:),duds(0:,:) 
      REAL,INTENT(IN)     :: vr(:,:,:) 
                                        !l-like charge                  
      REAL, INTENT(INOUT) :: qmtl(0:,:) 
      REAL,  INTENT(INOUT):: rho(:,0:,:,:) 
                                                                        
      !>                                                                
      !<--Locals                                                        
                                                                        
      INTEGER  :: l,lh,lm,lmp,lp,lv,j,jmem,lmxam,lmpmax,neqmax,bas1 
      INTEGER  :: m,mp,mv,na,itype,ns,nt,n,n_sph,lmax,lo,anglo,bas2 
      INTEGER  :: minat,maxat,lo1,lo2,locnum1,locnum2,locind1,locind2 
      INTEGER  :: nloitype 
      INTEGER,ALLOCATABLE:: basindex(:,:) 
                                 !result of sum over g,g'               
      COMPLEX  :: uu,ud,du,dd 
      COMPLEX  :: locorbuu,locorbdd,locorbdu,locorbud,locorblolo 
      COMPLEX  :: locorbulo,locorbdlo,locorblou,locorblod 
      COMPLEX  :: coef 
      REAL     :: sfp,gcr,gci,pi,norm 
      REAL,ALLOCATABLE    ::flo(:,:,:) 
      COMPLEX,ALLOCATABLE :: g_coef(:) 
      COMPLEX,ALLOCATABLE :: vecmat(:,:,:,:),leftlomat(:,:,:,:) 
      COMPLEX,ALLOCATABLE :: rightlomat(:,:,:,:) 
      REAL,ALLOCATABLE    :: s1(:),s2(:),s3(:),s4(:),s5(:),s6(:) 
      REAL,ALLOCATABLE    :: sfflo(:),sgflo(:),sfloflo(:) 
      REAL,ALLOCATABLE    :: sflof(:),sflog(:) 
      COMPLEX,ALLOCATABLE :: alocof(:,:,:),blocof(:,:,:),clocof(:,:,:) 
      INTEGER,ALLOCATABLE :: lh_sp(:) 
      LOGICAL             :: l_first,l_lolo,l_apwlo,l_loapw 
      REAL                :: uulon(atoms%nlod,atoms%ntype) 
      REAL                :: dulon(atoms%nlod,atoms%ntype) 
      REAL                :: uloulopn(atoms%nlod,atoms%nlod,atoms%ntype) 
      INTEGER             :: lnonsph(atoms%ntype) 
      COMPLEX CPP_BLAS_cdotc 
      EXTERNAL CPP_BLAS_cdotc 
                                                                        
                                                                        
                                                                        
!      print*,"lapw%nv(jspin)=",lapw%nv(jspin)                          
!      print*,"size(gmat,1)=",size(gmat,1)                              
!      print*,"size(gmat,2)=",size(gmat,2)                              
!      print*,"lapw%nmat=",lapw%nmat                                    
                                                                        
!      lnonsph(1:atoms%ntype)=atoms%lnonsph(1:atoms%ntype)              
      lnonsph(1:atoms%ntype)=atoms%lmax(1:atoms%ntype)
                                                                        
      !>                                                                
      sfp=2./SQRT(pimach())
      pi=pimach() 
      neqmax=maxval(atoms%neq(1:atoms%ntype)) 
      lmax=MAXVAL(atoms%lmax)
      lmxam=MAXVAL(lnonsph) 
      lmxam=lmxam*(lmxam+2) 
                                                                        
      ALLOCATE(g_coef((lmax**4+lmax**2)/8),lh_sp((lmax**4+lmax**2)/8)) 
      ALLOCATE(s5(maxval(atoms%jri))) 
      ALLOCATE(s6(maxval(atoms%jri))) 
      ALLOCATE(s1(maxval(atoms%jri))) 
      ALLOCATE(s2(maxval(atoms%jri))) 
      ALLOCATE(s3(maxval(atoms%jri))) 
      ALLOCATE(s4(maxval(atoms%jri))) 
                                                                        
                                                                        
      IF(atoms%nlotot>0)THEN 
          ALLOCATE(flo(maxval(atoms%jri),2,atoms%nlod)) 
          ALLOCATE(basindex(atoms%nat,atoms%nlod)) 
          ALLOCATE(alocof(atoms%nlotot,0:15,neqmax)) 
          ALLOCATE(blocof(atoms%nlotot,0:15,neqmax)) 
          ALLOCATE(clocof(atoms%nlotot,0:15,neqmax)) 
          ALLOCATE(sfflo(maxval(atoms%jri))) 
          ALLOCATE(sgflo(maxval(atoms%jri))) 
          ALLOCATE(sfloflo(maxval(atoms%jri))) 
          ALLOCATE(sflof(maxval(atoms%jri))) 
          ALLOCATE(sflog(maxval(atoms%jri))) 
      ENDIF 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
           !No of atom                                                  
      nt=0 
      !loop over atoms                                                  
      DO itype=1,atoms%ntype 
                               !ntypsy stores sym for atom not atomtype!
         ns=atoms%ntypsy(nt+1) 
         lmpmax=lnonsph(itype) 
         lmpmax=lmpmax*(lmpmax+2)+1 
         minat=nt+1 
         maxat=nt+atoms%neq(itype) 
         ALLOCATE(vecmat(lapw%nv_sphere(jspin),                                &
     &            0:lmpmax-1,1:atoms%neq(itype),1:2))                   
         CALL CPP_BLAS_cgemm('T','N',lapw%nv_sphere(jspin),                    &
     &              lmpmax*2*atoms%neq(itype),lapw%nv_sphere(jspin),           &
     &              cmplx(1.0,0.0),                                     &
     &              gmat,size(gmat,1),                                  &
     &              gf_ab_coef_matrix(lmpmax,minat,maxat,jspin),              &
     &              lapw%nv_sphere(jspin),                                    &
     &              cmplx(0.0,0.0),                                     &
     &              vecmat(:,0:lmpmax-1,1:atoms%neq(itype),1:2),        &
     &              lapw%nv_sphere(jspin))
                                                                        
         IF(atoms%nlo(itype)>=1)THEN 
            CALL gf_rholocorbs(atoms,itype,enpara,vr,jspin,             &
     &                 f(:,:,:,itype),g(:,:,:,itype),                   &
     &                 us,dus,uds,duds,ddn,lapw,cell,sym,bk,            &
     &                 flo,alocof,blocof,clocof,basindex,               &
     &                 uulon,dulon,uloulopn)                            
            nloitype=0 
            DO lo=1,atoms%nlo(itype) 
               l=1+2*atoms%llo(lo,itype) 
               nloitype=nloitype+l 
            ENDDO 
            ALLOCATE(  leftlomat(nloitype,                              &
     &                          0:lmpmax-1,1:atoms%neq(itype),1:2) )    
            ALLOCATE( rightlomat(nloitype,                              &
     &                          0:lmpmax-1,1:atoms%neq(itype),1:2) )    
            DO na=1,atoms%neq(itype) 
              CALL CPP_BLAS_cgemm('T','N',nloitype,                     &
     &              lmpmax*2,lapw%nv_sphere(jspin),                            &
     &              cmplx(1.0,0.0),                                     &
     &              gmat(1,lapw%nv_sphere(jspin)+basindex(nt+na,1)),           &
     &              size(gmat,1),                                       &
     &              gf_ab_coef_matrix(lmpmax,na+nt,na+nt,jspin),lapw%nv_sphere(jspin),    &
     &              cmplx(0.0,0.0),                                     &
     &              leftlomat(:,0:lmpmax-1,na,1:2),                     &
     &              nloitype)                                           
                                                                        
              CALL CPP_BLAS_cgemm('N','N',                              &
     &              nloitype,lmpmax*2,lapw%nv_sphere(jspin),                   &
     &              cmplx(1.0,0.0),                                     &
     &              gmat( (lapw%nv_sphere(jspin)+basindex(nt+na,1)):           &
     &                (lapw%nv_sphere(jspin)+basindex(nt+na,1)+nloitype-1),    &
     &                 1:lapw%nv_sphere(jspin)),                               &
     &              nloitype,                                           &
     &              conjg(gf_ab_coef_matrix(lmpmax,na+nt,na+nt,jspin)),             &
     &              lapw%nv_sphere(jspin),cmplx(0.0,0.0),                     &
     &              rightlomat(:,0:lmpmax-1,na,1:2),                    &
     &              nloitype)                                           
                  !na                                                   
            ENDDO 
               !nlo                                                     
         ENDIF 
                                                                        
                               !the first index of G                    
         DO l=0,lnonsph(itype) 
#ifdef CPP_LMPSYM                                                       
                     !second index of G                                 
          DO lp=0,l 
#else                                                                   
           DO lp=0,lnonsph(itype) 
#endif                                                                  
                            !l_first -> whether s has to be recalculated
            l_first=.TRUE. 
            l_lolo=.FALSE. 
            l_apwlo=.FALSE. 
            l_loapw=.FALSE. 
                                  !l,m                                  
            DO m=-l,l 
               lm=l*(l+1)+m 
                                         !lp,mp  (lp<=l)                
                  DO mp=-lp,lp 
                     lmp=lp*(lp+1)+mp 
#ifdef CPP_LMPSYM                                                       
                     IF (lmp>lm) CYCLE 
#endif                                                                  
                     !Generate a list of spherical harmonics            
                     n_sph=0 
                     DO lh=0,sphhar%nlh(ns) 
                        lv=sphhar%llh(lh,ns) 
                        IF (lv>lnonsph(itype)) CYCLE 
                                                       !l+lp+lv must    
                        IF (MOD(l+lp+lv,2) /= 0) CYCLE 
                                                       !be even,        
                                                       !cycle lh loop   
                        mv=m-mp 
                        !look if this spherical harmonics contains m    
                        jmem=0 
                        jloop:DO j=1,sphhar%nmem(lh,ns) 
                           IF (sphhar%mlh(j,lh,ns)==mv) THEN 
                              jmem=j 
                              EXIT jloop 
                           ENDIF 
                        ENDDO jloop 
                                           !m not found!                
                        IF (jmem==0) CYCLE 
                        coef=CONJG(sphhar%clnu(jmem,lh,ns))             &
     &                       *       fleur_gaunt(l,lv,lp,m,mv,mp        &
     &                       ,MAXVAL(atoms%lmax))
                                           !this test triangular eq.    
                        IF (coef==0) CYCLE 
                                           !found a spherical harmonic  
                        n_sph=n_sph+1 
#ifdef CPP_LMPSYM                                                       
                        IF (lmp<lm) coef=2.0*coef 
#endif                                                                  
                        g_coef(n_sph)=coef 
                        lh_sp(n_sph)=lh 
!                     WRITE(*,'(8i4,2f10.4)')itype,lh,l,m,lp,mp,lv      
!     +                    ,mv,coef                                     
                     ENDDO 
                                         !there are not sphhar for      
                     IF (n_sph==0) CYCLE 
                                         !this lm,lmp                   
                                                                        
                                                                        
                     !<-- multiply radial function                      
                     IF(l_first)THEN 
                      IF(atoms%nlo(itype)>0)THEN 
                        DO lo=1,atoms%nlo(itype) 
                           anglo=atoms%llo(lo,itype) 
                           IF(anglo==l)THEN 
                              IF(l_loapw) CALL juDFT_error("double-lo",calledby="gf_rhocdnmt.F90")
                              l_loapw=.TRUE. 
                              lo1=lo 
                           ENDIF 
                        ENDDO 
                        DO lo=1,atoms%nlo(itype) 
                           anglo=atoms%llo(lo,itype) 
                           IF(anglo==lp)THEN 
                              IF(l_apwlo) CALL juDFT_error("double-lo",calledby="gf_rhocdnmt.F90")
                              l_apwlo=.TRUE. 
                              lo2=lo 
                              IF(l_loapw)l_lolo=.TRUE. 
                           ENDIF 
                        ENDDO 
                            !local orbitals                             
                      ENDIF 
                      IF(l_lolo)THEN 
                       DO j = 1,atoms%jri(itype) 
                        s1(j)=f(j,1,l,itype)*f(j,1,lp,itype)+           &
     &                        f(j,2,l,itype)*f(j,2,lp,itype)            
                                                                        
                        s2(j)=g(j,1,l,itype)*g(j,1,lp,itype)+           &
     &                        g(j,2,l,itype)*g(j,2,lp,itype)            
                                                                        
                        s3(j)=f(j,1,l,itype)*g(j,1,lp,itype)+           &
     &                        f(j,2,l,itype)*g(j,2,lp,itype)            
                                                                        
                        s4(j)=g(j,1,l,itype)*f(j,1,lp,itype)+           &
     &                        g(j,2,l,itype)*f(j,2,lp,itype)            
                                                                        
                        sfflo(j)=f(j,1,l,itype)*flo(j,1,lo2)+           &
     &                          f(j,2,l,itype)*flo(j,2,lo2)             
                                                                        
                        sgflo(j)=g(j,1,l,itype)*flo(j,1,lo2)+           &
     &                           g(j,2,l,itype)*flo(j,2,lo2)            
                                                                        
                        sfloflo(j)=flo(j,1,lo1)*flo(j,1,lo2)+           &
     &                             flo(j,2,lo1)*flo(j,2,lo2)            
                                                                        
                        sflof(j)=f(j,1,lp,itype)*flo(j,1,lo1)+          &
     &                           f(j,2,lp,itype)*flo(j,2,lo1)           
                                                                        
                        sflog(j)=g(j,1,lp,itype)*flo(j,1,lo1)+          &
     &                           g(j,2,lp,itype)*flo(j,2,lo1)           
                       ENDDO 
                      ELSEIF(l_loapw)THEN 
                       DO j = 1,atoms%jri(itype) 
                        s1(j)=f(j,1,l,itype)*f(j,1,lp,itype)+           &
     &                        f(j,2,l,itype)*f(j,2,lp,itype)            
                                                                        
                        s2(j)=g(j,1,l,itype)*g(j,1,lp,itype)+           &
     &                        g(j,2,l,itype)*g(j,2,lp,itype)            
                                                                        
                        s3(j)=f(j,1,l,itype)*g(j,1,lp,itype)+           &
     &                        f(j,2,l,itype)*g(j,2,lp,itype)            
                                                                        
                        s4(j)=g(j,1,l,itype)*f(j,1,lp,itype)+           &
     &                        g(j,2,l,itype)*f(j,2,lp,itype)            
                                                                        
                        sflof(j)=f(j,1,lp,itype)*flo(j,1,lo1)+          &
     &                           f(j,2,lp,itype)*flo(j,2,lo1)           
                                                                        
                        sflog(j)=g(j,1,lp,itype)*flo(j,1,lo1)+          &
     &                           g(j,2,lp,itype)*flo(j,2,lo1)           
                       ENDDO 
                      ELSEIF(l_apwlo)THEN 
                       DO j = 1,atoms%jri(itype) 
                        s1(j)=f(j,1,l,itype)*f(j,1,lp,itype)+           &
     &                        f(j,2,l,itype)*f(j,2,lp,itype)            
                                                                        
                        s2(j)=g(j,1,l,itype)*g(j,1,lp,itype)+           &
     &                        g(j,2,l,itype)*g(j,2,lp,itype)            
                                                                        
                        s3(j)=f(j,1,l,itype)*g(j,1,lp,itype)+           &
     &                        f(j,2,l,itype)*g(j,2,lp,itype)            
                                                                        
                        s4(j)=g(j,1,l,itype)*f(j,1,lp,itype)+           &
     &                        g(j,2,l,itype)*f(j,2,lp,itype)            
                                                                        
                        sfflo(j)=f(j,1,l,itype)*flo(j,1,lo2)+           &
     &                           f(j,2,l,itype)*flo(j,2,lo2)            
                                                                        
                        sgflo(j)=g(j,1,l,itype)*flo(j,1,lo2)+           &
     &                           g(j,2,l,itype)*flo(j,2,lo2)            
                       ENDDO 
                      ELSE 
                       DO j = 1,atoms%jri(itype) 
                        s1(j)=f(j,1,l,itype)*f(j,1,lp,itype)+           &
     &                        f(j,2,l,itype)*f(j,2,lp,itype)            
                                                                        
                        s2(j)=g(j,1,l,itype)*g(j,1,lp,itype)+           &
     &                        g(j,2,l,itype)*g(j,2,lp,itype)            
                                                                        
                        s3(j)=f(j,1,l,itype)*g(j,1,lp,itype)+           &
     &                        f(j,2,l,itype)*g(j,2,lp,itype)            
                                                                        
                        s4(j)=g(j,1,l,itype)*f(j,1,lp,itype)+           &
     &                        g(j,2,l,itype)*f(j,2,lp,itype)            
                       ENDDO 
                            !l_lolo                                     
                      ENDIF 
                           !l_first                                     
                     ENDIF 
                     l_first=.FALSE. 
                                                                        
                                                                        
                     uu=0.0; dd=0.0 
                     ud=0.0; du=0.0 
                     !<--loop over equivalent atoms (for symmetry)      
                     DO na=1,atoms%neq(itype) 
                       n=nt+na 
                        !Use BLAS to construct the uu,ud,du,dd          
                       uu=uu+(CPP_BLAS_cdotc(lapw%nv_sphere(jspin),            &
     &                  gf_ab_coef_vector(lmp,n,1,jspin),1,vecmat(:,lm,na,1),1))
                       ud=ud+(CPP_BLAS_cdotc(lapw%nv_sphere(jspin),            &
     &                  gf_ab_coef_vector(lmp,n,2,jspin),1,vecmat(:,lm,na,1),1))
                       dd=dd+(CPP_BLAS_cdotc(lapw%nv_sphere(jspin),            &
     &                  gf_ab_coef_vector(lmp,n,2,jspin),1,vecmat(:,lm,na,2),1))
                       du=du+(CPP_BLAS_cdotc(lapw%nv_sphere(jspin),            &
     &                  gf_ab_coef_vector(lmp,n,1,jspin),1,vecmat(:,lm,na,2),1))
                                                                        
                     ENDDO 
                     !>                                                 
                                                                        
                     locorbuu=0.0;locorbdd=0.0 
                     locorbud=0.0;locorbdu=0.0 
                     locorbulo=0.0;locorbdlo=0.0 
                     locorblou=0.0;locorblod=0.0 
                     locorblolo=0.0; 
                                                                        
                     IF(l_lolo)THEN 
                      DO na=1,atoms%neq(itype) 
                       n=nt+na 
                       bas1=basindex(n,lo1) 
                       bas2=basindex(n,lo2) 
                       locnum1=2*atoms%llo(lo1,itype)+1 
                       locnum2=2*atoms%llo(lo2,itype)+1 
                       DO locind1=0,locnum1-1 
                        DO locind2=0,locnum2-1 
                                                                        
                         locorbuu=locorbuu+alocof(bas1+locind1,lm,na)*  &
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(alocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbdd=locorbdd+blocof(bas1+locind1,lm,na)*  &
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(blocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbud=locorbud+alocof(bas1+locind1,lm,na)*  &
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(blocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbdu=locorbdu+blocof(bas1+locind1,lm,na)*  &
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(alocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbulo=locorbulo+alocof(bas1+locind1,lm,na)*&
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(clocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbdlo=locorbdlo+blocof(bas1+locind1,lm,na)*&
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(clocof(bas2+locind2,lmp,na))      
                                                                        
                         locorblou=locorblou+clocof(bas1+locind1,lm,na)*&
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(alocof(bas2+locind2,lmp,na))      
                                                                        
                         locorblod=locorblod+clocof(bas1+locind1,lm,na)*&
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(blocof(bas2+locind2,lmp,na))      
                                                                        
                         locorblolo=locorblolo+                         &
     &                     clocof(bas1+locind1,lm,na)*                  &
     &                     gmat(lapw%nv_sphere(jspin)+locind1+bas1,            &
     &                          lapw%nv_sphere(jspin)+locind2+bas2)*           &
     &                          conjg(clocof(bas2+locind2,lmp,na))      
                                                                        
                              !locind2                                  
                        ENDDO 
                             !locind1                                   
                       ENDDO 
                            !na                                         
                      ENDDO 
                           !l_lolo                                      
                     ENDIF 
                                                                        
                     IF(l_loapw)THEN 
                      DO na=1,atoms%neq(itype) 
                       n=nt+na 
                       bas1=basindex(n,lo1) 
                       locnum1=2*atoms%llo(lo1,itype)+1 
                       DO locind1=0,locnum1-1 
                         locorbuu=locorbuu+alocof(bas1+locind1,lm,na)*  &
     &                          rightlomat(bas1+locind1+1-              &
     &                          basindex(n,1),lmp,na,1)                 
                                                                        
                         locorbdd=locorbdd+blocof(bas1+locind1,lm,na)*  &
     &                          rightlomat(bas1+locind1+1-              &
     &                          basindex(n,1),lmp,na,2)                 
                                                                        
                         locorbud=locorbud+alocof(bas1+locind1,lm,na)*  &
     &                          rightlomat(bas1+locind1+1-              &
     &                          basindex(n,1),lmp,na,2)                 
                                                                        
                         locorbdu=locorbdu+blocof(bas1+locind1,lm,na)*  &
     &                          rightlomat(bas1+locind1+1-              &
     &                          basindex(n,1),lmp,na,1)                 
                                                                        
                         locorblou=locorblou+clocof(bas1+locind1,lm,na)*&
     &                          rightlomat(bas1+locind1+1-              &
     &                          basindex(n,1),lmp,na,1)                 
                                                                        
                         locorblod=locorblod+clocof(bas1+locind1,lm,na)*&
     &                          rightlomat(bas1+locind1+1-              &
     &                          basindex(n,1),lmp,na,2)                 
                                                                        
                             !locind1                                   
                       ENDDO 
                            !na                                         
                      ENDDO 
                           !l_loapw                                     
                     ENDIF 
                                                                        
                     IF(l_apwlo)THEN 
                      DO na=1,atoms%neq(itype) 
                       n=nt+na 
                       bas2=basindex(n,lo2) 
                       locnum2=2*atoms%llo(lo2,itype)+1 
                       DO locind2=0,locnum2-1 
                         locorbuu=locorbuu+leftlomat( bas2+locind2+1-   &
     &                               basindex(n,1),lm,na,1 )*           &
     &                          conjg(alocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbdd=locorbdd+leftlomat( bas2+locind2+1-   &
     &                               basindex(n,1),lm,na,2 )*           &
     &                          conjg(blocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbud=locorbud+leftlomat( bas2+locind2+1-   &
     &                               basindex(n,1),lm,na,1 )*           &
     &                          conjg(blocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbdu=locorbdu+leftlomat( bas2+locind2+1-   &
     &                               basindex(n,1),lm,na,2 )*           &
     &                          conjg(alocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbulo=locorbulo+leftlomat( bas2+locind2+1- &
     &                               basindex(n,1),lm,na,1 )*           &
     &                          conjg(clocof(bas2+locind2,lmp,na))      
                                                                        
                         locorbdlo=locorbdlo+leftlomat( bas2+locind2+1- &
     &                               basindex(n,1),lm,na,2 )*           &
     &                          conjg(clocof(bas2+locind2,lmp,na))      
                             !locind2                                   
                       ENDDO 
                            !na                                         
                      ENDDO 
                           !l_apwlo                                     
                     ENDIF 
                                                                        
                     IF(l_lolo)THEN 
                       DO j=1,atoms%jri(itype) 
                          s5(j)=                                        &
     &                                 real(uu)*s1(j)+                  &
     &                                 real(dd)*s2(j)+                  &
     &                                 real(ud)*s3(j)+                  &
     &                                 real(du)*s4(j)+                  &
     &                                 real(locorbuu)*s1(j)+            &
     &                                 real(locorbdd)*s2(j)+            &
     &                                 real(locorbud)*s3(j)+            &
     &                                 real(locorbdu)*s4(j)+            &
     &                                 real(locorbulo)*sfflo(j)+        &
     &                                 real(locorbdlo)*sgflo(j)+        &
     &                                 real(locorblou)*sflof(j)+        &
     &                                 real(locorblod)*sflog(j)+        &
     &                                 real(locorblolo)*sfloflo(j)      
                          s6(j)=                                        &
     &                                 imag(uu)*s1(j)+                  &
     &                                 imag(dd)*s2(j)+                  &
     &                                 imag(ud)*s3(j)+                  &
     &                                 imag(du)*s4(j)+                  &
     &                                 imag(locorbuu)*s1(j)+            &
     &                                 imag(locorbdd)*s2(j)+            &
     &                                 imag(locorbud)*s3(j)+            &
     &                                 imag(locorbdu)*s4(j)+            &
     &                                 imag(locorbulo)*sfflo(j)+        &
     &                                 imag(locorbdlo)*sgflo(j)+        &
     &                                 imag(locorblou)*sflof(j)+        &
     &                                 imag(locorblod)*sflog(j)+        &
     &                                 imag(locorblolo)*sfloflo(j)      
                                                                        
                       ENDDO 
                     ELSEIF(l_loapw)THEN 
                       DO j=1,atoms%jri(itype) 
                          s5(j)=                                        &
     &                                 real(uu)*s1(j)+                  &
     &                                 real(dd)*s2(j)+                  &
     &                                 real(ud)*s3(j)+                  &
     &                                 real(du)*s4(j)+                  &
     &                                 real(locorbuu)*s1(j)+            &
     &                                 real(locorbdd)*s2(j)+            &
     &                                 real(locorbud)*s3(j)+            &
     &                                 real(locorbdu)*s4(j)+            &
     &                                 real(locorblou)*sflof(j)+        &
     &                                 real(locorblod)*sflog(j)         
                          s6(j)=                                        &
     &                                 imag(uu)*s1(j)+                  &
     &                                 imag(dd)*s2(j)+                  &
     &                                 imag(ud)*s3(j)+                  &
     &                                 imag(du)*s4(j)+                  &
     &                                 imag(locorbuu)*s1(j)+            &
     &                                 imag(locorbdd)*s2(j)+            &
     &                                 imag(locorbud)*s3(j)+            &
     &                                 imag(locorbdu)*s4(j)+            &
     &                                 imag(locorblou)*sflof(j)+        &
     &                                 imag(locorblod)*sflog(j)         
                                                                        
                       ENDDO 
                     ELSEIF(l_apwlo)THEN 
                       DO j=1,atoms%jri(itype) 
                          s5(j)=                                        &
     &                                 real(uu)*s1(j)+                  &
     &                                 real(dd)*s2(j)+                  &
     &                                 real(ud)*s3(j)+                  &
     &                                 real(du)*s4(j)+                  &
     &                                 real(locorbuu)*s1(j)+            &
     &                                 real(locorbdd)*s2(j)+            &
     &                                 real(locorbud)*s3(j)+            &
     &                                 real(locorbdu)*s4(j)+            &
     &                                 real(locorbulo)*sfflo(j)+        &
     &                                 real(locorbdlo)*sgflo(j)         
                          s6(j)=                                        &
     &                                 imag(uu)*s1(j)+                  &
     &                                 imag(dd)*s2(j)+                  &
     &                                 imag(ud)*s3(j)+                  &
     &                                 imag(du)*s4(j)+                  &
     &                                 imag(locorbuu)*s1(j)+            &
     &                                 imag(locorbdd)*s2(j)+            &
     &                                 imag(locorbud)*s3(j)+            &
     &                                 imag(locorbdu)*s4(j)+            &
     &                                 imag(locorbulo)*sfflo(j)+        &
     &                                 imag(locorbdlo)*sgflo(j)         
                                                                        
                       ENDDO 
                     ELSE 
                       DO j=1,atoms%jri(itype) 
                          s5(j)=                                        &
     &                                 real(uu)*s1(j)+                  &
     &                                 real(dd)*s2(j)+                  &
     &                                 real(ud)*s3(j)+                  &
     &                                 real(du)*s4(j)                   
                          s6(j)=                                        &
     &                                 aimag(uu)*s1(j)+                 &
     &                                 aimag(dd)*s2(j)+                 &
     &                                 aimag(ud)*s3(j)+                 &
     &                                 aimag(du)*s4(j)                  
                       ENDDO 
                     ENDIF 
                                                                        
                     !>                                                 
                     !<--Sum over sphhar again                          
                     norm=1.0/pi/atoms%neq(itype) 
                     DO n=1,n_sph 
                      lh=lh_sp(n) 
                      gci=aimag(g_coef(n))*norm 
                      gcr=real(g_coef(n))*norm 
                      DO j=1,atoms%jri(itype) 
                        rho(j,lh_sp(n),itype,jspin) = rho(j,lh_sp(n)    &
     &                       ,itype,jspin)+ (gcr*s6(j)+gci*s5(j))       
                                                                        
                      ENDDO 
                        !<--l-like charge                               
                      IF (lh==0) THEN 
                         IF(l_lolo)THEN 
                           qmtl(l,itype)=qmtl(l,itype)-                 &
     &                       AIMAG(g_coef(n)*(                          &
     &                               uu+                                &
     &                               dd*ddn(l,itype)+                   &
     &                               locorbuu+                          &
     &                               locorbdd*ddn(l,itype)+             &
     &                               locorbulo*uulon(lo1,itype)+        &
     &                               locorblou*uulon(lo1,itype)+        &
     &                               locorbdlo*dulon(lo1,itype)+        &
     &                               locorblod*dulon(lo1,itype)+        &
     &                               locorblolo*uloulopn(lo1,lo1,itype) &
     &                                   ))                             &
     &                       /atoms%neq(itype)*sfp                      
                                                                        
                         ELSE 
                           qmtl(l,itype)=qmtl(l,itype)-                 &
     &                          AIMAG(g_coef(n)*(uu+dd*ddn(l,itype)))   &
     &                          /atoms%neq(itype)*sfp                   
                                                                        
                         ENDIF 
                      ENDIF 
                        !>                                              
                     ENDDO 
                     !>                                                 
                        !mp-loop                                        
                  ENDDO 
                     !m-loop                                            
               ENDDO 
                 !lp-loop                                               
            ENDDO 
              !l-loop                                                   
         ENDDO 
         nt=nt+atoms%neq(itype) 
         DEALLOCATE(vecmat) 
         IF(atoms%nlo(itype)>=1)THEN 
            DEALLOCATE(leftlomat) 
            DEALLOCATE(rightlomat) 
         ENDIF 
           !itype                                                       
      ENDDO 
      DEALLOCATE (s1,s2,s3,s4,s5,s6,g_coef,lh_sp) 
                                                                        
      RETURN 
      END SUBROUTINE 
                                                                        
      END                                           
