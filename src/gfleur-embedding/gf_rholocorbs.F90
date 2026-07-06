!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_rholocorbs 
      use m_juDFT
          IMPLICIT NONE
!***************************************************************        
!     Make preparations for the calculation of the local-orbital        
!     contributions to the charge.                                      
!     Frank Freimuth, April 2007                                        
!***************************************************************        
      CONTAINS 
                                     !mpi,                              
      SUBROUTINE gf_rholocorbs(atoms,                                   &
     &           itype,enpara,vr,jspin,f,g,                             &
     &           us,dus,uds,duds,ddn,lapw,                              &
     &           cell,sym,bk,                                           &
     &           flo,alocof,blocof,clocof,basindex,                     &
     &           uulon,dulon,uloulopn)                                  
                                                                        
      USE m_radflo 
      USE m_ylm 
      USE m_constants ,ONLY: pimach 
      USE m_gf_types 
                                                                        
      IMPLICIT NONE 
      TYPE(t_cell),INTENT(IN)  :: cell 
      TYPE(t_atoms),INTENT(IN) :: atoms 
      TYPE(t_enpara),INTENT(IN):: enpara 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_sym),INTENT(IN)   :: sym 
      INTEGER,INTENT(IN)       :: jspin,itype 
      REAL, INTENT(IN)         :: vr(:,:,:) 
!      TYPE(t_mpi),INTENT(IN)      :: mpi                               
      REAL, INTENT(IN)         :: f(:,:,0:),g(:,:,0:) 
      REAL, INTENT(OUT)        :: flo(:,:,:) 
      INTEGER,INTENT(OUT)      :: basindex(:,:) 
                                                 ! Kpt                  
      REAL, INTENT(IN)         :: Bk(3) 
      COMPLEX,INTENT(OUT)      :: alocof(:,0:,:),blocof(:,0:,:) 
      COMPLEX,INTENT(OUT)      :: clocof(:,0:,:) 
      REAL,INTENT(IN)  :: us(0:,:),dus(0:,:),uds(0:,:),duds(0:,:) 
      REAL,INTENT(IN)  :: ddn(0:,:) 
      REAL,INTENT(OUT) :: uulon(atoms%nlod,atoms%ntype) 
      REAL,INTENT(OUT) :: dulon(atoms%nlod,atoms%ntype) 
      REAL,INTENT(OUT) :: uloulopn(atoms%nlod,atoms%nlod,atoms%ntype) 
                                                                        
      INTEGER jmtd,lmaxd,nlod,ntype,lo,l,bas,nalo,nat 
      INTEGER loc,naup,ind,natom,locnum,invloop,k 
      INTEGER nap,locvec,lm 
      REAL,ALLOCATABLE::ulos(:,:) 
      REAL,ALLOCATABLE::dulos(:,:) 
      REAL,ALLOCATABLE::uuilon(:,:),duilon(:,:) 
      REAL,ALLOCATABLE::ulouilopn(:,:,:) 
      REAL,ALLOCATABLE::alo1(:),blo1(:),clo1(:) 
      REAL ws,ka,kb,arg,tpi,const 
      COMPLEX phase 
      REAL v(3),fkr(3),fkp(3) 
      COMPLEX :: ylm(16),ci 
                                                                        
      lmaxd=maxval(atoms%lmax)
      jmtd=maxval(atoms%jri) 
      nlod=atoms%nlod 
      ntype=atoms%ntype 
      tpi=2*pimach() 
      const=4*pimach()/sqrt(cell%omtil) 
      ci=cmplx(0.0,1.0) 
                                                                        
      ALLOCATE(ulos(nlod,ntype)) 
      ALLOCATE(dulos(nlod,ntype)) 
      ALLOCATE(uuilon(nlod,ntype),duilon(nlod,ntype)) 
      ALLOCATE(ulouilopn(nlod,nlod,ntype)) 
                                                                        
      CALL radflo(                                                      &
     &        atoms%ntype,atoms%nlod,1,                                 &
     &        jmtd,                                                     &
     &        lmaxd,itype,jspin,                                        &
     &        enpara%ello(1:atoms%nlod,1:atoms%ntype,jspin),            &
     &        vr(1:jmtd,itype,jspin),atoms%jri(itype),                  &
     &        atoms%rmsh(1,itype),atoms%dx(itype),                      &
     &        f,                                                        &
     &        g,                                                        &
     &        atoms%llo,                                                &
     &        atoms%nlo,                                                &
     &        atoms%l_dulo(1:atoms%nlod,itype),                         &
     &        1,                                                        &
     &        atoms%ulo_der,                                            &
     &        ulos,dulos,uulon,dulon,uloulopn,                          &
     &        uuilon,duilon,ulouilopn,flo)                              
                                       !jspins,                         
                          !mpi%irank,                                   
                                                                        
      DEALLOCATE(uuilon,duilon,ulouilopn) 
                                                                        
      ALLOCATE(alo1(atoms%nlod)) 
      ALLOCATE(blo1(atoms%nlod)) 
      ALLOCATE(clo1(atoms%nlod)) 
                                                                        
      DO lo=1,atoms%nlo(itype) 
             l=atoms%llo(lo,itype) 
             ws = uds(l,itype)*dus(l,itype) - us(l,itype)*duds(l,itype) 
             ka = 1.0/ws* (duds(l,itype)*ulos(lo,itype)-                &
     &            uds(l,itype)*dulos(lo,itype))                         
             kb = 1.0/ws* (us(l,itype)*dulos(lo,itype)-                 &
     &            dus(l,itype)*ulos(lo,itype))                          
             clo1(lo) = 1.0/sqrt(ka**2+ (kb**2)*ddn(l,itype)+1.0+       &
     &                   2.0*ka*uulon(lo,itype)+2.0*kb*dulon(lo,itype)) 
             alo1(lo) = ka*clo1(lo) 
             blo1(lo) = kb*clo1(lo) 
      ENDDO 
                                                                        
      DEALLOCATE(ulos,dulos) 
                                                                        
      bas=1 
      nalo=1 
      DO ind=1,itype 
         naup=nalo+atoms%neq(ind)-1 
         DO natom=nalo,naup 
            IF(atoms%invsat(natom)==2)CYCLE 
            DO loc=1,atoms%nlo(ind) 
               basindex(natom,loc)=bas 
               locnum=(atoms%invsat(natom)+1)*(2*atoms%llo(loc,ind)+1) 
               bas=bas+locnum 
                  !loc                                                  
            ENDDO 
               !natom                                                   
         ENDDO 
         nalo=naup+1 
            !ind                                                        
      ENDDO 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
      nalo=naup-atoms%neq(itype)+1 
      DO natom=nalo,naup 
        IF(atoms%invsat(natom)==2)CYCLE 
        DO loc=1,atoms%nlo(itype) 
          l=atoms%llo(loc,itype) 
          bas=basindex(natom,loc) 
          locnum=2*atoms%llo(loc,itype)+1 
          DO invloop=0,atoms%invsat(natom) 
            DO ind=1,locnum 
              locvec=lapw%kveclo(bas) 
              v(1)=bk(1)+lapw%k%k1(locvec,jspin) 
              v(2)=bk(2)+lapw%k%k2(locvec,jspin) 
              v(3)=bk(3)+lapw%k%k3(locvec,jspin) 
              arg=    v(1)*atoms%taual(1,natom) 
              arg=arg+v(2)*atoms%taual(2,natom) 
              arg=arg+v(3)*atoms%taual(3,natom) 
              arg=arg*tpi 
              phase=cmplx(cos(arg),sin(arg)) 
              phase=phase*const*(atoms%rmt(itype))**2/2 
              phase=phase*ci**l 
              IF(invloop==1) CALL juDFT_error("not yet implemented",calledby="gf_rholocorbs.F90")
              nap=sym%invtab(atoms%ngopr(natom)) 
              fkr=MATMUL(v,real(sym%mrot(:,:,nap))) 
              fkp=MATMUL(fkr,cell%bmat) 
              CALL ylm4(3,fkp,ylm) 
              lm=l*(l+1) 
              alocof(bas,lm-l:lm+l,natom-nalo+1)=                       &
     &           alo1(loc)*phase*conjg(ylm(lm-l+1:lm+l+1))              
              blocof(bas,lm-l:lm+l,natom-nalo+1)=                       &
     &           blo1(loc)*phase*conjg(ylm(lm-l+1:lm+l+1))              
              clocof(bas,lm-l:lm+l,natom-nalo+1)=                       &
     &           clo1(loc)*phase*conjg(ylm(lm-l+1:lm+l+1))              
              bas=bas+1 
            ENDDO 
          ENDDO 
        ENDDO 
      ENDDO 
                                                                        
      DEALLOCATE(alo1,blo1,clo1) 
      END SUBROUTINE gf_rholocorbs 
      END                                           
