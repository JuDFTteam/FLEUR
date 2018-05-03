!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_mt_tofrom_grid
  USE m_types
  PRIVATE
  REAL,PARAMETER    :: d_15 = 1.e-15
  INTEGER,PARAMETER :: ndvgrd=6 ! this should be consistent across GGA derivative routines
  REAL, ALLOCATABLE :: ylh(:,:,:),ylht(:,:,:),ylhtt(:,:,:)
  REAL, ALLOCATABLE :: ylhf(:,:,:),ylhff(:,:,:),ylhtf(:,:,:)
  REAL, ALLOCATABLE :: wt(:),rx(:,:),thet(:)
  PUBLIC :: init_mt_grid,mt_to_grid,mt_from_grid,finish_mt_grid
CONTAINS
  SUBROUTINE init_mt_grid(nsp,jspins,atoms,sphhar,xcpot,sym,l_grad)
    USE m_gaussp
    USE m_lhglptg
    USE m_lhglpts
    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: nsp,jspins
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    CLASS(t_xcpot),INTENT(IN)   :: xcpot
    TYPE(t_sym),INTENT(IN)      :: sym
    LOGICAL,INTENT(IN)          :: l_grad
    
    ! generate nspd points on a sherical shell with radius 1.0
    ! angular mesh equidistant in phi,
    ! theta are zeros of the legendre polynomials
    ALLOCATE(wt(nsp),rx(3,nsp),thet(nsp))
    CALL gaussp(atoms%lmaxd, rx,wt)
    ! generate the lattice harmonics on the angular mesh
    ALLOCATE ( ylh(nsp,0:sphhar%nlhd,sphhar%ntypsd))
    IF (l_grad) ALLOCATE(ylht,ylhtt,ylhf,ylhff,ylhtf,MOLD=ylh )

    IF (l_grad) THEN
       CALL lhglptg(sphhar,atoms,rx,nsp,xcpot,sym,&
            ylh,thet,ylht,ylhtt,ylhf,ylhff,ylhtf)
    ELSE
       CALL lhglpts( sphhar,atoms, rx,nsp, sym, ylh)
    END IF
  END SUBROUTINE init_mt_grid
  
  SUBROUTINE mt_to_grid(atoms,sphhar,den,nsp,jspins,n,l_grad,ch,grad)
    USE m_grdchlh
    USE m_mkgylm
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_potden),INTENT(IN)   :: den
    INTEGER,INTENT(IN)          :: n,jspins,nsp
    LOGICAL,INTENT(IN)          :: l_grad
    REAL,INTENT(OUT)               :: ch(:,:)
    TYPE(t_gradients),INTENT(INOUT)::grad
    
    REAL, ALLOCATABLE :: chlh(:,:,:),chlhdr(:,:,:),chlhdrr(:,:,:)
    REAL, ALLOCATABLE :: chdr(:,:),chdt(:,:),chdf(:,:)
    REAL, ALLOCATABLE :: chdrr(:,:),chdtt(:,:),chdff(:,:),chdtf(:,:)
    REAL, ALLOCATABLE :: chdrt(:,:),chdrf(:,:)
    INTEGER:: nd,lh,js,jr,kt,k

    nd = atoms%ntypsy(SUM(atoms%neq(:n-1))+1)

    ALLOCATE ( chlh(atoms%jmtd,0:sphhar%nlhd,jspins))

    IF (l_grad)  THEN
       ALLOCATE(chdr(nsp,jspins),chdt(nsp,jspins),chdf(nsp,jspins),chdrr(nsp,jspins),&
            chdtt(nsp,jspins),chdff(nsp,jspins),chdtf(nsp,jspins),chdrt(nsp,jspins),&
            chdrf(nsp,jspins) )
       ALLOCATE (chlhdr(atoms%jmtd,0:sphhar%nlhd,jspins))
       ALLOCATE (chlhdrr(atoms%jmtd,0:sphhar%nlhd,jspins))
    ENDIF
    
    DO lh = 0,sphhar%nlh(nd)

       !         calculates gradients of radial charge densities of l=> 0.
       !         rho*ylh/r**2 is charge density. chlh=rho/r**2.
       !         charge density=sum(chlh*ylh).
       !         chlhdr=d(chlh)/dr, chlhdrr=dd(chlh)/drr.
       
       DO js = 1,jspins      
          DO jr = 1,atoms%jri(n)
             chlh(jr,lh,js) = den%mt(jr,lh,n,js)/(atoms%rmsh(jr,n)*atoms%rmsh(jr,n))
          ENDDO
          IF (l_grad) CALL grdchlh(1,1,atoms%jri(n),atoms%dx(n),atoms%rmsh(1,n),&
               chlh(1,lh,js),ndvgrd, chlhdr(1,lh,js),chlhdrr(1,lh,js))
          
       ENDDO ! js
    ENDDO   ! lh
    ch(:,:)    = 0.0     ! charge density (on extended grid for all jr)
    kt=0
    DO jr = 1,atoms%jri(n)
       !         following are at points on jr-th sphere.
       !  generate the densities on an angular mesh
       DO js = 1,jspins
          DO lh = 0,sphhar%nlh(nd)
             DO k = 1,nsp
                ch(kt+k,js) = ch(kt+k,js) + ylh(k,lh,nd)*chlh(jr,lh,js)
             ENDDO
          ENDDO
       ENDDO
       IF (l_grad) THEN
          chdr(:,:)  = 0.0     ! d(ch)/dr
          chdt(:,:)  = 0.0     ! d(ch)/dtheta
          chdf(:,:)  = 0.0     ! d(ch)/dfai
          chdrr(:,:) = 0.0     ! dd(ch)/drr
          chdtt(:,:) = 0.0     ! dd(ch)/dtt
          chdff(:,:) = 0.0     ! dd(ch)/dff
          chdtf(:,:) = 0.0     ! dd(ch)/dtf
          chdrt(:,:) = 0.0     ! d(d(ch)/dr)dt
          chdrf(:,:) = 0.0     ! d(d(ch)/dr)df
          !  generate the derivatives on an angular mesh
          DO js = 1,jspins
             DO lh = 0,sphhar%nlh(nd)
                ! 
                DO k = 1,nsp
                   chdr(k,js) =chdr(k,js)+ ylh(k,lh,nd)*chlhdr(jr,lh,js)
                   chdrr(k,js)=chdrr(k,js)+ylh(k,lh,nd)*chlhdrr(jr,lh,js)
                ENDDO
                
                DO k = 1,nsp
                   chdrt(k,js)=chdrt(k,js)+ylht(k,lh,nd)*chlhdr(jr,lh,js)
                   chdrf(k,js)=chdrf(k,js)+ylhf(k,lh,nd)*chlhdr(jr,lh,js)
                   chdt(k,js) =chdt(k,js) +ylht(k,lh,nd)*chlh(jr,lh,js)
                   chdf(k,js) =chdf(k,js) +ylhf(k,lh,nd)*chlh(jr,lh,js)
                   chdtt(k,js)=chdtt(k,js)+ylhtt(k,lh,nd)*chlh(jr,lh,js)
                   chdff(k,js)=chdff(k,js)+ylhff(k,lh,nd)*chlh(jr,lh,js)
                   chdtf(k,js)=chdtf(k,js)+ylhtf(k,lh,nd)*chlh(jr,lh,js)
                ENDDO
             ENDDO ! lh
          ENDDO   ! js
          
       
          CALL mkgylm(jspins,atoms%rmsh(jr,n),thet,nsp,&
               ch(kt+1:,:),chdr,chdt,chdf,chdrr,chdtt,chdff,chdtf,chdrt,chdrf,grad,kt)
       ENDIF
       kt=kt+nsp
    END DO
    !Set charge to minimum value
    ch(:,:)=MAX(ch(:,:),d_15)
    
  END SUBROUTINE mt_to_grid


  SUBROUTINE mt_from_grid(atoms,sphhar,nsp,n,jspins,v_in,vr)
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN) :: atoms
    TYPE(t_sphhar),INTENT(IN):: sphhar
    INTEGER,INTENT(IN)       :: nsp,jspins,n
    REAL,INTENT(IN)          :: v_in(:,:)
    REAL,INTENT(INOUT)       :: vr(:,0:,:,:)
    
    REAL    :: vpot(nsp),vlh
    INTEGER :: js,kt,lh,jr,nd
    nd=atoms%ntypsy(SUM(atoms%neq(:n-1))+1)
    DO js = 1,jspins
       !
       kt=0
       DO jr=1,atoms%jri(n)
          vpot=v_in(kt+1:kt+nsp,js)*wt(:)!  multiplicate v_in with the weights of the k-points
          
          DO lh = 0,sphhar%nlh(nd)
             !
             ! --->        determine the corresponding potential number
             !c            through gauss integration
             !
             vlh=dot_PRODUCT(vpot(:),ylh(:nsp,lh,nd))
             vr(jr,lh,n,js) = vr(jr,lh,n,js) + vlh
          ENDDO ! lh
          kt=kt+nsp
       ENDDO   ! jr
    ENDDO

  END SUBROUTINE mt_from_grid
    
  SUBROUTINE finish_mt_grid()
    DEALLOCATE(ylh,wt,rx,thet)
    IF (ALLOCATED(ylht)) DEALLOCATE(ylht,ylhtt,ylhf,ylhff,ylhtf)
  END SUBROUTINE finish_mt_grid

END MODULE m_mt_tofrom_grid
