!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_od_hsvac
  USE m_juDFT
CONTAINS
  SUBROUTINE od_hsvac(&
       vacuum,stars,DIMENSION, oneD,atoms, jsp,input,vxy,vz,evac,cell,&
       bkpt,lapw, MM,vM,m_cyl,n2d_1, n_size,n_rank,sym,noco,nv2,l_real,hamOvlp)

    !     subroutine for calculating the hamiltonian and overlap matrices in
    !     the vacuum in the case of 1-dimensional calculations
    !     Y. Mokrousov June 2002             

    USE m_cylbes
    USE m_dcylbs
    USE m_od_vacfun
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN)  :: DIMENSION
    TYPE(t_oneD),INTENT(IN)       :: oneD
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_vacuum),INTENT(IN)     :: vacuum
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_hamOvlp),INTENT(INOUT) :: hamOvlp
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: vM
    INTEGER, INTENT (IN) :: MM 
    INTEGER, INTENT (IN) :: jsp ,n_size,n_rank,n2d_1 
    INTEGER, INTENT (IN) :: m_cyl
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (INOUT) :: vxy(vacuum%nmzxyd,n2d_1-1,2)
    INTEGER, INTENT (OUT):: nv2(input%jspins)
    REAL,    INTENT (INOUT) :: vz(vacuum%nmz,2,4)
    REAL,    INTENT (IN) :: evac(2,input%jspins)
    REAL,    INTENT (IN) :: bkpt(3)

    LOGICAL, INTENT(IN)  :: l_real

    !     ..
    !     .. Local Scalars ..
    COMPLEX hij,sij,apw_lo,exp1,exp2,exp3,am,bm,ic
    REAL    d2,wronk,gr,gphi,qq,x,y
    INTEGER i,i2,ii,ik,j,jk,k,jspin,ipot,npot,ii0 ,l,i3,imz,m
    INTEGER jspin1,jspin2,jmax,irec2,irec3,ivac,ind1,gi
    INTEGER i_start,nc,nc_0,rotax,chiral,zi,m1,z,indm,indl
    !     ..
    !     .. Local Arrays ..

    INTEGER, ALLOCATABLE :: nvp(:,:),ind(:,:,:)
    INTEGER, ALLOCATABLE :: kvac3(:,:),map1(:,:)
    COMPLEX, ALLOCATABLE :: tddv(:,:,:,:)
    COMPLEX, ALLOCATABLE :: tduv(:,:,:,:)
    COMPLEX, ALLOCATABLE :: tudv(:,:,:,:)
    COMPLEX, ALLOCATABLE :: tuuv(:,:,:,:)
    COMPLEX, ALLOCATABLE ::  a(:,:,:),b(:,:,:)
    COMPLEX, ALLOCATABLE :: ai(:,:,:),bi(:,:,:)
    REAL, ALLOCATABLE :: bess(:),dbss(:),bess1(:)
    REAL, ALLOCATABLE :: ddnv(:,:,:),dudz(:,:,:)
    REAL, ALLOCATABLE :: duz(:,:,:)
    REAL, ALLOCATABLE :: udz(:,:,:),uz(:,:,:)
    !     ..
    ic  = CMPLX(0.,1.)
    d2 = SQRT(cell%omtil/cell%area)

 
    ALLOCATE (&
         ai(-vM:vM,DIMENSION%nv2d,DIMENSION%nvd),bi(-vM:vM,DIMENSION%nv2d,DIMENSION%nvd),&
         nvp(DIMENSION%nv2d,input%jspins),ind(stars%ng2,DIMENSION%nv2d,input%jspins),&
         kvac3(DIMENSION%nv2d,input%jspins),map1(DIMENSION%nvd,input%jspins),&
         tddv(-vM:vM,-vM:vM,DIMENSION%nv2d,DIMENSION%nv2d),&
         tduv(-vM:vM,-vM:vM,DIMENSION%nv2d,DIMENSION%nv2d),&
         tudv(-vM:vM,-vM:vM,DIMENSION%nv2d,DIMENSION%nv2d),&
         tuuv(-vM:vM,-vM:vM,DIMENSION%nv2d,DIMENSION%nv2d),&
         a(-vM:vM,DIMENSION%nvd,input%jspins),b(-vM:vM,DIMENSION%nvd,input%jspins),&
         bess(-vM:vM),dbss(-vM:vM),bess1(-vM:vM),&
         ddnv(-vM:vM,DIMENSION%nv2d,input%jspins),dudz(-vM:vM,DIMENSION%nv2d,input%jspins),&
         duz(-vM:vM,DIMENSION%nv2d,input%jspins),&
         udz(-vM:vM,DIMENSION%nv2d,input%jspins),uz(-vM:vM,DIMENSION%nv2d,input%jspins) )

    !--->     set up mapping function from 3d-->1d lapws
    !--->            creating arrays ind and nvp

    DO jspin = 1,input%jspins

       nv2(jspin) = 0
       k_loop:DO  k = 1,lapw%nv(jspin)
          DO  j = 1,nv2(jspin)
             IF (lapw%k3(k,jspin).EQ.kvac3(j,jspin)) THEN
                map1(k,jspin) = j
                CYCLE k_loop
             END IF
          ENDDO
          nv2(jspin) = nv2(jspin) + 1
          IF (nv2(jspin)>DIMENSION%nv2d)  CALL juDFT_error("dimension%nv2d",calledby ="od_hsvac")
          kvac3(nv2(jspin),jspin) = lapw%k3(k,jspin)
          map1(k,jspin) = nv2(jspin)
       END DO k_loop

       DO ik = 1,DIMENSION%nv2d
          nvp(ik,jspin) = 0
          DO i = 1,stars%ng2
             ind(i,ik,jspin) = 0
          END DO
       END DO

       DO k = 1,lapw%nv(jspin)
          ik = map1(k,jspin)
          nvp(ik,jspin) = nvp(ik,jspin) + 1
          ind(nvp(ik,jspin),ik,jspin) = k
       END DO

    ENDDO

    npot = 1      
    ivac = 1

    IF (noco%l_noco) THEN
       !--->         load the non-warping part of the potential
       READ (25)((vz(imz,ivac,ipot),imz=1,vacuum%nmz),ipot=1,4)
       npot = 3
    ENDIF

    DO ipot = 1,npot

       IF (noco%l_noco) THEN
          READ (25)((vxy(imz,k,ivac), imz=1,vacuum%nmzxy),k=1,n2d_1-1)
          !--->  l_J we want to average the diagonal elements of the pot. matrix
       ENDIF ! loco

       !     get the wavefunctions and set up the tuuv, etc matrices

       CALL od_vacfun(&
            m_cyl,cell,vacuum,DIMENSION,stars,&
            jsp,input,noco,ipot,oneD,n2d_1, ivac,evac(1,1),bkpt,MM,vM,&
            vxy(1,1,ivac),vz,kvac3,nv2, tuuv,tddv,tudv,tduv,uz,duz,udz,dudz,ddnv)

       IF (noco%l_noco) THEN

          DO jspin = 1,input%jspins

             DO k = 1,lapw%nv(jspin)
                irec3 = stars%ig(lapw%k1(k,jspin),lapw%k2(k,jspin),lapw%k3(k,jspin))
                IF (irec3.NE.0) THEN
                   irec2 = stars%ig2(irec3)
                   gr = stars%sk2(irec2)
                   gphi = stars%phi2(irec2)
                   i2 = map1(k,jspin)
                   qq = gr*cell%z1
                   CALL cylbes(vM,qq,bess)
                   CALL dcylbs(vM,qq,bess,dbss)
                   DO m = -vM,vM
                      wronk = uz(m,i2,jspin)*dudz(m,i2,jspin) - udz(m,i2,jspin)*duz(m,i2,jspin)
                      a(m,k,jspin)=EXP(-CMPLX(0.0,m*gphi))*(ic**m)*&
                           CMPLX(dudz(m,i2,jspin)*bess(m)- udz(m,i2,jspin)*gr*dbss(m),0.0) /(d2*wronk)
                      b(m,k,jspin)=EXP(-CMPLX(0.0,m*gphi))*(ic**m)*&
                           CMPLX(-duz(m,i2,jspin)*bess(m)+ uz(m,i2,jspin)*gr*dbss(m),0.0) /(d2*wronk)
                   END DO
                END IF
             ENDDO

          ENDDO  ! jspin

       ELSE 

          DO k = 1,lapw%nv(jsp)
             irec3 = stars%ig(lapw%k1(k,jsp),lapw%k2(k,jsp),lapw%k3(k,jsp))
             IF (irec3.NE.0) THEN
                irec2 = stars%ig2(irec3)
                gr = stars%sk2(irec2)
                gphi = stars%phi2(irec2)
                i2 = map1(k,jsp)
                qq = gr*cell%z1
                CALL cylbes(vM,qq,bess) 
                CALL dcylbs(vM,qq,bess,dbss)
                DO m = -vM,vM
                   wronk = uz(m,i2,jsp)*dudz(m,i2,jsp) - udz(m,i2,jsp)*duz(m,i2,jsp) 
                   a(m,k,1)=EXP(-CMPLX(0.0,m*gphi))*(ic**m)*&
                        CMPLX(dudz(m,i2,jsp)*bess(m)- udz(m,i2,jsp)*gr*dbss(m),0.0) /(d2*wronk)

                   b(m,k,1)=EXP(-CMPLX(0.0,m*gphi))*(ic**m)*&
                        CMPLX(-duz(m,i2,jsp)*bess(m)+ uz(m,i2,jsp)*gr*dbss(m),0.0) /(d2*wronk)

                END DO
             END IF
          ENDDO

       ENDIF ! loco
       !     update hamiltonian and overlap matrices

       IF (ipot.EQ.1 .OR. ipot.EQ.2) THEN
          jspin = ipot
          !+gb||
          IF (ipot.EQ.1) THEN
             nc = 0
             i_start = n_rank
          ELSE
             nc = nc + atoms%nlotot
             nc_0 = nc
             i_start = MOD(MOD(n_rank - (lapw%nv(1)+atoms%nlotot),n_size) + n_size,n_size)
          ENDIF

          DO  i = i_start+1,lapw%nv(jspin),n_size
             ik = map1(i,jspin)
             nc = nc + 1
             IF (ipot.EQ.1) THEN
                jspin = 1
                ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
             ELSEIF (ipot.EQ.2) THEN
                jspin = 2
                ii0=nc*(nc-1)/2*n_size-(nc-1)*(n_size-n_rank-1)+&
                     lapw%nv(1)+atoms%nlotot
             ENDIF
             jspin1 = jsp
             IF (noco%l_noco) jspin1 = jspin
             DO j = 1,i - 1
                ii = ii0 + j
                !     overlap: only  (g-g') parallel=0        
                IF (map1(j,jspin).EQ.ik) THEN
                   sij = (0.0,0.0)
                   DO m = -vM,vM
                      sij = sij + CONJG(a(m,i,jspin))*a(m,j,jspin) &
                           +CONJG(b(m,i,jspin))*b(m,j,jspin) *ddnv(m,ik,jspin1)
                   END DO
                   IF (l_real) THEN 
                      hamOvlp%b_r(ii) = hamOvlp%b_r(ii) + REAL(sij)
                   ELSE
                      hamOvlp%b_c(ii) = hamOvlp%b_c(ii) + sij
                   ENDIF
                END IF
             ENDDO
             ii = ii0 + i
             sij = (0.0,0.0)
             DO m = -vM,vM
                sij = sij + CONJG(a(m,i,jspin))*a(m,i,jspin)+ &
                     CONJG(b(m,i,jspin))*b(m,i,jspin)*ddnv(m,ik,jspin1)
             END DO

             IF (l_real) THEN
                hamOvlp%b_r(ii) = hamOvlp%b_r(ii) + REAL(sij)
             ELSE 
                hamOvlp%b_c(ii) = hamOvlp%b_c(ii) + sij
             ENDIF
          ENDDO
       ENDIF ! ipot.eq.1.or.2
       !   hamiltonian update 
       !   for the noncylindr. contributions we use the cutoff of m_cyl        
       IF (ipot.EQ.1) THEN
          jspin1 = 1
          jspin2 = 1
          nc = 0
          i_start = n_rank
       ELSEIF (ipot.EQ.2) THEN
          jspin1 = 2
          jspin2 = 2
          nc = nc_0
          i_start = MOD(MOD(n_rank - (lapw%nv(1)+atoms%nlotot),n_size) + n_size,n_size)
       ELSEIF (ipot.EQ.3) THEN
          jspin1 = 2
          jspin2 = 1
          nc = nc_0
          i_start = MOD(MOD(n_rank - (lapw%nv(1)+atoms%nlotot),n_size) + n_size,n_size)
       ENDIF

       ai(:,:,:) = CMPLX(0.,0.)
       bi(:,:,:) = CMPLX(0.,0.)

       DO ik = 1,nv2(jspin1)
          DO jk = 1,nv2(jspin2)
             i3 = kvac3(ik,jspin1) - kvac3(jk,jspin2) 
             DO l = -vM,vM
                DO m = -vM,vM
                   IF (l.EQ.m .OR. (iabs(m).LE.m_cyl .AND. iabs(l).LE.m_cyl)) THEN
                      ind1 = oneD%ig1(i3,m-l)
                      IF (ind1.NE.0) THEN
                         DO gi = 1,nvp(ik,jspin1)
                            i = ind(gi,ik,jspin1)
                            ai(l,jk,i) = ai(l,jk,i) +&
                                 CONJG(a(m,i,jspin1))*tuuv(m,l,ik,jk) + CONJG(b(m,i,jspin1))*tduv(m,l,ik,jk)
                            bi(l,jk,i) = bi(l,jk,i) +&
                                 CONJG(a(m,i,jspin1))*tudv(m,l,ik,jk) + CONJG(b(m,i,jspin1))*tddv(m,l,ik,jk)

                         END DO
                      END IF
                   END IF   ! noncyl. contributions
                END DO
             END DO
          END DO
       END DO

       DO i = i_start+1, lapw%nv(jspin1), n_size
          ik = map1(i,jspin1)
          nc = nc + 1
          IF (ipot.EQ.1) THEN
             ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
             jmax = i
          ELSEIF (ipot.EQ.2) THEN
             ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1) + lapw%nv(1)+atoms%nlotot
             jmax = i
          ELSEIF (ipot.EQ.3) THEN
             ii0 = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
             jmax = lapw%nv(jspin2)
          ENDIF
          DO j = 1,jmax
             ii = ii0 + j
             jk = map1(j,jspin2)
             hij = CMPLX(0.,0.)
             DO l = -vM,vM
                hij = hij + ai(l,jk,i)*a(l,j,jspin2) + bi(l,jk,i)*b(l,j,jspin2)
             END DO
             IF (l_real) THEN
                hamOvlp%a_r(ii) = hamOvlp%a_r(ii) + REAL(hij)
             ELSE 
                hamOvlp%a_c(ii) = hamOvlp%a_c(ii) + hij
             ENDIF
          END DO
       END DO

    ENDDO !ipot

    DEALLOCATE (ai,bi,nvp,ind,kvac3,map1, tddv,tduv,tudv,tuuv,a,b,bess,dbss,bess1, ddnv,dudz,duz,udz,uz )

    RETURN
  END SUBROUTINE od_hsvac
END MODULE m_od_hsvac
