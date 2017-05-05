!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt_overlap
  use m_juDFT
  implicit none
CONTAINS
  SUBROUTINE hsmt_overlap(input,atoms,n_size,n_rank,isp,l_socfirst,hlpmsize,ab_dim,&
       noco,cell,nintsp, lapw,usdus,gk,fj,gj,bb)
!Calculate overlap matrix
#include"cpp_double.h"
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_hsmt_spinor
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(INOUT)  :: lapw!lpaw%nv_tot is updated
    TYPE(t_usdus),INTENT(INOUT) :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: isp
    INTEGER, INTENT (IN) :: n_size,n_rank,ab_dim
    INTEGER, INTENT (IN) :: hlpmsize,nintsp
    LOGICAL, INTENT (IN) :: l_socfirst
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
#ifdef CPP_INVERSION
    REAL,    INTENT (INOUT) :: bb(:)!(matsize)
#else
    COMPLEX, INTENT (INOUT) :: bb(:)
#endif
    !     ..
    !     .. Local Scalars ..
    REAL tnn(3),fct,ski(3),fjkiln,gjkiln
    COMPLEX chi11,chi21,chi22
    INTEGER ii,iii,ij,k,ki,kj,l,n,n0,n1,nn,kjmax, nsp,iintsp,jintsp
    INTEGER nc ,kii

    !     ..
    !     .. Local Arrays ..
    COMPLEX chi(2,2),bbwa(maxval(lapw%nv))
    REAL qssbti(3),qssbtj(3)
    REAL fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)     
    REAL, ALLOCATABLE :: plegend(:,:)
    REAL, ALLOCATABLE :: rph(:,:),cph(:,:)
    COMPLEX, ALLOCATABLE :: bbhlp(:)


    IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
       !---> pk non-collinear
       !--->    initialize help array
       ALLOCATE ( bbhlp(hlpmsize) )
       bbhlp=cmplx(0.0,0.0)
    ENDIF
    !     ..
    DO l = 0,atoms%lmaxd
       fleg1(l) = real(l+l+1)/real(l+1)
       fleg2(l) = real(l)/real(l+1)
       fl2p1(l) = real(l+l+1)/fpi_const
    END DO
    !
    
    !$OMP PARALLEL  DEFAULT(shared)&
    !$OMP PRIVATE(kii,ki,nc,iii,kjmax,ski,kj,plegend,l,n1,n)&
    !$OMP PRIVATE(n0,rph,cph,nn,tnn,fjkiln,gjkiln)&
    !$OMP PRIVATE(fct,ij)&
    !$OMP PRIVATE(chi,chi11,chi21,chi22,nsp)&
    !$OMP PRIVATE(bbwa,ii) 
    ALLOCATE(rph(maxval(lapw%nv),ab_dim))
    ALLOCATE(cph(maxval(lapw%nv),ab_dim))
    ALLOCATE(plegend(maxval(lapw%nv),0:atoms%lmaxd))
   
    plegend=0.0
    DO iintsp = 1,nintsp
       IF (iintsp.EQ.1) THEN
          qssbti = - noco%qss/2
       ELSE
          qssbti = + noco%qss/2
       ENDIF
       DO jintsp = 1,iintsp
          IF (jintsp.EQ.1) THEN
             qssbtj = - noco%qss/2
          ELSE
             qssbtj = + noco%qss/2
          ENDIF

          nc = 0                                    ! maybe IF (iintsp.EQ
          IF ( noco%l_noco .AND. (n_size.GT.1) ) THEN
             !--->       for EV-parallelization & noco ( see comments at top )
             lapw%nv_tot = lapw%nv(1) + lapw%nv(2)
             IF (noco%l_ss)  CALL juDFT_error("ev-parallelization & spin-spiral !",calledby ="hsmt_overlap")
          ELSE
             lapw%nv_tot = lapw%nv(iintsp)
          ENDIF

          !
          !$OMP  DO SCHEDULE(DYNAMIC,1)
          DO  kii =  n_rank, lapw%nv_tot-1, n_size
             ki = mod(kii,lapw%nv(iintsp)) + 1
             nc = nc + 1
             !$ nc= 1+ (kii-n_rank)/n_size
             iii = nc*(nc-1)/2 * n_size - (nc-1)*(n_size - n_rank - 1)
             !  ii = nc*(nc+1)/2 * n_size -  nc   *(n_size - n_rank - 1)

             IF (noco%l_ss.OR.noco%l_constr.OR.l_socfirst) THEN
                kjmax = lapw%nv(jintsp)
             ELSE
                kjmax = ki
             ENDIF
             ski = (/lapw%k1(ki,iintsp),lapw%k2(ki,iintsp),lapw%k3(ki,iintsp)/) + qssbti
             !--->       legendre polynomials
             plegend(:,0)=1.0
             DO kj = 1,kjmax
                plegend(kj,1) = dot_product(gk(kj,:,jintsp),gk(ki,:,iintsp))
             END DO
             DO l = 1,maxval(atoms%lmax) - 1
                plegend(:,l+1) = fleg1(l)*plegend(:,1)*plegend(:,l) - fleg2(l)*plegend(:,l-1)
             END DO
             !include factor fl2p1 already here
             DO l=0,maxval(atoms%lmax)
                plegend(:kjmax,l)=plegend(:kjmax,l)*fl2p1(l)
             ENDDO
             !--->       loop over equivalent atoms
             n1 = 0
             DO n = 1,atoms%ntype

                IF (noco%l_noco) THEN
                   !--->          pk non-collinear
                   !--->             set up the spinors of this atom within global
                   !--->             spin-coordinateframe
                   call hsmt_spinor(isp,n, noco, input,chi, chi11, chi21, chi22)
                ENDIF
                !---> pk non-collinear
                n0 = n1 + 1
                n1 = n1 + atoms%neq(n)
                rph(:,1) = 0.0
                cph(:,1) = 0.0
                DO nn = n0,n1
                   tnn = tpi_const*atoms%taual(:,nn)
                   !--->             set up phase factors
                   DO kj = 1,kjmax
                      rph(kj,1) = rph(kj,1) +&
                           cos(dot_product(ski-(/lapw%k1(kj,jintsp),lapw%k2(kj,jintsp),lapw%k3(kj,jintsp)/)+qssbtj,tnn))

#ifndef CPP_INVERSION
                      !--->                if the system does not posses inversion symmetry
                      !--->                the complex part of the exponential is needed.
                      cph(kj,1) = cph(kj,1) +&
                           sin(dot_product((/lapw%k1(kj,jintsp),lapw%k2(kj,jintsp),lapw%k3(kj,jintsp)/)+qssbtj-ski,tnn))
#endif
                   END DO
                END DO

                !--->          update overlap and l-diagonal hamiltonian matrix
                DO  l = 0,atoms%lmax(n)
                   IF (noco%l_constr.or.l_socfirst) THEN
                      fjkiln = fj(ki,l,n,isp)
                      gjkiln = gj(ki,l,n,isp)*usdus%ddn(l,n,isp)
                   ELSE
                      fjkiln = fj(ki,l,n,iintsp)
                      gjkiln = gj(ki,l,n,iintsp)*usdus%ddn(l,n,isp)
                   ENDIF
                   !
                   !
                   
              
                   IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
                      !--->             pk non-collinear
#ifndef CPP_INVERSION
                      IF (noco%l_constr.or.l_socfirst) THEN
                         DO kj = 1,ki
                            fct  = plegend(kj,l)*&
                                 ( fjkiln*fj(kj,l,n,isp) + gjkiln*gj(kj,l,n,isp) )

                            bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                         ENDDO
                      ELSE
                         DO kj = 1,ki
                            fct  = plegend(kj,l)*&
                                 ( fjkiln*fj(kj,l,n,jintsp) + gjkiln*gj(kj,l,n,jintsp) )

                            bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                         ENDDO
                      ENDIF
                      !+||
                      IF ( kii+1.LE.lapw%nv(1) ) THEN
                         !--->                spin-up spin-up part
                         CALL CPP_BLAS_caxpy(ki,chi11,bbwa,1,bb(iii+1),1)
                         !--->                spin-down spin-up part, upper triangle.
                         !--->                the help array is used to allow vectorization on
                         !--->                Cray PVP systems. it is mapped onto the hamiltonian
                         !--->                matrix after the setup is complete.
                         DO kj = 1,ki - 1
                            ii = iii + kj
                            bbhlp(ii)=bbhlp(ii)+conjg(bbwa(kj))*chi21
                         END DO
                      ENDIF
                      IF ( (kii+1.GT.lapw%nv(1)).OR.(n_size.EQ.1) ) THEN
                         IF (n_size.EQ.1) THEN
                            ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2
                         ELSE
                            ii = iii
                         ENDIF
                         !--->                spin-down spin-up part, lower triangle
                         CALL CPP_BLAS_caxpy(ki,chi21,bbwa,1,bb(ii+1),1)
                         !--->                spin-down spin-down part
                         ii = ii + lapw%nv(1)+atoms%nlotot
                         CALL CPP_BLAS_caxpy(ki,chi22,bbwa,1,bb(ii+1),1)
                      ENDIF
                      !-||
                      !--->                when fj and gj are available for both local spins
                      !--->                (isp), add the contribution from the constraint
                      !--->                B-field.                      
                   ELSEIF ( noco%l_noco .AND. noco%l_ss ) THEN
                      IF ( iintsp.EQ.2 .AND. jintsp.EQ.1 ) THEN
                         kjmax = lapw%nv(1)
                      ELSE
                         kjmax = ki
                      ENDIF
                      DO kj = 1,kjmax
                         fct  = plegend(kj,l)* ( fjkiln*fj(kj,l,n,jintsp) +&
                              gjkiln*gj(kj,l,n,jintsp) )

                         bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                      ENDDO
                      IF ( iintsp.EQ.1 .AND. jintsp.EQ.1 ) THEN
                         !--->                   spin-up spin-up part
                         ii = (ki-1)*(ki)/2
                         DO kj = 1,ki
                            bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi11
                         ENDDO
                      ELSEIF ( iintsp.EQ.2 .AND. jintsp.EQ.2 ) THEN
                         !--->                   spin-down spin-down part
                         ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2 +&
                              lapw%nv(1)+atoms%nlotot
                         DO kj = 1,ki
                            bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi22
                         ENDDO
                      ELSE
                         !--->                   spin-down spin-up part
                         ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2
                         DO kj = 1,lapw%nv(1)
                            bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi21
                         ENDDO
                      ENDIF
#endif
                      !--->             pk non-collinear
                   ELSE
                      DO kj = 1,ki
                         fct  = plegend(kj,l)*( fjkiln*fj(kj,l,n,jintsp) + gjkiln*gj(kj,l,n,jintsp) )

                         ij = iii + kj
#ifdef CPP_INVERSION
                         bb(ij) = bb(ij) + rph(kj,1) * fct
                         !-APW
#else
                         bb(ij) = bb(ij) + cmplx(rph(kj,1),cph(kj,1))*fct
#endif
                      END DO
                   ENDIF
                   !--->          end loop over l
                enddo
                !--->       end loop over atom types (ntype)
             enddo
             !--->    end loop over ki
          enddo
          !$OMP END DO
          !---> end loops over interstitial spins
       ENDDO ! jintsp = 1,iintsp
    ENDDO   ! iintsp = 1,nintsp
    deallocate(plegend)
    deallocate(rph,cph)
    !$OMP END PARALLEL


    !---> pk non-collinear
    IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
       !--->    map the hamiltonian help array onto the hamitonian matrix
       IF (n_size.EQ.1) THEN
          DO ki = 1,lapw%nv(1)
             ii = (ki-1)*(ki)/2
             DO kj = 1,ki-1
                ij = (lapw%nv(1)+atoms%nlotot+kj-1)*(lapw%nv(1)+atoms%nlotot+kj)/2 + ki
                bb(ij) = bb(ij) + bbhlp(ii+kj)
             ENDDO
          ENDDO

       ELSE
#ifdef CPP_MPI
 !         CALL mingeselle(SUB_COMM,n_size,n_rank,lapw%nv, bbhlp, bb)
!
#endif
       ENDIF
    ENDIF

    RETURN
  END SUBROUTINE hsmt_overlap

  subroutine hsmt_overlap_zherk(input,sym,atoms,n_size,n_rank,isp,l_socfirst,hlpmsize,ab_dim,&
       noco,cell,nintsp, lapw,usdus,gk,vk,fj,gj,bb)
!Calculate overlap matrix
#include"cpp_double.h"
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_types
    USE m_ylm
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(INOUT)  :: lapw!lpaw%nv_tot is updated
    TYPE(t_usdus),INTENT(INOUT) :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: isp
    INTEGER, INTENT (IN) :: n_size,n_rank,ab_dim
    INTEGER, INTENT (IN) :: hlpmsize,nintsp
    LOGICAL, INTENT (IN) :: l_socfirst
    !     ..
    !     .. Array Arguments ..
    REAL,INTENT(IN) :: gk(:,:,:),vk(:,:,:)
    REAL,INTENT(IN) :: fj(:,0:,:,:),gj(:,0:,:,:)
#ifdef CPP_INVERSION
    REAL,    INTENT (INOUT) :: bb(:)!(matsize)
#else
    COMPLEX, INTENT (INOUT) :: bb(:)
#endif
    
    INTEGER:: n,nn,na,np,spin2,iintsp,k,l,ll1,m,i,ki,kj
    complex:: term
    real   :: th,v(3),bmrot(3,3),vmult(3)
    complex,allocatable:: c_ph(:,:),ylm(:),bb_tmp(:,:),a(:,:,:),b(:,:,:)
    real,allocatable   :: gkrot(:,:)
    
    ALLOCATE(c_ph(lapw%nv(1),1),ylm((atoms%lmaxd+1)**2),bb_tmp(lapw%nv(1),lapw%nv(1)))
    allocate(a(lapw%nv(1),(atoms%lmaxd+1)**2,1),b(lapw%nv(1),(atoms%lmaxd+1)**2,1))
    allocate(gkrot(lapw%nv(1),3))


    ntyploop: DO n=1,atoms%ntype
       DO nn = 1,atoms%neq(n)
	  a=0.0
          b=0.0
          na = SUM(atoms%neq(:n-1))+nn
          IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
             np = sym%invtab(atoms%ngopr(na))
             !--->       loop over interstitial spins
             DO iintsp = 1,nintsp
                IF (noco%l_constr.OR.l_socfirst) THEN
                   spin2=isp
                ELSE
                   spin2=iintsp
                ENDIF
                !--->          set up phase factors
                DO k = 1,lapw%nv(iintsp)
                   th= DOT_PRODUCT((/lapw%k1(k,iintsp),lapw%k2(k,iintsp),lapw%k3(k,iintsp)/)+(iintsp-1.5)*noco%qss,atoms%taual(:,na))
                   c_ph(k,iintsp) = CMPLX(COS(tpi_const*th),-SIN(tpi_const*th))
                END DO

                IF (np==1) THEN
                   gkrot( 1:lapw%nv(iintsp),:) = gk( 1:lapw%nv(iintsp),:,iintsp)
                ELSE
                   bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
                   DO k = 1,lapw%nv(iintsp)
                      !-->  apply the rotation that brings this atom into the
                      !-->  representative (this is the definition of ngopr(na)
                      !-->  and transform to cartesian coordinates
                      v(:) = vk(k,:,iintsp)
                      gkrot(k,:) = MATMUL(TRANSPOSE(bmrot),v)
                   END DO
                END IF
                DO k = 1,lapw%nv(1)
                   !-->    generate spherical harmonics
                   vmult(:) =  gkrot(k,:)
                   CALL ylm4(atoms%lmax(n),vmult,ylm)
                   !-->  synthesize the complex conjugates of a and b
                   DO l = 0,atoms%lmax(n)
                      ll1 = l* (l+1)
                      DO m = -l,l

                         term = c_ph(k,iintsp)*ylm(ll1+m+1)
                         a(k,ll1+m+1,iintsp) = fj(k,l,n,spin2)*term
                         b(k,ll1+m+1,iintsp) = gj(k,l,n,spin2)*term
                      END DO
                   END DO
                ENDDO !k-loop
                !--->       end loop over interstitial spin
                call timestart("zherk")
                CALL ZHERK("U","N",lapw%nv(iintsp),size(a,2),cmplx(1.,0),a(:,:,iintsp),size(a,1),cmplx(0.0,0.0),bb_tmp,size(bb_tmp,1))
                CALL ZHERK("U","N",lapw%nv(iintsp),size(a,2),cmplx(1.,0),b(:,:,iintsp),size(a,1),cmplx(1.0,0.0),bb_tmp,size(bb_tmp,1))
                call timestop("zherk")

                i=1
                DO ki=1,lapw%nv(1)
                   DO kj=ki,lapw%nv(1)
                      bb(i)=bb(i)+bb_tmp(ki,kj)
                      i=i+1
                   ENDDO
                ENDDO
             ENDDO
          end IF
       end DO
    end DO ntyploop
  end subroutine hsmt_overlap_zherk
END MODULE m_hsmt_overlap
