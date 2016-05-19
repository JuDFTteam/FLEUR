MODULE m_hsmt_sph
  use m_juDFT
  implicit none
CONTAINS
  SUBROUTINE hsmt_sph(dimension,atoms,SUB_COMM,n_size,n_rank,sphhar,isp,ab_dim,&
       input,hlpmsize,noco,l_socfirst,cell,nintsp, lapw,el,usdus,vr,gk,rsoc,isigma, aa,bb,fj,gj)

#include"cpp_double.h"
    USE m_constants, ONLY : fpi_const,tpi_const
    USE m_sphbes
    USE m_dsphbs
    USE m_ylm
    USE m_hsmt_socinit,ONLY:t_rsoc
    USE m_hsmt_spinor
    USE m_radovlp
#ifdef CPP_MPI
    USE m_mingeselle
#endif
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN):: dimension
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_sphhar),INTENT(IN)   :: sphhar
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(INOUT)  :: lapw!lpaw%nv_tot is updated
    TYPE(t_usdus),INTENT(INOUT) :: usdus
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: isp,ab_dim
    INTEGER, INTENT (IN) :: SUB_COMM,n_size,n_rank 
    INTEGER, INTENT (IN) :: hlpmsize,nintsp
    LOGICAL, INTENT (IN) :: l_socfirst
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntypd,dimension%jspd)
    REAL,    INTENT (IN) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,dimension%jspd)
    REAL,INTENT(IN)      :: gk(:,:,:)
    COMPLEX,INTENT(IN)   :: isigma(2,2,3)
    TYPE(t_rsoc),INTENT(IN) :: rsoc


#ifdef CPP_INVERSION
    REAL,    INTENT (INOUT) :: aa(:),bb(:)!(matsize)
#else
    COMPLEX, INTENT (INOUT) :: aa(:),bb(:)
#endif
    REAL,INTENT(OUT) :: fj(:,0:,:,:),gj(:,0:,:,:)
    !     ..
    !     .. Local Scalars ..
    REAL con1,ff,gg,gs,ws,tnn(3), elall,fct,fjkiln,gjkiln,ddnln,ski(3)
    REAL apw_lo1,apw_lo2,apw1,w1

    COMPLEX chi11,chi21,chi22,capw1
    INTEGER ii,iii,ij,k,ki,kj,l,lo,n,n0,n1,nn,kjmax, nsp,iintsp,jintsp
    INTEGER nc ,kii

    !     ..
    !     .. Local Arrays ..
    REAL fb(0:atoms%lmaxd),fleg1(0:atoms%lmaxd),fleg2(0:atoms%lmaxd),fl2p1(0:atoms%lmaxd)     
    REAL fl2p1bt(0:atoms%lmaxd),gb(0:atoms%lmaxd)
    REAL qssbti(3),qssbtj(3)


    REAL, ALLOCATABLE :: plegend(:,:)
    REAL, ALLOCATABLE :: rph(:,:),cph(:,:)
    REAL, ALLOCATABLE :: uun21(:,:),udn21(:,:),dun21(:,:),ddn21(:,:)


    COMPLEX chi(2,2),chj(2,2,2,atoms%ntypd),aawa(dimension%nvd),bbwa(dimension%nvd)
    COMPLEX, ALLOCATABLE :: aahlp(:),bbhlp(:)
    LOGICAL apw(0:atoms%lmaxd)


    ! for Spin-orbit...
    REAL, ALLOCATABLE :: dplegend(:,:)
    REAL, ALLOCATABLE :: cross_k(:,:)
    INTEGER :: j1,j2
    COMPLEX :: isigma_x(2,2),isigma_y(2,2),isigma_z(2,2)
    COMPLEX :: chi11so(2,2),chi21so(2,2),chi22so(2,2),angso(dimension%nvd,2,2)


    IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) ALLOCATE ( aahlp(hlpmsize),bbhlp(hlpmsize) )
    !     ..
    con1 = fpi_const/sqrt(cell%omtil)
    DO l = 0,atoms%lmaxd
       fleg1(l) = real(l+l+1)/real(l+1)
       fleg2(l) = real(l)/real(l+1)
       fl2p1(l) = real(l+l+1)/fpi_const
       fl2p1bt(l) = fl2p1(l)*0.5
    END DO
    !---> calculate the overlap of the radial basis functions with different
    !---> spin directions.
    IF (noco%l_constr) THEN
       ALLOCATE(uun21(0:atoms%lmaxd,atoms%ntypd),udn21(0:atoms%lmaxd,atoms%ntypd),&
            dun21(0:atoms%lmaxd,atoms%ntypd),ddn21(0:atoms%lmaxd,atoms%ntypd) )
       CALL rad_ovlp(atoms,usdus,input,vr,el(0:,:,:), uun21,udn21,dun21,ddn21)
    ENDIF
    !---> loop over each atom type

    DO iintsp = 1,nintsp
       !$OMP parallel do DEFAULT(SHARED)&
       !$OMP PRIVATE(n,l,apw,lo,k,gs,fb,gb,ws,ff,gg)
       DO n = 1,atoms%ntype

          DO l = 0,atoms%lmax(n)
             apw(l) = .false.
             DO lo = 1,atoms%nlo(n)
                IF (atoms%l_dulo(lo,n)) apw(l) = .true.
             ENDDO
#ifdef CPP_APW
             IF (atoms%lapw_l(n).GE.l) apw(l) = .false.
#endif
          ENDDO
          DO lo = 1,atoms%nlo(n)
             IF (atoms%l_dulo(lo,n)) apw(atoms%llo(lo,n)) = .true.
          ENDDO

          DO k = 1,lapw%nv(iintsp)
             gs = lapw%rk(k,iintsp)*atoms%rmt(n)
             CALL sphbes(atoms%lmax(n),gs, fb)
             CALL dsphbs(atoms%lmax(n),gs,fb, gb)
             DO l = 0,atoms%lmax(n)
                !---> set up wronskians for the matching conditions for each ntype
                ws = con1/(usdus%uds(l,n,isp)*usdus%dus(l,n,isp)&
                     - usdus%us(l,n,isp)*usdus%duds(l,n,isp))
                ff = fb(l)
                gg = lapw%rk(k,iintsp)*gb(l)
                IF ( apw(l) ) THEN
                   fj(k,l,n,iintsp) = 1.0*con1 * ff / usdus%us(l,n,isp)
                   gj(k,l,n,iintsp) = 0.0d0
                ELSE
                   IF (noco%l_constr.or.l_socfirst) THEN
                      !--->                 in a constrained calculation fj and gj are needed
                      !--->                 both local spin directions (isp) at the same time
                      fj(k,l,n,isp) = ws * ( usdus%uds(l,n,isp)*gg - usdus%duds(l,n,isp)*ff )
                      gj(k,l,n,isp) = ws * ( usdus%dus(l,n,isp)*ff - usdus%us(l,n,isp)*gg )
                   ELSE
                      !--->                 in a spin-spiral calculation fj and gj are needed
                      !--->                 both interstitial spin directions at the same time
                      !--->                 In all other cases iintsp runs from 1 to 1.
                      fj(k,l,n,iintsp) = ws * ( usdus%uds(l,n,isp)*gg - usdus%duds(l,n,isp)*ff )
                      gj(k,l,n,iintsp) = ws * ( usdus%dus(l,n,isp)*ff - usdus%us(l,n,isp)*gg )
                   ENDIF
                ENDIF
             ENDDO
          ENDDO ! k = 1, lapw%nv
       ENDDO    ! n = 1,atoms%ntype
       !$OMP end parallel do

    ENDDO       ! iintsp = 1,nintsp
    !
    IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
       !---> pk non-collinear
       !--->    initialize hamiltonian help array
       aahlp=cmplx(0.0,0.0)
       bbhlp=cmplx(0.0,0.0)
    ENDIF

    !$OMP PARALLEL  DEFAULT(shared)&
    !$OMP PRIVATE(kii,ki,nc,iii,kjmax,ski,kj,plegend,l,n1,n)&
    !$OMP PRIVATE(n0,rph,cph,nn,tnn,fjkiln,gjkiln)&
    !$OMP PRIVATE(w1,apw_lo1,apw_lo2,ddnln,elall,fct,ij,apw1)&
    !$OMP PRIVATE(cross_k,dplegend,chi,chi11,chi21,chi22,nsp,chj)&
    !$OMP PRIVATE(isigma_x,isigma_y,isigma_z,j1,j2,chi11so,chi21so,chi22so)&
    !$OMP PRIVATE(aawa,bbwa,capw1,ii) IF (.not.l_socfirst)
    !$     IF (l_socfirst) WRITE(*,*) "WARNING: first variation SOC does not work with OPENMP in hsmt_sph"
    !$     IF (l_socfirst) WRITE(*,*) "         switching off openmp parallelization"
    ALLOCATE(rph(dimension%nvd,ab_dim))
    ALLOCATE(cph(dimension%nvd,ab_dim))
    ALLOCATE(plegend(dimension%nvd,0:atoms%lmaxd))
    IF (l_socfirst)THEN
       ALLOCATE ( dplegend(dimension%nvd,0:atoms%lmaxd),cross_k(dimension%nvd,3))
       dplegend(:,0)=0.e0
       dplegend(:,1)=1.e0
    ENDIF

    plegend=0.0
    plegend(:,0)=1.0
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
             IF (noco%l_ss)  CALL juDFT_error("ev-|| & spin-spiral !",calledby ="hssphn")
          ELSE
             lapw%nv_tot = lapw%nv(iintsp)
          ENDIF

          !
          !$OMP  DO SCHEDULE(DYNAMIC,1)
          DO  kii =  n_rank, lapw%nv_tot-1, n_size
             ki = mod(kii,lapw%nv(iintsp)) + 1
             nc = nc + 1
             !$          nc= 1+ (kii-n_rank)/n_size
             iii = nc*(nc-1)/2 * n_size - (nc-1)*(n_size - n_rank - 1)
             !            ii = nc*(nc+1)/2 * n_size -  nc   *(n_size - n_rank - 1)

             IF (noco%l_ss.OR.noco%l_constr.OR.l_socfirst) THEN
                kjmax = lapw%nv(jintsp)
             ELSE
                kjmax = ki
             ENDIF
             ski = (/lapw%k1(ki,iintsp),lapw%k2(ki,iintsp),lapw%k3(ki,iintsp)/) + qssbti
             !--->       legendre polynomials
             DO kj = 1,kjmax
                plegend(kj,1) = dot_product(gk(kj,:,jintsp),gk(ki,:,iintsp))
                IF (l_socfirst) THEN
                   !#ifdef CPP_SOCFIRST
                   cross_k(kj,1)=gk(ki,2,jintsp)*gk(kj,3,iintsp)- gk(ki,3,jintsp)*gk(kj,2,iintsp)
                   cross_k(kj,2)=gk(ki,3,jintsp)*gk(kj,1,iintsp)- gk(ki,1,jintsp)*gk(kj,3,iintsp)
                   cross_k(kj,3)=gk(ki,1,jintsp)*gk(kj,2,iintsp)- gk(ki,2,jintsp)*gk(kj,1,iintsp)
                ENDIF
                !#endif
             END DO
             DO l = 1,maxval(atoms%lmax) - 1
                plegend(:,l+1) = fleg1(l)*plegend(:,1)*plegend(:,l) - fleg2(l)*plegend(:,l-1)
                IF (l_socfirst) THEN
                   !#ifdef CPP_SOCFIRST
                   dplegend(:,l+1)=REAL(l+1)*plegend(:,l)+ plegend(:,1)*dplegend(:,l)
                   !#endif
                ENDIF
             END DO
             !--->       loop over equivalent atoms
             n1 = 0
             DO n = 1,atoms%ntype

                IF (noco%l_noco) THEN
                   !--->          pk non-collinear
                   !--->             set up the spinors of this atom within global
                   !--->             spin-coordinateframe
                   call hsmt_spinor(isp,n, noco, input,chi, chi11, chi21, chi22,l_socfirst,&
                        isigma,isigma_x,isigma_y,isigma_z,chi11so,chi21so,chi22so,angso,chj,cross_k)
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
                      gjkiln = gj(ki,l,n,isp)
                   ELSE
                      fjkiln = fj(ki,l,n,iintsp)
                      gjkiln = gj(ki,l,n,iintsp)
                   ENDIF
                   !
                   w1 = 0.5 * ( usdus%uds(l,n,isp)*usdus%dus(l,n,isp) + &
                        usdus%us(l,n,isp)*usdus%duds(l,n,isp) )
                   apw_lo1 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( gjkiln * w1 +&
                        fjkiln * usdus%us(l,n,isp) * usdus%dus(l,n,isp) )
                   apw_lo2 = fl2p1(l) * 0.5 * atoms%rmt(n)**2 * ( fjkiln * w1 +&
                        gjkiln * usdus%uds(l,n,isp) * usdus%duds(l,n,isp) )
                   !
                   ddnln =  usdus%ddn(l,n,isp)
                   elall = el(l,n,isp)

                   IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
                      !--->             pk non-collinear
#ifndef CPP_INVERSION
                      IF (noco%l_constr.or.l_socfirst) THEN
                         DO kj = 1,ki
                            fct  = plegend(kj,l)*fl2p1(l)*&
                                 ( fjkiln*fj(kj,l,n,isp) + gjkiln*gj(kj,l,n,isp)*ddnln )

                            bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                            aawa(kj) = cmplx(rph(kj,1),cph(kj,1)) * ( &
                                 fct*elall + plegend(kj,l)*fl2p1bt(l)*&
                                 ( fjkiln*gj(kj,l,n,isp) + gjkiln*fj(kj,l,n,isp) ) )
#ifdef CPP_APW
                            capw1 = cmplx(rph(kj,1),cph(kj,1))*plegend(kj,l)&
                                 * ( apw_lo1 * fj(kj,l,n,isp) + apw_lo2 * gj(kj,l,n,isp) )
                            aawa(kj) = aawa(kj) + capw1
#endif
                         ENDDO
                      ELSE
                         DO kj = 1,ki
                            fct  = plegend(kj,l)*fl2p1(l)*&
                                 ( fjkiln*fj(kj,l,n,jintsp) + gjkiln*gj(kj,l,n,jintsp)*ddnln )

                            bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                            aawa(kj) = cmplx(rph(kj,1),cph(kj,1)) * (&
                                 fct*elall + plegend(kj,l)*fl2p1bt(l)*&
                                 ( fjkiln*gj(kj,l,n,jintsp) + gjkiln*fj(kj,l,n,jintsp) ) )
#ifdef CPP_APW
                            capw1 = cmplx(rph(kj,1),cph(kj,1))*plegend(kj,l)&
                                 * ( apw_lo1 * fj(kj,l,n,jintsp) + apw_lo2 * gj(kj,l,n,jintsp) )
                            aawa(kj) = aawa(kj) + capw1
#endif
                         ENDDO
                      ENDIF
                      !+||
                      IF ( kii+1.LE.lapw%nv(1) ) THEN
                         !--->                spin-up spin-up part
                         CALL CPP_BLAS_caxpy(ki,chi11,bbwa,1,bb(iii+1),1)
                         CALL CPP_BLAS_caxpy(ki,chi11,aawa,1,aa(iii+1),1)
                         !--->                spin-down spin-up part, upper triangle.
                         !--->                the help array is used to allow vectorization on
                         !--->                Cray PVP systems. it is mapped onto the hamiltonian
                         !--->                matrix after the setup is complete.
                         DO kj = 1,ki - 1
                            ii = iii + kj
                            aahlp(ii)=aahlp(ii)+conjg(aawa(kj))*chi21
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
                         CALL CPP_BLAS_caxpy(ki,chi21,aawa,1,aa(ii+1),1)
                         !--->                spin-down spin-down part
                         ii = ii + lapw%nv(1)+atoms%nlotot
                         CALL CPP_BLAS_caxpy(ki,chi22,bbwa,1,bb(ii+1),1)
                         CALL CPP_BLAS_caxpy(ki,chi22,aawa,1,aa(ii+1),1)
                      ENDIF
                      !-||
                      !--->                when fj and gj are available for both local spins
                      !--->                (isp), add the contribution from the constraint
                      !--->                B-field.
                      IF (noco%l_constr .AND. (isp .EQ. 2)) THEN
                         DO nsp = 1,input%jspins
                            IF (nsp.EQ.1) THEN
                               DO kj = 1,lapw%nv(1)
                                  aawa(kj) = (-0.5)*cmplx(noco%b_con(1,n),noco%b_con(2,n))*&
                                       cmplx(rph(kj,1),cph(kj,1))*&
                                       plegend(kj,l)*fl2p1(l)*(&
                                       + fj(ki,l,n,2)*fj(kj,l,n,1)*uun21(l,n)&
                                       + fj(ki,l,n,2)*gj(kj,l,n,1)*udn21(l,n)&
                                       + gj(ki,l,n,2)*fj(kj,l,n,1)*dun21(l,n)&
                                       + gj(ki,l,n,2)*gj(kj,l,n,1)*ddn21(l,n))
                               ENDDO
                            ELSE
                               DO kj = 1,lapw%nv(1)
                                  aawa(kj) = (-0.5)*cmplx(noco%b_con(1,n),-noco%b_con(2,n))*&
                                       cmplx(rph(kj,1),cph(kj,1))*&
                                       plegend(kj,l)*fl2p1(l)*(&
                                       + fj(ki,l,n,1)*fj(kj,l,n,2)*uun21(l,n)&
                                       + fj(ki,l,n,1)*gj(kj,l,n,2)*dun21(l,n)&
                                       + gj(ki,l,n,1)*fj(kj,l,n,2)*udn21(l,n)&
                                       + gj(ki,l,n,1)*gj(kj,l,n,2)*ddn21(l,n))
                               ENDDO
                            ENDIF
                            !--->                  spin-up spin-up part
                            ii = (ki-1)*(ki)/2
                            DO kj = 1,ki
                               aa(ii+kj) = aa(ii+kj) + aawa(kj)*chj(nsp,1,1,n)
                            ENDDO
                            !--->                  spin-down spin-down part
                            ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2 + &
                                 lapw%nv(1)+atoms%nlotot
                            DO kj = 1,ki
                               aa(ii+kj) = aa(ii+kj) + aawa(kj)*chj(nsp,2,2,n)
                            ENDDO
                            !--->                  spin-down spin-up part
                            ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2
                            DO kj = 1,lapw%nv(1)
                               aa(ii+kj) = aa(ii+kj) + aawa(kj)*chj(nsp,2,1,n)
                            ENDDO
                         ENDDO
                      ENDIF

                      !                 Add spin-orbit coupling
                      IF (isp.EQ.2) THEN
                         IF ((l.GT.0).AND.l_socfirst) THEN

                            DO j1 = 1,input%jspins
                               DO j2 = 1,input%jspins
                                  DO kj = 1,lapw%nv(1)
                                     aawa(kj)=cmplx(rph(kj,1),cph(kj,1))*(&
                                          + fj(ki,l,n,j1)*fj(kj,l,n,j2)*rsoc%rsopp(n,l,j1,j2)&
                                          + fj(ki,l,n,j1)*gj(kj,l,n,j2)*rsoc%rsopdp(n,l,j1,j2)&
                                          + gj(ki,l,n,j1)*fj(kj,l,n,j2)*rsoc%rsoppd(n,l,j1,j2)&
                                          + gj(ki,l,n,j1)*gj(kj,l,n,j2)*rsoc%rsopdpd(n,l,j1,j2))&
                                          *dplegend(kj,l)*fl2p1(l)*angso(kj,j1,j2)
                                  ENDDO
                                  IF (n_size.EQ.1) THEN

                                     !--->                  spin-up spin-up part
                                     ii = (ki-1)*(ki)/2
                                     DO kj = 1,ki
                                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi11so(j1,j2)
                                     ENDDO
                                     !--->                  spin-down spin-down part
                                     ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2 +&
                                          lapw%nv(1)+atoms%nlotot
                                     DO kj = 1,ki
                                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi22so(j1,j2)
                                     ENDDO
                                     !--->                  spin-down spin-up part
                                     ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2
                                     DO kj = 1,lapw%nv(1)
                                        aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi21so(j1,j2)
                                     ENDDO

                                  ELSE  ! eigenvalue parallelization

                                     IF ( kii+1.LE.lapw%nv(1) ) THEN
                                        !--->                    spin-up spin-up part
                                        CALL CPP_BLAS_caxpy(ki,chi11so(j1,j2),aawa,1, aa(iii+1),1)

                                        !--->                    spin-down spin-up part, upper triangle.
                                        DO kj = 1,ki - 1
                                           ii = iii + kj
                                           aahlp(ii) = aahlp(ii) + conjg(aawa(kj))*chi21so(j2,j1)
                                        END DO
                                     ENDIF
                                     IF  (kii+1.GT.lapw%nv(1)) THEN
                                        ii = iii
                                        !--->                    spin-down spin-up part, lower triangle
                                        CALL CPP_BLAS_caxpy(ki,chi21so(j1,j2),aawa,1, aa(ii+1),1)
                                        !--->                    spin-down spin-down part
                                        ii = ii + lapw%nv(1)+atoms%nlotot
                                        CALL CPP_BLAS_caxpy(ki,chi22so(j1,j2),aawa,1, aa(ii+1),1)
                                     ENDIF

                                  ENDIF ! eigenvalue par.

                               ENDDO ! j2
                            ENDDO    ! j1
                         ENDIF       ! ( l > 0 ) & socfirst
                      ENDIF          ! (isp = 2)
                      !               End spin-orbit
                      !#endif
                   ELSEIF ( noco%l_noco .AND. noco%l_ss ) THEN
                      IF ( iintsp.EQ.2 .AND. jintsp.EQ.1 ) THEN
                         kjmax = lapw%nv(1)
                      ELSE
                         kjmax = ki
                      ENDIF
                      DO kj = 1,kjmax
                         fct  = plegend(kj,l)*fl2p1(l)* ( fjkiln*fj(kj,l,n,jintsp) +&
                              gjkiln*gj(kj,l,n,jintsp)*ddnln )

                         bbwa(kj) = cmplx(rph(kj,1),cph(kj,1))*fct
                         aawa(kj) = cmplx(rph(kj,1),cph(kj,1)) * ( &
                              fct*elall + plegend(kj,l)*fl2p1bt(l)*&
                              ( fjkiln*gj(kj,l,n,jintsp) + gjkiln*fj(kj,l,n,jintsp) ) )
                      ENDDO
                      IF ( iintsp.EQ.1 .AND. jintsp.EQ.1 ) THEN
                         !--->                   spin-up spin-up part
                         ii = (ki-1)*(ki)/2
                         DO kj = 1,ki
                            bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi11
                            aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi11
                         ENDDO
                      ELSEIF ( iintsp.EQ.2 .AND. jintsp.EQ.2 ) THEN
                         !--->                   spin-down spin-down part
                         ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2 +&
                              lapw%nv(1)+atoms%nlotot
                         DO kj = 1,ki
                            bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi22
                            aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi22
                         ENDDO
                      ELSE
                         !--->                   spin-down spin-up part
                         ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2
                         DO kj = 1,lapw%nv(1)
                            bb(ii+kj) = bb(ii+kj) + bbwa(kj)*chi21
                            aa(ii+kj) = aa(ii+kj) + aawa(kj)*chi21
                         ENDDO
                      ENDIF
#endif
                      !--->             pk non-collinear
                   ELSE
                      DO kj = 1,ki
                         fct  = plegend(kj,l)*fl2p1(l)*&
                              ( fjkiln*fj(kj,l,n,jintsp) + gjkiln*gj(kj,l,n,jintsp)*ddnln )

                         ij = iii + kj
#ifdef CPP_INVERSION
                         bb(ij) = bb(ij) + rph(kj,1) * fct
                         aa(ij) = aa(ij) + rph(kj,1) * ( fct * elall + plegend(kj,l) * fl2p1bt(l) *&
                              ( fjkiln*gj(kj,l,n,jintsp) + gjkiln*fj(kj,l,n,jintsp) ) )
                         !+APW
#ifdef CPP_APW
                         apw1 = rph(kj,1) * plegend(kj,l)  * &
                              ( apw_lo1 * fj(kj,l,n,iintsp) + apw_lo2 * gj(kj,l,n,iintsp) )
                         aa(ij) = aa(ij) + apw1
#endif
                         !-APW
#else
                         bb(ij) = bb(ij) + cmplx(rph(kj,1),cph(kj,1))*fct
                         aa(ij) = aa(ij) + cmplx(rph(kj,1),cph(kj,1)) * (fct*elall + plegend(kj,l)*fl2p1bt(l) *&
                              ( fjkiln*gj(kj,l,n,jintsp) + gjkiln*fj(kj,l,n,jintsp) ) )
#ifdef CPP_APW
                         capw1 = cmplx(rph(kj,1),cph(kj,1))*plegend(kj,l)&
                              * ( apw_lo1 * fj(kj,l,n,iintsp) + apw_lo2 * gj(kj,l,n,iintsp) )
                         aa(ij) = aa(ij) + capw1
#endif
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
    IF (l_socfirst) DEALLOCATE(dplegend,cross_k)
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
                aa(ij) = aa(ij) + aahlp(ii+kj)
                bb(ij) = bb(ij) + bbhlp(ii+kj)
             ENDDO
          ENDDO

       ELSE
#ifdef CPP_MPI
          CALL mingeselle(SUB_COMM,n_size,n_rank,lapw%nv, aahlp, aa)
          CALL mingeselle(SUB_COMM,n_size,n_rank,lapw%nv, bbhlp, bb)
#endif
       ENDIF
    ENDIF



    RETURN
  END SUBROUTINE hsmt_sph
END MODULE m_hsmt_sph
