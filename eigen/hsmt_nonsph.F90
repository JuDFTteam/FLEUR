MODULE m_hsmt_nonsph
#define CPP_BLOCKSIZE 64
  !      USE m_juDFT
  !$     USE omp_lib

  !TODO:
  !  Check what can be done in l_noco=.true. case in terms of use of zgemm or aa_block
  !  Check what happens in case of CPP_INVERSION -> real matrix a

  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_nonsph(DIMENSION,atoms,sym,SUB_COMM, n_size,n_rank,input,isp,nintsp,&
       hlpmsize,noco,l_socfirst, lapw, cell,tlmplm, fj,gj,gk,vk,oneD,l_real,aa_r,aa_c)

#include"cpp_double.h"
    USE m_constants, ONLY : tpi_const
    USE m_ylm
    USE m_hsmt_spinor
    USE m_hsmt_hlptomat
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN):: DIMENSION
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(INOUT)  :: lapw !lapw%nv_tot is updated

    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nintsp,isp
    INTEGER, INTENT (IN) :: SUB_COMM,n_size,n_rank 
    INTEGER, INTENT (IN) :: hlpmsize
    LOGICAL, INTENT (IN) :: l_socfirst
    !     ..
    !     .. Array Arguments ..
    TYPE(t_tlmplm),INTENT(IN)::tlmplm
    REAL, INTENT(IN)     :: fj(:,0:,:,:),gj(:,0:,:,:)
    REAL,INTENT(IN)      :: gk(:,:,:),vk(:,:,:)
    !-odim
    !+odim
    LOGICAL, INTENT(IN)     :: l_real
    REAL,   ALLOCATABLE, INTENT (INOUT) :: aa_r(:)!(matsize)
    COMPLEX,ALLOCATABLE, INTENT (INOUT) :: aa_c(:)
    COMPLEX,PARAMETER :: one=CMPLX(1.0,0.0),zero=CMPLX(0.0,0.0)

    !     ..
    !     .. Local Scalars ..
    INTEGER :: i,iii,ii,ij,im,in,k,ki,kj,l,ll1,lm,lmp,lp,jd,m
    INTEGER :: mp,n,na,nn,np,kjmax,iintsp,jintsp
    INTEGER :: nc ,kii,spin2,ab_dim,lnonsphd,bsize,bsize2,kb
    REAL    :: th,invsfct
    COMPLEX :: term,chi11,chi21,chi22,chihlp


    !     ..
    !     .. Local Arrays ..
    COMPLEX,ALLOCATABLE :: aa_block(:,:)
    COMPLEX,ALLOCATABLE :: dtd(:,:),dtu(:,:),utd(:,:),utu(:,:)
    REAL   :: bmrot(3,3),gkrot(DIMENSION%nvd,3),vmult(3),v(3)
    COMPLEX:: ylm( (atoms%lmaxd+1)**2 ),chi(2,2)
    !     ..
    COMPLEX, ALLOCATABLE :: a(:,:,:),b(:,:,:),ax(:,:),bx(:,:)
    COMPLEX, ALLOCATABLE :: c_ph(:,:)
    COMPLEX,ALLOCATABLE :: aahlp(:),aa_tmphlp(:)
    INTEGER :: n_threads,thread,blocksize,maxloop
    INTEGER,ALLOCATABLE :: start_thread(:),stop_thread(:)


    
    lnonsphd=MAXVAL(atoms%lnonsph)*(MAXVAL(atoms%lnonsph)+2)
    ALLOCATE(dtd(0:lnonsphd,0:lnonsphd),utd(0:lnonsphd,0:lnonsphd),dtu(0:lnonsphd,0:lnonsphd),utu(0:lnonsphd,0:lnonsphd))
    !Decide how to distribute the work

    IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) ALLOCATE ( aahlp(hlpmsize),aa_tmphlp(hlpmsize) )

    ALLOCATE(aa_block(CPP_BLOCKSIZE,MAXVAL(lapw%nv)))

    ab_dim=1
    IF (noco%l_ss) ab_dim=2
    ALLOCATE(a(DIMENSION%nvd,0:DIMENSION%lmd,ab_dim),b(DIMENSION%nvd,0:DIMENSION%lmd,ab_dim))
    ALLOCATE(ax(DIMENSION%nvd,0:DIMENSION%lmd),bx(DIMENSION%nvd,0:DIMENSION%lmd))
    ALLOCATE(c_ph(DIMENSION%nvd,ab_dim))

    ntyploop: DO n=1,atoms%ntype
       IF (noco%l_noco) THEN
          IF (.NOT.noco%l_ss) aahlp=CMPLX(0.0,0.0)
          IF (.NOT.noco%l_ss) aa_tmphlp=CMPLX(0.0,0.0)
          CALL hsmt_spinor(isp,n, noco,input, chi, chi11, chi21, chi22)
       ENDIF
       DO nn = 1,atoms%neq(n)
	  a=0.0
	  b=0.0
          na = SUM(atoms%neq(:n-1))+nn
          IF (atoms%lnonsph(n)<0) CYCLE ntyploop
          IF ((atoms%invsat(na)==0) .OR. (atoms%invsat(na)==1)) THEN
             IF (atoms%invsat(na)==0) invsfct = 1
             IF (atoms%invsat(na)==1) invsfct = 2
             np = sym%invtab(atoms%ngopr(na))
             IF (oneD%odi%d1) np = oneD%ods%ngopr(na)
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
                   IF (oneD%odi%d1) THEN
                      bmrot=MATMUL(oneD%ods%mrot(:,:,np),cell%bmat)
                   ELSE
                      bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
                   END IF
                   DO k = 1,lapw%nv(iintsp)
                      !-->  apply the rotation that brings this atom into the
                      !-->  representative (this is the definition of ngopr(na)
                      !-->  and transform to cartesian coordinates
                      v(:) = vk(k,:,iintsp)
                      gkrot(k,:) = MATMUL(TRANSPOSE(bmrot),v)
                   END DO
                END IF
                DO k = 1,lapw%nv(iintsp)
                   !-->    generate spherical harmonics
                   vmult(:) =  gkrot(k,:)
                   CALL ylm4(atoms%lnonsph(n),vmult,ylm)
                   !-->  synthesize the complex conjugates of a and b
                   DO l = 0,atoms%lnonsph(n)
                      ll1 = l* (l+1)
                      DO m = -l,l
                         term = c_ph(k,iintsp)*ylm(ll1+m+1)
                         a(k,ll1+m,iintsp) = fj(k,l,n,spin2)*term
                         b(k,ll1+m,iintsp) = gj(k,l,n,spin2)*term
                      END DO
                   END DO
                ENDDO !k-loop
                !--->       end loop over interstitial spin
             ENDDO
             !--->       loops over the interstitial spin
             DO iintsp = 1,nintsp

                DO jintsp = 1,iintsp

                   jd = 1 ; IF (noco%l_noco) jd = isp
                   !--->       loop over l',m'
                   utu=0.0;utd=0.0;dtu=0.0;dtd=0.0
!!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(lp,mp,lmp,l,m,lm,in,utu,dtu,utd,dtd,im,k) &
!!$OMP SHARED(tlmplm,invsfct,lnonsph,nv,jintsp,jd,n)
                   DO lmp=0,atoms%lnonsph(n)*(atoms%lnonsph(n)+2)
                      lp=FLOOR(SQRT(1.0*lmp))
                      mp=lmp-lp*(lp+1)
                      IF (lp>atoms%lnonsph(n).OR.ABS(mp)>lp) STOP "BUG"
                      !--->             loop over l,m
                      DO l = 0,atoms%lnonsph(n)
                         DO m = -l,l
                            lm = l* (l+1) + m
                            in = tlmplm%ind(lmp,lm,n,jd)
                            IF (in/=-9999) THEN
                               IF (in>=0) THEN
                                  utu(lm,lmp) =CONJG(tlmplm%tuu(in,n,jd))*invsfct
                                  dtu(lm,lmp) =CONJG(tlmplm%tdu(in,n,jd))*invsfct
                                  utd(lm,lmp) =CONJG(tlmplm%tud(in,n,jd))*invsfct
                                  dtd(lm,lmp) =CONJG(tlmplm%tdd(in,n,jd))*invsfct
                               ELSE
                                  im = -in
                                  utu(lm,lmp) =tlmplm%tuu(im,n,jd)*invsfct
                                  dtu(lm,lmp) =tlmplm%tud(im,n,jd)*invsfct
                                  utd(lm,lmp) =tlmplm%tdu(im,n,jd)*invsfct
                                  dtd(lm,lmp) =tlmplm%tdd(im,n,jd)*invsfct
                               END IF
                               !--->    update ax, bx

                            END IF
                         END DO
                      END DO
                   ENDDO
!!$OMP END PARALLEL DO
                   lmp=atoms%lnonsph(n)*(atoms%lnonsph(n)+2)
                   !ax(:nv(jintsp),0:lmp)=(matmul(a(:nv(jintsp),0:lmp,jintsp),utu(0:lmp,0:lmp))+matmul(b(:nv(jintsp),0:lmp,jintsp),utd(0:lmp,0:lmp)))
                   !bx(:nv(jintsp),0:lmp)=(matmul(a(:nv(jintsp),0:lmp,jintsp),dtu(0:lmp,0:lmp))+matmul(b(:nv(jintsp),0:lmp,jintsp),dtd(0:lmp,0:lmp)))

                   CALL zgemm("N","N",lapw%nv(jintsp),lmp+1,lmp+1,one,a(1,0,jintsp),SIZE(a,1),utu(0,0),SIZE(utu,1),zero,ax,SIZE(ax,1))
                   CALL zgemm("N","N",lapw%nv(jintsp),lmp+1,lmp+1,one,b(1,0,jintsp),SIZE(a,1),utd(0,0),SIZE(utu,1),one,ax,SIZE(ax,1))

                   CALL zgemm("N","N",lapw%nv(jintsp),lmp+1,lmp+1,one,a(1,0,jintsp),SIZE(a,1),dtu(0,0),SIZE(utu,1),zero,bx,SIZE(ax,1))
                   CALL zgemm("N","N",lapw%nv(jintsp),lmp+1,lmp+1,one,b(1,0,jintsp),SIZE(a,1),dtd(0,0),SIZE(utu,1),one,bx,SIZE(ax,1))

                   !
                   !--->             update hamiltonian and overlap matrices
                   nc = 0
                   IF ( noco%l_noco .AND. (n_size>1) ) THEN
                      lapw%nv_tot = lapw%nv(1) + lapw%nv(2)
                   ELSE
                      lapw%nv_tot = lapw%nv(iintsp)
                   ENDIF
                   kii=n_rank
                   DO WHILE(kii<lapw%nv_tot)
                      !DO kii =  n_rank, nv_tot-1, n_size
                      ki = MOD(kii,lapw%nv(iintsp)) + 1
                      bsize=MIN(SIZE(aa_block,1),(lapw%nv(iintsp)-ki)/n_size+1) !Either use maximal blocksize or number of rows left to calculate
                      IF (bsize<1) EXIT !nothing more to do here
                      bsize2=bsize*n_size
                      bsize2=min(bsize2,lapw%nv(iintsp)-ki+1)
                      !aa_block(:bsize,:ki+bsize2-1)=matmul(a(ki:ki+bsize2-1:n_size,0:lmp,iintsp),conjg(transpose(ax(:ki+bsize2-1,0:lmp))))+ &
                      !                              matmul(b(ki:ki+bsize2-1:n_size,0:lmp,iintsp),conjg(transpose(bx(:ki+bsize2-1,0:lmp))))
                      IF (n_size==1) THEN !Make this a special case to avoid copy-in of a array
                         call zgemm("N","C",bsize,ki+bsize2-1,lmp+1,one,a(ki,0,iintsp),SIZE(a,1),ax(1,0),SIZE(ax,1),zero,aa_block,SIZE(aa_block,1))
                         call zgemm("N","C",bsize,ki+bsize2-1,lmp+1,one,b(ki,0,iintsp),SIZE(a,1),bx(1,0),SIZE(ax,1),one ,aa_block,SIZE(aa_block,1))
                      ELSE
                         CALL zgemm("N","C",bsize,ki+bsize2-1,lmp+1,one,a(ki:ki+bsize2-1:n_size,0:lmp,iintsp),SIZE(a(ki:ki+bsize2-1:n_size,0:lmp,iintsp),1),ax(1,0),SIZE(ax,1),zero,aa_block,SIZE(aa_block,1))
                         CALL zgemm("N","C",bsize,ki+bsize2-1,lmp+1,one,b(ki:ki+bsize2-1:n_size,0:lmp,iintsp),SIZE(a(ki:ki+bsize2-1:n_size,0:lmp,iintsp),1),bx(1,0),SIZE(ax,1),one,aa_block,SIZE(aa_block,1))
                      ENDIF
                      DO kb=1,bsize
                         IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
                            nc = 1+kii/n_size
                            ii = nc*(nc-1)/2*n_size-(nc-1)*(n_size-n_rank-1)
                            IF ( (n_size==1).OR.(kii+1<=lapw%nv(1)) ) THEN    !
                               aahlp(ii+1:ii+ki) = aahlp(ii+1:ii+ki)+MATMUL(CONJG(ax(:ki,:lmp)),a(ki,:,iintsp))+MATMUL(CONJG(bx(:ki,:lmp)),b(ki,:lmp,iintsp))
                            ELSE                    ! components for <2||2> block unused
                               aa_tmphlp(:ki) = MATMUL(CONJG(ax(:ki,:lmp)),a(ki,:lmp,iintsp))+MATMUL(CONJG(bx(:ki,:DIMENSION%lmd)),b(ki,:lmp,iintsp))
                               !--->                   spin-down spin-down part
                               ij = ii + lapw%nv(1)
                               aa_c(ij+1:ij+ki)=aa_c(ij+1:ij+ki)+chi22*aa_tmphlp(:ki)
                               !--->                   spin-down spin-up part, lower triangle
                               ij =  ii
                               aa_c(ij+1:ij+ki)=aa_c(ij+1:ij+ki)+chi21*aa_tmphlp(:ki)
                            ENDIF
                            !-||
                         ELSEIF ( noco%l_noco .AND. noco%l_ss ) THEN
                            IF ( iintsp==1 .AND. jintsp==1 ) THEN
                               !--->                      spin-up spin-up part
                               kjmax = ki
                               chihlp = chi11
                               ii = (ki-1)*(ki)/2
                            ELSEIF ( iintsp==2 .AND. jintsp==2 ) THEN
                               !--->                      spin-down spin-down part
                               kjmax = ki
                               chihlp = chi22
                               ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2+&
                                    lapw%nv(1)+atoms%nlotot
                            ELSE
                               !--->                      spin-down spin-up part
                               kjmax = lapw%nv(1)
                               chihlp = chi21
                               ii = (lapw%nv(1)+atoms%nlotot+ki-1)*(lapw%nv(1)+atoms%nlotot+ki)/2
                            ENDIF
                            aa_c(ii+1:ii+kjmax) = aa_c(ii+1:ii+kjmax) + chihlp*&
                                 (MATMUL(CONJG(ax(:kjmax,:lmp)),a(ki,:,iintsp))+MATMUL(CONJG(bx(:kjmax,:lmp)),b(ki,:lmp,iintsp)))
                         ELSE
                            nc = 1+kii/n_size
                            ii = nc*(nc-1)/2*n_size- (nc-1)*(n_size-n_rank-1)
                            if (l_real) THEN
                               aa_r(ii+1:ii+ki) = aa_r(ii+1:ii+ki) + aa_block(kb,:ki)
                            ELSE
                               aa_c(ii+1:ii+ki) = aa_c(ii+1:ii+ki) + aa_block(kb,:ki)
                            endif
                            !print*,ii,ki,kb
                            !                           IF (.not.apw(l)) THEN
                            !aa(ii+1:ii+ki) = aa(ii+1:ii+ki) + b(ki,lmp,iintsp)*bx(:ki)
                            !                           ENDIF
                         ENDIF
                         ki=ki+n_size
                         kii=kii+n_size
                      ENDDO
                      !--->             end loop over ki
                   END DO
                   !--->       end loops over interstitial spin
                ENDDO
             ENDDO
          ENDIF              ! atoms%invsat(na) = 0 or 1
          !--->    end loop over equivalent atoms
       END DO
       IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) CALL hsmt_hlptomat(atoms%nlotot,lapw%nv,sub_comm,chi11,chi21,chi22,aahlp,aa_c)
       !---> end loop over atom types (ntype)
    ENDDO ntyploop

    RETURN
  END SUBROUTINE hsmt_nonsph


END MODULE m_hsmt_nonsph
