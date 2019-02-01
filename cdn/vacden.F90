MODULE m_vacden
  USE m_juDFT
  !     *************************************************************
  !     determines the 2-d star function expansion coefficients of
  !     vacuum charge density. speed up by r. wu 1992
  !     *************************************************************
CONTAINS
  SUBROUTINE vacden(vacuum,DIMENSION,stars,oneD,kpts,input,sym,cell,atoms,noco,banddos,&
                    gVacMap,we,ikpt,jspin,vz,ne,lapw,evac,eig,den,zMat,dos)

    !***********************************************************************
    !     ****** change vacden(....,q) for vacuum density of states shz Jan.96
    !     ****** change vacden(......,dos%qstars) for starcoefficients, shz. Jan.99
    !     ****** changed for fleur dw
    !     In non-collinear calculations the density becomes a hermitian 2x2
    !     matrix. This subroutine generates this density matrix in the 
    !     vacuum region. The diagonal elements of this matrix (n_11 & n_22)
    !     are store in den%vacz and den%vacxy, while the real and imaginary part
    !     of the off-diagonal element are stored in den%vacz(:,:,3:4) and den%vacxy(:,:,:,3). 
    !
    !     Philipp Kurz 99/07
    !***********************************************************************

    !******** ABBREVIATIONS ************************************************
    !     qvac     : vacuum charge of each eigenstate, needed in in cdnval
    !                to determine the vacuum energy parameters
    !     vz       : non-warping part of the vacuum potential (matrix)
    !                collinear    : 2. index = ivac (# of vaccum)
    !                non-collinear: 2. index = ipot (comp. of pot. matr.)
    !     den%vacz : non-warping part of the vacuum density matrix,
    !                diagonal elements n_11 and n_22
    !     den%vacxy: warping part of the vacuum density matrix,
    !                diagonal elements n_11 and n_22
    !     den%vacz(:,:,3:4): non-warping part of the vacuum density matrix,
    !                off-diagonal elements n_21 (real part in (:,:,3), imaginary part in (:,:,4))
    !     den%vacxy(:,:,:,3): warping part of the vacuum density matrix,
    !                off-diagonal elements n_21
    !***********************************************************************
    !
    USE m_constants
    USE m_grdchlh
    USE m_qsf
    USE m_cylbes
    USE m_dcylbs
    USE m_od_abvac
    USE m_vacuz
    USE m_vacudz
    USE m_types
    IMPLICIT NONE
    TYPE(t_lapw),INTENT(INOUT)    :: lapw !for some reason the second spin data is reset in noco case
    TYPE(t_dimension),INTENT(IN)  :: DIMENSION
    TYPE(t_oneD),INTENT(IN)       :: oneD
    TYPE(t_banddos),INTENT(IN)    :: banddos
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_vacuum),INTENT(IN)     :: vacuum
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_kpts),INTENT(IN)       :: kpts
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_mat),INTENT(IN)        :: zMat
    TYPE(t_gVacMap),INTENT(IN)    :: gVacMap
    TYPE(t_potden),INTENT(INOUT)  :: den
    TYPE(t_dos),   INTENT(INOUT)  :: dos
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jspin      
    INTEGER, INTENT (IN) :: ne    
    INTEGER, INTENT (IN) :: ikpt
    INTEGER,PARAMETER    :: n2max=13
    REAL,PARAMETER        :: emax=2.0/hartree_to_ev_const
    !     .. Array Arguments ..
    REAL,    INTENT(IN)     :: evac(2,input%jspins)
    REAL,    INTENT(IN)     :: we(DIMENSION%neigd)
    REAL                    :: vz(vacuum%nmzd,2) ! Note this breaks the INTENT(IN) from cdnval. It may be read from a file in this subroutine.
    !     STM-Arguments
    REAL,    INTENT (IN)    :: eig(DIMENSION%neigd)
    !     local STM variables
    INTEGER nv2(input%jspins)
    INTEGER kvac1(DIMENSION%nv2d,input%jspins),kvac2(DIMENSION%nv2d,input%jspins),map2(DIMENSION%nvd,input%jspins)
    INTEGER kvac3(DIMENSION%nv2d,input%jspins),map1(DIMENSION%nvd,input%jspins)
    INTEGER mapg2k(DIMENSION%nv2d)
    !     .. Local Scalars ..
    COMPLEX aa,ab,av,ba,bb,bv,t1,aae,bbe,abe,bae,aaee,bbee,abee,baee,&
         &     factorx,factory,c_1,aa_1,ab_1,ba_1,bb_1,ic,av_1,bv_1,d,tempCmplx

    REAL arg,const,ddui,dduj,dduei,dduej,eps,ev,evacp,phs,phsp,qout,&
         &     scale,sign,uei,uej,ui,uj,wronk,zks,RESULT(1),ui2,uei2,&
         &     k_diff,k_d1,k_d2,ui_1,uei_1,uj_1,uej_1,wronk_1

    INTEGER i,ii,i1,i2,i3,ig3,ig3p,ik,ind2,ind2p,istar ,m1,m3,&
         &        ivac,j,jj,jz,k,l,ll,l1,n,n2,ispin,kspin,jsp_start,jsp_end,&
         &        ipot,ie,imz,isp_start,isp_end,&
         &        ind1,ind1p,irec2,irec3,m
    !
    !     .. Local Arrays ..
    REAL qssbti(3,2),qssbtii,vz0(2)
    REAL bess(-oneD%odi%mb:oneD%odi%mb),dbss(-oneD%odi%mb:oneD%odi%mb)
    COMPLEX, ALLOCATABLE :: ac(:,:,:),bc(:,:,:)
    REAL,    ALLOCATABLE :: dt(:),dte(:),du(:),ddu(:,:),due(:)
    REAL,    ALLOCATABLE :: ddue(:,:),t(:),te(:),tei(:,:),dummy(:)
    REAL,    ALLOCATABLE :: u(:,:,:),ue(:,:,:),v(:),yy(:)
    !-odim
    COMPLEX, ALLOCATABLE :: ac_1(:,:,:,:),bc_1(:,:,:,:)
    REAL,    ALLOCATABLE :: dt_1(:,:),dte_1(:,:)
    REAL,    ALLOCATABLE :: du_1(:,:),ddu_1(:,:,:),due_1(:,:)
    REAL,    ALLOCATABLE :: ddue_1(:,:,:),t_1(:,:)
    REAL,    ALLOCATABLE :: tei_1(:,:,:),te_1(:,:)
    REAL,    ALLOCATABLE :: u_1(:,:,:,:),ue_1(:,:,:,:)
    !+odim
    !     ..
    !     ..

    !     *******************************************************************************
    !

    !
    !    layers: no. of layers to be calculated (in vertical direction with z-values as given by izlay)
    !    izlay : defines vertical position of layers in delz (=0.1 a.u.) units from begining of vacuum region
    !    vacdos: =T: calculate vacuum dos in layers as given by the above
    !    integ : =T: integrate in vertical position between izlay(layer,1)..izlay(layer,2)
    !    nstm  : 0: s-Tip, 1: p_z-Tip, 2: d_z^2-Tip (following Chen's derivative rule) ->rhzgrd.f is used 
    !                 to calculate derivatives 
    !    tworkf: Workfunction of Tip (in hartree units) is needed for d_z^2-Orbital)
    !    starcoeff: =T: star coefficients are calculated at values of izlay for 0th (=q) to nstars-1. star 
    !                (dos%qstars(1..nstars-1))
    !    nstars: number of star functions to be used (0th star is given by value of q=charge integrated in 2D) 
    !
    !    further possibility: (readin of locx, locy has to be implemented in flapw7.f or they have to be set explicitly)
    !
    !     locx and locy can be used to calculate local DOS at a certain vertical position z (or integrated in z)
    !     within a restricted area of the 2D unit cell, the corners of this area is given by locx and locy
    !     they are defined in internal coordinates, i.e. \vec{r}_1=locx(1)*\vec{a}_1+locy(1)*\vec{a}_2
    !                                                    \vec{r}_2=locx(2)*\vec{a}_1+locy(2)*\vec{a}_2
    !                 \vec{a}_1,2 are the 2D lattice vectors           
    !  
    !     **************************************************************************************************

    CALL timestart("vacden")

    ALLOCATE ( ac(DIMENSION%nv2d,DIMENSION%neigd,input%jspins),bc(DIMENSION%nv2d,DIMENSION%neigd,input%jspins),dt(DIMENSION%nv2d),&
         &           dte(DIMENSION%nv2d),du(vacuum%nmzd),ddu(vacuum%nmzd,DIMENSION%nv2d),due(vacuum%nmzd),&
         &           ddue(vacuum%nmzd,DIMENSION%nv2d),t(DIMENSION%nv2d),te(DIMENSION%nv2d),&
         &           tei(DIMENSION%nv2d,input%jspins),u(vacuum%nmzd,DIMENSION%nv2d,input%jspins),ue(vacuum%nmzd,DIMENSION%nv2d,input%jspins),&
         &           v(3),yy(vacuum%nmzd))
    IF (oneD%odi%d1) THEN
       ALLOCATE (      ac_1(DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb,DIMENSION%neigd,input%jspins),&
            &                  bc_1(DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb,DIMENSION%neigd,input%jspins),&
            &                  dt_1(DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb),&
            &                 dte_1(DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb),&
            &                  du_1(vacuum%nmzd,-oneD%odi%mb:oneD%odi%mb),&
            &            ddu_1(vacuum%nmzd,DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb),&
            &                 due_1(vacuum%nmzd,-oneD%odi%mb:oneD%odi%mb),&
            &           ddue_1(vacuum%nmzd,DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb),&
            &                   t_1(DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb),&
            &                  te_1(DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb),&
            &                 tei_1(DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb,input%jspins),&
            &              u_1(vacuum%nmzd,DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb,input%jspins),&
            &             ue_1(vacuum%nmzd,DIMENSION%nv2d,-oneD%odi%mb:oneD%odi%mb,input%jspins) )
    END IF ! oneD%odi%d1
    !

    vz0(:) = vz(vacuum%nmz,:)
    eps=0.01
    ic = CMPLX(0.,1.)
    !    ------------------
   
    !     -----> set up mapping arrays
    IF (noco%l_ss) THEN
       jsp_start = 1
       jsp_end   = 2
    ELSE
       jsp_start = jspin
       jsp_end   = jspin
    ENDIF
    DO ispin = jsp_start,jsp_end
       IF (oneD%odi%d1) THEN
          n2 = 0
          k_loop:DO  k = 1,lapw%nv(ispin)
             DO  j = 1,n2
                IF (lapw%k3(k,ispin).EQ.kvac3(j,ispin)) THEN
                   map1(k,ispin) = j
                   CYCLE k_loop
                END IF
             ENDDO
             n2 = n2 + 1
             IF (n2>DIMENSION%nv2d)  CALL juDFT_error("vacden0",calledby ="vacden")
             kvac3(n2,ispin) =  lapw%k3(k,ispin)
             map1(k,ispin) = n2
          ENDDO k_loop
          nv2(ispin) = n2
       ELSE
          n2 = 0
          k_loop2:DO  k = 1,lapw%nv(ispin)
             DO  j = 1,n2
                IF ( lapw%gvec(1,k,ispin).EQ.kvac1(j,ispin) .AND.&
                     lapw%gvec(2,k,ispin).EQ.kvac2(j,ispin) ) THEN
                   map2(k,ispin) = j
                   CYCLE k_loop2
                END IF
             ENDDO
             n2 = n2 + 1
             IF (n2>DIMENSION%nv2d)  CALL juDFT_error("vacden0","vacden")
             kvac1(n2,ispin) = lapw%gvec(1,k,ispin)
             kvac2(n2,ispin) = lapw%gvec(2,k,ispin)
             map2(k,ispin) = n2
          ENDDO k_loop2
          nv2(ispin) = n2
       END IF
    ENDDO
    IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
       lapw%nv(2)  = lapw%nv(1)
       nv2(2) = nv2(1)
       DO k = 1,nv2(1)
          kvac1(k,2) = kvac1(k,1)
          kvac2(k,2) = kvac2(k,1)
          IF(oneD%odi%d1) kvac3(k,2) =  kvac3(k,1)
       ENDDO
       DO k = 1,lapw%nv(1)
          lapw%k3(k,2) = lapw%k3(k,1)
          map2(k,2) = map2(k,1)
          IF(oneD%odi%d1)THEN
             lapw%k1(k,2) = lapw%k1(k,1)
             lapw%k2(k,2) = lapw%k2(k,1)
             map1(k,2) = map1(k,1)
          ENDIF
       ENDDO
    ENDIF

    !+dw
    !    if tunneling current should be calculated we need to write out 
    !     info on electronic structure: --> mapping from kvac to gvac by mapg2k
    !                                             shz, Jan.99
    IF (vacuum%nstm.EQ.3) THEN
       DO j=1, n2max
          mapg2k(j)=j
          DO i=1, nv2(jspin)
             IF (kvac1(i,jspin).EQ.gVacMap%gvac1d(j).AND.kvac2(i,jspin).EQ.gVacMap%gvac2d(j)) mapg2k(j)=i
          END DO
       END DO
    END IF
    !
    !-dw
   

    wronk = 2.0
    const = 1.0 / ( SQRT(cell%omtil)*wronk )
    !-odim
    IF (oneD%odi%d1) THEN
       ac_1(:,:,:,:) = CMPLX(0.0,0.0)
       bc_1(:,:,:,:) = CMPLX(0.0,0.0)
    END IF
    !+odim
    DO  ivac = 1,vacuum%nvac
       ac(:,:,:) = CMPLX(0.0,0.0)
       bc(:,:,:) = CMPLX(0.0,0.0)
       sign = 3. - 2.*ivac
      
       IF (noco%l_noco) THEN
          !--->    In a non-collinear calculation vacden is only called once.
          !--->    Thus, the vaccum wavefunctions and the A- and B-coeff. (ac bc)
          !--->    have to be calculated for both spins on that call.
          !--->       setup the spin-spiral q-vector
          qssbti(1,1) = - noco%qss(1)/2
          qssbti(2,1) = - noco%qss(2)/2
          qssbti(1,2) = + noco%qss(1)/2
          qssbti(2,2) = + noco%qss(2)/2
          qssbti(3,1) = - noco%qss(3)/2
          qssbti(3,2) = + noco%qss(3)/2
          DO ispin = 1,input%jspins
             !     -----> set up vacuum wave functions
             IF (oneD%odi%d1) THEN
                CALL od_abvac(&
                     cell,vacuum,DIMENSION,stars,&
                     oneD,qssbti(3,ispin),&
                     oneD%odi%n2d,&
                     wronk,evacp,lapw%bkpt,oneD%odi%M,oneD%odi%mb,&
                     vz(1,ispin),kvac3(1,ispin),nv2(ispin),&
                     t_1(1,-oneD%odi%mb),dt_1(1,-oneD%odi%mb),u_1(1,1,-oneD%odi%mb,ispin),&
                     te_1(1,-oneD%odi%mb),dte_1(1,-oneD%odi%mb),&
                     tei_1(1,-oneD%odi%mb,ispin),&
                     ue_1(1,1,-oneD%odi%mb,ispin))
                DO k = 1,lapw%nv(ispin)
                   kspin = (lapw%nv(1)+atoms%nlotot)*(ispin-1) + k
                   l = map1(k,ispin) 
                   irec3 = stars%ig(lapw%gvec(1,k,ispin),lapw%gvec(2,k,ispin),lapw%gvec(3,k,ispin))
                   IF (irec3.NE.0) THEN
                      irec2 = stars%ig2(irec3)
                      zks = stars%sk2(irec2)*cell%z1
                      arg = stars%phi2(irec2)
                      CALL cylbes(oneD%odi%mb,zks,bess)
                      CALL dcylbs(oneD%odi%mb,zks,bess,dbss)
                      DO m = -oneD%odi%mb,oneD%odi%mb
                         wronk_1 = t_1(l,m)*dte_1(l,m) -&
                              te_1(l,m)*dt_1(l,m)
                         av_1 = EXP(-CMPLX(0.0,m*arg))*(ic**m)*&
                              CMPLX(dte_1(l,m)*bess(m) -&
                              te_1(l,m)*stars%sk2(irec2)*dbss(m),0.0)/&
                              ((wronk_1)*SQRT(cell%omtil))
                         bv_1 = EXP(-CMPLX(0.0,m*arg))*(ic**m)*&
                              CMPLX(-dt_1(l,m)*bess(m) +&
                              t_1(l,m)*stars%sk2(irec2)*dbss(m),0.0)/&
                              ((wronk_1)*SQRT(cell%omtil))
                         IF (zmat%l_real) THEN
                            ac_1(l,m,:ne,ispin) = ac_1(l,m,:ne,ispin) + zMat%data_r(kspin,:ne)*av_1
                            bc_1(l,m,:ne,ispin) = bc_1(l,m,:ne,ispin) + zMat%data_r(kspin,:ne)*bv_1
                         ELSE
                            ac_1(l,m,:ne,ispin) = ac_1(l,m,:ne,ispin) + zMat%data_c(kspin,:ne)*av_1
                            bc_1(l,m,:ne,ispin) = bc_1(l,m,:ne,ispin) + zMat%data_c(kspin,:ne)*bv_1
                         END IF
                      END DO      ! -mb:mb
                   END IF
                END DO
             ELSE ! 1-dimensional
                vz0(ispin) = vz(vacuum%nmz,ispin)
                evacp = evac(ivac,ispin)
                DO ik = 1,nv2(ispin)
                   v(1) = lapw%bkpt(1) + kvac1(ik,ispin) + qssbti(1,ispin)
                   v(2) = lapw%bkpt(2) + kvac2(ik,ispin) + qssbti(2,ispin)
                   v(3) = 0.
                   ev = evacp - 0.5*DOT_PRODUCT(v,MATMUL(v,cell%bbmat))
                   CALL vacuz(ev,vz(1,ispin),vz0(ispin),vacuum%nmz,vacuum%delz,t(ik),&
                        dt(ik),u(1,ik,ispin))
                   CALL vacudz(ev,vz(1,ispin),vz0(ispin),vacuum%nmz,vacuum%delz,te(ik),&
                        dte(ik),tei(ik,ispin),ue(1,ik,ispin),dt(ik),&
                        u(1,ik,ispin))
                   scale = wronk/ (te(ik)*dt(ik)-dte(ik)*t(ik))
                   te(ik) = scale*te(ik)
                   dte(ik) = scale*dte(ik)
                   tei(ik,ispin) = scale*tei(ik,ispin)
                   DO j = 1,vacuum%nmz
                      ue(j,ik,ispin) = scale*ue(j,ik,ispin)
                   ENDDO
                ENDDO
                !     -----> construct a and b coefficients
                DO k = 1,lapw%nv(ispin)
                   !--->          the coefficients of the spin-down basis functions are
                   !--->          stored in the second half of the eigenvector
                   kspin = (lapw%nv(1)+atoms%nlotot)*(ispin-1) + k
                   l = map2(k,ispin)
                   zks = lapw%k3(k,ispin)*cell%bmat(3,3)*sign
                   arg = zks*cell%z1
                   c_1 = CMPLX(COS(arg),SIN(arg)) * const
                   av = -c_1 * CMPLX( dte(l),zks*te(l) ) 
                   bv =  c_1 * CMPLX(  dt(l),zks* t(l) ) 
                   !     -----> loop over basis functions
                   IF (zmat%l_real) THEN
                      ac(l,:ne,ispin) = ac(l,:ne,ispin) + zMat%data_r(kspin,:ne)*av
                      bc(l,:ne,ispin) = bc(l,:ne,ispin) + zMat%data_r(kspin,:ne)*bv
                   ELSE
                      ac(l,:ne,ispin) = ac(l,:ne,ispin) + zMat%data_c(kspin,:ne)*av
                      bc(l,:ne,ispin) = bc(l,:ne,ispin) + zMat%data_c(kspin,:ne)*bv
                   ENDIF
                ENDDO
                !--->       end of spin loop
             ENDIF
             !---Y       end of geometry 1d/film
          ENDDO
          !--->       output for testing
          !            DO k = 1,10
          !               DO n = 1,5
          !                  DO ispin = 1,jspins
          !                     write(*,9000)k,n,ispin,ac(k,n,ispin),bc(k,n,ispin)
          !                  ENDDO
          !               ENDDO
          !            ENDDO
          ! 9000       FORMAT('k=',i3,' ie=',i3,' isp=',i3,
          !     +             ' ac= (',e12.6,',',e12.6,')',
          !     +             ' bc= (',e12.6,',',e12.6,')')
       ELSE
          !     -----> set up vacuum wave functions
          IF (oneD%odi%d1) THEN
             qssbtii = 0.
             evacp = evac(ivac,jspin)
             CALL od_abvac(&
                  &           cell,vacuum,DIMENSION,stars,&
                  &           oneD,qssbtii,&
                  &           oneD%odi%n2d,&
                  &           wronk,evacp,lapw%bkpt,oneD%odi%M,oneD%odi%mb,&
                  &           vz(1,ivac),kvac3(1,jspin),nv2(jspin),&
                  &           t_1(1,-oneD%odi%mb),dt_1(1,-oneD%odi%mb),u_1(1,1,-oneD%odi%mb,jspin),&
                  &           te_1(1,-oneD%odi%mb),dte_1(1,-oneD%odi%mb),&
                  &           tei_1(1,-oneD%odi%mb,jspin),&
                  &           ue_1(1,1,-oneD%odi%mb,jspin))
             DO k = 1,lapw%nv(jspin)
                l = map1(k,jspin)
                irec3 = stars%ig(lapw%gvec(1,k,jspin),lapw%gvec(2,k,jspin),lapw%gvec(3,k,jspin))
                IF (irec3.NE.0) THEN
                   irec2 = stars%ig2(irec3)
                   zks = stars%sk2(irec2)*cell%z1
                   arg = stars%phi2(irec2)
                   CALL cylbes(oneD%odi%mb,zks,bess)
                   CALL dcylbs(oneD%odi%mb,zks,bess,dbss)
                   DO m = -oneD%odi%mb,oneD%odi%mb
                      wronk_1 = t_1(l,m)*dte_1(l,m) -te_1(l,m)*dt_1(l,m)
                      av_1 = EXP(-CMPLX(0.0,m*arg))*(ic**m)*&
                           CMPLX(dte_1(l,m)*bess(m) -&
                           te_1(l,m)*stars%sk2(irec2)*dbss(m),0.0)/&
                           ((wronk_1)*SQRT(cell%omtil))
                      bv_1 = EXP(-CMPLX(0.0,m*arg))*(ic**m)*&
                           CMPLX(-dt_1(l,m)*bess(m) +&
                           t_1(l,m)*stars%sk2(irec2)*dbss(m),0.0)/&
                           ((wronk_1)*SQRT(cell%omtil))
                      IF (zmat%l_real) THEN
                         ac_1(l,m,:ne,jspin) = ac_1(l,m,:ne,jspin) + zMat%data_r(k,:ne)*av_1
                         bc_1(l,m,:ne,jspin) = bc_1(l,m,:ne,jspin) + zMat%data_r(k,:ne)*bv_1
                      ELSE
                         ac_1(l,m,:ne,jspin) = ac_1(l,m,:ne,jspin) + zMat%data_c(k,:ne)*av_1
                         bc_1(l,m,:ne,jspin) = bc_1(l,m,:ne,jspin) + zMat%data_c(k,:ne)*bv_1
                      ENDIF
                   END DO      ! -mb:mb
                END IF
             END DO         ! k = 1,lapw%nv
          ELSE     !oneD%odi%d1
             evacp = evac(ivac,jspin)
             DO ik = 1,nv2(jspin)
                v(1) = lapw%bkpt(1) + kvac1(ik,jspin)
                v(2) = lapw%bkpt(2) + kvac2(ik,jspin)
                v(3) = 0.
                ev = evacp - 0.5*DOT_PRODUCT(v,MATMUL(v,cell%bbmat))
                CALL vacuz(ev,vz(1,ivac),vz0(ivac),vacuum%nmz,vacuum%delz,t(ik),dt(ik),u(1,ik,jspin))
                CALL vacudz(ev,vz(1,ivac),vz0(ivac),vacuum%nmz,vacuum%delz,te(ik),&
                     &              dte(ik),tei(ik,jspin),ue(1,ik,jspin),dt(ik),&
                     &              u(1,ik,jspin))
                scale = wronk/ (te(ik)*dt(ik)-dte(ik)*t(ik))
                te(ik) = scale*te(ik)
                dte(ik) = scale*dte(ik)
                tei(ik,jspin) = scale*tei(ik,jspin)
                DO j = 1,vacuum%nmz
                   ue(j,ik,jspin) = scale*ue(j,ik,jspin)
                ENDDO
             ENDDO
             !     -----> construct a and b coefficients
             DO k = 1,lapw%nv(jspin)
                l = map2(k,jspin)
                zks = lapw%k3(k,jspin)*cell%bmat(3,3)*sign
                arg = zks*cell%z1
                c_1 = CMPLX(COS(arg),SIN(arg)) * const
                av = -c_1 * CMPLX( dte(l),zks*te(l) ) 
                bv =  c_1 * CMPLX(  dt(l),zks* t(l) ) 
                !     -----> loop over basis functions
                IF (zmat%l_real) THEN
                   ac(l,:ne,jspin) = ac(l,:ne,jspin) + zMat%data_r(k,:ne)*av
                   bc(l,:ne,jspin) = bc(l,:ne,jspin) + zMat%data_r(k,:ne)*bv
                ELSE
                   ac(l,:ne,jspin) = ac(l,:ne,jspin) + zMat%data_c(k,:ne)*av
                   bc(l,:ne,jspin) = bc(l,:ne,jspin) + zMat%data_c(k,:ne)*bv
                ENDIF
             ENDDO
          END IF ! D1
       ENDIF
       !
       !   ----> calculate first and second derivative of u,ue 
       !        in order to simulate p_z or d_z^2 Tip in Chen's model , shz. 97
       ! 
       IF (vacuum%nstm.GT.0) THEN
          DO  ik = 1,nv2(jspin)
             !               CALL rhzgrd(nmz,delz,u(1,ik,jspin),4,du,ddu(1,ik))
             !               CALL rhzgrd(nmz,delz,ue(1,ik,jspin),4,due,ddue(1,ik))   

             ALLOCATE ( dummy(vacuum%nmz) )
             CALL grdchlh(&
                  0,1,vacuum%nmz,vacuum%delz,dummy,u(1,ik,jspin),4,&
                  du,ddu(1,ik))
             CALL grdchlh(&
                  0,1,vacuum%nmz,vacuum%delz,dummy,ue(1,ik,jspin),4,&
                  due,ddue(1,ik))
             DEALLOCATE ( dummy )

             IF (vacuum%nstm.EQ.1) THEN
                u(:vacuum%nmz,ik,jspin)=du(:vacuum%nmz)
                ue(:vacuum%nmz,ik,jspin)=due(:vacuum%nmz)
             END IF
          ENDDO
       END IF

       !+dw

       !
       !       --> to calculate Tunneling Current between two systems
       !           within Bardeens Approach one needs ac(l,n), bc(l,n);
       !           they are written to the file vacwave 
       !                           IF nstm=3
       !                              tworkf is then the fermi energy (in hartree)
       !
       IF (vacuum%nstm.EQ.3) THEN
#ifdef CPP_MPI
          CALL judft_error("nstm==3 does not work in parallel",calledby="vacden")
#else
          i=0
          DO n = 1, ne
             IF (ABS(eig(n)-vacuum%tworkf).LE.emax) i=i+1
          END DO
          WRITE (87,FMT=990) lapw%bkpt(1),lapw%bkpt(2), i, n2max
          DO n = 1, ne
             IF (ABS(eig(n)-vacuum%tworkf).LE.emax) THEN
                WRITE (87,FMT=1000) eig(n)
                DO j=1,n2max
                   WRITE (87,FMT=1010) ac(mapg2k(j),n,jspin),&
                        bc(mapg2k(j),n,jspin)
                END DO
             END IF
          END DO
#endif



       END IF
990    FORMAT(2(f8.4,1x),i3,1x,i3)
1000   FORMAT(e10.4)
1010   FORMAT(2(2e20.8,1x))
       !
       !        ------------------------------------------------------------
       !-dw
       !
       !---->   non-warping part of the density (g=0 star)
       !
       IF (vacuum%nstm.EQ.2) THEN
          !
          !  ----> d_z^2-Tip needs: |d^2(psi)/dz^2 - kappa^2/3 psi|^2
          !
          DO l = 1,nv2(jspin)
             aa = 0.0
             bb = 0.0
             ba = 0.0
             ab = 0.0
             DO n = 1,ne
                aa = aa + we(n)*CONJG(ac(l,n,jspin))*ac(l,n,jspin)
                bb = bb + we(n)*CONJG(bc(l,n,jspin))*bc(l,n,jspin)
                ab = ab + we(n)*CONJG(ac(l,n,jspin))*bc(l,n,jspin)
                ba = ba + we(n)*CONJG(bc(l,n,jspin))*ac(l,n,jspin)
                qout = REAL(CONJG(ac(l,n,jspin))*ac(l,n,jspin)+tei(l,jspin)*CONJG(bc(l,n,jspin))*bc(l,n,jspin))
                dos%qvac(n,ivac,ikpt,jspin) = dos%qvac(n,ivac,ikpt,jspin) + qout*cell%area
             END DO
             aae=-aa*vacuum%tworkf*2/3
             bbe=-bb*vacuum%tworkf*2/3
             abe=-ab*vacuum%tworkf*2/3
             bae=-ba*vacuum%tworkf*2/3
             aaee=aa*vacuum%tworkf*vacuum%tworkf*4/9
             bbee=bb*vacuum%tworkf*vacuum%tworkf*4/9
             abee=ab*vacuum%tworkf*vacuum%tworkf*4/9
             baee=ba*vacuum%tworkf*vacuum%tworkf*4/9  
             DO  jz = 1,vacuum%nmz
                ui = u(jz,l,jspin)
                uei = ue(jz,l,jspin)
                ddui = ddu(jz,l) 
                dduei = ddue(jz,l)
                den%vacz(jz,ivac,jspin) = den%vacz(jz,ivac,jspin) +&
                     REAL(aaee*ui*ui+bbee*uei*uei+&
                     (abee+baee)*ui*uei+aa*ddui*ddui+&
                     bb*dduei*dduei+(ab+ba)*ddui*dduei+&
                     2*aae*ui*ddui+2*bbe*uei*dduei+&
                     (abe+bae)*(ui*dduei+uei*ddui))

             ENDDO
          END DO
          !
          !    -----> s-Tip: |psi|^2 and p-Tip: |d(psi)/dz|^2 
          !
       ELSE
          IF (noco%l_noco) THEN
             !--->          diagonal elements of the density matrix, n_11 and n_22
             !--->          the non-warping part of n_21 is calculated together with
             !--->          the warping part of n_21
             DO ispin = 1,input%jspins
                IF (oneD%odi%d1) THEN
                   DO  l = 1,nv2(ispin)
                      DO  m = -oneD%odi%mb,oneD%odi%mb
                         aa = CMPLX(0.0,0.0)
                         bb = CMPLX(0.0,0.0)
                         ba = CMPLX(0.0,0.0)
                         ab = CMPLX(0.0,0.0)
                         DO  n = 1,ne
                            aa = aa + we(n)*CONJG(ac_1(l,m,n,ispin))*ac_1(l,m,n,ispin)
                            bb = bb + we(n)*CONJG(bc_1(l,m,n,ispin))*bc_1(l,m,n,ispin)
                            ab = ab + we(n)*CONJG(ac_1(l,m,n,ispin))*bc_1(l,m,n,ispin)
                            ba = ba + we(n)*CONJG(bc_1(l,m,n,ispin))*ac_1(l,m,n,ispin)
                            qout = REAL(CONJG(ac_1(l,m,n,ispin))*ac_1(l,m,n,ispin) +&
                                 tei_1(l,m,ispin)*CONJG(bc_1(l,m,n,ispin))*&
                                 bc_1(l,m,n,ispin))
                            dos%qvac(n,ivac,ikpt,ispin) = dos%qvac(n,ivac,ikpt,ispin)+qout*cell%area
                         END DO
                         DO  jz = 1,vacuum%nmz
                            ui = u_1(jz,l,m,ispin)
                            uei = ue_1(jz,l,m,ispin)
                            den%vacz(jz,ivac,ispin) = den%vacz(jz,ivac,ispin) +REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
                         END DO
                      END DO
                   END DO
                ELSE 
                   DO l = 1,nv2(ispin)
                      aa = 0.0
                      bb = 0.0
                      ba = 0.0
                      ab = 0.0
                      DO n = 1,ne
                         aa=aa + we(n)*CONJG(ac(l,n,ispin))*ac(l,n,ispin)
                         bb=bb + we(n)*CONJG(bc(l,n,ispin))*bc(l,n,ispin)
                         ab=ab + we(n)*CONJG(ac(l,n,ispin))*bc(l,n,ispin)
                         ba=ba + we(n)*CONJG(bc(l,n,ispin))*ac(l,n,ispin)
                         qout = REAL(CONJG(ac(l,n,ispin))*ac(l,n,ispin)+tei(l,ispin)*CONJG(bc(l,n,ispin))*bc(l,n,ispin))
                         dos%qvac(n,ivac,ikpt,ispin) = dos%qvac(n,ivac,ikpt,ispin) + qout*cell%area
                      END DO
                      DO jz = 1,vacuum%nmz
                         ui = u(jz,l,ispin)
                         uei = ue(jz,l,ispin)
                         den%vacz(jz,ivac,ispin) = den%vacz(jz,ivac,ispin) +REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
                      ENDDO
                   ENDDO
                END IF ! one-dimensional
             ENDDO
          ELSE
             IF (oneD%odi%d1) THEN
                DO  l = 1,nv2(jspin)
                   DO  m = -oneD%odi%mb,oneD%odi%mb
                      aa = CMPLX(0.0,0.0)
                      bb = CMPLX(0.0,0.0)
                      ba = CMPLX(0.0,0.0)
                      ab = CMPLX(0.0,0.0)
                      DO  n = 1,ne
                         aa = aa + we(n)*CONJG(ac_1(l,m,n,jspin))*ac_1(l,m,n,jspin)
                         bb = bb + we(n)*CONJG(bc_1(l,m,n,jspin))*bc_1(l,m,n,jspin)
                         ab = ab + we(n)*CONJG(ac_1(l,m,n,jspin))*bc_1(l,m,n,jspin)
                         ba = ba + we(n)*CONJG(bc_1(l,m,n,jspin))*ac_1(l,m,n,jspin)
                         qout = REAL(CONJG(ac_1(l,m,n,jspin))*ac_1(l,m,n,jspin)+&
                              tei_1(l,m,jspin)*CONJG(bc_1(l,m,n,jspin))*bc_1(l,m,n,jspin))
                         dos%qvac(n,ivac,ikpt,jspin) = dos%qvac(n,ivac,ikpt,jspin)+qout*cell%area
                      END DO
                      DO  jz = 1,vacuum%nmz
                         ui = u_1(jz,l,m,jspin)
                         uei = ue_1(jz,l,m,jspin)
                         den%vacz(jz,ivac,jspin) = den%vacz(jz,ivac,jspin) +REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
                      END DO
                   END DO
                END DO
             ELSE          ! D1
                DO l = 1,nv2(jspin)
                   aa = CMPLX(0.0,0.0)
                   bb = CMPLX(0.0,0.0)
                   ba = CMPLX(0.0,0.0)
                   ab = CMPLX(0.0,0.0)
                   DO n = 1,ne
                      aa = aa + we(n)*CONJG(ac(l,n,jspin))*ac(l,n,jspin)
                      bb = bb + we(n)*CONJG(bc(l,n,jspin))*bc(l,n,jspin)
                      ab = ab + we(n)*CONJG(ac(l,n,jspin))*bc(l,n,jspin)
                      ba = ba + we(n)*CONJG(bc(l,n,jspin))*ac(l,n,jspin)
                      qout = REAL(CONJG(ac(l,n,jspin))*ac(l,n,jspin)+tei(l,jspin)*CONJG(bc(l,n,jspin))*bc(l,n,jspin))
                      dos%qvac(n,ivac,ikpt,jspin) = dos%qvac(n,ivac,ikpt,jspin) + qout*cell%area
                   END DO
                   DO  jz = 1,vacuum%nmz
                      ui = u(jz,l,jspin)
                      uei = ue(jz,l,jspin)
                      den%vacz(jz,ivac,jspin) = den%vacz(jz,ivac,jspin) +REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
                   ENDDO
                END DO
             END IF ! D1
          ENDIF
       END IF
       !
       !     ****************** change for vacuum density of states shz Jan.96 ***
       !
       IF (banddos%vacdos) THEN                       
          !
          !  ----> d_z^2-Tip needs: |d^2(psi)/dz^2 - kappa^2/3 psi|^2
          !
          IF (vacuum%nstm.EQ.2) THEN
             DO l=1,nv2(jspin)
                DO n = 1,ne
                   aa = CONJG(ac(l,n,jspin))*ac(l,n,jspin)
                   bb = CONJG(bc(l,n,jspin))*bc(l,n,jspin)
                   ab = CONJG(ac(l,n,jspin))*bc(l,n,jspin)
                   ba = CONJG(bc(l,n,jspin))*ac(l,n,jspin)
                   aae = -vacuum%tworkf*aa*2/3
                   bbe = -vacuum%tworkf*bb*2/3
                   abe = -vacuum%tworkf*ab*2/3
                   bae = -vacuum%tworkf*ba*2/3
                   aaee = aa*vacuum%tworkf*vacuum%tworkf*4/9
                   bbee = bb*vacuum%tworkf*vacuum%tworkf*4/9
                   abee = ab*vacuum%tworkf*vacuum%tworkf*4/9
                   baee = ba*vacuum%tworkf*vacuum%tworkf*4/9
                   DO jj = 1,vacuum%layers
                      !
                      !     ----> either integrated LDOS(z1,z2) or LDOS(z1)
                      !     
                      IF (input%integ) THEN
                         ll = 1
                         DO ii = vacuum%izlay(jj,1),vacuum%izlay(jj,2)
                            ui = u(ii,l,jspin)
                            uei = ue(ii,l,jspin)
                            ddui = ddu(ii,l)
                            dduei = ddue(ii,l)
                            yy(ll) = REAL(aaee*ui*ui+bbee*uei*uei+(abee+baee)*ui*uei+aa*ddui*ddui+bb*&
                                 dduei*dduei+(ab+ba)*ddui*dduei+2*aae*ui*ddui+2*bbe*uei*dduei+&
                                 (abe+bae)*(ui*dduei+uei*ddui))*cell%area
                            ll = ll+1
                         END DO
                         CALL qsf(vacuum%delz,yy,RESULT,ll-1,0)
                         dos%qvlay(n,jj,ivac,ikpt,jspin) = dos%qvlay(n,jj,ivac,ikpt,jspin) + RESULT(1)
                      ELSE
                         ui = u(vacuum%izlay(jj,1),l,jspin)       
                         uei = ue(vacuum%izlay(jj,1),l,jspin)
                         ddui = ddu(vacuum%izlay(jj,1),l)
                         dduei = ddue(vacuum%izlay(jj,1),l)
                         yy(1) = REAL(aaee*ui*ui+bbee*uei*uei+&
                              (abee+baee)*ui*uei+aa*ddui*ddui+&
                              bb*dduei*dduei+(ab+ba)*ddui*dduei+&
                              2*aae*ui*ddui+2*bbe*uei*dduei+&
                              (abe+bae)*(ui*dduei+uei*ddui))
                         dos%qvlay(n,jj,ivac,ikpt,jspin) = dos%qvlay(n,jj,ivac,ikpt,jspin) +yy (1)
                      END IF
                   END DO
                END DO
             END DO
             !     
             !     ----> s-Tip = calculate LDOS and(!) p_z-Tip (since u->du/dz, ue->due/dz)
             !     
          ELSE 
             IF (ABS(vacuum%locx(1)-vacuum%locx(2)).LE.eps) THEN
                !     
                !     ----> integrated over 2D-unit cell
                !
                IF (noco%l_noco) THEN
                   isp_start = 1
                   isp_end   = input%jspins
                ELSE
                   isp_start = jspin
                   isp_end   = jspin
                ENDIF
                DO ispin = isp_start, isp_end
                   DO l=1,nv2(ispin)
                      DO n = 1,ne
                         aa = CONJG(ac(l,n,ispin))*ac(l,n,ispin)
                         bb = CONJG(bc(l,n,ispin))*bc(l,n,ispin)
                         ab = CONJG(ac(l,n,ispin))*bc(l,n,ispin)
                         ba = CONJG(bc(l,n,ispin))*ac(l,n,ispin) 
                         DO jj = 1,vacuum%layers
                            !     
                            !     ---> either integrated (z1,z2) or slice (z1)
                            !     
                            IF (input%integ) THEN
                               ll = 1
                               DO ii = vacuum%izlay(jj,1),vacuum%izlay(jj,2)
                                  ui = u(ii,l,ispin)
                                  uei = ue(ii,l,ispin)
                                  yy(ll) = REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
                                  ll = ll+1
                               END DO
                               CALL qsf(vacuum%delz,yy,RESULT,ll-1,0)
                               dos%qvlay(n,jj,ivac,ikpt,ispin) = dos%qvlay(n,jj,ivac,ikpt,ispin) + RESULT(1)
                            ELSE
                               ui = u(vacuum%izlay(jj,1),l,ispin)       
                               uei = ue(vacuum%izlay(jj,1),l,ispin)
                               dos%qvlay(n,jj,ivac,ikpt,ispin) = dos%qvlay(n,jj,ivac,ikpt,ispin) + REAL(&
                                    aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)

                            END IF
                         END DO
                      END DO
                   END DO
                ENDDO
             ELSE
                !     
                !     ----> if LDOS should be calculated over restricted area of the 2D-unit cell
                !     lower left corner: (locx(1), locy(1))   }  in internal
                !     upper right corner: (locx(2), locy(2))  }  coordinates
                !     
                DO l=1, nv2(jspin)
                   DO l1=1, nv2(jspin)
                      IF (kvac1(l,jspin).EQ.kvac1(l1,jspin)) THEN
                         factorx = CMPLX((vacuum%locx(2)-vacuum%locx(1)), 0.)
                      ELSE
                         k_diff=tpi_const*(kvac1(l,jspin)-kvac1(l1,jspin))
                         k_d1 = k_diff*vacuum%locx(1)
                         k_d2 = k_diff*vacuum%locx(2)
                         factorx=( CMPLX( COS(k_d2), SIN(k_d2)) -&
                              CMPLX( COS(k_d1), SIN(k_d1)) ) /&
                              CMPLX( 0.,k_diff )
                      END IF
                      IF (kvac2(l,jspin).EQ.kvac2(l1,jspin)) THEN
                         factory = CMPLX((vacuum%locy(2)-vacuum%locy(1)), 0.)
                      ELSE
                         k_diff=tpi_const*(kvac2(l,jspin)-kvac2(l1,jspin))
                         k_d1 = k_diff*vacuum%locy(1)
                         k_d2 = k_diff*vacuum%locy(2)
                         factory=( CMPLX( COS(k_d2), SIN(k_d2)) -&
                              CMPLX( COS(k_d1), SIN(k_d1)) ) /&
                              CMPLX( 0.,k_diff )
                      END IF
                      DO n=1, ne
                         aa = CONJG(ac(l1,n,jspin))*ac(l,n,jspin)
                         bb = CONJG(bc(l1,n,jspin))*bc(l,n,jspin)
                         ab = CONJG(ac(l1,n,jspin))*bc(l,n,jspin)
                         ba = CONJG(bc(l1,n,jspin))*ac(l,n,jspin)   
                         DO jj = 1,vacuum%layers
                            !     
                            !     ---> either integrated (z1,z2) or slice (z1)
                            !     
                            IF (input%integ) THEN
                               ll = 1
                               DO ii = vacuum%izlay(jj,1), vacuum%izlay(jj,2)
                                  ui = u(ii,l,jspin)
                                  uei = ue(ii,l,jspin)
                                  uj = u(ii,l1,jspin)
                                  uej = ue(ii,l1,jspin)
                                  yy(ll) = REAL((aa*ui*uj+bb*uei*uej+ab*uei*uj+ba*ui*uej)*factorx*factory)
                                  ll = ll+1
                               END DO
                               CALL qsf(vacuum%delz,yy,RESULT,ll-1,0)
                               dos%qvlay(n,jj,ivac,ikpt,jspin) = dos%qvlay(n,jj,ivac,ikpt,jspin) + RESULT(1)
                            ELSE
                               ui = u(vacuum%izlay(jj,1),l,jspin)       
                               uei = ue(vacuum%izlay(jj,1),l,jspin)
                               uj = u(vacuum%izlay(jj,1),l1,jspin)
                               uej = ue(vacuum%izlay(jj,1),l1,jspin)
                               dos%qvlay(n,jj,ivac,ikpt,jspin) = REAL((aa*ui*uj + bb*uei*uej+ab*uei*uj+ba*ui**uej)*factorx*factory)
                            END IF
                         END DO
                      END DO
                   END DO
                END DO
             END IF
          END IF
       END IF

       !
       !     **********************************************************************
       !
       !--->    warping part of the density (g.ne.0 stars)
       !
       !   ---> d_z^2-Tip
       !
       IF (vacuum%nstm.EQ.2) THEN
          DO l = 1,nv2(jspin)
             DO  l1 = 1,l - 1
                i1 = kvac1(l,jspin) - kvac1(l1,jspin)
                i2 = kvac2(l,jspin) - kvac2(l1,jspin)
                i3 = 0
                IF (iabs(i1).GT.stars%mx1) CYCLE
                IF (iabs(i2).GT.stars%mx2) CYCLE
                ig3 = stars%ig(i1,i2,i3)
                IF (ig3.EQ.0)  CYCLE
                phs = stars%rgphs(i1,i2,i3)
                ig3p = stars%ig(-i1,-i2,i3)
                phsp = stars%rgphs(-i1,-i2,i3)
                ind2 = stars%ig2(ig3)
                ind2p = stars%ig2(ig3p)
                aa = 0.0
                bb = 0.0
                ba = 0.0
                ab = 0.0
                DO n = 1,ne
                   aa = aa + we(n)*CONJG(ac(l1,n,jspin))*ac(l,n,jspin)
                   bb = bb + we(n)*CONJG(bc(l1,n,jspin))*bc(l,n,jspin)
                   ab = ab + we(n)*CONJG(ac(l1,n,jspin))*bc(l,n,jspin)
                   ba = ba + we(n)*CONJG(bc(l1,n,jspin))*ac(l,n,jspin)
                END DO
                aae=-aa*2/3*vacuum%tworkf
                bbe=-bb*2/3*vacuum%tworkf
                abe=-ab*2/3*vacuum%tworkf
                bae=-ba*2/3*vacuum%tworkf
                aaee=aa*4/9*vacuum%tworkf*vacuum%tworkf
                bbee=bb*4/9*vacuum%tworkf*vacuum%tworkf
                abee=ab*4/9*vacuum%tworkf*vacuum%tworkf
                baee=ba*4/9*vacuum%tworkf*vacuum%tworkf
                DO  jz = 1,vacuum%nmzxy
                   ui = u(jz,l,jspin)
                   uj = u(jz,l1,jspin)
                   ddui = ddu(jz,l)
                   dduj = ddu(jz,l1)
                   uei = ue(jz,l,jspin)
                   uej = ue(jz,l1,jspin)
                   dduei = ddue(jz,l)
                   dduej = ddue(jz,l1)
                   t1=aaee*ui*uj+bbee*uei*uej+baee*ui*uej+abee*uei*uj&
                        + aae*(ui*dduj+uj*ddui)+bbe*(uei*dduej+uej*dduei)&
                        + abe*(ui*dduej+uj*dduei)+bae*(ddui*uej+dduj*uei)&
                        + aa*ddui*dduj+bb*dduei*dduej+ba*ddui*dduej&
                        + ab*dduei*dduj  
                   den%vacxy(jz,ind2-1,ivac,jspin) = den%vacxy(jz,ind2-1,ivac,jspin) + t1*phs/stars%nstr2(ind2)
                   den%vacxy(jz,ind2p-1,ivac,jspin)= den%vacxy(jz,ind2p-1,ivac,jspin) + CONJG(t1)*phsp/stars%nstr2(ind2p)
                ENDDO
             ENDDO
          END DO
          !
          ! ---> s-Tip and p_z-Tip
          !
       ELSE
          !=============================================================
          !           continuation of vacden....
          !=============================================================
          IF (noco%l_noco) THEN
             !--->       diagonal elements of the density matrix, n_11 and n_22
             DO ispin = 1,input%jspins
                IF (oneD%odi%d1) THEN
                   DO l = 1,nv2(ispin)
                      DO m = -oneD%odi%mb,oneD%odi%mb
                         lprimee: DO l1 = 1,l-1
                            mprimee: DO m1 = -oneD%odi%mb,m-1
                               i3 = kvac3(l,ispin) - kvac3(l1,ispin)
                               m3 = m-m1
                               IF (m3.EQ.0 .AND. i3.EQ.0) CYCLE mprimee
                               IF (iabs(m3).GT.oneD%odi%M) CYCLE mprimee
                               IF (iabs(i3).GT.stars%mx3) CYCLE lprimee
                               ind1 = oneD%odi%ig(i3,m3)
                               ind1p = oneD%odi%ig(-i3,-m3)
                               IF (ind1.NE.0 .OR. ind1p.NE.0) THEN
                                  aa = CMPLX(0.,0.)
                                  bb = CMPLX(0.,0.)
                                  ba = CMPLX(0.,0.)
                                  ab = CMPLX(0.,0.)
                                  DO n = 1,ne 
                                     aa=aa+we(n)*CONJG(ac_1(l1,m1,n,ispin))*ac_1(l,m,n,ispin)
                                     bb=bb+we(n)*CONJG(bc_1(l1,m1,n,ispin))*bc_1(l,m,n,ispin)
                                     ab=ab+we(n)*CONJG(ac_1(l1,m1,n,ispin))*bc_1(l,m,n,ispin)
                                     ba=ba+we(n)*CONJG(bc_1(l1,m1,n,ispin))*ac_1(l,m,n,ispin)
                                  END DO
                                  xys1: DO jz = 1,vacuum%nmzxy
                                     ui = u_1(jz,l,m,ispin)
                                     uj = u_1(jz,l1,m1,ispin)
                                     uei = ue_1(jz,l,m,ispin)
                                     uej = ue_1(jz,l1,m1,ispin)
                                     t1 = aa*ui*uj + bb*uei*uej + ba*ui*uej + ab*uei*uj
                                     IF (ind1.NE.0) THEN
                                        den%vacxy(jz,ind1-1,ivac,ispin) = den%vacxy(jz,ind1-1,ivac,ispin) + t1/ oneD%odi%nst2(ind1)
                                     END IF
                                     IF (ind1p.NE.0) THEN
                                        den%vacxy(jz,ind1p-1,ivac,ispin) = den%vacxy(jz,ind1p-1,ivac,ispin) + CONJG(t1)/ oneD%odi%nst2(ind1p)
                                     END IF

                                  END DO xys1
                               END IF   ! ind1 and ind1p =0
                            END DO mprimee
                         END DO lprimee
                      END DO  ! m
                   END DO   ! l
                ELSE
                   DO l = 1,nv2(ispin)
                      DO  l1 = 1,l - 1
                         i1 = kvac1(l,ispin) - kvac1(l1,ispin)
                         i2 = kvac2(l,ispin) - kvac2(l1,ispin)
                         i3 = 0
                         IF (iabs(i1).GT.stars%mx1) CYCLE
                         IF (iabs(i2).GT.stars%mx2) CYCLE
                         ig3 = stars%ig(i1,i2,i3)
                         IF (ig3.EQ.0)  CYCLE
                         phs = stars%rgphs(i1,i2,i3)
                         ig3p = stars%ig(-i1,-i2,i3)
                         phsp = stars%rgphs(-i1,-i2,i3)
                         ind2 = stars%ig2(ig3)
                         ind2p = stars%ig2(ig3p)
                         aa = 0.0
                         bb = 0.0
                         ba = 0.0
                         ab = 0.0
                         DO n = 1,ne
                            aa=aa+we(n)*CONJG(ac(l1,n,ispin))*ac(l,n,ispin)
                            bb=bb+we(n)*CONJG(bc(l1,n,ispin))*bc(l,n,ispin)
                            ab=ab+we(n)*CONJG(ac(l1,n,ispin))*bc(l,n,ispin)
                            ba=ba+we(n)*CONJG(bc(l1,n,ispin))*ac(l,n,ispin)
                         END DO
                         DO jz = 1,vacuum%nmzxy
                            ui = u(jz,l,ispin)
                            uj = u(jz,l1,ispin)
                            uei = ue(jz,l,ispin)
                            uej = ue(jz,l1,ispin)
                            t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                            den%vacxy(jz,ind2-1,ivac,ispin) = den%vacxy(jz,ind2-1,ivac,ispin) + t1*phs/stars%nstr2(ind2)
                            den%vacxy(jz,ind2p-1,ivac,ispin) = den%vacxy(jz,ind2p-1,ivac,ispin) + CONJG(t1)*phsp/stars%nstr2(ind2p)
                         ENDDO
                      ENDDO
                   ENDDO
                END IF ! oneD%odi%d1
             END DO
             !--->          off-diagonal element of the density matrix, n_21
             IF (oneD%odi%d1) THEN
                DO l = 1,nv2(1)
                   DO m = -oneD%odi%mb,oneD%odi%mb
                      lprimea: DO l1 = 1,nv2(2)
                         mprimea: DO m1 = -oneD%odi%mb,oneD%odi%mb
                            i3 = kvac3(l,1) - kvac3(l1,2)
                            m3 = m-m1
                            IF (iabs(m3).GT.oneD%odi%M) CYCLE mprimea
                            IF (iabs(i3).GT.stars%mx3) CYCLE lprimea
                            ind1 = oneD%odi%ig(i3,m3)
                            IF (ind1.NE.0) THEN
                               IF (m3.EQ.0 .AND. i3.EQ.0) THEN
                                  aa = CMPLX(0.,0.)
                                  bb = CMPLX(0.,0.)
                                  ba = CMPLX(0.,0.)
                                  ab = CMPLX(0.,0.)
                                  DO n = 1,ne
                                     aa=aa+we(n)*CONJG(ac_1(l1,m1,n,2))* ac_1(l,m,n,1)
                                     bb=bb+we(n)*CONJG(bc_1(l1,m1,n,2))* bc_1(l,m,n,1)
                                     ab=ab+we(n)*CONJG(ac_1(l1,m1,n,2))* bc_1(l,m,n,1)
                                     ba=ba+we(n)*CONJG(bc_1(l1,m1,n,2))* ac_1(l,m,n,1)
                                  END DO
                                  xys3: DO jz = 1,vacuum%nmzxy
                                     ui = u_1(jz,l,m,1)
                                     uj = u_1(jz,l1,m1,2)
                                     uei = ue_1(jz,l,m,1)
                                     uej = ue_1(jz,l1,m1,2)
                                     tempCmplx = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                                     den%vacz(jz,ivac,3) = den%vacz(jz,ivac,3) + REAL(tempCmplx)
                                     den%vacz(jz,ivac,4) = den%vacz(jz,ivac,4) + AIMAG(tempCmplx)
                                  END DO xys3
                               ELSE ! the warped part (ind1 > 1)
                                  aa = CMPLX(0.,0.)
                                  bb = CMPLX(0.,0.)
                                  ba = CMPLX(0.,0.)
                                  ab = CMPLX(0.,0.)
                                  DO n = 1,ne
                                     aa=aa+we(n)*CONJG(ac_1(l1,m1,n,2))* ac_1(l,m,n,1)
                                     bb=bb+we(n)*CONJG(bc_1(l1,m1,n,2))* bc_1(l,m,n,1)
                                     ab=ab+we(n)*CONJG(ac_1(l1,m1,n,2))* bc_1(l,m,n,1)
                                     ba=ba+we(n)*CONJG(bc_1(l1,m1,n,2))* ac_1(l,m,n,1)
                                  END DO
                                  xys2: DO jz = 1,vacuum%nmzxy
                                     ui = u_1(jz,l,m,1)
                                     uj = u_1(jz,l1,m1,2)
                                     uei = ue_1(jz,l,m,1)
                                     uej = ue_1(jz,l1,m1,2)
                                     t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                                     den%vacxy(jz,ind1-1,ivac,3) = den%vacxy(jz,ind1-1,ivac,3) + t1/oneD%odi%nst2(ind1)
                                  END DO xys2
                               END IF  ! the non-warped (ind1 = 1)
                            END IF   ! ind1.ne.0
                         END DO mprimea
                      END DO lprimea
                   END DO  ! m
                END DO   ! l
             ELSE ! oneD%odi%d1
                DO l = 1,nv2(1)
                   DO  l1 = 1,nv2(2)
                      i1 = kvac1(l,1) - kvac1(l1,2)
                      i2 = kvac2(l,1) - kvac2(l1,2)
                      i3 = 0
                      !--->                treat only the warping part
                      IF (iabs(i1).GT.stars%mx1) CYCLE
                      IF (iabs(i2).GT.stars%mx2) CYCLE
                      ig3 = stars%ig(i1,i2,i3)
                      IF (ig3.EQ.0)  CYCLE
                      phs = stars%rgphs(i1,i2,i3)
                      ind2 = stars%ig2(ig3)
                      IF ( ind2.EQ.1) THEN
                         !--->                non-warping part (1st star G=0)
                         aa = 0.0
                         bb = 0.0
                         ba = 0.0
                         ab = 0.0
                         DO ie = 1,ne
                            aa=aa+we(ie)*CONJG(ac(l1,ie,2))*ac(l,ie,1)
                            bb=bb+we(ie)*CONJG(bc(l1,ie,2))*bc(l,ie,1)
                            ab=ab+we(ie)*CONJG(ac(l1,ie,2))*bc(l,ie,1)
                            ba=ba+we(ie)*CONJG(bc(l1,ie,2))*ac(l,ie,1)
                         END DO
                         DO jz = 1,vacuum%nmz
                            ui = u(jz,l,1)
                            ui2 = u(jz,l1,2)
                            uei = ue(jz,l,1)
                            uei2 = ue(jz,l1,2)
                            tempCmplx = aa*ui2*ui + bb*uei2*uei + ab*ui2*uei + ba*uei2*ui
                            den%vacz(jz,ivac,3) = den%vacz(jz,ivac,3) + REAL(tempCmplx)
                            den%vacz(jz,ivac,4) = den%vacz(jz,ivac,4) + AIMAG(tempCmplx)
                         ENDDO
                      ELSE
                         !--->                warping part
                         aa = 0.0
                         bb = 0.0
                         ba = 0.0
                         ab = 0.0
                         DO ie = 1,ne
                            aa=aa + we(ie)*CONJG(ac(l1,ie,2))*ac(l,ie,1)
                            bb=bb + we(ie)*CONJG(bc(l1,ie,2))*bc(l,ie,1)
                            ab=ab + we(ie)*CONJG(ac(l1,ie,2))*bc(l,ie,1)
                            ba=ba + we(ie)*CONJG(bc(l1,ie,2))*ac(l,ie,1)
                         END DO
                         DO  jz = 1,vacuum%nmzxy
                            ui = u(jz,l,1)
                            uj = u(jz,l1,2)
                            uei = ue(jz,l,1)
                            uej = ue(jz,l1,2)
                            t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                            den%vacxy(jz,ind2-1,ivac,3) = den%vacxy(jz, ind2-1,ivac,3) + t1*phs/stars%nstr2(ind2)
                         ENDDO
                      ENDIF
                   ENDDO
                END DO
             END IF ! oneD%odi%d1
          ELSE                                ! collinear part

             IF (oneD%odi%d1) THEN
                DO l = 1,nv2(jspin)
                   DO m = -oneD%odi%mb,oneD%odi%mb
                      lprime: DO l1 = 1,l-1
                         mprime: DO m1 = -oneD%odi%mb,m-1
                            i3 = kvac3(l,jspin) - kvac3(l1,jspin)
                            m3 = m-m1
                            IF (m3.EQ.0 .AND. i3.EQ.0) CYCLE mprime
                            IF (iabs(m3).GT.oneD%odi%M) CYCLE mprime
                            IF (iabs(i3).GT.stars%mx3) CYCLE lprime
                            ind1 = oneD%odi%ig(i3,m3)
                            ind1p = oneD%odi%ig(-i3,-m3)
                            IF (ind1.NE.0 .OR. ind1p.NE.0) THEN
                               aa = CMPLX(0.,0.)
                               bb = CMPLX(0.,0.)
                               ba = CMPLX(0.,0.)
                               ab = CMPLX(0.,0.)
                               DO n = 1,ne
                                  aa=aa+we(n)*CONJG(ac_1(l1,m1,n,jspin))* ac_1(l,m,n,jspin)
                                  bb=bb+we(n)*CONJG(bc_1(l1,m1,n,jspin))* bc_1(l,m,n,jspin)
                                  ab=ab+we(n)*CONJG(ac_1(l1,m1,n,jspin))* bc_1(l,m,n,jspin)
                                  ba=ba+we(n)*CONJG(bc_1(l1,m1,n,jspin))* ac_1(l,m,n,jspin)
                               END DO
                               xys: DO jz = 1,vacuum%nmzxy
                                  ui = u_1(jz,l,m,jspin)
                                  uj = u_1(jz,l1,m1,jspin)
                                  uei = ue_1(jz,l,m,jspin)
                                  uej = ue_1(jz,l1,m1,jspin)
                                  t1 = aa*ui*uj + bb*uei*uej + ba*ui*uej + ab*uei*uj
                                  IF (ind1.NE.0) THEN
                                     den%vacxy(jz,ind1-1,ivac,jspin) = den%vacxy(jz,ind1-1,ivac,jspin) + t1/ oneD%odi%nst2(ind1)
                                  END IF
                                  IF (ind1p.NE.0) THEN
                                     den%vacxy(jz,ind1p-1,ivac,jspin) = den%vacxy(jz,ind1p-1,ivac,jspin) + CONJG(t1)/ oneD%odi%nst2(ind1p)
                                  END IF

                               END DO xys
                            END IF   ! ind1 and ind1p =0
                         END DO mprime
                      END DO lprime
                   END DO  ! m
                END DO   ! l

             ELSE         !D1

                DO l = 1,nv2(jspin)
                   DO  l1 = 1,l - 1
                      i1 = kvac1(l,jspin) - kvac1(l1,jspin)
                      i2 = kvac2(l,jspin) - kvac2(l1,jspin)
                      i3 = 0
                      IF (iabs(i1).GT.stars%mx1) CYCLE
                      IF (iabs(i2).GT.stars%mx2) CYCLE
                      ig3 = stars%ig(i1,i2,i3)
                      IF (ig3.EQ.0)  CYCLE
                      phs = stars%rgphs(i1,i2,i3)
                      ig3p = stars%ig(-i1,-i2,i3)
                      phsp = stars%rgphs(-i1,-i2,i3)
                      ind2 = stars%ig2(ig3)
                      ind2p = stars%ig2(ig3p)
                      aa = 0.0
                      bb = 0.0
                      ba = 0.0
                      ab = 0.0
                      DO n = 1,ne
                         aa=aa+we(n)*CONJG(ac(l1,n,jspin))*ac(l,n,jspin)
                         bb=bb+we(n)*CONJG(bc(l1,n,jspin))*bc(l,n,jspin)
                         ab=ab+we(n)*CONJG(ac(l1,n,jspin))*bc(l,n,jspin)
                         ba=ba+we(n)*CONJG(bc(l1,n,jspin))*ac(l,n,jspin)
                      END DO
                      DO  jz = 1,vacuum%nmzxy
                         ui = u(jz,l,jspin)
                         uj = u(jz,l1,jspin)
                         uei = ue(jz,l,jspin)
                         uej = ue(jz,l1,jspin)
                         t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                         den%vacxy(jz,ind2-1,ivac,jspin) = den%vacxy(jz,ind2-1, ivac,jspin) + t1*phs/stars%nstr2(ind2)
                         den%vacxy(jz,ind2p-1,ivac,jspin) = den%vacxy(jz,ind2p-1, ivac,jspin) + CONJG(t1)*phsp/stars%nstr2(ind2p)
                      ENDDO
                   ENDDO
                END DO
             END IF ! D1
          ENDIF
       END IF
       !=============================================================
       !
       !       calculate 1. to nstars. starcoefficient for each k and energy eigenvalue 
       !           to dos%qstars(ne,layer,ivac,ikpt) if starcoeff=T (the star coefficient values are written to vacdos)
       !
       IF (vacuum%starcoeff .AND. banddos%vacdos) THEN
          DO  n=1,ne
             DO l = 1,nv2(jspin)
                DO  l1 = 1,l - 1
                   i1 = kvac1(l,jspin) - kvac1(l1,jspin)
                   i2 = kvac2(l,jspin) - kvac2(l1,jspin)
                   i3 = 0
                   IF (iabs(i1).GT.stars%mx1) CYCLE
                   IF (iabs(i2).GT.stars%mx2) CYCLE
                   ig3 = stars%ig(i1,i2,i3)
                   IF (ig3.EQ.0)  CYCLE
                   ind2 = stars%ig2(ig3)
                   ig3p = stars%ig(-i1,-i2,i3)
                   ind2p = stars%ig2(ig3p)
                   IF ((ind2.GE.2.AND.ind2.LE.vacuum%nstars).OR.&
                        (ind2p.GE.2.AND.ind2p.LE.vacuum%nstars)) THEN
                      phs = stars%rgphs(i1,i2,i3)
                      phsp = stars%rgphs(-i1,-i2,i3)
                      aa = CONJG(ac(l1,n,jspin))*ac(l,n,jspin)
                      bb = CONJG(bc(l1,n,jspin))*bc(l,n,jspin)
                      ab = CONJG(ac(l1,n,jspin))*bc(l,n,jspin)
                      ba = CONJG(bc(l1,n,jspin))*ac(l,n,jspin)
                      DO jj = 1,vacuum%layers
                         ui = u(vacuum%izlay(jj,1),l,jspin)
                         uj = u(vacuum%izlay(jj,1),l1,jspin)
                         uei = ue(vacuum%izlay(jj,1),l,jspin)
                         uej = ue(vacuum%izlay(jj,1),l1,jspin)
                         t1 = aa*ui*uj + bb*uei*uej +ba*ui*uej + ab*uei*uj
                         IF (ind2.GE.2.AND.ind2.LE.vacuum%nstars) &
                              dos%qstars(ind2-1,n,jj,ivac,ikpt,jspin) = dos%qstars(ind2-1,n,jj,ivac,ikpt,jspin)+ t1*phs/stars%nstr2(ind2)
                         IF (ind2p.GE.2.AND.ind2p.LE.vacuum%nstars) &
                              dos%qstars(ind2p-1,n,jj,ivac,ikpt,jspin) = dos%qstars(ind2p-1,n,jj,ivac,ikpt,jspin) +CONJG(t1)*phs/stars%nstr2(ind2p)
                      END DO
                   END IF
                ENDDO
             END DO
          ENDDO
       END IF
    ENDDO
    DEALLOCATE (ac,bc,dt,dte,du,ddu,due,ddue,t,te,tei,u,ue,v,yy )

    IF (oneD%odi%d1) THEN
       DEALLOCATE (ac_1,bc_1,dt_1,dte_1,du_1,ddu_1,due_1,ddue_1)
       DEALLOCATE (t_1,te_1,tei_1,u_1,ue_1)
    END IF ! oneD%odi%d1

    IF(vacuum%nvac.EQ.1) THEN
       den%vacz(:,2,:) = den%vacz(:,1,:)
       IF (sym%invs) THEN
          den%vacxy(:,:,2,:) = CONJG(den%vacxy(:,:,1,:))
       ELSE
          den%vacxy(:,:,2,:) = den%vacxy(:,:,1,:)
       END IF
    END IF

    CALL timestop("vacden")

  END SUBROUTINE vacden
END MODULE m_vacden
