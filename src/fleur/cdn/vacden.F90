!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vacden
   USE m_juDFT
   !! Vacuum contribution to the valence charge density of a film.
   !!
   !! Determines the 2D star-function expansion coefficients of the vacuum charge
   !! density. Speed up by R. Wu 1992; vacuum DOS added by S. Heinze (Jan. 1996),
   !! star coefficients (Jan. 1999); non-collinear density matrix by P. Kurz (1999/07).
   !!
   !! In non-collinear calculations the density is a hermitian 2x2 matrix. Its
   !! diagonal elements n_11 and n_22 go into `den%vac(:,:,:,1)` and `den%vac(:,:,:,2)`,
   !! the (complex) off-diagonal element n_21 into `den%vac(:,:,:,3)`.
   implicit none
   PRIVATE
   PUBLIC :: vacden

   !******** ABBREVIATIONS ************************************************
   !     qvac     : vacuum charge of each eigenstate, needed in in cdnval
   !                to determine the vacuum energy parameters
   !     vz       : non-warping part of the vacuum potential (matrix)
   !                collinear    : 2. index = ivac (# of vaccum)
   !                non-collinear: 2. index = ipot (comp. of pot. matr.)
   !     den%vac  : vacuum density matrix, indexed (z, 2D star, ivac, spin).
   !                Star index 1 is the non-warping (G_|| = 0) part and is filled
   !                up to vacuum%nmz; all further stars are the warping part and
   !                are filled up to vacuum%nmzxy only.
   !                Spin slots 1 and 2 hold the diagonal elements n_11 and n_22,
   !                slot 3 (non-collinear only) the off-diagonal element n_21.
   !***********************************************************************
   !     *******************************************************************************
   !    layers: no. of layers to be calculated (in vertical direction with z-values as given by izlay)
   !    izlay : defines vertical position of layers in delz (=0.1 a.u.) units from begining of vacuum region
   !    vacdos: =T: calculate vacuum dos in layers as given by the above
   !    integ : =T: integrate in vertical position between izlay(layer,1)..izlay(layer,2)
   !    starcoeff: =T: star coefficients are calculated at values of izlay for 0th (=q) to nstars-1. star
   !                (vacdos%qstars(1..nstars-1))
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
CONTAINS

   SUBROUTINE vacden(vacuum,stars,input,cell,atoms,noco,nococonv,banddos,&
                     we,ikpt,jspin,vz,ne,ev_list,lapw,evac,den,zMat,vacdos,dos,lapwq,we1,zMat1)
      !! Calculates the vacuum part of the density and puts it into den%vac. The variable has
      !! four dimensions: The z, star, vacuum and spin index. Recent refactoring combined the
      !! real variable vacz and the complex star expansion vacxy into one. vacz is identical
      !! to the array for the first star (\(\boldsymbol{G}_{||}=0\)), which was not present in vacxy. The
      !! array is bounded by vacuum%nmzd in the first dimension, but only filled up to
      !! vacuum%nmzxyd for any star but the first.
      !!
      !! This routine is a driver: it sets up the vacuum basis (`t_vacbasis`) for each
      !! vacuum and then hands the accumulation over to the private workers below.
      USE m_types
      USE m_types_vacbasis
      USE m_types_vacdos
      USE m_types_dos

      IMPLICIT NONE

      TYPE(t_lapw),     INTENT(IN)    :: lapw !I just removed the assignment and made lapw INTENT(IN).
      TYPE(t_banddos),  INTENT(IN)    :: banddos
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_vacuum),   INTENT(IN)    :: vacuum
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_stars),    INTENT(IN)    :: stars
      TYPE(t_cell),     INTENT(IN)    :: cell
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_mat),      INTENT(IN)    :: zMat
      TYPE(t_potden),   INTENT(INOUT) :: den
      TYPE(t_vacdos),   INTENT(INOUT) :: vacdos
      TYPE(t_dos),      INTENT(INOUT) :: dos

      INTEGER,          INTENT (IN)   :: jspin
      INTEGER,          INTENT (IN)   :: ne
      INTEGER,          INTENT (IN)   :: ikpt

      INTEGER,          INTENT(IN)    :: ev_list(ne)
      REAL,             INTENT(IN)    :: evac(2,input%jspins)
      REAL,             INTENT(IN)    :: we(input%neig)
      REAL,             INTENT(IN)    :: vz(:,:,:) !(vacuum%nmzd,ivac,ispin)

      ! Optional DFPT variables
      TYPE(t_lapw), OPTIONAL, INTENT(IN) :: lapwq
      TYPE(t_mat),  OPTIONAL, INTENT(IN) :: zMat1
      REAL,         OPTIONAL, INTENT(IN) :: we1(:)

      !! Basis in the vacuum: vb at k, vbq at k+q (DFPT only).
      TYPE(t_vacbasis) :: vb, vbq

      REAL    :: const, zsign, wronk
      REAL    :: gshift(2,2), noshift(2,2)
      INTEGER :: ivac, ispin, k
      LOGICAL :: l_dfpt, l_center0

      REAL, ALLOCATABLE :: we1_l(:) !! local copy of we1, safe to name in an OMP clause

      CALL timestart("vacden")

      l_dfpt    = PRESENT(zMat1)
      l_center0 = norm2(stars%center) < 1e-8

      noshift = 0.0

      CALL vb%init(lapw, input%jspins, vacuum%nmzd, input%neig)
      IF (l_dfpt) THEN
         CALL vbq%init(lapwq, input%jspins, vacuum%nmzd, input%neig)
         we1_l = we1
      END IF

      IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) THEN
         ! Without a spin spiral both spins share one basis; mirror spin 1 onto spin 2.
         ! (lapw is INTENT(IN), so lapw%nv(2) and lapw%k3(:,2) are relied on to match.)
         vb%nv2(2) = vb%nv2(1)
         DO k = 1,vb%nv2(1)
            vb%kvac(:,k,2) = vb%kvac(:,k,1)
         END DO
         DO k = 1,lapw%nv(1)
            vb%map2(k,2) = vb%map2(k,1)
         END DO
      END IF

      wronk = vb%wronk
      const = 1.0 / ( SQRT(cell%omtil)*wronk )

      DO ivac = 1,vacuum%nvac
         vb%ac(:,:,:) = CMPLX(0.0,0.0)
         vb%bc(:,:,:) = CMPLX(0.0,0.0)
         IF (l_dfpt) THEN
            vbq%ac(:,:,:) = CMPLX(0.0,0.0)
            vbq%bc(:,:,:) = CMPLX(0.0,0.0)
         END IF
         zsign = 3. - 2.*ivac

         !---> vacuum wave functions and their A/B expansion coefficients
         IF (noco%l_noco) THEN
            !--->    In a non-collinear calculation vacden is only called once.
            !--->    Thus, the vaccum wavefunctions and the A- and B-coeff. (ac bc)
            !--->    have to be calculated for both spins on that call.
            !--->       setup the spin-spiral q-vector
            gshift(1:2,1) = - nococonv%qss(1:2)/2
            gshift(1:2,2) = + nococonv%qss(1:2)/2

            CALL vb%calc_radfun(vacuum, cell, lapw%bkpt, gshift, evac(ivac,:), vz(:,ivac,:), 1, input%jspins)
            DO ispin = 1,input%jspins
               !--->       the coefficients of the spin-down basis functions are
               !--->       stored in the second half of the eigenvector
               CALL vb%calc_abcoeff(lapw, cell, zMat, ne, ispin, zsign, const, &
                                    (lapw%nv(1)+atoms%nlotot)*(ispin-1), 1.0)
            END DO
         ELSE
            CALL vb%calc_radfun(vacuum, cell, lapw%bkpt, noshift, evac(ivac,:), vz(:,ivac,:), jspin, jspin)
            CALL vb%calc_abcoeff(lapw, cell, zMat, ne, jspin, zsign, const, 0, 1.0)
            IF (l_dfpt) THEN
               ! Same unperturbed potential and energy parameter, basis shifted by the phonon q.
               gshift(1,:) = lapwq%qphon(1)
               gshift(2,:) = lapwq%qphon(2)
               CALL vbq%calc_radfun(vacuum, cell, lapwq%bkpt, gshift, evac(ivac,:), vz(:,ivac,:), jspin, jspin)
               CALL vbq%calc_abcoeff(lapwq, cell, zMat1, ne, jspin, zsign, const, 0, 2.0)
            END IF
         END IF

         !---> non-warping part of the density (G_|| = 0 star) and the vacuum charges.
         !     For DFPT with a shifted star centre this term is produced by the warping
         !     loop instead, see priv_den_warp_col.
         IF (noco%l_noco) THEN
            CALL priv_den_g0_noco(vb, den, vacdos, dos, vacuum, cell, input, &
                                  we, ne, ev_list, ikpt, ivac)
         ELSE IF ((.NOT.l_dfpt).OR.l_center0) THEN
            CALL priv_den_g0_col(vb, vbq, den, vacdos, dos, vacuum, cell, &
                                 we, we1_l, ne, ev_list, ikpt, ivac, jspin, l_dfpt)
         END IF

         !---> layer-resolved vacuum DOS
         IF (banddos%vacdos) THEN
            CALL priv_dos_layers(vb, vacdos, banddos, input, vacuum, noco, &
                                 ne, ev_list, ikpt, ivac, jspin)
         END IF

         !---> warping part of the density (G_|| /= 0 stars)
         IF (noco%l_noco) THEN
            CALL timestart("vacden_warp_noco")
            CALL priv_den_warp_noco(vb, den, stars, vacuum, input, we, ne, ivac)
            CALL priv_den_offdiag(vb, den, stars, vacuum, we, ne, ivac)
            CALL timestop("vacden_warp_noco")
         ELSE
            CALL timestart("vacden_warp_col")
            CALL priv_den_warp_col(vb, vbq, den, stars, vacuum, we, we1_l, ne, ivac, jspin, &
                                   l_dfpt, l_center0)
            CALL timestop("vacden_warp_col")
         END IF

         !---> star coefficients 1..nstars for each k-point and eigenvalue
         IF (banddos%starcoeff .AND. banddos%vacdos) THEN
            CALL priv_starcoeff(vb, vacdos, banddos, stars, ne, ev_list, ikpt, ivac, jspin)
         END IF
      END DO ! vacuum loop

      ! For vacuum%nvac==1 the lower vacuum is filled from the upper one downstream
      ! via stars%fill_2nd_vac (see mix.F90 / cdn_io.F90).

      CALL timestop("vacden")

   END SUBROUTINE vacden

   !----------------------------------------------------------------------------------
   ! Private workers
   !----------------------------------------------------------------------------------

   PURE SUBROUTINE priv_pair_sum(acL, bcL, acR, bcR, we, ne, ip, ig, aa, bb, ab, ba)
      !! The weighted pair products \(\sum_n w_n \overline{X_L(p,n)} Y_R(g,n)\) that every
      !! density accumulation below is built from.
      !!
      !! The arrays are passed whole (a contiguous spin slice) rather than as `(:)`
      !! sections, so no temporary copy is made in the hot loop.
      IMPLICIT NONE
      COMPLEX, INTENT(IN)  :: acL(:,:), bcL(:,:), acR(:,:), bcR(:,:) !! (nv2d, neig)
      REAL,    INTENT(IN)  :: we(:)
      INTEGER, INTENT(IN)  :: ne, ip, ig
      COMPLEX, INTENT(OUT) :: aa, bb, ab, ba

      INTEGER :: n

      aa = 0.0
      bb = 0.0
      ba = 0.0
      ab = 0.0
      DO n = 1,ne
         aa = aa + we(n)*CONJG(acL(ip,n))*acR(ig,n)
         bb = bb + we(n)*CONJG(bcL(ip,n))*bcR(ig,n)
         ab = ab + we(n)*CONJG(acL(ip,n))*bcR(ig,n)
         ba = ba + we(n)*CONJG(bcL(ip,n))*acR(ig,n)
      END DO

   END SUBROUTINE priv_pair_sum

   SUBROUTINE priv_den_g0_noco(vb, den, vacdos, dos, vacuum, cell, input, &
                               we, ne, ev_list, ikpt, ivac)
      !! Non-warping (G_|| = 0) part of the diagonal density matrix elements n_11 and
      !! n_22, plus the vacuum charge of each eigenstate. The non-warping part of the
      !! off-diagonal n_21 is done together with its warping part in priv_den_offdiag.
      USE m_types
      USE m_types_vacbasis
      USE m_types_vacdos
      USE m_types_dos
      IMPLICIT NONE

      TYPE(t_vacbasis), INTENT(IN)    :: vb
      TYPE(t_potden),   INTENT(INOUT) :: den
      TYPE(t_vacdos),   INTENT(INOUT) :: vacdos
      TYPE(t_dos),      INTENT(INOUT) :: dos
      TYPE(t_vacuum),   INTENT(IN)    :: vacuum
      TYPE(t_cell),     INTENT(IN)    :: cell
      TYPE(t_input),    INTENT(IN)    :: input
      REAL,             INTENT(IN)    :: we(:)
      INTEGER,          INTENT(IN)    :: ne, ev_list(:), ikpt, ivac

      INTEGER :: ispin, ikG, n, jz
      REAL    :: qout, ui, uei
      COMPLEX :: aa, bb, ab, ba

      DO ispin = 1,input%jspins
         DO ikG = 1,vb%nv2(ispin)
            aa = 0.0
            bb = 0.0
            ba = 0.0
            ab = 0.0
            DO n = 1,ne
               aa=aa + we(n)*CONJG(vb%ac(ikG,n,ispin))*vb%ac(ikG,n,ispin)
               bb=bb + we(n)*CONJG(vb%bc(ikG,n,ispin))*vb%bc(ikG,n,ispin)
               ab=ab + we(n)*CONJG(vb%ac(ikG,n,ispin))*vb%bc(ikG,n,ispin)
               ba=ba + we(n)*CONJG(vb%bc(ikG,n,ispin))*vb%ac(ikG,n,ispin)
               qout = REAL(CONJG(vb%ac(ikG,n,ispin))*vb%ac(ikG,n,ispin)+vb%ddnv(ikG,ispin)*CONJG(vb%bc(ikG,n,ispin))*vb%bc(ikG,n,ispin))
               if (vacdos%l_initialized) vacdos%qvac(ev_list(n),ivac,ikpt,ispin) = vacdos%qvac(ev_list(n),ivac,ikpt,ispin) + qout*cell%area
               if (dos%l_initialized) dos%qTot(ev_list(n),ikpt,ispin) = dos%qTot(ev_list(n),ikpt,ispin) + qout*cell%area
            END DO
            DO jz = 1,vacuum%nmz
               ui = vb%u(jz,ikG,ispin)
               uei = vb%ue(jz,ikG,ispin)
               den%vac(jz,1,ivac,ispin) = den%vac(jz,1,ivac,ispin) + REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei) ! TODO: AN TB; sollte man das REAL killen?
            END DO
         END DO
      END DO

   END SUBROUTINE priv_den_g0_noco

   SUBROUTINE priv_den_g0_col(vb, vbq, den, vacdos, dos, vacuum, cell, &
                              we, we1_l, ne, ev_list, ikpt, ivac, jspin, l_dfpt)
      !! Non-warping (G_|| = 0) part of the collinear density plus the vacuum charge
      !! of each eigenstate. In the DFPT case the response density picks up both the
      !! occupation response (we1 with the unperturbed coefficients) and the
      !! coefficient response (we with vbq).
      USE m_types
      USE m_types_vacbasis
      USE m_types_vacdos
      USE m_types_dos
      IMPLICIT NONE

      TYPE(t_vacbasis), INTENT(IN)    :: vb, vbq
      TYPE(t_potden),   INTENT(INOUT) :: den
      TYPE(t_vacdos),   INTENT(INOUT) :: vacdos
      TYPE(t_dos),      INTENT(INOUT) :: dos
      TYPE(t_vacuum),   INTENT(IN)    :: vacuum
      TYPE(t_cell),     INTENT(IN)    :: cell
      REAL,             INTENT(IN)    :: we(:)
      REAL, ALLOCATABLE, INTENT(IN)   :: we1_l(:)
      INTEGER,          INTENT(IN)    :: ne, ev_list(:), ikpt, ivac, jspin
      LOGICAL,          INTENT(IN)    :: l_dfpt

      INTEGER :: ikG, n, jz
      REAL    :: qout, ui, uei
      COMPLEX :: aa, bb, ab, ba

      DO ikG = 1,vb%nv2(jspin)
         aa = CMPLX(0.0,0.0)
         bb = CMPLX(0.0,0.0)
         ba = CMPLX(0.0,0.0)
         ab = CMPLX(0.0,0.0)
         DO n = 1,ne
            IF (.NOT.l_dfpt) THEN
               aa = aa + we(n)*CONJG(vb%ac(ikG,n,jspin))*vb%ac(ikG,n,jspin)
               bb = bb + we(n)*CONJG(vb%bc(ikG,n,jspin))*vb%bc(ikG,n,jspin)
               ab = ab + we(n)*CONJG(vb%ac(ikG,n,jspin))*vb%bc(ikG,n,jspin)
               ba = ba + we(n)*CONJG(vb%bc(ikG,n,jspin))*vb%ac(ikG,n,jspin)
            ELSE
               aa = aa + we1_l(n)*CONJG(vb%ac(ikG,n,jspin))*vb%ac(ikG,n,jspin)
               bb = bb + we1_l(n)*CONJG(vb%bc(ikG,n,jspin))*vb%bc(ikG,n,jspin)
               ab = ab + we1_l(n)*CONJG(vb%ac(ikG,n,jspin))*vb%bc(ikG,n,jspin)
               ba = ba + we1_l(n)*CONJG(vb%bc(ikG,n,jspin))*vb%ac(ikG,n,jspin)
               aa = aa + we(n)*CONJG(vb%ac(ikG,n,jspin))*vbq%ac(ikG,n,jspin)
               bb = bb + we(n)*CONJG(vb%bc(ikG,n,jspin))*vbq%bc(ikG,n,jspin)
               ab = ab + we(n)*CONJG(vb%ac(ikG,n,jspin))*vbq%bc(ikG,n,jspin)
               ba = ba + we(n)*CONJG(vb%bc(ikG,n,jspin))*vbq%ac(ikG,n,jspin)
            END IF
            qout = REAL(CONJG(vb%ac(ikG,n,jspin))*vb%ac(ikG,n,jspin)+vb%ddnv(ikG,jspin)*CONJG(vb%bc(ikG,n,jspin))*vb%bc(ikG,n,jspin))
            if(vacdos%l_initialized) vacdos%qvac(ev_list(n),ivac,ikpt,jspin) = vacdos%qvac(ev_list(n),ivac,ikpt,jspin) + qout*cell%area
            if (dos%l_initialized) dos%qTot(ev_list(n),ikpt,jspin) = dos%qTot(ev_list(n),ikpt,jspin) + qout*cell%area
         END DO
         DO  jz = 1,vacuum%nmz
            ui = vb%u(jz,ikG,jspin)
            uei = vb%ue(jz,ikG,jspin)
            IF (.NOT.l_dfpt) THEN
               den%vac(jz,1,ivac,jspin) = den%vac(jz,1,ivac,jspin) +REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei) ! TODO: REAL weg?
            ELSE
               den%vac(jz,1,ivac,jspin) = den%vac(jz,1,ivac,jspin) + aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei
            END IF
         END DO
      END DO

   END SUBROUTINE priv_den_g0_col

   SUBROUTINE priv_dos_layers(vb, vacdos, banddos, input, vacuum, noco, &
                              ne, ev_list, ikpt, ivac, jspin)
      !! Layer-resolved vacuum DOS, either integrated over the whole 2D unit cell or,
      !! if locx/locy differ, over the rectangle spanned by
      !! (locx(1),locy(1)) .. (locx(2),locy(2)) in internal coordinates.
      USE m_constants
      USE m_qsf
      USE m_types
      USE m_types_vacbasis
      USE m_types_vacdos
      IMPLICIT NONE

      TYPE(t_vacbasis), INTENT(IN)    :: vb
      TYPE(t_vacdos),   INTENT(INOUT) :: vacdos
      TYPE(t_banddos),  INTENT(IN)    :: banddos
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_vacuum),   INTENT(IN)    :: vacuum
      TYPE(t_noco),     INTENT(IN)    :: noco
      INTEGER,          INTENT(IN)    :: ne, ev_list(:), ikpt, ivac, jspin

      INTEGER :: ispin, isp_start, isp_end, ikG, ikGPr, n, jj, ii, ll
      REAL    :: ui, uei, uj, uej, k_diff, k_d1, k_d2, qlay(1)
      REAL    :: yy(vacuum%nmzd)
      COMPLEX :: aa, bb, ab, ba, factorx, factory
      REAL, PARAMETER :: eps = 0.01

      IF (ABS(banddos%locx(1)-banddos%locx(2)).LE.eps) THEN
         !----> integrated over 2D-unit cell
         IF (noco%l_noco) THEN
            isp_start = 1
            isp_end   = input%jspins
         ELSE
            isp_start = jspin
            isp_end   = jspin
         END IF
         DO ispin = isp_start, isp_end
            DO ikG=1,vb%nv2(ispin)
               DO n = 1,ne
                  aa = CONJG(vb%ac(ikG,n,ispin))*vb%ac(ikG,n,ispin)
                  bb = CONJG(vb%bc(ikG,n,ispin))*vb%bc(ikG,n,ispin)
                  ab = CONJG(vb%ac(ikG,n,ispin))*vb%bc(ikG,n,ispin)
                  ba = CONJG(vb%bc(ikG,n,ispin))*vb%ac(ikG,n,ispin)
                  DO jj = 1,banddos%layers
                     !---> either integrated (z1,z2) or slice (z1)
                     IF (input%integ) THEN
                        ll = 1
                        DO ii = banddos%izlay(jj,1),banddos%izlay(jj,2)
                           ui = vb%u(ii,ikG,ispin)
                           uei = vb%ue(ii,ikG,ispin)
                           yy(ll) = REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)
                           ll = ll+1
                        END DO
                        CALL qsf(vacuum%delz,yy,qlay,ll-1,0)
                        vacdos%qvlay(ev_list(n),jj,ivac,ikpt,ispin) = vacdos%qvlay(ev_list(n),jj,ivac,ikpt,ispin) + qlay(1)
                     ELSE
                        ui = vb%u(banddos%izlay(jj,1),ikG,ispin)
                        uei = vb%ue(banddos%izlay(jj,1),ikG,ispin)
                        vacdos%qvlay(ev_list(n),jj,ivac,ikpt,ispin) = vacdos%qvlay(ev_list(n),jj,ivac,ikpt,ispin) &
                                                                    + REAL(aa*ui*ui+bb*uei*uei+(ab+ba)*ui*uei)

                     END IF
                  END DO
               END DO
            END DO
         END DO
      ELSE
         !----> if LDOS should be calculated over restricted area of the 2D-unit cell
         !     lower left corner: (locx(1), locy(1))   }  in internal
         !     upper right corner: (locx(2), locy(2))  }  coordinates
         !
         DO ikG=1, vb%nv2(jspin)
            DO ikGPr=1, vb%nv2(jspin)
               IF (vb%kvac(1,ikG,jspin).EQ.vb%kvac(1,ikGPr,jspin)) THEN
                  factorx = CMPLX((banddos%locx(2)-banddos%locx(1)), 0.)
               ELSE
                  k_diff=tpi_const*(vb%kvac(1,ikG,jspin)-vb%kvac(1,ikGPr,jspin))
                  k_d1 = k_diff*banddos%locx(1)
                  k_d2 = k_diff*banddos%locx(2)
                  factorx=(CMPLX( COS(k_d2), SIN(k_d2)) - &
                           CMPLX( COS(k_d1), SIN(k_d1)) ) / &
                           CMPLX( 0.,k_diff )
               END IF
               IF (vb%kvac(2,ikG,jspin).EQ.vb%kvac(2,ikGPr,jspin)) THEN
                  factory = CMPLX((banddos%locy(2)-banddos%locy(1)), 0.)
               ELSE
                  k_diff=tpi_const*(vb%kvac(2,ikG,jspin)-vb%kvac(2,ikGPr,jspin))
                  k_d1 = k_diff*banddos%locy(1)
                  k_d2 = k_diff*banddos%locy(2)
                  factory=(CMPLX( COS(k_d2), SIN(k_d2)) - &
                           CMPLX( COS(k_d1), SIN(k_d1)) ) / &
                           CMPLX( 0.,k_diff )
               END IF
               DO n=1, ne
                  aa = CONJG(vb%ac(ikGPr,n,jspin))*vb%ac(ikG,n,jspin)
                  bb = CONJG(vb%bc(ikGPr,n,jspin))*vb%bc(ikG,n,jspin)
                  ab = CONJG(vb%ac(ikGPr,n,jspin))*vb%bc(ikG,n,jspin)
                  ba = CONJG(vb%bc(ikGPr,n,jspin))*vb%ac(ikG,n,jspin)
                  DO jj = 1,banddos%layers
                     !---> either integrated (z1,z2) or slice (z1)
                     IF (input%integ) THEN
                        ll = 1
                        DO ii = banddos%izlay(jj,1), banddos%izlay(jj,2)
                           ui = vb%u(ii,ikG,jspin)
                           uei = vb%ue(ii,ikG,jspin)
                           uj = vb%u(ii,ikGPr,jspin)
                           uej = vb%ue(ii,ikGPr,jspin)
                           yy(ll) = REAL((aa*ui*uj+bb*uei*uej+ab*uei*uj+ba*ui*uej)*factorx*factory)
                           ll = ll+1
                        END DO
                        CALL qsf(vacuum%delz,yy,qlay,ll-1,0)
                        vacdos%qvlay(ev_list(n),jj,ivac,ikpt,jspin) = vacdos%qvlay(ev_list(n),jj,ivac,ikpt,jspin) + qlay(1)
                     ELSE
                        ui = vb%u(banddos%izlay(jj,1),ikG,jspin)
                        uei = vb%ue(banddos%izlay(jj,1),ikG,jspin)
                        uj = vb%u(banddos%izlay(jj,1),ikGPr,jspin)
                        uej = vb%ue(banddos%izlay(jj,1),ikGPr,jspin)
                        vacdos%qvlay(ev_list(n),jj,ivac,ikpt,jspin) = vacdos%qvlay(ev_list(n),jj,ivac,ikpt,jspin)&
                                                                  +REAL((aa*ui*uj + bb*uei*uej+ab*uei*uj+ba*ui*uej)*factorx*factory)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE priv_dos_layers

   SUBROUTINE priv_den_warp_noco(vb, den, stars, vacuum, input, we, ne, ivac)
      !! Warping part (G_|| /= 0) of the diagonal density matrix elements n_11 and n_22.
      USE m_types
      USE m_types_vacbasis
      IMPLICIT NONE

      TYPE(t_vacbasis), INTENT(IN)    :: vb
      TYPE(t_potden),   INTENT(INOUT) :: den
      TYPE(t_stars),    INTENT(IN)    :: stars
      TYPE(t_vacuum),   INTENT(IN)    :: vacuum
      TYPE(t_input),    INTENT(IN)    :: input
      REAL,             INTENT(IN)    :: we(:)
      INTEGER,          INTENT(IN)    :: ne, ivac

      INTEGER :: ispin, ikG, ikGPr, i1, i2, i3, ig3, ind2, ind2p, jz
      REAL    :: ui, uj, uei, uej
      COMPLEX :: aa, bb, ab, ba, t1, phs, phsp

      DO ispin = 1,input%jspins
         DO ikG = 1,vb%nv2(ispin)
            DO ikGPr = 1, ikG - 1
               i1 = vb%kvac(1,ikG,ispin) - vb%kvac(1,ikGPr,ispin)
               i2 = vb%kvac(2,ikG,ispin) - vb%kvac(2,ikGPr,ispin)
               i3 = 0
               IF (iabs(i1).GT.stars%mx1) CYCLE
               IF (iabs(i2).GT.stars%mx2) CYCLE
               ig3 = stars%ig(i1,i2,i3)
               IF (ig3.EQ.0)  CYCLE
               phs = CONJG(stars%r2gphs(i1,i2))
               phsp = CONJG(stars%r2gphs(-i1,-i2))
               ind2 = stars%i2g(i1,i2)
               ind2p = stars%i2g(-i1,-i2)
               CALL priv_pair_sum(vb%ac(:,:,ispin), vb%bc(:,:,ispin), &
                                  vb%ac(:,:,ispin), vb%bc(:,:,ispin), &
                                  we, ne, ikGPr, ikG, aa, bb, ab, ba)
               DO jz = 1,vacuum%nmzxy
                  ui = vb%u(jz,ikG,ispin)
                  uj = vb%u(jz,ikGPr,ispin)
                  uei = vb%ue(jz,ikG,ispin)
                  uej = vb%ue(jz,ikGPr,ispin)
                  t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                  den%vac(jz,ind2,ivac,ispin) = den%vac(jz,ind2,ivac,ispin) + t1*phs/stars%nstr2(ind2)
                  den%vac(jz,ind2p,ivac,ispin) = den%vac(jz,ind2p,ivac,ispin) + CONJG(t1)*phsp/stars%nstr2(ind2p)
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE priv_den_warp_noco

   SUBROUTINE priv_den_offdiag(vb, den, stars, vacuum, we, ne, ivac)
      !! Off-diagonal element n_21 of the non-collinear density matrix. Unlike the
      !! diagonal elements, its non-warping (1st star) part is accumulated here too,
      !! because both parts come from the same spin-1/spin-2 pair loop.
      USE m_types
      USE m_types_vacbasis
      IMPLICIT NONE

      TYPE(t_vacbasis), INTENT(IN)    :: vb
      TYPE(t_potden),   INTENT(INOUT) :: den
      TYPE(t_stars),    INTENT(IN)    :: stars
      TYPE(t_vacuum),   INTENT(IN)    :: vacuum
      REAL,             INTENT(IN)    :: we(:)
      INTEGER,          INTENT(IN)    :: ne, ivac

      INTEGER :: ikG, ikGPr, i1, i2, i3, ig3, ind2, jz
      REAL    :: ui, uj, uei, uej, ui2, uei2
      COMPLEX :: aa, bb, ab, ba, t1, phs, tempCmplx

      DO ikG = 1,vb%nv2(1)
         DO  ikGPr = 1,vb%nv2(2)
            i1 = vb%kvac(1,ikG,1) - vb%kvac(1,ikGPr,2)
            i2 = vb%kvac(2,ikG,1) - vb%kvac(2,ikGPr,2)
            i3 = 0
            !--->                treat only the warping part
            IF (iabs(i1).GT.stars%mx1) CYCLE
            IF (iabs(i2).GT.stars%mx2) CYCLE
            ig3 = stars%ig(i1,i2,i3)
            IF (ig3.EQ.0)  CYCLE
            phs = CONJG(stars%r2gphs(i1,i2)) ! 2D star phase, see above
            ind2 = stars%i2g(i1,i2)
            CALL priv_pair_sum(vb%ac(:,:,2), vb%bc(:,:,2), &
                               vb%ac(:,:,1), vb%bc(:,:,1), &
                               we, ne, ikGPr, ikG, aa, bb, ab, ba)
            IF ( ind2.EQ.1) THEN
               !--->                non-warping part (1st star G=0)
               DO jz = 1,vacuum%nmz
                  ui = vb%u(jz,ikG,1)
                  ui2 = vb%u(jz,ikGPr,2)
                  uei = vb%ue(jz,ikG,1)
                  uei2 = vb%ue(jz,ikGPr,2)
                  tempCmplx = aa*ui2*ui + bb*uei2*uei + ab*ui2*uei + ba*uei2*ui
                  den%vac(jz,1,ivac,3) = den%vac(jz,1,ivac,3) + tempCmplx
               END DO
            ELSE
               !--->                warping part
               DO jz = 1,vacuum%nmzxy
                  ui = vb%u(jz,ikG,1)
                  uj = vb%u(jz,ikGPr,2)
                  uei = vb%ue(jz,ikG,1)
                  uej = vb%ue(jz,ikGPr,2)
                  t1 = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
                  den%vac(jz,ind2,ivac,3) = den%vac(jz, ind2,ivac,3) + t1*phs/stars%nstr2(ind2)
               END DO
            END IF
         END DO
      END DO

   END SUBROUTINE priv_den_offdiag

   SUBROUTINE priv_den_warp_col(vb, vbq, den, stars, vacuum, we, we1_l, ne, ivac, jspin, &
                                l_dfpt, l_center0)
      !! Warping part (G_|| /= 0) of the collinear density.
      !!
      !! In the DFPT case the outer loop runs over the shifted basis at k+q and the
      !! pairing is not symmetric, so the inner loop covers all partners rather than
      !! only ikGPr < ikG. Note that for a shifted star centre (l_center0 false) the
      !! first star is *not* skipped here: priv_den_g0_col does not run in that case,
      !! so the G_|| = 0 term is produced by this loop instead.
      USE m_types
      USE m_types_vacbasis
      IMPLICIT NONE

      TYPE(t_vacbasis),  INTENT(IN)    :: vb, vbq
      TYPE(t_potden),    INTENT(INOUT) :: den
      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      REAL,              INTENT(IN)    :: we(:)
      REAL, ALLOCATABLE, INTENT(IN)    :: we1_l(:)
      INTEGER,           INTENT(IN)    :: ne, ivac, jspin
      LOGICAL,           INTENT(IN)    :: l_dfpt, l_center0

      INTEGER :: ikG, ikGPr, i1, i2, i3, ig3, ind2, ind2p, n, jz, nv2_outer
      REAL    :: ui, uj, uei, uej
      COMPLEX :: aa, bb, ab, ba, phs, phsp
      COMPLEX, ALLOCATABLE :: t1jz(:)

      IF (.NOT.l_dfpt) THEN
         nv2_outer = vb%nv2(jspin)
      ELSE
         nv2_outer = vbq%nv2(jspin)
      END IF

      ! we1_l rather than we1: naming an absent OPTIONAL dummy in a data-sharing
      ! clause is not standard-conforming.
      !$OMP PARALLEL DEFAULT(none) &
      !$OMP SHARED(vb,vbq,nv2_outer,jspin,stars,ne,we,we1_l,vacuum,den,ivac,l_dfpt,l_center0) &
      !$OMP PRIVATE(ikGPr,i1,i2,i3,ig3,phs,phsp,ind2,ind2p,n,jz,ui,uj,uei,uej)&
      !$OMP PRIVATE(aa,bb,ab,ba,t1jz,ikG)
      ALLOCATE(t1jz(vacuum%nmzxy))
      !$OMP DO SCHEDULE(dynamic,5)
      DO ikG = 1, nv2_outer
         IF (.NOT.l_dfpt) THEN
            DO ikGPr = 1, ikG - 1
               i1 = vb%kvac(1,ikG,jspin) - vb%kvac(1,ikGPr,jspin)
               i2 = vb%kvac(2,ikG,jspin) - vb%kvac(2,ikGPr,jspin)
               i3 = 0
               ! stars%ig/i2g are bounded by +-mx1,+-mx2: check before indexing
               IF (iabs(i1).GT.stars%mx1) CYCLE
               IF (iabs(i2).GT.stars%mx2) CYCLE
               ig3 = stars%ig(i1,i2,i3)
               IF (ig3.EQ.0)  CYCLE
               ind2 = stars%i2g(i1,i2)
               ind2p = stars%i2g(-i1,-i2)
               phs = conjg(stars%r2gphs(i1,i2))
               phsp = conjg(stars%r2gphs(-i1,-i2))
               CALL priv_pair_sum(vb%ac(:,:,jspin), vb%bc(:,:,jspin), &
                                  vb%ac(:,:,jspin), vb%bc(:,:,jspin), &
                                  we, ne, ikGPr, ikG, aa, bb, ab, ba)
               DO  jz = 1,vacuum%nmzxy
                  ui = vb%u(jz,ikG,jspin)
                  uj = vb%u(jz,ikGPr,jspin)
                  uei = vb%ue(jz,ikG,jspin)
                  uej = vb%ue(jz,ikGPr,jspin)
                  t1jz(jz) = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
               END DO
               !$OMP CRITICAL ! (denvacxy,denvac)
               den%vac(:vacuum%nmzxy,ind2,ivac,jspin) = den%vac(:vacuum%nmzxy,ind2, ivac,jspin) &
                                                      + t1jz(:vacuum%nmzxy)*phs/stars%nstr2(ind2)
               den%vac(:vacuum%nmzxy,ind2p,ivac,jspin) = den%vac(:vacuum%nmzxy,ind2p,ivac,jspin) &
                                                      + CONJG(t1jz(:vacuum%nmzxy))*phsp/stars%nstr2(ind2p)
               !$OMP END CRITICAL ! (denvacxy,denvac)
            END DO
         ELSE
            DO ikGPr = 1, vb%nv2(jspin)
               i1 = vbq%kvac(1,ikG,jspin) - vb%kvac(1,ikGPr,jspin)
               i2 = vbq%kvac(2,ikG,jspin) - vb%kvac(2,ikGPr,jspin)
               i3 = 0
               ! stars%ig/i2g are bounded by +-mx1,+-mx2: check before indexing
               IF (iabs(i1).GT.stars%mx1) CYCLE
               IF (iabs(i2).GT.stars%mx2) CYCLE
               ig3 = stars%ig(i1,i2,i3)
               IF (ig3.EQ.0) CYCLE
               ind2 = stars%i2g(i1,i2)
               IF ((ind2==1).AND.l_center0) CYCLE
               phs = CONJG(stars%r2gphs(i1,i2)) ! 2D star phase, see above
               aa = 0.0
               bb = 0.0
               ba = 0.0
               ab = 0.0
               DO n = 1,ne
                  aa=aa+we(n)*CONJG(vb%ac(ikGPr,n,jspin))*vbq%ac(ikG,n,jspin)
                  bb=bb+we(n)*CONJG(vb%bc(ikGPr,n,jspin))*vbq%bc(ikG,n,jspin)
                  ab=ab+we(n)*CONJG(vb%ac(ikGPr,n,jspin))*vbq%bc(ikG,n,jspin)
                  ba=ba+we(n)*CONJG(vb%bc(ikGPr,n,jspin))*vbq%ac(ikG,n,jspin)
                  IF (l_center0) THEN
                     aa=aa+we1_l(n)*CONJG(vb%ac(ikGPr,n,jspin))*vb%ac(ikG,n,jspin)
                     bb=bb+we1_l(n)*CONJG(vb%bc(ikGPr,n,jspin))*vb%bc(ikG,n,jspin)
                     ab=ab+we1_l(n)*CONJG(vb%ac(ikGPr,n,jspin))*vb%bc(ikG,n,jspin)
                     ba=ba+we1_l(n)*CONJG(vb%bc(ikGPr,n,jspin))*vb%ac(ikG,n,jspin)
                  END IF
               END DO
               DO  jz = 1,vacuum%nmzxy
                  ui = vbq%u(jz,ikG,jspin)
                  uj = vb%u(jz,ikGPr,jspin)
                  uei = vbq%ue(jz,ikG,jspin)
                  uej = vb%ue(jz,ikGPr,jspin)
                  t1jz(jz) = aa*ui*uj+bb*uei*uej+ba*ui*uej+ab*uei*uj
               END DO
               !$OMP CRITICAL ! (denvacxy,denvac)
               den%vac(:vacuum%nmzxy,ind2,ivac,jspin) = den%vac(:vacuum%nmzxy,ind2, ivac,jspin) &
                                                      + t1jz(:vacuum%nmzxy)*phs/stars%nstr2(ind2)
               !$OMP END CRITICAL ! (denvacxy,denvac)
            END DO
         END IF
      END DO
      !$OMP END DO
      DEALLOCATE(t1jz)
      !$OMP END PARALLEL

   END SUBROUTINE priv_den_warp_col

   SUBROUTINE priv_starcoeff(vb, vacdos, banddos, stars, ne, ev_list, ikpt, ivac, jspin)
      !! Star coefficients 1 .. nstars for each k-point and eigenvalue, written to
      !! vacdos%qstars. The 0th star is the charge integrated over the 2D cell and is
      !! already covered by qvac.
      USE m_types
      USE m_types_vacbasis
      USE m_types_vacdos
      IMPLICIT NONE

      TYPE(t_vacbasis), INTENT(IN)    :: vb
      TYPE(t_vacdos),   INTENT(INOUT) :: vacdos
      TYPE(t_banddos),  INTENT(IN)    :: banddos
      TYPE(t_stars),    INTENT(IN)    :: stars
      INTEGER,          INTENT(IN)    :: ne, ev_list(:), ikpt, ivac, jspin

      INTEGER :: n, ikG, ikGPr, i1, i2, i3, ig3, ind2, ind2p, jj
      REAL    :: ui, uj, uei, uej
      COMPLEX :: aa, bb, ab, ba, t1, phs, phsp

      DO  n=1,ne
         DO ikG = 1,vb%nv2(jspin)
            DO ikGPr = 1, ikG - 1
               i1 = vb%kvac(1,ikG,jspin) - vb%kvac(1,ikGPr,jspin)
               i2 = vb%kvac(2,ikG,jspin) - vb%kvac(2,ikGPr,jspin)
               i3 = 0
               IF (iabs(i1).GT.stars%mx1) CYCLE
               IF (iabs(i2).GT.stars%mx2) CYCLE
               ig3 = stars%ig(i1,i2,i3)
               IF (ig3.EQ.0)  CYCLE
               ind2 = stars%i2g(i1,i2)
               ind2p = stars%i2g(-i1,-i2)
               IF ((ind2.GE.2.AND.ind2.LE.banddos%nstars).OR.&
                  (ind2p.GE.2.AND.ind2p.LE.banddos%nstars)) THEN
                  phs = CONJG(stars%r2gphs(i1,i2))   ! 2D star phase, see above
                  phsp = CONJG(stars%r2gphs(-i1,-i2))
                  aa = CONJG(vb%ac(ikGPr,n,jspin))*vb%ac(ikG,n,jspin)
                  bb = CONJG(vb%bc(ikGPr,n,jspin))*vb%bc(ikG,n,jspin)
                  ab = CONJG(vb%ac(ikGPr,n,jspin))*vb%bc(ikG,n,jspin)
                  ba = CONJG(vb%bc(ikGPr,n,jspin))*vb%ac(ikG,n,jspin)
                  DO jj = 1,banddos%layers
                     ui = vb%u(banddos%izlay(jj,1),ikG,jspin)
                     uj = vb%u(banddos%izlay(jj,1),ikGPr,jspin)
                     uei = vb%ue(banddos%izlay(jj,1),ikG,jspin)
                     uej = vb%ue(banddos%izlay(jj,1),ikGPr,jspin)
                     t1 = aa*ui*uj + bb*uei*uej +ba*ui*uej + ab*uei*uj
                     IF (ind2.GE.2.AND.ind2.LE.banddos%nstars) &
                         vacdos%qstars(ind2-1,ev_list(n),jj,ivac,ikpt,jspin) = vacdos%qstars(ind2-1,ev_list(n),jj,ivac,ikpt,jspin)+ t1*phs/stars%nstr2(ind2)
                     IF (ind2p.GE.2.AND.ind2p.LE.banddos%nstars) &
                         vacdos%qstars(ind2p-1,ev_list(n),jj,ivac,ikpt,jspin) = vacdos%qstars(ind2p-1,ev_list(n),jj,ivac,ikpt,jspin) +CONJG(t1)*phsp/stars%nstr2(ind2p)
                  END DO
               END IF
            END DO
         END DO
      END DO

   END SUBROUTINE priv_starcoeff

END MODULE m_vacden
