!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vgen_finalize
   USE m_juDFT
   USE m_xcBfield
   USE m_plot
   USE m_constants
   USE m_lattHarmsSphHarmsConv

CONTAINS

   SUBROUTINE vgen_finalize(mpi,oneD,field,cell,atoms,stars,vacuum,sym,noco,nococonv,input,xcpot,sphhar,vTot,vCoul,denRot,sliceplot)
      !--------------------------------------------------------------------------
      ! FLAPW potential generator (finalization)
      !
      ! Non-noco: Some rescaling is done here.
      !
      ! Noco: rotate_int_den_from_local is called to generate 2x2 interstitial V matrix.
      !
      ! Fully fully noco: rotate_mt_den_from_local does so for the Muffin Tins.
      !
      ! Sourcefree: The xc-B-field is scaled up an source terms are purged out.
      !--------------------------------------------------------------------------
      USE m_constants
      USE m_rotate_int_den_tofrom_local
      USE m_types
      USE m_rotate_mt_den_tofrom_local
      USE m_magnMomFromDen
      USE m_pw_tofrom_grid

      IMPLICIT NONE

      TYPE(t_mpi),      INTENT(IN)    :: mpi
      TYPE(t_oneD),     INTENT(IN)    :: oneD
      TYPE(t_field),    INTENT(IN)    :: field
      TYPE(t_cell),     INTENT(IN)    :: cell
      TYPE(t_vacuum),   INTENT(IN)    :: vacuum
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_sym),      INTENT(IN)    :: sym
      TYPE(t_stars),    INTENT(IN)    :: stars
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_input),    INTENT(IN)    :: input
      CLASS(t_xcpot),   INTENT(IN)    :: xcpot
      TYPE(t_sphhar),   INTENT(IN)    :: sphhar
      TYPE(t_potden),   INTENT(INOUT) :: vTot, vCoul, denRot
      TYPE(t_sliceplot), INTENT(IN)    :: sliceplot

      TYPE(t_potden)                  :: div, phi, checkdiv, vScal, vCorr, v2, vdiff
      TYPE(t_potden), DIMENSION(3)    :: cvec, corrB, bxc
      TYPE(t_gradients)               :: tmp_grad

      INTEGER                         :: i, js, n, lh, nat, nd, indmax
      REAL                            :: sfscale, r2(atoms%jmtd)
      REAL                            :: b(3,atoms%ntype), dummy1(atoms%ntype), dummy2(atoms%ntype)
      REAL, ALLOCATABLE               :: intden(:,:)

COMPLEX, ALLOCATABLE :: flm(:,:,:,:)

      indmax=(atoms%lmaxd+1)**2

      ALLOCATE(flm(atoms%jmtd,indmax,atoms%ntype,4))

      IF (.NOT.noco%l_noco) THEN
         ! Rescale vTot%pw_w with number of stars:
         DO js=1,SIZE(vtot%pw_w,2)
            DO i=1,stars%ng3
               vTot%pw_w(i,js)=vtot%pw_w(i,js) / stars%nstr(i)
            END DO
         END DO
      ELSE IF(noco%l_noco) THEN
         ! Rotate interstital potential back to global frame:
         CALL rotate_int_den_from_local(stars,atoms,vacuum,sym,input,denRot,vTot)
         IF (noco%l_mtnocoPot) THEN
            ! Rotate Muffin Tin potential back to global frame:
            CALL rotate_mt_den_from_local(atoms,sphhar,sym,denRot,noco,vtot)
         END IF
      END IF

      IF (noco%l_mtnocoPot.AND.noco%l_scaleMag) THEN
         sfscale=noco%mag_scale
         CALL vTot%SpinsToChargeAndMagnetisation()
         vTot%mt(:,0:,:,  2:4) = sfscale*vTot%mt(:,0:,:,2:4)
         vTot%pw(:,       2:3) = sfscale*vTot%pw(:,     2:3)
         vTot%vacz(:,:,   2:4) = sfscale*vTot%vacz(:,:, 2:4)
         vTot%vacxy(:,:,:,2:3) = sfscale*vTot%vacxy(:,:,:,2:3)
         CALL vTot%ChargeAndMagnetisationToSpins()
      END IF

      IF (noco%l_mtnocoPot.AND.noco%l_sourceFree) THEN

         CALL magnMomFromDen(input,atoms,noco,vTot,b,dummy1,dummy2)
         DO i=1,atoms%ntype
            WRITE  (6,8025) i,b(1,i),b(2,i),b(3,i),SQRT(b(1,i)**2+b(2,i)**2+b(3,i)**2)
            8025 FORMAT(2x,'Bfield before SF [local frame, atom ',i2,']: ','Bx=',f9.5,' By=',f9.5,' Bz=',f9.5,' |B|=',f9.5)
         END DO

         CALL timestart("Purging source terms in B-field")

         CALL timestart("Building B")
         CALL makeVectorField(sym,stars,atoms,sphhar,vacuum,input,noco,nococonv,vTot,2.0,vScal,bxc,cell)
         CALL timestop("Building B")

         CALL timestart("SF subroutine")
         CALL sourcefree(mpi,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,bxc,vScal,div,phi,vCorr,cvec,corrB,checkdiv)
         CALL timestop("SF subroutine")

         !CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, nococonv, &
         !         .FALSE., .FALSE., 'div                 ', div)

         !CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, nococonv, &
         !            .FALSE., .TRUE., 'phiDiv              ', phi)

         !CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, nococonv, &
         !            .FALSE., .FALSE., 'gradPhiDiv          ', cvec(1), cvec(1), cvec(2), cvec(3))

         !CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, nococonv, &
         !            .FALSE., .FALSE., 'bCorrected          ', corrB(1), corrB(1), corrB(2), corrB(3))

         CALL div%resetPotDen()
         CALL checkdiv%resetPotDen()
         CALL phi%resetPotDen()

         DO i=1,3
            CALL corrB(i)%resetPotDen()
         END DO

         CALL timestart("Correcting vTot")

         !ALLOCATE (vCorr%pw_w,  mold=vCorr%pw)
         !vTot%pw_w=CMPLX(0.0,0.0)
         !vCorr%pw_w=CMPLX(0.0,0.0)

         !CALL correctPot(vTot,cvec)

         DO js=1,4
            DO i=1,atoms%ntype
               DO lh=0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:i - 1)) + 1))
                  r2=atoms%rmsh(:,i)**2
                  vCorr%mt(:,lh,i,js) = vCorr%mt(:,lh,i,js)/r2
               END DO !lh
            END DO !i
         END DO !js

         !CALL init_pw_grid(.FALSE.,stars,sym,cell)
         !CALL pw_to_grid(.FALSE.,3,.FALSE.,stars,cell,vCorr%pw,tmp_grad,rho=intden)
         !CALL pw_from_grid(.FALSE.,stars,.TRUE.,intden,vCorr%pw,vCorr%pw_w)
         !intden=0.0
         !CALL pw_to_grid(.FALSE.,3,.FALSE.,stars,cell,vTot%pw,tmp_grad,rho=intden)
         !CALL pw_from_grid(.FALSE.,stars,.TRUE.,intden,vTot%pw,vTot%pw_w)
         !CALL finish_pw_grid()

         vTot=vCorr

         CALL timestop("Correcting vTot")

         CALL timestop("Purging source terms in B-field")

         CALL magnMomFromDen(input,atoms,noco,vTot,b,dummy1,dummy2)
         DO i=1,atoms%ntype
            WRITE  (6,8026) i,b(1,i),b(2,i),b(3,i),SQRT(b(1,i)**2+b(2,i)**2+b(3,i)**2)
            8026 FORMAT(2x,'Bfield after SF [local frame, atom ',i2,']: ','Bx=',f9.5,' By=',f9.5,' Bz=',f9.5,' |B|=',f9.5)
         END DO
      END IF

      IF ((sliceplot%iplot.NE.0 )) THEN
         CALL makeplots(stars, atoms, sphhar, vacuum, input, mpi,oneD, sym, cell, &
                        noco,nococonv, vTot, PLOT_POT_TOT, sliceplot)
         !CALL makeplots(fi%sym,stars,fi%vacuum,fi%atoms,sphhar,fi%input,fi%cell,fi%oneD,fi%noco,fi%sliceplot,vCoul,PLOT_POT_COU)
         !CALL subPotDen(vxcForPlotting,vTot,vCoul)
         !CALL makeplots(fi%sym,stars,fi%vacuum,fi%atoms,sphhar,fi%input,fi%cell,fi%oneD,fi%noco,fi%sliceplot,vxcForPlotting,PLOT_POT_VXC)
      END IF

      ! Store vTot(L=0) component as r*vTot(L=0)/sqrt(4*pi):
      ! (Used input%jspins instead of SIZE(vtot%mt,4) since
      ! the off diagonal part of VTot%mt needs no rescaling!)
      DO js = 1, input%jspins
         DO n = 1, atoms%ntype
            vTot%mt(:atoms%jri(n),0,n,js)  = atoms%rmsh(:atoms%jri(n),n)*vTot%mt(:atoms%jri(n),0,n,js)/sfp_const
         END DO ! n
      END DO ! js

      ! Rescale vCoul%pw_w with number of stars:
      ! (This normalization is needed for gw!)
      DO js = 1, SIZE(vCoul%pw_w,2)
         DO i = 1, stars%ng3
            vcoul%pw_w(i,js) = vcoul%pw_w(i,js) / stars%nstr(i)
         END DO
      END DO

      ! Copy first vacuum into second vacuum if this was not calculated before:
      IF (vacuum%nvac==1) THEN
         vTot%vacz(:,2,:)  = vTot%vacz(:,1,:)
         IF (sym%invs) THEN
            vTot%vacxy(:,:,2,:)  = CMPLX(vTot%vacxy(:,:,1,:))
         ELSE
            vTot%vacxy(:,:,2,:)  = vTot%vacxy(:,:,1,:)
         END IF
      END IF

   END SUBROUTINE vgen_finalize

END MODULE m_vgen_finalize
