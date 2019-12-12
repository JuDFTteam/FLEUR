!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_xcBfield
   USE m_types
   USE m_constants
   USE m_plot
   USE m_divergence

   IMPLICIT NONE

   !-----------------------------------------------------------------------------
   ! This module contains all the operations on exchange-correlation B-fields
   ! that are necessary to project out source terms. This way, the whole modifi-
   ! cation towards source-free fields can be done by one call, either as a post-
   ! process test or in the scf-loop to achieve said fields self-consistently.
   !-----------------------------------------------------------------------------

   PUBLIC :: makeVectorField, sourcefree, correctPot

CONTAINS
   SUBROUTINE makeVectorField(sym,stars,atoms,sphhar,vacuum,input,noco,denmat,factor,aVec,icut)

      ! Contructs the exchange-correlation magnetic field from the total poten-
      ! tial matrix or the magnetic density for the density matrix. Only used for
      ! the implementation of source free fields and therefore only applicable in
      ! a (fully) non-collinear description of magnetism.
      !
      ! Assumes the following form for the density/potential matrix:
      ! rho_mat = n*Id_(2x2) + conj(sigma_vec)*m_vec
      ! V_mat   = V*Id_(2x2) + sigma_vec*B_vec
      !
      ! A_vec is saved as a density type with an additional r^2-factor.
      !
      ! TODO: How do we constuct B when not only it is saved to the density
      !       matrix (SOC, LDA+U etc.)?
      TYPE(t_sym),                  INTENT(IN)  :: sym
      TYPE(t_stars),                INTENT(IN)  :: stars
      TYPE(t_atoms),                INTENT(IN)  :: atoms
      TYPE(t_sphhar),               INTENT(IN)  :: sphhar
      TYPE(t_vacuum),               INTENT(IN)  :: vacuum
      TYPE(t_input),                INTENT(IN)  :: input
      TYPE(t_noco),                 INTENT(IN)  :: noco
      TYPE(t_potden),               INTENT(IN)  :: denmat
      REAL,                         INTENT(IN)  :: factor
      TYPE(t_potden), DIMENSION(3), INTENT(OUT) :: aVec

      TYPE(t_potden), DIMENSION(4)              :: dummyDen
      INTEGER                                   :: i, itype, ir, lh
      REAL                                      :: r2(atoms%jmtd), fcut

      fcut=1.e-12

      ! Initialize and fill a dummy density array, that takes the initial result
      ! of matrixsplit.
      DO i=1, 4
         CALL dummyDen(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd, &
                          atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE., &
                          POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(dummyDen(i)%pw_w,mold=dummyDen(i)%pw)
      ENDDO

      CALL matrixsplit(sym,stars,atoms,sphhar,vacuum,input,noco,factor,denmat, &
                       dummyDen(1),dummyDen(2),dummyDen(3),dummyDen(4))

      r2=1.0

      ! Initialize and fill the vector field.
      DO i=1,3
         CALL aVec(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd, &
                     atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE., &
                     POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(aVec(i)%pw_w,mold=aVec(i)%pw)
         aVec(i)%mt(:,:,:,:) = dummyDen(i+1)%mt(:,:,:,:)
         DO itype=1,atoms%ntype
            DO lh=0, sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:itype - 1)) + 1))
               IF (factor==2.0) THEN
                  r2=atoms%rmsh(:,itype)**2
               END IF
               aVec(i)%mt(:,lh,itype,1) = aVec(i)%mt(:,lh,itype,1)*r2
            END DO !lh
         END DO !itype
         aVec(i)%pw(1:,:)          = dummyDen(i+1)%pw(1:,:)
         aVec(i)%vacz(1:,1:,:)     = dummyDen(i+1)%vacz(1:,1:,:)
         aVec(i)%vacxy(1:,1:,1:,:) = dummyDen(i+1)%vacxy(1:,1:,1:,:)
      END DO !i

   END SUBROUTINE makeVectorField

   SUBROUTINE sourcefree(mpi,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,aVec,icut,div,phi,cvec,corrB,checkdiv)
      USE m_vgen_coulomb
      USE m_gradYlm
      USE m_grdchlh
      USE m_sphpts
      USE m_checkdop

      ! Takes a vectorial quantity, i.e. a t_potden variable of dimension 3, and
      ! makes it into a source free vector field as follows:
      !
      ! a) Build the divergence d of the vector field A_vec as d=nabla_vec*A_vec.
      ! b) Solve the Poisson equation (nabla_vec*nabla_vec)phi=-4*pi*d for phi.
      ! c) Construct an auxiliary vector field C_vec=(nabla_vec phi)/(4*pi).
      ! d) Build A_vec_sf=A_vec+C_vec, which is source free by construction.
      !
      ! Note: A_vec is assumed to be a density with an additional factor of r^2.
      !       A field of the same form will also be calculated.

      TYPE(t_mpi),                  INTENT(IN)     :: mpi
      TYPE(t_field),                INTENT(INOUT)  :: field
      TYPE(t_stars),                INTENT(IN)     :: stars
      TYPE(t_atoms),                INTENT(IN)     :: atoms
      TYPE(t_sphhar),               INTENT(IN)     :: sphhar
      TYPE(t_vacuum),               INTENT(IN)     :: vacuum
      TYPE(t_input),                INTENT(IN)     :: input
      TYPE(t_oneD),                 INTENT(IN)     :: oneD
      TYPE(t_sym),                  INTENT(IN)     :: sym
      TYPE(t_cell),                 INTENT(IN)     :: cell
      TYPE(t_noco),                 INTENT(IN)     :: noco
      TYPE(t_potden), DIMENSION(3), INTENT(INOUT)  :: aVec
      TYPE(t_potden),               INTENT(OUT)    :: div, phi, checkdiv
      TYPE(t_potden), DIMENSION(3), INTENT(OUT)    :: cvec, corrB

      TYPE(t_potden)                               :: divloc
      TYPE(t_atoms)                                :: atloc
      INTEGER                                      :: n, jr, lh, lhmax, jcut, nat
      REAL                                         :: xp(3,(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1))

      CALL div%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype, &
                                  atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN, &
                                  vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(div%pw_w,mold=div%pw)

      CALL divergence2(stars,atoms,sphhar,vacuum,sym,cell,noco,aVec,div)

      ! Local atoms variable with no charges;
      ! needed for the potential generation from the divergence.
      atloc=atoms
      atloc%zatom=0.0
      CALL phi%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_POTCOUL,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(phi%pw_w(SIZE(phi%pw,1),size(phi%pw,2)))
      phi%pw_w = CMPLX(0.0,0.0)

      CALL vgen_coulomb(1,mpi,oneD,input,field,vacuum,sym,stars,cell,sphhar,atloc,.TRUE.,div,phi)

      DO i=1,3
         CALL cvec(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(cvec(i)%pw_w,mold=cvec(i)%pw)
      ENDDO

      CALL divpotgrad2(stars,atoms,sphhar,vacuum,sym,cell,noco,phi,cvec)

      DO i=1,3
         CALL corrB(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(corrB(i)%pw_w,mold=corrB(i)%pw)
         CALL corrB(i)%addPotDen(aVec(i),cvec(i))
      ENDDO

      CALL checkdiv%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype, &
                                  atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN, &
                                  vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(checkdiv%pw_w,mold=checkdiv%pw)

      CALL divergence2(stars,atoms,sphhar,vacuum,sym,cell,noco,corrB,checkdiv)

      DO n=1,atoms%ntype
         lhmax=sphhar%nlh(sym%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
         DO lh=0, lhmax
            cvec(1)%mt(:,lh,n,1)=cvec(1)%mt(:,lh,n,1)/(atoms%rmsh(:, n)**2)
            cvec(2)%mt(:,lh,n,1)=cvec(2)%mt(:,lh,n,1)/(atoms%rmsh(:, n)**2)
            cvec(3)%mt(:,lh,n,1)=cvec(3)%mt(:,lh,n,1)/(atoms%rmsh(:, n)**2)
         END DO
      END DO

   END SUBROUTINE sourcefree

   SUBROUTINE correctPot(vTot,c)
      USE m_types

      ! Takes a vectorial quantity c and saves its components into the appro
      ! priate components of the potential matrix V.
      !
      ! An initial V_mat = V*Id_(2x2) + sigma_vec*B_vec will become
      !            V_mat = V*Id_(2x2) + sigma_vec*(B_vec+C_vec)
      !
      ! TODO: Both quantities are assumed to be in the global frame of refe-
      !       rence. Make sure this is true. Also: consider SOC, LDA+U etc.

      TYPE(t_potden),               INTENT(INOUT) :: vTot
      TYPE(t_potden), DIMENSION(3), INTENT(IN)    :: c

      REAL :: pwr(SIZE(vTot%pw(1:,3))), pwi(SIZE(vTot%pw(1:,3)))

      vTot%mt(:,0:,:,1)=vTot%mt(:,0:,:,1)+c(3)%mt(:,0:,:,1)
      vTot%mt(:,0:,:,2)=vTot%mt(:,0:,:,2)-c(3)%mt(:,0:,:,1)
      vTot%mt(:,0:,:,3)=vTot%mt(:,0:,:,3)+c(1)%mt(:,0:,:,1)
      vTot%mt(:,0:,:,4)=vTot%mt(:,0:,:,4)+c(2)%mt(:,0:,:,1)

      vTot%pw(1:,1)=vTot%pw(1:,1)+c(3)%pw(1:,1)
      vTot%pw(1:,2)=vTot%pw(1:,2)-c(3)%pw(1:,1)
      pwr=REAL(vTot%pw(1:,3)) + REAL(c(1)%pw(1:,1)) - AIMAG(c(2)%pw(1:,1))
      pwi=AIMAG(vTot%pw(1:,3)) + AIMAG(c(1)%pw(1:,1)) + REAL(c(2)%pw(1:,1))
      vTot%pw(1:,3)=CMPLX(pwr,pwi)

      vTot%pw_w(1:,1)=vTot%pw_w(1:,1)+c(3)%pw_w(1:,1)
      vTot%pw_w(1:,2)=vTot%pw_w(1:,2)-c(3)%pw_w(1:,1)
      pwr=REAL(vTot%pw_w(1:,3)) + REAL(c(1)%pw_w(1:,1)) - AIMAG(c(2)%pw_w(1:,1))
      pwi=AIMAG(vTot%pw_w(1:,3)) + AIMAG(c(1)%pw_w(1:,1)) + REAL(c(2)%pw_w(1:,1))
      vTot%pw_w(1:,3)=CMPLX(pwr,pwi)

      !vTot%vacz(1:,1:,1)=vTot%vacz(1:,1:,1)+c(3)%vacz(1:,1:,1)
      !vTot%vacz(1:,1:,2)=vTot%vacz(1:,1:,2)-c(3)%vacz(1:,1:,1)
      !vTot%vacz(1:,1:,3)=vTot%vacz(1:,1:,3)+c(1)%vacz(1:,1:,1)
      !vTot%vacz(1:,1:,4)=vTot%vacz(1:,1:,4)+c(2)%vacz(1:,1:,1)

      !vTot%vacxy(1:,1:,1:,1)=vTot%vacxy(1:,1:,1:,1)+c(3)%vacxy(1:,1:,1:,1)
      !vTot%vacxy(1:,1:,1:,2)=vTot%vacxy(1:,1:,1:,2)-c(3)%vacxy(1:,1:,1:,1)
      !vTot%vacxy(1:,1:,1:,3)=vTot%vacxy(1:,1:,1:,3)+c(1)%vacxy(1:,1:,1:,1)
      !vTot%vacxy(1:,1:,1:,4)=vTot%vacxy(1:,1:,1:,4)+c(2)%vacxy(1:,1:,1:,1)

   END SUBROUTINE correctPot

END MODULE m_xcBfield
