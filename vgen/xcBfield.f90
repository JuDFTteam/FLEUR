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
   USE m_juDFT

   IMPLICIT NONE

   !-----------------------------------------------------------------------------
   ! This module contains all the operations on exchange-correlation B-fields
   ! that are necessary to project out source terms. This way, the whole modifi-
   ! cation towards source-free fields can be done by one call, either as a post-
   ! process test or in the scf-loop to achieve said fields self-consistently.
   !-----------------------------------------------------------------------------

   PUBLIC :: makeVectorField, sourcefree, correctPot

CONTAINS
   SUBROUTINE makeVectorField(sym,stars,atoms,sphhar,vacuum,input,noco,nococonv,denmat,factor,vScal,aVec,cell)
      
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
      TYPE(t_nococonv),             INTENT(IN)  :: nococonv
      TYPE(t_potden),               INTENT(IN)  :: denmat
      REAL,                         INTENT(IN)  :: factor
      TYPE(t_potden), DIMENSION(3), INTENT(OUT) :: aVec
      TYPE(t_potden),               INTENT(OUT) :: vScal
      TYPE(t_cell),                 INTENT(IN)  :: cell

      TYPE(t_gradients)                         :: tmp_grad
      TYPE(t_potden), DIMENSION(4)              :: dummyDen
      INTEGER                                   :: i, itype, ir, lh
      REAL                                      :: r2(atoms%jmtd), fcut
      REAL, ALLOCATABLE                         :: intden(:,:)

      fcut=1.e-12

      ! Initialize and fill a dummy density array, that takes the initial result
      ! of matrixsplit.
      DO i=1, 4
         CALL dummyDen(i)%init_potden_simple(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd, &
                          atoms%ntype,atoms%n_u,atoms%n_vPairs,1,.FALSE.,.FALSE., &
                          POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(dummyDen(i)%pw_w,mold=dummyDen(i)%pw)
      ENDDO

      CALL matrixsplit(sym,stars,atoms,sphhar,vacuum,input,noco,nococonv,factor,denmat, &
                       dummyDen(1),dummyDen(2),dummyDen(3),dummyDen(4))

      vScal=dummyDen(1)

      r2=1.0

      ! Initialize and fill the vector field.
      DO i=1,3
         CALL aVec(i)%init_potden_simple(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd, &
                     atoms%ntype,atoms%n_u,atoms%n_vPairs,1,.FALSE.,.FALSE., &
                     POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(aVec(i)%pw_w,mold=aVec(i)%pw)
         aVec(i)%mt(:,:,:,:) = dummyDen(i+1)%mt(:,:,:,:)
         DO itype=1,atoms%ntype
            DO lh=0, sphhar%nlh(sym%ntypsy(atoms%firstAtom(itype)))
               IF (factor==2.0) THEN
                  r2=atoms%rmsh(:,itype)**2
               END IF
               IF (i==1) THEN
                  vScal%mt(:,lh,itype,1) = vScal%mt(:,lh,itype,1)*r2
               END IF
               aVec(i)%mt(:,lh,itype,1) = aVec(i)%mt(:,lh,itype,1)*r2
            END DO !lh
         END DO !itype
         aVec(i)%pw(1:,:)          = dummyDen(i+1)%pw(1:,:)
         aVec(i)%pw_w(1:,:)        = dummyDen(i+1)%pw_w(1:,:)
         !aVec(i)%vacz(1:,1:,:)     = dummyDen(i+1)%vacz(1:,1:,:)
         !aVec(i)%vacxy(1:,1:,1:,:) = dummyDen(i+1)%vacxy(1:,1:,1:,:)
         aVec(i)%vac(1:,1:,1:,:) = dummyDen(i+1)%vac(1:,1:,1:,:)
      END DO !i

   END SUBROUTINE makeVectorField

   SUBROUTINE sourcefree(fmpi,field,stars,atoms,sphhar,vacuum,input ,sym,cell,noco,aVec,vScal,vCorr)
      USE m_vgen_coulomb
      USE m_gradYlm
      USE m_grdchlh
      USE m_sphpts
      USE m_checkdop
      USE m_BfieldtoVmat
      USE m_pw_tofrom_grid

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

      TYPE(t_mpi),                  INTENT(IN)     :: fmpi
      TYPE(t_field),                INTENT(IN)     :: field
      TYPE(t_stars),                INTENT(IN)     :: stars
      TYPE(t_atoms),                INTENT(IN)     :: atoms
      TYPE(t_sphhar),               INTENT(IN)     :: sphhar
      TYPE(t_vacuum),               INTENT(IN)     :: vacuum
      TYPE(t_input),                INTENT(IN)     :: input
       
      TYPE(t_sym),                  INTENT(IN)     :: sym
      TYPE(t_cell),                 INTENT(IN)     :: cell
      TYPE(t_noco),                 INTENT(IN)     :: noco
      TYPE(t_potden), DIMENSION(3), INTENT(INOUT)  :: aVec
      TYPE(t_potden),               INTENT(IN)     :: vScal
      TYPE(t_potden),               INTENT(OUT)    :: vCorr

      TYPE(t_potden)                               :: div, phi, divloc
      TYPE(t_potden), DIMENSION(3)                 :: cvec, corrB
      TYPE(t_atoms)                                :: atloc
      INTEGER                                      :: n, jr, lh, lhmax, jcut, nat ,i
      REAL                                         :: xp(3,(atoms%lmaxd+1+mod(atoms%lmaxd+1,2))*(2*atoms%lmaxd+1))
      REAL, ALLOCATABLE                            :: intden(:,:)
      TYPE(t_gradients)               :: tmp_grad
      complex                          :: sigma_loc(2)

      CALL div%init_potden_simple(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd,atoms%ntype, &
                                  atoms%n_u,atoms%n_vPairs,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN, &
                                  vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(div%pw_w,mold=div%pw)
      div%pw_w = CMPLX(0.0,0.0)

      CALL timestart("Building divergence")
      CALL divergence(input,stars,atoms,sphhar,vacuum,sym,cell,noco,aVec,div)
      CALL timestop("Building divergence")

      ! Local atoms variable with no charges;
      ! needed for the potential generation from the divergence.
      atloc=atoms
      atloc%zatom=0.0
      CALL phi%init_potden_simple(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd,atoms%ntype,atoms%n_u,atoms%n_vPairs,1,.FALSE.,.FALSE.,POTDEN_TYPE_POTCOUL,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(phi%pw_w(SIZE(phi%pw,1),size(phi%pw,2)))
      phi%pw_w = CMPLX(0.0,0.0)

      CALL timestart("Building potential")
      sigma_loc = cmplx(0.0,0.0)
      CALL vgen_coulomb(1,fmpi ,input,field,vacuum,sym,stars,cell,sphhar,atloc,.TRUE.,div,phi,sigma_loc)
      CALL timestop("Building potential")

      DO i=1,3
         CALL cvec(i)%init_potden_simple(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd,atoms%ntype,atoms%n_u,atoms%n_vPairs,1,.FALSE.,.FALSE.,POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(cvec(i)%pw_w,mold=cvec(i)%pw)
         cvec(i)%pw_w=CMPLX(0.0,0.0)
      ENDDO

      CALL timestart("Building correction field")
      CALL divpotgrad(input,stars,atloc,sphhar,vacuum,sym,cell,noco,phi,cvec)
      CALL timestop("Building correction field")

      CALL init_pw_grid(stars,sym,cell)
      DO i=1,3
         CALL pw_to_grid(.FALSE.,1,.FALSE.,stars,cell,cvec(i)%pw,tmp_grad,rho=intden)
         cvec(i)%pw=CMPLX(0.0,0.0)
         cvec(i)%pw_w=CMPLX(0.0,0.0)
         CALL pw_from_grid(stars,intden,cvec(i)%pw,cvec(i)%pw_w)
         DEALLOCATE(intden)
      END DO !i
      CALL finish_pw_grid()

      DO i=1,3
         CALL corrB(i)%init_potden_simple(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd,atoms%ntype,atoms%n_u,atoms%n_vPairs,1,.FALSE.,.FALSE.,&
                                          POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(corrB(i)%pw_w,mold=corrB(i)%pw)
         corrB(i)%pw_w=CMPLX(0.0,0.0)
         CALL corrB(i)%addPotDen(aVec(i),cvec(i))
      ENDDO

      CALL BfieldtoVmat(sym, stars, atoms, sphhar, vacuum, vScal, corrB(1), corrB(2), corrB(3), vCorr)

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

      !vtot%mt(:,0:,:,1)=(vtot%mt(:,0:,:,1)+vtot%mt(:,0:,:,2))/2
      !vtot%mt(:,0:,:,2)=vtot%mt(:,0:,:,1)
      !vtot%mt(:,0:,:,3:4)=0.0

      !vTot%pw(1:,1)=(vTot%pw(1:,1)+vTot%pw(1:,2))/2
      !vTot%pw(1:,2)=vTot%pw(1:,1)
      !vTot%pw(1:,3)=CMPLX(0.0,0.0)

      !vTot%pw_w(1:,1)=(vTot%pw_w(1:,1)+vTot%pw_w(1:,2))/2
      !vTot%pw_w(1:,2)=vTot%pw_w(1:,1)
      !vTot%pw_w(1:,3)=CMPLX(0.0,0.0)

      !vTot%mt(:,0:,:,1)=vTot%mt(:,0:,:,1)+b(3)%mt(:,0:,:,1)/2
      !vTot%mt(:,0:,:,2)=vTot%mt(:,0:,:,2)-b(3)%mt(:,0:,:,1)/2
      !vTot%mt(:,0:,:,3)=b(1)%mt(:,0:,:,1)/2
      !vTot%mt(:,0:,:,4)=b(2)%mt(:,0:,:,1)/2

      !vTot%pw(1:,1)=vTot%pw(1:,1)+b(3)%pw(1:,1)/2
      !vTot%pw(1:,2)=vTot%pw(1:,2)-b(3)%pw(1:,1)/2
      !vTot%pw(1:,3)=(b(1)%pw(1:,1)+ImagUnit*b(2)%pw(1:,1))/2

      !vTot%pw_w(1:,1)=vTot%pw_w(1:,1)+b(3)%pw_w(1:,1)/2
      !vTot%pw_w(1:,2)=vTot%pw_w(1:,2)-b(3)%pw_w(1:,1)/2
      !vTot%pw_w(1:,3)=(b(1)%pw_w(1:,1)+ImagUnit*b(2)%pw_w(1:,1))/2

      vTot%mt(:,0:,:,1)=c(3)%mt(:,0:,:,1)+vTot%mt(:,0:,:,1)
      vTot%mt(:,0:,:,2)=-c(3)%mt(:,0:,:,1)+vTot%mt(:,0:,:,2)
      vTot%mt(:,0:,:,3)=c(1)%mt(:,0:,:,1)+vTot%mt(:,0:,:,3)
      vTot%mt(:,0:,:,4)=c(2)%mt(:,0:,:,1)+vTot%mt(:,0:,:,4)

      vTot%pw(1:,1)=vTot%pw(1:,1)+c(3)%pw(1:,1)
      vTot%pw(1:,2)=vTot%pw(1:,2)-c(3)%pw(1:,1)
      pwr=REAL(vTot%pw(1:,3)) + REAL(c(1)%pw(1:,1)) - AIMAG(c(2)%pw(1:,1))
      pwi=AIMAG(vTot%pw(1:,3)) + AIMAG(c(1)%pw(1:,1)) + REAL(c(2)%pw(1:,1))
      vTot%pw(1:,3)=CMPLX(pwr,pwi)

      !vTot%pw_w(1:,1)=vTot%pw_w(1:,1)+c(3)%pw_w(1:,1)
      !vTot%pw_w(1:,2)=vTot%pw_w(1:,2)-c(3)%pw_w(1:,1)
      !pwr=REAL(vTot%pw_w(1:,3)) + REAL(c(1)%pw_w(1:,1)) - AIMAG(c(2)%pw_w(1:,1))
      !pwi=AIMAG(vTot%pw_w(1:,3)) + AIMAG(c(1)%pw_w(1:,1)) + REAL(c(2)%pw_w(1:,1))
      !vTot%pw_w(1:,3)=CMPLX(pwr,pwi)

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
