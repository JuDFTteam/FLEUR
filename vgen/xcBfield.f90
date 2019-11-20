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
   
   !-----------------------------------------------------------------------------
   ! When finished, this module will contain all the operations on exchange-cor-
   ! relation B-fields, that are currenlty done in fleur.F90 after the scf-loop.
   ! This way, the whole modification towards source-free fields can be done by
   ! one call, either as a postprocess test or in the scf-loop to achieve self-
   ! consistent source-free fields.
   !-----------------------------------------------------------------------------

   PUBLIC :: makeBxc, sourcefree ! TODO: Build a routine to pack A_vec 
                                 !       back into the matrix correctly. 

CONTAINS
   SUBROUTINE makeVectorField(stars,atoms,sphhar,vacuum,input,noco,denmat,factor,aVec)
      
      ! Contructs the exchange-correlation magnetic field from the total poten-
      ! tial matrix or the magnetic density for the density matrix. Only used for
      ! the implementation of source free fields and therefore only applicable in
      ! a fully non-collinear description of magnetism.
      ! 
      ! Assumes the following form for the density/potential matrix:
      ! rho_mat = n*Id_(2x2) + conj(sigma_vec)*m_vec
      ! V_mat   = V*Id_(2x2) + sigma_vec*B_vec
      ! 
      ! A_vec is saved as a density type with an additional r^2-factor.

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
      INTEGER                                   :: i, itype, ir
      REAL                                      :: r2
      
      ! Initialize and fill a dummy density array, that takes the initial result
      ! of matrixsplit.
      DO i=1, 4
         CALL dummyDen(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd, &
                          atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE., &
                          POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(dummyDen(i)%pw_w,mold=dummyDen(i)%pw)
      ENDDO

      CALL matrixsplit(stars,atoms,sphhar,vacuum,input,noco,factor,denmat, &
                       dummyDen(1),dummyDen(2),dummyDen(3),dummyDen(4))

      r2=1.0
      
      ! Initialize and fill the vector field.
      DO i=1,3
         CALL aVec(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd, &
                     atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE., &
                     POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(aVec(i)%pw_w,mold=aVec(i)%pw)
         DO itype=1,atoms%ntype
            DO ir=1,atoms%jri(itype)
               IF (factor==2.0) THEN 
                  r2=atoms%rmsh(ir,itype)**2
               END IF 
               aVec(i)%mt(ir,0:,itype,:) = dummyDen(i+1)%mt(ir,0:,itype,:)*r2
            END DO !ir
         END DO !itype
         aVec(i)%pw(1:,:)          = dummyDen(i+1)%pw(1:,:)
         aVec(i)%vacz(1:,1:,:)     = dummyDen(i+1)%vacz(1:,1:,:)
         aVec(i)%vacxy(1:,1:,1:,:) = dummyDen(i+1)%vacxy(1:,1:,1:,:)
      END DO !i
     
   END SUBROUTINE makeVectorField

   SUBROUTINE sourcefree(mpi,dimension,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,bxc,div,phi,cvec,corrB,checkdiv)
      USE m_vgen_coulomb

      ! Takes a vectorial quantity, i.e. a t_potden variable of dimension 3, and
      ! makes it into a source free vector field as follows:
      ! 
      ! a) Build the divergence d of the vector field A_vec as d=nabla_vec*A_vec.
      ! b) Solve the Poisson equation (nabla_vec*nabla_vec)phi=-4*pi*d.
      ! c) Construct an auxiliary vector field C_vec=(nabla_vec phi)/(4*pi).
      ! d) Build A_vec_sf=A_vec+C_vec, which is source free by construction.
      !       
      ! Note: A_vec is assumed to be a density with an additional factor of r^2.
      !       A field of the same form will also be calculated.  

      TYPE(t_mpi),                  INTENT(IN)     :: mpi
      TYPE(t_dimension),            INTENT(IN)     :: dimension
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
      TYPE(t_potden), DIMENSION(3), INTENT(INOUT)  :: bxc
      TYPE(t_potden),               INTENT(OUT)    :: div, phi, checkdiv
      TYPE(t_potden), DIMENSION(3), INTENT(OUT)    :: cvec, corrB
      
      TYPE(t_atoms)                                :: atloc
      INTEGER                                      :: n, jr, lh, lhmax

      CALL div%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype, &
                                  atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN, &
                                  vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(div%pw_w,mold=div%pw)
      
      !CALL divergence(stars,atoms,sphhar,vacuum,sym,cell,noco,bxc,div)
      CALL divergence2(stars,atoms,sphhar,vacuum,sym,cell,noco,bxc,div)
 
      atloc=atoms
      atloc%zatom=0.0 !Local atoms variable with no charges; needed for the potential generation.
      eps=1.e-10

      CALL phi%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_POTCOUL,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(phi%pw_w(SIZE(phi%pw,1),size(phi%pw,2)))
      phi%pw_w = CMPLX(0.0,0.0)

      CALL vgen_coulomb(1,mpi,dimension,oneD,input,field,vacuum,sym,stars,cell,sphhar,atloc,.TRUE.,div,phi)
   
      !DO n=1,atoms%ntype   
       !  lhmax=sphhar%nlh(atoms%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
        ! DO jr = 1, atoms%jri(n)
         !   DO lh=0, lhmax
          !     !IF (ABS(phi%mt(jr,lh,n,1))<eps) THEN
           !    IF (lh/=1) THEN
            !      phi%mt(jr,lh,n,:)=0.0
             !  END IF
            !END DO
        ! END DO
     ! END DO

      DO i=1,3
         CALL cvec(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(cvec(i)%pw_w,mold=cvec(i)%pw)
      ENDDO

      !CALL divpotgrad(stars,atoms,sphhar,vacuum,sym,cell,noco,phi,cvec)
      CALL divpotgrad2(stars,atoms,sphhar,vacuum,sym,cell,noco,phi,cvec)

      DO i=1,3
         CALL corrB(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(corrB(i)%pw_w,mold=corrB(i)%pw)
         CALL corrB(i)%addPotDen(bxc(i),cvec(i))
      ENDDO
      
      CALL checkdiv%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype, &
                                  atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN, &
                                  vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(checkdiv%pw_w,mold=checkdiv%pw)
      
      CALL divergence2(stars,atoms,sphhar,vacuum,sym,cell,noco,corrB,checkdiv)

      !checkdiv%mt(:,2:,:,:)=0.0
      !checkdiv%mt(:,0,:,:)=0.0

   END SUBROUTINE sourcefree

END MODULE m_xcBfield
