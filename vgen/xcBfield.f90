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

   PUBLIC :: makeBxc, sourcefree

CONTAINS
   SUBROUTINE makeBxc(stars,atoms,sphhar,vacuum,input,noco,vTot,bxc)
      
      ! Contructs the exchange-correlation magnetic field from the total poten-
      ! tial matrix. Only used for the implementation of source free fields and   
      ! therefore only applicable in a fully non-collinear description of magne-
      ! tism.
      ! 
      ! Assumes the following form for the potential matrix:
      ! V_mat = V*Id_(2x2) + sigma_vec*B_vec
      ! 
      ! B_vec is saved as a density type with no additional r^2-factor.

      TYPE(t_stars),                INTENT(IN)  :: stars
      TYPE(t_atoms),                INTENT(IN)  :: atoms
      TYPE(t_sphhar),               INTENT(IN)  :: sphhar
      TYPE(t_vacuum),               INTENT(IN)  :: vacuum
      TYPE(t_input),                INTENT(IN)  :: input
      TYPE(t_noco),                 INTENT(IN)  :: noco
      TYPE(t_potden),               INTENT(IN)  :: vTot
      TYPE(t_potden), DIMENSION(3), INTENT(OUT) :: bxc

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

      CALL matrixsplit(stars,atoms,sphhar,vacuum,input,noco,2.0,vtot, &
                       dummyDen(1),dummyDen(2),dummyDen(3),dummyDen(4))
      
      ! Initialize and fill the magnetic field.
      DO itype=1,atoms%ntype
         DO ir=1,atoms%jri(itype)
            r2=atoms%rmsh(ir,itype)**2
            DO i=1,3
               CALL bxc(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd, &
                           atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE., &
                           POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
               ALLOCATE(bxc(i)%pw_w,mold=bxc(i)%pw)

               bxc(i)%mt(:,0:,:,:)      = dummyDen(i+1)%mt(:,0:,:,:)/r2
               bxc(i)%pw(1:,:)          = dummyDen(i+1)%pw(1:,:)
               bxc(i)%vacz(1:,1:,:)     = dummyDen(i+1)%vacz(1:,1:,:)
               bxc(i)%vacxy(1:,1:,1:,:) = dummyDen(i+1)%vacxy(1:,1:,1:,:)
            ENDDO
         END DO
      END DO
     
   END SUBROUTINE makeBxc

   SUBROUTINE sourcefree(stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,bxc)

      ! Takes a vectorial quantity, i.e. a t_potden variable of dimension 3, and
      ! makes it into a source free vector field as follows:
      ! 
      ! a) Build the divergence d of the vector field A_vec as d=nabla_vec*A_vec.
      ! b) Solve the Poisson equation (nabla_vec*nabla_vec)phi=-4*pi*d.
      ! c) Construct an auxiliary vector field C_vec=(nabla phi)/(4*pi).
      ! d) Build A_vec_sf=A_vec+C_vec, which is source free by construction.

      TYPE(t_stars),                INTENT(IN)     :: stars
      TYPE(t_atoms),                INTENT(IN)     :: atoms
      TYPE(t_sphhar),               INTENT(IN)     :: sphhar
      TYPE(t_vacuum),               INTENT(IN)     :: vacuum
      TYPE(t_input),                INTENT(IN)     :: input
      TYPE(t_oneD),                 INTENT(IN)     :: oneD
      TYPE(t_sym),                  INTENT(IN)     :: sym
      TYPE(t_cell),                 INTENT(IN)     :: cell
      TYPE(t_noco),                 INTENT(IN)     :: noco
      TYPE(t_potden), DIMENSION(3), INTENT(IN)  :: bxc
      
      TYPE(t_potden)                               :: div!,phi
      !TYPE(t_potden), DIMENSION(3)                 :: cvec
      INTEGER                                      :: n

      CALL div%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype, &
                                  atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN, &
                                  vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      ALLOCATE(div%pw_w,mold=div%pw)
    
      DO n=1,atoms%ntype
         CALL divergence(n,stars,atoms,sphhar,vacuum,sym,cell,noco,bxc,div)
      END DO

      !CALL vDiv%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_POTCOUL,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      !ALLOCATE(vDiv%pw_w(SIZE(vDiv%pw,1),size(vDiv%pw,2)))
      !vDiv%pw_w = CMPLX(0.0,0.0)

      !CALL vgen_coulomb(1,mpi,dimension,oneD,input,field,vacuum,sym,stars,cell,sphhar,atoms,divB,vDiv)
    
      !DO i=1,3
      !   CALL graddiv(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      !   ALLOCATE(graddiv(i)%pw_w,mold=graddiv(i)%pw)
      !   CALL corrB(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
      !   ALLOCATE(corrB(i)%pw_w,mold=corrB(i)%pw)
      !ENDDO

      !DO i=1,atoms%ntype
      !   CALL divpotgrad(input%jspins,i,stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft,atoms,sphhar,sym,stars,cell,vacuum,noco,vDiv,graddiv)
      !END DO

      !DO i=1,3
      !   CALL corrB(i)%addPotDen(xcB(i), graddiv(i))
      !END DO

      CALL plotBtest(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                     noco, div)!, vDiv, graddiv(1), graddiv(2), graddiv(3), &
                     !corrB(1), corrB(2), corrB(3))

   END SUBROUTINE sourcefree

   !SUBROUTINE builddivtest()

   !   DO i=1,3
   !      CALL bxc(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd, &
   !                  atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE., &
   !                  POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
   !      ALLOCATE(bxc(i)%pw_w,mold=bxc(i)%pw)

   !      bxc(i)%mt(:,0:,:,:)      = dummyDen(i+1)%mt(:,0:,:,:)/r2
   !      bxc(i)%pw(1:,:)          = dummyDen(i+1)%pw(1:,:)
   !      bxc(i)%vacz(1:,1:,:)     = dummyDen(i+1)%vacz(1:,1:,:)
   !      bxc(i)%vacxy(1:,1:,1:,:) = dummyDen(i+1)%vacxy(1:,1:,1:,:)
   !   ENDDO

   !END SUBROUTINE builddivtest

END MODULE m_xcBfield
