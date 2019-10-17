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

   PUBLIC :: makeBxc, sourcefree, builddivtest

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

   SUBROUTINE builddivtest(stars,atoms,sphhar,vacuum,sym,cell,itest,Avec)
      USE m_mt_tofrom_grid
      USE m_pw_tofrom_grid

      IMPLICIT NONE

      TYPE(t_stars),                INTENT(IN)     :: stars
      TYPE(t_atoms),                INTENT(IN)     :: atoms
      TYPE(t_sphhar),               INTENT(IN)     :: sphhar
      TYPE(t_vacuum),               INTENT(IN)     :: vacuum
      TYPE(t_sym),                  INTENT(IN)     :: sym
      TYPE(t_cell),                 INTENT(IN)     :: cell
      INTEGER,                      INTENT(IN)     :: itest
      TYPE(t_potden), DIMENSION(3), INTENT(OUT)    :: Avec

      INTEGER                                      :: nsp, n, kt, ir, k, i, ifftxc3
      REAL                                         :: r, th, ph
      REAL, ALLOCATABLE                            :: thet(:), phi(:), A_temp(:,:)!space grid, index

      IF (itest.EQ.0) THEN 
         RETURN
      END IF

      nsp = atoms%nsp()
      ifftxc3=stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft

      DO i=1,3
         CALL Avec(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd, &
                       atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE., &
                       POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(Avec(i)%pw_w,mold=Avec(i)%pw)
      END DO

      ALLOCATE (thet(atoms%nsp()),phi(atoms%nsp()))

      CALL init_mt_grid(1, atoms, sphhar, .FALSE., sym, thet, phi)
      CALL init_pw_grid(.FALSE.,stars,sym,cell)
      !--------------------------------------------------------------------------
      ! Test case 1.
      ! In MT: radial function f(r)=r*R_MT in direction e_r.
      ! In interstitial: f(r_vec)=r_vec.
      ! So by construction, the divergence for a set of atoms with R_MT=rmt(n)
      ! for every n should be an isotropic 3*R_MT.
      ! TODO: Add interstitial part and if-test (1---> this test).

      DO n=1,atoms%ntype
         ALLOCATE (A_temp(atoms%jri(n)*nsp,3))
         kt = 0
         DO ir = 1, atoms%jri(n)
            r = atoms%rmsh(ir, n)
            DO k = 1, nsp
               th = thet(k)
               ph = phi(k)
               A_temp(kt+nsp,1)=COS(ph)*SIN(th)*atoms%rmt(n)*atoms%rmsh(ir,n)
               A_temp(kt+nsp,2)=SIN(ph)*SIN(th)*atoms%rmt(n)*atoms%rmsh(ir,n)
               A_temp(kt+nsp,3)=        COS(th)*atoms%rmt(n)*atoms%rmsh(ir,n)
            ENDDO ! k
            kt = kt + nsp
         END DO ! ir

         CALL mt_from_grid(atoms, sphhar, n, 1, A_temp(:,:), Avec(1)%mt(:,0:,n,:))
         CALL mt_from_grid(atoms, sphhar, n, 1, A_temp(:,:), Avec(2)%mt(:,0:,n,:))
         CALL mt_from_grid(atoms, sphhar, n, 1, A_temp(:,:), Avec(3)%mt(:,0:,n,:))
         DEALLOCATE (A_temp)
      END DO ! n

      ALLOCATE (A_temp(ifftxc3,3))

      !DO i = 1, ifftxc3
      !   div_temp(i,1)=
      !   div_temp(i,2)=
      !   div_temp(i,3)=
      !ENDDO ! i

      ! End test case 1.
      !--------------------------------------------------------------------------
      CALL finish_mt_grid
      CALL finish_pw_grid()

   END SUBROUTINE builddivtest

END MODULE m_xcBfield
