!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_usetup
   !-------------------------------------------------------------------+
   ! Sets up the quantities needed for the LDA+U subroutines:          |                                         |
   !     potential matrix: vs_mmp                                      |
   !     total energy contribution: e_ldau                             |
   !                                                  G.B. Oct. 2000   |
   !                                                                   |
   !     Extension to multiple U per atom type  G.M. 2017              |
   !-------------------------------------------------------------------+
   USE m_juDFT
   USE m_umtx
   USE m_uj2f
   USE m_nmat_rot
   USE m_vmmp
   USE m_vmmp21
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE u_setup(atoms,input,noco,mpi,hub1inp,inDen,pot,results)

      TYPE(t_mpi),    INTENT(IN)    :: mpi
      TYPE(t_input),  INTENT(IN)    :: input
      TYPE(t_atoms),  INTENT(IN)    :: atoms
      TYPE(t_noco),   INTENT(IN)    :: noco
      TYPE(t_hub1inp),INTENT(IN)    :: hub1inp
      TYPE(t_potden), INTENT(IN)    :: inDen
      TYPE(t_potden), INTENT(INOUT) :: pot
      TYPE(t_results),INTENT(INOUT) :: results


      INTEGER :: itype,ispin,j,k,l,i_u,n_u
      REAL    :: e_off
      CHARACTER(len=2) :: l_type
      CHARACTER(len=9) :: l_form

      REAL :: zero(atoms%n_u+atoms%n_hia)
      REAL :: f0(atoms%n_u+atoms%n_hia),f2(atoms%n_u+atoms%n_hia)
      REAL :: f4(atoms%n_u+atoms%n_hia),f6(atoms%n_u+atoms%n_hia)
      REAL,    ALLOCATABLE :: u(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: n_mmp(:,:,:,:)

      ! look, whether density matrix exists already:
      IF (ANY(ABS(inDen%mmpMat).GT.1e-12)) THEN

         n_u = atoms%n_u+atoms%n_hia

         ! calculate slater integrals from u and j
         CALL uj2f(input%jspins,atoms%lda_u(:),n_u,f0,f2,f4,f6)

         ! set up e-e- interaction matrix
         ALLOCATE ( u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
                      -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,n_u)), source=0.0 )
         CALL umtx(atoms%lda_u(:),n_u,f0,f2,f4,f6,u)

         ALLOCATE ( n_mmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,n_u),&
                          input%jspins), source=cmplx_0 )
         n_mmp = inDen%mmpMat(:,:,:,1:input%jspins)

         ! check for possible rotation of n_mmp
         zero = 0.0
         CALL nmat_rot(atoms%lda_u(:)%phi,atoms%lda_u(:)%theta,zero,3,n_u,input%jspins,atoms%lda_u(:)%l,n_mmp)

         ! calculate potential matrix and total energy correction
         CALL v_mmp(atoms,input%jspins,hub1inp%l_dftspinpol,n_mmp,u,f0,f2,pot%mmpMat,results%e_ldau)

         !spin off-diagonal elements (no rotation yet)
         IF(noco%l_mtNocoPot) THEN
            IF(ANY(atoms%lda_u(:)%phi.NE.0.0).OR.ANY(atoms%lda_u(:)%theta.NE.0.0)) CALL juDFT_error("vmmp21+Rot not implemented", calledby="u_setup")
            CALL v_mmp_21(atoms,inDen%mmpMat(:,:,:,3),u,pot%mmpMat(:,:,:,3),e_off)
            results%e_ldau = results%e_ldau + e_off
         ELSE IF(noco%l_mperp.AND.mpi%irank.EQ.0) THEN
            WRITE(*,*) "Offdiagonal LDA+U ignored"
            WRITE(oUnit,*) "The offdiagonal contributions to LDA+U only enters into the Hamiltonian if l_mtNocoPot is also true"
         ENDIF

         IF (mpi%irank.EQ.0) THEN
            DO ispin = 1,SIZE(pot%mmpMat,4)
               WRITE (oUnit,'(a7,i3)') 'spin #',ispin
               DO i_u = 1, n_u
                  itype = atoms%lda_u(i_u)%atomType
                  l     = atoms%lda_u(i_u)%l

                  WRITE (l_type,'(i2)') 2*(2*l+1)
                  l_form = '('//l_type//'f12.7)'
                  WRITE (oUnit,'(a20,i3)') 'n-matrix for atom # ',itype
                  IF(i_u > atoms%n_u) WRITE(oUnit,"(A)") 'n-matrix calculated with DFT+Hubbard-1'
                  IF(ispin < 3) THEN
                     WRITE (oUnit,l_form) ((n_mmp(k,j,i_u,ispin),k=-l,l),j=-l,l)
                  ELSE
                     WRITE (oUnit,l_form) ((inDen%mmpMat(k,j,i_u,ispin),k=-l,l),j=-l,l)
                  ENDIF
                  WRITE (oUnit,'(a20,i3)') 'V-matrix for atom # ',itype
                  IF (atoms%lda_u(i_u)%l_amf) THEN
                     WRITE (oUnit,*) 'using the around-mean-field limit '
                  ELSE
                     WRITE (oUnit,*) 'using the atomic limit of LDA+U '
                  ENDIF
                  WRITE (oUnit,l_form) ((pot%mmpMat(k,j,i_u,ispin),k=-l,l),j=-l,l)
               END DO
            END DO
            WRITE (oUnit,*) results%e_ldau
            IF(noco%l_mtNocoPot) WRITE(oUnit,*) e_off
         ENDIF
         DEALLOCATE (u,n_mmp)
      ELSE
         IF (mpi%irank.EQ.0) THEN
            WRITE (*,*) 'no density matrix found ... skipping LDA+U'
         ENDIF
         pot%mmpMat = cmplx_0
         results%e_ldau = 0.0
      ENDIF

   END SUBROUTINE u_setup
END MODULE m_usetup
