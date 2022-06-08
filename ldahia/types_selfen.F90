MODULE m_types_selfen

   !------------------------------------------------------------------------
   !This type contains the array for the selfenergy from the impurity solver
   !We separate the individual elements because the number of energy points
   !can differ massively and would lead to wasted storage
   !------------------------------------------------------------------------

   USE m_constants

   IMPLICIT NONE

   PRIVATE

   TYPE t_selfen

      INTEGER :: l = -1
      REAL,ALLOCATABLE    :: muMatch(:)

      COMPLEX, ALLOCATABLE :: data(:,:,:,:)

      CONTAINS
         PROCEDURE, PASS :: init    => init_selfen
         PROCEDURE       :: collect => collect_selfen
         PROCEDURE       :: postProcess => postProcess_selfen

   END TYPE t_selfen

   PUBLIC t_selfen

   CONTAINS

      SUBROUTINE init_selfen(this,l,nz,jspins,l_fullMatch)

         CLASS(t_selfen), INTENT(INOUT) :: this
         INTEGER,         INTENT(IN)    :: l
         INTEGER,         INTENT(IN)    :: nz
         INTEGER,         INTENT(IN)    :: jspins
         LOGICAL,         INTENT(IN)    :: l_fullMatch

         this%l = l
         ALLOCATE(this%muMatch(MERGE(1,jspins,l_fullMatch)),source=0.0)
         ALLOCATE(this%data(2*(2*l+1),2*(2*l+1),nz,2),source = cmplx_0)

      END SUBROUTINE init_selfen

      SUBROUTINE collect_selfen(this,mpi_communicator)

#ifdef CPP_MPI
         USE mpi
#endif

         CLASS(t_selfen),     INTENT(INOUT) :: this
         INTEGER,             INTENT(IN)    :: mpi_communicator
#ifdef CPP_MPI

         INTEGER:: ierr,irank,n
         COMPLEX,ALLOCATABLE::ctmp(:)
         REAL, ALLOCATABLE :: rtmp(:)

         CALL MPI_COMM_RANK(mpi_communicator,irank,ierr)

         n = SIZE(this%muMatch)
         ALLOCATE(rtmp(n))
         CALL MPI_REDUCE(this%muMatch,rtmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_communicator,ierr)
         IF(irank.EQ.0) this%muMatch = reshape(rtmp,[n])
         DEALLOCATE(rtmp)

         n = SIZE(this%data)
         ALLOCATE(ctmp(n))
         CALL MPI_REDUCE(this%data,ctmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_communicator,ierr)
         IF(irank.EQ.0) CALL zcopy(n,ctmp,1,this%data,1)
         DEALLOCATE(ctmp)
#endif

      END SUBROUTINE collect_selfen

      SUBROUTINE postProcess_selfen(this,noco,nococonv,atomType,l,jspins,vmmp)

         USE m_types_noco
         USE m_types_nococonv
         USE m_rotMMPmat

         CLASS(t_selfen), INTENT(INOUT) :: this
         TYPE(t_noco),    INTENT(IN)    :: noco
         TYPE(t_nococonv),INTENT(IN)    :: nococonv
         INTEGER,         INTENT(IN)    :: atomType,l
         INTEGER,         INTENT(IN)    :: jspins
         COMPLEX,         INTENT(IN)    :: vmmp(-lmaxU_const:,-lmaxU_const:,:)

         INTEGER :: i,j,iz,ipm,m,mp,ispin,ns
         COMPLEX,ALLOCATABLE :: swapMat(:,:)
         COMPLEX,ALLOCATABLE :: vmmp_local(:,:,:)

         ns = 2*this%l+1

         ALLOCATE(swapMat(2*ns,2*ns),source=cmplx_0)
         ALLOCATE(vmmp_local(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(vmmp,3)),source=cmplx_0)

         !Transformation matrix is a Block matrix of form
         ! | 0  I |
         ! | I  0 |
         !to swap the order of the spins
         swapMat = 0.0
         DO i = 1, ns
            swapMat(i,ns+i) = 1.0
            swapMat(ns+i,i) = 1.0
         ENDDO

         !The DFT+U correction is in the global frame of real space
         !For the calculation of the impurity greens function we shift into the local frame
         IF(noco%l_noco) THEN
            vmmp_local = rotMMPmat(vmmp,0.0,-nococonv%beta(atomType),-nococonv%alph(atomType),l)
         ELSE IF(noco%l_soc) THEN
            vmmp_local = rotMMPmat(vmmp,0.0,-nococonv%theta,-nococonv%phi,l)
         ELSE
            vmmp_local = vmmp
         ENDIF

         DO ipm = 1, 2
            DO iz = 1, SIZE(this%data,3)
               !---------------------------------------------
               ! Convert the selfenergy to hartree
               !---------------------------------------------
               this%data(:,:,iz,ipm) = this%data(:,:,iz,ipm)/hartree_to_ev_const
               !---------------------------------------------
               ! The order of spins is reversed in the Solver (transformation matrix is symmetric)
               !---------------------------------------------
               this%data(:,:,iz,ipm) = matmul(this%data(:,:,iz,ipm),swapMat)
               this%data(:,:,iz,ipm) = matmul(swapMat,this%data(:,:,iz,ipm))
               !---------------------------------------------------------------------
               ! The DFT green's function also includes the previous DFT+U correction
               ! This is removed by substracting it from the selfenergy
               !---------------------------------------------------------------------
               DO i = 1, ns
                  m  = i-1-this%l
                  DO j = 1, ns
                     mp = j-1-this%l
                     DO ispin = 1, SIZE(vmmp_local,3)
                        IF(ispin < 3) THEN
                           this%data(i+(ispin-1)*ns,j+(ispin-1)*ns,iz,ipm) = this%data(i+(ispin-1)*ns,j+(ispin-1)*ns,iz,ipm) &
                                                                             - vmmp_local(m,mp,ispin)/(3.0-jspins)
                           IF(jspins.EQ.1) this%data(i+ns,j+ns,iz,ipm) = this%data(i+ns,j+ns,iz,ipm) - vmmp_local(-m,-mp,ispin)/(3.0-jspins)
                        ELSE
                           !----------------------------------------------------------------------------
                           ! The offdiagonal elements only have to be removed if they are actually added
                           ! to the hamiltonian (so noco%l_mperp and any(noco%l_unrestrictMT))
                           !----------------------------------------------------------------------------
                           this%data(i+ns,j,iz,ipm) = this%data(i+ns,j,iz,ipm) - vmmp_local(m,mp,ispin)
                           this%data(i,j+ns,iz,ipm) = this%data(i,j+ns,iz,ipm) - conjg(vmmp_local(mp,m,ispin))
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      END SUBROUTINE postProcess_selfen

END MODULE m_types_selfen
