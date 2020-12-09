!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_greensf

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_types_greensf
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION:
   !>  Contains a type for onsite and intersite green's functions in the mt-sphere
   !>  It stores the energy contour in the complex plane and the corresponding
   !>  matrix elements of the green's function
   !------------------------------------------------------------------------------

   USE m_juDFT
   USE m_constants
   USE m_types_setup
   USE m_types_greensfContourData

   IMPLICIT NONE

   PRIVATE

   TYPE t_greensf

      LOGICAL :: l_calc = .FALSE.

      !Energy contour parameters
      TYPE(t_greensfContourData) :: contour

      !Pointer to the element type in gfinp
      TYPE(t_gfelementtype), POINTER :: elem => NULL()

      !Arrays for Green's function
      COMPLEX, ALLOCATABLE :: gmmpMat(:,:,:,:,:)

      !for radial dependence
      COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: du(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:)

      COMPLEX, ALLOCATABLE :: uulo(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: ulou(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: dulo(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: ulod(:,:,:,:,:,:)

      COMPLEX, ALLOCATABLE :: uloulop(:,:,:,:,:,:,:)

      CONTAINS
         PROCEDURE, PASS :: init                => init_greensf
         PROCEDURE       :: mpi_bc              => mpi_bc_greensf
         PROCEDURE       :: collect             => collect_greensf
         PROCEDURE       :: get                 => get_gf
         PROCEDURE       :: getRadial           => getRadial_gf
         PROCEDURE       :: getRadialSpin       => getRadialSpin_gf
         PROCEDURE       :: getRadialRadial     => getRadialRadial_gf!(Full Radial dependence for intersite)
         PROCEDURE       :: getRadialRadialSpin => getRadialRadialSpin_gf
         PROCEDURE       :: integrateOverMT     => integrateOverMT_greensf
         PROCEDURE       :: set                 => set_gf
         PROCEDURE       :: rotate              => rotate_gf
         PROCEDURE       :: reset               => reset_gf
         PROCEDURE       :: resetSingleElem     => resetSingleElem_gf
         PROCEDURE       :: checkEmpty          => checkEmpty_greensf
   END TYPE t_greensf

   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE init_greensf(this,gfelem,gfinp,atoms,input,contour_in,l_sphavg_in)

         CLASS(t_greensf),                     INTENT(INOUT)  :: this
         TYPE(t_gfelementtype),      TARGET,   INTENT(IN)     :: gfelem
         TYPE(t_gfinp),                        INTENT(IN)     :: gfinp
         TYPE(t_atoms),                        INTENT(IN)     :: atoms
         TYPE(t_input),                        INTENT(IN)     :: input
         !Pass a already calculated energy contour to the type
         TYPE(t_greensfContourData), OPTIONAL, INTENT(IN)     :: contour_in
         LOGICAL,                    OPTIONAL, INTENT(IN)     :: l_sphavg_in !To overwrite the allocation for integrateOverMT_greensf

         INTEGER spin_dim,lmax,nLO
         LOGICAL l_sphavg

         this%elem => gfelem
         this%l_calc = .FALSE.

         nLO = this%elem%countLOs(atoms)

         !Initialize the contour
         CALL this%contour%init(gfinp%contour(this%elem%iContour),contour_in=contour_in)

         spin_dim = MERGE(3,input%jspins,gfinp%l_mperp)
         lmax = lmaxU_const

         l_sphavg = this%elem%l_sphavg
         IF(PRESENT(l_sphavg_in)) l_sphavg = l_sphavg_in

         IF(l_sphavg) THEN
            ALLOCATE(this%gmmpMat(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
         ELSE
            ALLOCATE(this%uu(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
            ALLOCATE(this%dd(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
            ALLOCATE(this%du(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
            ALLOCATE(this%ud(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)

            IF(nLO>0) THEN
               ALLOCATE(this%uulo(this%contour%nz,-lmax:lmax,-lmax:lmax,nLO,spin_dim,2),source=cmplx_0)
               ALLOCATE(this%ulou(this%contour%nz,-lmax:lmax,-lmax:lmax,nLO,spin_dim,2),source=cmplx_0)
               ALLOCATE(this%dulo(this%contour%nz,-lmax:lmax,-lmax:lmax,nLO,spin_dim,2),source=cmplx_0)
               ALLOCATE(this%ulod(this%contour%nz,-lmax:lmax,-lmax:lmax,nLO,spin_dim,2),source=cmplx_0)

               ALLOCATE(this%uloulop(this%contour%nz,-lmax:lmax,-lmax:lmax,nLO,nLO,spin_dim,2),source=cmplx_0)
            ENDIF
         ENDIF

      END SUBROUTINE init_greensf

      SUBROUTINE mpi_bc_greensf(this,mpi_comm,irank)
         USE m_mpi_bc_tool
         CLASS(t_greensf), INTENT(INOUT)::this
         INTEGER, INTENT(IN):: mpi_comm
         INTEGER, INTENT(IN), OPTIONAL::irank
         INTEGER ::rank
         IF (PRESENT(irank)) THEN
            rank = irank
         ELSE
            rank = 0
         END IF

         CALL mpi_bc(this%l_calc,rank,mpi_comm)

         CALL this%contour%mpi_bc(mpi_comm,rank)

         IF(ALLOCATED(this%gmmpMat)) CALL mpi_bc(this%gmmpMat,rank,mpi_comm)
         IF(ALLOCATED(this%uu)) CALL mpi_bc(this%uu,rank,mpi_comm)
         IF(ALLOCATED(this%ud)) CALL mpi_bc(this%ud,rank,mpi_comm)
         IF(ALLOCATED(this%du)) CALL mpi_bc(this%du,rank,mpi_comm)
         IF(ALLOCATED(this%dd)) CALL mpi_bc(this%dd,rank,mpi_comm)
         IF(ALLOCATED(this%uulo)) CALL mpi_bc(this%uulo,rank,mpi_comm)
         IF(ALLOCATED(this%ulou)) CALL mpi_bc(this%ulou,rank,mpi_comm)
         IF(ALLOCATED(this%dulo)) CALL mpi_bc(this%dulo,rank,mpi_comm)
         IF(ALLOCATED(this%ulod)) CALL mpi_bc(this%ulod,rank,mpi_comm)
         IF(ALLOCATED(this%uloulop)) CALL mpi_bc(this%uloulop,rank,mpi_comm)

      END SUBROUTINE mpi_bc_greensf

      SUBROUTINE collect_greensf(this,mpi_communicator)

#ifdef CPP_MPI
         USE mpi
#endif

         CLASS(t_greensf),     INTENT(INOUT) :: this
         INTEGER,              INTENT(IN)    :: mpi_communicator
#ifdef CPP_MPI
#include"cpp_double.h"
         INTEGER:: ierr,n
         COMPLEX,ALLOCATABLE::ctmp(:)

         IF(ALLOCATED(this%gmmpMat)) THEN
            n = SIZE(this%gmmpMat)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%gmmpMat,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%gmmpMat,1)
            DEALLOCATE(ctmp)
         ELSE
            n = SIZE(this%uu)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%uu,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%uu,1)
            CALL MPI_ALLREDUCE(this%ud,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%ud,1)
            CALL MPI_ALLREDUCE(this%du,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%du,1)
            CALL MPI_ALLREDUCE(this%dd,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%dd,1)
            DEALLOCATE(ctmp)

            IF(ALLOCATED(this%uulo)) THEN
               n = SIZE(this%uulo)
               ALLOCATE(ctmp(n))
               CALL MPI_ALLREDUCE(this%uulo,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
               CALL CPP_BLAS_ccopy(n,ctmp,1,this%uulo,1)
               CALL MPI_ALLREDUCE(this%ulou,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
               CALL CPP_BLAS_ccopy(n,ctmp,1,this%ulou,1)
               CALL MPI_ALLREDUCE(this%dulo,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
               CALL CPP_BLAS_ccopy(n,ctmp,1,this%dulo,1)
               CALL MPI_ALLREDUCE(this%ulod,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
               CALL CPP_BLAS_ccopy(n,ctmp,1,this%ulod,1)
               DEALLOCATE(ctmp)

               n = SIZE(this%uloulop)
               ALLOCATE(ctmp(n))
               CALL MPI_ALLREDUCE(this%uloulop,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
               CALL CPP_BLAS_ccopy(n,ctmp,1,this%uloulop,1)
               DEALLOCATE(ctmp)
            ENDIF
         ENDIF
#endif
      END SUBROUTINE collect_greensf

      !----------------------------------------------------------------------------------
      ! Following this comment there are multiple definitions for functions
      ! to access the data in the greensFunction Type:
      !     get_gf -> Get the (m,mp) matrix of the spherically averaged GF
      !               at a certain energy point. If the correct scalar products
      !               are provided, the radial dependent GF can also be recombined here
      !     getRadial_gf -> Returns the radial and energy dependent GF for a certain spin
      !                     and m,mp pair
      !     getRadialSpin_gf -> Returns the radial and energy dependent GF for a certain
      !                         m,mp pair. Also returns the 2x2 spin matrix at that point
      !     set_gf -> Set the value of the (m,mp) Matrix at a
      !               certain energy point with an input matrix
      !----------------------------------------------------------------------------------

      SUBROUTINE get_gf(this,atoms,iz,l_conjg,gmat,spin,usdus,denCoeffsOffDiag,scalarGF)

         USE m_types_mat
         USE m_types_usdus
         USE m_types_denCoeffsOffDiag
         USE m_types_scalarGF

         !Returns the matrix belonging to energy point iz with l,lp,nType,nTypep
         !can also return the spherically averaged GF with the given scalar products

         CLASS(t_greensf),                   INTENT(IN)     :: this
         TYPE(t_atoms),                      INTENT(IN)     :: atoms
         INTEGER,                            INTENT(IN)     :: iz
         LOGICAL,                            INTENT(IN)     :: l_conjg
         TYPE(t_mat),                        INTENT(INOUT)  :: gmat !Return matrix
         INTEGER,                  OPTIONAL, INTENT(IN)     :: spin
         TYPE(t_usdus),            OPTIONAL, INTENT(IN)     :: usdus
         TYPE(t_denCoeffsOffDiag), OPTIONAL, INTENT(IN)     :: denCoeffsOffDiag
         TYPE(t_scalarGF),         OPTIONAL, INTENT(IN)     :: scalarGF

         INTEGER matsize1,matsize2,ind1,ind2,ind1_start,ind2_start
         INTEGER m,mp,spin1,spin2,ipm,ispin,spin_start,spin_end,spin_ind,m_ind,mp_ind
         INTEGER l,lp,atomType,atomTypep,nspins,ilo,ilop,iLO_ind,iLOp_ind
         LOGICAL l_full,l_scalar,l_scalarGF

         REAL :: uun,dun,udn,ddn
         REAL :: uulon(atoms%nlod),dulon(atoms%nlod),ulodn(atoms%nlod),uloun(atoms%nlod)
         REAL :: uloulopn(atoms%nlod,atoms%nlod)

         IF(.NOT.this%l_calc) THEN
            CALL juDFT_error("The requested Green's Function element was not calculated", calledby="get_gf")
         ENDIF

         l  = this%elem%l
         lp = this%elem%lp
         atomType  = this%elem%atomType
         atomTypep = this%elem%atomTypep

         IF(ALLOCATED(this%gmmpMat)) THEN
            nspins = SIZE(this%gmmpMat,4)
         ELSE
            nspins = SIZE(this%uu,4)
         ENDIF

         l_scalarGF = PRESENT(scalarGF)
         IF(l_scalarGF) l_scalarGF = scalarGF%done.AND.this%elem%isOffDiag()
         l_scalar = PRESENT(usdus).OR.PRESENT(denCoeffsOffDiag)
         IF(l_scalar.AND.nspins==3) THEN
            IF(.NOT.PRESENT(denCoeffsOffDiag)) THEN
                  CALL juDFT_error("Offdiagonal Scalar products missing", calledby="get_gf")
            ENDIF
         ENDIF

         IF(l_scalar.AND. .NOT.ALLOCATED(this%uu)) THEN
            CALL juDFT_error("l_scalar/l_radial only without l_sphavg", calledby="get_gf")
         ENDIF

         IF(PRESENT(spin)) THEN
            IF(spin.GT.4 .OR. spin.LT.1) THEN
               CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
            ENDIF
         END IF

         !Determine matsize for the result gmat (if spin is given only return this diagonal element)
         l_full = .NOT.PRESENT(spin)
         matsize1 = (2*l+1) * MERGE(2,1,l_full)
         matsize2 = (2*lp+1) * MERGE(2,1,l_full)

         IF(.NOT.ALLOCATED(gmat%data_c)) THEN
            CALL gmat%init(.FALSE.,matsize1,matsize2)
         ELSE IF(matsize1.NE.gmat%matsize1.OR.matsize2.NE.gmat%matsize2) THEN
            CALL juDFT_error("Mismatch in matsizes", calledby="get_gf")
         ENDIF

         ipm = MERGE(2,1,l_conjg)

         gmat%data_c = cmplx_0

         IF(l_full) THEN
            spin_start = 1
            spin_end   = MERGE(4,2,nspins.EQ.3)
         ELSE
            spin_start = spin
            spin_end   = spin
         ENDIF

         DO ispin = spin_start, spin_end
            !Find the corresponding physical spin indices
            IF(ispin < 3) THEN
               spin1 = ispin
               spin2 = ispin
            ELSE IF(ispin.EQ.3) THEN
               spin1 = 2
               spin2 = 1
            ELSE
               spin1 = 1
               spin2 = 2
            ENDIF
            !Find the correct spin index in gmmpMat arrays
            spin_ind = MERGE(1,ispin,nspins.EQ.1)
            spin_ind = MERGE(3,spin_ind,ispin.EQ.4)
            !Find the right quadrant in gmat
            IF(l_full) THEN
               ind1_start = (spin2-1)*(2*l+1)
               ind2_start = (spin1-1)*(2*lp+1)
            ELSE
               ind1_start = 0
               ind2_start = 0
            ENDIF

            IF(l_scalar.OR.l_scalarGF) THEN
               !Select the correct scalar products or integrals (So we do not have to repeat the actual combination)
               IF(l_scalarGF) THEN
                  uun = scalarGF%uun(spin1,spin2)
                  dun = scalarGF%dun(spin1,spin2)
                  udn = scalarGF%udn(spin1,spin2)
                  ddn = scalarGF%ddn(spin1,spin2)

                  IF(ALLOCATED(this%uulo)) THEN
                     uulon(:) = scalarGF%uulon(:,spin1,spin2)
                     uloun(:) = scalarGF%uloun(:,spin1,spin2)
                     dulon(:) = scalarGF%dulon(:,spin1,spin2)
                     ulodn(:) = scalarGF%ulodn(:,spin1,spin2)

                     uloulopn(:,:) = scalarGF%uloulopn(:,:,spin1,spin2)
                  ENDIF
               ELSE IF(spin_ind<3) THEN
                  uun = 1.0
                  dun = 0.0
                  udn = 0.0
                  ddn = usdus%ddn(l,atomType,spin_ind)

                  IF(ALLOCATED(this%uulo)) THEN
                     uulon(:) = usdus%uulon(:,atomType,spin_ind)
                     uloun(:) = usdus%uulon(:,atomType,spin_ind)
                     dulon(:) = usdus%dulon(:,atomType,spin_ind)
                     ulodn(:) = usdus%dulon(:,atomType,spin_ind)

                     uloulopn(:,:) = usdus%uloulopn(:,:,atomType,spin_ind)
                  ENDIF
               ELSE
                  uun = denCoeffsOffDiag%uu21n(l,atomType)
                  dun = denCoeffsOffDiag%du21n(l,atomType)
                  udn = denCoeffsOffDiag%ud21n(l,atomType)
                  ddn = denCoeffsOffDiag%dd21n(l,atomType)

                  IF(ALLOCATED(this%uulo)) THEN
                     uulon(:) = denCoeffsOffDiag%uulo21n(:,atomType)
                     uloun(:) = denCoeffsOffDiag%ulou21n(:,atomType)
                     dulon(:) = denCoeffsOffDiag%dulo21n(:,atomType)
                     ulodn(:) = denCoeffsOffDiag%ulod21n(:,atomType)

                     uloulopn(:,:) = denCoeffsOffDiag%uloulop21n(:,:,atomType)
                  ENDIF
               ENDIF

            ENDIF

            ind1 = ind1_start
            DO m = -l,l
               ind1 = ind1 + 1
               ind2 = ind2_start
               DO mp = -lp,lp
                  ind2 = ind2 + 1

                  !-------------------------------------------------------------------
                  ! Check wether we need to do some operation on the indices m and mp
                  !-------------------------------------------------------------------
                  IF(ispin.EQ.2 .AND. nspins.EQ.1) THEN
                     !For a non-spin-polarized calculation we might still want the full
                     !matrix. Then we need to reverse the order (SOC prop m*s_z)
                     m_ind  = -m
                     mp_ind = -mp
                  ELSE IF(ispin.EQ.4) THEN
                     !We only calculate spin21. spin12 is obtained as hermitian conjugate
                     !(Complex conjugation happens afterwards)
                     m_ind  = mp
                     mp_ind = m
                  ELSE
                     !Do nothing
                     m_ind  = m
                     mp_ind = mp
                  ENDIF
                  !-------------------
                  ! Fetch the values
                  !-------------------
                  IF(l_scalar.OR.l_scalarGF) THEN
                     gmat%data_c(ind1,ind2) =   this%uu(iz,m_ind,mp_ind,spin_ind,ipm) * uun &
                                              + this%dd(iz,m_ind,mp_ind,spin_ind,ipm) * ddn &
                                              + this%du(iz,m_ind,mp_ind,spin_ind,ipm) * dun &
                                              + this%ud(iz,m_ind,mp_ind,spin_ind,ipm) * udn
                     IF(ALLOCATED(this%uulo)) THEN
                        iLO_ind = 0
                        DO ilo = 1, atoms%nlo(atomType)
                           IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                           iLO_ind = iLO_ind + 1
                           gmat%data_c(ind1,ind2) =   gmat%data_c(ind1,ind2) &
                                                    + this%uulo(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * uulon(ilo) &
                                                    + this%dulo(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * dulon(ilo)
                        ENDDO
                        iLO_ind = 0
                        DO ilo = 1, atoms%nlo(atomTypep)
                           IF(atoms%llo(ilo,atomTypep).NE.lp) CYCLE
                           iLO_ind = iLO_ind + 1
                           gmat%data_c(ind1,ind2) =   gmat%data_c(ind1,ind2) &
                                                    + this%ulou(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * uloun(ilo) &
                                                    + this%dulo(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ulodn(ilo)
                        ENDDO
                        iLO_ind = 0
                        DO ilo = 1,atoms%nlo(atomType)
                           IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                           iLO_ind = iLO_ind + 1
                           iLOp_ind = 0
                           DO ilop = 1, atoms%nlo(atomTypep)
                              IF(atoms%llo(ilop,atomTypep).NE.lp) CYCLE
                              iLOp_ind = iLOp_ind + 1
                              gmat%data_c(ind1,ind2) =   gmat%data_c(ind1,ind2) &
                                                       + this%uloulop(iz,m_ind,mp_ind,iLO_ind,iLOp_ind,spin_ind,ipm) &
                                                        * uloulopn(ilo,ilop)
                           ENDDO
                        ENDDO
                     ENDIF
                  ELSE
                     gmat%data_c(ind1,ind2) = this%gmmpMat(iz,m_ind,mp_ind,spin_ind,ipm)
                  ENDIF
                  !------------------------
                  ! Additional operations
                  !------------------------
                  !Spin-degeneracy when using a full matrix and having input%jspins.EQ.1
                  IF(l_full) gmat%data_c(ind1,ind2) = gmat%data_c(ind1,ind2) * MERGE(0.5,1.0,nspins.EQ.1)
                  !Complex conjugate for spin 4
                  IF(ispin.EQ.4) gmat%data_c(ind1,ind2) = conjg(gmat%data_c(ind1,ind2))

               ENDDO!mp
            ENDDO!m
         ENDDO!ispin

      END SUBROUTINE get_gf

      SUBROUTINE getRadial_gf(this,atoms,m,mp,l_conjg,spin,f,g,flo,gmat)

         !Returns the green's function on the radial and energy mesh
         !for a certain m,mp,spin combination. Attention: The correct radial functions have to be provided

         CLASS(t_greensf),    INTENT(IN)     :: this
         TYPE(t_atoms),       INTENT(IN)     :: atoms
         INTEGER,             INTENT(IN)     :: m,mp
         LOGICAL,             INTENT(IN)     :: l_conjg
         INTEGER,             INTENT(IN)     :: spin
         REAL   ,             INTENT(IN)     :: f(:,:,0:,:,:)
         REAL   ,             INTENT(IN)     :: g(:,:,0:,:,:)
         REAL   ,             INTENT(IN)     :: flo(:,:,:,:,:)
         COMPLEX, ALLOCATABLE,INTENT(INOUT)  :: gmat(:,:) !Return matrix

         INTEGER spin1,spin2,ipm,spin_ind,m_ind,mp_ind,ilo,ilop,iLO_ind,iLOp_ind
         INTEGER l,lp,atomType,atomTypep,nspins,iz

         IF(.NOT.this%l_calc) THEN
            CALL juDFT_error("The requested Green's Function element was not calculated", calledby="get_gf")
         ENDIF

         l  = this%elem%l
         lp = this%elem%lp
         atomType  = this%elem%atomType
         atomTypep = this%elem%atomTypep

         IF(ALLOCATED(this%gmmpMat)) THEN
            CALL juDFT_error("Green's function not calculated for radial dependence", calledby="get_gf")
         ENDIF

         nspins = SIZE(this%uu,4)

         IF(spin.GT.4 .OR. spin.LT.1) THEN
            CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
         ENDIF

         ipm = MERGE(2,1,l_conjg)

         IF(.NOT.ALLOCATED(gmat)) ALLOCATE(gmat(SIZE(f,1),this%contour%nz),source=cmplx_0)
         gmat = cmplx_0

         IF(spin < 3) THEN
            spin1 = spin
            spin2 = spin
         ELSE IF(spin.EQ.3) THEN
            spin1 = 2
            spin2 = 1
         ELSE
            spin1 = 1
            spin2 = 2
         ENDIF
         !Find the correct spin index in gmmpMat arrays
         spin_ind = MERGE(1,spin,nspins.EQ.1)
         spin_ind = MERGE(3,spin_ind,spin.EQ.4)

         !-------------------------------------------------------------------
         ! Check wether we need to do some operation on the indices m and mp
         !-------------------------------------------------------------------
         IF(spin.EQ.2 .AND. nspins.EQ.1) THEN
            !For a non-spin-polarized calculation we might still want the full
            !matrix. Then we need to reverse the order (SOC prop m*s_z)
            m_ind  = -m
            mp_ind = -mp
         ELSE IF(spin.EQ.4) THEN
            !We only calculate spin21. spin12 is obtained as hermitian conjugate
            !(Complex conjugation happens afterwards)
            m_ind  = mp
            mp_ind = m
         ELSE
            !Do nothing
            m_ind  = m
            mp_ind = mp
         ENDIF
         !-------------------
         ! Fetch the values
         !-------------------
         DO iz = 1, this%contour%nz
            gmat(:,iz) =   this%uu(iz,m_ind,mp_ind,spin_ind,ipm) * ( f(:,1,l,spin2,atomType) * f(:,1,lp,spin1,atomTypep) &
                                                                    +f(:,2,l,spin2,atomType) * f(:,2,lp,spin1,atomTypep))&
                         + this%dd(iz,m_ind,mp_ind,spin_ind,ipm) * ( g(:,1,l,spin2,atomType) * g(:,1,lp,spin1,atomTypep) &
                                                                    +g(:,2,l,spin2,atomType) * g(:,2,lp,spin1,atomTypep))&
                         + this%ud(iz,m_ind,mp_ind,spin_ind,ipm) * ( g(:,1,l,spin2,atomType) * f(:,1,lp,spin1,atomTypep) &
                                                                    +g(:,2,l,spin2,atomType) * f(:,2,lp,spin1,atomTypep))&
                         + this%du(iz,m_ind,mp_ind,spin_ind,ipm) * ( f(:,1,l,spin2,atomType) * g(:,1,lp,spin1,atomTypep) &
                                                                    +f(:,2,l,spin2,atomType) * g(:,2,lp,spin1,atomTypep))

            IF(ALLOCATED(this%uulo)) THEN
               iLO_ind = 0
               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  iLO_ind = iLO_ind + 1
                  gmat(:,iz) = gmat(:,iz) &
                              + this%uulo(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ( f(:,1,lp,spin1,atomTypep) *flo(:,1,ilo,spin2,atomType) &
                                                                                   +f(:,2,lp,spin1,atomTypep) *flo(:,2,ilo,spin2,atomType))&
                              + this%dulo(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ( g(:,1,lp,spin1,atomTypep) *flo(:,1,ilo,spin2,atomType) &
                                                                                   +g(:,2,lp,spin1,atomTypep) *flo(:,2,ilo,spin2,atomType))
               ENDDO
               iLO_ind = 0
               DO ilo = 1, atoms%nlo(atomTypep)
                  IF(atoms%llo(ilo,atomTypep).NE.lp) CYCLE
                  iLO_ind = iLO_ind + 1
                  gmat(:,iz) = gmat(:,iz) &
                              + this%ulou(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ( flo(:,1,ilo,spin1,atomTypep)*f(:,1,l,spin2,atomType) &
                                                                                   +flo(:,2,ilo,spin1,atomTypep)*f(:,2,l,spin2,atomType))&
                              + this%ulod(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ( flo(:,1,ilo,spin1,atomTypep)*g(:,1,l,spin2,atomType) &
                                                                                   +flo(:,2,ilo,spin1,atomTypep)*g(:,2,l,spin2,atomType))
               ENDDO
               iLO_ind = 0
               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  iLOp_ind = 0
                  DO ilop = 1, atoms%nlo(atomTypep)
                     IF(atoms%llo(ilop,atomTypep).NE.lp) CYCLE
                     iLOp_ind = iLOp_ind + 1
                     gmat(:,iz) = gmat(:,iz) &
                                 + this%uloulop(iz,m_ind,mp_ind,iLO_ind,iLOp_ind,spin_ind,ipm) *( flo(:,1,ilo,spin2,atomType)*flo(:,1,ilop,spin1,atomTypep) &
                                                                                                 +flo(:,2,ilo,spin2,atomType)*flo(:,2,ilop,spin1,atomTypep))
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         !------------------------
         ! Additional operations
         !------------------------
         !Complex conjugate for spin 4
         IF(spin.EQ.4) gmat = conjg(gmat)

      END SUBROUTINE getRadial_gf

      SUBROUTINE getRadialRadial_gf(this,atoms,iz,m,mp,l_conjg,spin,f,g,flo,gmat)

         !Returns the green's function on the radial and energy mesh (r/=r')
         !for a certain m,mp,spin combination. Attention: The correct radial functions have to be provided

         CLASS(t_greensf),    INTENT(IN)     :: this
         TYPE(t_atoms),       INTENT(IN)     :: atoms
         INTEGER,             INTENT(IN)     :: iz
         INTEGER,             INTENT(IN)     :: m,mp
         LOGICAL,             INTENT(IN)     :: l_conjg
         INTEGER,             INTENT(IN)     :: spin
         REAL   ,             INTENT(IN)     :: f(:,:,0:,:,:)
         REAL   ,             INTENT(IN)     :: g(:,:,0:,:,:)
         REAL   ,             INTENT(IN)     :: flo(:,:,:,:,:)
         COMPLEX, ALLOCATABLE,INTENT(INOUT)  :: gmat(:,:) !Return matrix

         INTEGER spin1,spin2,ipm,spin_ind,m_ind,mp_ind,ilo,ilop,iLO_ind,iLOp_ind
         INTEGER l,lp,atomType,atomTypep,nspins,jr,jrp

         IF(.NOT.this%l_calc) THEN
            CALL juDFT_error("The requested Green's Function element was not calculated", calledby="get_gf")
         ENDIF

         l  = this%elem%l
         lp = this%elem%lp
         atomType  = this%elem%atomType
         atomTypep = this%elem%atomTypep

         IF(ALLOCATED(this%gmmpMat)) THEN
            CALL juDFT_error("Green's function not calculated for radial dependence", calledby="get_gf")
         ENDIF

         nspins = SIZE(this%uu,4)

         IF(spin.GT.4 .OR. spin.LT.1) THEN
            CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
         ENDIF

         ipm = MERGE(2,1,l_conjg)

         IF(.NOT.ALLOCATED(gmat)) ALLOCATE(gmat(SIZE(f,1),SIZE(f,1)),source=cmplx_0)

         IF(spin < 3) THEN
            spin1 = spin
            spin2 = spin
         ELSE IF(spin.EQ.3) THEN
            spin1 = 2
            spin2 = 1
         ELSE
            spin1 = 1
            spin2 = 2
         ENDIF
         !Find the correct spin index in gmmpMat arrays
         spin_ind = MERGE(1,spin,nspins.EQ.1)
         spin_ind = MERGE(3,spin_ind,spin.EQ.4)

         !-------------------------------------------------------------------
         ! Check wether we need to do some operation on the indices m and mp
         !-------------------------------------------------------------------
         IF(spin.EQ.2 .AND. nspins.EQ.1) THEN
            !For a non-spin-polarized calculation we might still want the full
            !matrix. Then we need to reverse the order (SOC prop m*s_z)
            m_ind  = -m
            mp_ind = -mp
         ELSE IF(spin.EQ.4) THEN
            !We only calculate spin21. spin12 is obtained as hermitian conjugate
            !(Complex conjugation happens afterwards)
            m_ind  = mp
            mp_ind = m
         ELSE
            !Do nothing
            m_ind  = m
            mp_ind = mp
         ENDIF
         !-------------------
         ! Fetch the values
         !-------------------
         !$OMP parallel do default(none) &
         !$OMP shared(this,atoms,f,g,flo,gmat,m_ind,mp_ind,spin_ind,ipm,atomType,atomTypep,spin1,spin2,l,lp,iz) &
         !$OMP private(jr,iLO_ind,iLOp_ind,ilo,ilop)
         DO jr = 1, atoms%jri(atomType)
            gmat(:,jr) =  this%uu(iz,m_ind,mp_ind,spin_ind,ipm) * ( f(jr,1,l,spin2,atomType) * f(:,1,lp,spin1,atomTypep) &
                                                                   +f(jr,2,l,spin2,atomType) * f(:,2,lp,spin1,atomTypep))&
                        + this%dd(iz,m_ind,mp_ind,spin_ind,ipm) * ( g(jr,1,l,spin2,atomType) * g(:,1,lp,spin1,atomTypep) &
                                                                   +g(jr,2,l,spin2,atomType) * g(:,2,lp,spin1,atomTypep))&
                        + this%ud(iz,m_ind,mp_ind,spin_ind,ipm) * ( g(jr,1,l,spin2,atomType) * f(:,1,lp,spin1,atomTypep) &
                                                                   +g(jr,2,l,spin2,atomType) * f(:,2,lp,spin1,atomTypep))&
                        + this%du(iz,m_ind,mp_ind,spin_ind,ipm) * ( f(jr,1,l,spin2,atomType) * g(:,1,lp,spin1,atomTypep) &
                                                                   +f(jr,2,l,spin2,atomType) * g(:,2,lp,spin1,atomTypep))

            IF(ALLOCATED(this%uulo)) THEN
               iLO_ind = 0
               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  iLO_ind = iLO_ind + 1
                  gmat(:,jr) = gmat(:,jr) &
                              + this%uulo(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ( f(:,1,lp,spin1,atomTypep) *flo(jr,1,ilo,spin2,atomType) &
                                                                                   +f(:,2,lp,spin1,atomTypep) *flo(jr,2,ilo,spin2,atomType))&
                              + this%dulo(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ( g(:,1,lp,spin1,atomTypep) *flo(jr,1,ilo,spin2,atomType) &
                                                                                   +g(:,2,lp,spin1,atomTypep) *flo(jr,2,ilo,spin2,atomType))
               ENDDO
               iLO_ind = 0
               DO ilo = 1, atoms%nlo(atomTypep)
                  IF(atoms%llo(ilo,atomTypep).NE.lp) CYCLE
                  iLO_ind = iLO_ind + 1
                  gmat(:,jr) = gmat(:,jr) &
                              + this%ulou(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ( flo(:,1,ilo,spin1,atomTypep)*f(jr,1,l,spin2,atomType) &
                                                                                   +flo(:,2,ilo,spin1,atomTypep)*f(jr,2,l,spin2,atomType))&
                              + this%ulod(iz,m_ind,mp_ind,iLO_ind,spin_ind,ipm) * ( flo(:,1,ilo,spin1,atomTypep)*g(jr,1,l,spin2,atomType) &
                                                                                   +flo(:,2,ilo,spin1,atomTypep)*g(jr,2,l,spin2,atomType))
               ENDDO
               iLO_ind = 0
               DO ilo = 1, atoms%nlo(atomType)
                  IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                  iLOp_ind = 0
                  DO ilop = 1, atoms%nlo(atomTypep)
                     IF(atoms%llo(ilop,atomType).NE.lp) CYCLE
                     iLOp_ind = iLOp_ind + 1
                     gmat(:,jr) = gmat(:,jr) &
                                 + this%uloulop(iz,m_ind,mp_ind,iLO_ind,iLOp_ind,spin_ind,ipm) *( flo(jr,1,ilo,spin2,atomType)*flo(:,1,ilop,spin1,atomTypep) &
                                                                                                 +flo(jr,2,ilo,spin2,atomType)*flo(:,2,ilop,spin1,atomTypep))
                  ENDDO
               ENDDO
            ENDIF
            !Get the right normalization
            gmat(:,jr) = gmat(:,jr) * atoms%rmsh(jr,atomType) * atoms%rmsh(:atoms%jri(atomTypep),atomTypep)
         ENDDO
         !$OMP end parallel do

         !------------------------
         ! Additional operations
         !------------------------
         !Complex conjugate for spin 4
         IF(spin.EQ.4) gmat = conjg(gmat)

      END SUBROUTINE getRadialRadial_gf

      SUBROUTINE getRadialSpin_gf(this,atoms,m,mp,l_conjg,f,g,flo,gmat)

         !Returns the green's function on the radial and energy mesh and in a 2x2 spin matrix
         !for a certain m,mp,spin combination. Attention: The correct radial functions have to be provided

         CLASS(t_greensf),    INTENT(IN)     :: this
         TYPE(t_atoms),       INTENT(IN)     :: atoms
         INTEGER,             INTENT(IN)     :: m,mp
         LOGICAL,             INTENT(IN)     :: l_conjg
         REAL   ,             INTENT(IN)     :: f(:,:,0:,:,:)
         REAL   ,             INTENT(IN)     :: g(:,:,0:,:,:)
         REAL   ,             INTENT(IN)     :: flo(:,:,:,:,:)
         COMPLEX, ALLOCATABLE,INTENT(INOUT)  :: gmat(:,:,:,:) !Return matrix

         INTEGER :: spin,spin1,spin2
         COMPLEX,ALLOCATABLE :: temp(:,:)

         IF(.NOT.ALLOCATED(gmat)) ALLOCATE(gmat(2,2,SIZE(f,1),this%contour%nz),source=cmplx_0)

         DO spin = 1, 4
            IF(spin < 3) THEN
               spin1 = spin
               spin2 = spin
            ELSE IF(spin.EQ.3) THEN
               spin1 = 2
               spin2 = 1
            ELSE
               spin1 = 1
               spin2 = 2
            ENDIF
            IF(spin>=3 .AND.SIZE(this%uu,4)<3) THEN
               gmat(spin1,spin2,:,:) = cmplx_0
            ELSE
               CALL this%getRadial(atoms,m,mp,l_conjg,spin,f,g,flo,temp)
               gmat(spin1,spin2,:,:) = temp(:,:)
            ENDIF
         ENDDO

      END SUBROUTINE getRadialSpin_gf

      SUBROUTINE getRadialRadialSpin_gf(this,atoms,iz,m,mp,l_conjg,f,g,flo,gmat)

         !Returns the green's function on the radial and energy mesh and in a 2x2 spin matrix
         !for a certain m,mp,spin combination. Attention: The correct radial functions have to be provided

         CLASS(t_greensf),    INTENT(IN)     :: this
         TYPE(t_atoms),       INTENT(IN)     :: atoms
         INTEGER,             INTENT(IN)     :: iz
         INTEGER,             INTENT(IN)     :: m,mp
         LOGICAL,             INTENT(IN)     :: l_conjg
         REAL   ,             INTENT(IN)     :: f(:,:,0:,:,:)
         REAL   ,             INTENT(IN)     :: g(:,:,0:,:,:)
         REAL   ,             INTENT(IN)     :: flo(:,:,:,:,:)
         COMPLEX, ALLOCATABLE,INTENT(INOUT)  :: gmat(:,:,:,:) !Return matrix

         INTEGER :: spin,spin1,spin2
         COMPLEX,ALLOCATABLE :: temp(:,:)

         IF(.NOT.ALLOCATED(gmat)) ALLOCATE(gmat(2,2,SIZE(f,1),SIZE(f,1)),source=cmplx_0)

         DO spin = 1, 4
            IF(spin < 3) THEN
               spin1 = spin
               spin2 = spin
            ELSE IF(spin.EQ.3) THEN
               spin1 = 2
               spin2 = 1
            ELSE
               spin1 = 1
               spin2 = 2
            ENDIF
            IF(spin>=3 .AND.SIZE(this%uu,4)<3) THEN
               gmat(spin1,spin2,:,:) = cmplx_0
            ELSE
               CALL this%getRadialRadial(atoms,iz,m,mp,l_conjg,spin,f,g,flo,temp)
               gmat(spin1,spin2,:,:) = temp(:,:)
            ENDIF
         ENDDO

      END SUBROUTINE getRadialRadialSpin_gf

      SUBROUTINE set_gf(this,iz,l_conjg,gmat,spin)

         USE m_types_mat

         !Sets the spherically averaged greens function matrix belonging to energy point iz with l,lp,nType,nTypep
         !equal to gmat

         CLASS(t_greensf),    INTENT(INOUT)  :: this
         INTEGER,             INTENT(IN)     :: iz
         LOGICAL,             INTENT(IN)     :: l_conjg
         TYPE(t_mat),         INTENT(IN)     :: gmat
         INTEGER, OPTIONAL,   INTENT(IN)     :: spin

         INTEGER matsize1,matsize2,ind1,ind2,ind1_start,ind2_start
         INTEGER l,lp,atomType,atomTypep,m,mp,spin1,spin2,ipm,ispin,spin_start,spin_end
         LOGICAL l_full


         this%l_calc = .TRUE. !If its set it counts as calculated

         l  = this%elem%l
         lp = this%elem%lp
         atomType  = this%elem%atomType
         atomTypep = this%elem%atomTypep


         IF(ALLOCATED(this%uu)) THEN
            CALL juDFT_error("Can only set spherically averaged Green's Function", calledby="set_gf")
         ENDIF

         IF(PRESENT(spin)) THEN
            IF(spin.GT.4 .OR. spin.LT.1) THEN
               CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
            ENDIF
         ENDIF

         l_full = .NOT.PRESENT(spin)
         !Determine matsize for the result gmat (if spin is given only return this digonal element)
         matsize1 = (2*l+1) * MERGE(2,1,l_full)
         matsize2 = (2*lp+1) * MERGE(2,1,l_full)

         !Check the expected matsizes against the actual
         IF(matsize1.NE.gmat%matsize1.OR.matsize2.NE.gmat%matsize2) THEN
            CALL juDFT_error("Mismatch in matsizes", calledby="set_gf")
         ENDIF

         ipm = MERGE(2,1,l_conjg)

         IF(l_full) THEN
            spin_start = 1
            spin_end   = SIZE(this%gmmpMat,4)
         ELSE
            spin_start = spin
            spin_end   = spin
         ENDIF

         DO ispin = spin_start, spin_end
            !Find the right quadrant in gmat according to the spin index
            IF(ispin.EQ.2 .AND.SIZE(this%gmmpMat,4).EQ.1) CYCLE
            IF(l_full) THEN
               IF(ispin < 3) THEN
                  spin1 = ispin
                  spin2 = ispin
               ELSE IF(ispin.EQ.3) THEN
                  spin1 = 2
                  spin2 = 1
               ELSE
                  spin1 = 1
                  spin2 = 2
               ENDIF
               ind1_start = (spin2-1)*(2*l+1)
               ind2_start = (spin1-1)*(2*lp+1)
            ELSE
               ind1_start = 0
               ind2_start = 0
            ENDIF
            ind1 = ind1_start
            DO m = -l,l
               ind1 = ind1 + 1
               ind2 = ind2_start
               DO mp = -lp,lp
                  ind2 = ind2 + 1
                  this%gmmpMat(iz,m,mp,ispin,ipm) = gmat%data_c(ind1,ind2)
                  IF(l_full) this%gmmpMat(iz,m,mp,ispin,ipm) = this%gmmpMat(iz,m,mp,ispin,ipm) &
                                                               * MERGE(2.0,1.0,SIZE(this%gmmpMat,4).EQ.1)
               ENDDO
            ENDDO
         ENDDO

      END SUBROUTINE set_gf

      SUBROUTINE rotate_gf(this,sym,atoms)

         !Applies the given symmetry operation to the greens function
         CLASS(t_greensf),    INTENT(INOUT)  :: this
         TYPE(t_sym),         INTENT(IN)     :: sym
         TYPE(t_atoms),       INTENT(IN)     :: atoms

         INTEGER :: l,lp,iop,atomType,atomTypep,nspins
         INTEGER :: ipm,ispin,iz,ilo,ilop,iLO_ind,iLOp_ind

         IF(this%elem%representative_elem<0) RETURN !Nothing to be done

         CALL timestart("Green's Function: Rotate")
         l  = this%elem%l
         lp = this%elem%lp
         atomType = this%elem%atomType
         atomTypep = this%elem%atomTypep
         iop = this%elem%representative_op

         nspins = 0
         IF(ALLOCATED(this%gmmpMat)) nspins = SIZE(this%gmmpMat,4)
         IF(ALLOCATED(this%uu)) nspins = SIZE(this%uu,4)

         DO ipm = 1, 2
            DO ispin = 1, nspins
               DO iz = 1, this%contour%nz
                  IF(ALLOCATED(this%gmmpMat)) THEN
                     this%gmmpMat(iz,:,:,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                             this%gmmpMat(iz,:,:,ispin,ipm))
                     this%gmmpMat(iz,:,:,ispin,ipm) = matmul(this%gmmpMat(iz,:,:,ispin,ipm),&
                                                             sym%d_wgn(:,:,lp,iop))
                  ELSE IF(ALLOCATED(this%uu)) THEN
                     this%uu(iz,:,:,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                        this%uu(iz,:,:,ispin,ipm))
                     this%uu(iz,:,:,ispin,ipm) = matmul(this%uu(iz,:,:,ispin,ipm),&
                                                        sym%d_wgn(:,:,lp,iop))
                     this%dd(iz,:,:,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                        this%dd(iz,:,:,ispin,ipm))
                     this%dd(iz,:,:,ispin,ipm) = matmul(this%dd(iz,:,:,ispin,ipm),&
                                                        sym%d_wgn(:,:,lp,iop))
                     this%ud(iz,:,:,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                        this%ud(iz,:,:,ispin,ipm))
                     this%ud(iz,:,:,ispin,ipm) = matmul(this%ud(iz,:,:,ispin,ipm),&
                                                        sym%d_wgn(:,:,lp,iop))
                     this%du(iz,:,:,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                        this%du(iz,:,:,ispin,ipm))
                     this%du(iz,:,:,ispin,ipm) = matmul(this%du(iz,:,:,ispin,ipm),&
                                                        sym%d_wgn(:,:,lp,iop))

                     IF(ALLOCATED(this%uulo)) THEN
                        iLO_ind = 0
                        DO ilo = 1, atoms%nlo(atomType)
                           IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                           iLO_ind = iLO_ind + 1
                           this%uulo(iz,:,:,iLO_ind,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                                        this%uulo(iz,:,:,iLO_ind,ispin,ipm))
                           this%uulo(iz,:,:,iLO_ind,ispin,ipm) = matmul(this%uulo(iz,:,:,iLO_ind,ispin,ipm),&
                                                                        sym%d_wgn(:,:,lp,iop))
                           this%dulo(iz,:,:,iLO_ind,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                                        this%dulo(iz,:,:,iLO_ind,ispin,ipm))
                           this%dulo(iz,:,:,iLO_ind,ispin,ipm) = matmul(this%dulo(iz,:,:,iLO_ind,ispin,ipm),&
                                                                        sym%d_wgn(:,:,lp,iop))
                        ENDDO
                        iLO_ind = 0
                        DO ilo = 1, atoms%nlo(atomTypep)
                           IF(atoms%llo(ilo,atomTypep).NE.lp) CYCLE
                           iLO_ind = iLO_ind + 1
                           this%ulou(iz,:,:,iLO_ind,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                                        this%ulou(iz,:,:,iLO_ind,ispin,ipm))
                           this%ulou(iz,:,:,iLO_ind,ispin,ipm) = matmul(this%ulou(iz,:,:,iLO_ind,ispin,ipm),&
                                                                        sym%d_wgn(:,:,lp,iop))
                           this%ulod(iz,:,:,iLO_ind,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                                        this%ulod(iz,:,:,iLO_ind,ispin,ipm))
                           this%ulod(iz,:,:,iLO_ind,ispin,ipm) = matmul(this%ulod(iz,:,:,iLO_ind,ispin,ipm),&
                                                                        sym%d_wgn(:,:,lp,iop))
                        ENDDO
                        iLO_ind = 0
                        DO ilo = 1, atoms%nlo(atomType)
                           IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                           iLOp_ind = 0
                           DO ilop = 1, atoms%nlo(atomTypep)
                              IF(atoms%llo(ilop,atomType).NE.lp) CYCLE
                              iLOp_ind = iLOp_ind + 1
                              this%uloulop(iz,:,:,iLO_ind,iLOp_ind,ispin,ipm) = matmul(conjg(transpose(sym%d_wgn(:,:,l,iop))),&
                                                                                       this%uloulop(iz,:,:,iLO_ind,iLOp_ind,ispin,ipm))
                              this%uloulop(iz,:,:,iLO_ind,iLOp_ind,ispin,ipm) = matmul(this%uloulop(iz,:,:,iLO_ind,iLOp_ind,ispin,ipm) ,&
                                                                                       sym%d_wgn(:,:,lp,iop))
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         CALL timestop("Green's Function: Rotate")
      END SUBROUTINE rotate_gf

      SUBROUTINE reset_gf(this)

         !---------------------------------------------------
         ! Sets all gmmpMat arrays back to 0
         !---------------------------------------------------

         CLASS(t_greensf),       INTENT(INOUT)  :: this

         IF(ALLOCATED(this%gmmpMat)) this%gmmpMat = cmplx_0
         IF(ALLOCATED(this%uu)) THEN
            this%uu = cmplx_0
            this%ud = cmplx_0
            this%du = cmplx_0
            this%dd = cmplx_0
         ENDIF
         IF(ALLOCATED(this%uulo)) THEN
            this%uulo = cmplx_0
            this%ulou = cmplx_0
            this%dulo = cmplx_0
            this%ulod = cmplx_0

            this%uloulop = cmplx_0
         ENDIF

      END SUBROUTINE reset_gf

      PURE FUNCTION checkEmpty_greensf(this,m,mp,spin,ipm) Result(l_empty)

         CLASS(t_greensf),         INTENT(IN)   :: this
         INTEGER,                  INTENT(IN)   :: m
         INTEGER,                  INTENT(IN)   :: mp
         INTEGER,                  INTENT(IN)   :: spin
         INTEGER,                  INTENT(IN)   :: ipm

         LOGICAL :: l_empty

         IF(ALLOCATED(this%gmmpMat)) THEN
            l_empty = ALL(ABS(this%gmmpMat(:,m,mp,spin,ipm)).LT.1e-12)
         ELSE
            l_empty =      ALL(ABS(this%uu(:,m,mp,spin,ipm)).LT.1e-12) &
                     .AND. ALL(ABS(this%dd(:,m,mp,spin,ipm)).LT.1e-12) &
                     .AND. ALL(ABS(this%ud(:,m,mp,spin,ipm)).LT.1e-12) &
                     .AND. ALL(ABS(this%du(:,m,mp,spin,ipm)).LT.1e-12)
            IF(ALLOCATED(this%uulo)) THEN
               l_empty = l_empty .AND. ALL(ABS(this%uulo(:,m,mp,:,spin,ipm)).LT.1e-12) &
                                 .AND. ALL(ABS(this%ulou(:,m,mp,:,spin,ipm)).LT.1e-12) &
                                 .AND. ALL(ABS(this%ulod(:,m,mp,:,spin,ipm)).LT.1e-12) &
                                 .AND. ALL(ABS(this%dulo(:,m,mp,:,spin,ipm)).LT.1e-12) &
                                 .AND. ALL(ABS(this%uloulop(:,m,mp,:,:,spin,ipm)).LT.1e-12)
            ENDIF
         ENDIF

      END FUNCTION checkEmpty_greensf

      SUBROUTINE resetSingleElem_gf(this,m,mp,spin,ipm)

         !---------------------------------------------------
         ! Sets one Element in gmmpMat arrays back to 0
         !---------------------------------------------------

         CLASS(t_greensf),       INTENT(INOUT)  :: this
         INTEGER,                INTENT(IN)     :: m
         INTEGER,                INTENT(IN)     :: mp
         INTEGER,                INTENT(IN)     :: spin
         INTEGER,                INTENT(IN)     :: ipm

         IF(ALLOCATED(this%gmmpMat)) this%gmmpMat(:,m,mp,spin,ipm) = cmplx_0
         IF(ALLOCATED(this%uu)) THEN
            this%uu(:,m,mp,spin,ipm) = cmplx_0
            this%ud(:,m,mp,spin,ipm) = cmplx_0
            this%du(:,m,mp,spin,ipm) = cmplx_0
            this%dd(:,m,mp,spin,ipm) = cmplx_0
         ENDIF
         IF(ALLOCATED(this%uulo)) THEN
            this%uulo(:,m,mp,:,spin,ipm) = cmplx_0
            this%ulou(:,m,mp,:,spin,ipm) = cmplx_0
            this%dulo(:,m,mp,:,spin,ipm) = cmplx_0
            this%ulod(:,m,mp,:,spin,ipm) = cmplx_0

            this%uloulop(:,m,mp,:,:,spin,ipm) = cmplx_0
         ENDIF

      END SUBROUTINE resetSingleElem_gf

      FUNCTION integrateOverMT_greensf(this,atoms,input,gfinp,f,g,flo,usdus,denCoeffsOffDiag,scalarGF,l_fullRadial) Result(gIntegrated)

         USE m_intgr
         USE m_types_scalarGF
         USE m_types_usdus
         USE m_types_denCoeffsOffDiag
         USE m_types_mat

         CLASS(t_greensf),                   INTENT(IN) :: this
         TYPE(t_atoms),                      INTENT(IN) :: atoms
         TYPE(t_input),                      INTENT(IN) :: input
         TYPE(t_gfinp),                      INTENT(IN) :: gfinp
         REAL,                               INTENT(IN) :: f(:,:,0:,:,:)
         REAL,                               INTENT(IN) :: g(:,:,0:,:,:)
         REAL,                               INTENT(IN) :: flo(:,:,:,:,:)
         TYPE(t_usdus),            OPTIONAL, INTENT(IN) :: usdus
         TYPE(t_denCoeffsOffDiag), OPTIONAL, INTENT(IN) :: denCoeffsOffDiag
         TYPE(t_scalarGF),         OPTIONAL, INTENT(IN) :: scalarGF
         LOGICAL,                  OPTIONAL, INTENT(IN) :: l_fullRadial

         TYPE(t_greensf) :: gIntegrated

         LOGICAL :: l_fullRadialArg,l_explicit
         INTEGER :: l,lp,atomType,atomTypep,ipm,spin,m,mp,iz,jr,jrp
         REAL    :: realPart, imagPart, atomDiff(3)
         COMPLEX, ALLOCATABLE :: gmatR(:,:)
         COMPLEX :: gmat(atoms%jmtd)
         TYPE(t_mat) :: gmatTmp

         l_fullRadialArg = .FALSE.
         IF(PRESENT(l_fullRadial)) l_fullRadialArg = l_fullRadial

         IF(this%elem%l_sphavg) CALL juDFT_error("GF has to be provided with radial dependence",&
                                                 calledby="integrateOverMT_greensf")


         CALL timestart("Green's Function: Average over MT")
         CALL gIntegrated%init(this%elem,gfinp,atoms,input,contour_in=this%contour,l_sphavg_in=.TRUE.)
         gIntegrated%l_calc = .TRUE.
         l  = this%elem%l
         lp = this%elem%lp
         atomType  = this%elem%atomType
         atomTypep = this%elem%atomTypep
         atomDiff  = this%elem%atomDiff
         !Do we have the offdiagonal scalar products
         l_explicit = .TRUE.
         IF(PRESENT(scalarGF)) THEN
            IF(scalarGF%done.AND.this%elem%isOffDiag()) l_explicit = .FALSE.
         ENDIF
         IF(PRESENT(usdus).OR.PRESENT(denCoeffsOffDiag).AND..NOT.this%elem%isOffDiag()) l_explicit = .FALSE.

         !only intersite arguments have independent radial arguments ??
         l_fullRadialArg = l_fullRadialArg.AND.(atomType.NE.atomTypep.OR.ANY(ABS(atomDiff).GT.1e-12))

         DO ipm = 1, 2
            DO spin = 1 , SIZE(this%uu,4)
               IF(.NOT.l_explicit) THEN
                  DO iz = 1, this%contour%nz
                     CALL this%get(atoms,iz,ipm==2,gmatTmp,spin=spin,usdus=usdus,&
                                   denCoeffsOffDiag=denCoeffsOffDiag,scalarGF=scalarGF)
                     CALL gIntegrated%set(iz,ipm==2,gmatTmp,spin=spin)
                  ENDDO
               ELSE
                  DO mp = -lp, lp
                     DO m = -l, l
                        IF(this%checkEmpty(m,mp,spin,ipm)) CYCLE
                        IF(.NOT.l_fullRadialArg) THEN
                           CALL this%getRadial(atoms,m,mp,ipm==2,spin,f,g,flo,gmatR)
                        ENDIF
                        DO iz = 1, this%contour%nz
                           IF(l_fullRadialArg) THEN
                              CALL this%getRadialRadial(atoms,iz,m,mp,ipm==2,spin,f,g,flo,gmatR)
                              gmat = cmplx_0
                              DO jr = 1, SIZE(gmat)
                                 CALL intgr3(REAL(gmatR(:,jr)),atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep),realPart)
                                 CALL intgr3(AIMAG(gmatR(:,jr)),atoms%rmsh(:,atomTypep),atoms%dx(atomTypep),atoms%jri(atomTypep),imagPart)

                                 gmat(jr) = realPart + ImagUnit * imagPart
                              ENDDO
                           ELSE
                              gmat = gmatR(:,iz)
                           ENDIF
                           CALL intgr3(REAL(gmat),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realPart)
                           CALL intgr3(AIMAG(gmat),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),imagPart)

                           gIntegrated%gmmpMat(iz,m,mp,spin,ipm) = realPart + ImagUnit * imagPart
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CALL timestop("Green's Function: Average over MT")

      END FUNCTION integrateOverMT_greensf

END MODULE m_types_greensf
