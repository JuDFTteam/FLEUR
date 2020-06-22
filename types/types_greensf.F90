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
   !>  We have the following cases
   !>    -onsite
   !>       -we look at l=l' but m\=m'
   !>       -we treat non-magnetic/collinear and noco (not tested)
   !>       -we look at r=r' and spherically averaged gf
   !>    -intersite
   !>       -l\=l' and m\=m'
   !>       -r\=r' (not stored we calculate the gf by calling calc_intersite in m_intersite for specific r and r')
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

      CONTAINS
         PROCEDURE, PASS :: init    => init_greensf
         PROCEDURE       :: mpi_bc  => mpi_bc_greensf
         PROCEDURE       :: collect => collect_greensf
         PROCEDURE       :: get     => get_gf
         PROCEDURE       :: set     => set_gf
         PROCEDURE       :: reset   => reset_gf
   END TYPE t_greensf

   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE init_greensf(this,gfelem,gfinp,input,contour_in)

         CLASS(t_greensf),             INTENT(INOUT)  :: this
         TYPE(t_gfelementtype), TARGET,INTENT(IN)     :: gfelem
         TYPE(t_gfinp),                INTENT(IN)     :: gfinp
         TYPE(t_input),                INTENT(IN)     :: input
         !Pass a already calculated energy contour to the type
         TYPE(t_greensfContourData), OPTIONAL, INTENT(IN)   :: contour_in

         INTEGER spin_dim,lmax

         this%elem => gfelem

         !Initialize the contour
         CALL this%contour%init(gfinp%contour(this%elem%iContour),contour_in=contour_in)

         spin_dim = MERGE(3,input%jspins,gfinp%l_mperp)
         lmax = lmaxU_const

         IF(gfinp%l_sphavg) THEN
            ALLOCATE(this%gmmpMat(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
         ELSE
            ALLOCATE(this%uu(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
            ALLOCATE(this%dd(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
            ALLOCATE(this%du(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
            ALLOCATE(this%ud(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
         ENDIF

      END SUBROUTINE init_greensf

      SUBROUTINE mpi_bc_greensf(this,mpi_comm,irank)
         USE m_mpi_bc_tool
         CLASS(t_greensf), INTENT(INOUT)::this
         INTEGER, INTENT(IN):: mpi_comm
         INTEGER, INTENT(IN), OPTIONAL::irank
         INTEGER ::rank,myrank,n,ierr
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

      END SUBROUTINE mpi_bc_greensf

      SUBROUTINE collect_greensf(this,mpi_comm)

         CLASS(t_greensf),     INTENT(INOUT) :: this
         INTEGER,              INTENT(IN)    :: mpi_comm
#ifdef CPP_MPI
         include 'mpif.h'
#include"cpp_double.h"
         INTEGER:: ierr,n
         COMPLEX,ALLOCATABLE::ctmp(:)

         IF(ALLOCATED(this%gmmpMat)) THEN
            n = SIZE(this%gmmpMat)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%gmmpMat,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_comm,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%gmmpMat,1)
            DEALLOCATE(ctmp)
         ELSE
            n = SIZE(this%uu)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%uu,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_comm,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%uu,1)
            CALL MPI_ALLREDUCE(this%ud,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_comm,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%ud,1)
            CALL MPI_ALLREDUCE(this%du,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_comm,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%du,1)
            CALL MPI_ALLREDUCE(this%dd,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_comm,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%dd,1)
            DEALLOCATE(ctmp)
         ENDIF
#endif
      END SUBROUTINE collect_greensf

      SUBROUTINE get_gf(this,iz,l_conjg,gmat,spin,u,udot)

         USE m_types_mat

         !Returns the matrix belonging to energy point iz with l,lp,nType,nTypep
         !when jr (and jrp) are given return for that radial point

         CLASS(t_greensf),    INTENT(IN)     :: this
         INTEGER,             INTENT(IN)     :: iz
         LOGICAL,             INTENT(IN)     :: l_conjg
         TYPE(t_mat),         INTENT(INOUT)  :: gmat !Return matrix
         INTEGER, OPTIONAL,   INTENT(IN)     :: spin
         REAL   , OPTIONAL,   INTENT(IN)     :: u(:,:)       !Radial functions at the point where you want to evaluate the greens function
         REAL   , OPTIONAL,   INTENT(IN)     :: udot(:,:)

         INTEGER matsize1,matsize2,i,j,ind1,ind2,ind1_start,ind2_start
         INTEGER m,mp,spin1,spin2,ipm,ispin,spin_start,spin_end,spin_ind,m_ind,mp_ind
         INTEGER l,lp,atomType,atomTypep
         LOGICAL l_radial,l_full

         IF(.NOT.this%l_calc) THEN
            CALL juDFT_error("The requested Green's Function element was not calculated", calledby="get_gf")
         ENDIF

         l  = this%elem%l
         lp = this%elem%lp
         atomType  = this%elem%atomType
         atomTypep = this%elem%atomTypep

         IF(PRESENT(u).OR.PRESENT(udot).AND.ALLOCATED(this%gmmpMat)) THEN
            CALL juDFT_error("Green's function not calculated for radial dependence", calledby="get_gf")
         ENDIF

         IF((PRESENT(u).AND..NOT.PRESENT(udot)).OR.&
            (PRESENT(udot).AND..NOT.PRESENT(u))) THEN
            CALL juDFT_error("Not a valid input: Either provide both u and udot or neither of them", calledby="get_gf")
         ENDIF

         l_radial = PRESENT(u).AND.PRESENT(udot)

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
            spin_end   = MERGE(4,2,SIZE(this%gmmpMat,4).EQ.3)
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
            spin_ind = MERGE(1,ispin,SIZE(this%gmmpMat,4).EQ.1)
            spin_ind = MERGE(3,spin_ind,ispin.EQ.4)
            !Find the right quadrant in gmat
            IF(l_full) THEN
               ind1_start = (spin1-1)*(2*l+1)
               ind2_start = (spin2-1)*(2*lp+1)
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

                  !-------------------------------------------------------------------
                  ! Check wether we need to do some operation on the indices m and mp
                  !-------------------------------------------------------------------
                  IF(ispin.EQ.2 .AND. SIZE(this%gmmpMat,4).EQ.1) THEN
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
                  IF(l_radial) THEN
                     gmat%data_c(ind1,ind2) = this%uu(iz,m_ind,mp_ind,spin_ind,ipm) * u(1,spin1)    * u(2,spin2)     + &
                                              this%dd(iz,m_ind,mp_ind,spin_ind,ipm) * udot(1,spin1) * udot(2,spin2)  + &
                                              this%du(iz,m_ind,mp_ind,spin_ind,ipm) * udot(1,spin1) * u(2,spin2)     + &
                                              this%ud(iz,m_ind,mp_ind,spin_ind,ipm) * u(1,spin1)    * udot(2,spin2)
                  ELSE
                     gmat%data_c(ind1,ind2) = this%gmmpMat(iz,m_ind,mp_ind,spin_ind,ipm)
                  ENDIF
                  !------------------------
                  ! Additional operations
                  !------------------------
                  !Spin-degeneracy when using a full matrix and having input%jspins.EQ.1
                  IF(l_full) gmat%data_c(ind1,ind2) = gmat%data_c(ind1,ind2) * MERGE(0.5,1.0,SIZE(this%gmmpMat,4).EQ.1)
                  !Complex conjugate for spin 4
                  IF(ispin.EQ.4) gmat%data_c(ind1,ind2) = conjg(gmat%data_c(ind1,ind2))

               ENDDO!mp
            ENDDO!m
         ENDDO!ispin

      END SUBROUTINE get_gf

      SUBROUTINE set_gf(this,iz,l_conjg,gmat,spin)

         USE m_types_mat

         !Sets the spherically averaged greens function matrix belonging to energy point iz with l,lp,nType,nTypep
         !equal to gmat

         CLASS(t_greensf),    INTENT(INOUT)  :: this
         INTEGER,             INTENT(IN)     :: iz
         LOGICAL,             INTENT(IN)     :: l_conjg
         TYPE(t_mat),         INTENT(IN)     :: gmat
         INTEGER, OPTIONAL,   INTENT(IN)     :: spin

         INTEGER matsize1,matsize2,i,j,ind1,ind2,ind1_start,ind2_start
         INTEGER l,lp,atomType,atomTypep,m,mp,spin1,spin2,ipm,ispin,spin_start,spin_end
         LOGICAL l_full


         this%l_calc = .TRUE. !If its set it counts as calculated

         l  = this%elem%l
         lp = this%elem%lp
         atomType  = this%elem%atomType
         atomTypep = this%elem%atomTypep

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
               ind1_start = (spin1-1)*(2*l+1)
               ind2_start = (spin2-1)*(2*lp+1)
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

      END SUBROUTINE reset_gf

END MODULE m_types_greensf
