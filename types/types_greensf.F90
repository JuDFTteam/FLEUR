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
   USE m_types_scalarGF
   USE m_types_nococonv

   IMPLICIT NONE

   PRIVATE

   TYPE t_greensf

      LOGICAL :: l_calc = .FALSE.
      LOGICAL :: l_sphavg
      LOGICAL :: l_kresolved

      !Energy contour parameters
      TYPE(t_greensfContourData) :: contour

      !Pointer to the element type in gfinp
      TYPE(t_gfelementtype), POINTER :: elem => NULL()
      TYPE(t_scalarGF) :: scalarProducts

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

      COMPLEX, ALLOCATABLE :: gmmpMat_k(:,:,:,:,:,:)

      CONTAINS
         PROCEDURE, PASS :: init                   => init_greensf
         PROCEDURE       :: mpi_bc                 => mpi_bc_greensf
         PROCEDURE       :: collect                => collect_greensf
         PROCEDURE       :: get                    => get_gf
         PROCEDURE       :: occmtx_greensf_spin
         PROCEDURE       :: occmtx_greensf_complete
         GENERIC         :: occupationMatrix       => occmtx_greensf_spin, occmtx_greensf_complete
         PROCEDURE       :: getFullMatrix          => getFullMatrix_gf
         PROCEDURE       :: getRadial              => getRadial_gf
         PROCEDURE       :: getRadialRadial        => getRadialRadial_gf!(Full Radial dependence for intersite)
         PROCEDURE       :: integrateOverMT        => integrateOverMT_greensf
         PROCEDURE       :: set                    => set_gf
         PROCEDURE       :: set_gfdata             => set_gfdata
         PROCEDURE       :: rotate                 => rotate_gf
         PROCEDURE       :: rotate_euler_angles    => rotate_euler_angles_gf
         PROCEDURE       :: reset                  => reset_gf
         PROCEDURE       :: resetSingleElem        => resetSingleElem_gf
         PROCEDURE       :: checkEmpty             => checkEmpty_greensf
   END TYPE t_greensf

   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE init_greensf(this,gfelem,gfinp,atoms,input,contour_in,l_sphavg_in,l_kresolved_in,nkpt)

         CLASS(t_greensf),                     INTENT(INOUT)  :: this
         TYPE(t_gfelementtype),      TARGET,   INTENT(IN)     :: gfelem
         TYPE(t_gfinp),                        INTENT(IN)     :: gfinp
         TYPE(t_atoms),                        INTENT(IN)     :: atoms
         TYPE(t_input),                        INTENT(IN)     :: input
         !Pass a already calculated energy contour to the type
         TYPE(t_greensfContourData), OPTIONAL, INTENT(IN)     :: contour_in
         LOGICAL,                    OPTIONAL, INTENT(IN)     :: l_sphavg_in !To overwrite the allocation for integrateOverMT_greensf
         LOGICAL,                    OPTIONAL, INTENT(IN)     :: l_kresolved_in !To overwrite the allocation for integrateOverBZ_greensf
         INTEGER,                    OPTIONAL, INTENT(IN)     :: nkpt

         INTEGER spin_dim,lmax,nLO
         LOGICAL l_sphavg

         this%elem => gfelem
         this%l_calc = .FALSE.

         nLO = this%elem%countLOs(atoms)

         !Initialize the contour
         CALL this%contour%init(gfinp%contour(this%elem%iContour),contour_in=contour_in)

         spin_dim = MERGE(3,input%jspins,gfinp%l_mperp)
         lmax = lmaxU_const

         this%l_sphavg = this%elem%l_sphavg
         IF(PRESENT(l_sphavg_in)) this%l_sphavg = l_sphavg_in
      
         this%l_kresolved = this%elem%l_kresolved
         IF(PRESENT(l_kresolved_in)) this%l_kresolved = l_kresolved_in

         IF(this%l_kresolved.AND. .NOT.PRESENT(nkpt)) CALL juDFT_error("Missing argument nkpt for k-resolved green's function", calledby='init_greensf')

         IF(this%l_sphavg) THEN
            IF(this%l_kresolved) THEN
               ALLOCATE(this%gmmpMat_k(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2,nkpt),source=cmplx_0)
            ELSE
               ALLOCATE(this%gmmpMat(this%contour%nz,-lmax:lmax,-lmax:lmax,spin_dim,2),source=cmplx_0)
            ENDIF
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
         CALL mpi_bc(this%l_sphavg,rank,mpi_comm)
         CALL mpi_bc(this%l_kresolved,rank,mpi_comm)

         CALL this%contour%mpi_bc(mpi_comm,rank)

         IF(ALLOCATED(this%gmmpMat)) CALL mpi_bc(this%gmmpMat,rank,mpi_comm)
         IF(ALLOCATED(this%gmmpMat_k)) CALL mpi_bc(this%gmmpMat_k,rank,mpi_comm)
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
         ELSE IF(ALLOCATED(this%gmmpMat_k)) THEN
            n = SIZE(this%gmmpMat_k)
            ALLOCATE(ctmp(n))
            CALL MPI_ALLREDUCE(this%gmmpMat_k,ctmp,n,CPP_MPI_COMPLEX,MPI_SUM,mpi_communicator,ierr)
            CALL CPP_BLAS_ccopy(n,ctmp,1,this%gmmpMat_k,1)
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


      FUNCTION occmtx_greensf_complete(this,gfinp,input,atoms,noco,nococonv,l_write,check,occError) Result(occmtx)

         !calculates the occupation of the greens function
         !The Greens-function should already be prepared on a energy contour ending at e_fermi
         !The occupation is calculated with:
         !
         ! n^sigma_mm' = -1/2pi int^Ef dz (G^+(z)^sigma_mm'-G^-(z)^sigma_mm')
         !
         ! If l_write is given the density matrix together with the spin up/down trace is written to the out files

         CLASS(t_greensf),                 INTENT(IN)    :: this
         TYPE(t_gfinp),                    INTENT(IN)    :: gfinp
         TYPE(t_input),                    INTENT(IN)    :: input
         TYPE(t_atoms),                    INTENT(IN)    :: atoms
         TYPE(t_noco),                     INTENT(IN)    :: noco
         TYPE(t_nococonv),                 INTENT(IN)    :: nococonv
         LOGICAL,                 OPTIONAL,INTENT(IN)    :: l_write !write the occupation matrix to out file
         LOGICAL,                 OPTIONAL,INTENT(IN)    :: check
         LOGICAL,                 OPTIONAL,INTENT(INOUT) :: occError

         COMPLEX, ALLOCATABLE :: occmtx(:,:,:)

         INTEGER :: nspins, ispin, m, mp
         LOGICAL :: all_occError, single_occError
         REAL    :: nup, ndwn
         COMPLEX :: offd
         CHARACTER(len=2) :: l_type
         CHARACTER(len=8) :: l_form
         TYPE(t_contourInp) :: contourInp

         contourInp = gfinp%contour(this%elem%iContour)

         IF(ALLOCATED(this%gmmpMat)) nspins = SIZE(this%gmmpMat,4)
         IF(ALLOCATED(this%uu)) nspins = SIZE(this%uu,4)

         ALLOCATE(occmtx(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const, nspins), source=cmplx_0)

         all_occError = .FALSE.
         DO ispin = 1, nspins
            occmtx(:,:,ispin) = this%occmtx_greensf_spin(ispin,gfinp,input,atoms,noco,nococonv,&
                                                         check=check,occError=single_occError)
            all_occError = all_occError.OR.single_occError
         ENDDO
         IF(PRESENT(occError)) occError = all_occError

         !Io-part (ATM this subroutine is only called from rank 0)
         IF(PRESENT(l_write)) THEN
            IF(l_write) THEN
               !Write to file
               WRITE (l_type,'(i2)') 2*(2*this%elem%l+1)
               l_form = '('//l_type//'f8.4)'
               IF(this%elem%isIntersite()) THEN
8998              FORMAT(/,"Occupation matrix obtained from the green's function for atom: ",I3," atomp: ",I3," atomDiff: ", 3f8.4," l: ",I3," lp: ",I3)
                  WRITE(oUnit,8998) this%elem%atomType, this%elem%atomTypep, this%elem%atomDiff, this%elem%l, this%elem%lp
               ELSE IF (this%elem%isOffDiag()) THEN
8999              FORMAT(/,"Occupation matrix obtained from the green's function for atom: ",I3," l: ",I3," lp: ",I3)
                  WRITE(oUnit,8999) this%elem%atomType, this%elem%l, this%elem%lp
               ELSE
9000              FORMAT(/,"Occupation matrix obtained from the green's function for atom: ",I3," l: ",I3)
                  WRITE(oUnit,9000) this%elem%atomType, this%elem%l
               ENDIF
               WRITE(oUnit,"(A)") "In the |L,S> basis:"
               DO ispin = 1, MERGE(3, input%jspins, gfinp%l_mperp)
                  WRITE(oUnit,'(A,I0)') "Spin: ", ispin
                  WRITE(oUnit,l_form) ((occmtx(m,mp,ispin),m=-this%elem%l, this%elem%l),mp=-this%elem%lp, this%elem%lp)
               ENDDO

               IF(this%elem%l.EQ.this%elem%lp) THEN
                  nup = 0.0
                  DO m = -this%elem%l, this%elem%l
                     nup = nup + REAL(occmtx(m,m,1))
                  ENDDO
                  WRITE(oUnit,'(/,1x,A,I0,A,A,A,f8.4)') "l--> ",this%elem%l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Up trace: ", nup

                  IF(input%jspins.EQ.2) THEN
                     ndwn = 0.0
                     DO m = -this%elem%l, this%elem%l
                        ndwn = ndwn + REAL(occmtx(m,m,2))
                     ENDDO
                     WRITE(oUnit,'(1x,A,I0,A,A,A,f8.4)') "l--> ",this%elem%l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Down trace: ", ndwn
                  ENDIF

                  IF(gfinp%l_mperp) THEN
                     offd = cmplx_0
                     DO m = -this%elem%l, this%elem%l
                        offd = offd + occmtx(m,m,3)
                     ENDDO
                     WRITE(oUnit,'(1x,A,I0,A,A,A,f8.4)') "l--> ",this%elem%l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Offd trace (x): ", REAL(offd)
                     WRITE(oUnit,'(1x,A,I0,A,A,A,f8.4)') "l--> ",this%elem%l, " Contour(",TRIM(ADJUSTL(contourInp%label)),")    Spin-Offd trace (y): ", AIMAG(offd)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

      END FUNCTION occmtx_greensf_complete

      FUNCTION occmtx_greensf_spin(this,spin,gfinp,input,atoms,noco,nococonv,check,occError) Result(occmtx)

         USE m_rotMMPmat
         USE m_types_mat

         !calculates the occupation of the greens function for a given spin
         !The Greens-function should already be prepared on a energy contour ending at e_fermi
         !The occupation is calculated with:
         !
         ! n^sigma_mm' = -1/2pi int^Ef dz (G^+(z)^sigma_mm'-G^-(z)^sigma_mm')

         CLASS(t_greensf),                 INTENT(IN)    :: this
         INTEGER,                          INTENT(IN)    :: spin
         TYPE(t_gfinp),                    INTENT(IN)    :: gfinp
         TYPE(t_input),                    INTENT(IN)    :: input
         TYPE(t_atoms),                    INTENT(IN)    :: atoms
         TYPE(t_noco),                     INTENT(IN)    :: noco
         TYPE(t_nococonv),                 INTENT(IN)    :: nococonv
         LOGICAL,                 OPTIONAL,INTENT(IN)    :: check
         LOGICAL,                 OPTIONAL,INTENT(INOUT) :: occError

         COMPLEX :: occmtx(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

         INTEGER :: ind1,ind2,ipm,iz
         INTEGER :: l,lp,atomType,atomTypep,m,mp
         REAL    :: tr
         COMPLEX :: weight
         TYPE(t_mat) :: gmat
         CHARACTER(len=300) :: message
         TYPE(t_contourInp) :: contourInp

         contourInp = gfinp%contour(this%elem%iContour)
         IF(PRESENT(occError)) occError = .FALSE.

         l = this%elem%l
         lp = this%elem%lp
         atomType = this%elem%atomType
         atomTypep = this%elem%atomTypep

         !Check for Contours not reproducing occupations
         IF(contourInp%shape.EQ.CONTOUR_SEMICIRCLE_CONST.AND.ABS(contourInp%et).GT.1e-12) &
            WRITE(oUnit,*) "Energy contour not ending at efermi: These are not the actual occupations"
         IF(contourInp%shape.EQ.CONTOUR_DOS_CONST.AND..NOT.contourInp%l_dosfermi) &
            WRITE(oUnit,*) "Energy contour not weighted for occupations: These are not the actual occupations"

         occmtx = cmplx_0

         DO ipm = 1, 2
            !Integrate over the contour:
            DO iz = 1, this%contour%nz
               !get the corresponding gf-matrix
               weight = MERGE(this%contour%de(iz),conjg(this%contour%de(iz)),ipm.EQ.1)
               CALL this%get(atoms,iz,ipm.EQ.2,spin,gmat)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     occmtx(m,mp) = occmtx(m,mp) + ImagUnit/tpi_const * (-1)**(ipm-1) * gmat%data_c(ind1,ind2) &
                                                  * weight
                  ENDDO
               ENDDO
            ENDDO
            !For the contour 3 (real Axis just shifted with sigma) we can add the tails on both ends
            IF(contourInp%shape.EQ.CONTOUR_DOS_CONST.AND.contourInp%l_anacont) THEN
               !left tail
               weight = MERGE(this%contour%de(1),conjg(this%contour%de(1)),ipm.EQ.1)
               CALL this%get(atoms,1,ipm.EQ.2,spin,gmat)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     occmtx(m,mp) = occmtx(m,mp) - 1/tpi_const * gmat%data_c(ind1,ind2) * weight
                  ENDDO
               ENDDO
               !right tail
               weight = MERGE(this%contour%de(this%contour%nz),conjg(this%contour%de(this%contour%nz)),ipm.EQ.1)
               CALL this%get(atoms,this%contour%nz,ipm.EQ.2,spin,gmat)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0
                  DO mp = -lp,lp
                     ind2 = ind2 + 1
                     occmtx(m,mp) = occmtx(m,mp) + 1/tpi_const * gmat%data_c(ind1,ind2) * weight
                  ENDDO
               ENDDO
            ENDIF
         ENDDO

         !Rotate the occupation matrix into the global frame in real-space
         IF(noco%l_noco) THEN
            occmtx = rotMMPmat(occmtx,nococonv%alph(atomType),nococonv%beta(atomType),0.0,l)
         ELSE IF(noco%l_soc) THEN
            occmtx = rotMMPmat(occmtx,nococonv%phi,nococonv%theta,0.0,l)
         ENDIF

         !Sanity check are the occupations reasonable?
         IF(PRESENT(check)) THEN
            IF(check) THEN
               IF(spin<=input%jspins) THEN !Only the spin-diagonal parts
                  tr = 0.0
                  DO m = -l,l
                     tr = tr + REAL(occmtx(m,m))/(3.0-input%jspins)
                     IF(REAL(occmtx(m,m))/(3.0-input%jspins).GT. 1.05 .OR.&
                        REAL(occmtx(m,m))/(3.0-input%jspins).LT.-0.01) THEN

                        IF(PRESENT(occError)) THEN
                           occError = .TRUE.
                        ELSE
                           WRITE(message,9100) spin,m,REAL(occmtx(m,m))
9100                       FORMAT("Invalid element in mmpmat (spin ",I1,",m ",I2"): ",f14.8)
                           CALL juDFT_warn(TRIM(ADJUSTL(message)),calledby="occupationMatrix")
                        ENDIF
                     ENDIF
                  ENDDO
                  IF(tr.LT.-0.01.OR.tr.GT.2*l+1.1) THEN
                     IF(PRESENT(occError)) THEN
                        occError = .TRUE.
                     ELSE
                        WRITE(message,9110) spin,tr
9110                    FORMAT("Invalid occupation for spin ",I1,": ",f14.8)
                        CALL juDFT_warn(TRIM(ADJUSTL(message)),calledby="occupationMatrix")
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

      END FUNCTION occmtx_greensf_spin

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

      SUBROUTINE get_gf(this,atoms,iz,l_conjg,spin,gmat)

         USE m_types_mat

         !Returns the matrix belonging to energy point iz with l,lp,nType,nTypep
         !can also return the spherically averaged GF with the given scalar products

         CLASS(t_greensf),                   INTENT(IN)     :: this
         TYPE(t_atoms),                      INTENT(IN)     :: atoms
         INTEGER,                            INTENT(IN)     :: iz
         LOGICAL,                            INTENT(IN)     :: l_conjg
         INTEGER,                            INTENT(IN)     :: spin
         TYPE(t_mat),                        INTENT(INOUT)  :: gmat !Return matrix

         INTEGER matsize1,matsize2,ind1,ind2
         INTEGER m,mp,spin1,spin2,ipm,spin_start,spin_end,spin_ind,m_ind,mp_ind
         INTEGER l,lp,atomType,atomTypep,nspins,ilo,ilop,iLO_ind,iLOp_ind

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

         IF(spin.GT.3 .OR. spin.LT.1) THEN
            CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
         ENDIF

         matsize1 = 2*l+1
         matsize2 = 2*lp+1

         IF(.NOT.ALLOCATED(gmat%data_c)) THEN
            CALL gmat%init(.FALSE.,matsize1,matsize2)
         ELSE IF(matsize1.NE.gmat%matsize1.OR.matsize2.NE.gmat%matsize2) THEN
            CALL juDFT_error("Mismatch in matsizes", calledby="get_gf")
         ENDIF

         ipm = MERGE(2,1,l_conjg)

         gmat%data_c = cmplx_0

         IF(spin < 3) THEN
            spin1 = spin
            spin2 = spin
         ELSE
            spin1 = 2
            spin2 = 1
         ENDIF

         !Find the correct spin index in gmmpMat arrays
         spin_ind = MERGE(1,spin,nspins.EQ.1)

         IF(.NOT.this%l_sphavg) THEN
            IF(.NOT.ALLOCATED(this%scalarProducts%uun)) THEN
               CALL juDFT_error('Scalar products not available')
            ENDIF
            uun = this%scalarProducts%uun(spin1,spin2)
            dun = this%scalarProducts%dun(spin1,spin2)
            udn = this%scalarProducts%udn(spin1,spin2)
            ddn = this%scalarProducts%ddn(spin1,spin2)

            IF(ALLOCATED(this%uulo)) THEN
               uulon(:) = this%scalarProducts%uulon(:,spin1,spin2)
               uloun(:) = this%scalarProducts%uloun(:,spin1,spin2)
               dulon(:) = this%scalarProducts%dulon(:,spin1,spin2)
               ulodn(:) = this%scalarProducts%ulodn(:,spin1,spin2)

               uloulopn(:,:) = this%scalarProducts%uloulopn(:,:,spin1,spin2)
            ENDIF
         ENDIF

         ind1 = 0
         DO m = -l,l
            ind1 = ind1 + 1
            ind2 = 0
            DO mp = -lp,lp
               ind2 = ind2 + 1

               !-------------------------------------------------------------------
               ! Check wether we need to do some operation on the indices m and mp
               !-------------------------------------------------------------------
               IF(spin.EQ.2 .AND. nspins.EQ.1) THEN
                  !For a non-spin-polarized calculation we might still want the full
                  !matrix. Then we need to reverse the order (SOC prop m*s_z)
                  m_ind  = -m
                  mp_ind = -mp
               ELSE
                  !Do nothing
                  m_ind  = m
                  mp_ind = mp
               ENDIF
               !-------------------
               ! Fetch the values
               !-------------------
               IF(.NOT.this%l_sphavg) THEN
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

            ENDDO!mp
         ENDDO!m

      END SUBROUTINE get_gf

      SUBROUTINE getFullMatrix_gf(this,atoms,iz,l_conjg,gmat)

         USE m_types_mat

         !Return the full matrix with all spin blocks for the given energy point

         CLASS(t_greensf),                   INTENT(IN)     :: this
         TYPE(t_atoms),                      INTENT(IN)     :: atoms
         INTEGER,                            INTENT(IN)     :: iz
         LOGICAL,                            INTENT(IN)     :: l_conjg
         TYPE(t_mat),                        INTENT(INOUT)  :: gmat !Return matrix

         INTEGER :: matsize1, matsize2, nspins, ispin

         TYPE(t_mat) :: gmat_spin

         matsize1  = 2*this%elem%l+1
         matsize2 = 2*this%elem%l+1

         IF(.NOT.ALLOCATED(gmat%data_c)) THEN
            CALL gmat%init(.FALSE.,2*matsize1,2*matsize2)
         ELSE IF(2*matsize1.NE.gmat%matsize1.OR.2*matsize2.NE.gmat%matsize2) THEN
            CALL juDFT_error("Mismatch in matsizes", calledby="getFullMatrix_gf")
         ENDIF

         IF(ALLOCATED(this%gmmpMat)) THEN
            nspins = SIZE(this%gmmpMat,4)
         ELSE
            nspins = SIZE(this%uu,4)
         ENDIF

         DO ispin = 1, MAX(nspins,2)
            CALL this%get(atoms,iz,l_conjg,ispin,gmat_spin)

            IF(ispin<3) THEN
               gmat%data_c((ispin-1)*matsize1+1:ispin*matsize1,(ispin-1)*matsize2+1:ispin*matsize2) = gmat_spin%data_c
            ELSE IF(ispin.EQ.3) THEN
               gmat%data_c(1:matsize1,matsize2+1:2*matsize2) = gmat_spin%data_c
               gmat%data_c(matsize1+1:2*matsize1,1:matsize2) = conjg(transpose(gmat_spin%data_c))
            ELSE
               gmat%data_c(1:matsize1,matsize2+1:2*matsize2) = ImagUnit*conjg(gmat_spin%data_c)
               gmat%data_c(matsize1+1:2*matsize1,1:matsize2) = -ImagUnit*transpose(gmat_spin%data_c)
            ENDIF
         ENDDO

         IF(nspins==1) gmat%data_c = gmat%data_c * 0.5


      END SUBROUTINE getFullMatrix_gf

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

         IF(spin.GT.3 .OR. spin.LT.1) THEN
            CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
         ENDIF

         ipm = MERGE(2,1,l_conjg)

         IF(.NOT.ALLOCATED(gmat)) ALLOCATE(gmat(SIZE(f,1),this%contour%nz),source=cmplx_0)
         gmat = cmplx_0

         IF(spin < 3) THEN
            spin1 = spin
            spin2 = spin
         ELSE
            spin1 = 2
            spin2 = 1
         ENDIF
         !Find the correct spin index in gmmpMat arrays
         spin_ind = MERGE(1,spin,nspins.EQ.1)

         !-------------------------------------------------------------------
         ! Check wether we need to do some operation on the indices m and mp
         !-------------------------------------------------------------------
         IF(spin.EQ.2 .AND. nspins.EQ.1) THEN
            !For a non-spin-polarized calculation we might still want the full
            !matrix. Then we need to reverse the order (SOC prop m*s_z)
            m_ind  = -m
            mp_ind = -mp
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
                  iLO_ind = iLO_ind + 1
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

         IF(spin.GT.3 .OR. spin.LT.1) THEN
            CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
         ENDIF

         ipm = MERGE(2,1,l_conjg)

         IF(.NOT.ALLOCATED(gmat)) ALLOCATE(gmat(SIZE(f,1),SIZE(f,1)),source=cmplx_0)

         IF(spin < 3) THEN
            spin1 = spin
            spin2 = spin
         ELSE
            spin1 = 2
            spin2 = 1
         ENDIF
         !Find the correct spin index in gmmpMat arrays
         spin_ind = MERGE(1,spin,nspins.EQ.1)

         !-------------------------------------------------------------------
         ! Check wether we need to do some operation on the indices m and mp
         !-------------------------------------------------------------------
         IF(spin.EQ.2 .AND. nspins.EQ.1) THEN
            !For a non-spin-polarized calculation we might still want the full
            !matrix. Then we need to reverse the order (SOC prop m*s_z)
            m_ind  = -m
            mp_ind = -mp
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
                  iLO_ind = iLO_ind + 1
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

      END SUBROUTINE getRadialRadial_gf

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
            IF(spin.GT.3 .OR. spin.LT.1) THEN
               CALL juDFT_error("Invalid argument for spin",calledby="set_gf")
            ENDIF
         ENDIF

         l_full = .NOT.PRESENT(spin)
         IF(l_full.AND.SIZE(this%gmmpMat,4)>=3) CALL juDFT_error("Not implemented", calledby="set_gf")
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

      SUBROUTINE set_gfdata(this,repr_gf)

         !Sets all arrays equal to the data from another greensfunction

         CLASS(t_greensf),    INTENT(INOUT)  :: this
         TYPE(t_greensf),     INTENT(IN)     :: repr_gf

         !TODO check array sizes before

         IF(ALLOCATED(repr_gf%gmmpMat)) this%gmmpMat = repr_gf%gmmpMat
         IF(ALLOCATED(repr_gf%gmmpMat_k)) this%gmmpMat_k = repr_gf%gmmpMat_k
         IF(ALLOCATED(repr_gf%uu)) this%uu = repr_gf%uu
         IF(ALLOCATED(repr_gf%ud)) this%ud = repr_gf%ud
         IF(ALLOCATED(repr_gf%du)) this%du = repr_gf%du
         IF(ALLOCATED(repr_gf%dd)) this%dd = repr_gf%dd
         IF(ALLOCATED(repr_gf%uulo)) this%uulo = repr_gf%uulo
         IF(ALLOCATED(repr_gf%ulou)) this%ulou = repr_gf%ulou
         IF(ALLOCATED(repr_gf%ulod)) this%ulod = repr_gf%ulod
         IF(ALLOCATED(repr_gf%dulo)) this%dulo = repr_gf%dulo
         IF(ALLOCATED(repr_gf%uloulop)) this%uloulop = repr_gf%uloulop

      END SUBROUTINE set_gfdata

      SUBROUTINE rotate_gf(this,sym,atoms)

         USE m_rotMMPmat

         !Applies the given symmetry operation to the greens function
         CLASS(t_greensf),    INTENT(INOUT)  :: this
         TYPE(t_sym),         INTENT(IN)     :: sym
         TYPE(t_atoms),       INTENT(IN)     :: atoms

         INTEGER :: l,lp,iop,atomType,atomTypep,ikpt
         INTEGER :: ipm,ispin,iz,ilo,ilop,iLO_ind,iLOp_ind

         IF(this%elem%representative_elem<0) RETURN !Nothing to be done

         CALL timestart("Green's Function: Rotate (Symmetry)")
         l  = this%elem%l
         lp = this%elem%lp
         atomType = this%elem%atomType
         atomTypep = this%elem%atomTypep
         iop = this%elem%representative_op

         DO ipm = 1, 2
            DO iz = 1, this%contour%nz
               IF(ALLOCATED(this%gmmpMat)) THEN
                  this%gmmpMat(iz,:,:,:,ipm) = rotMMPmat(this%gmmpMat(iz,:,:,:,ipm),sym,iop,l,lp=lp)
               ELSE IF(ALLOCATED(this%uu)) THEN
                  this%uu(iz,:,:,:,ipm) = rotMMPmat(this%uu(iz,:,:,:,ipm),sym,iop,l,lp=lp)
                  this%dd(iz,:,:,:,ipm) = rotMMPmat(this%dd(iz,:,:,:,ipm),sym,iop,l,lp=lp)
                  this%ud(iz,:,:,:,ipm) = rotMMPmat(this%ud(iz,:,:,:,ipm),sym,iop,l,lp=lp)
                  this%du(iz,:,:,:,ipm) = rotMMPmat(this%du(iz,:,:,:,ipm),sym,iop,l,lp=lp)

                  IF(ALLOCATED(this%uulo)) THEN
                     iLO_ind = 0
                     DO ilo = 1, atoms%nlo(atomType)
                        IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                        iLO_ind = iLO_ind + 1
                        this%uulo(iz,:,:,iLO_ind,:,ipm) = rotMMPmat(this%uulo(iz,:,:,iLO_ind,:,ipm),sym,iop,l,lp=lp)
                        this%dulo(iz,:,:,iLO_ind,:,ipm) = rotMMPmat(this%dulo(iz,:,:,iLO_ind,:,ipm),sym,iop,l,lp=lp)
                     ENDDO
                     iLO_ind = 0
                     DO ilo = 1, atoms%nlo(atomTypep)
                        IF(atoms%llo(ilo,atomTypep).NE.lp) CYCLE
                        iLO_ind = iLO_ind + 1
                        this%ulou(iz,:,:,iLO_ind,:,ipm) = rotMMPmat(this%ulou(iz,:,:,iLO_ind,:,ipm),sym,iop,l,lp=lp)
                        this%ulod(iz,:,:,iLO_ind,:,ipm) = rotMMPmat(this%ulod(iz,:,:,iLO_ind,:,ipm),sym,iop,l,lp=lp)
                     ENDDO
                     iLO_ind = 0
                     DO ilo = 1, atoms%nlo(atomType)
                        IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                        iLOp_ind = 0
                        iLO_ind = iLO_ind + 1
                        DO ilop = 1, atoms%nlo(atomTypep)
                           IF(atoms%llo(ilop,atomType).NE.lp) CYCLE
                           iLOp_ind = iLOp_ind + 1
                           this%uloulop(iz,:,:,iLO_ind,iLOp_ind,:,ipm) = rotMMPmat(this%uloulop(iz,:,:,iLO_ind,iLOp_ind,:,ipm),sym,iop,l,lp=lp)
                        ENDDO
                     ENDDO
                  ENDIF
               ELSE IF(ALLOCATED(this%gmmpMat_k)) THEN
                  DO ikpt = 1, SIZE(this%gmmpMat_k,6)
                     this%gmmpMat_k(iz,:,:,:,ipm,ikpt) = rotMMPmat(this%gmmpMat_k(iz,:,:,:,ipm,ikpt),sym,iop,l,lp=lp)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CALL timestop("Green's Function: Rotate (Symmetry)")
      END SUBROUTINE rotate_gf

      SUBROUTINE rotate_euler_angles_gf(this,atoms,alpha,beta,gamma,spin_rotation,real_space_rotation)

         USE m_rotMMPmat

         !Applies the given symmetry operation to the greens function
         CLASS(t_greensf),    INTENT(INOUT)  :: this
         TYPE(t_atoms),       INTENT(IN)     :: atoms
         REAL,                INTENT(IN)     :: alpha
         REAL,                INTENT(IN)     :: beta
         REAL,                INTENT(IN)     :: gamma
         LOGICAL, OPTIONAL,   INTENT(IN)     :: spin_rotation
         LOGICAL, OPTIONAL,   INTENT(IN)     :: real_space_rotation

         INTEGER :: l,lp,atomType,atomTypep,ikpt
         INTEGER :: ipm,ispin,iz,ilo,ilop,iLO_ind,iLOp_ind

         IF(ABS(alpha).LT.1e-12.AND.ABS(beta).LT.1e-12.AND.ABS(gamma).LT.1e-12) RETURN !Nothing to be done

         CALL timestart("Green's Function: Rotate (Angles)")
         l  = this%elem%l
         lp = this%elem%lp
         atomType = this%elem%atomType
         atomTypep = this%elem%atomTypep

         DO ipm = 1, 2
            DO iz = 1, this%contour%nz
               IF(ALLOCATED(this%gmmpMat)) THEN
                  this%gmmpMat(iz,:,:,:,ipm) = rotMMPmat(this%gmmpMat(iz,:,:,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                         spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
               ELSE IF(ALLOCATED(this%uu)) THEN
                  this%uu(iz,:,:,:,ipm) = rotMMPmat(this%uu(iz,:,:,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                    spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                  this%ud(iz,:,:,:,ipm) = rotMMPmat(this%ud(iz,:,:,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                    spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                  this%du(iz,:,:,:,ipm) = rotMMPmat(this%du(iz,:,:,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                    spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                  this%dd(iz,:,:,:,ipm) = rotMMPmat(this%dd(iz,:,:,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                    spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)

                  IF(ALLOCATED(this%uulo)) THEN
                     iLO_ind = 0
                     DO ilo = 1, atoms%nlo(atomType)
                        IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                        iLO_ind = iLO_ind + 1
                        this%uulo(iz,:,:,iLO_ind,:,ipm) = rotMMPmat(this%uulo(iz,:,:,iLO_ind,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                                    spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                        this%dulo(iz,:,:,iLO_ind,:,ipm) = rotMMPmat(this%dulo(iz,:,:,iLO_ind,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                                    spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                     ENDDO
                     iLO_ind = 0
                     DO ilo = 1, atoms%nlo(atomTypep)
                        IF(atoms%llo(ilo,atomTypep).NE.lp) CYCLE
                        iLO_ind = iLO_ind + 1
                        this%ulou(iz,:,:,iLO_ind,:,ipm) = rotMMPmat(this%ulou(iz,:,:,iLO_ind,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                                    spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                        this%ulod(iz,:,:,iLO_ind,:,ipm) = rotMMPmat(this%ulod(iz,:,:,iLO_ind,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                                    spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                     ENDDO
                     iLO_ind = 0
                     DO ilo = 1, atoms%nlo(atomType)
                        IF(atoms%llo(ilo,atomType).NE.l) CYCLE
                        iLOp_ind = 0
                        iLO_ind = iLO_ind + 1
                        DO ilop = 1, atoms%nlo(atomTypep)
                           IF(atoms%llo(ilop,atomTypep).NE.lp) CYCLE
                           iLOp_ind = iLOp_ind + 1
                           this%uloulop(iz,:,:,iLO_ind,iLOp_ind,:,ipm) = rotMMPmat(this%uloulop(iz,:,:,iLO_ind,iLOp_ind,:,ipm),alpha,beta,gamma,l,lp=lp,&
                                                                                   spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                        ENDDO
                     ENDDO
                  ENDIF
               ELSE IF(ALLOCATED(this%gmmpMat_k)) THEN
                  DO ikpt = 1, SIZE(this%gmmpMat_k,6)
                     this%gmmpMat_k(iz,:,:,:,ipm,ikpt) = rotMMPmat(this%gmmpMat_k(iz,:,:,:,ipm,ikpt),alpha,beta,gamma,l,lp=lp,&
                                                                   spin_rotation=spin_rotation,real_space_rotation=real_space_rotation)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         CALL timestop("Green's Function: Rotate (Angles)")
      END SUBROUTINE rotate_euler_angles_gf

      PURE FUNCTION checkEmpty_greensf(this,m,mp,spin,ipm) Result(l_empty)

         CLASS(t_greensf),         INTENT(IN)   :: this
         INTEGER,                  INTENT(IN)   :: m
         INTEGER,                  INTENT(IN)   :: mp
         INTEGER,                  INTENT(IN)   :: spin
         INTEGER,                  INTENT(IN)   :: ipm

         LOGICAL :: l_empty

         IF(ALLOCATED(this%gmmpMat)) THEN
            l_empty = ALL(ABS(this%gmmpMat(:,m,mp,spin,ipm)).LT.1e-12)
         ELSE IF(ALLOCATED(this%uu)) THEN
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
         ELSE IF(ALLOCATED(this%gmmpMat_k)) THEN
            l_empty = ALL(ABS(this%gmmpMat_k(:,m,mp,spin,ipm,:)).LT.1e-12)
         ENDIF

      END FUNCTION checkEmpty_greensf

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
         IF(ALLOCATED(this%gmmpMat_k)) this%gmmpMat_k = cmplx_0

      END SUBROUTINE reset_gf


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
         IF(ALLOCATED(this%gmmpMat_k)) this%gmmpMat_k(:,m,mp,spin,ipm,:) = cmplx_0

      END SUBROUTINE resetSingleElem_gf

      FUNCTION integrateOverMT_greensf(this,atoms,input,gfinp,f,g,flo,l_fullRadial) Result(gIntegrated)

         USE m_intgr
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
         LOGICAL,                  OPTIONAL, INTENT(IN) :: l_fullRadial

         TYPE(t_greensf) :: gIntegrated

         LOGICAL :: l_fullRadialArg,l_explicit
         INTEGER :: l,lp,atomType,atomTypep,ipm,spin,m,mp,iz,jr,jrp
         REAL    :: realPart, imagPart
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
         !Do we have the offdiagonal scalar products
         l_explicit = .TRUE.
         IF(ALLOCATED(this%scalarProducts%uun)) l_explicit = .FALSE.

         !only intersite arguments have independent radial arguments ??
         l_fullRadialArg = l_fullRadialArg.AND.this%elem%isIntersite()

         DO ipm = 1, 2
            DO spin = 1 , SIZE(this%uu,4)
               IF(.NOT.l_explicit) THEN
                  DO iz = 1, this%contour%nz
                     CALL this%get(atoms,iz,ipm==2,spin,gmatTmp)
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
