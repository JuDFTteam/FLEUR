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

   IMPLICIT NONE

   PRIVATE

   TYPE t_greensf

      !Energy contour parameters
      INTEGER  :: nz = 0   !number of points in the contour
      COMPLEX, ALLOCATABLE  :: e(:)  !energy points
      COMPLEX, ALLOCATABLE  :: de(:) !weights for integration

      !Arrays for Green's function
      COMPLEX, ALLOCATABLE :: gmmpMat(:,:,:,:,:,:)

      !for radial dependence
      COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: du(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:,:)

      CONTAINS
         PROCEDURE, PASS :: init    => greensf_init
         PROCEDURE       :: get     => get_gf
         PROCEDURE       :: set     => set_gf
         PROCEDURE       :: reset   => reset_gf
   END TYPE t_greensf

   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE greensf_init(thisGREENSF,gfinp,input,noco,nz_in,e_in,de_in)

         CLASS(t_greensf),    INTENT(INOUT)  :: thisGREENSF
         TYPE(t_gfinp),       INTENT(IN)     :: gfinp
         TYPE(t_input),       INTENT(IN)     :: input
         TYPE(t_noco),        INTENT(IN)     :: noco !Stays for now until everyting clear
         !Pass a already calculated energy contour to the type
         INTEGER, OPTIONAL,   INTENT(IN)     :: nz_in
         COMPLEX, OPTIONAL,   INTENT(IN)     :: e_in(:)
         COMPLEX, OPTIONAL,   INTENT(IN)     :: de_in(:)

         INTEGER spin_dim,lmax

         !
         !Setting up parameters for the energy contour if one was passed
         !
         IF(PRESENT(nz_in)) THEN
            thisGREENSF%nz = nz_in

            ALLOCATE(thisGREENSF%e(thisGREENSF%nz),source=cmplx_0)
            ALLOCATE(thisGREENSF%de(thisGREENSF%nz),source=cmplx_0)

            IF(PRESENT(e_in).AND.PRESENT(de_in)) THEN
               thisGREENSF%e(:) = e_in(:)
               thisGREENSF%de(:)= de_in(:)
            ELSE
               CALL juDFT_error("e_in and de_in arrays expected, but not found",&
                               hint="This is a bug in FLEUR, please report",&
                               calledby="greensf_init")
            ENDIF
         ELSE
            IF(gfinp%mode.EQ.1) THEN
               thisGREENSF%nz = gfinp%n1+gfinp%n2+gfinp%n3+gfinp%nmatsub
            ELSE IF(gfinp%mode.EQ.2) THEN
               thisGREENSF%nz = gfinp%ncirc
            ELSE IF(gfinp%mode.EQ.3) THEN
               thisGREENSF%nz = gfinp%nDOS
            ELSE
               CALL juDFT_error("No valid energy contour mode",calledby="greensf_init")
            ENDIF

            ALLOCATE(thisGREENSF%e(thisGREENSF%nz),source=cmplx_0)
            ALLOCATE(thisGREENSF%de(thisGREENSF%nz),source=cmplx_0)
         ENDIF

         spin_dim = MERGE(3,input%jspins,gfinp%l_mperp)
         lmax = lmaxU_const
         IF(gfinp%l_sphavg) THEN
            ALLOCATE(thisGREENSF%gmmpMat(thisGREENSF%nz,-lmax:lmax,-lmax:lmax,spin_dim,2,MAX(1,gfinp%n)),source=cmplx_0)
         ELSE
            ALLOCATE(thisGREENSF%uu(thisGREENSF%nz,-lmax:lmax,-lmax:lmax,spin_dim,2,MAX(1,gfinp%n)),source=cmplx_0)
            ALLOCATE(thisGREENSF%dd(thisGREENSF%nz,-lmax:lmax,-lmax:lmax,spin_dim,2,MAX(1,gfinp%n)),source=cmplx_0)
            ALLOCATE(thisGREENSF%du(thisGREENSF%nz,-lmax:lmax,-lmax:lmax,spin_dim,2,MAX(1,gfinp%n)),source=cmplx_0)
            ALLOCATE(thisGREENSF%ud(thisGREENSF%nz,-lmax:lmax,-lmax:lmax,spin_dim,2,MAX(1,gfinp%n)),source=cmplx_0)
         ENDIF

      END SUBROUTINE greensf_init

      SUBROUTINE get_gf(this,gmat,gfinp,input,iz,l,nType,l_conjg,spin,lp,nTypep,u,udot)

         USE m_types_mat

         !Returns the matrix belonging to energy point iz with l,lp,nType,nTypep
         !when jr (and jrp) are given return for that radial point

         CLASS(t_greensf),    INTENT(IN)  :: this
         TYPE(t_input),       INTENT(IN)  :: input
         TYPE(t_gfinp),       INTENT(IN)  :: gfinp
         TYPE(t_mat),         INTENT(OUT) :: gmat !Return matrix

         INTEGER,             INTENT(IN)  :: iz
         INTEGER,             INTENT(IN)  :: nType
         INTEGER,             INTENT(IN)  :: l
         LOGICAL,             INTENT(IN)  :: l_conjg
         INTEGER, OPTIONAL,   INTENT(IN)  :: spin
         INTEGER, OPTIONAL,   INTENT(IN)  :: nTypep
         INTEGER, OPTIONAL,   INTENT(IN)  :: lp
         REAL   , OPTIONAL,   INTENT(IN)  :: u(:,:)       !Radial functions at the point where you want to evaluate the greens function
         REAL   , OPTIONAL,   INTENT(IN)  :: udot(:,:)

         INTEGER matsize1,matsize2,i_gf,i,j,ind1,ind2,ind1_start,ind2_start
         INTEGER m,mp,spin1,spin2,ipm,ispin,ispin_end,spin_ind,m_ind,mp_ind
         INTEGER lp_loop
         LOGICAL l_radial,l_full

         IF(PRESENT(u).OR.PRESENT(udot).AND.gfinp%l_sphavg) THEN
            CALL juDFT_error("Greens function not calculated for radial dependence", calledby="get_gf")
         ENDIF

         IF((PRESENT(u).AND..NOT.PRESENT(udot)).OR.&
            (PRESENT(udot).AND..NOT.PRESENT(u))) THEN
            CALL juDFT_error("Not a valid input: Either provide both u and udot or neither of them", calledby="get_gf")
         ENDIF

         l_radial = PRESENT(u).AND.PRESENT(udot)

         IF(PRESENT(spin)) THEN
            IF(spin.GT.4.OR.spin.LT.1) THEN
               CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
            ENDIF
         END IF

         !Determine matsize for the result gmat (if spin is given only return this diagonal element)
         l_full = .NOT.PRESENT(spin)
         matsize1 = (2*l+1) * MERGE(2,1,l_full)

         IF(PRESENT(lp)) THEN
            matsize2 = (2*lp+1) * MERGE(2,1,l_full)
         ELSE
            matsize2 = matsize1
         ENDIF

         IF(.NOT.ALLOCATED(gmat%data_c)) THEN
            CALL gmat%init(.FALSE.,matsize1,matsize2)
         ELSE IF(matsize1.NE.gmat%matsize1.OR.matsize2.NE.gmat%matsize2) THEN
            CALL juDFT_error("Mismatch in matsizes", calledby="get_gf")
         ENDIF

         IF(.NOT.PRESENT(lp)) THEN
            lp_loop = l
         ELSE
            lp_loop = lp
         ENDIF

         !Find the index i_gf corresponding to l,lp,nType,nTypep
         i_gf = gfinp%find(l,nType,lp=lp,nTypep=nTypep)
         ipm = MERGE(2,1,l_conjg)

         gmat%data_c = cmplx_0
         ispin_end = MERGE(4,2,gfinp%l_mperp)

         DO ispin = MERGE(1,spin,l_full), MERGE(ispin_end,spin,l_full)
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
            spin_ind = MERGE(ispin,1,input%jspins.EQ.2)
            spin_ind = MERGE(3,spin_ind,ispin.EQ.4)
            !Find the right quadrant in gmat
            IF(l_full) THEN
               ind1_start = (spin1-1)*(2*l+1)
               ind2_start = (spin2-1)*(2*lp_loop+1)
            ELSE
               ind1_start = 0
               ind2_start = 0
            ENDIF

            ind1 = ind1_start
            DO m = -l,l
               ind1 = ind1 + 1
               ind2 = ind2_start
               DO mp = -lp_loop,lp_loop
                  ind2 = ind2 + 1

                  !-------------------------------------------------------------------
                  ! Check wether we need to do some operation on the indices m and mp
                  !-------------------------------------------------------------------
                  IF(ispin.EQ.2.AND.input%jspins.EQ.1) THEN
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
                  !-------------------_ind
                  IF(l_radial) THEN
                     gmat%data_c(ind1,ind2) = this%uu(iz,i_gf,m_ind,mp_ind,spin_ind,ipm) * u(1,spin1)    * u(2,spin2)     + &
                                              this%dd(iz,i_gf,m_ind,mp_ind,spin_ind,ipm) * udot(1,spin1) * udot(2,spin2)  + &
                                              this%du(iz,i_gf,m_ind,mp_ind,spin_ind,ipm) * udot(1,spin1) * u(2,spin2)     + &
                                              this%ud(iz,i_gf,m_ind,mp_ind,spin_ind,ipm) * u(1,spin1)    * udot(2,spin2)
                  ELSE
                     gmat%data_c(ind1,ind2) = this%gmmpMat(iz,m_ind,mp_ind,spin_ind,ipm,i_gf)
                  ENDIF
                  !------------------------
                  ! Additional operations
                  !------------------------
                  !Spin-degeneracy when using a full matrix and having input%jspins.EQ.1
                  IF(l_full) gmat%data_c(ind1,ind2) = gmat%data_c(ind1,ind2)/(3.0-input%jspins)
                  !Complex conjugate for spin 4
                  IF(ispin.EQ.4) gmat%data_c(ind1,ind2) = conjg(gmat%data_c(ind1,ind2))

               ENDDO!mp
            ENDDO!m
         ENDDO!ispin

      END SUBROUTINE get_gf

      SUBROUTINE set_gf(this,gmat,gfinp,input,iz,l,nType,l_conjg,spin,lp,nTypep)

         USE m_types_mat

         !Sets the spherically averaged greens function matrix belonging to energy point iz with l,lp,nType,nTypep
         !equal to gmat

         CLASS(t_greensf),    INTENT(INOUT)  :: this
         TYPE(t_gfinp),       INTENT(IN)     :: gfinp
         TYPE(t_input),       INTENT(IN)     :: input
         TYPE(t_mat),         INTENT(IN)     :: gmat

         INTEGER,             INTENT(IN)     :: iz
         INTEGER,             INTENT(IN)     :: nType
         INTEGER,             INTENT(IN)     :: l
         LOGICAL,             INTENT(IN)     :: l_conjg
         INTEGER, OPTIONAL,   INTENT(IN)     :: spin
         INTEGER, OPTIONAL,   INTENT(IN)     :: nTypep
         INTEGER, OPTIONAL,   INTENT(IN)     :: lp

         INTEGER matsize1,matsize2,i_gf,i,j,ind1,ind2,ind1_start,ind2_start
         INTEGER m,mp,spin1,spin2,ipm,ispin,ispin_end
         INTEGER lp_loop

         IF(PRESENT(spin)) THEN
            IF(spin.GT.4.OR.spin.LT.1) THEN
               CALL juDFT_error("Invalid argument for spin",calledby="get_gf")
            ENDIF
         ENDIF

         !Determine matsize for the result gmat (if spin is given only return this digonal element)
         matsize1 = (2*l+1) * MERGE(1,2,PRESENT(spin))
         IF(PRESENT(lp)) THEN
            matsize2 = (2*lp+1) * MERGE(1,2,PRESENT(spin))
         ELSE
            matsize2 = matsize1
         ENDIF

         !Check the expected matsizes against the actual
         IF(matsize1.NE.gmat%matsize1.OR.matsize2.NE.gmat%matsize2) THEN
            CALL juDFT_error("Mismatch in matsizes", calledby="set_gf")
         ENDIF

         IF(.NOT.PRESENT(lp)) THEN
            lp_loop = l
         ELSE
            lp_loop = lp
         ENDIF

         !Find the index i_gf corresponding to l,lp,nType,nTypep
         i_gf = gfinp%find(l,nType,lp=lp,nTypep=nTypep)
         ipm = MERGE(2,1,l_conjg)

         ispin_end = MERGE(3,input%jspins,gfinp%l_mperp)

         DO ispin = MERGE(spin,1,PRESENT(spin)), MERGE(spin,ispin_end,PRESENT(spin))
            !Find the right quadrant in gmat according to the spin index
            IF(ispin.EQ.2.AND.input%jspins.EQ.1) CYCLE
            IF(.NOT.PRESENT(spin)) THEN
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
               ind2_start = (spin2-1)*(2*lp_loop+1)
            ELSE
               ind1_start = 0
               ind2_start = 0
            ENDIF
            ind1 = ind1_start
            DO m = -l,l
               ind1 = ind1 + 1
               ind2 = ind2_start
               DO mp = -lp_loop,lp_loop
                  ind2 = ind2 + 1
                  this%gmmpMat(iz,m,mp,ispin,ipm,i_gf) = gmat%data_c(ind1,ind2)*MERGE(1.0,2.0/input%jspins,PRESENT(spin))
               ENDDO
            ENDDO
         ENDDO

      END SUBROUTINE set_gf

      SUBROUTINE reset_gf(this,gfinp,l,nType,lp,ntypep)

         !---------------------------------------------------
         ! Sets all gmmpMat arrays back to 0
         ! if no element is specified all elements are reset
         !---------------------------------------------------

         CLASS(t_greensf),       INTENT(INOUT)  :: this
         TYPE(t_gfinp),          INTENT(IN)     :: gfinp
         INTEGER, OPTIONAL,      INTENT(IN)     :: l, lp
         INTEGER, OPTIONAL,      INTENT(IN)     :: nType, nTypep

         INTEGER :: i_gf,i_reset
         LOGICAL :: lFullreset

         lFullreset = .TRUE.

         IF(PRESENT(l).OR.PRESENT(nType)) THEN
            IF(PRESENT(l).AND.PRESENT(nType)) THEN
               i_reset = gfinp%find(l,nType,lp=lp,nTypep=nTypep)
               lFullreset = .FALSE.
            ELSE
               CALL juDFT_error("Please specify both l and nType or none", &
                     hint="This is a bug in FLEUR please report", calledby="reset_gf")
            ENDIF
         ENDIF

         DO i_gf = MERGE(1,i_reset,lFullreset), MERGE(gfinp%n,i_reset,lFullreset)
            this%gmmpMat(:,:,:,:,:,i_gf) = cmplx_0
            IF(.NOT.gfinp%l_sphavg) THEN
               this%uu(:,:,:,:,:,i_gf) = cmplx_0
               this%ud(:,:,:,:,:,i_gf) = cmplx_0
               this%du(:,:,:,:,:,i_gf) = cmplx_0
               this%dd(:,:,:,:,:,i_gf) = cmplx_0
            ENDIF
         ENDDO

      END SUBROUTINE reset_gf

END MODULE m_types_greensf
