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

   IMPLICIT NONE

   PRIVATE

      TYPE t_greensf

         LOGICAL  :: l_onsite !This switch determines wether we look at a intersite or an onsite gf

         !Energy contour parameters
         INTEGER  :: mode  !Determines the shape of the contour (more information in kkintgr.f90)
         INTEGER  :: nz    !number of points in the contour
         INTEGER  :: nmatsub

         !array for energy contour
         COMPLEX, ALLOCATABLE  :: e(:)  !energy points
         COMPLEX, ALLOCATABLE  :: de(:) !weights for integration

         !Arrays for Green's function
         COMPLEX, ALLOCATABLE :: gmmpMat(:,:,:,:,:,:) 
         !Off-diagonal elements for noco calculations
         COMPLEX, ALLOCATABLE :: gmmpMat21(:,:,:,:,:)
         !Arrays for intersite Greens-functions argument order (E,n,n',L,L',spin) n is the site index
         !We store the radial function and the coefficients for those to obtain the imaginary part as the torage demands 
         !to store the whole green's function is too big
         !Because of that we need information on the energy grid in this case
         !Energy grid for Imaginary part
         INTEGER  :: ne       !number of energy grid points for imaginary part calculations
         REAL     :: e_top    !Cutoff energies
         REAL     :: e_bot
         REAL     :: del
         REAL     :: sigma       !Smoothing parameter(not used at the moment)


         !for radial dependence
         COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: du(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:,:)

         

         CONTAINS
            PROCEDURE, PASS :: init => greensf_init
            PROCEDURE       :: init_e_contour
      END TYPE t_greensf


   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE greensf_init(thisGREENSF,input,lmax,atoms,l_onsite,noco,nz_in,e_in,de_in,matsub_in)

         USE m_juDFT
         USE m_types_setup
         USE m_constants, only : lmaxU_const

         CLASS(t_greensf),    INTENT(INOUT)  :: thisGREENSF
         TYPE(t_atoms),       INTENT(IN)     :: atoms
         TYPE(t_input),       INTENT(IN)     :: input
         INTEGER,             INTENT(IN)     :: lmax
         TYPE(t_noco),        INTENT(IN)     :: noco
         LOGICAL,             INTENT(IN)     :: l_onsite
         !Pass a already calculated energy contour to the type (not used)
         INTEGER, OPTIONAL,   INTENT(IN)     :: nz_in
         INTEGER, OPTIONAL,   INTENT(IN)     :: matsub_in
         COMPLEX, OPTIONAL,   INTENT(IN)     :: e_in(:)
         COMPLEX, OPTIONAL,   INTENT(IN)     :: de_in(:)

         INTEGER i,j,r_dim,l_dim
         REAL    tol,n
         LOGICAL l_new

         thisGREENSF%l_onsite = l_onsite

         IF(.NOT.l_onsite.AND.noco%l_mperp) CALL juDFT_error("NOCO + intersite gf not implemented",calledby="greensf_init")
         IF(l_onsite.AND.noco%l_mperp) CALL juDFT_error("NOCO + onsite gf not implemented",calledby="greensf_init")

         !
         !Set up general parameters for the Green's function (intersite and onsite)
         !
         !
         !Setting up parameters for the energy contour
         !
         thisGREENSF%mode     = input%onsite_mode
         IF(PRESENT(nz_in)) THEN
            thisGREENSF%nz = nz_in 
            thisGREENSF%nmatsub = matsub_in
         ELSE
            !Parameters for the energy contour in the complex plane

            IF(thisGREENSF%mode.EQ.1) THEN
                  thisGREENSF%nz = input%onsite_nz+input%onsite_nmatsub
                  thisGREENSF%nmatsub = input%onsite_nmatsub
            ELSE IF(thisGREENSF%mode.EQ.2) THEN
               thisGREENSF%nz = input%onsite_nz+input%onsite_nmatsub
               thisGREENSF%nmatsub = input%onsite_nmatsub
            ENDIF
         END IF

         ALLOCATE (thisGREENSF%e(thisGREENSF%nz))
         ALLOCATE (thisGREENSF%de(thisGREENSF%nz))

         IF(PRESENT(e_in)) THEN
            thisGREENSF%e(:) = e_in(:)
            thisGREENSF%de(:)= de_in(:)
         ELSE
            !If no energy contour is given it is set up to zero
            thisGREENSF%e(:) = CMPLX(0.0,0.0)
            thisGREENSF%de(:)= CMPLX(0.0,0.0)
         END IF


         IF(thisGREENSF%l_onsite) THEN
            IF(atoms%n_gf.GT.0) THEN !Are there Green's functions to be calculated?

               IF(input%onsite_sphavg) THEN
                  ALLOCATE (thisGREENSF%gmmpMat(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins,2))
                  thisGREENSF%gmmpMat = 0.0
               ELSE
                  ALLOCATE (thisGREENSF%uu(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins,2))
                  ALLOCATE (thisGREENSF%dd(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins,2))
                  ALLOCATE (thisGREENSF%du(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins,2))
                  ALLOCATE (thisGREENSF%ud(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins,2))

                  thisGREENSF%uu = 0.0
                  thisGREENSF%dd = 0.0
                  thisGREENSF%du = 0.0
                  thisGREENSF%ud = 0.0
               ENDIF


               !Allocate arrays for non-colinear part
               IF(noco%l_mperp) THEN
                  ALLOCATE (thisGREENSF%gmmpMat21(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
                  thisGREENSF%gmmpMat21 = CMPLX(0.0,0.0)
               ENDIF
            ENDIF
         ELSE
            !
            !intersite case; Here we look at l/=l' r/=r' and multiple sites
            !
            !We cannot store the green's function for every radial point r,r' because it takes up way too much space
            r_dim = MAXVAL(atoms%jri(:)) 
            l_dim = lmax**2 + lmax

         ENDIF 

      END SUBROUTINE greensf_init

      SUBROUTINE init_e_contour(this,eb,et,ef,sigma,n1,n2,n3,nmatsub,beta)

         ! calculates the energy contour where the greens function is calculated
         ! mode determines the kind of contour between e_bot and the fermi energy
         ! mode = 1 gives a equidistant contour with imaginary part g%sigma with g%nz points

         ! mode = 2 gives a half circle with 2**g%nz points

         USE m_constants
         USE m_juDFT
         USE m_grule
         IMPLICIT NONE

         CLASS(t_greensf),  INTENT(INOUT)  :: this
         REAL,              INTENT(IN)     :: eb  
         REAL,              INTENT(IN)     :: et
         REAL,              INTENT(IN)     :: ef
         REAL, OPTIONAL,    INTENT(IN)     :: sigma
         INTEGER, OPTIONAL, INTENT(IN)     :: n1
         INTEGER, OPTIONAL, INTENT(IN)     :: n2
         INTEGER, OPTIONAL, INTENT(IN)     :: n3
         INTEGER, OPTIONAL, INTENT(IN)     :: nmatsub
         REAL, OPTIONAL,    INTENT(IN)     :: beta

         INTEGER i, j, iz, np,imatsub

         REAL e1, e2, del
         REAL r, xr, xm, c, s, a, b
         REAL x(this%nz), w(this%nz)



         IF(this%mode.EQ.1) THEN

            !Using a rectangular contour ending at efermi and including N_matsub matsubara frequencies 
            !at the moment we use a equidistant mesh to make interfacing with the hubbard 1 solver easier

            e1 = -1.0+ef
            e2 = 1.0+ef

            !Determine the imaginary part of the rectangle

            this%sigma = 2 * pi_const  * 1./beta/hartree_to_ev_const
            this%nmatsub = nmatsub

            del = (e2-e1)/REAL(n2-1)

            !ATM we only use an equidistant mesh at the imaginary part 2*pi*n_matsub*kT and 
            !the matsubara frequencies as thats the only thing we can comfortably get from the solver

            DO iz = 1, n2
               IF(iz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
               this%e(iz) = (iz-1) * del + e1 + ImagUnit * this%sigma
               this%de(iz) = del
            ENDDO
            iz = n2
            DO imatsub = nmatsub , 1 , -1
               iz = iz + 1 
               IF(iz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
               this%e(iz)  = e2 + ImagUnit * ((2*imatsub-1)*pi_const*1./beta)
               this%de(iz) =  ImagUnit * 2.*pi_const*1./beta
            ENDDO 


         ELSE IF(this%mode.EQ.2) THEN

            !In this mode we use a ellipsoid form for our energy contour with 2**n_in points
            !Further we use four-point gaussian quadrature
            !The method is based on an old kkr version 

            e1 = eb
            e2 = ef
            this%nmatsub = 0
            !Radius
            r  = (e2-e1)*0.5
            !midpoint
            xr = (e2+e1)*0.5

            CALL grule(this%nz,x(1:(this%nz)/2),w(1:(this%nz)/2))

            DO i = 1, this%nz/2
               x(this%nz-i+1) = -x(i)
               w(this%nz-i+1) =  w(i)
            ENDDO
            DO i = 1, this%nz
               this%e(i) = xr + ImagUnit * r * EXP(-ImagUnit*pi_const/2.0 * x(this%nz-i+1))
               this%de(i) = pi_const/2.0 * r * EXP(-ImagUnit*pi_const/2.0 * x(this%nz-i+1)) * w(this%nz-i+1)
            ENDDO

         ELSE

            CALL juDFT_error("Invalid mode for energy contour in Green's function calculation", calledby="init_e_contour")

         END IF



      END SUBROUTINE init_e_contour


END MODULE m_types_greensf