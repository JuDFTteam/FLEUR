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
         INTEGER, ALLOCATABLE :: nr(:) !dimension(atoms%n_gf) number of radial points
                                       !in case of spherical average nr(:) = 1

         !Energy contour parameters
         INTEGER  :: mode  !Determines the shape of the contour (more information in kkintgr.f90)
         INTEGER  :: nz    !number of points in the contour

         !array for energy contour
         COMPLEX, ALLOCATABLE  :: e(:)  !energy points
         COMPLEX, ALLOCATABLE  :: de(:) !weights for integration

         !Arrays for Green's function
         COMPLEX, ALLOCATABLE :: gmmpMat(:,:,:,:,:,:,:) 
         !Off-diagonal elements for noco calculations
         COMPLEX, ALLOCATABLE :: gmmpMat21(:,:,:,:,:,:)
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

         REAL, ALLOCATABLE :: uu_int(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: dd_int(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: du_int(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: ud_int(:,:,:,:,:,:)
         REAL, ALLOCATABLE :: R(:,:,:,:,:,:)

         

         CONTAINS
            PROCEDURE, PASS :: init => greensf_init
            PROCEDURE       :: init_e_contour
      END TYPE t_greensf


   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE greensf_init(thisGREENSF,input,lmax,atoms,kpts,noco,l_onsite,nz_in,e_in,de_in)

         USE m_juDFT
         USE m_types_setup
         USE m_types_kpts
         USE m_constants, only : lmaxU_const

         CLASS(t_greensf),       INTENT(INOUT)  :: thisGREENSF
         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_input),          INTENT(IN)     :: input
         INTEGER,                INTENT(IN)     :: lmax
         TYPE(t_kpts), OPTIONAL, INTENT(IN)     :: kpts
         TYPE(t_noco), OPTIONAL, INTENT(IN)     :: noco
         LOGICAL,                INTENT(IN)     :: l_onsite
         !Pass a already calculated energy contour to the type (not used)
         INTEGER, OPTIONAL,      INTENT(IN)     :: nz_in
         COMPLEX, OPTIONAL,      INTENT(IN)     :: e_in(:)
         COMPLEX, OPTIONAL,      INTENT(IN)     :: de_in(:)

         INTEGER i,j,r_dim,l_dim
         REAL    tol,n
         LOGICAL l_new

         thisGREENSF%l_onsite = l_onsite

         IF(.NOT.l_onsite.AND.noco%l_mperp) CALL juDFT_error("NOCO + intersite gf not implented",calledby="greensf_init")

         !
         !Set up general parameters for the Green's function (intersite and onsite)
         !
         !
         !Setting up parameters for the energy contour
         !
         IF(PRESENT(nz_in)) THEN
            thisGREENSF%nz = nz_in
         ELSE
            !Parameters for the energy contour in the complex plane
            thisGREENSF%mode     = input%onsite_mode

            IF(thisGREENSF%mode.EQ.1) THEN
               thisGREENSF%nz = input%onsite_nz
            ELSE IF(thisGREENSF%mode.EQ.2) THEN
               !We want a power of 2 as number of points
               n = LOG(REAL(input%onsite_nz))/LOG(2.0)
               IF(n.NE.AINT(n)) THEN
                  WRITE(*,*) "This mode for the energy contour uses 2^n number of points."
                  WRITE(*,*) "Setting nz = ", 2**AINT(n) 
               END IF
               thisGREENSF%nz = 2**AINT(n)
            END IF
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
               !Set number of radial points
               ALLOCATE(thisGREENSF%nr(MAX(1,atoms%n_gf)))

               IF(input%onsite_sphavg) THEN
                  thisGREENSF%nr(:) = 1
               ELSE
                  DO i = 1, atoms%n_gf
                     thisGREENSF%nr(i) = atoms%jri(atoms%onsiteGF(i)%atomType)
                  ENDDO
               END IF

               ALLOCATE (thisGREENSF%gmmpMat(MAXVAL(thisGREENSF%nr(:)),thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins,2))
               thisGREENSF%gmmpMat = CMPLX(0.0,0.0)

               !Allocate arrays for non-colinear part
               IF(noco%l_mperp) THEN
                  ALLOCATE (thisGREENSF%gmmpMat21(MAXVAL(thisGREENSF%nr(:)),thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
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

            ALLOCATE (thisGREENSF%uu_int(thisGREENSF%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))
            ALLOCATE (thisGREENSF%dd_int(thisGREENSF%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))
            ALLOCATE (thisGREENSF%du_int(thisGREENSF%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))
            ALLOCATE (thisGREENSF%ud_int(thisGREENSF%ne,atoms%nat,atoms%nat,0:l_dim,0:l_dim,input%jspins))

            ALLOCATE (thisGREENSF%R(thisGREENSF%ne,r_dim,2,atoms%nat,0:l_dim,input%jspins))

            thisGREENSF%uu_int      = 0.0
            thisGREENSF%dd_int      = 0.0
            thisGREENSF%du_int      = 0.0
            thisGREENSF%ud_int      = 0.0

            thisGREENSF%R           = 0.0

         ENDIF 

      END SUBROUTINE greensf_init

      SUBROUTINE init_e_contour(this,eb,ef,sigma)

         ! calculates the energy contour where the greens function is calculated
         ! mode determines the kind of contour between e_bot and the fermi energy
         ! mode = 1 gives a equidistant contour with imaginary part g%sigma with g%nz points

         ! mode = 2 gives a half circle with 2**g%nz points

         USE m_constants
         USE m_juDFT

         IMPLICIT NONE

         CLASS(t_greensf),  INTENT(INOUT)  :: this
         REAL,              INTENT(IN)     :: eb  
         REAL,              INTENT(IN)     :: ef
         REAL, OPTIONAL,    INTENT(IN)     :: sigma


         INTEGER i, j, iz, np

         REAL e1, e2, del
         REAL psi(4), wpsi(4), r, xr, xm, c, s, a, b



         IF(this%mode.EQ.1) THEN

            e1 = eb
            e2 = ef

            del = (e2-e1)/REAL(this%nz-1)

            DO i = 1, this%nz
               IF(PRESENT(sigma)) THEN
                  this%e(i) = (i-1)*del + e1 + ImagUnit * sigma
               ELSE
                  CALL juDFT_error("Sigma not given for energy contour",calledby="init_e_contour")
               ENDIF
            ENDDO

            this%de(:) = del

         ELSE IF(this%mode.EQ.2) THEN

            !In this mode we use a ellipsoid form for our energy contour with 2**n_in points
            !Further we use four-point gaussian quadrature
            !The method is based on an old kkr version 

            np = INT(this%nz/4.)

            e1 = eb
            e2 = ef

            !Radius
            r  = (e2-e1)*0.5
            !midpoint
            xr = (e2+e1)*0.5

            !supports for four-point gaussian quadrature
            a = 0.43056815579702629
            b = 0.16999052179242813

            psi(1) =    a/np
            psi(2) =    b/np
            psi(3) =   -b/np
            psi(4) =   -a/np

            !weights for four-point gaussian quadrature
            a = 0.17392742256872693
            b = 0.32607257743127307

            wpsi(1) =   a/np
            wpsi(2) =   b/np
            wpsi(3) =   b/np
            wpsi(4) =   a/np

            iz = 1

            DO i = 1, np

               !midpoint for the current interval in terms of angle
               xm = (np-i+0.5)/np

               DO j = 1, 4

                  !the squaring moves the points closer to the right end of the contour where the fermi energy is located

                  c = cos((psi(j)+xm)**2*pi_const)
                  s = sin((psi(j)+xm)**2*pi_const)

                  !TODO: implement sigma to ensure the integral can be calculated with finite sigma (look at weights)

                  this%e(iz) = CMPLX(xr+r*c, r*s*0.25)

                  this%de(iz) = pi_const * CMPLX((psi(j)+xm)*r*s*wpsi(j)*2.0,&
                                             -(psi(j)+xm)*r*c*wpsi(j)*0.5)

                  iz = iz+1

               ENDDO

            ENDDO

         ELSE

            CALL juDFT_error("Invalid mode for energy contour in Green's function calculation", calledby="init_e_contour")

         END IF



      END SUBROUTINE init_e_contour

END MODULE m_types_greensf