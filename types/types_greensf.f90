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

         !for radial dependence
         COMPLEX, ALLOCATABLE :: uu(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: dd(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: du(:,:,:,:,:,:)
         COMPLEX, ALLOCATABLE :: ud(:,:,:,:,:,:)

         

         CONTAINS
            PROCEDURE, PASS :: init => greensf_init
            PROCEDURE       :: e_contour
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

         INTEGER i,j,r_dim,l_dim,spin_dim
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
         thisGREENSF%mode     = input%gf_mode
         IF(PRESENT(nz_in)) THEN
            thisGREENSF%nz = nz_in 
            thisGREENSF%nmatsub = matsub_in
         ELSE
            !Parameters for the energy contour in the complex plane

            IF(thisGREENSF%mode.EQ.1) THEN
                  thisGREENSF%nz = input%gf_n1+input%gf_n2+input%gf_n3+input%gf_nmatsub
                  thisGREENSF%nmatsub = input%gf_nmatsub
            ELSE IF(thisGREENSF%mode.EQ.2) THEN
               thisGREENSF%nz = input%gf_n
               thisGREENSF%nmatsub = 0
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

               IF(input%l_gfsphavg) THEN
                  spin_dim = MERGE(3,input%jspins,input%l_gfmperp)
                  ALLOCATE (thisGREENSF%gmmpMat(thisGREENSF%nz,MAX(1,atoms%n_gf),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,spin_dim,2))
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
            ENDIF
         ELSE
            !
            !intersite case; Here we look at l/=l' r/=r' and multiple sites
            r_dim = MAXVAL(atoms%jri(:)) 
            l_dim = lmax**2 + lmax

         ENDIF 

      END SUBROUTINE greensf_init

      SUBROUTINE e_contour(this,input,eb,et,ef)


         USE m_types_setup
         USE m_constants
         USE m_juDFT
         USE m_grule
         USE m_ExpSave

         IMPLICIT NONE

         CLASS(t_greensf),  INTENT(INOUT)  :: this
         TYPE(t_input),     INTENT(IN)     :: input
         REAL,              INTENT(IN)     :: eb  
         REAL,              INTENT(IN)     :: et
         REAL,              INTENT(IN)     :: ef
         
         INTEGER i, j, iz,nz

         REAL e1, e2, del, sigma
         COMPLEX de
         REAL r, xr, xm, c, s, a, b
         REAL x(this%nz), w(this%nz)



         IF(this%mode.EQ.1) THEN
            
            sigma = input%gf_sigma * pi_const

            IF(this%nmatsub.EQ.0) THEN
               !Equidistant contour (without vertical edges)
      
               de = (et-eb)/REAL(this%nz-1)
               DO iz = 1, this%nz
                  this%e(iz) = (iz-1) * del + eb + ImagUnit * sigma
                  this%de(iz) = de
               ENDDO

            ELSE IF(this%nmatsub > 0) THEN

               nz = 0
               !Left Vertical part (eb,0) -> (eb,sigma)
               de = this%nmatsub * CMPLX(0.0,sigma)
               CALL grule(input%gf_n1,x(1:(input%gf_n1)/2),w(1:(input%gf_n1)/2))
               x = -x
               DO i = 1, (input%gf_n1+3)/2-1
                  x(input%gf_n1-i+1) = -x(i)
                  w(input%gf_n1-i+1) =  w(i)
               ENDDO
               DO iz = 1, input%gf_n1
                  nz = nz + 1
                  IF(nz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
                  this%e(nz) = eb + de + de * x(iz) 
                  this%de(nz) = w(iz)*de
               ENDDO

               !Horizontal Part (eb,sigma) -> (et,sigma)
               de = (ef-30*input%gf_sigma-eb)/2.0
               CALL grule(input%gf_n2,x(1:(input%gf_n2)/2),w(1:(input%gf_n2)/2))
               x = -x
               DO i = 1, (input%gf_n2+3)/2-1
                  x(input%gf_n2-i+1) = -x(i)
                  w(input%gf_n2-i+1) =  w(i)
               ENDDO
               DO iz = 1, input%gf_n2
                  nz = nz + 1
                  IF(nz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
                  this%e(nz) = de*x(iz) + de + eb + 2 * this%nmatsub * ImagUnit * sigma
                  this%de(nz) = de*w(iz)
               ENDDO

               !Right Vertical part (et,sigma) -> infty
               CALL grule(input%gf_n3,x(1:(input%gf_n3)/2),w(1:(input%gf_n3)/2))
               x = -x
               DO i = 1, (input%gf_n3+3)/2-1
                  x(input%gf_n3-i+1) = -x(i)
                  w(input%gf_n3-i+1) =  w(i)
               ENDDO
               de = 30*input%tkb
               DO iz = 1, input%gf_n3
                  nz = nz + 1
                  IF(nz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
                  this%e(nz) = de*x(iz)+ef +  2 * this%nmatsub * ImagUnit * sigma
                  this%de(nz) = w(iz)*de/(1.0+exp_save((REAL(this%e(nz))-ef)/input%gf_sigma))
               ENDDO

               !Matsubara frequencies
               DO iz = this%nmatsub , 1, -1 
                  nz = nz + 1 
                  IF(nz.GT.this%nz) CALL juDFT_error("Dimension error in energy mesh",calledby="init_e_contour")
                  this%e(nz)  = ef + (2*iz-1) * ImagUnit *sigma 
                  this%de(nz) =  -2 *ImagUnit * sigma
               ENDDO 
               WRITE(*,1000) this%nz, this%nmatsub,input%gf_n1,input%gf_n2,input%gf_n3
            ELSE
               !DOES NOTHING ATM
            ENDIF
1000     FORMAT("Energy Contour for Green's Function Integration with nz: ", I5.1,"; nmatsub: ", I5.1,"; n1: ", I5.1,"; n2: ", I5.1,"; n3: ", I5.1)
         ELSE IF(this%mode.EQ.2) THEN

            !Semicircle
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

      END SUBROUTINE e_contour


END MODULE m_types_greensf