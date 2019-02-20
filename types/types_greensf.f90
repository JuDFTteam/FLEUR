
! This module defines a type for the Greens-functions used in the LDA+HIA formalism
! The Greens function is an on-Site Green' function which is stored in the matrix gmmpMat 
! and only contains Blocks with l = lprime in the MT-sphere

MODULE m_types_greensf

   IMPLICIT NONE

   PRIVATE

      TYPE t_greensf 

         !Energy grid for Imaginary part
         INTEGER  :: ne       !number of energy grid points for imaginary part calculations
         REAL     :: e_top    !Cutoff energies
         REAL     :: e_bot
         REAL     :: del


         REAL     :: sigma       !Smoothing parameter(not used at the moment)
     
         LOGICAL  :: l_tetra  !Determines wether to use the tetrahedron method for Brillouin-Zone integration
         LOGICAL  :: l_ef     !This switch determines wether the energy contour ends at efermi
         INTEGER, ALLOCATABLE :: nr(:) !dimension(atoms%n_hia) number of radial points
                                       !in case of spherical average nr(:) = 1

         !Energy contour parameters
         INTEGER  :: mode  !Determines the shape of the contour (more information in kkintgr.f90)
         INTEGER  :: nz    !number of points in the contour

         !array for energy contour
         COMPLEX, ALLOCATABLE  :: e(:)  !energy points
         COMPLEX, ALLOCATABLE  :: de(:) !weights for integration

         !Arrays for Green's function
         COMPLEX, ALLOCATABLE :: gmmpMat(:,:,:,:,:,:) 

         !These arrays are needed to calculate the onsite-Green's function from the eigenstates and eigenvalues of the DFT calculation
         REAL, ALLOCATABLE :: im_gmmpMat(:,:,:,:,:,:)

         ! These arrays are only used in the case we want the green's function with radial dependence
         REAL, ALLOCATABLE :: uu(:,:,:,:,:)
         REAL, ALLOCATABLE :: dd(:,:,:,:,:)
         REAL, ALLOCATABLE :: du(:,:,:,:,:)

         CONTAINS
            PROCEDURE, PASS :: init => greensf_init
            PROCEDURE       :: init_e_contour
            PROCEDURE       :: calc_mmpmat
      END TYPE t_greensf


   PUBLIC t_greensf

   CONTAINS

      SUBROUTINE greensf_init(thisGREENSF,input,atoms,kpts,dimension,l_onsite,nz_in,e_in,de_in)

         USE m_types_setup
         USE m_types_kpts
         USE m_constants, only : lmaxU_const

         CLASS(t_greensf),       INTENT(INOUT)  :: thisGREENSF
         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_input),          INTENT(IN)     :: input
         TYPE(t_kpts),           INTENT(IN)     :: kpts
         TYPE(t_dimension),      INTENT(IN)     :: dimension
         LOGICAL,                INTENT(IN)     :: l_onsite !this switch determines wether we want to use this type to calculate the on-site green's function
         !Pass a already calculated energy contour to the type
         INTEGER, OPTIONAL,      INTENT(IN)     :: nz_in
         COMPLEX, OPTIONAL,      INTENT(IN)     :: e_in(:)
         COMPLEX, OPTIONAL,      INTENT(IN)     :: de_in(:)

         INTEGER n, i_hia

         !Parameters for calculation of the imaginary part
         thisGREENSF%ne       = input%ldahia_ne
         thisGREENSF%e_top    = input%ldahia_etop
         thisGREENSF%e_bot    = input%ldahia_ebot
         thisGREENSF%sigma    = input%ldahia_sigma

         thisGREENSF%l_tetra  = input%ldahia_tetra
         thisGREENSF%l_ef     = .false.

         !set up energy grid for imaginary part
         thisGREENSF%del      = (thisGREENSF%e_top-thisGREENSF%e_bot)/REAL(thisGREENSF%ne-1)

         !Set number of radial points
         ALLOCATE(thisGREENSF%nr(atoms%n_hia))

         IF(input%ldahia_sphavg) THEN
            thisGREENSF%nr(:) = 1
         ELSE
            DO i_hia = 1, atoms%n_hia
               n = atoms%lda_hia(i_hia)%atomType  
   
               thisGREENSF%nr(i_hia) = atoms%jri(n)
            ENDDO
         END IF

         IF(PRESENT(nz_in)) THEN
            thisGREENSF%nz = nz_in
         ELSE
            !Parameters for the energy contour in the complex plane
            thisGREENSF%mode     = input%ldahia_mode

            IF(thisGREENSF%mode.EQ.1) THEN
               thisGREENSF%nz = input%ldahia_nin
            ELSE IF(thisGREENSF%mode.EQ.2) THEN
               n = input%ldahia_nin
               !ensure that we don't flood the memory accidentally
               IF(n.LT.2) n = 2
               IF(n.GT.7) n = 7
               thisGREENSF%nz = 2**n
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

         IF(l_onsite) THEN
            IF (.NOT.input%ldahia_sphavg) THEN 
               ALLOCATE (thisGREENSF%uu(thisGREENSF%ne,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
               ALLOCATE (thisGREENSF%dd(thisGREENSF%ne,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
               ALLOCATE (thisGREENSF%du(thisGREENSF%ne,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
              
               thisGREENSF%uu      = 0.0
               thisGREENSF%dd      = 0.0
               thisGREENSF%du      = 0.0

            END IF


            ALLOCATE (thisGREENSF%im_gmmpMat(MAXVAL(thisGREENSF%nr(:)),thisGREENSF%ne,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))
            thisGREENSF%im_gmmpMat  = 0.0
         END IF

         ALLOCATE (thisGREENSF%gmmpMat(MAXVAL(thisGREENSF%nr(:)),thisGREENSF%nz,atoms%n_hia,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins))

         thisGREENSF%gmmpMat = CMPLX(0.0,0.0)

      END SUBROUTINE greensf_init

      SUBROUTINE init_e_contour(this,ef)

         ! calculates the energy contour where the greens function is calculated
         ! mode determines the kind of contour between e_bot and the fermi energy (if l_ef = .true.) 
         ! mode = 1 gives a equidistant contour with imaginary part g%sigma with g%nz points

         ! mode = 2 gives a half circle with 2**g%nz points

         USE m_constants
         USE m_juDFT

         IMPLICIT NONE

         CLASS(t_greensf),  INTENT(INOUT)  :: this
         REAL,   OPTIONAL, INTENT(IN)     :: ef


         INTEGER i, j, iz, np

         REAL e1, e2, del
         REAL psi(4), wpsi(4), r, xr, xm, c, s, a, b



         IF(this%mode.EQ.1) THEN

            e1 = this%e_bot
            e2 = this%e_top

            IF(PRESENT(ef)) THEN
               e2 = ef
               this%l_ef = .true.
            ENDIF

            del = (e2-e1)/REAL(this%nz-1)

            DO i = 1, this%nz

               this%e(i) = (i-1)*del + e1 + ImagUnit * this%sigma

            ENDDO

            this%de(:) = del

         ELSE IF(this%mode.EQ.2) THEN

            !In this mode we use a ellipsoid form for our energy contour with 2**n_in points
            !Further we use four-point gaussian quadrature
            !The method is based on an old kkr version 

            np = INT(this%nz/4.)

            e1 = this%e_bot
            e2 = this%e_top

            IF(PRESENT(ef)) THEN
               e2 = ef
               this%l_ef = .true.
            ENDIF


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

      SUBROUTINE calc_mmpmat(this,atoms,sym,jspins,mmpMat)


         !calculates the occupation of a orbital treated with DFT+HIA from the related greens function

         !The Greens-function should already be prepared on a energy contour ending at e_fermi

         USE m_types_setup
         USE m_constants
         USE m_juDFT
         USE m_intgr

         IMPLICIT NONE

         CLASS(t_greensf),       INTENT(IN)     :: this
         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_sym),            INTENT(IN)     :: sym
         COMPLEX,                INTENT(INOUT)  :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,atoms%n_hia),jspins)

         INTEGER,                INTENT(IN)     :: jspins

         INTEGER i, m,mp, l, i_hia, ispin, n, it,is, isi, natom, nn
         REAL imag, re, fac, n_l
         COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),nr_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
         COMPLEX n1_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const), d_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)


         mmpMat(:,:,:,:) = CMPLX(0.0,0.0)

         IF(this%l_ef) THEN
            DO i_hia = 1, atoms%n_hia
               n_l = 0.0
               l = atoms%lda_hia(i_hia)%l
               n = atoms%lda_hia(i_hia)%atomType  
  

               DO ispin = 1, jspins
                  n_tmp(:,:) = CMPLX(0.0,0.0)
                  DO m = -l, l
                     DO mp = -l, l
                        DO i = 1, this%nz
                           IF(this%nr(i_hia).NE.1) THEN

                              CALL intgr3(REAL(this%gmmpMat(:,i,i_hia,m,mp,ispin)),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),re)
                              CALL intgr3(AIMAG(this%gmmpMat(:,i,i_hia,m,mp,ispin)),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imag)

                              n_tmp(m,mp) = n_tmp(m,mp) + AIMAG((re+ImagUnit*imag)*this%de(i))

                           ELSE

                              n_tmp(m,mp) = n_tmp(m,mp) + AIMAG(this%gmmpMat(1,i,i_hia,m,mp,ispin)*this%de(i))
                           
                           END IF
                        ENDDO

                        mmpMat(m,mp,i_hia,ispin) = -1/pi_const * n_tmp(m,mp)
                        IF(m.EQ.mp) n_l = n_l + mmpMat(m,mp,i_hia,ispin)

                     ENDDO
                  ENDDO
               ENDDO
               WRITE(*,*) "OCCUPATION: ", n_l
            ENDDO 

         ELSE

            CALL juDFT_error("Green's Function is calculated on a contour not ending at Efermi", calledby="calc_mmpmat")
         
         ENDIF

      END SUBROUTINE calc_mmpmat

END MODULE m_types_greensf