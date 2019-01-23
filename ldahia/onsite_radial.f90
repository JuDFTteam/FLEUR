MODULE m_gOnsite_radial

   CONTAINS

      SUBROUTINE onsite_coeffs(atoms,ispin,noccbd,ikpt,eig,eigVecCoeffs,gOnsite)

         !This Subroutine calculates the Coefficients of the imaginary part of the
         !crystal greens function

         !We enforce l = lprime

         USE m_types
         USE m_constants

         IMPLICIT NONE

         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_eigVecCoeffs),   INTENT(IN)     :: eigVecCoeffs
         TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite

         INTEGER,                INTENT(IN)     :: ispin
         INTEGER,                INTENT(IN)     :: noccbd
         INTEGER,                INTENT(IN)     :: ikpt
         REAL,                   INTENT(IN)     :: eig(noccbd)

         INTEGER i_hia, l,n,natom,nn,i,m,lm,mp,lmp
         REAL fac


         DO i_hia = 1, atoms%n_hia

            l = atoms%lda_hia(i_hia)%l
            n = atoms%lda_hia(i_hia)%atomType

            natom = 0
            DO i = 1, n-1
               natom = natom + atoms%neq(i)
            ENDDO


            DO nn = 1, atoms%neq(n)
               natom = natom +1
               fac = 1./atoms%neq(n)
               !calculate the coefficients for the radial functions
               DO i = 1,noccbd
                  IF(eig(i).LT.gonsite%e_top) THEN
                     DO m = -l, l
                        lm = l*(l+1)+m
                        DO mp = -l, l
                           lmp = l*(l+1)+mp    
                           gOnsite%uu(i,ikpt,i_hia,m,mp,ispin) = gOnsite%uu(i,ikpt,i_hia,m,mp,ispin) - pi_const * fac *&
                                                    conjg(eigVecCoeffs%acof(i,lm,natom,ispin))*eigVecCoeffs%acof(i,lmp,natom,ispin)      
                           gOnsite%dd(i,ikpt,i_hia,m,mp,ispin) = gOnsite%dd(i,ikpt,i_hia,m,mp,ispin) - pi_const * fac *&
                                                    conjg(eigVecCoeffs%bcof(i,lm,natom,ispin))*eigVecCoeffs%bcof(i,lmp,natom,ispin) 
                           gOnsite%du(i,ikpt,i_hia,m,mp,ispin) = gOnsite%du(i,ikpt,i_hia,m,mp,ispin) - pi_const * fac *&
                                                    (conjg(eigVecCoeffs%acof(i,lm,natom,ispin))*eigVecCoeffs%bcof(i,lmp,natom,ispin)+&
                                                    conjg(eigVecCoeffs%bcof(i,lm,natom,ispin))*eigVecCoeffs%acof(i,lmp,natom,ispin))  
                        ENDDO
                     ENDDO
                     !
                     !  TODO: add local orbital contribution; treat local symmetries
                     !
                  ENDIF
               ENDDO

            ENDDO
         ENDDO
      END SUBROUTINE onsite_coeffs

      SUBROUTINE calc_onsite_radial(atoms,enpara,vr,jspin,jspins,neigd,ntetra,nkpt,itetra,voltet,nevk,eigv,gOnsite,ef,sym)

         USE m_types
         USE m_constants
         USE m_juDFT
         USE m_tetra
         USE m_smooth
         USE m_kkintgr
         USE m_radfun

         IMPLICIT NONE

         TYPE(t_atoms),          INTENT(IN)     :: atoms
         TYPE(t_enpara),         INTENT(IN)     :: enpara
         TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite
         TYPE(t_sym),            INTENT(IN)     :: sym

         INTEGER,                INTENT(IN)     :: neigd,jspin,jspins
         INTEGER,                INTENT(IN)     :: ntetra,nkpt


         INTEGER,                INTENT(IN)     :: itetra(4,6*nkpt),nevk(nkpt)
         REAL,                   INTENT(IN)     :: voltet(6*nkpt)
         REAL,                   INTENT(INOUT)  :: eigv(neigd,nkpt)

         REAL,                   INTENT(IN)     :: ef
         REAL,                   INTENT(IN)     :: vr(atoms%jmtd,atoms%ntype,jspins)


         TYPE(t_usdus) usdus
         INTEGER i, i_hia, l, m, mp, jr, noded, nodeu, n, j
         REAL wronk


         REAL,    ALLOCATABLE :: f(:,:,:,:),g(:,:,:,:)

         ALLOCATE ( f(atoms%jmtd,2,0:atoms%lmaxd,jspins) )
         ALLOCATE ( g(atoms%jmtd,2,0:atoms%lmaxd,jspins) )

         CALL usdus%init(atoms,jspins)



         DO i_hia = 1,atoms%n_hia
            l = atoms%lda_hia(i_hia)%l
            n = atoms%lda_hia(i_hia)%atomType

            CALL radfun(l,n,jspin,enpara%el0(l,n,jspin),vr(:,n,jspin),atoms,&
                  f(:,:,l,jspin),g(:,:,l,jspin),usdus, nodeu,noded,wronk)

            DO m = -l, l
               DO mp = -l,l
                  !calculate the imaginary part if we use the tetrahedron method
                  IF(gOnsite%l_tetra) THEN
                     CALL timestart("On-Site: Tetrahedron method")
                     CALL int_tetra(nkpt,ntetra,itetra,voltet,neigd,nevk,gOnsite%uu(:,:,i_hia,m,mp,jspin)&
                                                ,eigv,gOnsite%ne,gOnsite%im_uu(:,i_hia,m,mp,jspin),gOnsite%e_top,gOnsite%e_bot)
                     CALL int_tetra(nkpt,ntetra,itetra,voltet,neigd,nevk,gOnsite%dd(:,:,i_hia,m,mp,jspin)&
                                                ,eigv,gOnsite%ne,gOnsite%im_dd(:,i_hia,m,mp,jspin),gOnsite%e_top,gOnsite%e_bot)
                     CALL int_tetra(nkpt,ntetra,itetra,voltet,neigd,nevk,gOnsite%du(:,:,i_hia,m,mp,jspin)&
                                                ,eigv,gOnsite%ne,gOnsite%im_du(:,i_hia,m,mp,jspin),gOnsite%e_top,gOnsite%e_bot)
                     IF(jspins.EQ.1) THEN
                        gOnsite%im_uu(:,i_hia,m,mp,1) = 2.0 * gOnsite%im_uu(:,i_hia,m,mp,1)
                        gOnsite%im_dd(:,i_hia,m,mp,1) = 2.0 * gOnsite%im_dd(:,i_hia,m,mp,1)
                        gOnsite%im_du(:,i_hia,m,mp,1) = 2.0 * gOnsite%im_du(:,i_hia,m,mp,1)
                     END IF
                     CALL timestop("On-Site: Tetrahedron method")
                  ENDIF

                  DO jr = 1, atoms%jri(n)
                     DO j = 1, gOnsite%ne
                        gOnsite%im_gmmpmat(jr,j,i_hia,m,mp,jspin) = gOnsite%im_gmmpmat(jr,j,i_hia,m,mp,jspin) + &
                                          gOnsite%im_uu(j,i_hia,m,mp,jspin) * (f(jr,1,l,jspin)*f(jr,1,l,jspin)+f(jr,2,l,jspin)*f(jr,2,l,jspin)) +&
                                          gOnsite%im_dd(j,i_hia,m,mp,jspin) * (g(jr,1,l,jspin)*g(jr,1,l,jspin)+g(jr,2,l,jspin)*g(jr,2,l,jspin)) +&
                                          gOnsite%im_du(j,i_hia,m,mp,jspin) * (f(jr,1,l,jspin)*g(jr,1,l,jspin)+f(jr,2,l,jspin)*g(jr,2,l,jspin))

                     ENDDO
                  ENDDO  
               ENDDO
            ENDDO
            !
            !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration 
            !
            CALL greensf_cutoff_radial(gOnsite,atoms,jspins)

            CALL gOnsite%init_e_contour(ef)

            CALL timestart("On-Site: Kramer-Kronigs-Integration")
            DO m= -l,l
               DO mp= -l,l
                  DO jr = 1, gOnsite%nr(i_hia)
                     IF(gOnsite%mode.EQ.1) THEN
                        CALL kkintgr_real(gOnsite%nz,gOnsite%e(:),gOnsite%ne,gOnsite%sigma,gOnsite%del,gOnsite%e_bot,&
                                          gOnsite%im_gmmpMat(jr,:,i_hia,m,mp,jspin),gOnsite%gmmpMat(jr,:,i_hia,m,mp,jspin))
                     ELSE IF(gOnsite%mode.EQ.2) THEN
                        CALL kkintgr_complex(gOnsite%nz,gOnsite%e(:),gOnsite%ne,gOnsite%sigma,gOnsite%del,gOnsite%e_bot,&
                                             gOnsite%im_gmmpMat(jr,:,i_hia,m,mp,jspin),gOnsite%gmmpMat(jr,:,i_hia,m,mp,jspin))
                     END IF
                  ENDDO
               ENDDO
            ENDDO
            CALL timestop("On-Site: Kramer-Kronigs-Integration")

            CALL gOnsite%calc_mmpmat(atoms,jspins)


         ENDDO

      END SUBROUTINE calc_onsite_radial


      SUBROUTINE greensf_cutoff_radial(gOnsite,atoms,jspins)
         !This Subroutine detemines the cutoff energy for the kramers-kronig-integration
         !This cutoff energy is defined so that the integral over the fDOS up to this cutoff 
         !is equal to 14 (the number of states in the 4f shell) or not to small


         USE m_types
         USE m_intgr
         USE m_juDFT
         USE m_constants
         USE m_kkintgr
         
         IMPLICIT NONE

         TYPE(t_greensf),     INTENT(INOUT)  :: gOnsite
         TYPE(t_atoms),       INTENT(IN)     :: atoms
         INTEGER,             INTENT(IN)     :: jspins


         REAL, ALLOCATABLE :: fDOS(:)

         INTEGER i_hia, i,l, m,mp, ispin, j, n_c, kkintgr_cut, jr, n

         REAL integral

         REAL a,b, imag
         LOGICAL l_write

         ALLOCATE(fDOS(gOnsite%ne))

         l_write=.true.

         DO i_hia = 1, atoms%n_hia

            fDOS(:) = 0.0
            l = atoms%lda_hia(i_hia)%l   
            n = atoms%lda_hia(i_hia)%atomType      

            !Calculate the trace over m,mp of the Greens-function matrix to obtain the fDOS 

            !n_f(e) = -1/pi * TR[Im(G_f(e))]
            DO ispin = 1, jspins
               DO m = -l , l
                  DO j = 1, gOnsite%ne
                     IF(gOnsite%nr(i_hia).NE.1) THEN
                        CALL intgr3(gOnsite%im_gmmpMat(:,j,i_hia,m,m,ispin),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imag)
                     ELSE
                        imag = gOnsite%im_gmmpMat(1,j,i_hia,m,m,ispin)
                     END IF
                     fDOS(j) = fDOS(j) + imag
                  ENDDO
               ENDDO
            ENDDO

            fDOS(:) = -1/pi_const * fDOS(:)


            IF(l_write) THEN

               open(unit=1337,file="fDOS.txt",status="replace", action="write")

               DO j = 1, gOnsite%ne

                  WRITE(1337,*) (j-1) * gOnsite%del + gOnsite%e_bot, fDOS(j)

               ENDDO

               close(unit=1337)

            END IF       



            CALL trapz(fDOS(:), gOnsite%del, gOnsite%ne, integral)

            WRITE(*,*) "Integral over fDOS: ", integral

            kkintgr_cut = gOnsite%ne

            IF(integral.LT.13.5) THEN
               ! If the integral is to small we stop here to avoid problems
               CALL juDFT_error("fDOS-integral too small: make sure numbands is big enough", calledby="greensf_cutoff")
               
            ELSE IF((integral.GT.14).AND.((integral-14).GT.0.001)) THEN
               !IF the integral is bigger than 14, search for the cutoff using the bisection method   

               a = gOnsite%e_bot
               b = gOnsite%e_top


               DO

                  n_c = INT(((a+b)/2.0-gOnsite%e_bot)/gOnsite%del)+1
                  
                  CALL trapz(fDOS(:),gOnsite%del,n_c,integral)

                  IF((ABS(integral-14).LT.0.001).OR.(ABS(a-b)/2.0.LT.gonsite%del)) THEN

                     kkintgr_cut = INT(((a+b)/2.0-gOnsite%e_bot)/gOnsite%del)+1
                     EXIT

                  ELSE IF((integral-14).LT.0) THEN
                     a = (a+b)/2.0
                  ELSE IF((integral-14).GT.0) THEN
                     b = (a+b)/2.0
                  END IF

               ENDDO

               CALL trapz(fDOS(:),gOnsite%del,kkintgr_cut,integral)

               WRITE(*,*) "CALCULATED CUTOFF: ", kkintgr_cut
               WRITE(*,*) "INTEGRAL OVER fDOS with cutoff: ", integral
            END IF


            !Now we set the imaginary part of the greens function to zero above this cutoff
            DO ispin = 1, jspins
               DO i = kkintgr_cut+1, gOnsite%ne
                  DO m = -l, l
                     DO mp = -l,l
                        DO jr = 1, gOnsite%nr(i_hia)
                           gOnsite%im_gmmpMat(jr,i,i_hia,m,mp,ispin) = 0.0e0
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

         ENDDO



      END SUBROUTINE greensf_cutoff_radial

END MODULE m_gOnsite_radial