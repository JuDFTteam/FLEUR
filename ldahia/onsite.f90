MODULE m_gOnsite

CONTAINS

SUBROUTINE calc_qalmmpMat(atoms,sym,ispin,noccbd,ikpt,usdus,eigVecCoeffs,gOnsite)

   !This Subroutine is used if we want to use the tetrahedron method 
   !It calculates the weights for the individual k-points for the integration,
   !which is performed in the subroutine calc_onsite

   !It is essentially the f-density of states in a (m,mp) matrix with an additional factor - pi

   USE m_types
   USE m_constants

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_sym),            INTENT(IN)     :: sym
   TYPE(t_usdus),          INTENT(IN)     :: usdus
   TYPE(t_eigVecCoeffs),   INTENT(IN)     :: eigVecCoeffs
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite

   INTEGER,                INTENT(IN)     :: noccbd
   INTEGER,                INTENT(IN)     :: ikpt
   INTEGER,                INTENT(IN)     :: ispin

   INTEGER i_hia, n, l, natom, i,nn,m,lm,mp,lmp, it,is, isi
   REAL fac

   COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),nr_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
   COMPLEX n1_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const), d_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)



!One HIA per atom at the moment

   DO i_hia = 1, atoms%n_hia
      n = atoms%lda_hia(i_hia)%atomType
      l = atoms%lda_hia(i_hia)%l

      !finding the right starting index
      natom = 0
      DO i = 1, n-1
         natom = natom +atoms%neq(i)
      ENDDO

      DO nn = 1, atoms%neq(n)
         natom = natom +1
         DO i = 1, noccbd
            n_tmp(:,:) = cmplx(0.0,0.0)
            DO m = -l, l
               lm = l*(l+1)+m
               DO mp = -l,l
                  lmp = l*(l+1)+mp

                     n_tmp(m,mp) = n_tmp(m,mp) -  pi_const*&
                                  (conjg(eigVecCoeffs%acof(i,lmp,natom,ispin))*eigVecCoeffs%acof(i,lm,natom,ispin) +&
                                   conjg(eigVecCoeffs%bcof(i,lmp,natom,ispin))*eigVecCoeffs%bcof(i,lm,natom,ispin) *&
                                   usdus%ddn(l,n,ispin))
               ENDDO
            ENDDO

            !TODO: add local orbital contribution

            !
            !  n_mmp should be rotated by D_mm' ; compare force_a21; taken from n_mat.f90
            !
            DO it = 1, sym%invarind(natom)

                fac = 1.0  /  ( sym%invarind(natom) * atoms%neq(n) )
                is = sym%invarop(natom,it)
                isi = sym%invtab(is)
                d_tmp(:,:) = cmplx(0.0,0.0)
                DO m = -l,l
                   DO mp = -l,l
                      d_tmp(m,mp) = sym%d_wgn(m,mp,l,isi)
                   ENDDO
                ENDDO
                nr_tmp = matmul( transpose( conjg(d_tmp) ) , n_tmp)
                n1_tmp =  matmul( nr_tmp, d_tmp )


                DO m = -l,l
                   DO mp = -l,l
                      gOnsite%qalmmpMat(i,ikpt,i_hia,m,mp,ispin) = gOnsite%qalmmpMat(i,ikpt,i_hia,m,mp,ispin) + conjg(n1_tmp(m,mp)) * fac
                   ENDDO
                ENDDO

             ENDDO

         ENDDO !loop over bands
      ENDDO !loop over equivalent atoms
   ENDDO !loop over number of DFT+HIAs



END SUBROUTINE calc_qalmmpMat


SUBROUTINE im_gmmpMathist(atoms,sym,ispin,jspins,noccbd,wtkpt,eig,usdus,eigVecCoeffs,gOnsite)

   !This Subroutine calculates the imaginary part of the Matrix elements G^[n \sigma]_{Lm Lm'}(E+i*sigma)
   !at the current k-Point (it is called in cdnval) inside the MT-sphere (averaged over the radial component)
   !and sums over the Brillouin-Zone using the histogram method

   !It is essentially the f-density of states in a (m,mp) matrix with an additional factor - pi

   USE m_types
   USE m_constants

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_sym),            INTENT(IN)     :: sym
   TYPE(t_eigVecCoeffs),   INTENT(IN)     :: eigVecCoeffs
   TYPE(t_usdus),          INTENT(IN)     :: usdus
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite

   INTEGER,                INTENT(IN)     :: ispin
   INTEGER,                INTENT(IN)     :: jspins
   INTEGER,                INTENT(IN)     :: noccbd

   REAL,                   INTENT(IN)     :: wtkpt
   REAL,                   INTENT(IN)     :: eig(noccbd)

   INTEGER i_hia, i, j, n, nn, natom, l, m, mp, lm, lmp, it,is, isi
   REAL fac, wk, del

   COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),nr_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
   COMPLEX n1_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const), d_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

   del = (gOnsite%e_top - gOnsite%e_bot)/REAL(gOnsite%ne-1)
   wk = 2.*wtkpt/REAL(jspins)/del

   !One HIA per atom at the moment

   DO i_hia = 1, atoms%n_hia
      n = atoms%lda_hia(i_hia)%atomType
      l = atoms%lda_hia(i_hia)%l

      !finding the right starting index
      natom = 0
      DO i = 1, n-1
         natom = natom +atoms%neq(i)
      ENDDO

      DO nn = 1, atoms%neq(n)
         natom = natom +1
         DO i = 1, noccbd
            j = NINT((eig(i)-gOnsite%e_bot)/del)+1
            IF( (j.LE.gOnsite%ne).AND.(j.GE.1) ) THEN
               n_tmp(:,:) = cmplx(0.0,0.0)
               DO m = -l, l
                  lm = l*(l+1)+m
                  DO mp = -l,l
                     lmp = l*(l+1)+mp

                        n_tmp(m,mp) = n_tmp(m,mp) -  pi_const*&
                                     (conjg(eigVecCoeffs%acof(i,lmp,natom,ispin))*eigVecCoeffs%acof(i,lm,natom,ispin) +&
                                      conjg(eigVecCoeffs%bcof(i,lmp,natom,ispin))*eigVecCoeffs%bcof(i,lm,natom,ispin) *&
                                      usdus%ddn(l,n,ispin))
                  ENDDO
               ENDDO

               !TODO: add local orbital contribution

               !
               !  n_mmp should be rotated by D_mm' ; compare force_a21; taken from n_mat.f90
               !
               DO it = 1, sym%invarind(natom)

                   fac = 1.0  /  ( sym%invarind(natom) * atoms%neq(n) )
                   is = sym%invarop(natom,it)
                   isi = sym%invtab(is)
                   d_tmp(:,:) = cmplx(0.0,0.0)
                   DO m = -l,l
                      DO mp = -l,l
                         d_tmp(m,mp) = sym%d_wgn(m,mp,l,isi)
                      ENDDO
                   ENDDO
                   nr_tmp = matmul( transpose( conjg(d_tmp) ) , n_tmp)
                   n1_tmp =  matmul( nr_tmp, d_tmp )


                   DO m = -l,l
                      DO mp = -l,l
                         gOnsite%im_gmmpMat(j,i_hia,m,mp,ispin) = gOnsite%im_gmmpMat(j,i_hia,m,mp,ispin) + conjg(n1_tmp(m,mp)) * fac * wk
                      ENDDO
                   ENDDO

                ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE im_gmmpMathist

SUBROUTINE calc_onsite(atoms,jspin,jspins,neigd,ntetra,nkpt,itetra,voltet,nevk,eigv,gOnsite,ef,sym)

   USE m_types
   USE m_constants
   USE m_tetra
   USE m_smooth
   USE m_kkintgr
   USE m_occupation

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite
   TYPE(t_sym),            INTENT(IN)     :: sym

   INTEGER,                INTENT(IN)     :: neigd,jspin,jspins
   INTEGER,                INTENT(IN)     :: ntetra,nkpt


   INTEGER,                INTENT(IN)     :: itetra(4,6*nkpt),nevk(nkpt)
   REAL,                   INTENT(IN)     :: voltet(6*nkpt)
   REAL,                   INTENT(INOUT)  :: eigv(neigd,nkpt)

   REAL,                   INTENT(IN)     :: ef

   INTEGER i, i_hia, l, m, mp
   REAL del

   REAL, ALLOCATABLE :: e(:), im(:)
   !construct an energy grid for smoothing (required by the function in m_smooth)

   ALLOCATE (e(gOnsite%ne))
   ALLOCATE (im(gOnsite%ne))
   
   del = (gOnsite%e_top - gOnsite%e_bot)/REAL(gOnsite%ne-1)
   DO i = 1, gOnsite%ne
      e(i) = del * (i-1) + gOnsite%e_bot
   ENDDO

   DO i_hia = 1,atoms%n_hia
      l = atoms%lda_hia(i_hia)%l
      DO m = -l, l
         DO mp = -l,l

            !calculate the imaginary part if we use the tetrahedron method
            IF(gOnsite%l_tetra) THEN
               CALL int_tetra(nkpt,ntetra,itetra,voltet,neigd,nevk,gOnsite%qalmmpMat(:,:,i_hia,m,mp,jspin)&
                                          ,eigv,gOnsite%ne,gOnsite%im_gmmpMat(:,i_hia,m,mp,jspin),gOnsite%e_top,gOnsite%e_bot)
               IF(jspins.EQ.1) gOnsite%im_gmmpMat(:,i_hia,m,mp,1) = 2.0 * gOnsite%im_gmmpMat(:,i_hia,m,mp,1)
            END IF
            
            !smooth the imaginary part using gaussian broadening 
            IF(gOnsite%sigma.NE.0.0) CALL smooth(e(:),gOnsite%im_gmmpMat(:,i_hia,m,mp,jspin),gOnsite%sigma,gOnsite%ne)

         ENDDO
      ENDDO


      !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration 

      CALL greensf_cutoff(gOnsite,atoms,jspins)

      CALL energy_contour(gOnsite,ef)

      DO m= -l,l
         DO mp= -l,l
            IF(gOnsite%mode.EQ.1) THEN
               CALL kkintgr_real(gOnsite%nz,gOnsite%e(:),gOnsite%ne,gOnsite%sigma,del,gOnsite%e_bot,&
                                 gOnsite%im_gmmpMat(:,i_hia,m,mp,jspin),gOnsite%gmmpMat(:,i_hia,m,mp,jspin))
            ELSE IF(gOnsite%mode.EQ.2) THEN
               CALL kkintgr_complex(gOnsite%nz,gOnsite%e(:),gOnsite%ne,gOnsite%sigma,del,gOnsite%e_bot,&
                                    gOnsite%im_gmmpMat(:,i_hia,m,mp,jspin),gOnsite%gmmpMat(:,i_hia,m,mp,jspin))
            END IF
         ENDDO
      ENDDO

      CALL calc_occ_from_g(gOnsite,atoms,jspins,ef,sym)


   ENDDO

END SUBROUTINE calc_onsite


SUBROUTINE greensf_cutoff(gOnsite,atoms,jspins)
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

   INTEGER i_hia, i,l, m,mp, ispin, j, n_c, kkintgr_cut

   REAL integral, del

   REAL a,b
   LOGICAL l_write

   ALLOCATE(fDOS(gOnsite%ne))

   l_write=.true.

   DO i_hia = 1, atoms%n_hia

      fDOS(:) = 0.0
      l = atoms%lda_hia(i_hia)%l      

      !Calculate the trace over m,mp of the Greens-function matrix to obtain the fDOS 

      !n_f(e) = -1/pi * TR[Im(G_f(e))]
      DO ispin = 1, jspins
         DO m = -l , l
            DO j = 1, gOnsite%ne
               fDOS(j) = fDOS(j) + gOnsite%im_gmmpMat(j,i_hia,m,m,ispin)
            ENDDO
         ENDDO
      ENDDO

      fDOS(:) = -1/pi_const * fDOS(:)

      del = (gOnsite%e_top-gOnsite%e_bot)/REAL(gOnsite%ne-1)

      IF(l_write) THEN

         open(unit=1337,file="fDOS.txt",status="replace", action="write")

         DO j = 1, gOnsite%ne

            WRITE(1337,*) (j-1) * del + gOnsite%e_bot, fDOS(j)

         ENDDO

         close(unit=1337)

      END IF       



      CALL trapz(fDOS(:), del, gOnsite%ne, integral)

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

            n_c = INT(((a+b)/2.0-gOnsite%e_bot)/del)+1
            
            CALL trapz(fDOS(:),del,n_c,integral)

            IF((ABS(integral-14).LT.0.001).OR.(ABS(a-b)/2.0.LT.del)) THEN

               kkintgr_cut = INT(((a+b)/2.0-gOnsite%e_bot)/del)+1
               EXIT

            ELSE IF((integral-14).LT.0) THEN
               a = (a+b)/2.0
            ELSE IF((integral-14).GT.0) THEN
               b = (a+b)/2.0
            END IF

         ENDDO

         CALL trapz(fDOS(:),del,kkintgr_cut,integral)

         WRITE(*,*) "CALCULATED CUTOFF: ", kkintgr_cut
         WRITE(*,*) "INTEGRAL OVER fDOS with cutoff: ", integral
      END IF


      !Now we set the imaginary part of the greens function to zero above this cutoff
      DO ispin = 1, jspins
         DO i = kkintgr_cut+1, gOnsite%ne
            DO m = -l, l
               DO mp = -l,l
                  gOnsite%im_gmmpMat(i,i_hia,m,mp,ispin) = 0.0e0
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   ENDDO



END SUBROUTINE greensf_cutoff


END MODULE m_gOnsite