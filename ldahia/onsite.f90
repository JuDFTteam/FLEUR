MODULE m_gOnsite

CONTAINS

SUBROUTINE calc_qalmmpMat(atoms,ispin,noccbd,ikpt,usdus,eigVecCoeffs,gOnsite)

   !This Subroutine is used if we want to use the tetrahedron method 
   !It calcualtes the weights for the individual k-points for the integration,
   !which is performed in the subroutine calc_onsite

   USE m_types

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_usdus),          INTENT(IN)     :: usdus
   TYPE(t_eigVecCoeffs),   INTENT(IN)     :: eigVecCoeffs
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite

   INTEGER,                INTENT(IN)     :: noccbd
   INTEGER,                INTENT(IN)     :: ikpt
   INTEGER,                INTENT(IN)     :: ispin

   INTEGER i_hia, n, l, natom, i,nn,m,lm,mp,lmp
   REAL fac




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
         fac = 1./atoms%neq(n)
         DO m = -l, l
            lm = l*(l+1)+m
            DO mp = -l,l
               lmp = l*(l+1)+mp
               DO i = 1,noccbd
                  gOnsite%qalmmpMat(i,ikpt,i_hia,m,mp,ispin) = gOnsite%qalmmpMat(i,ikpt,i_hia,m,mp,ispin) - fac *&
                                                          (conjg(eigVecCoeffs%acof(i,lmp,natom,ispin))*eigVecCoeffs%acof(i,lm,natom,ispin) +&
                                                           conjg(eigVecCoeffs%bcof(i,lmp,natom,ispin))*eigVecCoeffs%bcof(i,lm,natom,ispin) *&
                                                           usdus%ddn(l,n,ispin))
               ENDDO 
            ENDDO
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE calc_qalmmpMat


SUBROUTINE im_gmmpMathist(atoms,ispin,jspins,noccbd,wtkpt,eig,usdus,eigVecCoeffs,gOnsite)

   !This Subroutine calculates the imaginary part of the Matrix elements G^[n \sigma]_{Lm Lm'}(E+i*sigma)
   !at the current k-Point (it is called in cdnval) inside the MT-sphere (averaged over the radial component)
   !and sums over the Brillouin-Zone using the histogram method

   USE m_types
   USE m_constants

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_eigVecCoeffs),   INTENT(IN)     :: eigVecCoeffs
   TYPE(t_usdus),          INTENT(IN)     :: usdus
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite

   INTEGER,                INTENT(IN)     :: ispin
   INTEGER,                INTENT(IN)     :: jspins
   INTEGER,                INTENT(IN)     :: noccbd

   REAL,                   INTENT(IN)     :: wtkpt
   REAL,                   INTENT(IN)     :: eig(noccbd)

   INTEGER i_hia, i, j, n, nn, natom, l, m, mp, lm, lmp
   REAL fac, wk, del

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
         fac = 1./atoms%neq(n)
         DO m = -l, l
            lm = l*(l+1)+m
            DO mp = -l,l
               lmp = l*(l+1)+mp
               DO i = 1,noccbd

                  j = NINT((eig(i)-gOnsite%e_bot)/del)+1
                  IF( (j.LE.gOnsite%ne).AND.(j.GE.1) ) THEN
                     gOnsite%gmmpMat(j,i_hia,m,mp,ispin) = gOnsite%gmmpMat(j,i_hia,m,mp,ispin) - ImagUnit * wk * fac *&
                                                          (conjg(eigVecCoeffs%acof(i,lmp,natom,ispin))*eigVecCoeffs%acof(i,lm,natom,ispin) +&
                                                           conjg(eigVecCoeffs%bcof(i,lmp,natom,ispin))*eigVecCoeffs%bcof(i,lm,natom,ispin) *&
                                                           usdus%ddn(l,n,ispin))
                  ENDIF

               ENDDO 
            ENDDO
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE im_gmmpMathist

SUBROUTINE calc_onsite(atoms,jspin,jspins,neigd,ntetra,nkpt,itetra,voltet,nevk,eigv,gOnsite)

   USE m_types
   USE m_constants
   USE m_tetra
   USE m_smooth
   USE m_kkintgr

   IMPLICIT NONE

   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite

   INTEGER,                INTENT(IN)     :: neigd,jspin,jspins
   INTEGER,                INTENT(IN)     :: ntetra,nkpt


   INTEGER,                INTENT(IN)     :: itetra(4,6*nkpt),nevk(nkpt)
   REAL,                   INTENT(IN)     :: voltet(6*nkpt)
   REAL,                   INTENT(INOUT)  :: eigv(neigd,nkpt)

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
                                          ,eigv,gOnsite%ne,gOnsite%gmmpMat(:,i_hia,m,mp,jspin),gOnsite%e_top,gOnsite%e_bot)
               IF(jspins.EQ.1) gOnsite%gmmpMat(:,i_hia,m,mp,1) = 2.0 * gOnsite%gmmpMat(:,i_hia,m,mp,1)
            END IF
            
            !smooth the imaginary part using gaussian broadening

            im = AIMAG(gOnsite%gmmpMat(:,i_hia,m,mp,jspin))
            CALL smooth(e(:),im(:),gOnsite%sigma,gOnsite%ne)
            gOnsite%gmmpMat(:,i_hia,m,mp,jspin) = ImagUnit * im(:)


         ENDDO
      ENDDO

      !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration 

      CALL greensf_cutoff(gOnsite,atoms)

      DO m= -l,l
         DO mp= -l,l
            im = AIMAG(gOnsite%gmmpMat(:,i_hia,m,mp,jspin))
            CALL kkintgr(gOnsite%ne,gOnsite%sigma,del,im(:),gOnsite%gmmpMat(:,i_hia,m,mp,jspin))
         ENDDO
      ENDDO

   ENDDO

END SUBROUTINE calc_onsite


SUBROUTINE greensf_cutoff(gOnsite,atoms)
   !This Subroutine detemines the cutoff energy for the kramers-kronig-integration
   !This cutoff energy is defined so that the integral over the fDOS up to this cutoff 
   !is equal to 14 (the number of states in the 4f shell) or not to small


   USE m_types
   USE m_intgr
   USE m_juDFT
   
   IMPLICIT NONE

   TYPE(t_greensf),     INTENT(INOUT)  :: gOnsite
   TYPE(t_atoms),       INTENT(IN)     :: atoms


   REAL, ALLOCATABLE :: fDOS(:)

   INTEGER i, l, m, jspin, j, n_c

   REAL integral, del

   REAL a,b

   ALLOCATE(fDOS(gOnsite%ne))
   fDOS(:) = 0.0

   DO i = 1, atoms%n_hia

      l = atoms%lda_hia(i)%l      


      jspin = 1 !TEMPORARY!!!
      !Calculate the trace over m,mp of the Greens-function matrix to obtain the fDOS
      DO m = -l , l
         DO j = 1, gOnsite%ne
            fDOS(j) = fDOS(j) - AIMAG(gOnsite%gmmpMat(j,i,m,m,jspin))
         ENDDO
      ENDDO

      del = (gOnsite%e_top-gOnsite%e_bot)/REAL(gOnsite%ne-1)

      CALL intgz0(fDOS(:), del, gOnsite%ne, integral,.false.)

      WRITE(*,*) "Integral over fDOS: ", integral

      gOnsite%kkintgr_cut = gOnsite%e_top

      IF(integral.LT.13.5) THEN
         ! If the integral is to small we stop here to avoid problems
         CALL juDFT_error("fDOS-integral too small", calledby="greensf_cutoff")
         
      ELSE IF((integral.GT.14).AND.((integral-14).GT.0.1)) THEN
         !IF the integral is bigger than 14, search for the cutoff using the bisection method   

         a = gOnsite%e_bot
         b = gOnsite%e_top


         DO

            n_c = INT(((a+b)/2.0-gOnsite%e_bot)/del)+1
            
            CALL intgz0(fDOS(:),del,n_c,integral,.false.)

            IF((ABS(integral-14).LT.0.01).OR.(ABS(a-b)/2.0.LT.del)) THEN

               gOnsite%kkintgr_cut = (a+b)/2.0
               EXIT

            ELSE IF((integral-14).LT.0) THEN
               a = (a+b)/2.0
            ELSE IF((integral-14).GT.0) THEN
               b = (a+b)/2.0
            END IF

         ENDDO

         n_c = INT((gOnsite%kkintgr_cut-gOnsite%e_bot)/del)+1

         CALL intgz0(fDOS(:),del,n_c,integral,.false.)

         WRITE(*,*) "CALCULATED CUTOFF: ", gOnsite%kkintgr_cut
         WRITE(*,*) "INTEGRAL OVER fDOS with cutoff: ", integral
         WRITE(*,*) "fDOS at the cutoff", fDOS(n_c)
   END IF
   ENDDO



END SUBROUTINE greensf_cutoff


END MODULE m_gOnsite