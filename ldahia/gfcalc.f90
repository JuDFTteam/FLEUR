MODULE m_gfcalc

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_gfcalc
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  This module contains various subroutines calculating properties from the
   !>  green's functions:
   !>       -calculates the effective exchange interaction from the onsite
   !>          green's function according to Condens. Matter 26 (2014) 476003 EQ.1
   !>       -calculates the occuaption matrix 
   !
   ! REVISION HISTORY:
   ! February 2019 - Initial Version
   ! March    2019 - Changed calculation of the onsite exchange matrix
   !------------------------------------------------------------------------------

   USE m_juDFT
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE eff_excinteraction(g0,atoms,input,j0,onsite_exc_split)

      TYPE(t_greensf),        INTENT(IN)  :: g0
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      COMPLEX,                INTENT(OUT) :: j0
      TYPE(t_input),          INTENT(IN)  :: input
      REAL,                   INTENT(IN)  :: onsite_exc_split

      COMPLEX tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),integrand
      INTEGER i,iz,m,l,mp,ispin,n,i_gf,matsize,ipm
      LOGICAL l_matinv

      TYPE(t_mat) :: calcup,calcdwn
      TYPE(t_mat) :: delta
      TYPE(t_mat) :: calc

      l_matinv = .true. !Determines how the onsite exchange splitting is calculated
      
      DO i_gf = 1, atoms%n_gf
         j0 = 0.0
         l = atoms%onsiteGF(i_gf)%l
         n = atoms%onsiteGF(i_gf)%atomType
         matsize = 2*l+1
         CALL calcup%init(.false.,matsize,matsize)
         CALL calcdwn%init(.false.,matsize,matsize)
         CALL delta%init(.false.,matsize,matsize)
         CALL calc%init(.false.,matsize,matsize)
         IF(.NOT.l_matinv) THEN
            DO i = 1, matsize
               delta%data_c(i,i) = onsite_exc_split
            ENDDO
         ENDIF

         DO iz = 1, g0%nz
            !
            !calculate the onsite exchange matrix if we use matrix inversion
            !
            IF(l_matinv) THEN
               !---------------------------------------------
               !\Delta = (G_up)^-1-(G_down)^-1
               !---------------------------------------------
               !Symmetrize the green's function for up/down 
               !spin with respect to the complex plane
               !Here we assume that the onsite Hamiltonian
               !is real
               !---------------------------------------------
               !G^(up/down)^-1 = 1/2 * (G+^(up/down)^-1 + G-^(up/down)^-1)
               !---------------------------------------------
               delta%data_c = 0.0
               DO ispin = 1, input%jspins
                  DO ipm = 1, 2
                     CALL to_tmat(calc,g0%gmmpMat(1,iz,i_gf,:,:,ispin,ipm),1,l)
                     CALL calc%inverse()
                     delta%data_c = delta%data_c + 1/2.0 * (-1)**(ispin+ipm-2) * calc%data_c
                  ENDDO
               ENDDO
            ENDIF
            !
            !  Tr[\Delta (G_up-G-down) + \Delta G_up \Delta G-down]
            !
            DO ipm = 1, 2
               CALL to_tmat(calcup,g0%gmmpMat(1,iz,i_gf,:,:,1,ipm),1,l)
               CALL to_tmat(calcdwn,g0%gmmpMat(1,iz,i_gf,:,:,1,ipm),1,l)
               calcup%data_c = matmul(delta%data_c,calcup%data_c)
               calcdwn%data_c = matmul(delta%data_c,calcdwn%data_c)
               
               calc%data_c = calcup%data_c - calcdwn%data_c + matmul(calcup%data_c,calcdwn%data_c)

               !Calculate the trace
               integrand = 0.0
               DO i = 1,matsize
                  integrand = integrand + calc%data_c(i,i)
               ENDDO
               j0 = j0 + (-1)**(ipm-1) * integrand*(REAL(g0%de(iz)) + (-1)**(ipm-1) * ImagUnit * AIMAG(g0%de(iz)))
            ENDDO
         ENDDO

         
         j0 = -1/(2.0*fpi_const)*hartree_to_ev_const * AIMAG(j0)
         WRITE(*,*)  "Eff. Exchange Interaction for atom", n, ": ", j0, "eV"

         CALL calcup%free()
         CALL calcdwn%free()
         CALL calc%free()
         CALL delta%free()
      ENDDO
   END SUBROUTINE eff_excinteraction

   SUBROUTINE occmtx(g,i_gf,atoms,sym,jspins,beta,mmpMat)

      !calculates the occupation of a orbital treated with DFT+HIA from the related greens function
      !The Greens-function should already be prepared on a energy contour ending at e_fermi
      !The occupation is calculated with:
      !
      ! n^sigma_mm' = -1/2pi int^Ef dz (G^+(z)^sigma_mm'-G^-(z)^sigma_mm')

      IMPLICIT NONE

      TYPE(t_greensf),        INTENT(IN)  :: g
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      TYPE(t_sym),            INTENT(IN)  :: sym
      COMPLEX,                INTENT(OUT) :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
      INTEGER,                INTENT(IN)  :: i_gf
      REAL,                   INTENT(IN)  :: beta
      INTEGER,                INTENT(IN)  :: jspins

      INTEGER i, m,mp, l, ispin, n,ipm,iz
      LOGICAL l_vertcorr
      COMPLEX :: g_int(g%nz)

      l_vertcorr = .true. !Enables/Disables a correction for the vertical parts of the rectangular contour

      mmpMat(:,:,:) = CMPLX(0.0,0.0)
      l = atoms%onsiteGF(i_gf)%l
      n = atoms%onsiteGF(i_gf)%atomType 

      DO ispin = 1, jspins
         DO m = -l, l
            DO mp = -l, l
               DO ipm = 1, 2
                  !If necessary here the radial averaging is performed
                  CALL int_sph(g%gmmpMat(:,:,i_gf,m,mp,ispin,ipm),atoms,n,g%nz,g%nr(i_gf),g_int)
                  !APPROXIMATION FOR THE VERTICAL PARTS OF THE CONTOUR:
                  IF(g%mode.EQ.1.AND.l_vertcorr) THEN
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + ImagUnit * REAL(g_int(1)) * AIMAG(g%e(1))
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) - ImagUnit * REAL(g_int(g%nz-g%nmatsub)) * AIMAG(g%e(1))
                  ENDIF
                  !Integrate over the contour:
                  DO iz = 1, g%nz-g%nmatsub
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + (-1)**(ipm-1) * g_int(iz)*(REAL(g%de(iz)) + (-1)**(ipm-1) * ImagUnit * AIMAG(g%de(iz)))
                  ENDDO
                  !NOT WORKING: MATSUBARA FREQ
                  !iz = g%nz-g%nmatsub
                  !DO i = g%nmatsub, 1, -1
                  !   iz = iz + 1
                  !   mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + (-1)**(ipm-1) * g%gmmpMat(1,iz,i_gf,m,mp,ispin,ipm)*g%de(iz)
                  !ENDDO
               ENDDO
               mmpMat(m,mp,ispin) = -1/(2.0 * pi_const) * AIMAG(mmpMat(m,mp,ispin))
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE occmtx

   SUBROUTINE int_sph(g,atoms,nType,nz,nr,g_int)

      USE m_types
      USE m_intgr
      USE m_constants
      !Performs the radial averaging if the green's function is calculated with radial dependence

      TYPE(t_atoms), INTENT(IN)  :: atoms
      COMPLEX,       INTENT(IN)  :: g(nr,nz)
      COMPLEX,       INTENT(OUT) :: g_int(nz)
      INTEGER,       INTENT(IN)  :: nr 
      INTEGER,       INTENT(IN)  :: nz
      INTEGER,       INTENT(IN)  :: nType

      INTEGER iz 
      REAL    re,imag

      DO iz = 1, nz
         IF(nr.EQ.1) THEN
            !Green's function was already calculated with spherical average
            g_int(iz) = g(1,iz)
         ELSE
            CALL intgr3( REAL(g(:,iz)),atoms%rmsh(:,nType),atoms%dx(nType),atoms%jri(nType),re)
            CALL intgr3(AIMAG(g(:,iz)),atoms%rmsh(:,nType),atoms%dx(nType),atoms%jri(nType),imag)
            g_int(iz) = re + ImagUnit * imag
         ENDIF
      ENDDO

   END SUBROUTINE int_sph



   SUBROUTINE ldosmtx(app,g,i_gf,atoms,sym,jspins)

      !calculates the l-dos from the onsite green's function 
      !If mode = 2 this may not make sense

      USE m_intgr

      CHARACTER(len=*),       INTENT(IN)  :: app
      TYPE(t_greensf),        INTENT(IN)  :: g
      TYPE(t_atoms),          INTENT(IN)  :: atoms 
      TYPE(t_sym),            INTENT(IN)  :: sym 
      INTEGER,                INTENT(IN)  :: i_gf
      INTEGER,                INTENT(IN)  :: jspins
      REAL :: dos(g%nz,jspins,2)

      INTEGER ispin,m,j,l,n,i,ipm
      REAL    imag

      dos = 0.0

     
      !n(E) = 1/2pii * Tr(G^+-G^-)
      l = atoms%onsiteGF(i_gf)%l
      n = atoms%onsiteGF(i_gf)%atomType
      DO ispin = 1, jspins
         DO m = -l , l
            DO ipm = 1, 2
               DO j = 1, g%nz
                  IF(g%nr(i_gf).EQ.1) THEN
                     imag = AIMAG(g%gmmpMat(1,j,i_gf,m,m,ispin,ipm))
                  ELSE
                     CALL intgr3(AIMAG(g%gmmpMat(:,j,i_gf,m,m,ispin,ipm)),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imag)
                  END IF
                  dos(j,ispin,ipm) = dos(j,ispin,ipm) - 1/(pi_const) * imag
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      OPEN(1337,file="lDOS_up_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

      DO i = 1, g%nz
         WRITE(1337,*) REAL(g%e(i)), dos(i,1,1), dos(i,1,2)
      ENDDO

      CLOSE(unit = 1337)
      IF(jspins.EQ.2) THEN
         OPEN(1337,file="lDOS_dwn_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

         DO i = 1, g%nz
            WRITE(1337,*) REAL(g%e(i)), dos(i,2,1), dos(i,2,2)
         ENDDO

         CLOSE(unit = 1337)
      ENDIF

   END SUBROUTINE

   SUBROUTINE indexgf(atoms,l,n,ind)

      !Find the index of the greens function associated with this l,n 

      USE m_types

      !Finds the corresponding entry in gmmpMat for given atomType and l

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      INTEGER,             INTENT(IN)  :: l,n
      INTEGER,             INTENT(OUT) :: ind

      ind = 0
      DO 
         ind = ind + 1
         IF(atoms%onsiteGF(ind)%atomType.EQ.n.AND.atoms%onsiteGF(ind)%l.EQ.l) THEN
            EXIT
         ENDIF
         IF(ind.EQ.atoms%n_gf) CALL juDFT_error("Green's function element not found", hint="This is a bug in FLEUR, please report")
      ENDDO
   END SUBROUTINE

      !Conversion between TYPE(t_mat) and array in TYPE(t_greensf) to make use of inversion routines

   SUBROUTINE to_tmat(gmat,g,jspins,l)

      USE m_types
      USE m_constants
      !Writes the array consisting of (m,mp,ispin) into a 2d matrix of type(t_mat)

      TYPE(t_mat),   INTENT(INOUT)  :: gmat 
      COMPLEX,       INTENT(IN)     :: g(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
      INTEGER,       INTENT(IN)     :: jspins
      INTEGER,       INTENT(IN)     :: l

      INTEGER ns,ind1,ind2,i,j,spin

      ns = 2*l + 1

      gmat%data_c = 0.0
      DO spin = 1, jspins
         DO i = 1, ns
            DO j = 1, ns
               ind1 = i + (spin-1)*ns
               ind2 = j + (spin-1)*ns
               IF(jspins.EQ.1) THEN
                  gmat%data_c(ind1,ind2) = 0.5*g(i-l-1,j-l-1,1)
               ELSE
                  gmat%data_c(ind1,ind2) = g(i-l-1,j-l-1,spin)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   
   END SUBROUTINE to_tmat

   SUBROUTINE to_g(gmat,g,jspins,l)

      USE m_types
      USE m_constants
      !Writes the array consisting of (m,mp,ispin) into a 2d matrix of type(t_mat)

      TYPE(t_mat),   INTENT(IN)  :: gmat 
      COMPLEX,       INTENT(INOUT) :: g(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
      INTEGER,       INTENT(IN)  :: jspins
      INTEGER,       INTENT(IN)  :: l

      INTEGER ns,ind1,ind2,i,j,spin

      ns = 2*l + 1

      g = 0.0
      DO spin = 1, 2
         DO i = 1, ns
            DO j = 1, ns
               ind1 = i + (spin-1)*ns
               ind2 = j + (spin-1)*ns
               IF(jspins.EQ.1) THEN
                  g(i-l-1,j-l-1,1) = g(i-l-1,j-l-1,1) + gmat%data_c(ind1,ind2)
               ELSE
                  g(i-l-1,j-l-1,spin) = gmat%data_c(ind1,ind2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
   
   END SUBROUTINE to_g



END MODULE m_gfcalc