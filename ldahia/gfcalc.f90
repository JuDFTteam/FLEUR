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
   !>       -calculates the occupation matrix 
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

   SUBROUTINE eff_excinteraction(g0,atoms,input,j0,ef,onsite_exc_split)

      USE m_ExpSave
      TYPE(t_greensf),        INTENT(IN)  :: g0
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      REAL,                INTENT(OUT) :: j0
      REAL,                INTENT(IN)  :: ef
      TYPE(t_input),          INTENT(IN)  :: input
      REAL,                   INTENT(IN)  :: onsite_exc_split

      COMPLEX integrand
      INTEGER i,iz,m,l,mp,ispin,n,i_gf,matsize,ipm
      REAL beta
      LOGICAL l_matinv
      TYPE(t_mat) :: calcup,calcdwn
      TYPE(t_mat) :: delta
      TYPE(t_mat) :: calc

      beta = input%onsite_beta *hartree_to_ev_const
      l_matinv = .false. !Determines how the onsite exchange splitting is calculated
      !WRITE(*,*) onsite_exc_split

      DO i_gf = 1, atoms%n_gf
         j0 = 0.0
         l = atoms%onsiteGF(i_gf)%l
         !We want to calculate j0 from the d-bands
         IF(l.NE.2) CYCLE
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
         OPEN(unit=1337,file="j0",status="replace",action="write")
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
                     CALL to_tmat(calc,g0%gmmpMat(1,iz,i_gf,:,:,ispin,ipm),1,1,l)
                     CALL calc%inverse()
                     delta%data_c = delta%data_c + 1/2.0 * (-1)**(ispin-1) * calc%data_c
                  ENDDO
               ENDDO
            ENDIF
            !
            !  Tr[\Delta (G_up-G_down) + \Delta G_up \Delta G_down]
            !
            ! calculated for G^+/- and then substract to obtain imaginary part
            integrand = 0.0
            DO ipm = 1, 2
               CALL to_tmat(calcup,g0%gmmpMat(1,iz,i_gf,:,:,1,ipm),1,1,l)
               CALL to_tmat(calcdwn,g0%gmmpMat(1,iz,i_gf,:,:,2,ipm),1,1,l)
               
               calcup%data_c  = matmul(delta%data_c,calcup%data_c)
               calcdwn%data_c = matmul(delta%data_c,calcdwn%data_c)

               calc%data_c = calcup%data_c - calcdwn%data_c + matmul(calcup%data_c,calcdwn%data_c)

               !Calculate the trace
               DO i = 1,matsize
                  integrand = integrand + (-1)**(ipm-1) * calc%data_c(i,i)
               ENDDO
            ENDDO
            WRITE(1337,"(2f14.8)") REAL(g0%e(iz)), AIMAG(integrand/(1.0+exp_save(beta*(real(g0%e(iz))-ef))))
            j0 = j0 + AIMAG(integrand*(REAL(g0%de(iz)))/(1.0+exp_save(beta*(real(g0%e(iz))-ef))))
         ENDDO
         CLOSE(1337)
         j0 = -1/(2.0*fpi_const)*hartree_to_ev_const * j0
         WRITE(*,*)  "Eff. Exchange Interaction for atom", n, ": ", j0, "eV"

         CALL calcup%free()
         CALL calcdwn%free()
         CALL calc%free()
         CALL delta%free()
      ENDDO
   END SUBROUTINE eff_excinteraction

   SUBROUTINE occmtx(g,i_gf,atoms,sym,jspins,beta,ef,mmpMat)


      USE m_ExpSave
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
      REAL,                   INTENT(IN)  :: ef
      INTEGER,                INTENT(IN)  :: jspins

      INTEGER i, m,mp, l, ispin, n,ipm,iz
      LOGICAL l_vertcorr
      COMPLEX :: g_int(g%nz)

      l_vertcorr = .false. !Enables/Disables a correction for the vertical parts of the rectangular contour

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
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + (-1)**(ipm-1) * ImagUnit * REAL(g_int(1)) * AIMAG(g%e(1))
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) - (-1)**(ipm-1) * ImagUnit * REAL(g_int(g%nz-g%nmatsub)) * AIMAG(g%e(1))
                  ENDIF
                  !Integrate over the contour:
                  DO iz = 1, g%nz-g%nmatsub
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + (-1)**(ipm-1) * AIMAG(g_int(iz))*REAL(g%de(iz))/(1.0+exp_save(beta*real(g%e(iz)-ef)))
                  ENDDO
                  !NOT WORKING: MATSUBARA FREQ
                  !iz = g%nz-g%nmatsub
                  !DO i = g%nmatsub, 1, -1
                  !   iz = iz + 1
                  !   mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) + (-1)**(ipm-1) * g%gmmpMat(1,iz,i_gf,m,mp,ispin,ipm)*g%de(iz)
                  !ENDDO
               ENDDO
               mmpMat(m,mp,ispin) = -1/(2.0 * pi_const) * mmpMat(m,mp,ispin)
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
      REAL :: dos(-lmaxU_const:lmaxU_const,g%nz,jspins)
      REAL :: re(-lmaxU_const:lmaxU_const,g%nz,jspins)
      INTEGER ispin,m,j,l,n,i,ipm
      REAL    imag

      dos = 0.0
      re = 0.0

     
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
                  dos(m,j,ispin) = dos(m,j,ispin) - (-1)**(ipm-1)  * imag
                  re(m,j,ispin) = re(m,j,ispin)  +  REAL(g%gmmpMat(1,j,i_gf,m,m,ispin,ipm))
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      OPEN(1337,file="lDOS_up_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

      DO i = 1, g%nz
         WRITE(1337,"(9f14.8)") REAL(g%e(i)), SUM(dos(-l:l,i,1)), (dos(m,i,1), m = -lmaxU_const, lmaxU_const)
      ENDDO

      CLOSE(unit = 1337)
      IF(jspins.EQ.2) THEN
         OPEN(1337,file="lDOS_dwn_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

         DO i = 1, g%nz
            WRITE(1337,"(9f14.8)") REAL(g%e(i)),SUM(dos(-l:l,i,2)), (dos(m,i,2), m = -lmaxU_const, lmaxU_const)
         ENDDO

         CLOSE(unit = 1337)
      ENDIF

      OPEN(1337,file="Re_up_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

      DO i = 1, g%nz
         WRITE(1337,"(9f15.8)") REAL(g%e(i)), SUM(re(-l:l,i,1)), (re(m,i,1), m = -lmaxU_const, lmaxU_const)
      ENDDO

      CLOSE(unit = 1337)
      IF(jspins.EQ.2) THEN
         OPEN(1337,file="Re_dwn_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

         DO i = 1, g%nz
            WRITE(1337,"(9f15.8)") REAL(g%e(i)),SUM(re(-l:l,i,2)), (re(m,i,2), m = -lmaxU_const, lmaxU_const)
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

   SUBROUTINE to_tmat(gmat,g,jsp_in,jsp_out,l)

      USE m_types
      USE m_constants
      !Writes the array consisting of (m,mp,ispin) into a 2d matrix of type(t_mat)

      TYPE(t_mat),   INTENT(INOUT)  :: gmat 
      COMPLEX,       INTENT(IN)     :: g(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jsp_in)
      INTEGER,       INTENT(IN)     :: jsp_in
      INTEGER,       INTENT(IN)     :: jsp_out
      INTEGER,       INTENT(IN)     :: l

      INTEGER ns,ind1,ind2,i,j,spin

      ns = 2*l + 1

      gmat%data_c = 0.0
      DO spin = 1, jsp_out
         DO i = 1, ns
            DO j = 1, ns
               ind1 = i + (spin-1)*ns
               ind2 = j + (spin-1)*ns
               IF(jsp_in.EQ.1.AND.jsp_out.EQ.2) THEN
                  gmat%data_c(ind1,ind2) = 0.5*g(i-l-1,j-l-1,1)
               ELSE
                  gmat%data_c(ind1,ind2) = g(i-l-1,j-l-1,spin)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   
   END SUBROUTINE to_tmat

   SUBROUTINE to_g(gmat,g,jsp_in,jsp_out,l)

      USE m_types
      USE m_constants
      !Writes the array consisting of (m,mp,ispin) into a 2d matrix of type(t_mat)

      TYPE(t_mat),   INTENT(IN)  :: gmat 
      COMPLEX,       INTENT(INOUT) :: g(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jsp_out)
      INTEGER,       INTENT(IN)     :: jsp_in
      INTEGER,       INTENT(IN)     :: jsp_out
      INTEGER,       INTENT(IN)  :: l

      INTEGER ns,ind1,ind2,i,j,spin

      ns = 2*l + 1

      g = 0.0
      DO spin = 1, jsp_in
         DO i = 1, ns
            DO j = 1, ns
               ind1 = i + (spin-1)*ns
               ind2 = j + (spin-1)*ns
               IF(jsp_in.EQ.2.AND.jsp_out.EQ.1) THEN
                  g(i-l-1,j-l-1,1) = g(i-l-1,j-l-1,1) + gmat%data_c(ind1,ind2)
               ELSE
                  g(i-l-1,j-l-1,spin) = gmat%data_c(ind1,ind2)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      
   
   END SUBROUTINE to_g

   SUBROUTINE local_sym(mat,atoms,sym,l,n)

      USE m_types
      !Make sure that the elements that should be zero by symmetry are actually zero
   
      COMPLEX, INTENT(INOUT) :: mat(-l:l,-l:l)
      TYPE(t_sym), INTENT(IN) :: sym 
      TYPE(t_atoms),INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: l
      INTEGER, INTENT(IN) :: n
   
      INTEGER natom,it,is,isi,m,mp,n_op
      COMPLEX mat_tmp(-l:l,-l:l)
      COMPLEX mat_op(-l:l,-l:l),d_tmp(-l:l,-l:l)
      COMPLEX sym_diag(-l:l)
   
      !Store the original matrix
      mat_tmp(-l:l,-l:l) = mat(-l:l,-l:l)
      mat = 0.0
      n_op =0
      !Is the loop over equivalent atoms necessary ??
      DO natom = SUM(atoms%neq(:n-1)) + 1, SUM(atoms%neq(:n))
         IF(sym%invarind(natom).EQ.0) CALL juDFT_error("No local symmetries found for projected DOS",&
                                       hint="This is a bug in FLEUR, please report",calledby="local_sym")
         DO it = 1, sym%invarind(natom)
            is = sym%invarop(natom,it)
            isi = sym%invtab(is)
            d_tmp(:,:) = cmplx(0.0,0.0)
            DO m = -l,l
               DO mp = -l,l
                  d_tmp(m,mp) = sym%d_wgn(m,mp,l,isi)
               ENDDO
            ENDDO
            DO m = -l,l
               sym_diag(m) = d_tmp(m,m)
               d_tmp(m,m) = 0.0
            ENDDO
            !Exclude all symmetries that would prevent splitting of the levels
            !IF(ANY(d_tmp(:,:).NE.0.0)) CYCLE
            n_op = n_op + 1
            DO m = -l,l
               d_tmp(m,m) = sym_diag(m)
            ENDDO
            mat_op = matmul( transpose( conjg(d_tmp) ) , mat_tmp)
            mat_op =  matmul( mat_op, d_tmp )
            DO m = -l,l
               DO mp = -l,l
                  mat(m,mp) = mat(m,mp) + conjg(mat_op(m,mp))
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   
      mat = mat * 1.0/(REAL(n_op))
   
   END SUBROUTINE local_sym


END MODULE m_gfcalc