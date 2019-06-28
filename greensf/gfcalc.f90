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
   !>       -calculates the crystal-field-contrbution for the local hamiltonian
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

   SUBROUTINE occmtx(g,i_gf,atoms,sym,input,ef,mmpMat,el0,vr)


      USE m_ExpSave
      USE m_radfun
      USE m_intgr
      !calculates the occupation of a orbital treated with DFT+HIA from the related greens function
      !The Greens-function should already be prepared on a energy contour ending at e_fermi
      !The occupation is calculated with:
      !
      ! n^sigma_mm' = -1/2pi int^Ef dz (G^+(z)^sigma_mm'-G^-(z)^sigma_mm')

      IMPLICIT NONE

      TYPE(t_greensf),        INTENT(IN)  :: g
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      TYPE(t_sym),            INTENT(IN)  :: sym
      TYPE(t_input),          INTENT(IN)  :: input
      COMPLEX,                INTENT(OUT) :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,input%jspins)
      INTEGER,                INTENT(IN)  :: i_gf
      REAL,                   INTENT(IN)  :: ef
      REAL, OPTIONAL,         INTENT(IN)  :: el0(input%jspins)  !energy parameter for radial functions
      REAL, OPTIONAL,         INTENT(IN)  :: vr(atoms%jmtd,input%jspins)

      INTEGER i, m,mp, l, ispin, nType,ipm,iz
      INTEGER nodeu,noded
      LOGICAL l_vertcorr,l_hia
      COMPLEX g_int,nmmp,weight
      REAL    wronk,re,imag,beta
      TYPE(t_usdus) :: usdus

      REAL vrav(atoms%jmtd)
      REAL, ALLOCATABLE :: u(:,:),udot(:,:)
      COMPLEX, ALLOCATABLE :: gr(:)

      IF(.NOT.input%l_gfsphavg.AND.(.NOT.PRESENT(el0).OR..NOT.PRESENT(vr))) THEN
         CALL juDFT_error("Cannot calculate radial dependence for green's function", calledby="occmtx")
      ENDIF
      IF(.NOT.input%l_gfsphavg) THEN
         ALLOCATE(u(atoms%jmtd,2))
         ALLOCATE(udot(atoms%jmtd,2))
         ALLOCATE(gr(atoms%jmtd))
         u = 0.0
         udot = 0.0
         gr = 0.0
      ENDIF 
      l_vertcorr = .false. !Enables/Disables a correction for the vertical parts of the rectangular contour

      mmpMat(:,:,:) = CMPLX(0.0,0.0)
      l = atoms%onsiteGF(i_gf)%l
      nType = atoms%onsiteGF(i_gf)%atomType 


      DO ispin = 1, input%jspins
         IF(.NOT.input%l_gfsphavg) THEN
            !If we have the radial dependence of the greens function calculate the radial functions here
            !Check if the orbital is to be treated with Hubbard 1
            l_hia=.FALSE.
            DO i = 1, atoms%n_hia
               IF(atoms%lda_hia(i)%atomType.EQ.nType.AND.atoms%lda_hia(i)%l.EQ.l) THEN
                  l_hia=.TRUE.
               ENDIF
            ENDDO
            IF(l_hia.AND.input%jspins.EQ.2) THEN
               !In the case of a spin-polarized calculation with Hubbard 1 we want to treat 
               !the correlated orbitals with a non-spin-polarized basis
              
               vrav(1:atoms%jri(nType)) = (vr(1:atoms%jri(nType),1) + vr(1:atoms%jri(nType),2))/2.0
               CALL radfun(l,nType,ispin,el0(ispin),vrav,atoms,u,udot,usdus, nodeu,noded,wronk)  
            ELSE
               CALL radfun(l,nType,ispin,el0(ispin),vr(:,ispin),atoms,u,udot,usdus,nodeu,noded,wronk)
            ENDIF
         ENDIF
         DO m = -l, l
            DO mp = -l, l
               nmmp = 0.0
               DO ipm = 1, 2
                  !Integrate over the contour:
                  DO iz = 1, g%nz
                     !weight for the energy integration
                     !If necessary here the radial averaging is performed
                     IF(input%l_gfsphavg) THEN
                        g_int = g%gmmpMat(iz,i_gf,m,mp,ispin,ipm)
                     ELSE
                        gr = onsite_radial(g%uu(iz,i_gf,m,mp,ispin,ipm),g%dd(iz,i_gf,m,mp,ispin,ipm),g%du(iz,i_gf,m,mp,ispin,ipm),&
                                           g%ud(iz,i_gf,m,mp,ispin,ipm),u,udot,atoms%jmtd)
                        CALL intgr3(REAL(gr),atoms%rmsh(:,nType),atoms%dx(nType),atoms%jri(nType),re)
                        CALL intgr3(AIMAG(gr),atoms%rmsh(:,nType),atoms%dx(nType),atoms%jri(nType),imag)
                        g_int = re + ImagUnit * imag
                     ENDIF
                     IF((iz.EQ.1.OR.iz.EQ.g%nz-g%nmatsub).AND.(g%mode.EQ.1.AND.l_vertcorr)) THEN
                        !APPROXIMATION FOR THE VERTICAL PARTS OF THE CONTOUR:
                        nmmp = nmmp + (-1)**(ipm-1) * ImagUnit * REAL(g_int) * AIMAG(g%e(1))
                     ENDIF
                     nmmp = nmmp + (-1)**(ipm-1) * g_int * MERGE(g%de(iz),conjg(g%de(iz)),ipm.EQ.1)
                  ENDDO
               ENDDO
               mmpMat(m,mp,ispin) = -1/(2.0 * pi_const) * AIMAG(nmmp)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE occmtx


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

   SUBROUTINE crystal_field(atoms,input,greensfCoeffs,hub1,vu)

      !calculates the crystal-field matrix for the local hamiltonian
      !In addition we calculate the onsite exchange splitting for the calculation of j0


      USE m_kkintgr

      IMPLICIT NONE

      !-Type Arguments
      TYPE(t_greensfCoeffs), INTENT(IN)    :: greensfCoeffs
      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_input),         INTENT(IN)    :: input
      TYPE(t_hub1ham),       INTENT(INOUT) :: hub1

      !-Array Arguments
      COMPLEX,               INTENT(IN)    :: vu(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins) !LDA+U potential (should be removed from h_loc)

      !-Local Scalars
      INTEGER i_gf,l,nType,jspin,m,mp,ie,i_hia

      !-Local Arrays
      REAL :: h_loc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_gf,input%jspins)
      REAL :: integrand(greensfCoeffs%ne), norm(input%jspins),tr_hloc(atoms%n_gf,input%jspins)

      h_loc = 0.0
      tr_hloc = 0.0
      DO i_hia = 1, atoms%n_hia

         l     = atoms%lda_hia(i_hia)%l
         nType = atoms%lda_hia(i_hia)%atomType
         CALL indexgf(atoms,l,nType,i_gf)
         !Perform the integration 
         !
         ! \int_{E_b}^{E_c} dE E * N_LL'(E)
         !
         norm = 0.0
         DO jspin = 1, input%jspins
            DO m = -l, l
               DO mp = -l, l
                  integrand = 0.0
                  DO ie = greensfCoeffs%kkintgr_cutoff(i_gf,jspin,1), greensfCoeffs%kkintgr_cutoff(i_gf,jspin,2)
                     integrand(ie) = ((ie-1) * greensfCoeffs%del+greensfCoeffs%e_bot) * greensfCoeffs%projdos(ie,i_gf,m,mp,jspin)
                  ENDDO
                  h_loc(m,mp,i_gf,jspin) = trapz(integrand,greensfCoeffs%del,greensfCoeffs%ne)
               ENDDO
               !trace of the integrated E*projdos
               tr_hloc(i_gf,jspin) = tr_hloc(i_gf,jspin) + h_loc(m,m,i_gf,jspin) - REAL(vu(m,m,i_hia,jspin))
            ENDDO
         ENDDO

         !Average over spins
         hub1%ccfmat(i_hia,:,:) = 0.0
         DO m = -l, l
            DO mp = -l, l
               hub1%ccfmat(i_hia,m,mp) = SUM(h_loc(m,mp,i_gf,:))/REAL(input%jspins) &
                                       - SUM(vu(m,mp,i_hia,:))/REAL(input%jspins)
            ENDDO
            hub1%ccfmat(i_hia,m,m) = hub1%ccfmat(i_hia,m,m) - SUM(tr_hloc(i_gf,:))/REAL(input%jspins*(2*l+1))
         ENDDO
      ENDDO


   END SUBROUTINE crystal_field

   FUNCTION onsite_radial(uu,dd,du,ud,u,udot,nr)

      !Return the radial dependence of the onsite green's function for one specific matrix elements
   
      IMPLICIT NONE

      COMPLEX,        INTENT(IN) :: uu,dd,du,ud 
      INTEGER,        INTENT(IN) :: nr
      REAL, OPTIONAL, INTENT(IN) :: u(nr,2),udot(nr,2)
      COMPLEX :: onsite_radial(nr)

      INTEGER jr 

      onsite_radial = 0.0
      DO jr = 1, nr
         onsite_radial(jr) = uu * (u(jr,1)*u(jr,1)+u(jr,2)*u(jr,2)) + &
                             dd * (udot(jr,1)*udot(jr,1)+udot(jr,2)*udot(jr,2))+&
                             ud * (u(jr,1)*udot(jr,1)+u(jr,2)*udot(jr,2))+&
                             du * (udot(jr,1)*u(jr,1)+udot(jr,2)*u(jr,2))
      ENDDO
   
   END FUNCTION onsite_radial


      SUBROUTINE ldosmtx(app,g,i_gf,atoms,sym,input)

         !calculates the l-dos from the onsite green's function 
         !If mode = 2 this may not make sense

      USE m_intgr
      USE m_ExpSave
      CHARACTER(len=*),       INTENT(IN)  :: app
      TYPE(t_greensf),        INTENT(IN)  :: g
      TYPE(t_atoms),          INTENT(IN)  :: atoms 
      TYPE(t_sym),            INTENT(IN)  :: sym 
      TYPE(t_input),          INTENT(IN)  :: input
      INTEGER,                INTENT(IN)  :: i_gf
      REAL :: dos(-lmaxU_const:lmaxU_const,g%nz,input%jspins)
      REAL :: re(-lmaxU_const:lmaxU_const,g%nz,input%jspins)
      INTEGER ispin,m,j,l,n,i,ipm
      REAL    imag

      dos = 0.0
      re = 0.0

      IF(.NOT.input%l_gfsphavg) CALL juDFT_warn("Not implemented for radial dependence",calledby="ldosmtx")
      !n(E) = 1/2pii * Tr(G^+-G^-)
      l = atoms%onsiteGF(i_gf)%l
      n = atoms%onsiteGF(i_gf)%atomType
      DO ispin = 1, input%jspins
         DO m = -l , l
            DO ipm = 1, 2
               DO j = 1, g%nz

                  imag = AIMAG(g%gmmpMat(j,i_gf,m,m,ispin,ipm))
                  dos(m,j,ispin) = dos(m,j,ispin) + (-1)**(ipm-1)  * imag
                  re(m,j,ispin) = re(m,j,ispin)  +  REAL(g%gmmpMat(j,i_gf,m,m,ispin,ipm))
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      OPEN(1337,file="lDOS_up_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

      DO i = 1, g%nz
         WRITE(1337,"(9f15.8)") REAL(g%e(i)), SUM(dos(-l:l,i,1)), (dos(m,i,1), m = -lmaxU_const, lmaxU_const)
      ENDDO

      CLOSE(unit = 1337)
      IF(input%jspins.EQ.2) THEN
         OPEN(1337,file="lDOS_dwn_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

         DO i = 1, g%nz
            WRITE(1337,"(9f15.8)") REAL(g%e(i)),SUM(dos(-l:l,i,2)), (dos(m,i,2), m = -lmaxU_const, lmaxU_const)
         ENDDO

         CLOSE(unit = 1337)
      ENDIF

      OPEN(1337,file="Re_up_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

      DO i = 1, g%nz
         WRITE(1337,"(9f15.8)") REAL(g%e(i)), SUM(re(-l:l,i,1)), (re(m,i,1), m = -lmaxU_const, lmaxU_const)
      ENDDO

      CLOSE(unit = 1337)
      IF(input%jspins.EQ.2) THEN
         OPEN(1337,file="Re_dwn_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

         DO i = 1, g%nz
            WRITE(1337,"(9f15.8)") REAL(g%e(i)),SUM(re(-l:l,i,2)), (re(m,i,2), m = -lmaxU_const, lmaxU_const)
         ENDDO

         CLOSE(unit = 1337)
      ENDIF

   END SUBROUTINE

END MODULE m_gfcalc