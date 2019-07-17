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

   SUBROUTINE occmtx(g,l,nType,atoms,sym,input,mmpMat,lp,nTypep)

      USE m_ind_greensf
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
      COMPLEX,                INTENT(OUT) :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MERGE(3,input%jspins,input%l_gfmperp))
      INTEGER,                INTENT(IN)  :: l
      INTEGER,                INTENT(IN)  :: nType
      INTEGER, OPTIONAL,      INTENT(IN)  :: lp
      INTEGER, OPTIONAL,      INTENT(IN)  :: nTypep

      INTEGER ind1,ind2,ipm,iz,ispin,m,mp,lp_loop
      LOGICAL l_vertcorr
      REAL    re,imag
      TYPE(t_mat) :: gmat

      l_vertcorr = .false. !Enables/Disables a correction for the vertical parts of the rectangular contour

      mmpMat(:,:,:) = CMPLX(0.0,0.0)

      IF(.NOT.PRESENT(lp)) THEN
         lp_loop = l 
      ELSE 
         lp_loop = lp 
      ENDIF

      !REPLACE: input%jspins --> MERGE(3,input%jspins,input%l_gfmperp)
      DO ispin = 1, MERGE(3,input%jspins,input%l_gfmperp)
         DO ipm = 1, 2
            !Integrate over the contour:
            DO iz = 1, g%nz
               !get the corresponding gf-matrix
               CALL g%get_gf(gmat,atoms,input,iz,l,nType,ipm.EQ.2,spin=ispin,lp=lp,nTypep=nTypep)
               ind1 = 0
               DO m = -l, l
                  ind1 = ind1 + 1
                  ind2 = 0 
                  DO mp = -lp_loop,lp_loop
                     ind2 = ind2 + 1 
                     mmpMat(m,mp,ispin) = mmpMat(m,mp,ispin) - 1/(2.0*pi_const*ImagUnit) * (-1)**(ipm-1) * gmat%data_c(ind1,ind2) &
                                                             * MERGE(g%de(iz),conjg(g%de(iz)),ipm.EQ.1)
                  ENDDO
               ENDDO
               CALL gmat%free()
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE occmtx

   SUBROUTINE crystal_field(atoms,input,greensfCoeffs,hub1,vu)

      !calculates the crystal-field matrix for the local hamiltonian
      !In addition we calculate the onsite exchange splitting for the calculation of j0


      USE m_kkintgr
      USE m_ind_greensf

      IMPLICIT NONE

      !-Type Arguments
      TYPE(t_greensfCoeffs), INTENT(IN)    :: greensfCoeffs
      TYPE(t_atoms),         INTENT(IN)    :: atoms
      TYPE(t_input),         INTENT(IN)    :: input
      TYPE(t_hub1ham),       INTENT(INOUT) :: hub1

      !-Array Arguments
      COMPLEX,               INTENT(IN)    :: vu(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins) !LDA+U potential (should be removed from h_loc)

      !-Local Scalars
      INTEGER i_gf,l,nType,jspin,m,mp,ie,i_hia,kkcut,spin_cut

      !-Local Arrays
      REAL :: h_loc(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_hia,input%jspins)
      REAL :: integrand(greensfCoeffs%ne), norm(input%jspins),tr_hloc(atoms%n_hia,input%jspins)

      h_loc = 0.0
      tr_hloc = 0.0
      DO i_hia = 1, atoms%n_hia

         l     = atoms%lda_u(atoms%n_u+i_hia)%l
         nType = atoms%lda_u(atoms%n_u+i_hia)%atomType
         i_gf = ind_greensf(atoms,l,nType)
         !Perform the integration 
         !
         ! \int_{E_b}^{E_c} dE E * N_LL'(E)
         !
         norm = 0.0
         DO jspin = 1, input%jspins
            spin_cut = MERGE(1,jspin,jspin.GT.2)
            kkcut = greensfCoeffs%kkintgr_cutoff(i_gf,spin_cut,2)
            DO m = -l, l
               DO mp = -l, l
                  integrand = 0.0
                  DO ie = 1, kkcut
                     integrand(ie) = ((ie-1) * greensfCoeffs%del+greensfCoeffs%e_bot) * greensfCoeffs%projdos(ie,i_gf,m,mp,jspin)
                  ENDDO
                  h_loc(m,mp,i_hia,jspin) = trapz(integrand,greensfCoeffs%del,greensfCoeffs%ne)
               ENDDO
               !trace of the integrated E*projdos
               tr_hloc(i_hia,jspin) = tr_hloc(i_hia,jspin) + h_loc(m,m,i_hia,jspin) - REAL(vu(m,m,i_hia,jspin))
            ENDDO
         ENDDO

         !Average over spins
         hub1%ccfmat(i_hia,:,:) = 0.0
         DO m = -l, l
            DO mp = -l, l
               hub1%ccfmat(i_hia,m,mp) = SUM(h_loc(m,mp,i_hia,:))/REAL(input%jspins) &
                                       - SUM(vu(m,mp,i_hia,:))/REAL(input%jspins)
            ENDDO
            hub1%ccfmat(i_hia,m,m) = hub1%ccfmat(i_hia,m,m) - SUM(tr_hloc(i_hia,:))/REAL(input%jspins*(2*l+1))
         ENDDO
      ENDDO


   END SUBROUTINE crystal_field


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
      l = atoms%gfelem(i_gf)%l
      n = atoms%gfelem(i_gf)%atomType
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