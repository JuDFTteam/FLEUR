MODULE m_cfmat

   !------------------------------------------------------------------------------
   !
   ! MODULE:  m_cfmat
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION: 
   !>  calculates the crytal field contribution to the atomic hamiltonian
   !>  in DFT+Hubbard 1
   !
   ! REVISION HISTORY:
   ! March 2019 - Initial Version
   !------------------------------------------------------------------------------

   USE m_juDFT
   USE m_types
   USE m_gfcalc
   USE m_constants
   USE m_kkintgr

   IMPLICIT NONE

   CONTAINS 

   SUBROUTINE cfcontrib(projDOS,l,n,nb,nc,ne,del,jspins,cfmat,pot)

      REAL,                INTENT(IN)  :: projDOS(ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
      INTEGER,             INTENT(IN)  :: l
      INTEGER,             INTENT(IN)  :: n
      INTEGER,             INTENT(IN)  :: nb
      INTEGER,             INTENT(IN)  :: nc
      INTEGER,             INTENT(IN)  :: ne
      REAL,                INTENT(IN)  :: del 
      INTEGER,             INTENT(IN)  :: jspins
      REAL,                INTENT(OUT) :: cfmat(-l:l,-l:l)
      !Stuff to be removed from the local hamiltonian to obtain cfmat
      REAL,                INTENT(IN)  :: pot(-l:l,-l:l,jspins)
      
      REAL :: cfmat_sp(-l:l,-l:l,jspins)
      INTEGER ispin,m,mp,i
      REAL trac

      REAL,ALLOCATABLE :: integrand(:)
      !Perform the integration 
      !
      ! \int_{E_b}^{E_c} dE E * N_LL'(E)
      !
      ALLOCATE(integrand(nc-nb))
      integrand = 0.0
      cfmat_sp = 0.0
      DO ispin = 1, jspins
         DO m = -l, l
            DO mp = -l, l
               DO i = 1, nc-nb
                  integrand(i) = (i-1+nb) * del * projDOS(i+nb,m,mp,ispin)
               ENDDO
               CALL trapz(integrand,del,nc-nb,cfmat_sp(m,mp,ispin))
                !Remove the LDA+U potential
               cfmat_sp(m,mp,ispin) = cfmat_sp(m,mp,ispin) - pot(m,mp,ispin) 
            ENDDO
         ENDDO
      ENDDO
      

      !Average over spins
      cfmat = 0.0
      DO m = -l, l
         DO mp = -l, l
            cfmat(m,mp) = SUM(cfmat_sp(m,mp,:))/REAL(jspins)
         ENDDO
      ENDDO
      !Make the matrix traceless
      trac = 0.0
      DO m = -l, l
         trac = trac + cfmat(m,m)
      ENDDO
      !WRITE(*,*) trac/(2*l+1)
      DO m = -l, l
         cfmat(m,m) = cfmat(m,m) - trac/(2*l+1)
      ENDDO


   END SUBROUTINE cfcontrib



END MODULE m_cfmat