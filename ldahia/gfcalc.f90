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

   SUBROUTINE eff_excinteraction(gOnsite,atoms,input,j0,onsite_exc_split)

      TYPE(t_greensf),        INTENT(IN)  :: gOnsite
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      REAL,                   INTENT(OUT) :: j0
      TYPE(t_input),          INTENT(IN)  :: input
      REAL,                   INTENT(IN)  :: onsite_exc_split

      COMPLEX tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),integrand(2)
      INTEGER iz,m,l,mp,ispin,n,i_gf,matsize,i
      LOGICAL l_matinv



      COMPLEX :: delta(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX :: g_up(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX :: g_dwn(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      INTEGER              :: info
      INTEGER, ALLOCATABLE :: ipiv(:)
      COMPLEX, ALLOCATABLE :: work(:)
      COMPLEX, ALLOCATABLE :: invup(:,:)
      COMPLEX, ALLOCATABLE :: invdwn(:,:)
      COMPLEX, ALLOCATABLE :: inv(:,:)

      l_matinv = .true. !Determines how the onsite exchange splitting is calculated
      
      DO i_gf = 1, atoms%n_gf
         j0 = 0.0
         l = atoms%onsiteGF(i_gf)%l
         n = atoms%onsiteGF(i_gf)%atomType
         IF(l_matinv) THEN
            matsize = 2*l+1
            ALLOCATE(work(matsize))
            ALLOCATE(invup(matsize,matsize))
            ALLOCATE(invdwn(matsize,matsize))
            ALLOCATE(inv(matsize,matsize))
            ALLOCATE(ipiv(matsize))
         ENDIF 

         DO iz = 1, gOnsite%nef
            !
            !calculate the onsite exchange matrix
            !
            IF(l_matinv) THEN
               !First Way: Matrix Inversion
               !---------------------------------------------
               !\Delta = (G_up)^-1-(G_down)^-1
               !---------------------------------------------
               !Symmetrize the green's function for up/down 
               !spin with respect to the complex plane
               !Here we assume that the onsite Hamiltonian
               !is real
               !---------------------------------------------
               !G^(up/down)^-1 = 1/2 * (G+^(up/down) + G-^(up/down))
               !---------------------------------------------
               !TODO: Replace Matrix inversion with calls to tmat%inverse
               inv(1:matsize,1:matsize) = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,1,1)
               CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               invup(1:matsize,1:matsize) = inv(1:matsize,1:matsize)
               inv(1:matsize,1:matsize) = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,1,2)
               CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               invup(1:matsize,1:matsize) = 1/2.0*(invup(1:matsize,1:matsize) +inv(1:matsize,1:matsize))


               inv(1:matsize,1:matsize) = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,2,1)
               CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               invdwn(1:matsize,1:matsize) = inv(1:matsize,1:matsize)
               inv(1:matsize,1:matsize) = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,2,2)
               CALL zgetrf(matsize,matsize,inv,matsize,ipiv,info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               CALL zgetri(matsize,inv,matsize,ipiv,work,size(work),info)
               IF(info.NE.0) CALL judft_error("Failed to invert matrix: dpotrf failed.",calledby="eff_excinteraction")
               invdwn(1:matsize,1:matsize) = 1/2.0*(invdwn(1:matsize,1:matsize) +inv(1:matsize,1:matsize))

               delta(-l:l,-l:l) = invup(1:matsize,1:matsize) - invdwn(1:matsize,1:matsize)
            ELSE
               !Second Way: onsite_exc_split is the difference in the center of gravity of the up/down bands
               delta = 0.0
               DO m = -l, l
                  delta(m,m) = onsite_exc_split
               ENDDO
            ENDIF
            !
            !  Tr[\Delta (G_up-G-down) + \Delta G_up \Delta G-down]
            !
            integrand = 0.0
            DO i = 1, 2
               g_up(-l:l,-l:l)   = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,1,i)
               g_dwn(-l:l,-l:l)  = gOnsite%gmmpMat(1,iz,i_gf,-l:l,-l:l,2,i)
               tmp(-l:l,-l:l)    = g_up(-l:l,-l:l)-g_dwn(-l:l,-l:l)
               tmp(-l:l,-l:l)    = matmul(delta(-l:l,-l:l),tmp(-l:l,-l:l))
               tmp(-l:l,-l:l)    = tmp(-l:l,-l:l) + matmul(matmul(delta(-l:l,-l:l),g_up(-l:l,-l:l)),&
                                 matmul(delta(-l:l,-l:l),g_dwn(-l:l,-l:l)))
               !Trace over m
               DO m = -l,l
                  integrand(i) = integrand(i) + tmp(m,m)
               ENDDO
            ENDDO
            !WRITE(*,*) gOnsite%e(iz), gOnsite%de(iz),  1/2.0 * AIMAG(integrand(1)*gOnsite%de(iz)-integrand(2)*conjg(gOnsite%de(iz)))
            WRITE(*,*) j0
            j0 = j0 + 1/2.0 * integrand(1)*gOnsite%de(iz)-integrand(2)*conjg(gOnsite%de(iz))
         ENDDO

         
         j0 = -j0*1/fpi_const*hartree_to_ev_const
         WRITE(*,*)  "Eff. Exchange Interaction for atom", n, ": ", j0, "eV"
         IF(ALLOCATED(work)) DEALLOCATE(work)
         IF(ALLOCATED(ipiv)) DEALLOCATE(ipiv)
         IF(ALLOCATED(invup)) DEALLOCATE(invup)
         IF(ALLOCATED(invdwn)) DEALLOCATE(invdwn)
      ENDDO
   END SUBROUTINE eff_excinteraction

   SUBROUTINE occmtx(g,i_gf,atoms,sym,jspins,mmpMat)

      USE m_intgr

      !calculates the occupation of a orbital treated with DFT+HIA from the related greens function
      !The Greens-function should already be prepared on a energy contour ending at e_fermi

      IMPLICIT NONE

      TYPE(t_greensf),        INTENT(IN)  :: g
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      TYPE(t_sym),            INTENT(IN)  :: sym
      COMPLEX,                INTENT(OUT) :: mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
      INTEGER,                INTENT(IN)  :: i_gf
      INTEGER,                INTENT(IN)  :: jspins

      INTEGER i, m,mp, l, ispin, n, it,is, isi, natom, nn
      REAL imag, re, fac
      COMPLEX n_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),nr_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX n1_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const), d_tmp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)


      mmpMat(:,:,:) = CMPLX(0.0,0.0)
      l = atoms%onsiteGF(i_gf)%l
      n = atoms%onsiteGF(i_gf)%atomType 

      DO ispin = 1, jspins
         n_tmp(:,:) = CMPLX(0.0,0.0)
         DO m = -l, l
            DO mp = -l, l
               DO i = 1, g%nef
                  IF(g%nr(i_gf).NE.1) THEN
                     CALL intgr3(REAL(g%gmmpMat(:,i,i_gf,m,mp,ispin,1)-g%gmmpMat(:,i,i_gf,m,mp,ispin,2)),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),re)
                     CALL intgr3(AIMAG(g%gmmpMat(:,i,i_gf,m,mp,ispin,1)-g%gmmpMat(:,i,i_gf,m,mp,ispin,2)),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imag)

                     n_tmp(m,mp) = n_tmp(m,mp) + (re+ImagUnit*imag)*g%de(i)

                  ELSE  
                     n_tmp(m,mp) = n_tmp(m,mp) + 1/2.0 * AIMAG(g%gmmpMat(1,i,i_gf,m,mp,ispin,1)*g%de(i)-g%gmmpMat(1,i,i_gf,m,mp,ispin,2)*conjg(g%de(i)))
                  END IF
               ENDDO

               mmpMat(m,mp,ispin) = -1/pi_const * n_tmp(m,mp)
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE occmtx

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
      REAL :: dos(g%nz,jspins)

      INTEGER ispin,m,j,l,n,i
      REAL    imagpl,imagmi

      dos = 0.0

      IF(g%mode.EQ.2) CALL juDFT_warn("The green's function is calculated on a semicircle in the complex plane; The ldos might not make sense", calledby="ldosmtx")
      !n(E) = 1/2pii * Tr(G^+-G^-)
      l = atoms%onsiteGF(i_gf)%l
      n = atoms%onsiteGF(i_gf)%atomType
      DO ispin = 1, jspins
         DO m = -l , l
            DO j = 1, g%nz
               IF(g%nr(i_gf).EQ.1) THEN
                  imagpl = AIMAG(g%gmmpMat(1,j,i_gf,m,m,ispin,1))
                  imagmi = AIMAG(g%gmmpMat(1,j,i_gf,m,m,ispin,2))
               ELSE
                  CALL intgr3(AIMAG(g%gmmpMat(:,j,i_gf,m,m,ispin,1)),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imagpl)
                  CALL intgr3(AIMAG(g%gmmpMat(:,j,i_gf,m,m,ispin,2)),atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n),imagmi)
               END IF
               dos(j,ispin) = dos(j,ispin) - 1/(2*pi_const) * (imagpl-imagmi)
            ENDDO
         ENDDO
      ENDDO

      OPEN(1337,file="lDOS_up_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

      DO i = 1, g%nz
         WRITE(1337,*) REAL(g%e(i)), dos(i,1)
      ENDDO

      CLOSE(unit = 1337)
      IF(jspins.EQ.2) THEN
         OPEN(1337,file="lDOS_dwn_" // TRIM(ADJUSTL(app)) // ".txt",action="write",status="replace")

         DO i = 1, g%nz
            WRITE(1337,*) REAL(g%e(i)), dos(i,2)
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



END MODULE m_gfcalc