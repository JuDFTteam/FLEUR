MODULE m_onsite

!------------------------------------------------------------------------------
!
! MODULE: m_onsite
!
!> @author
!> Henning Janßen
!
! DESCRIPTION: 
!>  This module contains the functions to calculate the imaginary part of the 
!>  onsite GF with and without radial dependence
!>  Further we can transform this imaginary part to obtain the onsite GF
!>  using the Kramer Kronig Transformation
!
!------------------------------------------------------------------------------

USE m_juDFT
USE m_types
USE m_constants

LOGICAL, PARAMETER :: l_debug = .FALSE.

CONTAINS

SUBROUTINE onsite_coeffs(atoms,input,ispin,nbands,tetweights,ind,wtkpt,eig,usdus,eigVecCoeffs,greensfCoeffs)

   !This Subroutine calculates the contribution to the imaginary part of the Matrix elements G^[n \sigma]_{Lm Lm'}(E+i*sigma)
   !of the current k-Point (it is called in cdnval) inside the MT-sphere 
   !and sums over the Brillouin-Zone using the histogram method or linear tetrahedron method
   !It is essentially the l-density of states in a (m,mp) matrix with an additional factor - pi

   IMPLICIT NONE

   !-Type Arguments
   TYPE(t_atoms),         INTENT(IN)    :: atoms
   TYPE(t_eigVecCoeffs),  INTENT(IN)    :: eigVecCoeffs
   TYPE(t_usdus),         INTENT(IN)    :: usdus
   TYPE(t_greensfCoeffs), INTENT(INOUT) :: greensfCoeffs
   TYPE(t_input),         INTENT(IN)    :: input 

   !-Scalar Arguments 
   INTEGER,               INTENT(IN)    :: ispin  !Current spin index
   INTEGER,               INTENT(IN)    :: nbands !Number of bands to be considered
   REAL,                  INTENT(IN)    :: wtkpt  !Weight of the current k-point

   !-Array Arguments
   REAL,                  INTENT(IN)    :: tetweights(greensfCoeffs%ne,nbands) !Precalculated tetrahedron weights for the current k-point
   INTEGER,               INTENT(IN)    :: ind(nbands,2)                       !Gives the range where the tetrahedron weights are non-zero
   REAL,                  INTENT(IN)    :: eig(nbands)                         !Eigenvalues for the current k-point

   !-Local Scalars
   LOGICAL l_zero
   INTEGER i_gf,ib,ie,j,nType,natom,l,m,mp,lm,lmp,ilo,ilop
   REAL    weight


   !Loop through the gf elements to be calculated
   DO i_gf = 1, atoms%n_gf

      l     = atoms%onsiteGF(i_gf)%l
      nType = atoms%onsiteGF(i_gf)%atomType

      !Loop through equivalent atoms
      DO natom = SUM(atoms%neq(:nType-1)) + 1, SUM(atoms%neq(:nType))
         !Loop through bands
         DO ib = 1, nbands

            !Check wether there is a non-zero weight for the energy window
            l_zero = .true.
            IF(input%tria) THEN
               !TETRAHEDRON METHOD: check if the weight for this eigenvalue is non zero
               IF(ANY(tetweights(ind(ib,1):ind(ib,2),ib).NE.0.0)) l_zero = .false.
            ELSE
               !HISTOGRAM METHOD: check if eigenvalue is inside the energy range
               j = NINT((eig(ib)-greensfCoeffs%e_bot)/greensfCoeffs%del)+1
               IF( (j.LE.greensfCoeffs%ne).AND.(j.GE.1) )         l_zero = .false.
            END IF

            IF(l_zero) CYCLE 

            !$OMP PARALLEL DEFAULT(none) &
            !$OMP SHARED(j,ib,natom,l,nType,ispin,wtkpt,i_gf) &
            !$OMP SHARED(atoms,input,eigVecCoeffs,usdus,greensfCoeffs,eig,tetweights,ind) &
            !$OMP PRIVATE(ie,m,mp,lm,lmp,ilo,ilop,weight)

            !$OMP DO
            DO m = -l, l
               lm = l*(l+1)+m
               DO mp = -l,l
                  lmp = l*(l+1)+mp
                  !Choose the relevant energy points depending on the bz-integration method
                  DO ie = MERGE(ind(ib,1),j,input%tria), MERGE(ind(ib,2),j,input%tria)
                     !weight for the bz-integration including spin-degeneracy
                     weight = 2.0/input%jspins * 1.0/atoms%neq(nType) * MERGE(tetweights(ie,ib),wtkpt/greensfCoeffs%del,input%tria)  
                     !
                     !Contribution from states
                     !
                     !CHANGE: We drop the dependency of u/u_dot on spin and average over the two 
                     greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                               REAL(conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) +&
                                                                     conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin) *&
                                                                     SUM(usdus%ddn(l,nType,:)/input%jspins))
                     IF(.NOT.input%onsite_sphavg) THEN
                        greensfCoeffs%uu(ie,i_gf,m,mp,ispin) = greensfCoeffs%uu(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                               conjg(eigVecCoeffs%acof(ib,lm,natom,ispin))*eigVecCoeffs%acof(ib,lmp,natom,ispin)
                        greensfCoeffs%dd(ie,i_gf,m,mp,ispin) = greensfCoeffs%dd(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                               conjg(eigVecCoeffs%bcof(ib,lm,natom,ispin))*eigVecCoeffs%bcof(ib,lmp,natom,ispin)
                        greensfCoeffs%ud(ie,i_gf,m,mp,ispin) = greensfCoeffs%ud(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                               conjg(eigVecCoeffs%acof(ib,lm,natom,ispin))*eigVecCoeffs%bcof(ib,lmp,natom,ispin)
                        greensfCoeffs%du(ie,i_gf,m,mp,ispin) = greensfCoeffs%uu(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                               conjg(eigVecCoeffs%bcof(ib,lm,natom,ispin))*eigVecCoeffs%acof(ib,lmp,natom,ispin)
                     END IF
                     !
                     ! add local orbital contribution (not implemented for radial dependence yet and not tested for average)
                     !
                     DO ilo = 1, atoms%nlo(nType)
                        IF(atoms%llo(ilo,nType).EQ.l) THEN
                           greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) &
                                                               - pi_const * weight * (  usdus%uulon(ilo,nType,ispin) * (&
                                                               conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                                               conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) )&
                                                               + usdus%dulon(ilo,nType,ispin) * (&
                                                               conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%ccof(m,ib,ilo,natom,ispin) +&
                                                               conjg(eigVecCoeffs%ccof(mp,ib,ilo,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin)))
                        ENDIF
                     ENDDO
                     DO ilop = 1, atoms%nlo(nType)
                        IF (atoms%llo(ilop,nType).EQ.l) THEN
                           greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) &
                                                               - pi_const * weight * usdus%uloulopn(ilo,ilop,nType,ispin) *&
                                                               conjg(eigVecCoeffs%ccof(mp,ib,ilop,natom,ispin)) *eigVecCoeffs%ccof(m,ib,ilo,natom,ispin)

                        ENDIF
                     ENDDO
                  ENDDO! ie
               ENDDO !mp
            ENDDO !m
            !$OMP END DO
            !$OMP END PARALLEL
         ENDDO !ib
      ENDDO !natom
   ENDDO !i_gf

END SUBROUTINE onsite_coeffs

SUBROUTINE calc_onsite(atoms,input,noco,ef,greensfCoeffs,gOnsite,sym)

   USE m_kkintgr
   USE m_onsite21
   USE m_gfcalc

   IMPLICIT NONE

   !-Type Arguments
   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_greensfCoeffs),  INTENT(INOUT)  :: greensfCoeffs     !This is INTENT(INOUT) because the projected dos is useful for other things 
   TYPE(t_greensf),        INTENT(INOUT)  :: gOnsite
   TYPE(t_sym),            INTENT(IN)     :: sym
   TYPE(t_noco),           INTENT(IN)     :: noco
   REAL,                   INTENT(IN)     :: ef
   TYPE(t_input),          INTENT(IN)     :: input

   !-Local Scalars
   INTEGER i_gf,ie,l,m,mp,nType,jspin,ipm
   REAL    fac

   COMPLEX mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,atoms%n_gf,input%jspins)

   DO i_gf = 1, atoms%n_gf
      l =     atoms%onsiteGF(i_gf)%l
      nType = atoms%onsiteGF(i_gf)%atomType
      !
      !Enforcing that the projected density of states follows the local symmetries
      !
      DO ie = 1, greensfCoeffs%ne
         DO jspin = 1, input%jspins
            CALL local_sym(greensfCoeffs%projdos(ie,i_gf,-l:l,-l:l,jspin),l,nType,sym,atoms)
            IF(.NOT.input%onsite_sphavg) THEN
               CALL local_sym(greensfCoeffs%uu(ie,i_gf,-l:l,-l:l,jspin),l,nType,sym,atoms)
               CALL local_sym(greensfCoeffs%dd(ie,i_gf,-l:l,-l:l,jspin),l,nType,sym,atoms)
               CALL local_sym(greensfCoeffs%du(ie,i_gf,-l:l,-l:l,jspin),l,nType,sym,atoms)
               CALL local_sym(greensfCoeffs%ud(ie,i_gf,-l:l,-l:l,jspin),l,nType,sym,atoms)
            ENDIF
         ENDDO
      ENDDO     
      !
      !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration 
      !
      CALL greensf_cutoff(greensfCoeffs%projdos(:,i_gf,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,:),atoms,l,input%jspins,greensfCoeffs%ne,greensfCoeffs%del&
                          ,greensfCoeffs%e_bot,greensfCoeffs%e_top,greensfCoeffs%kkintgr_cutoff(i_gf,:))
      !
      ! Set the imaginary part to 0 outside the energy cutoffs
      !
      DO ie = 1, greensfCoeffs%ne
         IF( ie.GE.greensfCoeffs%kkintgr_cutoff(i_gf,1).AND.ie.LE.greensfCoeffs%kkintgr_cutoff(i_gf,2) ) CYCLE
         DO jspin = 1, input%jspins
            DO m= -l,l
               DO mp= -l,l
                  IF(input%onsite_sphavg) THEN
                     greensfCoeffs%projdos(ie,i_gf,m,mp,jspin) = 0.0
                  ELSE
                     greensfCoeffs%uu(ie,i_gf,m,mp,jspin) = 0.0
                     greensfCoeffs%dd(ie,i_gf,m,mp,jspin) = 0.0
                     greensfCoeffs%du(ie,i_gf,m,mp,jspin) = 0.0
                     greensfCoeffs%ud(ie,i_gf,m,mp,jspin) = 0.0
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      !
      !Perform the Kramers-Kronig-Integration
      !
      CALL timestart("On-Site: Kramer-Kronigs-Integration")
      DO jspin = 1, input%jspins
         DO m= -l,l
            DO mp= -l,l
               DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))
                  IF(input%onsite_sphavg) THEN
                     CALL kkintgr(greensfCoeffs%projdos(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%gmmpMat(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,1)
                  ELSE
                  ! In the case of radial dependence we perform the kramers-kronig-integration seperately for uu,dd,etc.
                  ! We can do this because the radial functions are independent of E
                     CALL kkintgr(greensfCoeffs%uu(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%uu(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,1)
                     CALL kkintgr(greensfCoeffs%dd(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%dd(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,1)
                     CALL kkintgr(greensfCoeffs%du(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%du(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,1)
                     CALL kkintgr(greensfCoeffs%ud(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%ud(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,1)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL timestop("On-Site: Kramer-Kronigs-Integration")
      !TODO: Add checks for the Green's function
      IF(l_debug) THEN
         CALL occmtx(gOnsite,i_gf,atoms,sym,input,ef,mmpMat(:,:,i_gf,:))
         DO jspin = 1, input%jspins
            DO m = -l, l
               WRITE(*,*) jspin, m, REAL(mmpMat(m,m,i_gf,jspin))
            ENDDO
         ENDDO
         CALL ldosmtx("g",gOnsite,1,atoms,sym,input)
      ENDIF
   ENDDO

   !In the noco case we need to rotate into the global frame
   IF(input%onsite_sphavg.AND.noco%l_mperp) CALL rot_onsite(atoms,noco,gOnsite)

END SUBROUTINE calc_onsite


SUBROUTINE greensf_cutoff(im,atoms,l,jspins,ne,del,e_bot,e_top,cutoff)
   !This Subroutine determines the cutoff energy for the kramers-kronig-integration
   !This cutoff energy is defined so that the integral over the fDOS up to this cutoff 
   !is equal to 2*(2l+1) (the number of states in the correlated shell) or not to small
   
   USE m_kkintgr
   
   IMPLICIT NONE

   REAL,                INTENT(INOUT)  :: im(ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
   TYPE(t_atoms),       INTENT(IN)     :: atoms

   INTEGER,             INTENT(IN)     :: l
   INTEGER,             INTENT(IN)     :: jspins
   INTEGER,             INTENT(IN)     :: ne
   REAL,                INTENT(IN)     :: del
   REAL,                INTENT(IN)     :: e_bot
   REAL,                INTENT(IN)     :: e_top

   INTEGER,             INTENT(OUT)    :: cutoff(2)


   INTEGER i,m,n_c,ispin
   REAL integral
   REAL a,b, n_states
   LOGICAL l_write

   REAL :: fDOS(ne,jspins)

   l_write=.false.  !Debugging Output
   fDOS = 0.0

   !Calculate the trace over m,mp of the Greens-function matrix to obtain the fDOS 
   !n_f(e) = -1/pi * TR[Im(G_f(e))]
   DO ispin = 1, jspins
      DO m = -l , l
         DO i = 1, ne
            fDOS(i,ispin) = fDOS(i,ispin) + im(i,m,m,ispin)
         ENDDO
      ENDDO
   ENDDO
   fDOS(:,:) = -1/pi_const*fDOS(:,:)

   !For Debugging:
   IF(l_write) THEN
      
      OPEN(1337,file="fDOS_up.txt",action="write",status="replace")

      DO i = 1, ne
         WRITE(1337,*) ((i-1)*del + e_bot)*hartree_to_ev_const, fDOS(i,1)/hartree_to_ev_const
      ENDDO

      CLOSE(unit = 1337)

      IF(jspins.EQ.2) THEN
         OPEN(1337,file="fDOS_dwn.txt",action="write",status="replace")

         DO i = 1, ne
            WRITE(1337,*) ((i-1)*del + e_bot)*hartree_to_ev_const, -fDOS(i,2)/hartree_to_ev_const
         ENDDO

         CLOSE(unit = 1337)
      ENDIF
   ENDIF
   
   IF(jspins.EQ.2) fDOS(:,1) = fDOS(:,1) + fDOS(:,2)

   integral =  trapz(fDOS(1:ne,1), del, ne)

   n_states = 2*(2*l+1)
   
   IF(l_write) WRITE(*,*) "Integral over DOS: ", integral

   cutoff(1) = 1   !at the moment we don't modify the lower bound
   cutoff(2) = ne

   IF(integral.LT.n_states-0.1) THEN
      ! If the integral is to small we stop here to avoid problems
      CALL juDFT_error("integral over DOS too small", calledby="greensf_cutoff")
      
   ELSE IF((integral.GT.n_states).AND.((integral-n_states).GT.0.001)) THEN
      !IF the integral is bigger than 14, search for the cutoff using the bisection method   

      a = e_bot
      b = e_top

      DO

         cutoff(2) = INT(((a+b)/2.0-e_bot)/del)+1
         integral =  trapz(fDOS(1:cutoff(2),1),del,cutoff(2))

         IF((ABS(integral-n_states).LT.0.001).OR.(ABS(a-b)/2.0.LT.del)) THEN
            !The integral is inside the desired accuracy
            EXIT
         ELSE IF((integral-n_states).LT.0) THEN
            !integral to small -> choose the right interval
            a = (a+b)/2.0
         ELSE IF((integral-n_states).GT.0) THEN
            !integral to big   -> choose the left interval
            b = (a+b)/2.0
         END IF

      ENDDO

      IF(l_write) THEN
         WRITE(*,*) "CALCULATED CUTOFF: ", cutoff(2)
         WRITE(*,*) "INTEGRAL OVER fDOS with cutoff: ", integral
      ENDIF
   ENDIF

END SUBROUTINE greensf_cutoff

SUBROUTINE local_sym(mat,l,nType,sym,atoms)

   IMPLICIT NONE

   TYPE(t_sym),   INTENT(IN)    :: sym 
   TYPE(t_atoms), INTENT(IN)    :: atoms
   INTEGER,       INTENT(IN)    :: l 
   INTEGER,       INTENT(IN)    :: nType
   REAL,          INTENT(INOUT) :: mat(-l:l,-l:l)

   !-Local Scalars
   INTEGER natom,it,is,isi,m,mp,nop
   REAL    fac

   !-Local Arrays
   COMPLEX orig_mat(-l:l,-l:l), calc_mat(-l:l,-l:l), d_mat(-l:l,-l:l), diag(-l:l)

   orig_mat(-l:l,-l:l) = mat(-l:l,-l:l)

   mat = 0.0
   nop = 0
   DO natom = SUM(atoms%neq(:nType-1)) + 1, SUM(atoms%neq(:nType))
      DO it = 1, sym%invarind(natom)
         is = sym%invarop(natom,it)
         isi = sym%invtab(is)
         d_mat(:,:) = cmplx(0.0,0.0)
         DO m = -l,l
            DO mp = -l,l
               d_mat(m,mp) = sym%d_wgn(m,mp,l,isi)
            ENDDO
         ENDDO
         DO m = -l,l
            diag(m) = d_mat(m,m)
            d_mat(m,m) = 0.0
         ENDDO
         !Exclude all symmetries that would prevent splitting of the levels
         IF(ANY(d_mat(:,:).NE.0.0)) CYCLE
         nop = nop + 1
         DO m = -l,l
            d_mat(m,m) = diag(m)
         ENDDO
         calc_mat = matmul( transpose( conjg(d_mat) ) , orig_mat)
         calc_mat =  matmul( calc_mat, d_mat )
         DO m = -l,l
            DO mp = -l,l
               mat(m,mp) = mat(m,mp) + conjg(calc_mat(m,mp))
            ENDDO
         ENDDO
      ENDDO
   ENDDO
   IF(nop.EQ.0) CALL juDFT_error("No symmetry operations found",calledby="local_sym")
   fac = 1.0/(REAL(nop))
   mat = mat * fac

END SUBROUTINE local_sym

END MODULE m_onsite