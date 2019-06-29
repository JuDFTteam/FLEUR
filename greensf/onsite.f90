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

LOGICAL, PARAMETER :: l_debug = .TRUE.
INTEGER, PARAMETER :: int_method(3) = (/3,3,3/)

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
                     greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) -  pi_const * weight *&
                                                               REAL(conjg(eigVecCoeffs%acof(ib,lmp,natom,ispin))*eigVecCoeffs%acof(ib,lm,natom,ispin) +&
                                                                     conjg(eigVecCoeffs%bcof(ib,lmp,natom,ispin))*eigVecCoeffs%bcof(ib,lm,natom,ispin) *&
                                                                     usdus%ddn(l,nType,ispin))
                     IF(.NOT.input%l_gfsphavg) THEN
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
                        DO ilop = 1, atoms%nlo(nType)
                           IF (atoms%llo(ilop,nType).EQ.l) THEN
                              greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) = greensfCoeffs%projdos(ie,i_gf,m,mp,ispin) &
                                                                  - pi_const * weight * usdus%uloulopn(ilo,ilop,nType,ispin) *&
                                                                  conjg(eigVecCoeffs%ccof(mp,ib,ilop,natom,ispin)) *eigVecCoeffs%ccof(m,ib,ilo,natom,ispin)
   
                           ENDIF
                        ENDDO
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
   USE m_gfcalc
   USE m_kk_cutoff

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
   REAL    fac, n_occ

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
            IF(.NOT.input%l_gfsphavg) THEN
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
      CALL kk_cutoff(greensfCoeffs%projdos(:,i_gf,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,:),atoms,noco,&
                     l,input%jspins,greensfCoeffs%ne,greensfCoeffs%del,greensfCoeffs%e_bot,greensfCoeffs%e_top,&
                     greensfCoeffs%kkintgr_cutoff(i_gf,:,:))
      !
      ! Set the imaginary part to 0 outside the energy cutoffs
      !
     
      DO jspin = 1, input%jspins
         DO ie = 1, greensfCoeffs%ne
            IF( ie.GE.greensfCoeffs%kkintgr_cutoff(i_gf,jspin,1).AND.ie.LE.greensfCoeffs%kkintgr_cutoff(i_gf,jspin,2) ) CYCLE
            DO m= -l,l
               DO mp= -l,l
                  IF(input%l_gfsphavg) THEN
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
                  IF(input%l_gfsphavg) THEN
                     CALL kkintgr(greensfCoeffs%projdos(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%gmmpMat(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,int_method(gOnsite%mode))
                  ELSE
                  ! In the case of radial dependence we perform the kramers-kronig-integration seperately for uu,dd,etc.
                  ! We can do this because the radial functions are independent of E
                     CALL kkintgr(greensfCoeffs%uu(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%uu(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,int_method(gOnsite%mode))
                     CALL kkintgr(greensfCoeffs%dd(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%dd(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,int_method(gOnsite%mode))
                     CALL kkintgr(greensfCoeffs%du(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%du(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,int_method(gOnsite%mode))
                     CALL kkintgr(greensfCoeffs%ud(:,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,greensfCoeffs%ne,&
                                 gOnsite%ud(:,i_gf,m,mp,jspin,ipm),gOnsite%e,(ipm.EQ.2),gOnsite%mode,gOnsite%nz,int_method(gOnsite%mode))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL timestop("On-Site: Kramer-Kronigs-Integration")
      !TODO: Add checks for the Green's function
   ENDDO

END SUBROUTINE calc_onsite

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
      fac = 1.0/(sym%invarind(natom)*atoms%neq(nType))
      IF(sym%invarind(natom).EQ.0) CALL juDFT_error("No symmetry operations",calledby="local_sym")
      DO it = 1, sym%invarind(natom)
         is = sym%invarop(natom,it)
         isi = sym%invtab(is)
         d_mat(:,:) = cmplx(0.0,0.0)
         DO m = -l,l
            DO mp = -l,l
               d_mat(m,mp) = sym%d_wgn(m,mp,l,isi)
            ENDDO
         ENDDO
         calc_mat = matmul( transpose( conjg(d_mat) ) , orig_mat)
         calc_mat =  matmul( calc_mat, d_mat )
         DO m = -l,l
            DO mp = -l,l
               mat(m,mp) = mat(m,mp) + fac * conjg(calc_mat(m,mp))
            ENDDO
         ENDDO
      ENDDO
   ENDDO

END SUBROUTINE local_sym

END MODULE m_onsite