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

SUBROUTINE calc_onsite(atoms,input,noco,greensfCoeffs,g,sym)

   USE m_kkintgr
   USE m_gfcalc
   USE m_kk_cutoff

   IMPLICIT NONE

   !-Type Arguments
   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_greensfCoeffs),  INTENT(INOUT)  :: greensfCoeffs     !This is INTENT(INOUT) because the projected dos is useful for other things 
   TYPE(t_greensf),        INTENT(INOUT)  :: g
   TYPE(t_sym),            INTENT(IN)     :: sym
   TYPE(t_noco),           INTENT(IN)     :: noco
   TYPE(t_input),          INTENT(IN)     :: input

   !-Local Scalars
   INTEGER i_gf,ie,l,m,mp,nType,jspin,ipm,kkcut,lp,nTypep,spin_cut
   REAL    fac

   DO i_gf = 1, atoms%n_gf
      l =     atoms%gfelem(i_gf)%l
      lp =     atoms%gfelem(i_gf)%lp
      nType = atoms%gfelem(i_gf)%atomType
      nTypep = atoms%gfelem(i_gf)%atomTypep
      !
      !Enforcing that the projected density of states follows the local symmetries
      !
      IF(nType.EQ.nTypep.AND.l.EQ.lp) THEN   
         !
         !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration 
         !
         CALL kk_cutoff(greensfCoeffs%projdos(:,i_gf,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,:),atoms,noco,&
                        l,input%jspins,greensfCoeffs%ne,greensfCoeffs%del,greensfCoeffs%e_bot,greensfCoeffs%e_top,&
                        greensfCoeffs%kkintgr_cutoff(i_gf,:,:))
      ELSE 
         !For all other elements we just use ef+elup as a hard cutoff
         greensfCoeffs%kkintgr_cutoff(i_gf,:,1) = 1
         greensfCoeffs%kkintgr_cutoff(i_gf,:,2) = greensfCoeffs%ne
      ENDIF
      !
      !Perform the Kramers-Kronig-Integration
      !
      CALL timestart("On-Site: Kramer-Kronigs-Integration")
      DO jspin = 1, MERGE(3,input%jspins,input%l_gfmperp)
         spin_cut = MERGE(1,jspin,jspin.GT.2)
         kkcut = greensfCoeffs%kkintgr_cutoff(i_gf,spin_cut,2)
         DO m= -l,l
            DO mp= -lp,lp
               DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))
                  CALL kkintgr(greensfCoeffs%projdos(1:kkcut,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,kkcut,&
                              g%gmmpMat(:,i_gf,m,mp,jspin,ipm),g%e,(ipm.EQ.2),g%mode,g%nz,int_method(g%mode))
                  IF(.NOT.input%l_gfsphavg) THEN
                     ! In the case of radial dependence we perform the kramers-kronig-integration seperately for uu,dd,etc.
                     ! We can do this because the radial functions are independent of E
                     CALL kkintgr(greensfCoeffs%uu(1:kkcut,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,kkcut,&
                                 g%uu(:,i_gf,m,mp,jspin,ipm),g%e,(ipm.EQ.2),g%mode,g%nz,int_method(g%mode))
                     CALL kkintgr(greensfCoeffs%dd(1:kkcut,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,kkcut,&
                                 g%dd(:,i_gf,m,mp,jspin,ipm),g%e,(ipm.EQ.2),g%mode,g%nz,int_method(g%mode))
                     CALL kkintgr(greensfCoeffs%du(1:kkcut,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,kkcut,&
                                 g%du(:,i_gf,m,mp,jspin,ipm),g%e,(ipm.EQ.2),g%mode,g%nz,int_method(g%mode))
                     CALL kkintgr(greensfCoeffs%ud(1:kkcut,i_gf,m,mp,jspin),greensfCoeffs%e_bot,greensfCoeffs%del,kkcut,&
                                 g%ud(:,i_gf,m,mp,jspin,ipm),g%e,(ipm.EQ.2),g%mode,g%nz,int_method(g%mode))
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      CALL timestop("On-Site: Kramer-Kronigs-Integration")
      !IF(input%l_gfmperp) THEN
      !   CALL rot_gf_mat(g,noco)
      !ENDIF
   ENDDO

END SUBROUTINE calc_onsite

!SUBROUTINE rot_gf_mat(g,noco,i_gf,nType)
!
!   USE m_rotdenmat
!
!   IMPLICIT NONE 
!
!   TYPE(t_greensf),     INTENT(INOUT)  :: g 
!   TYPE(t_noco),        INTENT(IN)     :: noco
!   INTEGER,             INTENT(IN)     :: i_gf  
!   INTEGER,             INTENT(IN)     :: nType
!
!   CALL rot_den_mat(noco%alph(nType),noco%beta(nType),&
!        g%gmmpMat(iz,i_gf,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,1,ipm),&
!        g%gmmpMat(iz,i_gf,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,2,ipm),&
!        g%gmmpMat(iz,i_gf,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,3,ipm))
!
!END SUBROUTINE

END MODULE m_onsite