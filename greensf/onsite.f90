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
INTEGER, PARAMETER :: int_method(3) = (/3,3,1/)

CONTAINS

SUBROUTINE calc_onsite(atoms,input,sym,noco,angle,greensfCoeffs,g)

   USE m_kkintgr
   USE m_kk_cutoff

   IMPLICIT NONE

   !-Type Arguments
   TYPE(t_atoms),          INTENT(IN)     :: atoms
   TYPE(t_greensfCoeffs),  INTENT(INOUT)  :: greensfCoeffs     !This is INTENT(INOUT) because the projected dos is useful for other things 
   TYPE(t_greensf),        INTENT(INOUT)  :: g
   TYPE(t_sym),            INTENT(IN)     :: sym
   TYPE(t_noco),           INTENT(IN)     :: noco
   TYPE(t_input),          INTENT(IN)     :: input
   REAL,                   INTENT(IN)     :: angle(sym%nop)

   !-Local Scalars
   INTEGER i_gf,ie,l,m,mp,nType,jspin,ipm,kkcut,lp,nTypep,spin_cut,nn,natom
   REAL    fac
   INTEGER it,is,isi
   COMPLEX phase
   COMPLEX d_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const),calc_mat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
   COMPLEX g21(g%nz,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

   DO i_gf = 1, atoms%n_gf
      l =     atoms%gfelem(i_gf)%l
      lp =     atoms%gfelem(i_gf)%lp
      nType = atoms%gfelem(i_gf)%atomType
      nTypep = atoms%gfelem(i_gf)%atomTypep
      !
      !Enforcing that the projected density of states follows the local symmetries
      !
      CALL timestart("On-Site: Integration Cutoff")
      IF(nType.EQ.nTypep.AND.l.EQ.lp) THEN
         !
         !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration
         !
         CALL kk_cutoff(greensfCoeffs%projdos(:,i_gf,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,1:input%jspins),atoms,noco,&
                        l,input%jspins,greensfCoeffs%ne,greensfCoeffs%del,greensfCoeffs%e_bot,greensfCoeffs%e_top,&
                        greensfCoeffs%kkintgr_cutoff(i_gf,:,:))
      ELSE
         !For all other elements we just use ef+elup as a hard cutoff
         greensfCoeffs%kkintgr_cutoff(i_gf,:,1) = 1
         greensfCoeffs%kkintgr_cutoff(i_gf,:,2) = greensfCoeffs%ne
      ENDIF
      CALL timestop("On-Site: Integration Cutoff")
      !
      !Perform the Kramers-Kronig-Integration if not already calculated
      !
      CALL timestart("On-Site: Kramer-Kronigs-Integration")
      DO jspin = 1, input%jspins
         spin_cut = MERGE(1,jspin,jspin.GT.2)
         kkcut = greensfCoeffs%kkintgr_cutoff(i_gf,spin_cut,2)
         DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))
            DO m= -l,l
               DO mp= -lp,lp
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
      IF(input%l_gfmperp) THEN
         kkcut = greensfCoeffs%kkintgr_cutoff(i_gf,1,2)
         DO nn = 1, atoms%neq(nType)
            natom = SUM(atoms%neq(:nType-1)) + nn
            DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))
               DO m= -l,l
                  DO mp= -lp,lp
                     CALL kkintgr(greensfCoeffs%projdos21(1:kkcut,i_gf,nn,m,mp),greensfCoeffs%e_bot,greensfCoeffs%del,kkcut,&
                                 g21(:,m,mp),g%e,(ipm.EQ.2),g%mode,g%nz,int_method(g%mode))
                  ENDDO
               ENDDO
               fac = 1.0/(sym%invarind(natom)*atoms%neq(nType))
               IF(sym%invarind(natom).EQ.0) CALL juDFT_error("No symmetry operations available",calledby="greensfImag")
               DO it = 1, sym%invarind(natom)
                  is = sym%invarop(natom,it)
                  isi = sym%invtab(is)
                  d_mat(:,:) = cmplx(0.0,0.0)
                  DO m = -l,l
                     DO mp = -l,l
                        d_mat(m,mp) = sym%d_wgn(m,mp,l,isi)
                     ENDDO
                  ENDDO
                  DO ie = 1, g%nz
                     calc_mat = matmul( transpose( conjg(d_mat) ) , &
                                 g21(ie,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
                     calc_mat =  matmul( calc_mat, d_mat )
                     phase = exp(ImagUnit*angle(isi))
                     DO m = -l,l
                        DO mp = -l,l
                           g%gmmpMat(ie,i_gf,m,mp,3,ipm) =&
                           g%gmmpMat(ie,i_gf,m,mp,3,ipm) + phase * fac * conjg(calc_mat(m,mp))
                        ENDDO
                     ENDDO
                  ENDDO!ie
               ENDDO!it
            ENDDO
         ENDDO!natom
      ENDIF

      CALL timestop("On-Site: Kramer-Kronigs-Integration")
      !IF(input%l_gfmperp) THEN
      !   CALL rot_gf_mat(g,noco)
      !ENDIF
   ENDDO

END SUBROUTINE calc_onsite

END MODULE m_onsite