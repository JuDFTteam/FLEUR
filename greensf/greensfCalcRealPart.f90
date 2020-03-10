MODULE m_greensfCalcRealPart

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
   USE m_kkintgr
   USE m_kk_cutoff

   IMPLICIT NONE

   INTEGER, PARAMETER :: int_method(3) = (/3,3,1/)

   CONTAINS

   SUBROUTINE greensfCalcRealPart(atoms,gfinp,input,sym,noco,ef,greensfImagPart,g)

      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_gfinp),          INTENT(IN)     :: gfinp
      TYPE(t_greensfImagPart),INTENT(INOUT)  :: greensfImagPart     !This is INTENT(INOUT) because the projected dos is useful for other things 
      TYPE(t_greensf),        INTENT(INOUT)  :: g
      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_input),          INTENT(IN)     :: input
      REAL,                   INTENT(IN)     :: ef

      INTEGER i_gf,ie,l,m,mp,nType,jspin,ipm,kkcut,lp,nTypep,spin_cut,nn,natom
      REAL    fac,del,eb,et

      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(ef,del,eb,et)

      DO i_gf = 1, gfinp%n
         l =      gfinp%elem(i_gf)%l
         lp =     gfinp%elem(i_gf)%lp
         nType =  gfinp%elem(i_gf)%atomType
         nTypep = gfinp%elem(i_gf)%atomTypep
         CALL timestart("On-Site: Integration Cutoff")
         IF(nType.EQ.nTypep.AND.l.EQ.lp) THEN
            !
            !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration
            !
            CALL kk_cutoff(greensfImagPart%sphavg(:,:,:,i_gf,:),noco,l,input%jspins,&
                           gfinp%ne,del,eb,et,greensfImagPart%kkintgr_cutoff(i_gf,:,:))
         ELSE
            !For all other elements we just use ef+elup as a hard cutoff
            greensfImagPart%kkintgr_cutoff(i_gf,:,1) = 1
            greensfImagPart%kkintgr_cutoff(i_gf,:,2) = gfinp%ne
         ENDIF
         CALL timestop("On-Site: Integration Cutoff")
         !
         !Perform the Kramers-Kronig-Integration if not already calculated
         !
         CALL timestart("On-Site: Kramer-Kronigs-Integration")
         DO jspin = 1, MERGE(3,input%jspins,gfinp%l_mperp)
            spin_cut = MERGE(1,jspin,jspin.GT.2)
            kkcut = greensfImagPart%kkintgr_cutoff(i_gf,spin_cut,2)
            DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))
               DO m= -l,l
                  DO mp= -lp,lp
                     CALL kkintgr(greensfImagPart%sphavg(:,m,mp,i_gf,jspin),eb,del,kkcut,&
                                 g%gmmpMat(:,m,mp,jspin,ipm,i_gf),g%e,(ipm.EQ.2),gfinp%mode,g%nz,int_method(gfinp%mode))
                     IF(.NOT.gfinp%l_sphavg) THEN
                        ! In the case of radial dependence we perform the kramers-kronig-integration seperately for uu,dd,etc.
                        ! We can do this because the radial functions are independent of E
                        CALL kkintgr(greensfImagPart%uu(:,m,mp,i_gf,jspin),eb,del,kkcut,&
                                    g%uu(:,m,mp,jspin,ipm,i_gf),g%e,(ipm.EQ.2),gfinp%mode,g%nz,int_method(gfinp%mode))
                        CALL kkintgr(greensfImagPart%dd(:,m,mp,i_gf,jspin),eb,del,kkcut,&
                                    g%dd(:,m,mp,jspin,ipm,i_gf),g%e,(ipm.EQ.2),gfinp%mode,g%nz,int_method(gfinp%mode))
                        CALL kkintgr(greensfImagPart%du(:,m,mp,i_gf,jspin),eb,del,kkcut,&
                                    g%du(:,m,mp,jspin,ipm,i_gf),g%e,(ipm.EQ.2),gfinp%mode,g%nz,int_method(gfinp%mode))
                        CALL kkintgr(greensfImagPart%ud(:,m,mp,i_gf,jspin),eb,del,kkcut,&
                                    g%ud(:,m,mp,jspin,ipm,i_gf),g%e,(ipm.EQ.2),gfinp%mode,g%nz,int_method(gfinp%mode))
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL timestop("On-Site: Kramer-Kronigs-Integration")
      ENDDO

   END SUBROUTINE greensfCalcRealPart
END MODULE m_greensfCalcRealPart