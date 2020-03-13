MODULE m_greensfCalcRealPart

   !------------------------------------------------------------------------------
   !
   ! MODULE: m_greensfCalcRealPart
   !
   !> @author
   !> Henning Janßen
   !
   ! DESCRIPTION:
   !>  This module contains the functions to calculate the imaginary part of the
   !>  onsite GF with and without radial dependence
   !>  Further we can transform this imaginary part to obtain the Green's Function
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
      TYPE(t_greensf),        INTENT(INOUT)  :: g(:)
      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_input),          INTENT(IN)     :: input
      REAL,                   INTENT(IN)     :: ef

      INTEGER i_gf,i_elem,ie,l,m,mp,nType,jspin,ipm,kkcut,lp,nTypep,spin_cut,nn,natom,contourShape,dummy
      REAL    fac,del,eb,et

      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(ef,del,eb,et)

      DO i_gf = 1, gfinp%n
         l =      gfinp%elem(i_gf)%l
         lp =     gfinp%elem(i_gf)%lp
         nType =  gfinp%elem(i_gf)%atomType
         nTypep = gfinp%elem(i_gf)%atomTypep
         contourShape = gfinp%contour(gfinp%elem(i_gf)%iContour)%shape

         dummy = gfinp%uniqueElements(ind=i_gf,indUnique=i_elem)

         CALL timestart("On-Site: Integration Cutoff")
         IF(nType.EQ.nTypep.AND.l.EQ.lp.AND.gfinp%l_sphavg) THEN
            !
            !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration
            !
            CALL kk_cutoff(greensfImagPart%sphavg(:,:,:,i_elem,:),noco,l,input%jspins,&
                           gfinp%ne,del,eb,et,greensfImagPart%kkintgr_cutoff(i_gf,:,:))
         ELSE
            !For all other elements we just use ef+elup as a hard cutoff
            !(maybe give option to specify outside of changing the realAxis grid)
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
                     IF(gfinp%l_sphavg) THEN
                        CALL kkintgr(greensfImagPart%sphavg(:,m,mp,i_elem,jspin),eb,del,kkcut,&
                                     g(i_gf)%gmmpMat(:,m,mp,jspin,ipm),g(i_gf)%contour%e,(ipm.EQ.2),contourShape,g(i_gf)%contour%nz,int_method(contourShape))
                     ELSE
                        ! In the case of radial dependence we perform the kramers-kronig-integration seperately for uu,dd,etc.
                        ! We can do this because the radial functions are independent of E
                        CALL kkintgr(greensfImagPart%uu(:,m,mp,i_elem,jspin),eb,del,kkcut,&
                                     g(i_gf)%uu(:,m,mp,jspin,ipm),g(i_gf)%contour%e,(ipm.EQ.2),contourShape,g(i_gf)%contour%nz,int_method(contourShape))
                        CALL kkintgr(greensfImagPart%dd(:,m,mp,i_elem,jspin),eb,del,kkcut,&
                                     g(i_gf)%dd(:,m,mp,jspin,ipm),g(i_gf)%contour%e,(ipm.EQ.2),contourShape,g(i_gf)%contour%nz,int_method(contourShape))
                        CALL kkintgr(greensfImagPart%du(:,m,mp,i_elem,jspin),eb,del,kkcut,&
                                     g(i_gf)%du(:,m,mp,jspin,ipm),g(i_gf)%contour%e,(ipm.EQ.2),contourShape,g(i_gf)%contour%nz,int_method(contourShape))
                        CALL kkintgr(greensfImagPart%ud(:,m,mp,i_elem,jspin),eb,del,kkcut,&
                                     g(i_gf)%ud(:,m,mp,jspin,ipm),g(i_gf)%contour%e,(ipm.EQ.2),contourShape,g(i_gf)%contour%nz,int_method(contourShape))
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL timestop("On-Site: Kramer-Kronigs-Integration")
      ENDDO

   END SUBROUTINE greensfCalcRealPart
END MODULE m_greensfCalcRealPart