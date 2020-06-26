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

   INTEGER, PARAMETER :: int_method(3) = [method_direct,method_direct,method_maclaurin]

   CONTAINS

   SUBROUTINE greensfCalcRealPart(atoms,gfinp,input,sym,noco,vTot,enpara,fmpi,hub1inp,ef,greensfImagPart,g)

      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_gfinp),          INTENT(IN)     :: gfinp
      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_input),          INTENT(IN)     :: input
      TYPE(t_potden),         INTENT(IN)     :: vTot
      TYPE(t_enpara),         INTENT(IN)     :: enpara
      TYPE(t_hub1inp),        INTENT(IN)     :: hub1inp
      TYPE(t_mpi),            INTENT(IN)     :: fmpi
      REAL,                   INTENT(IN)     :: ef
      TYPE(t_greensfImagPart),INTENT(INOUT)  :: greensfImagPart
      TYPE(t_greensf),        INTENT(INOUT)  :: g(:)

      INTEGER :: i_gf,i_elem,ie,l,m,mp,nType,indUnique
      INTEGER :: jspin,nspins,ipm,kkcut,lp,nTypep
      INTEGER :: spin_cut,contourShape
      INTEGER :: i_gf_start,i_gf_end,spin_start,spin_end
      INTEGER :: n_gf_task,extra
      LOGICAL :: l_onsite,l_fixedCutoffset,l_skip
      REAL    :: fac,del,eb,et,fixedCutoff
      REAL    :: scalingFactor(input%jspins)
      REAL,    ALLOCATABLE :: eMesh(:)

      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(ef,del,eb,eMesh=eMesh)

      nspins = MERGE(3,input%jspins,gfinp%l_mperp)

      IF(fmpi%irank.EQ.0) THEN
         CALL timestart("Green's Function: Integration Cutoff")
         DO i_gf = 1, gfinp%n

            !Get the information of ith current element
            l  = g(i_gf)%elem%l
            lp = g(i_gf)%elem%lp
            nType  = g(i_gf)%elem%atomType
            nTypep = g(i_gf)%elem%atomTypep
            l_fixedCutoffset = g(i_gf)%elem%l_fixedCutoffset
            fixedCutoff      = g(i_gf)%elem%fixedCutoff

            i_elem = gfinp%uniqueElements(ind=i_gf,indUnique=indUnique)

            IF(i_gf /= indUnique) THEN
               !This cutoff was already calculated
               greensfImagPart%kkintgr_cutoff(i_gf,:,:) = greensfImagPart%kkintgr_cutoff(indUnique,:,:)
            ELSE
               !Is the current element suitable for automatic finding of the cutoff
               l_onsite = nType.EQ.nTypep.AND.l.EQ.lp
               IF(l_onsite.AND..NOT.l_fixedCutoffset) THEN
                  !
                  !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration
                  !
                  IF(gfinp%l_sphavg) THEN
                     CALL kk_cutoff(greensfImagPart%sphavg(:,:,:,i_elem,:),noco,gfinp%l_mperp,l,input%jspins,&
                                    eMesh,greensfImagPart%kkintgr_cutoff(i_gf,:,:),scalingFactor)
                  ELSE
                     !Onsite element with radial dependence
                     CALL kk_cutoffRadial(greensfImagPart%uu(:,:,:,i_elem,:),greensfImagPart%ud(:,:,:,i_elem,:),&
                                          greensfImagPart%du(:,:,:,i_elem,:),greensfImagPart%dd(:,:,:,i_elem,:),&
                                          noco,atoms,vTot,enpara,fmpi,hub1inp,gfinp%l_mperp,l,nType,input,eMesh,&
                                          greensfImagPart%kkintgr_cutoff(i_gf,:,:),scalingFactor)
                  ENDIF
               ELSE IF (l_fixedCutoffset) THEN
                  greensfImagPart%kkintgr_cutoff(i_gf,:,1) = 1
                  greensfImagPart%kkintgr_cutoff(i_gf,:,2) = INT((fixedCutoff+ef-eb)/del)+1
               ELSE
                  !For all other elements we just use ef+elup as a hard cutoff
                  greensfImagPart%kkintgr_cutoff(i_gf,:,1) = 1
                  greensfImagPart%kkintgr_cutoff(i_gf,:,2) = gfinp%ne
               ENDIF
            ENDIF

            DO jspin = 1, nspins
               spin_cut = MERGE(1,jspin,jspin.GT.2)
               kkcut = greensfImagPart%kkintgr_cutoff(i_gf,spin_cut,2)
               !------------------------------------------------------------
               ! Set everything above the cutoff in the imaginary part to 0
               ! We do this explicitely because when we just use the hard cutoff index
               ! Things might get lost when the imaginary part is smoothed explicitely
               !------------------------------------------------------------
               IF(kkcut.ne.SIZE(eMesh)) THEN
                  IF(gfinp%l_sphavg) THEN
                     greensfImagPart%sphavg(kkcut+1:,-l:l,-l:l,i_elem,jspin) = 0.0
                  ELSE
                     greensfImagPart%uu(kkcut+1:,-l:l,-l:l,i_elem,jspin) = 0.0
                     greensfImagPart%ud(kkcut+1:,-l:l,-l:l,i_elem,jspin) = 0.0
                     greensfImagPart%du(kkcut+1:,-l:l,-l:l,i_elem,jspin) = 0.0
                     greensfImagPart%dd(kkcut+1:,-l:l,-l:l,i_elem,jspin) = 0.0
                  ENDIF
               ENDIF
               IF(nspins == 2) THEN
                  IF(gfinp%l_sphavg) THEN
                     greensfImagPart%sphavg(:kkcut,-l:l,-l:l,i_elem,jspin) = greensfImagPart%sphavg(:kkcut,-l:l,-l:l,i_elem,jspin) &
                                                                            * scalingFactor(jspin)
                  ELSE
                     greensfImagPart%uu(:kkcut,-l:l,-l:l,i_elem,jspin) = greensfImagPart%uu(:kkcut,-l:l,-l:l,i_elem,jspin) &
                                                                        * scalingFactor(jspin)
                     greensfImagPart%ud(:kkcut,-l:l,-l:l,i_elem,jspin) = greensfImagPart%ud(:kkcut,-l:l,-l:l,i_elem,jspin) &
                                                                        * scalingFactor(jspin)
                     greensfImagPart%du(:kkcut,-l:l,-l:l,i_elem,jspin) = greensfImagPart%du(:kkcut,-l:l,-l:l,i_elem,jspin) &
                                                                        * scalingFactor(jspin)
                     greensfImagPart%dd(:kkcut,-l:l,-l:l,i_elem,jspin) = greensfImagPart%dd(:kkcut,-l:l,-l:l,i_elem,jspin) &
                                                                        * scalingFactor(jspin)
                  ENDIF
               ENDIF
            ENDDO

         ENDDO
         CALL timestop("Green's Function: Integration Cutoff")
      ENDIF

      !Broadcast cutoffs and modified imaginary parts
      CALL greensfImagPart%mpi_bc(fmpi%mpi_comm)


      !Distribute the Calculations
#ifdef CPP_MPI
      IF(fmpi%isize>1) THEN
         IF(gfinp%n>=fmpi%isize) THEN
            !Just distribute the individual gf elements over the ranks
            n_gf_task = FLOOR(REAL(gfinp%n)/(fmpi%isize))
            extra = gfinp%n - n_gf_task*fmpi%isize
            i_gf_start = fmpi%irank*n_gf_task + 1 + extra
            i_gf_end = (fmpi%irank+1)*n_gf_task   + extra
            IF(fmpi%irank < extra) THEN
               i_gf_start = i_gf_start - (extra - fmpi%irank)
               i_gf_end = i_gf_end - (extra - fmpi%irank - 1)
            ENDIF
            spin_start = 1
            spin_end   = nspins
         ELSE IF(gfinp%n*nspins>fmpi%isize) THEN
            !Just fill up the ranks
            i_gf_start = fmpi%irank + 1
            i_gf_end   = fmpi%irank + 1
            spin_start = 1
            spin_end   = nspins
         ELSE
            !If there are few enough gf elements then distribute the spins
            spin_start = MOD(fmpi%irank,nspins) + 1
            spin_end   = MOD(fmpi%irank,nspins) + 1
            i_gf_start = 1 + FLOOR(REAL(fmpi%irank)/nspins)
            i_gf_end   = 1 + FLOOR(REAL(fmpi%irank)/nspins)
         ENDIF
      ELSE
         !Distribute nothing
         i_gf_start = 1
         i_gf_end = gfinp%n
         spin_start = 1
         spin_end   = nspins
      ENDIF
#else
      i_gf_start = 1
      i_gf_end = gfinp%n
      spin_start = 1
      spin_end   = nspins
#endif

      DO i_gf = i_gf_start, i_gf_end

         IF(i_gf.LT.1 .OR. i_gf.GT.gfinp%n) CYCLE !Make sure to not produce segfaults with mpi

         !Get the information of ith current element
         l  = g(i_gf)%elem%l
         lp = g(i_gf)%elem%lp
         nType  = g(i_gf)%elem%atomType
         nTypep = g(i_gf)%elem%atomTypep
         contourShape = gfinp%contour(g(i_gf)%elem%iContour)%shape

         i_elem = gfinp%uniqueElements(ind=i_gf)

         CALL timestart("Green's Function: Kramer-Kronigs-Integration")
         DO jspin = spin_start, spin_end
            DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))
               DO m= -l,l
                  DO mp= -lp,lp

                     !Don't waste time on empty elements
                     l_skip = .FALSE.
                     DO ie = 1, SIZE(eMesh)
                        IF(gfinp%l_sphavg) THEN
                           IF(ABS(greensfImagPart%sphavg(ie,m,mp,i_elem,jspin)).GT.1e-12) EXIT
                        ELSE
                           IF(ABS(greensfImagPart%uu(ie,m,mp,i_elem,jspin)).GT.1e-12) EXIT
                           IF(ABS(greensfImagPart%ud(ie,m,mp,i_elem,jspin)).GT.1e-12) EXIT
                           IF(ABS(greensfImagPart%du(ie,m,mp,i_elem,jspin)).GT.1e-12) EXIT
                           IF(ABS(greensfImagPart%dd(ie,m,mp,i_elem,jspin)).GT.1e-12) EXIT
                        ENDIF
                        IF(ie==SIZE(eMesh)) l_skip = .TRUE.
                     ENDDO
                     IF(l_skip) THEN
                        IF(gfinp%l_sphavg) THEN
                           g(i_gf)%gmmpMat(:,m,mp,jspin,ipm) = cmplx_0
                        ELSE
                           g(i_gf)%uu(:,m,mp,jspin,ipm) = cmplx_0
                           g(i_gf)%ud(:,m,mp,jspin,ipm) = cmplx_0
                           g(i_gf)%du(:,m,mp,jspin,ipm) = cmplx_0
                           g(i_gf)%dd(:,m,mp,jspin,ipm) = cmplx_0
                        ENDIF
                        CYCLE
                     ENDIF

                     IF(gfinp%l_sphavg) THEN
                        CALL kkintgr(greensfImagPart%sphavg(:,m,mp,i_elem,jspin),eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%gmmpMat(:,m,mp,jspin,ipm),int_method(contourShape))
                     ELSE
                        ! In the case of radial dependence we perform the kramers-kronig-integration seperately for uu,dd,etc.
                        ! We can do this because the radial functions are independent of E
                        CALL kkintgr(greensfImagPart%uu(:,m,mp,i_elem,jspin),eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%uu(:,m,mp,jspin,ipm),int_method(contourShape))
                        CALL kkintgr(greensfImagPart%ud(:,m,mp,i_elem,jspin),eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%ud(:,m,mp,jspin,ipm),int_method(contourShape))
                        CALL kkintgr(greensfImagPart%du(:,m,mp,i_elem,jspin),eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%du(:,m,mp,jspin,ipm),int_method(contourShape))
                        CALL kkintgr(greensfImagPart%dd(:,m,mp,i_elem,jspin),eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%dd(:,m,mp,jspin,ipm),int_method(contourShape))
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL timestop("Green's Function: Kramer-Kronigs-Integration")
      ENDDO

#ifdef CPP_MPI
      CALL timestart("Green's Function: Collect")
      !Collect all the greensFuntions
      DO i_gf = 1, gfinp%n
         CALL g(i_gf)%collect(fmpi%mpi_comm)
      ENDDO
      CALL timestop("Green's Function: Collect")
#endif

   END SUBROUTINE greensfCalcRealPart
END MODULE m_greensfCalcRealPart