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

   SUBROUTINE greensfCalcRealPart(atoms,gfinp,sym,input,noco,usdus,denCoeffsOffDiag,fmpi,ef,greensfImagPart,g)

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_usdus),             INTENT(IN)     :: usdus
      TYPE(t_denCoeffsOffDiag),  INTENT(IN)     :: denCoeffsOffDiag
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_mpi),               INTENT(IN)     :: fmpi
      REAL,                      INTENT(IN)     :: ef
      TYPE(t_greensfImagPart),   INTENT(INOUT)  :: greensfImagPart
      TYPE(t_greensf),           INTENT(INOUT)  :: g(:)

      INTEGER :: i_gf,i_elem,l,m,mp,nType,indUnique,nLO,iLO,iLOp,i_elemLO
      INTEGER :: jspin,nspins,ipm,lp,nTypep,refCutoff
      INTEGER :: contourShape
      INTEGER :: i_gf_start,i_gf_end,spin_start,spin_end
      INTEGER :: n_gf_task,extra
      LOGICAL :: l_onsite,l_fixedCutoffset,l_sphavg
      REAL    :: del,eb,fixedCutoff,atomDiff(3)
      REAL,    ALLOCATABLE :: eMesh(:),imag(:)

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
            l_sphavg = g(i_gf)%elem%l_sphavg
            l_fixedCutoffset = g(i_gf)%elem%l_fixedCutoffset
            fixedCutoff      = g(i_gf)%elem%fixedCutoff
            refCutoff        = g(i_gf)%elem%refCutoff
            atomDiff(:) = g(i_gf)%elem%atomDiff(:)

            i_elem = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,indUnique=indUnique)

            IF(i_gf /= indUnique.AND..NOT.l_fixedCutoffset.AND.refCutoff==-1&
               .AND..NOT.g(indUnique)%elem%l_fixedCutoffset.AND.g(indUnique)%elem%refCutoff==-1) THEN
               !This cutoff was already calculated
               greensfImagPart%kkintgr_cutoff(i_gf,:,:) = greensfImagPart%kkintgr_cutoff(indUnique,:,:)
            ELSE
               !Is the current element suitable for automatic finding of the cutoff
               l_onsite = nType.EQ.nTypep.AND.l.EQ.lp.AND.ALL(ABS(atomDiff(:)).LT.1e-12)
               IF(l_onsite.AND..NOT.l_fixedCutoffset.AND.refCutoff==-1 .AND. g(i_gf)%elem%countLOs(atoms)==0) THEN
                  !
                  !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration
                  ! with LOs I just use a fixed cutoff or reference otherwise I would need to check whether
                  ! the LO lies in the energy boundary and raise the expected number of states accordingly
                  IF(l_sphavg) THEN
                     CALL kk_cutoff(greensfImagPart%sphavg(:,:,:,i_elem,:),noco,gfinp%l_mperp,l,input%jspins,&
                                    eMesh,greensfImagPart%kkintgr_cutoff(i_gf,:,:),greensfImagPart%scalingFactorSphavg(i_elem,:))
                  ELSE
                     !Onsite element with radial dependence
                     CALL kk_cutoffRadial(greensfImagPart%uu(:,:,:,i_elem,:),greensfImagPart%ud(:,:,:,i_elem,:),&
                                          greensfImagPart%du(:,:,:,i_elem,:),greensfImagPart%dd(:,:,:,i_elem,:),&
                                          noco,usdus,denCoeffsOffDiag,gfinp%l_mperp,l,nType,input,eMesh,&
                                          greensfImagPart%kkintgr_cutoff(i_gf,:,:),greensfImagPart%scalingFactorRadial(i_elem,:))
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

         ENDDO

         !Getting reference Cutoffs and perform scaling
         DO i_gf = 1, gfinp%n
            l  = g(i_gf)%elem%l
            lp = g(i_gf)%elem%lp
            nType  = g(i_gf)%elem%atomType
            nTypep = g(i_gf)%elem%atomTypep
            l_fixedCutoffset = g(i_gf)%elem%l_fixedCutoffset
            fixedCutoff      = g(i_gf)%elem%fixedCutoff
            refCutoff        = g(i_gf)%elem%refCutoff
            l_sphavg = g(i_gf)%elem%l_sphavg
            i_elem = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,indUnique=indUnique)
            i_elemLO = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,indUnique=indUnique,lo=.TRUE.)
            nLO = g(i_gf)%elem%countLOs(atoms)

            IF(refCutoff/=-1) THEN
               !Overwrite cutoff with reference from other elements
               greensfImagPart%kkintgr_cutoff(i_gf,:,:) = greensfImagPart%kkintgr_cutoff(refCutoff,:,:)
            ENDIF
            CALL greensfImagPart%scale(i_elem,i_elemLO,l_sphavg,nLO)
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
         l_sphavg = g(i_gf)%elem%l_sphavg
         contourShape = gfinp%contour(g(i_gf)%elem%iContour)%shape
         nLO = g(i_gf)%elem%countLOs(atoms)
         IF(g(i_gf)%elem%representative_elem > 0) CYCLE

         i_elem = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg)
         i_elemLO = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,lo=.TRUE.)

         CALL timestart("Green's Function: Kramer-Kronigs-Integration")
         DO jspin = spin_start, spin_end
            DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))
               DO m= -l,l
                  DO mp= -lp,lp

                     IF(greensfImagPart%checkEmpty(i_elem,i_elemLO,nLO,m,mp,jspin,l_sphavg)) THEN
                        CALL g(i_gf)%resetSingleElem(m,mp,jspin,ipm)
                        CYCLE
                     ENDIF

                     IF(l_sphavg) THEN
                        imag = greensfImagPart%applyCutoff(i_elem,i_gf,m,mp,jspin,l_sphavg)
                        CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%gmmpMat(:,m,mp,jspin,ipm),int_method(contourShape))
                     ELSE
                        ! In the case of radial dependence we perform the kramers-kronig-integration seperately for uu,dd,etc.
                        ! We can do this because the radial functions are independent of E
                        imag = greensfImagPart%applyCutoff(i_elem,i_gf,m,mp,jspin,l_sphavg,imat=1)
                        CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%uu(:,m,mp,jspin,ipm),int_method(contourShape))
                        imag = greensfImagPart%applyCutoff(i_elem,i_gf,m,mp,jspin,l_sphavg,imat=2)
                        CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%dd(:,m,mp,jspin,ipm),int_method(contourShape))
                        imag = greensfImagPart%applyCutoff(i_elem,i_gf,m,mp,jspin,l_sphavg,imat=3)
                        CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%ud(:,m,mp,jspin,ipm),int_method(contourShape))
                        imag = greensfImagPart%applyCutoff(i_elem,i_gf,m,mp,jspin,l_sphavg,imat=4)
                        CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                     g(i_gf)%du(:,m,mp,jspin,ipm),int_method(contourShape))

                        !KKT for LOs
                        IF(nLO>0) THEN
                           DO iLO = 1, nLO
                              imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,m,mp,jspin,l_sphavg,imat=1,iLO=iLO)
                              CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                           g(i_gf)%uulo(:,m,mp,iLO,jspin,ipm),int_method(contourShape))
                              imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,m,mp,jspin,l_sphavg,imat=2,iLO=iLO)
                              CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                           g(i_gf)%ulou(:,m,mp,iLO,jspin,ipm),int_method(contourShape))
                              imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,m,mp,jspin,l_sphavg,imat=3,iLO=iLO)
                              CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                           g(i_gf)%dulo(:,m,mp,iLO,jspin,ipm),int_method(contourShape))
                              imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,m,mp,jspin,l_sphavg,imat=4,iLO=iLO)
                              CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                           g(i_gf)%ulod(:,m,mp,iLO,jspin,ipm),int_method(contourShape))

                              DO iLOp = 1, nLO
                                 imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,m,mp,jspin,l_sphavg,iLO=iLO,iLOp=iLop)
                                 CALL kkintgr(imag,eMesh,g(i_gf)%contour%e,(ipm.EQ.2),&
                                              g(i_gf)%uloulop(:,m,mp,iLO,iLOp,jspin,ipm),int_method(contourShape))
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL timestop("Green's Function: Kramer-Kronigs-Integration")
      ENDDO


      !perform rotations for intersite elements
      DO i_gf = i_gf_start, i_gf_end
         IF(i_gf.LT.1 .OR. i_gf.GT.gfinp%n) CYCLE !Make sure to not produce segfaults with mpi
         IF(g(i_gf)%elem%representative_elem <= 0) CYCLE
         g(i_gf) = g(g(i_gf)%elem%representative_elem)
         CALL g(i_gf)%rotate(sym,atoms)
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