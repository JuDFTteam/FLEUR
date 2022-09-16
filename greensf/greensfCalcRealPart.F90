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

   private
   public greensfCalcRealPart

   CONTAINS

   SUBROUTINE greensfCalcRealPart(atoms,gfinp,sym,input,enpara,noco,kpts,fmpi,ef,greensfImagPart,g)

      TYPE(t_atoms),             INTENT(IN)     :: atoms
      TYPE(t_gfinp),             INTENT(IN)     :: gfinp
      TYPE(t_sym),               INTENT(IN)     :: sym
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_kpts),              INTENT(IN)     :: kpts
      TYPE(t_input),             INTENT(IN)     :: input
      TYPE(t_mpi),               INTENT(IN)     :: fmpi
      type(t_enpara),            intent(in)     :: enpara
      REAL,                      INTENT(IN)     :: ef
      TYPE(t_greensfImagPart),   INTENT(INOUT)  :: greensfImagPart
      TYPE(t_greensf),           INTENT(INOUT)  :: g(:)

      INTEGER :: i_gf,i_elem,l,m,mp,indUnique,nLO,iLO,iLOp,i_elemLO
      INTEGER :: jspin,nspins,ipm,lp,refCutoff
      INTEGER :: contourShape, iContour
      INTEGER :: i_gf_start,i_gf_end,spin_start,spin_end
      INTEGER :: ikpt, ikpt_i, refElem
      LOGICAL :: l_fixedCutoffset,l_sphavg,l_kresolved_int,l_kresolved
      REAL    :: del,eb,fixedCutoff,bk(3)
      REAL,    ALLOCATABLE :: eMesh(:)
      COMPLEX, ALLOCATABLE :: gmat(:,:,:),imag(:,:,:)

      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(ef,del,eb,eMesh=eMesh)

      nspins = MERGE(3,input%jspins,gfinp%l_mperp)

      IF(fmpi%irank.EQ.0) THEN
         CALL timestart("Green's Function: Integration Cutoff")
         DO i_gf = 1, gfinp%n

            !Get the information of ith current element
            l  = g(i_gf)%elem%l
            lp = g(i_gf)%elem%lp
            l_sphavg = g(i_gf)%elem%l_sphavg
            l_fixedCutoffset = g(i_gf)%elem%l_fixedCutoffset
            fixedCutoff      = g(i_gf)%elem%fixedCutoff
            refCutoff        = g(i_gf)%elem%refCutoff
            l_kresolved_int = g(i_gf)%elem%l_kresolved_int

            IF(refCutoff /= -1) CYCLE

            IF(l_fixedCutoffset) THEN
               greensfImagPart%kkintgr_cutoff(i_gf,:,1) = 1
               greensfImagPart%kkintgr_cutoff(i_gf,:,2) = INT((fixedCutoff+ef-eb)/del)+1
               CYCLE
            ENDIF

            IF(.NOT.gfinp%isUnique(i_gf,distinct_kresolved_int=.TRUE.)) THEN
               indUnique = gfinp%getUniqueElement(i_gf,distinct_kresolved_int=.TRUE.)
               !This cutoff was already calculated
               greensfImagPart%kkintgr_cutoff(i_gf,:,:) = greensfImagPart%kkintgr_cutoff(indUnique,:,:)
            ELSE
               i_elem = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,l_kresolved_int=l_kresolved_int)
               IF(.not. g(i_gf)%elem%isOffDiag() .and. &
                  .not. has_sclos(g(i_gf)%elem, atoms,enpara) .and. &
                  .not. l_kresolved_int) THEN
                  !
                  !Check the integral over the PDOS to define a cutoff for the Kramer-Kronigs-Integration
                  ! with SCLOs I just use a fixed cutoff or reference otherwise I would need to check whether
                  ! the SCLO lies in the energy boundary and raise the expected number of states accordingly
                  IF(l_sphavg) THEN
                     CALL kk_cutoff(greensfImagPart%sphavg(:,:,:,i_elem,:),noco,gfinp%l_mperp,l,input%jspins,&
                                    eMesh,greensfImagPart%kkintgr_cutoff(i_gf,:,:),greensfImagPart%scalingFactorSphavg(i_elem,:))
                  ELSE IF(g(i_gf)%elem%countLOs(atoms)==0) then
                     !Onsite element with radial dependence
                     CALL kk_cutoffRadial(greensfImagPart%uu(:,:,:,i_elem,:),greensfImagPart%ud(:,:,:,i_elem,:),&
                                          greensfImagPart%du(:,:,:,i_elem,:),greensfImagPart%dd(:,:,:,i_elem,:),&
                                          noco,g(i_gf)%scalarProducts,gfinp%l_mperp,l,input,eMesh,&
                                          greensfImagPart%kkintgr_cutoff(i_gf,:,:),greensfImagPart%scalingFactorRadial(i_elem,:))
                  ELSE
                     CALL kk_cutoffRadialLO(greensfImagPart%uu(:,:,:,i_elem,:),greensfImagPart%ud(:,:,:,i_elem,:),&
                                          greensfImagPart%du(:,:,:,i_elem,:),greensfImagPart%dd(:,:,:,i_elem,:),&
                                          greensfImagPart%uulo(:,:,:,:,i_elem,:),greensfImagPart%dulo(:,:,:,:,i_elem,:),&
                                          greensfImagPart%ulou(:,:,:,:,i_elem,:),greensfImagPart%ulod(:,:,:,:,i_elem,:),&
                                          greensfImagPart%uloulop(:,:,:,:,:,i_elem,:),&
                                          atoms,noco,g(i_gf)%scalarProducts,gfinp%l_mperp,l,g(i_gf)%elem%atomType,input,eMesh,&
                                          greensfImagPart%kkintgr_cutoff(i_gf,:,:),greensfImagPart%scalingFactorRadial(i_elem,:))
                  ENDIF
               ELSE
                  !For all other elements we just use ef+elup as a hard cutoff
                  greensfImagPart%kkintgr_cutoff(i_gf,:,1) = 1
                  greensfImagPart%kkintgr_cutoff(i_gf,:,2) = gfinp%ne
               ENDIF
            ENDIF

         ENDDO

         !Getting reference Cutoffs and perform scaling
         DO i_gf = 1, gfinp%n
            refCutoff       = g(i_gf)%elem%refCutoff
            l_kresolved_int = g(i_gf)%elem%l_kresolved_int
            l_sphavg = g(i_gf)%elem%l_sphavg
            i_elem = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,l_kresolved_int=l_kresolved_int)
            i_elemLO = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,l_kresolved_int=l_kresolved_int,lo=.TRUE.)
            nLO = g(i_gf)%elem%countLOs(atoms)

            IF(refCutoff/=-1) THEN
               !Overwrite cutoff with reference from other elements
               greensfImagPart%kkintgr_cutoff(i_gf,:,:) = greensfImagPart%kkintgr_cutoff(refCutoff,:,:)
               
               ! refElem = gfinp%getUniqueElement(refCutoff,distinct_kresolved_int=.TRUE.)
               ! if (l_sphavg .and. .not. l_kresolved_int) then
               !    greensfImagPart%scalingFactorSphavg(i_elem,:) = greensfImagPart%scalingFactorSphavg(refElem,:)
               ! else if(.not.l_sphavg) then
               !    greensfImagPart%scalingFactorRadial(i_elem,:) = greensfImagPart%scalingFactorRadial(refElem,:)
               ! else
               !    greensfImagPart%scalingFactorSphavgKres(i_elem,:) = greensfImagPart%scalingFactorSphavgKres(refElem,:)
               ! endif
            ENDIF
            if(.NOT.gfinp%isUnique(i_gf,distinct_kresolved_int=.TRUE.)) cycle
            CALL greensfImagPart%scale(i_elem,i_elemLO,l_sphavg,nLO,k_resolved=l_kresolved_int)
         ENDDO
         CALL timestop("Green's Function: Integration Cutoff")
      ENDIF

      !Broadcast cutoffs and modified imaginary parts
      CALL greensfImagPart%mpi_bc(fmpi%mpi_comm)

      !Distribute the Calculations
      CALL gfinp%distribute_elements(fmpi%irank, fmpi%isize, nspins, i_gf_start, i_gf_end, spin_start, spin_end)

      !Initialize kkintgr_module variables
      DO i_gf = i_gf_start, i_gf_end
         IF(i_gf.LT.1 .OR. i_gf.GT.gfinp%n) CYCLE !Make sure to not produce segfaults with mpi
         contourShape = gfinp%contour(g(i_gf)%elem%iContour)%shape

         CALL kkintgr_init(eMesh,g(i_gf)%contour%e,g(i_gf)%elem%iContour,gfinp%numberContours, contourShape, gfinp%additional_smearing)
      ENDDO

      DO i_gf = i_gf_start, i_gf_end

         IF(i_gf.LT.1 .OR. i_gf.GT.gfinp%n) CYCLE !Make sure to not produce segfaults with mpi

         !Get the information of ith current element
         l  = g(i_gf)%elem%l
         lp = g(i_gf)%elem%lp
         l_sphavg = g(i_gf)%elem%l_sphavg
         iContour = g(i_gf)%elem%iContour
         nLO = g(i_gf)%elem%countLOs(atoms)
         IF(g(i_gf)%elem%representative_elem > 0) CYCLE
         IF(g(i_gf)%elem%l_kresolved_int) CYCLE

         i_elem = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,l_kresolved_int=.FALSE.)
         i_elemLO = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,lo=.TRUE.,l_kresolved_int=.FALSE.)

         CALL timestart("Green's Function: Kramer-Kronigs-Integration")
         DO jspin = spin_start, spin_end
            DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))

               IF(l_sphavg) THEN
                  imag = greensfImagPart%applyCutoff(i_elem,i_gf,jspin,l_sphavg)
                  CALL kkintgr(imag,ipm==2,g(i_gf)%gmmpMat(:,:,:,jspin,ipm),iContour)
               ELSE
                  ! In the case of radial dependence we perform the kramers-kronig-integration seperately for uu,dd,etc.
                  ! We can do this because the radial functions are independent of E
                  imag = greensfImagPart%applyCutoff(i_elem,i_gf,jspin,l_sphavg,imat=1)
                  CALL kkintgr(imag,ipm==2,g(i_gf)%uu(:,:,:,jspin,ipm),iContour)
                  imag = greensfImagPart%applyCutoff(i_elem,i_gf,jspin,l_sphavg,imat=2)
                  CALL kkintgr(imag,ipm==2,g(i_gf)%dd(:,:,:,jspin,ipm),iContour)
                  imag = greensfImagPart%applyCutoff(i_elem,i_gf,jspin,l_sphavg,imat=3)
                  CALL kkintgr(imag,ipm==2,g(i_gf)%ud(:,:,:,jspin,ipm),iContour)
                  imag = greensfImagPart%applyCutoff(i_elem,i_gf,jspin,l_sphavg,imat=4)
                  CALL kkintgr(imag,ipm==2,g(i_gf)%du(:,:,:,jspin,ipm),iContour)

                  !KKT for LOs
                  IF(nLO>0) THEN
                     DO iLO = 1, nLO
                        imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,jspin,l_sphavg,imat=1,iLO=iLO)
                        CALL kkintgr(imag,ipm==2,g(i_gf)%uulo(:,:,:,iLO,jspin,ipm),iContour)
                        imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,jspin,l_sphavg,imat=2,iLO=iLO)
                        CALL kkintgr(imag,ipm==2,g(i_gf)%ulou(:,:,:,iLO,jspin,ipm),iContour)
                        imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,jspin,l_sphavg,imat=3,iLO=iLO)
                        CALL kkintgr(imag,ipm==2,g(i_gf)%dulo(:,:,:,iLO,jspin,ipm),iContour)
                        imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,jspin,l_sphavg,imat=4,iLO=iLO)
                        CALL kkintgr(imag,ipm==2,g(i_gf)%ulod(:,:,:,iLO,jspin,ipm),iContour)

                        DO iLOp = 1, nLO
                           imag = greensfImagPart%applyCutoff(i_elemLO,i_gf,jspin,l_sphavg,iLO=iLO,iLOp=iLop)
                           CALL kkintgr(imag,ipm==2,g(i_gf)%uloulop(:,:,:,iLO,iLOp,jspin,ipm),iContour)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         CALL timestop("Green's Function: Kramer-Kronigs-Integration")
      ENDDO
      CALL kkintgr_free()

      IF(ANY(gfinp%elem(:)%l_kresolved_int)) THEN
         CALL gfinp%distribute_elements(fmpi%n_rank, fmpi%n_size, nspins, i_gf_start, i_gf_end, spin_start, spin_end, k_resolved=.TRUE.)
         !Initialize kkintgr_module variables
         DO i_gf = i_gf_start, i_gf_end
            IF(i_gf.LT.1 .OR. i_gf.GT.gfinp%n) CYCLE !Make sure to not produce segfaults with mpi
            contourShape = gfinp%contour(g(i_gf)%elem%iContour)%shape
            CALL kkintgr_init(eMesh,g(i_gf)%contour%e,g(i_gf)%elem%iContour,gfinp%numberContours, contourShape, gfinp%additional_smearing)
         ENDDO
         CALL timestart("Green's Function: K-Resolved Kramer-Kronigs-Integration")
         DO ikpt_i = 1, SIZE(fmpi%k_list)
            ikpt = fmpi%k_list(ikpt_i)
            bk = kpts%bk(:,ikpt)
            DO i_gf = i_gf_start, i_gf_end

               IF(i_gf.LT.1 .OR. i_gf.GT.gfinp%n) CYCLE !Make sure to not produce segfaults with mpi

               !Get the information of ith current element
               l  = g(i_gf)%elem%l
               lp = g(i_gf)%elem%lp
               l_sphavg = g(i_gf)%elem%l_sphavg
               l_kresolved = g(i_gf)%elem%l_kresolved
               iContour = g(i_gf)%elem%iContour
               nLO = g(i_gf)%elem%countLOs(atoms)
               IF(g(i_gf)%elem%representative_elem > 0) CYCLE
               IF(.NOT.g(i_gf)%elem%l_kresolved_int) CYCLE

               i_elem = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,l_kresolved_int=.TRUE.)
               i_elemLO = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,lo=.TRUE.,l_kresolved_int=.TRUE.)

               IF(.NOT.l_kresolved) ALLOCATE(gmat(SIZE(g(i_gf)%contour%e),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const), source=cmplx_0)

               CALL timestart("Green's Function: Kramer-Kronigs-Integration")
               DO jspin = spin_start, spin_end
                  DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))

                     IF(l_sphavg) THEN
                        imag = greensfImagPart%applyCutoff(i_elem,i_gf,jspin,l_sphavg,ikpt=ikpt_i)
                        IF(l_kresolved) THEN
                           CALL kkintgr(imag,ipm==2,g(i_gf)%gmmpMat_k(:,:,:,jspin,ipm,ikpt),iContour)
                        ELSE
                           CALL kkintgr(imag,ipm==2,gmat,iContour)
                           g(i_gf)%gmmpMat(:,:,:,jspin,ipm) = g(i_gf)%gmmpMat(:,:,:,jspin,ipm) + gmat
                        ENDIF
                     ELSE
                        CALL juDFT_error("No Green's function with k-resolution and radial dependence implemented")
                     ENDIF

                  ENDDO
               ENDDO
               CALL timestop("Green's Function: Kramer-Kronigs-Integration")
               IF(.NOT.l_kresolved) DEALLOCATE(gmat)
            ENDDO
         ENDDO
         CALL timestop("Green's Function: K-Resolved Kramer-Kronigs-Integration")
         CALL kkintgr_free()
      ENDIF


#ifdef CPP_MPI
      CALL timestart("Green's Function: Collect")
      !Collect all the greensFuntions
      DO i_gf = 1, gfinp%n
         CALL g(i_gf)%collect(fmpi%mpi_comm)
      ENDDO
      CALL timestop("Green's Function: Collect")
#endif

      IF(fmpi%irank.EQ.0) THEN
         !perform rotations for intersite elements
         DO i_gf = 1, gfinp%n
            IF(g(i_gf)%elem%representative_elem <= 0) CYCLE
            CALL g(i_gf)%set_gfdata(g(g(i_gf)%elem%representative_elem))
            CALL g(i_gf)%rotate(sym,atoms)
         ENDDO
      ENDIF

      DO i_gf = 1, gfinp%n
         CALL g(i_gf)%mpi_bc(fmpi%mpi_comm)
      ENDDO

   END SUBROUTINE greensfCalcRealPart

   logical function has_sclos(elem,atoms, enpara)

      type(t_gfelementtype), intent(in) :: elem
      type(t_atoms),         intent(in) :: atoms
      type(t_enpara),        intent(in) :: enpara

      integer :: ilo

      has_sclos = .false.
      DO ilo = 1, atoms%nlo(elem%atomType)
         if(atoms%llo(ilo,elem%atomType).NE.elem%l) cycle
         !Check for HELO (negative energy parameter) and HDLO (eDeriv > 0)
         if(all(enpara%qn_ello(ilo,elem%atomType,:)<0).or.atoms%ulo_der(ilo, elem%atomType)>=1) cycle
         has_sclos = .true.
      ENDDO

      IF(elem%l.NE.elem%lp.OR.elem%atomType.NE.elem%atomTypep) THEN
         DO ilo = 1, atoms%nlo(elem%atomTypep)
            if(atoms%llo(ilo,elem%atomTypep).NE.elem%lp) cycle
            !Check for HELO (negative energy parameter) and HDLO (eDeriv > 0)
            if(all(enpara%qn_ello(ilo,elem%atomTypep,:)<0).or.atoms%ulo_der(ilo, elem%atomTypep)>=1) cycle
            has_sclos = .true.
         ENDDO
      ENDIF

   end function

END MODULE m_greensfCalcRealPart
