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

#ifdef CPP_MPI
   INCLUDE 'mpif.h'
#endif

   INTEGER, PARAMETER :: int_method(3) = [method_direct,method_direct,method_maclaurin]

   CONTAINS

   SUBROUTINE greensfCalcRealPart(atoms,gfinp,input,sym,noco,mpi,ef,greensfImagPart,g)

      TYPE(t_atoms),          INTENT(IN)     :: atoms
      TYPE(t_gfinp),          INTENT(IN)     :: gfinp
      TYPE(t_sym),            INTENT(IN)     :: sym
      TYPE(t_noco),           INTENT(IN)     :: noco
      TYPE(t_input),          INTENT(IN)     :: input
      TYPE(t_mpi),            INTENT(IN)     :: mpi
      REAL,                   INTENT(IN)     :: ef
      TYPE(t_greensfImagPart),INTENT(INOUT)  :: greensfImagPart
      TYPE(t_greensf),        INTENT(INOUT)  :: g(:)

      INTEGER :: i_gf,i_elem,ie,l,m,mp,nType
      INTEGER :: jspin,nspins,ipm,kkcut,lp,nTypep
      INTEGER :: spin_cut,nn,natom,contourShape,dummy
      INTEGER :: i_gf_start,i_gf_end,spin_start,spin_end
      INTEGER :: n_gf_task,extra,n,ierr
      LOGICAL :: l_onsite,l_fixedCutoffset,l_skip
      REAL    :: fac,del,eb,et,fixedCutoff
      REAL,    ALLOCATABLE :: eMesh(:)
      INTEGER, ALLOCATABLE :: itmp(:)

      !Get the information on the real axis energy mesh
      CALL gfinp%eMesh(ef,del,eb,eMesh=eMesh)

      nspins = MERGE(3,input%jspins,gfinp%l_mperp)

      !Distribute the Calculations
#ifdef CPP_MPI
      IF(mpi%isize>1) THEN
         IF(gfinp%n>=mpi%isize) THEN
            !Just distribute the individual gf elements over the ranks
            n_gf_task = FLOOR(REAL(gfinp%n)/(mpi%isize))
            extra = gfinp%n - n_gf_task*mpi%isize
            i_gf_start = mpi%irank*n_gf_task + 1 + extra
            i_gf_end = (mpi%irank+1)*n_gf_task   + extra
            IF(mpi%irank < extra) THEN
               i_gf_start = i_gf_start - (extra - mpi%irank)
               i_gf_end = i_gf_end - (extra - mpi%irank - 1)
            ENDIF
            spin_start = 1
            spin_end   = nspins
         ELSE IF(gfinp%n*nspins>mpi%isize) THEN
            !Just fill up the ranks
            i_gf_start = mpi%irank + 1
            i_gf_end   = mpi%irank + 1
            spin_start = 1
            spin_end   = nspins
         ELSE
            !If there are few enough gf elements then distribute the spins
            spin_start = MOD(mpi%irank,nspins) + 1
            spin_end   = MOD(mpi%irank,nspins) + 1
            i_gf_start = 1 + FLOOR(REAL(mpi%irank)/nspins)
            i_gf_end   = 1 + FLOOR(REAL(mpi%irank)/nspins)
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
         contourShape     = gfinp%contour(g(i_gf)%elem%iContour)%shape
         l_fixedCutoffset = g(i_gf)%elem%l_fixedCutoffset
         fixedCutoff      = g(i_gf)%elem%fixedCutoff

         CALL uniqueElements_gfinp(gfinp,dummy,ind=i_gf,indUnique=i_elem)

         !Is the current element suitable for automatic finding of the cutoff
         l_onsite = nType.EQ.nTypep.AND.l.EQ.lp.AND.gfinp%l_sphavg
         IF(l_onsite.AND..NOT.l_fixedCutoffset) THEN
            !
            !Check the integral over the fDOS to define a cutoff for the Kramer-Kronigs-Integration
            !
            CALL timestart("On-Site: Integration Cutoff")
            CALL kk_cutoff(greensfImagPart%sphavg(:,:,:,i_elem,:),noco,gfinp%l_mperp,l,input%jspins,&
                           eMesh,greensfImagPart%kkintgr_cutoff(i_gf,:,:))
            CALL timestop("On-Site: Integration Cutoff")
         ELSE IF (l_fixedCutoffset) THEN
            greensfImagPart%kkintgr_cutoff(i_gf,:,1) = 1
            greensfImagPart%kkintgr_cutoff(i_gf,:,2) = INT((fixedCutoff+ef-eb)/del)+1
         ELSE
            !For all other elements we just use ef+elup as a hard cutoff
            greensfImagPart%kkintgr_cutoff(i_gf,:,1) = 1
            greensfImagPart%kkintgr_cutoff(i_gf,:,2) = gfinp%ne
         ENDIF
         !
         !Perform the Kramers-Kronig-Integration if not already calculated
         !
         CALL timestart("Green's Function: Kramer-Kronigs-Integration")
         DO jspin = spin_start, spin_end
            spin_cut = MERGE(1,jspin,jspin.GT.2)
            kkcut = greensfImagPart%kkintgr_cutoff(i_gf,spin_cut,2)
            !------------------------------------------------------------
            ! Set everything above the cutoff in the imaginary part to 0
            ! We do this explicitely because when we just use the hard cutoff index
            ! Things might get lost when the imaginary part is smoothed explicitely
            !------------------------------------------------------------
            IF(kkcut.ne.SIZE(eMesh)) THEN
               greensfImagPart%sphavg(kkcut+1:,-l:l,-l:l,i_elem,jspin) = 0.0
            ENDIF
            DO ipm = 1, 2 !upper or lower half of the complex plane (G(E \pm i delta))
               DO m= -l,l
                  DO mp= -lp,lp

                     !Don't waste time on empty elements
                     l_skip = .FALSE.
                     DO ie = 1, SIZE(eMesh)
                        IF(ABS(greensfImagPart%sphavg(ie,m,mp,i_elem,jspin)).GT.1e-12) EXIT
                        IF(ie==SIZE(eMesh)) l_skip = .TRUE.
                     ENDDO
                     IF(l_skip) THEN
                        g(i_gf)%gmmpMat(:,m,mp,jspin,ipm) = cmplx_0
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

      !Collect all the greensFuntions
      DO i_gf = 1, gfinp%n
         CALL g(i_gf)%collect(mpi%mpi_comm)
      ENDDO

#ifdef CPP_MPI
      !Collect all cutoffs
      n = SIZE(greensfImagPart%kkintgr_cutoff)
      ALLOCATE(itmp(n))
      CALL MPI_REDUCE(greensfImagPart%kkintgr_cutoff,itmp,n,MPI_INTEGER,MPI_SUM,0,mpi%mpi_comm,ierr)
      if (mpi%irank==0) greensfImagPart%kkintgr_cutoff = reshape(itmp,shape(greensfImagPart%kkintgr_cutoff))
      DEALLOCATE(itmp)
#endif

      IF(mpi%irank.EQ.0) THEN
         DO i_gf = 1, gfinp%n
            l  = g(i_gf)%elem%l
            CALL uniqueElements_gfinp(gfinp,dummy,ind=i_gf,indUnique=i_elem)

            DO jspin = 1, nspins
               spin_cut = MERGE(1,jspin,jspin.GT.2)
               kkcut = greensfImagPart%kkintgr_cutoff(i_gf,spin_cut,2)
               !------------------------------------------------------------
               ! Set everything above the cutoff in the imaginary part to 0
               ! We do this explicitely because when we just use the hard cutoff index
               ! Things might get lost when the imaginary part is smoothed explicitely
               !------------------------------------------------------------
               IF(kkcut.ne.SIZE(eMesh)) THEN
                  greensfImagPart%sphavg(kkcut+1:,-l:l,-l:l,i_elem,jspin) = 0.0
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      CALL greensfImagPart%mpi_bc(mpi%mpi_comm)

   END SUBROUTINE greensfCalcRealPart
END MODULE m_greensfCalcRealPart