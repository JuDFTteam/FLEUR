MODULE m_writeCFOutput

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_intgr
   USE m_mt_tofrom_grid

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE writeCFOutput(atoms,input,sym,sphhar,noco,vTot,hub1data,pot)

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      TYPE(t_input),       INTENT(IN)  :: input
      TYPE(t_sym),         INTENT(IN)  :: sym
      TYPE(t_sphhar),      INTENT(IN)  :: sphhar
      TYPE(t_noco),        INTENT(IN)  :: noco
      TYPE(t_potden),      INTENT(IN)  :: vTot
      TYPE(t_hub1data),    INTENT(IN)  :: hub1data
      LOGICAL, OPTIONAL,   INTENT(IN)  :: pot

      INTEGER, PARAMETER :: lcf = 3

      INTEGER :: iType,l,m,lm,io_error,iGrid,nd, lh,lv
      REAL    :: n_0Norm
      LOGICAL :: processPot

      COMPLEX, ALLOCATABLE :: vlm(:,:,:)
      REAL, ALLOCATABLE :: vTotch(:,:)
      REAL :: n_0(atoms%jmtd)

      TYPE(t_gradients) :: grad

      processPot = .FALSE.
      IF(PRESENT(pot)) processPot = pot

      IF(processPot) THEN
         ALLOCATE(vlm(atoms%jmtd,0:MAXVAL(sphhar%llh)*(MAXVAL(sphhar%llh)+2),input%jspins),source=cmplx_0)
         CALL init_mt_grid(input%jspins, atoms, sphhar, .FALSE., sym, l_mdependency=.TRUE.)
      ENDIF
      DO iType = 1, atoms%ntype

         IF(atoms%l_outputCFcdn(iType)) THEN
            !Calculate n_4f^0(r) (normed spherical part of the 4f charge density)
            n_0 = hub1data%cdn_spherical(:,lcf,iType)
            !Norm to int r^2 n_4f(r) dr = 1
            CALL intgr3(n_0,atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),n_0Norm)
            n_0 = n_0/n_0Norm

            OPEN(unit=29,file='n4f.'//int2str(iType)//'.dat',status='replace',action='write',iostat=io_error)
            IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
            DO iGrid = 1, atoms%jri(iType)
               WRITE(29,'(2e20.8)') atoms%rmsh(iGrid,iType), n_0(iGrid)
            ENDDO
            CLOSE(unit=29,iostat=io_error)
            IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")

         ENDIF

         IF(atoms%l_outputCFpot(iType).AND.processPot) THEN

            !                          sigma
            !Decompose potential into V(r)
            !                          lm
            !IF(ALLOCATED(vTotch)) DEALLOCATE(vTotch)
            !ALLOCATE(vTotch(atoms%nsp()*atoms%jri(iType),input%jspins))
            !CALL mt_to_grid(.FALSE., input%jspins, atoms,sym,sphhar,.True.,vTot%mt(:,0:,iType,:),iType,noco,grad,vTotch)
            !modified mt_from_grid with lm index
            !vlm = cmplx_0
            !CALL mt_from_gridlm(atoms, sym, sphhar, iType, input%jspins, vTotch, vlm)

            !Missing: only write out relevant components

            !DO l = 2, 6, 2
            !   DO m = -l, l
            !      lm = l*(l+1) + m
            nd = sym%ntypsy(SUM(atoms%neq(:iType - 1)) + 1)
            DO lh = 1,sphhar%nlh(nd)
               lv = sphhar%llh(lh,nd)
               IF(sphhar%nmem(lh,nd) == 1) THEN
                  OPEN(unit=29,file='V_'//int2str(lv)//int2str(0)//'.'//int2str(iType)//'.dat',status='replace',&
                       action='write',iostat=io_error)
                  IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
                  DO iGrid = 1, atoms%jri(iType)
                     WRITE(29,'(5e20.8)') atoms%rmsh(iGrid,iType), vTot%mt(iGrid,lh,iType,1), vTot%mt(iGrid,lh,iType,input%jspins)!vlm(iGrid,lm,1), vlm(iGrid,lm,input%jspins)
                  ENDDO
                  CLOSE(unit=29,iostat=io_error)
                  IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
               ENDIF
            ENDDO

         ENDIF

      ENDDO
      IF(processPot) CALL finish_mt_grid()

   END SUBROUTINE writeCFOutput

END MODULE m_writeCFOutput