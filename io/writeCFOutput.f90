MODULE m_writeCFOutput

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_intgr
   USE m_mt_tofrom_grid

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE writeCFOutput(atoms,input,sym,sphhar,noco,vTot,hub1data)

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      TYPE(t_input),       INTENT(IN)  :: input
      TYPE(t_sym),         INTENT(IN)  :: sym
      TYPE(t_sphhar),      INTENT(IN)  :: sphhar
      TYPE(t_noco),        INTENT(IN)  :: noco
      TYPE(t_potden),      INTENT(IN)  :: vTot
      TYPE(t_hub1data),    INTENT(IN)  :: hub1data

      INTEGER, PARAMETER :: lcf = 3

      INTEGER :: iType,l,m,lm,io_error,iGrid
      REAL    :: n_0Norm

      COMPLEX, ALLOCATABLE :: vlm(:,:,:)
      REAL, ALLOCATABLE :: vTotch(:,:)
      REAL :: n_0(atoms%jmtd)

      TYPE(t_gradients) :: grad

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
               WRITE(29,'(2e14.8)') atoms%rmsh(iGrid,iType), n_0(iGrid)
            ENDDO
            CLOSE(unit=29,iostat=io_error)
            IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")

         ENDIF

         IF(atoms%l_outputCFpot(iType)) THEN

            !                          sigma
            !Decompose potential into V(r)
            !                          lm
            IF(ALLOCATED(vTotch)) DEALLOCATE(vTotch)
            ALLOCATE(vTotch(atoms%nsp()*atoms%jri(iType),input%jspins))
            CALL mt_to_grid(.FALSE., input%jspins, atoms,sym,sphhar,.True.,vTot%mt(:,0:,iType,:),iType,noco,grad,vTotch)
            !modified mt_from_grid with lm index
            vlm = cmplx_0
            CALL mt_from_gridlm(atoms, sym, sphhar, iType, input%jspins, vTotch, vlm)

            !Missing: only write out relevant components

            DO l = 2, 6, 2
               DO m = -l, l
                  lm = l*(l+1) + m
                  OPEN(unit=29,file='V_'//int2str(l)//int2str(m)//'.'//int2str(iType)//'.dat',status='replace',&
                       action='write',iostat=io_error)
                  IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
                  DO iGrid = 1, atoms%jri(iType)
                     WRITE(29,'(5e14.8)') atoms%rmsh(iGrid,iType), vlm(iGrid,lm,1), vlm(iGrid,lm,input%jspins)
                  ENDDO
                  CLOSE(unit=29,iostat=io_error)
                  IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
               ENDDO
            ENDDO

         ENDIF

      ENDDO

   END SUBROUTINE writeCFOutput

END MODULE m_writeCFOutput