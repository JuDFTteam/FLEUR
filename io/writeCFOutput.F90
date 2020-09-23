MODULE m_writeCFOutput

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_mt_tofrom_grid
   USE m_cfOutput_hdf
   USE m_genMTBasis

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE writeCFOutput(atoms,input,sym,sphhar,noco,vTot,hub1data,enpara,fmpi)

      TYPE(t_atoms),       INTENT(IN)  :: atoms
      TYPE(t_input),       INTENT(IN)  :: input
      TYPE(t_sym),         INTENT(IN)  :: sym
      TYPE(t_sphhar),      INTENT(IN)  :: sphhar
      TYPE(t_noco),        INTENT(IN)  :: noco
      TYPE(t_potden),      INTENT(IN)  :: vTot
      TYPE(t_hub1data),    INTENT(IN)  :: hub1data
      TYPE(t_enpara),      INTENT(IN)  :: enpara
      TYPE(t_mpi),         INTENT(IN)  :: fmpi

      INTEGER, PARAMETER :: lcf = 3

      INTEGER :: iType,l,m,lm,io_error,iGrid,ispin
      LOGICAL :: processPot

      COMPLEX, ALLOCATABLE :: vlm(:,:,:)
      REAL, ALLOCATABLE :: vTotch(:,:)
      REAL, ALLOCATABLE :: f(:,:,:),g(:,:,:),flo(:,:,:)
      REAL :: n_0(atoms%jmtd)

#ifdef CPP_HDF
      INTEGER(HID_T) :: cfFileID
#endif

      TYPE(t_gradients) :: grad
      TYPE(t_potden)    :: vTotProcess
      TYPE(t_usdus)     :: usdus

      vTotProcess = vTot
      ALLOCATE(vlm(atoms%jmtd,0:MAXVAL(sphhar%llh)*(MAXVAL(sphhar%llh)+2),input%jspins),source=cmplx_0)
      CALL init_mt_grid(input%jspins, atoms, sphhar, .FALSE., sym, l_mdependency=.TRUE.)

      ALLOCATE (f(atoms%jmtd,2,0:atoms%lmaxd),source=0.0)
      ALLOCATE (g(atoms%jmtd,2,0:atoms%lmaxd),source=0.0)
      ALLOCATE (flo(atoms%jmtd,2,atoms%nlod),source=0.0)
      CALL usdus%init(atoms,input%jspins)

#ifdef CPP_HDF
      CALL opencfFile(cfFileID, atoms, l_create = .TRUE.)
#endif
      DO iType = 1, atoms%ntype

         IF(atoms%l_outputCFcdn(iType)) THEN
            IF(vTot%potdenType.EQ.POTDEN_TYPE_CRYSTALFIELD) THEN
               CALL juDFT_error("Simultaneous calculation of cf potential and charge density not supported yet",&
                                calledby="writeCFOutput")
            ENDIF
            n_0 = 0.0
            DO ispin = 1, input%jspins
               CALL genMTBasis(atoms,enpara,vTot,fmpi,iType,ispin,usdus,f,g,flo,hub1data,.FALSE.)
               n_0(:) = n_0(:) + f(:,1,lcf)*f(:,1,lcf) + f(:,2,lcf)*f(:,2,lcf)
            ENDDO

#ifdef CPP_HDF
            CALL writeCFcdn(cfFileID, atoms, iType, n_0)
#else
            !Stupid text output
            OPEN(unit=29,file='n4f.'//int2str(iType)//'.dat',status='replace',action='write',iostat=io_error)
            IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
            DO iGrid = 1, atoms%jri(iType)
               WRITE(29,'(2e20.8)') atoms%rmsh(iGrid,iType), n_0(iGrid)
            ENDDO
            CLOSE(unit=29,iostat=io_error)
            IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
#endif

         ENDIF

         IF(atoms%l_outputCFpot(iType).AND.processPot) THEN
            IF(vTot%potdenType.NE.POTDEN_TYPE_CRYSTALFIELD) THEN
               CALL juDFT_error("Wrong potential type for crystalfield",calledby="writeCFOutput")
            ENDIF
            !                          sigma
            !Decompose potential into V(r)
            !                          lm
            DO ispin =1, input%jspins
               DO iGrid=1,atoms%jri(iType)
                  vTotProcess%mt(iGrid,:,iType,ispin)=vTotProcess%mt(iGrid,:,iType,ispin)*atoms%rmsh(iGrid,iType)**2
               END DO
            ENDDO
            IF(ALLOCATED(vTotch)) DEALLOCATE(vTotch)
            ALLOCATE(vTotch(atoms%nsp()*atoms%jri(iType),input%jspins))
            CALL mt_to_grid(.FALSE., input%jspins, atoms,sym,sphhar,.True.,vTotProcess%mt(:,0:,iType,:),iType,noco,grad,vTotch)
            !modified mt_from_grid with lm index
            vlm = cmplx_0
            CALL mt_from_gridlm(atoms, sym, sphhar, iType, input%jspins, vTotch, vlm)

            !Missing: only write out relevant components
#ifdef CPP_HDF
            CALL writeCFpot(cfFileID, atoms, input, iType, vlm)
#else
            !Stupid text output
            DO l = 2, 6, 2
               DO m = -l, l
                  lm = l*(l+1) + m
                  OPEN(unit=29,file='V_'//int2str(l)//int2str(m)//'.'//int2str(iType)//'.dat',status='replace',&
                       action='write',iostat=io_error)
                  IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
                  DO iGrid = 1, atoms%jri(iType)
                     WRITE(29,'(5e20.8)') atoms%rmsh(iGrid,iType), vlm(iGrid,lm,1), vlm(iGrid,lm,input%jspins)
                  ENDDO
                  CLOSE(unit=29,iostat=io_error)
                  IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
               ENDDO
            ENDDO
#endif
         ENDIF

      ENDDO

#ifdef CPP_HDF
      CALL closecfFile(cfFileID)
#endif
      IF(processPot) CALL finish_mt_grid()

   END SUBROUTINE writeCFOutput

END MODULE m_writeCFOutput