MODULE m_writeCFOutput

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_lattHarmsSphHarmsConv
   USE m_cfOutput_hdf
   USE m_vgen
   USE m_intgr
   USE m_mpi_bc_tool

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE writeCFOutput(fi,stars,hybdat,sphhar,xcpot,EnergyDen,inDen,hub1data,nococonv,enpara,fmpi)

      TYPE(t_fleurinput),  INTENT(IN)  :: fi
      TYPE(t_stars),       INTENT(IN)  :: stars
      TYPE(t_hybdat),      INTENT(IN)  :: hybdat
      TYPE(t_sphhar),      INTENT(IN)  :: sphhar
      CLASS(t_xcpot),      INTENT(IN)  :: xcpot
      TYPE(t_potden),      INTENT(IN)  :: EnergyDen
      TYPE(t_potden),      INTENT(IN)  :: inDen
      TYPE(t_hub1data),    INTENT(IN)  :: hub1data
      TYPE(t_nococonv),    INTENT(IN)  :: nococonv
      TYPE(t_enpara),      INTENT(IN)  :: enpara
      TYPE(t_mpi),         INTENT(IN)  :: fmpi

      INTEGER, PARAMETER :: lcf = 3
#ifdef CPP_HDF
      INTEGER(HID_T) :: cfFileID
#endif

      INTEGER :: iType,l,m,lm,io_error,iGrid,ispin
      REAL    :: n_0Norm
      COMPLEX, ALLOCATABLE :: vlm(:,:,:)
      REAL,    ALLOCATABLE :: f(:,:,:),g(:,:,:),flo(:,:,:)
      REAL :: n_0(fi%atoms%jmtd)

      !Dummy variables to avoid accidental changes to them in vgen
      TYPE(t_results)   :: results_dummy
      TYPE(t_nococonv)  :: nococonv_dummy
      TYPE(t_atoms)     :: atoms_dummy

      !Modified densities and potentials for crystalfield
      TYPE(t_potden)    :: inDenCF
      TYPE(t_potden)    :: vCF,vCoul,vx,vxc,exc

      CALL timestart("Crystal Field Output")

      ALLOCATE(vlm(fi%atoms%jmtd,fi%atoms%lmaxd*(fi%atoms%lmaxd+2)+1,fi%input%jspins),source=cmplx_0)

      !POTDEN_TYPE_CRYSTALFIELD excludes the external potential in the coulomb potential
      CALL vCF%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL vCoul%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_CRYSTALFIELD)
      CALL vx%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTCOUL)
      CALL vxc%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT)
      CALL exc%init(stars, fi%atoms, sphhar, fi%vacuum, fi%noco, fi%input%jspins, POTDEN_TYPE_POTTOT)

      CALL results_dummy%init(fi%input,fi%atoms,fi%kpts,fi%noco)

      ALLOCATE (f(fi%atoms%jmtd,2,0:fi%atoms%lmaxd),source=0.0)
      ALLOCATE (g(fi%atoms%jmtd,2,0:fi%atoms%lmaxd),source=0.0)
      ALLOCATE (flo(fi%atoms%jmtd,2,fi%atoms%nlod),source=0.0)

#ifdef CPP_HDF
      IF(fmpi%irank==0) CALL opencfFile(cfFileID, fi%atoms, fi%cell, l_create = .TRUE.)
#endif
      DO iType = 1, fi%atoms%ntype

         IF(fi%atoms%l_outputCFcdn(iType)) THEN
            n_0 = 0.0
            DO ispin = 1, fi%input%jspins
               n_0(:) = n_0(:) + hub1data%cdn_atomic(:,lcf,iType,ispin)
            ENDDO
            CALL intgr3(n_0,fi%atoms%rmsh(:,iType),fi%atoms%dx(iType),fi%atoms%jri(iType),n_0Norm)
            n_0 = n_0/n_0Norm

            IF(fmpi%irank==0) THEN
#ifdef CPP_HDF
               CALL writeCFcdn(cfFileID, fi%atoms, iType, n_0)
#else
               !Stupid text output
               OPEN(unit=29,file='n4f.'//int2str(iType)//'.dat',status='replace',&
                    action='write',iostat=io_error)
               IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
               DO iGrid = 1, fi%atoms%jri(iType)
                  WRITE(29,'(2e20.8)') fi%atoms%rmsh(iGrid,iType), n_0(iGrid)
               ENDDO
               CLOSE(unit=29,iostat=io_error)
               IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
#endif
            ENDIF

         ENDIF

         IF(fi%atoms%l_outputCFpot(iType)) THEN
            !Run vgen again to obtain the right potential (without external and 4f)
            inDenCF = inDen
            atoms_dummy = fi%atoms

            IF(fi%atoms%l_outputCFcdn(iType).AND.fi%atoms%l_outputCFremove4f(iType)) THEN
               !Remove atomic 4f density before vgen
               DO ispin = 1, fi%input%jspins
                  inDenCF%mt(:,0,iType,ispin) = inDenCF%mt(:,0,iType,ispin) - hub1data%cdn_atomic(:,lcf,iType,ispin)
               ENDDO
               !Remove the same amount of protons from the core to keep everything charge neutral for vgen
               atoms_dummy%zatom(iType) = atoms_dummy%zatom(iType) - n_0Norm*sfp_const*atoms_dummy%neq(iType)
            ENDIF

            nococonv_dummy = nococonv
            CALL vgen(hybdat, fi%field, fi%input, xcpot, atoms_dummy, sphhar, stars, fi%vacuum, fi%sym, &
                      fi%cell,  fi%sliceplot, fmpi, results_dummy, fi%noco, nococonv_dummy,&
                      EnergyDen, inDenCF, vCF, vx, vCoul, vxc, exc)


            IF(fmpi%irank==0) THEN
               !                          sigma
               !Decompose potential into V(r)
               !                          lm
               vlm = cmplx_0
               DO ispin = 1, fi%input%jspins
                  CALL lattHarmsRepToSphHarms(fi%sym, fi%atoms, sphhar, iType, vCF%mt(:,0:,iType,ispin), vlm(:,:,ispin))
               ENDDO

               !Missing: only write out relevant components
#ifdef CPP_HDF
               CALL writeCFpot(cfFileID, fi%atoms, fi%input, iType, vlm)
#else
               !Stupid text output
               DO l = 2, 6, 2
                  DO m = -l, l
                     lm = l*(l+1) + m + 1
                     OPEN(unit=29,file='V_'//int2str(l)//int2str(m)//'.'//int2str(iType)//'.dat',status='replace',&
                          action='write',iostat=io_error)
                     IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
                     DO iGrid = 1, fi%atoms%jri(iType)
                        WRITE(29,'(5e20.8)') fi%atoms%rmsh(iGrid,iType), vlm(iGrid,lm,1), vlm(iGrid,lm,fi%input%jspins)
                     ENDDO
                     CLOSE(unit=29,iostat=io_error)
                     IF(io_error/=0) CALL juDFT_error("IO error", calledby="writeCFOutput")
                  ENDDO
               ENDDO
#endif
            ENDIF
         ENDIF

      ENDDO

#ifdef CPP_HDF
      IF(fmpi%irank==0) CALL closecfFile(cfFileID)
#endif
      CALL timestop("Crystal Field Output")

   END SUBROUTINE writeCFOutput

END MODULE m_writeCFOutput
