MODULE m_greensfTorque

   USE m_types
   USE m_juDFT
   USE m_constants
   USE m_intgr
   USE m_gaunt
   USE m_xmlOutput
   USE m_lattHarmsSphHarmsConv
#ifdef CPP_MPI
   USE mpi
#endif

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE greensfTorque(greensFunction,gfinp,fmpi,sphhar,atoms,sym,noco,nococonv,input,f,g,flo,vTot)

      !--------------------------------------------------------------------------
      ! This Subroutine implements the formula:
      !   alpha     1                  ->    Ef        i ->    alpha     ->->
      !  J      = - -  Im Tr      int dx  int   dE    B (x) sigma    G  (x,x,E)
      !   i         pi      sigma           -infinity  xc             ii
      !
      ! For the evaluation of the torque at site i
      !--------------------------------------------------------------------------

      TYPE(t_greensf),        INTENT(IN)  :: greensFunction(:)
      TYPE(t_gfinp),          INTENT(IN)  :: gfinp
      TYPE(t_mpi),            INTENT(IN)  :: fmpi
      TYPE(t_sphhar),         INTENT(IN)  :: sphhar
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      TYPE(t_sym),            INTENT(IN)  :: sym
      TYPE(t_noco),           INTENT(IN)  :: noco
      TYPE(t_nococonv),       INTENT(IN)  :: nococonv
      TYPE(t_input),          INTENT(IN)  :: input
      REAL,                   INTENT(IN)  :: f(:,:,0:,:,:)
      REAL,                   INTENT(IN)  :: g(:,:,0:,:,:)
      REAL,                   INTENT(IN)  :: flo(:,:,:,:,:)
      TYPE(t_potden),         INTENT(IN)  :: vTot

      INTEGER :: l,lp,iContour,iGrid,ispin,iTorque,atomType,index_task,extra,ierr
      INTEGER :: lh,mu,m,mp,iz,ipm,jr,alpha,lhmu,index,index_start,index_end,n,i_gf
      COMPLEX :: phaseFactor, weight
      REAL    :: realIntegral
      CHARACTER(LEN=20) :: attributes(5)

      TYPE(t_greensf) :: gf_rot

      REAL,    ALLOCATABLE :: torque(:,:),rtmp(:)
      COMPLEX, ALLOCATABLE :: bxc(:,:,:)
      COMPLEX, ALLOCATABLE :: integrand(:,:)
      COMPLEX, ALLOCATABLE :: g_ii(:,:), mag_ii(:,:)
      COMPLEX, ALLOCATABLE :: vlm(:,:,:)
      INTEGER, ALLOCATABLE :: gf_indices(:)

      ALLOCATE(gf_indices(SUM(gfinp%numTorqueElems(:))), source=-1)
      ALLOCATE(bxc(atoms%jmtd,atoms%lmaxd*(atoms%lmaxd+2)+1,atoms%ntype), source=cmplx_0)
      CALL timestart("Green's Function Torque: init")

      DO atomType = 1, atoms%ntype
         IF(gfinp%numTorqueElems(atomType)==0) CYCLE

         iContour = -1
         DO iTorque = 1, gfinp%numTorqueElems(atomType)
            gf_indices(SUM(gfinp%numTorqueElems(:atomType-1))+iTorque) = gfinp%torqueElem(atomType,iTorque)

            !Check that its actually right
            IF(greensFunction(gfinp%torqueElem(atomType,iTorque))%elem%atomType.NE.atomType.OR.&
               greensFunction(gfinp%torqueElem(atomType,iTorque))%elem%atomTypep.NE.atomType) THEN
               CALL juDFT_error("Provided greensFunction for wrong atomType", calledby="greensFunctionTorque")
            ENDIF

            IF(iContour == -1) THEN
               iContour = greensFunction(gfinp%torqueElem(atomType,iTorque))%elem%iContour
            ELSE IF(greensFunction(gfinp%torqueElem(atomType,iTorque))%elem%iContour/=iContour) THEN
               CALL juDFT_error("Provided different energy contours", calledby="greensFunctionTorque")
            ENDIF
         ENDDO
         !Get Bxc from the total potential (local frame)
         !TODO: FFN components
         ALLOCATE(vlm(atoms%jmtd,atoms%lmaxd*(atoms%lmaxd+2)+1,input%jspins),source=cmplx_0)
         vlm = cmplx_0
         DO ispin = 1, input%jspins
            CALL lattHarmsRepToSphHarms(sym, atoms, sphhar, atomType, vTot%mt(:,0:,atomType,ispin), vlm(:,:,ispin))
         ENDDO
         !Get the Bxc part of the potential
         bxc(:,:,atomType) = (vlm(:,:,1) - vlm(:,:,2))/2.0
         DEALLOCATE(vlm)

         !L=0 of potential has an additional rescaling of r/sqrt(4pi)
         bxc(:atoms%jri(atomType),1,atomType) = bxc(:atoms%jri(atomType),1,atomType) &
                                               * sfp_const/atoms%rmsh(:atoms%jri(atomType),atomType)
      ENDDO

      IF(ANY(gf_indices<1) .OR. ANY(gf_indices>SIZE(greensFunction))) THEN
         CALL juDFT_error("Invalid index in greensFunction mapping array", calledby="greensFunctionTorque")
      ENDIF

#ifdef CPP_MPI
      IF(fmpi%isize > 1) THEN
         !Just distribute the individual gf elements over the ranks
         index_task = FLOOR(REAL(SIZE(gf_indices))/(fmpi%isize))
         extra = SIZE(gf_indices) - index_task*fmpi%isize
         index_start = fmpi%irank*index_task + 1 + extra
         index_end = (fmpi%irank+1)*index_task   + extra
         IF(fmpi%irank < extra) THEN
            index_start = index_start - (extra - fmpi%irank)
            index_end = index_end - (extra - fmpi%irank - 1)
         ENDIF
      ELSE
         index_start = 1
         index_end = SIZE(gf_indices)
      ENDIF
#else
      index_start = 1
      index_end = SIZE(gf_indices)
#endif


      CALL timestop("Green's Function Torque: init")
      CALL timestart("Green's Function Torque: Integration")

      ALLOCATE(torque(3,atoms%ntype), source=0.0)
      !$ phaseFactor = gaunt1(0,0,0,0,0,0,atoms%lmaxd)
      DO index = index_start, index_end
         IF(index.LT.1 .OR. index.GT.SIZE(gf_indices)) CYCLE
         i_gf = gf_indices(index)

         gf_rot = greensFunction(i_gf)
         l  = gf_rot%elem%l
         lp = gf_rot%elem%lp
         atomType = gf_rot%elem%atomType

         !Rotate the greens function into the global real space frame
         IF(noco%l_noco) THEN
            CALL gf_rot%rotate_euler_angles(atoms,nococonv%alph(atomType),nococonv%beta(atomType),0.0)
         ELSE IF(noco%l_soc) THEN
            CALL gf_rot%rotate_euler_angles(atoms,nococonv%phi,nococonv%theta,0.0)
         ENDIF

#ifndef CPP_NOTYPEPROCINOMP
         !$OMP parallel default(none) &
         !$OMP shared(sphhar,atoms,input,gf_rot,f,g,flo,bxc) &
         !$OMP shared(l,lp,atomType,torque) &
         !$OMP private(lh,m,mu,mp,lhmu,phaseFactor,weight,ispin,ipm,iz,alpha,jr) &
         !$OMP private(realIntegral,integrand,g_ii,mag_ii)
#endif
         ALLOCATE(integrand(atoms%jmtd,3),source=cmplx_0)
         ALLOCATE(g_ii(atoms%jmtd,gf_rot%contour%nz),source=cmplx_0)
         ALLOCATE(mag_ii(atoms%jmtd,gf_rot%contour%nz),source=cmplx_0)
#ifndef CPP_NOTYPEPROCINOMP
         !$OMP do collapse(2)
#endif
         DO lh = 0, atoms%lmaxd
            DO m = -l, l
               IF(MOD(lh+l+lp,2) .NE. 0) CYCLE
               IF(lh.GT.l+lp) CYCLE
               IF(lh.LT.abs(l-lp)) CYCLE
               DO mu = -lh, lh
                  lhmu = lh * (lh+1) + mu + 1
                  mp = m + mu
                  IF(ABS(mp).GT.lp) CYCLE
                  phaseFactor = gaunt1(lp,lh,l,mp,mu,m,atoms%lmaxd)
                  IF(ABS(phaseFactor).LT.1e-12) CYCLE
                  integrand = cmplx_0
                  DO ipm = 1, 2
                     DO alpha = 1, 3 !(x,y,z)
                        IF (alpha.EQ.1) THEN
                           !magnetization in x-direction
                           CALL gf_rot%getRadial(atoms,m,mp,ipm==2,3,f,g,flo,mag_ii)
                           CALL gf_rot%getRadial(atoms,mp,m,ipm==2,3,f,g,flo,g_ii)
                           mag_ii = mag_ii + conjg(g_ii)
                        ELSE IF (alpha.EQ.2) THEN
                           !magnetization in y-direction
                           CALL gf_rot%getRadial(atoms,m,mp,ipm==2,3,f,g,flo,mag_ii)
                           CALL gf_rot%getRadial(atoms,mp,m,ipm==2,3,f,g,flo,g_ii)
                           mag_ii = ImagUnit * (mag_ii - conjg(g_ii))
                        ELSE
                           !magnetization in z-direction
                           mag_ii = cmplx_0
                           DO ispin = 1, input%jspins
                              CALL gf_rot%getRadial(atoms,m,mp,ipm==2,ispin,f,g,flo,g_ii)
                              mag_ii = mag_ii + (-1)**(ispin-1) * g_ii
                           ENDDO
                        ENDIF

                        DO iz = 1, SIZE(mag_ii,2)
                           weight = gf_rot%contour%de(iz) * phaseFactor
                           weight = MERGE(weight, conjg(weight), ipm==1)

                           DO jr = 1, atoms%jri(atomType)
                              integrand(jr,alpha) = integrand(jr,alpha) + ImagUnit/tpi_const * (-1)**(ipm-1) * mag_ii(jr,iz) &
                                                                         * MERGE(bxc(jr,lhmu,atomType),conjg(bxc(jr,lhmu,atomType)), ipm==1) * weight
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO

                  DO alpha = 1, 3 !(x,y,z)
                     CALL intgr3(REAL(integrand(:,alpha)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realIntegral)
#ifndef CPP_NOTYPEPROCINOMP
                     !$OMP critical
                     torque(alpha,atomType) = torque(alpha,atomType) + realIntegral
                     !$OMP end critical
#else
                     torque(alpha,atomType) = torque(alpha,atomType) + realIntegral
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
#ifndef CPP_NOTYPEPROCINOMP
         !$OMP end do
         DEALLOCATE(integrand,g_ii,mag_ii)
         !$OMP end parallel
#else
         DEALLOCATE(integrand,g_ii,mag_ii)
#endif

      ENDDO
      CALL timestop("Green's Function Torque: Integration")

#ifdef CPP_MPI
      !Collect the torque to rank 0
      n = SIZE(torque)
      ALLOCATE(rtmp(n))
      CALL MPI_REDUCE(torque,rtmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF(fmpi%irank.EQ.0) CALL dcopy(n,rtmp,1,torque,1)
      DEALLOCATE(rtmp)
#endif

      IF(fmpi%irank.EQ.0) THEN
         CALL openXMLElementNoAttributes('noncollinearTorque')
         WRITE(oUnit,'(/,A)') 'Torque Calculation (noco):'
         WRITE(oUnit,'(/,A)') '---------------------------'
         DO atomType = 1, atoms%ntype
            IF(gfinp%numTorqueElems(atomType)==0) CYCLE
            WRITE(oUnit,'(A,I4,A,3f16.8,A)') '  atom: ', atomType, '   torque: ', torque(:,atomType) * hartree_to_ev_const * 1000, ' meV'

            attributes = ''
            WRITE(attributes(1),'(i0)') atomType
            WRITE(attributes(2),'(f14.8)') torque(1,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(3),'(f14.8)') torque(2,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(4),'(f14.8)') torque(3,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(5),'(a3)') 'meV'
            CALL writeXMLElementForm('torque',['atomType','sigma_x ','sigma_y ','sigma_z ','units   '],&
                                     attributes,reshape([8,7,7,7,4,6,14,14,14,3],[5,2]))
         ENDDO
         CALL closeXMLElement('noncollinearTorque')
      ENDIF

   END SUBROUTINE greensfTorque

END MODULE m_greensfTorque
