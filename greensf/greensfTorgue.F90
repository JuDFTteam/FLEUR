MODULE m_greensfTorgue

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

#include"cpp_double.h"

   CONTAINS

   SUBROUTINE greensfTorgue(greensFunction,gfinp,fmpi,sphhar,atoms,sym,noco,nococonv,input,f,g,flo,vTot)

      !--------------------------------------------------------------------------
      ! This Subroutine implements the formula:
      !   alpha     1                  ->    Ef        i ->    alpha     ->->
      !  J      = - -  Im Tr      int dx  int   dE    B (x) sigma    G  (x,x,E)
      !   i         pi      sigma           -infinity  xc             ii
      !
      ! For the evaluation of the torgue at site i
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

      INTEGER :: l,lp,iContour,iGrid,ispin,iTorgue,atomType,index_task,extra,ierr
      INTEGER :: lh,mu,m,mp,iz,ipm,jr,alpha,lhmu,index,index_start,index_end,n,i_gf
      COMPLEX :: phaseFactor, weight
      REAL    :: realIntegral
      COMPLEX :: sigma(2,2,3),g_Spin(2,2)
      CHARACTER(LEN=20) :: attributes(5)

      REAL,    ALLOCATABLE :: torgue(:,:),rtmp(:)
      COMPLEX, ALLOCATABLE :: bxc(:,:,:)
      COMPLEX, ALLOCATABLE :: integrand(:,:)
      COMPLEX, ALLOCATABLE :: g_ii(:,:,:,:)
      COMPLEX, ALLOCATABLE :: vlm(:,:,:)
      INTEGER, ALLOCATABLE :: gf_indices(:)

      ALLOCATE(gf_indices(SUM(gfinp%numTorgueElems(:))), source=-1)
      ALLOCATE(bxc(atoms%jmtd,atoms%lmaxd*(atoms%lmaxd+2)+1,atoms%ntype), source=cmplx_0)
      CALL timestart("Green's Function Torgue: init")

      DO atomType = 1, atoms%ntype
         IF(gfinp%numTorgueElems(atomType)==0) CYCLE

         iContour = -1
         DO iTorgue = 1, gfinp%numTorgueElems(atomType)
            gf_indices(SUM(gfinp%numTorgueElems(:atomType-1))+iTorgue) = gfinp%torgueElem(atomType,iTorgue)

            !Check that its actually right
            IF(greensFunction(gfinp%torgueElem(atomType,iTorgue))%elem%atomType.NE.atomType.OR.&
               greensFunction(gfinp%torgueElem(atomType,iTorgue))%elem%atomTypep.NE.atomType) THEN
               CALL juDFT_error("Provided greensFunction for wrong atomType", calledby="greensFunctionTorgue")
            ENDIF

            IF(iContour == -1) THEN
               iContour = greensFunction(gfinp%torgueElem(atomType,iTorgue))%elem%iContour
            ELSE IF(greensFunction(gfinp%torgueElem(atomType,iTorgue))%elem%iContour/=iContour) THEN
               CALL juDFT_error("Provided different energy contours", calledby="greensFunctionTorgue")
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
         CALL juDFT_error("Invalid index in greensFunction mapping array", calledby="greensFunctionTorgue")
      ENDIF

      ! sigma are the Pauli matrices
      sigma=cmplx_0
      sigma(1,2,1)=CMPLX(1.0,0.0)
      sigma(2,1,1)=CMPLX(1.0,0.0)
      sigma(1,2,2)=CMPLX(0.0,-1.0)
      sigma(2,1,2)=CMPLX(0.0,1.0)
      sigma(1,1,3)=CMPLX(1.0,0.0)
      sigma(2,2,3)=CMPLX(-1.0,0.0)

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


      CALL timestop("Green's Function Torgue: init")
      CALL timestart("Green's Function Torgue: Integration")

      ALLOCATE(torgue(3,atoms%ntype), source=0.0)
      DO index = index_start, index_end
         IF(index.LT.1 .OR. index.GT.SIZE(gf_indices)) CYCLE
         i_gf = gf_indices(index)

         l  = greensFunction(i_gf)%elem%l
         lp = greensFunction(i_gf)%elem%lp
         atomType = greensFunction(i_gf)%elem%atomType

#ifndef CPP_NOTYPEPROCINOMP
         !$OMP parallel default(none) &
         !$OMP shared(sphhar,atoms,greensFunction,i_gf,f,g,flo,sigma,bxc) &
         !$OMP shared(l,lp,atomType,torgue) &
         !$OMP private(lh,m,mu,mp,lhmu,phaseFactor,weight,ipm,iz,alpha,jr) &
         !$OMP private(realIntegral,integrand,g_ii,g_Spin)
#endif
         ALLOCATE(integrand(atoms%jmtd,3),source=cmplx_0)
         ALLOCATE(g_ii(2,2,atoms%jmtd,greensFunction(i_gf)%contour%nz),source=cmplx_0)
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
                     CALL greensFunction(i_gf)%getRadialSpin(atoms,m,mp,ipm==2,f,g,flo,g_ii)
                     DO iz = 1, SIZE(g_ii,4)
                        weight = greensFunction(i_gf)%contour%de(iz) * phaseFactor

                        IF(ipm == 1) THEN
                           DO alpha = 1, 3 !(x,y,z)
                              DO jr = 1, atoms%jri(atomType)
                                 g_Spin = matmul(sigma(:,:,alpha),g_ii(:,:,jr,iz))
                                 integrand(jr,alpha) = integrand(jr,alpha) + ImagUnit/tpi_const * (g_Spin(1,1) + g_Spin(2,2)) &
                                                                            * bxc(jr,lhmu,atomType) * weight
                              ENDDO
                           ENDDO
                        ELSE
                           DO alpha = 1, 3 !(x,y,z)
                              DO jr = 1, atoms%jri(atomType)
                                 g_Spin = matmul(conjg(sigma(:,:,alpha)),g_ii(:,:,jr,iz))
                                 integrand(jr,alpha) = integrand(jr,alpha) - ImagUnit/tpi_const * (g_Spin(1,1) + g_Spin(2,2)) &
                                                                            * conjg(bxc(jr,lhmu,atomType) * weight)
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO

                  DO alpha = 1, 3 !(x,y,z)
                     CALL intgr3(REAL(integrand(:,alpha)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realIntegral)
#ifndef CPP_NOTYPEPROCINOMP
                     !$OMP critical
                     torgue(alpha,atomType) = torgue(alpha,atomType) + realIntegral
                     !$OMP end critical
#else
                     torgue(alpha,atomType) = torgue(alpha,atomType) + realIntegral
#endif
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
#ifndef CPP_NOTYPEPROCINOMP
         !$OMP end do
         DEALLOCATE(integrand,g_ii)
         !$OMP end parallel
#else
         DEALLOCATE(integrand,g_ii)
#endif

      ENDDO
      CALL timestop("Green's Function Torgue: Integration")

#ifdef CPP_MPI
      !Collect the torgue to rank 0
      n = SIZE(torgue)
      ALLOCATE(rtmp(n))
      CALL MPI_REDUCE(torgue,rtmp,n,CPP_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF(fmpi%irank.EQ.0) CALL CPP_BLAS_scopy(n,rtmp,1,torgue,1)
      DEALLOCATE(rtmp)
#endif

      IF(fmpi%irank.EQ.0) THEN
         CALL openXMLElementNoAttributes('noncollinearTorgue')
         WRITE(oUnit,'(/,A)') 'Torgue Calculation (noco):'
         WRITE(oUnit,'(/,A)') '---------------------------'
         DO atomType = 1, atoms%ntype
            IF(gfinp%numTorgueElems(atomType)==0) CYCLE
            WRITE(oUnit,'(A,I4,A,3f14.8,A)') '  atom: ', atomType, '   torgue: ', torgue(:,atomType) * hartree_to_ev_const * 1000, ' meV'

            attributes = ''
            WRITE(attributes(1),'(i0)') atomType
            WRITE(attributes(2),'(f14.8)') torgue(1,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(3),'(f14.8)') torgue(2,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(4),'(f14.8)') torgue(3,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(5),'(a3)') 'meV'
            CALL writeXMLElementForm('torgue',['atomType','sigma_x ','sigma_y ','sigma_z ','unit    '],&
                                     attributes,reshape([8,7,7,7,4,6,14,14,14,3],[5,2]))
         ENDDO
         CALL closeXMLElement('noncollinearTorgue')
      ENDIF

   END SUBROUTINE greensfTorgue

   SUBROUTINE greensfSOTorgue(greensFunction,gfinp,fmpi,sphhar,atoms,sym,noco,nococonv,input,enpara,f,g,flo,vTot)

      USE m_sointg
      USE m_spnorb
      USE m_fourProduct

      !--------------------------------------------------------------------------
      ! This Subroutine implements the formula:
      !   alpha     1                  ->    Ef        i ->    alpha     ->->
      !  J      = - -  Im Tr      int dx  int   dE    V (x) sigma    G  (x,x,E)
      !   i         pi      sigma           -infinity  SO             ii
      !
      ! For the evaluation of the torgue at site i (argument atomType) caused by Spin-orbit-coupling
      !--------------------------------------------------------------------------

      TYPE(t_greensf),        INTENT(IN)  :: greensFunction(:)
      TYPE(t_gfinp),          INTENT(IN)  :: gfinp
      TYPE(t_mpi),            INTENT(IN)  :: fmpi
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      TYPE(t_sphhar),         INTENT(IN)  :: sphhar
      TYPE(t_sym),            INTENT(IN)  :: sym
      TYPE(t_noco),           INTENT(IN)  :: noco
      TYPE(t_nococonv),       INTENT(IN)  :: nococonv
      TYPE(t_input),          INTENT(IN)  :: input
      TYPE(t_enpara),         INTENT(IN)  :: enpara
      REAL,                   INTENT(IN)  :: f(:,:,0:,:,:)
      REAL,                   INTENT(IN)  :: g(:,:,0:,:,:)
      REAL,                   INTENT(IN)  :: flo(:,:,:,:,:)
      TYPE(t_potden),         INTENT(IN)  :: vTot

      INTEGER :: jspin,na,nsym,nh,i_gf,l,lp,spin,iContour,i,n,iTorgue
      INTEGER :: index_start,index_end, index_task, extra, index,ierr
      INTEGER :: lh,mh,mhp,m,mp,iz,ipm,jr,alpha,jspin1,jspin2,atomType
      COMPLEX :: phaseFactor,weight
      REAL    :: realIntegral, imagIntegral, e
      COMPLEX :: sigma(2,2,3),chi(2,2),g_Spin(2,2)
      CHARACTER(LEN=20) :: attributes(5)

      REAL :: v0(atoms%jmtd),vso_tmp(atoms%jmtd,2)
      COMPLEX, ALLOCATABLE :: integrand(:,:),g_ii(:,:,:,:)
      COMPLEX,ALLOCATABLE :: soangl(:,:,:,:,:,:)
      COMPLEX,ALLOCATABLE :: vso(:,:,:,:,:,:,:)
      INTEGER,ALLOCATABLE :: gf_indices(:)
      REAL, ALLOCATABLE   :: torgue(:,:)
      REAL, ALLOCATABLE   :: rtmp(:)


      ALLOCATE(gf_indices(SUM(gfinp%numTorgueElems(:))), source=-1)
      ALLOCATE(vso(atoms%jmtd,2,2,-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,0:atoms%lmaxd,atoms%ntype), source=cmplx_0)
      CALL timestart("Green's Function SOTorgue: init")

      DO atomType = 1, atoms%ntype
         IF(gfinp%numTorgueElems(atomType)==0) CYCLE

         iContour = -1
         DO iTorgue = 1, gfinp%numTorgueElems(atomType)
            gf_indices(SUM(gfinp%numTorgueElems(:atomType-1))+iTorgue) = gfinp%torgueElem(atomType,iTorgue)

            !Check that its actually right
            IF(greensFunction(gfinp%torgueElem(atomType,iTorgue))%elem%atomType.NE.atomType.OR.&
               greensFunction(gfinp%torgueElem(atomType,iTorgue))%elem%atomTypep.NE.atomType) THEN
               CALL juDFT_error("Provided greensFunction for wrong atomType", calledby="greensFunctionTorgue")
            ENDIF

            IF(iContour == -1) THEN
               iContour = greensFunction(gfinp%torgueElem(atomType,iTorgue))%elem%iContour
            ELSE IF(greensFunction(gfinp%torgueElem(atomType,iTorgue))%elem%iContour/=iContour) THEN
               CALL juDFT_error("Provided different energy contours", calledby="greensFunctionTorgue")
            ENDIF
         ENDDO
         !
         !---> common spin-orbit integrant V   (average spin directions)
         !                                  SO
         ALLOCATE(soangl(atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2,atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,2))
         IF(.NOT.noco%l_noco) THEN
            CALL spnorb_angles(atoms,fmpi,nococonv%beta(atomType),nococonv%alph(atomType),soangl)
         ELSE
            CALL spnorb_angles(atoms,fmpi,nococonv%theta,nococonv%phi,soangl)
         ENDIF

         DO l = 0, atoms%lmaxd
            v0(:) = 0.0
            DO i = 1,atoms%jri(atomType)
               v0(i) = (vtot%mt(i,0,atomType,1)+vtot%mt(i,0,atomType,input%jspins))/2.
            END DO
            e = (enpara%el0(l,atomType,1)+enpara%el0(l,atomType,input%jspins))/2.

            CALL sointg(atomType,e,vtot%mt(:,0,atomType,:),v0,atoms,input,vso_tmp)
            IF (.TRUE.) THEN
               DO i= 1,atoms%jmtd
                  vso_tmp(i,1)= (vso_tmp(i,1)+vso_tmp(i,2))/2.
                  vso_tmp(i,2)= vso_tmp(i,1)
               ENDDO
            ENDIF
            DO m = -l, l
               DO mp = -l,l
                  DO jspin1 = 1,2
                     DO jspin2 = 1,2
                        vso(:,jspin1,jspin2,m,mp,l,atomType) = vso_tmp(:,1) * soangl(l,m,jspin1,l,mp,jspin2)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(soangl)
      ENDDO

      IF(ANY(gf_indices<1) .OR. ANY(gf_indices>SIZE(greensFunction))) THEN
         CALL juDFT_error("Invalid index in greensFunction mapping array", calledby="greensFunctionTorgue")
      ENDIF

      ! sigma are the Pauli matrices
      sigma=cmplx_0
      sigma(1,2,1)=CMPLX(1.0,0.0)
      sigma(2,1,1)=CMPLX(1.0,0.0)
      sigma(1,2,2)=CMPLX(0.0,-1.0)
      sigma(2,1,2)=CMPLX(0.0,1.0)
      sigma(1,1,3)=CMPLX(1.0,0.0)
      sigma(2,2,3)=CMPLX(-1.0,0.0)

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

      CALL timestop("Green's Function SOTorgue: init")
      CALL timestart("Green's Function SOTorgue: Integration")

      ALLOCATE(torgue(3,atoms%ntype), source=0.0)
      DO index = index_start, index_end
         IF(index.LT.1 .OR. index.GT.SIZE(gf_indices)) CYCLE
         i_gf = gf_indices(index)

         l  = greensFunction(i_gf)%elem%l
         lp = greensFunction(i_gf)%elem%lp
         atomType = greensFunction(i_gf)%elem%atomType

#ifndef CPP_NOTYPEPROCINOMP
         !$OMP parallel default(none) &
         !$OMP shared(sphhar,atoms,greensFunction,i_gf,f,g,flo,sigma,vso) &
         !$OMP shared(l,lp,atomType,torgue) &
         !$OMP private(lh,m,mh,mhp,mp,phaseFactor,weight,ipm,iz,alpha,jr) &
         !$OMP private(realIntegral,integrand,g_ii,g_Spin)
#endif
         ALLOCATE(integrand(atoms%jmtd,3),source=cmplx_0)
         ALLOCATE(g_ii(2,2,atoms%jmtd,greensFunction(i_gf)%contour%nz),source=cmplx_0)
#ifndef CPP_NOTYPEPROCINOMP
         !$OMP do collapse(3)
#endif
         DO lh = 0, SIZE(vso,6)-1
            DO m = -l, l
               DO mp = -lp, lp
                  DO mh = -lh,lh
                     DO mhp = -lh, lh
                        IF(MAXVAL(ABS(vso(:,:,:,mh,mhp,lh,atomType))).LT.1e-12) CYCLE
                        phaseFactor = fourProduct(lp,lh,l,lh,mp,mh,m,mhp,atoms%lmaxd)
                        IF(ABS(phaseFactor).LT.1e-12) CYCLE
                        integrand = cmplx_0
                        DO ipm = 1, 2
                           CALL greensFunction(i_gf)%getRadialSpin(atoms,m,mp,ipm==2,f,g,flo,g_ii)
                           DO iz = 1, SIZE(g_ii,4)
                              weight = greensFunction(i_gf)%contour%de(iz) * phaseFactor

                              IF(ipm == 1) THEN
                                 DO alpha = 1, 3 !(x,y,z)
                                    DO jr = 1, atoms%jri(atomType)
                                       g_Spin = matmul(sigma(:,:,alpha),g_ii(:,:,jr,iz))
                                       g_Spin = matmul(vso(jr,:,:,mh,mhp,lh,atomType),g_Spin)
                                       integrand(jr,alpha) = integrand(jr,alpha) + ImagUnit/tpi_const * (g_Spin(1,1) + g_Spin(2,2)) &
                                                                                  * weight
                                    ENDDO
                                 ENDDO
                              ELSE
                                 DO alpha = 1, 3 !(x,y,z)
                                    DO jr = 1, atoms%jri(atomType)
                                       g_Spin = matmul(conjg(sigma(:,:,alpha)),g_ii(:,:,jr,iz))
                                       g_Spin = matmul(conjg(vso(jr,:,:,mh,mhp,lh,atomType)),g_Spin)
                                       integrand(jr,alpha) = integrand(jr,alpha) - ImagUnit/tpi_const * (g_Spin(1,1) + g_Spin(2,2)) &
                                                                                  * conjg(weight)
                                    ENDDO
                                 ENDDO
                              ENDIF
                           ENDDO
                        ENDDO

                        DO alpha = 1, 3 !(x,y,z)
                           CALL intgr3(REAL(integrand(:,alpha)),atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),realIntegral)
#ifndef CPP_NOTYPEPROCINOMP
                           !$OMP critical
                           torgue(alpha,atomType) = torgue(alpha,atomType) + realIntegral
                           !$OMP end critical
#else
                           torgue(alpha,atomType) = torgue(alpha,atomType) + realIntegral
#endif
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
#ifndef CPP_NOTYPEPROCINOMP
         !$OMP end do
         DEALLOCATE(integrand,g_ii)
         !$OMP end parallel
#else
         DEALLOCATE(integrand,g_ii)
#endif

      ENDDO
      CALL timestop("Green's Function SOTorgue: Integration")

#ifdef CPP_MPI
      !Collect the torgue to rank 0
      n = SIZE(torgue)
      ALLOCATE(rtmp(n))
      CALL MPI_REDUCE(torgue,rtmp,n,CPP_MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      IF(fmpi%irank.EQ.0) CALL CPP_BLAS_scopy(n,rtmp,1,torgue,1)
      DEALLOCATE(rtmp)
#endif

      IF(fmpi%irank.EQ.0) THEN
         CALL openXMLElementNoAttributes('spinorbitTorgue')
         WRITE(oUnit,'(/,A)') 'Torgue Calculation (spin-orbit):'
         WRITE(oUnit,'(/,A)') '---------------------------'
         DO atomType = 1, atoms%ntype
            IF(gfinp%numTorgueElems(atomType)==0) CYCLE
            WRITE(oUnit,'(A,I4,A,3f14.8,A)') '  atom: ', atomType, '   torgue: ', torgue(:,atomType) * hartree_to_ev_const * 1000, ' meV'

            attributes = ''
            WRITE(attributes(1),'(i0)') atomType
            WRITE(attributes(2),'(f14.8)') torgue(1,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(3),'(f14.8)') torgue(2,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(4),'(f14.8)') torgue(3,atomType) * hartree_to_ev_const * 1000
            WRITE(attributes(5),'(a3)') 'meV'
            CALL writeXMLElementForm('torgue',['atomType','sigma_x ','sigma_y ','sigma_z ','unit    '],&
                                     attributes,reshape([8,7,7,7,4,6,14,14,14,3],[5,2]))
         ENDDO
         CALL closeXMLElement('spinorbitTorgue')
      ENDIF


   END SUBROUTINE greensfSOTorgue

END MODULE m_greensfTorgue
