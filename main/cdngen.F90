!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_cdngen

USE m_juDFT

CONTAINS

SUBROUTINE cdngen(eig_id,mpi,input,banddos,sliceplot,vacuum,&
                  dimension,kpts,atoms,sphhar,stars,sym,obsolete,&
                  enpara,cell,noco,jij,vTot,results,oneD,coreSpecInput,&
                  inDen,outDen)

   !*****************************************************
   !    Charge density generator
   !    calls cdnval to generate the valence charge and the
   !    core routines for the core contribution
   !*****************************************************

   USE m_constants
   USE m_prpqfftmap
   USE m_cdnval
   USE m_cdn_io
   USE m_wrtdop
   USE m_cdntot
   USE m_cdnovlp
   USE m_qfix
   USE m_genNewNocoInp
   USE m_types
   USE m_xmlOutput
   USE m_magMoms
   USE m_orbMagMoms
   USE m_cdncore
#ifdef CPP_MPI
   USE m_mpi_bc_potden
#endif

   IMPLICIT NONE

   ! Type instance arguments
   TYPE(t_results),INTENT(INOUT)    :: results
   TYPE(t_mpi),INTENT(IN)           :: mpi
   TYPE(t_dimension),INTENT(IN)     :: dimension
   TYPE(t_oneD),INTENT(IN)          :: oneD
   TYPE(t_enpara),INTENT(INOUT)     :: enpara
   TYPE(t_obsolete),INTENT(IN)      :: obsolete
   TYPE(t_banddos),INTENT(IN)       :: banddos
   TYPE(t_sliceplot),INTENT(IN)     :: sliceplot
   TYPE(t_input),INTENT(IN)         :: input
   TYPE(t_vacuum),INTENT(IN)        :: vacuum
   TYPE(t_noco),INTENT(IN)          :: noco
   TYPE(t_jij),INTENT(IN)           :: jij
   TYPE(t_sym),INTENT(IN)           :: sym
   TYPE(t_stars),INTENT(IN)         :: stars
   TYPE(t_cell),INTENT(IN)          :: cell
   TYPE(t_kpts),INTENT(IN)          :: kpts
   TYPE(t_sphhar),INTENT(IN)        :: sphhar
   TYPE(t_atoms),INTENT(IN)         :: atoms
   TYPE(t_coreSpecInput),INTENT(IN) :: coreSpecInput
   TYPE(t_potden),INTENT(IN)        :: vTot
   TYPE(t_potden),INTENT(INOUT)     :: inDen,outDen

   !Scalar Arguments
   INTEGER, INTENT (IN)             :: eig_id

   ! Local type instances
   TYPE(t_noco) :: noco_new

   !Local Scalars
   REAL fix,qtot,scor,stot,sval,dummy
   INTEGER ivac,j,jspin,jspmax,k,iType
   LOGICAL l_enpara

   !Local Arrays
   REAL stdn(atoms%ntype,dimension%jspd),svdn(atoms%ntype,dimension%jspd),alpha_l(atoms%ntype)
   REAL chmom(atoms%ntype,dimension%jspd),clmom(3,atoms%ntype,dimension%jspd)
   INTEGER,ALLOCATABLE :: igq_fft(:)
   REAL   ,ALLOCATABLE :: qvac(:,:,:,:),qvlay(:,:,:,:,:)

   !pk non-collinear (start)
   INTEGER igq2_fft(0:stars%kq1_fft*stars%kq2_fft-1)
   COMPLEX,ALLOCATABLE :: qa21(:), cdomvz(:,:)
   !pk non-collinear (end)

   CALL outDen%init(stars,atoms,sphhar,vacuum,noco,oneD,input%jspins,noco%l_noco,POTDEN_TYPE_DEN)

   IF (mpi%irank.EQ.0) THEN
      INQUIRE(file='enpara',exist=l_enpara)
      IF (l_enpara) OPEN (40,file ='enpara',form = 'formatted',status ='unknown')
   ENDIF
   ALLOCATE (qa21(atoms%ntype))
   ALLOCATE (qvac(dimension%neigd,2,kpts%nkpt,dimension%jspd))
   ALLOCATE (qvlay(dimension%neigd,vacuum%layerd,2,kpts%nkpt,dimension%jspd))
   ALLOCATE (igq_fft(0:stars%kq1_fft*stars%kq2_fft*stars%kq3_fft-1))

   !initialize density arrays with zero
   qa21(:) = cmplx(0.0,0.0)
   qvac(:,:,:,:) = 0.0 
   qvlay(:,:,:,:,:) = 0.0

   outDen%iter = inDen%iter
        
   !Set up pointer for backtransformation of from g-vector in
   !positive domain fof carge density fftibox into stars
   !In principle this can also be done in main program once.
   !It is done here to save memory.
   CALL prp_qfft_map(stars,sym, input, igq2_fft,igq_fft)

   !in a non-collinear calcuation where the off-diagonal part of
   !density matrix in the muffin-tins is calculated, the a- and
   !b-coef. for both spins are needed at once. Thus, cdnval is only
   !called once and both spin directions are calculated in a single
   !go.
   IF (mpi%irank.EQ.0) CALL openXMLElementNoAttributes('valenceDensity')

   jspmax = input%jspins
   IF (noco%l_mperp) jspmax = 1
   DO jspin = 1,jspmax
      CALL timestart("cdngen: cdnval")
      CALL cdnval(eig_id,&
                  mpi,kpts,jspin,sliceplot,noco, input,banddos,cell,atoms,enpara,stars, vacuum,dimension,&
                  sphhar,sym,obsolete,igq_fft,vTot,oneD,coreSpecInput,&
                  outDen,results,qvac,qvlay,qa21, chmom,clmom)
      CALL timestop("cdngen: cdnval")
   END DO

   IF (mpi%irank.EQ.0) THEN
      IF (l_enpara) CLOSE (40)
      CALL cdntot(stars,atoms,sym, vacuum,input,cell,oneD, outDen%pw,outDen%mt,outDen%vacz,.TRUE., qtot,dummy)
      CALL closeXMLElement('valenceDensity')
   END IF ! mpi%irank = 0

   CALL cdncore(results,mpi,dimension,oneD,sliceplot,input,vacuum,noco,sym,&
                stars,cell,sphhar,atoms,vTot,outDen,stdn,svdn)

   IF (mpi%irank.EQ.0) THEN
      !block 2 unnecessary for slicing: begin
      IF (.NOT.sliceplot%slice) THEN
         CALL openXMLElementNoAttributes('allElectronCharges')
         CALL qfix(stars,atoms,sym,vacuum, sphhar,input,cell,oneD,outDen,.TRUE.,.true.,fix)
         CALL closeXMLElement('allElectronCharges')
         !pk non-collinear (start)
         IF (noco%l_noco) THEN
            !fix also the off-diagonal part of the density matrix
            outDen%pw(:stars%ng3,3) = fix*outDen%pw(:stars%ng3,3)
            IF (input%film) THEN
               outDen%vacz(:,:,3:4) = fix*outDen%vacz(:,:,3:4)
               outDen%vacxy(:,:,:,3) = fix*outDen%vacxy(:,:,:,3)
            END IF
         END IF
         !pk non-collinear (end)

         !spin densities at the nucleus
         !and magnetic moment in the spheres
         IF (input%jspins.EQ.2) THEN
            WRITE (6,FMT=8000)
            WRITE (16,FMT=8000)
            DO iType = 1,atoms%ntype
               sval = svdn(iType,1) - svdn(iType,input%jspins)
               stot = stdn(iType,1) - stdn(iType,input%jspins)
               scor = stot - sval
               WRITE (6,FMT=8010) iType,stot,sval,scor,svdn(iType,1),stdn(iType,1)
               WRITE (16,FMT=8010) iType,stot,sval,scor,svdn(iType,1),stdn(iType,1)
            END DO

            noco_new = noco

            CALL magMoms(dimension,input,atoms,noco_new,vTot,chmom,qa21)

            !Generate and save the new nocoinp file if the directions of the local
            !moments are relaxed or a constraint B-field is calculated.
            IF (ANY(noco%l_relax(:atoms%ntype)).OR.noco%l_constr) THEN
               CALL genNewNocoInp(input,atoms,jij,noco,noco_new)
            END IF

            IF (noco%l_soc) CALL orbMagMoms(dimension,atoms,noco,clmom)
         END IF
      !block 2 unnecessary for slicing: end
      END IF ! .NOT.sliceplot%slice

      8000 FORMAT (/,/,10x,'spin density at the nucleus:',/,10x,'type',t25,&
                   'input%total',t42,'valence',t65,'core',t90,&
                   'majority valence and input%total density',/)
      8010 FORMAT (i13,2x,3e20.8,5x,2e20.8)

      IF (sliceplot%slice) THEN
         OPEN (20,file='cdn_slice',form='unformatted',status='unknown')
         CALL wrtdop(stars,vacuum,atoms,sphhar, input,sym, 20, outDen%iter,outDen%mt,outDen%pw,outDen%vacz,outDen%vacxy)
         IF (noco%l_noco) THEN
            WRITE (20) (outDen%pw(k,3),k=1,stars%ng3)
            IF (input%film) THEN
               ALLOCATE(cdomvz(vacuum%nmz,vacuum%nvac))
               DO ivac = 1, vacuum%nvac
                  DO j = 1, vacuum%nmz
                     cdomvz(j,ivac) = CMPLX(outDen%vacz(j,ivac,3),outDen%vacz(j,ivac,4))
                  END DO
               END DO
               WRITE (20) ((cdomvz(j,ivac),j=1,vacuum%nmz),ivac=1,vacuum%nvac)
               WRITE (20) (((outDen%vacxy(j,k-1,ivac,3),j=1,vacuum%nmzxy),k=2,oneD%odi%nq2) ,ivac=1,vacuum%nvac)
               DEALLOCATE(cdomvz)
            END IF
         END IF
         CLOSE(20) 
         CALL juDFT_end("slice OK")
      END IF

      IF(vacuum%nvac.EQ.1) THEN
         outDen%vacz(:,2,:) = outDen%vacz(:,1,:)
         IF (sym%invs) THEN
            outDen%vacxy(:,:,2,:) = CONJG(outDen%vacxy(:,:,1,:))
         ELSE
            outDen%vacxy(:,:,2,:) = outDen%vacxy(:,:,1,:)
         END IF
      END IF

   ENDIF ! mpi%irank.EQ.0

#ifdef CPP_MPI
   CALL mpi_bc_potden(mpi,stars,sphhar,atoms,input,vacuum,oneD,noco,outDen)
#endif

   DEALLOCATE (qvac,qvlay,qa21)
   DEALLOCATE (igq_fft)

   IF (sliceplot%slice) CALL juDFT_end("sliceplot%slice OK",mpi%irank)

   RETURN

END SUBROUTINE cdngen

END MODULE m_cdngen
