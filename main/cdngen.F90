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
                  inIter,inDen,outDen)

   !*****************************************************
   !    Charge density generator
   !    calls cdnval to generate the valence charge and the
   !    core routines for the core contribution
   !*****************************************************

   USE m_constants
   USE m_umix
   USE m_prpqfftmap
   USE m_cdnval
   USE m_cdn_io
   USE m_wrtdop
   USE m_cdntot
   USE m_cdnovlp
   USE m_qfix
   USE m_rwnoco
   USE m_cored
   USE m_coredr
   USE m_m_perp
   USE m_types
   USE m_xmlOutput
#ifdef CPP_MPI
   USE m_mpi_bc_potden
   USE m_mpi_bc_coreden
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
   INTEGER, INTENT (IN)             :: inIter

   ! Local type instances
   TYPE(t_noco) :: noco_new

   !Local Scalars
   REAL fix,qtot,scor,seig,smom,stot,sval,dummy
   REAL slmom,slxmom,slymom,sum,thetai,phii,fermiEnergyTemp
   INTEGER iter,ivac,j,jspin,jspmax,k,n,nt,ieig,ikpt
   INTEGER  ityp,ilayer,urec,itype,iatom,archiveType
   LOGICAL l_relax_any,exst,n_exist,l_qfix, l_enpara
   LOGICAL, PARAMETER :: l_st=.FALSE.

   !Local Arrays
   REAL stdn(atoms%ntype,dimension%jspd),svdn(atoms%ntype,dimension%jspd),alpha_l(atoms%ntype)
   REAL rh(dimension%msh,atoms%ntype,dimension%jspd),qint(atoms%ntype,dimension%jspd)
   REAL tec(atoms%ntype,DIMENSION%jspd),rhTemp(dimension%msh,atoms%ntype,dimension%jspd)
   REAL chmom(atoms%ntype,dimension%jspd),clmom(3,atoms%ntype,dimension%jspd)
   INTEGER,ALLOCATABLE :: igq_fft(:)
   REAL   ,ALLOCATABLE :: qvac(:,:,:,:),qvlay(:,:,:,:,:)
   CHARACTER(LEN=20)   :: attributes(4)

   !pk non-collinear (start)
   REAL    rhoint,momint,alphdiff(atoms%ntype)
   INTEGER igq2_fft(0:stars%kq1_fft*stars%kq2_fft-1)
   COMPLEX,ALLOCATABLE :: qa21(:)
   !pk non-collinear (end)

   iter = inIter
   CALL outDen%init(stars,atoms,sphhar,vacuum,oneD,input%jspins,.FALSE.,POTDEN_TYPE_DEN)
   archiveType = CDN_ARCHIVE_TYPE_CDN1_const
   IF(noco%l_noco) archiveType = CDN_ARCHIVE_TYPE_NOCO_const

   IF (mpi%irank.EQ.0) THEN
      INQUIRE(file='enpara',exist=l_enpara)
      IF (l_enpara) OPEN (40,file ='enpara',form = 'formatted',status ='unknown')
   ENDIF
   ALLOCATE (outDen%cdom(stars%ng3),outDen%cdomvz(vacuum%nmzd,2))
   ALLOCATE (outDen%cdomvxy(vacuum%nmzxyd,oneD%odi%n2d-1,2))
   ALLOCATE (qa21(atoms%ntype))
   ALLOCATE (qvac(dimension%neigd,2,kpts%nkpt,dimension%jspd))
   ALLOCATE (qvlay(dimension%neigd,vacuum%layerd,2,kpts%nkpt,dimension%jspd))
   ALLOCATE (igq_fft(0:stars%kq1_fft*stars%kq2_fft*stars%kq3_fft-1))

   !initialize density arrays with zero
   qa21(:) = cmplx(0.0,0.0)
   qvac(:,:,:,:) = 0.0 
   qvlay(:,:,:,:,:) = 0.0
   outDen%mt(:,:,:,:) = 0.0
   outDen%pw(:,:) = cmplx(0.0,0.0)
   outDen%cdom(:) =  cmplx(0.0,0.0)
   IF (input%film) THEN
      outDen%vacz(:,:,:) = 0.0
      outDen%cdomvz(:,:) = cmplx(0.0,0.0)
      outDen%vacxy(:,:,:,:) = cmplx(0.0,0.0)
      outDen%cdomvxy(:,:,:) = cmplx(0.0,0.0)
   END IF
   outDen%iter = iter
        
   !Set up pointer for backtransformation of from g-vector in
   !positive domain fof carge density fftibox into stars
   !In principle this can also be done in main program once.
   !It is done here to save memory.
   CALL prp_qfft_map(stars,sym, input, igq2_fft,igq_fft)

   !LDA+U: initialise density-matrix if needed
   ALLOCATE (outDen%mmpMat(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,MAX(1,atoms%n_u),input%jspins))
   outDen%mmpMat(:,:,:,:) = CMPLX(0.0,0.0)


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
                  sphhar,sym,obsolete,igq_fft,vTot%mt,vTot%vacz(:,:,jspin),oneD,coreSpecInput,&
                  outDen%mmpMat(-lmaxU_const:,-lmaxU_const:,:,jspin),results, outDen%pw,outDen%vacxy,outDen%mt,outDen%vacz,&
                  outDen%cdom,outDen%cdomvz,outDen%cdomvxy,qvac,qvlay,qa21, chmom,clmom)
      CALL timestop("cdngen: cdnval")
!-fo
   END DO

   ! lda+u
   IF ((atoms%n_u.GT.0).and.(mpi%irank.EQ.0)) CALL u_mix(input,atoms,inDen%mmpMat,outDen%mmpMat)

!+t3e
   IF (mpi%irank.EQ.0) THEN
!-t3e
      IF (l_enpara) CLOSE (40)
      CALL cdntot(stars,atoms,sym, vacuum,input,cell,oneD, outDen%pw,outDen%mt,outDen%vacz,.TRUE., qtot,dummy)
      CALL closeXMLElement('valenceDensity')
!---> changes
   END IF ! mpi%irank = 0
   IF (input%kcrel.EQ.0) THEN
      results%seigc = 0.

      ! Generate input file ecore for subsequent GW calculation
      ! 11.2.2004 Arno Schindlmayr
      IF ((input%gw.eq.1 .or. input%gw.eq.3).AND.(mpi%irank.EQ.0)) THEN
         OPEN (15,file='ecore',status='unknown', action='write',form='unformatted')
      END IF

      rh = 0.0
      tec = 0.0
      qint = 0.0
      IF (input%frcor) THEN
         IF (mpi%irank.EQ.0) THEN
            CALL readCoreDensity(input,atoms,dimension,rh,tec,qint)
         END IF
#ifdef CPP_MPI
         CALL mpi_bc_coreDen(mpi,atoms,input,dimension,rh,tec,qint)
#endif
      END IF

      DO jspin = 1,input%jspins
         IF ((input%jspins.EQ.2).AND.(mpi%irank.EQ.0)) THEN
            DO n = 1,atoms%ntype
               svdn(n,jspin) = outDen%mt(1,0,n,jspin)/ (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
            END DO
         END IF

         !block 1 unnecessary for slicing: begin
         IF (.NOT.sliceplot%slice) THEN

            !add in core density
            IF (mpi%irank.EQ.0) THEN
               CALL cored(input,jspin,atoms,outDen%mt,dimension,sphhar,vTot%mt(:,0,:,jspin), qint,rh,tec,seig)
               rhTemp(:,:,jspin) = rh(:,:,jspin)
               results%seigc = results%seigc + seig
               IF (input%jspins.EQ.2) THEN
                  DO n = 1,atoms%ntype
                     stdn(n,jspin) = outDen%mt(1,0,n,jspin)/ (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
                  END DO
               END IF
            END IF  ! mpi%irank = 0

            !add core tail charge to outDen%pw
            IF ((noco%l_noco).AND.(mpi%irank.EQ.0)) THEN
               !pk non-collinear (start)
               !add the coretail-charge to the constant interstitial
               !charge (star 0), taking into account the direction of
               !magnetisation of this atom
               IF (jspin .EQ. 2) THEN
                  DO ityp = 1,atoms%ntype
                     rhoint = (qint(ityp,1) + qint(ityp,2)) /cell%volint/input%jspins/2.0
                     momint = (qint(ityp,1) - qint(ityp,2)) /cell%volint/input%jspins/2.0
                     !rho_11
                     outDen%pw(1,1) = outDen%pw(1,1) + rhoint + momint*cos(noco%beta(ityp))
                     !rho_22
                     outDen%pw(1,2) = outDen%pw(1,2) + rhoint - momint*cos(noco%beta(ityp))
                     !real part rho_21
                     outDen%cdom(1) = outDen%cdom(1) + cmplx(0.5*momint *cos(noco%alph(ityp))*sin(noco%beta(ityp)),0.0)
                     !imaginary part rho_21
                     outDen%cdom(1) = outDen%cdom(1) + cmplx(0.0,-0.5*momint *sin(noco%alph(ityp))*sin(noco%beta(ityp)))
                  END DO
               END IF
               !pk non-collinear (end)

            ELSE IF (input%ctail) THEN
               CALL cdnovlp(mpi,sphhar,stars,atoms,sym, dimension,vacuum,&
                            cell, input,oneD,l_st, jspin,rh(:,:,jspin),&
                            outDen%pw,outDen%vacxy,outDen%mt,outDen%vacz)
            ELSE IF (mpi%irank.EQ.0) THEN
               DO ityp = 1,atoms%ntype
                  outDen%pw(1,jspin) = outDen%pw(1,jspin) + qint(ityp,jspin)/input%jspins/cell%volint
               END DO
            END IF
         !block 1 unnecessary for slicing: end
         END IF ! .NOT.sliceplot%slice

      END DO ! loop over spins
      IF (mpi%irank.EQ.0) THEN
         CALL writeCoreDensity(input,atoms,dimension,rhTemp,tec,qint)
      END IF
      IF ((input%gw.eq.1 .or. input%gw.eq.3).AND.(mpi%irank.EQ.0)) CLOSE(15)
   ELSE ! input%kcrel.EQ.0
      !relativistic core implementation : kcrel.eq.1
      results%seigc = 0.
      IF ((input%jspins.EQ.2).AND.(mpi%irank.EQ.0)) THEN
         DO jspin = 1,input%jspins
            DO n = 1,atoms%ntype
               svdn(n,jspin) = outDen%mt(1,0,n,jspin)/ (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
            END DO
         END DO
      END IF
      !block 1 unnecessary for slicing: begin
      IF (.NOT.sliceplot%slice) THEN
         !add in core density
         IF (mpi%irank.EQ.0) THEN
            CALL coredr(input,atoms,seig, outDen%mt,dimension,sphhar,vTot%mt(:,0,:,:),qint,rh)
            results%seigc = results%seigc + seig
            IF (input%jspins.EQ.2) THEN
               DO jspin = 1,input%jspins
                  DO n = 1,atoms%ntype
                     stdn(n,jspin) = outDen%mt(1,0,n,jspin)/ (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
                  END DO
               END DO
            END IF
         END IF

         IF ((noco%l_noco).AND.(mpi%irank.EQ.0)) THEN
            !pk non-collinear (start)
            !add the coretail-charge to the constant interstitial
            !charge (star 0), taking into account the direction of
            !magnetisation of this atom
            DO ityp = 1,atoms%ntype
               rhoint = (qint(ityp,1) + qint(ityp,2)) /cell%volint/input%jspins/2.0
               momint = (qint(ityp,1) - qint(ityp,2)) /cell%volint/input%jspins/2.0
               !rho_11
               outDen%pw(1,1) = outDen%pw(1,1) + rhoint + momint*cos(noco%beta(ityp))
               !rho_22
               outDen%pw(1,2) = outDen%pw(1,2) + rhoint - momint*cos(noco%beta(ityp))
               !real part rho_21
               outDen%cdom(1) = outDen%cdom(1) + cmplx(0.5*momint *cos(noco%alph(ityp))*sin(noco%beta(ityp)),0.0)
               !imaginary part rho_21
               outDen%cdom(1) = outDen%cdom(1) + cmplx(0.0,-0.5*momint *sin(noco%alph(ityp))*sin(noco%beta(ityp)))
            END DO
            !pk non-collinear (end)
         ELSE
            DO jspin = 1,input%jspins
               IF (input%ctail) THEN
!+gu hope this works as well
                  CALL cdnovlp(mpi,sphhar,stars,atoms,sym,dimension,vacuum,&
                               cell,input,oneD,l_st,jspin,rh(1,1,jspin),&
                               outDen%pw,outDen%vacxy,outDen%mt,outDen%vacz)
               ELSE IF (mpi%irank.EQ.0) THEN
                  DO ityp = 1,atoms%ntype
                     outDen%pw(1,jspin) = outDen%pw(1,jspin) + qint(ityp,jspin)/input%jspins/cell%volint
                  END DO
               END IF
            END DO
         END IF
      !block 1 unnecessary for slicing: end
      END IF ! .NOT.sliceplot%slice
! end relativistic core
   END IF ! input%kcrel.EQ.0

   IF (mpi%irank.EQ.0) THEN
      !block 2 unnecessary for slicing: begin
      IF (.NOT.sliceplot%slice) THEN
         CALL openXMLElementNoAttributes('allElectronCharges')
         CALL qfix(stars,atoms,sym,vacuum, sphhar,input,cell,oneD,&
                   outDen%pw,outDen%vacxy,outDen%mt,outDen%vacz,.TRUE.,.true.,fix)
         CALL closeXMLElement('allElectronCharges')
         !pk non-collinear (start)
         IF (noco%l_noco) THEN
            !fix also the off-diagonal part of the density matrix
            outDen%cdom(:stars%ng3) = fix*outDen%cdom(:stars%ng3)
            IF (input%film) THEN
               outDen%cdomvz(:,:) = fix*outDen%cdomvz(:,:)
               outDen%cdomvxy(:,:,:) = fix*outDen%cdomvxy(:,:,:)
            END IF
         END IF
         !pk non-collinear (end)

         !spin densities at the nucleus
         !and magnetic moment in the spheres
         IF (input%jspins.EQ.2) THEN
            WRITE (6,FMT=8000)
            WRITE (16,FMT=8000)
            DO n = 1,atoms%ntype
               sval = svdn(n,1) - svdn(n,input%jspins)
               stot = stdn(n,1) - stdn(n,input%jspins)
               scor = stot - sval
               WRITE (6,FMT=8010) n,stot,sval,scor,svdn(n,1),stdn(n,1)
               WRITE (16,FMT=8010) n,stot,sval,scor,svdn(n,1),stdn(n,1)
            END DO
            IF (noco%l_mperp) THEN
               ! angles in nocoinp file are (alph-alphdiff)
               iatom = 1
               DO n = 1,atoms%ntype
                  IF (noco%l_ss) THEN
                     alphdiff(n)= 2.*pi_const*(noco%qss(1)*atoms%taual(1,iatom) + &
                                               noco%qss(2)*atoms%taual(2,iatom) + &
                                               noco%qss(3)*atoms%taual(3,iatom) )
                  ELSE
                     alphdiff(n)= 0.
                  END IF
                  iatom= iatom + atoms%neq(n)
               END DO
            END IF
            WRITE (6,FMT=8020)
            WRITE (16,FMT=8020)
            noco_new = noco
            CALL openXMLElement('magneticMomentsInMTSpheres',(/'units'/),(/'muBohr'/))
            DO n = 1, atoms%ntype
               smom = chmom(n,1) - chmom(n,input%jspins)
               WRITE (6,FMT=8030) n,smom, (chmom(n,j),j=1,input%jspins)
               WRITE (16,FMT=8030) n,smom, (chmom(n,j),j=1,input%jspins)
               attributes = ''
               WRITE(attributes(1),'(i0)') n
               WRITE(attributes(2),'(f15.10)') smom
               WRITE(attributes(3),'(f15.10)') chmom(n,1)
               WRITE(attributes(4),'(f15.10)') chmom(n,2)
               CALL writeXMLElementFormPoly('magneticMoment',(/'atomType      ','moment        ','spinUpCharge  ',&
                                                               'spinDownCharge'/),&
                                            attributes,reshape((/8,6,12,14,6,15,15,15/),(/4,2/)))
               IF (noco%l_mperp) THEN
                  !calculate the perpendicular part of the local moment
                  !and relax the angle of the local moment or calculate
                  !the constraint B-field.
                  CALL m_perp(atoms,n,noco_new,vTot%mt(:,0,:,:),chmom,qa21,alphdiff)
               END IF
            END DO
            CALL closeXMLElement('magneticMomentsInMTSpheres')

            !save the new nocoinp file if the dierctions of the local
            !moments are relaxed or a constraint B-field is calculated.
            l_relax_any = .false.
            iatom = 1
            DO itype = 1,atoms%ntype
               l_relax_any = l_relax_any.OR.noco%l_relax(itype)
            END DO
            IF (l_relax_any.OR.noco%l_constr) THEN
               IF (.not. noco%l_mperp) THEN
                  CALL juDFT_error ("(l_relax_any.OR.noco).AND.(.NOT. )" ,calledby ="cdngen")
               END IF
               DO itype = 1, atoms%ntype
                  IF (noco%l_ss) THEN
                     noco_new%alph(itype) = noco%alph(itype) - alphdiff(itype)
                     DO WHILE (noco_new%alph(n) > +pi_const)
                        noco_new%alph(n)= noco_new%alph(n) - 2.*pi_const
                     END DO
                     DO WHILE (noco_new%alph(n) < -pi_const)
                        noco_new%alph(n)= noco_new%alph(n) + 2.*pi_const
                     END DO
                  ELSE
                     noco_new%alph(itype) = noco%alph(itype)
                  END IF
               END DO

               OPEN (24,file='nocoinp',form='formatted', status='old')
               REWIND (24)
               CALL rw_noco_write(atoms,jij,noco_new, input)
               CLOSE (24)
            END IF

            IF (noco%l_soc) THEN
               thetai = noco%theta
               phii   = noco%phi
               WRITE (6,FMT=9020)
               WRITE (16,FMT=9020)
               CALL openXMLElement('orbitalMagneticMomentsInMTSpheres',(/'units'/),(/'muBohr'/))
               DO n = 1, atoms%ntype
                  IF (noco%l_noco) THEN
                     thetai = noco%beta(n)
                     phii =   noco%alph(n)
                  END IF

                  ! magn. moment(-)
                  slxmom = clmom(1,n,1)+clmom(1,n,2)
                  slymom = clmom(2,n,1)+clmom(2,n,2)
                  slmom =  clmom(3,n,1)+clmom(3,n,2)

                  ! rotation: orbital moment || spin moment (extended to incude phi - hopefully)
                  slmom   = cos(thetai)*slmom + sin(thetai)* (cos(phii)*slxmom + sin(phii)*slymom)
                  clmom(3,n,1) = cos(thetai)*clmom(3,n,1) + &
                                 sin(thetai)*(cos(phii)*clmom(1,n,1) + sin(phii)*clmom(2,n,1))
                  clmom(3,n,2) = cos(thetai)*clmom(3,n,2) + &
                                 sin(thetai)*(cos(phii)*clmom(1,n,2) + sin(phii)*clmom(2,n,2))

                  WRITE (6,FMT=8030) n,slmom,(clmom(3,n,j),j=1,2)
                  WRITE (16,FMT=8030) n,slmom,(clmom(3,n,j),j=1,2)
                  attributes = ''
                  WRITE(attributes(1),'(i0)') n
                  WRITE(attributes(2),'(f15.10)') slmom
                  WRITE(attributes(3),'(f15.10)') clmom(3,n,1)
                  WRITE(attributes(4),'(f15.10)') clmom(3,n,2)
                  CALL writeXMLElementFormPoly('orbMagMoment',(/'atomType      ','moment        ','spinUpCharge  ',&
                                                                'spinDownCharge'/),&
                                               attributes,reshape((/8,6,12,14,6,15,15,15/),(/4,2/)))
               END DO
               CALL closeXMLElement('orbitalMagneticMomentsInMTSpheres')
            END IF
         END IF
      !block 2 unnecessary for slicing: end
      END IF ! .NOT.sliceplot%slice

      9020 FORMAT (/,/,10x,'orb. magnetic moments in the spheres:',/,10x,&
                   'type',t22,'moment',t33,'spin-up',t43,'spin-down')
      8000 FORMAT (/,/,10x,'spin density at the nucleus:',/,10x,'type',t25,&
                   'input%total',t42,'valence',t65,'core',t90,&
                   'majority valence and input%total density',/)
      8010 FORMAT (i13,2x,3e20.8,5x,2e20.8)
      8020 FORMAT (/,/,2x,'-->  magnetic moments in the spheres:',/,2x,&
                   'mm -->   type',t22,'moment',t33,'spin-up',t43,'spin-down')
      8030 FORMAT (2x,'--> mm',i8,2x,3f12.5)

      IF (sliceplot%slice) THEN
         OPEN (20,file='cdn_slice',form='unformatted',status='unknown')
         CALL wrtdop(stars,vacuum,atoms,sphhar, input,sym, 20, outDen%iter,outDen%mt,outDen%pw,outDen%vacz,outDen%vacxy)
         IF (noco%l_noco) THEN
            WRITE (20) (outDen%cdom(k),k=1,stars%ng3)
            IF (input%film) THEN
               WRITE (20) ((outDen%cdomvz(j,ivac),j=1,vacuum%nmz),ivac=1,vacuum%nvac)
               WRITE (20) (((outDen%cdomvxy(j,k-1,ivac),j=1,vacuum%nmzxy),k=2,oneD%odi%nq2) ,ivac=1,vacuum%nvac)
            END IF
         END IF
         CLOSE(20) 
         CALL juDFT_end("slice OK")
      END IF

      CALL writeDensity(stars,vacuum,atoms,cell,sphhar,input,sym,oneD,archiveType,&
                        CDN_OUTPUT_DEN_const,0,results%last_distance,results%ef,.FALSE.,outDen)
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
