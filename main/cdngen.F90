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
   USE m_types
   USE m_xmlOutput
   USE m_magMoms
   USE m_orbMagMoms
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
   REAL fix,qtot,scor,seig,stot,sval,dummy
   REAL sum,fermiEnergyTemp
   INTEGER iter,ivac,j,jspin,jspmax,k,n,nt,ieig,ikpt
   INTEGER  ityp,ilayer,urec,itype,iatom
   LOGICAL l_relax_any,exst,n_exist,l_qfix, l_enpara
   LOGICAL, PARAMETER :: l_st=.FALSE.

   !Local Arrays
   REAL stdn(atoms%ntype,dimension%jspd),svdn(atoms%ntype,dimension%jspd),alpha_l(atoms%ntype)
   REAL rh(dimension%msh,atoms%ntype,dimension%jspd),qint(atoms%ntype,dimension%jspd)
   REAL tec(atoms%ntype,DIMENSION%jspd),rhTemp(dimension%msh,atoms%ntype,dimension%jspd)
   REAL chmom(atoms%ntype,dimension%jspd),clmom(3,atoms%ntype,dimension%jspd)
   INTEGER,ALLOCATABLE :: igq_fft(:)
   REAL   ,ALLOCATABLE :: qvac(:,:,:,:),qvlay(:,:,:,:,:)

   !pk non-collinear (start)
   REAL    rhoint,momint,alphdiff
   INTEGER igq2_fft(0:stars%kq1_fft*stars%kq2_fft-1)
   COMPLEX,ALLOCATABLE :: qa21(:), cdomvz(:,:)
   !pk non-collinear (end)

   iter = inIter
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

   outDen%iter = iter
        
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
                  sphhar,sym,obsolete,igq_fft,vTot%mt,vTot%vacz(:,:,jspin),oneD,coreSpecInput,&
                  outDen,results,qvac,qvlay,qa21, chmom,clmom)
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

   results%seigc = 0.0
   IF (mpi%irank.EQ.0) THEN
      DO jspin = 1,input%jspins
         DO n = 1,atoms%ntype
            svdn(n,jspin) = outDen%mt(1,0,n,jspin)/ (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
         END DO
      END DO
   END IF

   IF (input%kcrel.EQ.0) THEN
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
   END IF

   IF (.NOT.sliceplot%slice) THEN
      !add in core density
      IF (mpi%irank.EQ.0) THEN
         IF (input%kcrel.EQ.0) THEN
            DO jspin = 1,input%jspins
               CALL cored(input,jspin,atoms,outDen%mt,dimension,sphhar,vTot%mt(:,0,:,jspin), qint,rh,tec,seig)
               rhTemp(:,:,jspin) = rh(:,:,jspin)
               results%seigc = results%seigc + seig
            END DO
         ELSE
            CALL coredr(input,atoms,seig, outDen%mt,dimension,sphhar,vTot%mt(:,0,:,:),qint,rh)
            results%seigc = results%seigc + seig
         END IF
      END IF
      DO jspin = 1,input%jspins
         IF (mpi%irank.EQ.0) THEN
            DO n = 1,atoms%ntype
               stdn(n,jspin) = outDen%mt(1,0,n,jspin)/ (sfp_const*atoms%rmsh(1,n)*atoms%rmsh(1,n))
            END DO
         END IF
         IF ((noco%l_noco).AND.(mpi%irank.EQ.0)) THEN
            IF (jspin.EQ.2) THEN
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
                  outDen%pw(1,3) = outDen%pw(1,3) + cmplx(0.5*momint *cos(noco%alph(ityp))*sin(noco%beta(ityp)),0.0)
                  !imaginary part rho_21
                  outDen%pw(1,3) = outDen%pw(1,3) + cmplx(0.0,-0.5*momint *sin(noco%alph(ityp))*sin(noco%beta(ityp)))
               END DO
               !pk non-collinear (end)
            END IF
         ELSE
            IF (input%ctail) THEN
               !+gu hope this works as well
               CALL cdnovlp(mpi,sphhar,stars,atoms,sym,dimension,vacuum,&
                            cell,input,oneD,l_st,jspin,rh(:,:,jspin),&
                            outDen%pw,outDen%vacxy,outDen%mt,outDen%vacz)
            ELSE IF (mpi%irank.EQ.0) THEN
               DO ityp = 1,atoms%ntype
                  outDen%pw(1,jspin) = outDen%pw(1,jspin) + qint(ityp,jspin)/input%jspins/cell%volint
               END DO
            END IF
         END IF
      END DO
   END IF

   IF (input%kcrel.EQ.0) THEN
      IF (mpi%irank.EQ.0) THEN
         CALL writeCoreDensity(input,atoms,dimension,rhTemp,tec,qint)
      END IF
      IF ((input%gw.eq.1 .or. input%gw.eq.3).AND.(mpi%irank.EQ.0)) CLOSE(15)
   END IF

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
            DO n = 1,atoms%ntype
               sval = svdn(n,1) - svdn(n,input%jspins)
               stot = stdn(n,1) - stdn(n,input%jspins)
               scor = stot - sval
               WRITE (6,FMT=8010) n,stot,sval,scor,svdn(n,1),stdn(n,1)
               WRITE (16,FMT=8010) n,stot,sval,scor,svdn(n,1),stdn(n,1)
            END DO

            noco_new = noco

            CALL magMoms(dimension,input,atoms,noco_new,vTot,chmom,qa21)

            !save the new nocoinp file if the dierctions of the local
            !moments are relaxed or a constraint B-field is calculated.
            l_relax_any = .false.
            DO itype = 1,atoms%ntype
               l_relax_any = l_relax_any.OR.noco%l_relax(itype)
            END DO
            IF (l_relax_any.OR.noco%l_constr) THEN
               IF (.not. noco%l_mperp) THEN
                  CALL juDFT_error ("(l_relax_any.OR.noco).AND.(.NOT. )" ,calledby ="cdngen")
               END IF
               iatom = 1
               DO itype = 1, atoms%ntype
                  IF (noco%l_ss) THEN
                     alphdiff = 2.0*pi_const*(noco%qss(1)*atoms%taual(1,iatom) + &
                                              noco%qss(2)*atoms%taual(2,iatom) + &
                                              noco%qss(3)*atoms%taual(3,iatom) )
                     noco_new%alph(itype) = noco%alph(itype) - alphdiff
                     DO WHILE (noco_new%alph(n) > +pi_const)
                        noco_new%alph(n)= noco_new%alph(n) - 2.*pi_const
                     END DO
                     DO WHILE (noco_new%alph(n) < -pi_const)
                        noco_new%alph(n)= noco_new%alph(n) + 2.*pi_const
                     END DO
                  ELSE
                     noco_new%alph(itype) = noco%alph(itype)
                  END IF
                  iatom= iatom + atoms%neq(n)
               END DO

               OPEN (24,file='nocoinp',form='formatted', status='old')
               REWIND (24)
               CALL rw_noco_write(atoms,jij,noco_new, input)
               CLOSE (24)
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
