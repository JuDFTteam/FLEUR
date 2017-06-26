!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     This module generates the cmt coefficients and eigenvectors z   !
!     at all kpoints nkpt from the irreducible kpoints nkpti          !
!     and writes them out in cmt and z, respectively.                 !
!                                                 M.Betzinger(09/07)  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      MODULE m_gen_wavf

      CONTAINS

      SUBROUTINE gen_wavf (&
     &          nkpti,kpts,it,sym,&
     &          atoms,el_eig,ello_eig,cell,&
     &          dimension,hybrid,vr0,&
     &          hybdat,&
     &          noco,oneD,mpi,irank2,&
     &          input,jsp,&
     &          zmat)


      ! nkpti      ::     number of irreducible k-points
      ! nkpt       ::     number of all k-points 

      USE m_apws
      USE m_radfun
      USE m_radflo
      USE m_abcof
      USE m_trafo     ,ONLY: waveftrafo_genwavf
      USE m_util      ,ONLY: modulo1
      USE m_setabc1locdn
      USE m_olap
      USE m_types
      USE m_abcrot
      IMPLICIT NONE

      TYPE(t_hybdat),INTENT(INOUT)   :: hybdat
      TYPE(t_mpi),INTENT(IN)         :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_oneD),INTENT(IN)        :: oneD
      TYPE(t_hybrid),INTENT(IN)      :: hybrid
      TYPE(t_input),INTENT(IN)       :: input
      TYPE(t_noco),INTENT(IN)        :: noco
      TYPE(t_sym),INTENT(IN)         :: sym
      TYPE(t_cell),INTENT(IN)        :: cell
      TYPE(t_kpts),INTENT(IN)        :: kpts
      TYPE(t_atoms),INTENT(IN)       :: atoms

!     - - scalars - -
      INTEGER,INTENT(IN)      ::  nkpti ,it    

!     scalars for abcoff

!     scalars for apws
      INTEGER,INTENT(IN)      :: jsp

!     - - arrays - -
      INTEGER,INTENT(IN)      ::  irank2(nkpti)


      REAL,INTENT(IN)         ::  vr0(:,:,:)!(jmtd,ntype,jspd)
      
      REAL,INTENT(IN)         ::  el_eig(0:atoms%lmaxd,atoms%ntype)
      REAL,INTENT(IN)         ::  ello_eig(atoms%nlod,atoms%ntype)
      TYPE(t_zmat),INTENT(IN) :: zmat(:) !for all kpoints 

  !     - - local scalars - - 
      INTEGER                 ::  ilo,idum ,m,maxlmindx
      REAL                    ::  rdum,merror
      COMPLEX                 ::  cdum,cdum1,cdum2

!     local scalars for apws
      INTEGER                 ::  nred
      INTEGER                 ::  ikpt0,ikpt,ikpt1,iband,itype,iop,&
     &                            ispin,ieq,ic,indx,iatom
      INTEGER                 ::  i,j,l ,ll,lm,nrkpt,ng,ok
      REAL                    ::  time1,time2,time3,time4
      COMPLEX                 ::  img=(0d0,1d0)

!     local scalars for radfun
      INTEGER                 ::  nodem,noded
      REAL                    ::  wronk

!     reduced dimension for parallel calculation
      INTEGER                 ::  lower, upper
      LOGICAL                 ::  found

!     - - local arrays - -
      INTEGER                 ::  rrot(3,3,sym%nsym)
      INTEGER                 ::  map_lo(atoms%nlod)
      INTEGER                 ::  iarr(0:atoms%lmaxd,atoms%ntype)
      REAL                    ::  rotkpt(3)
      COMPLEX,ALLOCATABLE     ::  acof(:,:,:),bcof(:,:,:),ccof(:,:,:,:)
      REAL,ALLOCATABLE        ::  zhlp_r(:,:)
      COMPLEX,ALLOCATABLE     ::  zhlp_c(:,:)

      COMPLEX,ALLOCATABLE     ::  cmt(:,:,:),cmthlp(:,:,:)


!     local arrays for abcof1
!       COMPLEX                 ::  a(nvd,0:lmd,natd,nkpti),b(nvd,0:lmd,natd,nkpti)

!     local arrays for radfun
      REAL                    ::  vr(atoms%jmtd,atoms%ntype,dimension%jspd)
      REAL,ALLOCATABLE        ::  f(:,:,:),df(:,:,:)


!     local arrays for radflo
      REAL                    ::  flo(atoms%jmtd,2,atoms%nlod)
      REAL                    ::  uloulopn(atoms%nlod,atoms%nlod,atoms%ntype)
      REAL                    ::  uuilon(atoms%nlod,atoms%ntype),duilon(atoms%nlod,atoms%ntype)
      REAL                    ::  ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)

!     local arrays for setabc1locdn
      INTEGER                 ::  nbasf0(atoms%nlod,atoms%nat),nkvec(atoms%nlod,atoms%nat)
      INTEGER                 ::  kvec(2*(2*atoms%llod+1)  )
      LOGICAL                 ::  enough(atoms%nat)


!     local arrays for apws
      INTEGER                 ::  matind(dimension%nbasfcn,2)
      INTEGER                 ::  gpthlp1(3,dimension%nvd,dimension%jspd),nvhlp1(dimension%jspd)
      INTEGER                 ::  gpthlp2(3,dimension%nvd,dimension%jspd),nvhlp2(dimension%jspd)
      INTEGER                 ::  k1hlp(dimension%nvd,dimension%jspd),k2hlp(dimension%nvd,dimension%jspd),&
     &                            k3hlp(dimension%nvd,dimension%jspd),nvhlp(dimension%jspd)
      REAL                    :: bkpt(3)

      REAL,ALLOCATABLE        ::  olapmt(:,:,:,:)
#if ( defined(CPP_INVERSION) )
      REAL,ALLOCATABLE        ::  olappw(:,:)
#else
      COMPLEX,ALLOCATABLE     ::  olappw(:,:)
#endif
      INTEGER                 ::  irecl_cmt,irecl_z

      INTEGER                 ::  gpt(3,dimension%nvd,dimension%jspd,kpts%nkpt),ngpt(dimension%jspd,kpts%nkpt)

      TYPE(t_lapw) :: lapw
      TYPE(t_usdus):: usdus
      
      CALL CPU_TIME(time1)

      ! setup rotations in reciprocal space
      DO iop=1,sym%nsym
        IF( iop .le. sym%nop ) THEN
          rrot(:,:,iop) = transpose( sym%mrot(:,:,sym%invtab(iop)) )
        ELSE
          rrot(:,:,iop) = -rrot(:,:,iop-sym%nop)
        END IF
      END DO

      ! generate G-vectors, which fulfill |k+G|<rkmax
      ! for all k-points
      DO ikpt=1,kpts%nkpt
        CALL apws(dimension,input,noco,&
     &            kpts,ikpt,cell,sym%zrfs,&
     &            1,jsp,bkpt,lapw,matind,nred)

      END DO

      ! read in spherical component of the potential from
      ! the previous iteration
      ! eigen_hf writes spherical component in file vr0
      ! if it = 1 vr0 is identical to the potential of the
      ! previous calculation

      IF( it .ne. 1) THEN
        OPEN(unit=220,file='vr0',form='unformatted')
        DO ispin=1,dimension%jspd
          DO itype=1,atoms%ntype
            DO i=1,atoms%jmtd
              READ(220) vr(i,itype,ispin)
            END DO
          END DO
        END DO
        CLOSE(220)
      ELSE
        vr = vr0
      END IF


!       ALLOCATE ( z_out(nbasfcn,neigd,nkpti),stat=ok )
!       IF ( ok .ne. 0) STOP 'gen_wavf: failure allocation z'
!       z_out = 0
!       z_out(:,:,:nkpti) = z_in


      ! calculate radial basis functions belonging to the
      ! potential vr stored in bas1 and bas2
      ! bas1 denotes the large component
      ! bas2    "     "  small component

      ALLOCATE( f(atoms%jmtd,2,0:atoms%lmaxd), df(atoms%jmtd,2,0:atoms%lmaxd) )
      f    = 0
      df   = 0
      iarr = 2
      DO itype=1,atoms%ntype
        if ( mpi%irank == 0 ) WRITE (6,FMT=8000) itype
        ng = atoms%jri(itype)
        DO l=0,atoms%lmax(itype)
           CALL radfun(l,itype,1,el_eig(l,itype),vr(:,itype,jsp),&
                atoms,f(:,:,l),df(:,:,l),usdus,nodem,noded,wronk)
          IF ( mpi%irank == 0 ) WRITE (6,FMT=8010) l,el_eig(l,itype),&
               usdus%us(l,itype,1),usdus%dus(l,itype,1),nodem,usdus%uds(l,itype,1),&
               usdus%duds(l,itype,1),noded,usdus%ddn(l,itype,1),wronk

          hybdat%bas1(1:ng,1,l,itype) =  f(1:ng,1,l)
          hybdat%bas2(1:ng,1,l,itype) =  f(1:ng,2,l)
          hybdat%bas1(1:ng,2,l,itype) = df(1:ng,1,l)
          hybdat%bas2(1:ng,2,l,itype) = df(1:ng,2,l)

            hybdat%bas1_MT(1,l,itype) =   usdus%us(l,itype,1)
          hybdat%drbas1_MT(1,l,itype) =  usdus%dus(l,itype,1)
            hybdat%bas1_MT(2,l,itype) =  usdus%uds(l,itype,1)
          hybdat%drbas1_MT(2,l,itype) = usdus%duds(l,itype,1)
        END DO

        IF (atoms%nlo(itype).GE.1) THEN
          CALL radflo( atoms,itype,jsp,&
     &                 ello_eig,vr(:,itype,jsp),&
     &                 f,df,mpi,&
     &                 usdus,&
     &                 uuilon,duilon,ulouilopn,flo)

          DO ilo=1,atoms%nlo(itype)
            iarr(atoms%llo(ilo,itype),itype) = iarr(atoms%llo(ilo,itype),itype) + 1
            hybdat%bas1(1:ng,iarr(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),itype)&
     &                                                = flo(1:ng,1,ilo)
            hybdat%bas2(1:ng,iarr(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),itype)&
     &                                                = flo(1:ng,2,ilo)

              hybdat%bas1_MT(iarr(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),itype)&
     &                                                = usdus%ulos(ilo,itype,1)
            hybdat%drbas1_MT(iarr(atoms%llo(ilo,itype),itype),atoms%llo(ilo,itype),itype)&
     &                                                = usdus%dulos(ilo,itype,1)
          END DO

        END IF
      END DO
      DEALLOCATE (f,df)
#if CPP_DEBUG
      ! consistency check
      IF( .not. all(iarr .eq. hybrid%nindx ) ) STOP 'gen_wavf: counting error'
#endif

 8000    FORMAT (1x,/,/,' wavefunction parameters for atom type',i3,':',&
     &          /,t32,'radial function',t79,'energy derivative',/,t3,&
     &          'l',t8,'energy',t26,'value',t39,'derivative',t53,&
     &          'nodes',t68,'value',t81,'derivative',t95,'nodes',t107,&
     &          'norm',t119,'wronskian')
 8010    FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)


      ! determine boundaries for parallel calculations
      lower = 1
      upper = nkpti
      found = .false.
#     ifdef CPP_MPI
        DO ikpt = 1, nkpti
          IF ( irank2(ikpt) == 0 .AND. .not. found ) THEN
            lower = ikpt
            found = .true.
          ELSE IF ( irank2(ikpt) /= 0 .AND. found ) THEN
            upper = ikpt-1
            EXIT
          END IF
        END DO
#     else
        found = .true.
#     endif
      IF ( .not. found ) THEN
        upper = 0
      END IF


      maxlmindx = maxval(&
     &       (/ ( sum( (/ (hybrid%nindx(l,itype)*(2*l+1),l=0,atoms%lmax(itype)) /) ),&
     &                                              itype=1,atoms%ntype) /) )

      ! calculate wavefunction expansion in the the MT region
      ! (acof,bcof,ccof) and APW-basis coefficients
      ! (a,b,bascofold_lo) at irred. kpoints

!       CALL cpu_time(time2)
!       WRITE(*,*) 'time for generating radial functions',time2-time1

      ALLOCATE( acof(dimension%neigd,0:dimension%lmd,atoms%nat),stat=ok )
      IF( ok .ne. 0 ) STOP 'gen_wavf: failure allocation acof'
      ALLOCATE( bcof(dimension%neigd,0:dimension%lmd,atoms%nat),stat=ok )
      IF( ok .ne. 0 ) STOP 'gen_wavf: failure allocation bcof'
      ALLOCATE( ccof(-atoms%llod:atoms%llod,dimension%neigd,atoms%nlod,atoms%nat),stat=ok )
      IF( ok .ne. 0 ) STOP 'gen_wavf: failure allocation ccof'
      ALLOCATE ( cmt(dimension%neigd,maxlmindx,atoms%nat), stat=ok)
      IF(  ok .ne. 0 ) STOP 'gen_wavf: Failure allocation cmt'
      ALLOCATE ( cmthlp(dimension%neigd,maxlmindx,atoms%nat), stat=ok)
      IF( ok .ne. 0) STOP 'gen_wavf: failure allocation cmthlp'
      if (zmat(1)%l_real) THEN
         ALLOCATE ( zhlp_r(dimension%nbasfcn,dimension%neigd), stat=ok)
      ELSE
         ALLOCATE ( zhlp_c(dimension%nbasfcn,dimension%neigd), stat=ok)
      ENDIF
      IF( ok .ne. 0) STOP 'gen_wavf: failure allocation zhlp'

      irecl_cmt = dimension%neigd*maxlmindx*atoms%nat*16
      OPEN(unit=777,file='cmt',form='unformatted',access='direct',&
     &     recl=irecl_cmt)

#     ifdef CPP_INVERSION
        irecl_z   =  dimension%nbasfcn*dimension%neigd*8
#     else
        irecl_z   =  dimension%nbasfcn*dimension%neigd*16
#     endif
      OPEN(unit=778,file='z',form='unformatted',access='direct',&
     &     recl=irecl_z)

      DO ikpt0 = lower, upper

        acof = 0; bcof = 0; ccof = 0

        ! abcof calculates the wavefunction coefficients
        ! stored in acof,bcof,ccof
        CALL abcof(&
              input,atoms,hybdat%nbands(ikpt0),sym, cell, Kpts%bk(:,ikpt0), lapw, &
              ngpt(jsp,ikpt0)+hybdat%nbands(ikpt0),usdus,noco,jsp,hybdat%kveclo_eig(:,ikpt0),&
              oneD,acof(: hybdat%nbands(ikpt0),:,:),bcof(: hybdat%nbands(ikpt0),:,:),ccof(:,: hybdat%nbands(ikpt0),:,:),&
              zmat(ikpt))
        

! call was ...
          ! gpt(1,:,:,ikpt0),gpt(2,:,:,ikpt0),&
          ! gpt(3,:,:,ikpt0),ngpt(:,ikpt0),&!k1hlp,k2hlp,k3hlp,nvhlp,&
          !    ngpt(jsp,ikpt0)+nbands(ikpt0),z(:,:,ikpt0),&!nvhlp(jsp)+ &
          !   &usdus,&
          !    noco,&
          !    jsp,kveclo_eig(:ikpt0),oneD,oneD,&
          !    acof(:nbands(ikpt0),:,:),&
          !    bcof(:nbands(ikpt0),:,:),ccof(:,:nbands(ikpt0),:,:) )

        ! MT wavefunction coefficients are calculated in a local coordinate system
        ! rotate them in the global one

        CALL abcrot(&
                atoms,hybdat%nbands(ikpt0),&
                 sym,&
                cell,oneD,&
                acof(: hybdat%nbands(ikpt0),:,:),bcof(: hybdat%nbands(ikpt0),:,:),&
                ccof(:,: hybdat%nbands(ikpt0),:,:) )

!       CALL cpu_time(time3)

!       WRITE(*,*) 'time for abcoff and abcrot',time3-time2


!       CALL cpu_time(time2)


        ! decorate acof, bcof, ccof with coefficient i**l and store them
        ! in the field cmt(neigd,nkpt,maxlmindx,nat), i.e.
        ! where maxlmindx subsumes l,m and nindx

        cmt   = 0
        iatom = 0
        DO itype=1,atoms%ntype
          DO ieq=1,atoms%neq(itype)
            iatom   = iatom + 1
            indx    = 0
            DO l=0,atoms%lmax(itype)
              ll   = l*(l+1)
              cdum = img**l

              ! determine number of local orbitals with quantum number l
              ! map returns the number of the local orbital of quantum
              ! number l in the list of all local orbitals of the atom type
              idum   = 0
              map_lo = 0
              IF( hybrid%nindx(l,itype) .gt.2) THEN
                DO j=1,atoms%nlo(itype)
                  IF( atoms%llo(j,itype) .eq. l) THEN
                    idum = idum + 1
                    map_lo(idum) = j
                  END IF
                END DO
              END IF



              DO M=-l,l
                lm=ll+M
                DO i=1,hybrid%nindx(l,itype)
                  indx = indx + 1
                  IF( i .eq. 1 ) THEN
                    cmt(:,indx,iatom) =cdum*acof(:,lm,iatom)
                  ELSEIF ( i .eq. 2 ) THEN
                    cmt(:,indx,iatom) =cdum*bcof(:,lm,iatom)
                  ELSE
                    idum = i - 2
                    cmt(:,indx,iatom) =cdum*ccof(M,:,map_lo(idum),iatom)
                  END IF
                END DO
              END DO
            END DO

          END DO
        END DO

!       CALL cpu_time(time3)

        ! write cmt at irreducible k-points in direct-access file cmt
        WRITE(777,rec=ikpt0) cmt

        ! write z at irreducible k-points in direct-access file z
        IF (zmat(1)%l_real) THEN
           WRITE(778,rec=ikpt0) zmat(ikpt0)%z_r(:,:)
        ELSE
           WRITE(778,rec=ikpt0) zmat(ikpt0)%z_c(:,:)
        ENDIF
       

        ! generate wavefunctions coefficients at all k-points from
        ! irreducible k-points

        DO ikpt=1,kpts%nkpt
          IF ( kpts%bkp(ikpt) .eq. ikpt0 .and. ikpt0 .ne. ikpt ) THEN
            iop = kpts%bksym(ikpt)


 
            
              CALL waveftrafo_genwavf( cmthlp,zhlp_r,zhlp_c,&
     &                 cmt(:,:,:),zmat(1)%l_real,zmat(ikpt0)%z_r(:,:),zmat(ikpt0)%z_c(:,:),ikpt0,iop,atoms,&
     &                 hybrid,kpts,sym,&
     &                 jsp,dimension,hybdat%nbands(ikpt0),&
     &                 cell,gpt(:,:ngpt(jsp,ikpt0),jsp,ikpt0),&
     &                 ngpt(:,ikpt0),gpt(:,:ngpt(jsp,ikpt),jsp,ikpt),&
     &                 ngpt(:,ikpt),.true.)

              WRITE(777,rec=ikpt) cmthlp
              IF (zmat(1)%l_real) THEN
                 WRITE(778,rec=ikpt) zhlp_r
              ELSE
                 WRITE(778,rec=ikpt) zhlp_c
              ENDIF
          END IF
        END DO  !ikpt
      END DO !ikpt0

      DEALLOCATE( acof,bcof,ccof )
      DEALLOCATE( cmt,cmthlp)

      !close file cmt and z
      CLOSE(777)
      CLOSE(778)

      END SUBROUTINE gen_wavf

      END MODULE m_gen_wavf
