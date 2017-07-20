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
      USE m_olap
      USE m_types
      USE m_abcrot
      USE m_io_hybrid
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
      INTEGER                 ::  ilo,idum ,m
      COMPLEX                 ::  cdum
      TYPE(t_mat)             :: zhlp
!     local scalars for apws
      INTEGER                 ::  nred
      INTEGER                 ::  ikpt0,ikpt,itype,iop,ispin,ieq,indx,iatom
      INTEGER                 ::  i,j,l ,ll,lm,ng,ok
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
      COMPLEX,ALLOCATABLE     ::  acof(:,:,:),bcof(:,:,:),ccof(:,:,:,:)
      
      COMPLEX,ALLOCATABLE     ::  cmt(:,:,:),cmthlp(:,:,:)


!     local arrays for abcof1
!       COMPLEX                 ::  a(nvd,0:lmd,natd,nkpti),b(nvd,0:lmd,natd,nkpti)

!     local arrays for radfun
      REAL                    ::  vr(atoms%jmtd,atoms%ntype,dimension%jspd)
      REAL,ALLOCATABLE        ::  f(:,:,:),df(:,:,:)


!     local arrays for radflo
      REAL                    ::  flo(atoms%jmtd,2,atoms%nlod)
      REAL                    ::  uuilon(atoms%nlod,atoms%ntype),duilon(atoms%nlod,atoms%ntype)
      REAL                    ::  ulouilopn(atoms%nlod,atoms%nlod,atoms%ntype)
   

!     local arrays for apws
      INTEGER                 ::  matind(dimension%nbasfcn,2)
      REAL                    :: bkpt(3) 

      INTEGER                 ::  irecl_cmt,irecl_z


      TYPE(t_lapw) :: lapw(kpts%nkptf)
      TYPE(t_usdus):: usdus

      !CALL CPU_TIME(time1)
      call usdus%init(atoms,dimension%jspd)
      call zhlp%alloc(zmat(1)%l_real,zmat(1)%nbasfcn,zmat(1)%nbands)

      
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
      DO ikpt=1,kpts%nkptf
        CALL apws(dimension,input,noco,&
     &            kpts,ikpt,cell,sym%zrfs,&
     &            1,jsp,bkpt,lapw(ikpt),matind,nred)

      END DO

      ! read in spherical component of the potential from
      ! the previous iteration
      ! eigen_hf writes spherical component in file vr0
      ! if it = 1 vr0 is identical to the potential of the
      ! previous calculation
      vr = vr0


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
      ALLOCATE ( cmt(dimension%neigd,hybrid%maxlmindx,atoms%nat), stat=ok)
      IF(  ok .ne. 0 ) STOP 'gen_wavf: Failure allocation cmt'
      ALLOCATE ( cmthlp(dimension%neigd,hybrid%maxlmindx,atoms%nat), stat=ok)
      IF( ok .ne. 0) STOP 'gen_wavf: failure allocation cmthlp'
    
   

      DO ikpt0 = lower, upper

        acof = 0; bcof = 0; ccof = 0

        ! abcof calculates the wavefunction coefficients
        ! stored in acof,bcof,ccof
        lapw(ikpt0)%nmat=lapw(ikpt0)%nv(jsp)+atoms%nlotot
        CALL abcof(&
              input,atoms,hybrid%nbands(ikpt0),sym, cell, Kpts%bk(:,ikpt0), lapw(ikpt0), &
              hybrid%nbands(ikpt0),usdus,noco,jsp,hybdat%kveclo_eig(:,ikpt0),&
              oneD,acof(: hybrid%nbands(ikpt0),:,:),bcof(: hybrid%nbands(ikpt0),:,:),ccof(:,: hybrid%nbands(ikpt0),:,:),&
              zmat(ikpt0))
        

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
                hybrid,atoms,hybrid%nbands(ikpt0),&
                 sym,&
                cell,oneD,&
                acof(: hybrid%nbands(ikpt0),:,:),bcof(: hybrid%nbands(ikpt0),:,:),&
                ccof(:,: hybrid%nbands(ikpt0),:,:) )

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
        call write_cmt(cmt,ikpt0)
        call zhlp%alloc(zmat(1)%l_real,zmat(1)%nbasfcn,zmat(1)%nbands)
        
        IF (zhlp%l_real) THEN
           zhlp%data_r=zmat(ikpt0)%z_r
        ELSE
           zhlp%data_c=zmat(ikpt0)%z_c
        end IF
        call write_z(zhlp,ikpt0)
       

        ! generate wavefunctions coefficients at all k-points from
        ! irreducible k-points
        
        DO ikpt=1,kpts%nkptf
          IF ( kpts%bkp(ikpt) .eq. ikpt0 .and. ikpt0 .ne. ikpt ) THEN
            iop = kpts%bksym(ikpt)
              CALL waveftrafo_genwavf( cmthlp,zhlp%data_r,zhlp%data_c,&
     &                 cmt(:,:,:),zmat(1)%l_real,zmat(ikpt0)%z_r(:,:),zmat(ikpt0)%z_c(:,:),ikpt0,iop,atoms,&
     &                 hybrid,kpts,sym,&
     &                 jsp,dimension,hybrid%nbands(ikpt0),&
     &                 cell,lapw(ikpt0),lapw(ikpt),.true.)

              call write_cmt(cmthlp,ikpt)
              call write_z(zhlp,ikpt)
          END IF
        END DO  !ikpt
      END DO !ikpt0

      DEALLOCATE( acof,bcof,ccof )
      DEALLOCATE( cmt,cmthlp)

     
      END SUBROUTINE gen_wavf

      END MODULE m_gen_wavf
