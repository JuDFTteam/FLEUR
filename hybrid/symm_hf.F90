!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   This module generates the little group of k and the extended irr. !
!   BZ. Furthermore it calculates the irr. representation             !
!                                                                     !
!   P(R,T)\phi_n,k = \sum_{n'} rep_v(n',n) *\phi_n',k        !
!   where                                                             !
!         P  is an element of the little group of k                   !
!         n' runs over group of degenerat states belonging to n.      !
!                                                                     !
!                                             M.Betzinger (09/07)     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MODULE m_symm_hf

#define irreps .false.

      CONTAINS

      SUBROUTINE symm_hf(kpts,nkpti,nk,sym,&
     &                   dimension,ne_eig,eig_irr,nbands,&
     &                   atoms,hybrid,bas1,bas2,cell,&
     &                   lapw,jsp,&
     &                   gpt,&
     &                   lmaxcd,&
     &                   mpi,irank2,&
     &                   nsymop,psym,nkpt_EIBZ,n_q,parent,&
     &                   symop,degenerat,pointer_EIBZ,maxndb,nddb,&
     &                   nsest,indx_sest,rep_c ) 


      USE m_constants
      USE m_util   ,ONLY: modulo1,intgrf,intgrf_init
      USE m_olap   ,ONLY: wfolap,wfolap1,wfolap_init
      USE m_trafo  ,ONLY: waveftrafo_symm
    USE m_types
      IMPLICIT NONE

      TYPE(t_mpi),INTENT(IN)   :: mpi
      TYPE(t_dimension),INTENT(IN)   :: dimension
      TYPE(t_hybrid),INTENT(IN)   :: hybrid
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_cell),INTENT(IN)   :: cell
      TYPE(t_kpts),INTENT(IN)   :: kpts
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_lapw),INTENT(IN)   :: lapw

!     - scalars -
      INTEGER,INTENT(IN)              :: nkpti  ,nk,ne_eig,nbands
      INTEGER,INTENT(IN)              :: jsp
      INTEGER,INTENT(IN)              :: lmaxcd
      INTEGER,INTENT(IN)              :: irank2
      INTEGER,INTENT(OUT)             :: nkpt_EIBZ
      INTEGER,INTENT(OUT)             :: nsymop
      INTEGER,INTENT(OUT)             :: maxndb,nddb

!     - arrays -
      INTEGER,INTENT(IN)              :: gpt(3,lapw%nv(jsp))
      INTEGER,INTENT(OUT)             :: parent(kpts%nkpt)
      INTEGER,INTENT(OUT)             :: symop(kpts%nkpt)
      INTEGER,INTENT(INOUT)           :: degenerat(ne_eig)
      INTEGER,INTENT(OUT)             :: nsest(nbands), indx_sest(nbands,nbands)  
      INTEGER,ALLOCATABLE,INTENT(OUT) :: pointer_EIBZ(:)
      INTEGER,ALLOCATABLE,INTENT(OUT) :: psym(:),n_q(:)      

      REAL,INTENT(IN)                 :: eig_irr(dimension%neigd,nkpti)
      REAL,INTENT(IN)               :: bas1(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype),&
     &                                 bas2(atoms%jmtd,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype)
      COMPLEX,ALLOCATABLE,INTENT(OUT) :: rep_c(:,:,:,:,:)

!     - local scalars -
      INTEGER                         :: ikpt,ikpt1,iop,isym,iisym,m
      INTEGER                         :: itype,ieq,iatom,ratom
      INTEGER                         :: iband,iband1,iband2,iatom0
      INTEGER                         :: i,j,ic,ic1,ic2
      INTEGER                         :: irecl_cmt,irecl_z
      INTEGER                         :: ok
      INTEGER                         :: l,lm
      INTEGER                         :: n1,n2,nn
      INTEGER                         :: ndb,ndb1,ndb2
      INTEGER                         :: nrkpt

      REAL                            :: tolerance,pi

      COMPLEX                         :: cdum
      COMPLEX , PARAMETER             :: img = (0d0,1d0)

!     - local arrays -
      INTEGER                         :: rrot(3,3,sym%nsym)
      INTEGER                         :: neqvkpt(kpts%nkpt)
      INTEGER                         :: list(kpts%nkpt)
      INTEGER,ALLOCATABLE             :: help(:)
      
      REAL                            :: rotkpt(3),g(3)
      REAL   ,ALLOCATABLE             :: olapmt(:,:,:,:)
      REAL   ,ALLOCATABLE             :: gridf(:,:)

      COMPLEX                         :: cmt(dimension%neigd,hybrid%maxlmindx,atoms%nat)
      COMPLEX                         :: carr1(hybrid%maxlmindx,atoms%nat)
      COMPLEX,ALLOCATABLE             :: carr(:),wavefolap(:,:)
      COMPLEX,ALLOCATABLE             :: cmthlp(:,:,:)
      COMPLEX,ALLOCATABLE             :: cpwhlp(:,:)
      COMPLEX,ALLOCATABLE             :: trace(:,:)

#if ( defined(CPP_INVERSION) )
      REAL                            :: carr2(lapw%nv(jsp))
      REAL                            :: z(dimension%nbasfcn,dimension%neigd)
      REAL,ALLOCATABLE                :: olappw(:,:)
#else
      COMPLEX                         :: carr2(lapw%nv(jsp))
      COMPLEX                         :: z(dimension%nbasfcn,dimension%neigd)
      COMPLEX,ALLOCATABLE             :: olappw(:,:)
#endif 
      COMPLEX,ALLOCATABLE             :: rep_d(:,:,:)
      LOGICAL,ALLOCATABLE             :: symequivalent(:,:)

      IF ( irank2 == 0 ) THEN
        WRITE(6,'(A)') new_line('n') // new_line('n') // '### subroutine: symm ###'
      END IF

      ! calculate rotations in reciprocal space
      DO i = 1,sym%nsym
        IF( i .le. sym%nop ) THEN
          rrot(:,:,i) = transpose(sym%mrot(:,:,sym%invtab(i)))
        ELSE
          rrot(:,:,i) = -rrot(:,:,i-sym%nop)
        END IF
      END DO


      ! determine little group of k., i.e. those symmetry operations
      ! which keep bk(:,nk) invariant
      ! nsymop :: number of such symmetry-operations
      ! psym   :: points to the symmetry-operation

      ic = 0
      ALLOCATE(psym(sym%nsym))
      
      DO iop=1,sym%nsym

        rotkpt = matmul( rrot(:,:,iop), kpts%bk(:,nk) )

        !transfer rotkpt into BZ
        rotkpt = modulo1(rotkpt,kpts%nkpt3)

        !check if rotkpt is identical to bk(:,nk)
        IF( maxval( abs( rotkpt - kpts%bk(:,nk) ) ) .le. 1E-07) THEN
          ic = ic + 1
          psym(ic) = iop
        END IF
      END DO
      nsymop = ic
     

      IF ( irank2 == 0 ) THEN
        WRITE(6,'(A,i3)') ' nk',nk
        WRITE(6,'(A,3f10.5)') ' kpts%bk(:,nk):',kpts%bk(:,nk)
        WRITE(6,'(A,i3)') ' Number of elements in the little group:',nsymop
      END IF

      ! reallocate psym
      ALLOCATE(help(ic))
      help = psym(1:ic)
      DEALLOCATE(psym)
      ALLOCATE(psym(ic))
      psym =help
      DEALLOCATE(help)

      ! determine extented irreducible BZ of k ( EIBZ(k) ), i.e.
      ! those k-points, which can generate the whole BZ by
      ! applying the symmetry operations of the little group of k

      neqvkpt = 0

      DO i=1,kpts%nkpt
        list(i) = i-1
      END DO

      symop = 0
      DO ikpt=2,kpts%nkpt
        DO iop=1,nsymop

          rotkpt = matmul( rrot(:,:,psym(iop)), kpts%bk(:,ikpt) )

          !transfer rotkpt into BZ
          rotkpt = modulo1(rotkpt,kpts%nkpt3)

          !determine number of rotkpt
          nrkpt = 0
          DO ikpt1=1,kpts%nkpt
            IF ( maxval( abs( rotkpt - kpts%bk(:,ikpt1) ) ) <= 1E-06 ) THEN
              nrkpt = ikpt1
              EXIT
            END IF
          END DO
          IF( nrkpt .eq.0 ) STOP 'symm: Difference vector not found !'

          IF( list(nrkpt) .ne. 0) THEN
            list(nrkpt)   = 0
            neqvkpt(ikpt) = neqvkpt(ikpt) + 1
            parent(nrkpt) = ikpt
            symop(nrkpt)  = psym(iop)
          END IF
          IF ( all(list .eq. 0) ) EXIT

        END DO
      END DO

      ! for the Gamma-point holds:
      parent(1)  = 1
      neqvkpt(1) = 1

#ifdef CPP_DEBUG
      IF( sum(neqvkpt(:)) .ne. kpts%nkpt) THEN
        IF ( mpi%irank == 0 ) WRITE(6,'(A,i3)') ' Check EIBZ(k) by summing&&
     & up equivalent k-points:',sum(neqvkpt(:))
        STOP 'symm: neqvkpt not identical to nkpt'
      END IF
#endif

      ! determine number of members in the EIBZ(k)
      ic = 0
      DO ikpt=1,kpts%nkpt
        IF(parent(ikpt) .eq. ikpt) ic = ic + 1
      END DO
      nkpt_EIBZ = ic

      ALLOCATE( pointer_EIBZ(nkpt_EIBZ) )
      ic = 0
      DO ikpt=1,kpts%nkpt
        IF(parent(ikpt) .eq. ikpt) THEN
          ic = ic + 1
          pointer_EIBZ(ic) = ikpt
        END IF
      END DO

      IF ( irank2 == 0 ) THEN
        WRITE(6,'(A,i5)') ' Number of k-points in the EIBZ',nkpt_EIBZ
      END IF


      ! determine the factor n_q, that means the number of symmetrie operations of the little group of bk(:,nk)
      ! which keep q (in EIBZ) invariant

      ALLOCATE( n_q(nkpt_EIBZ) )

      ic  = 0
      n_q = 0
      DO ikpt = 1,kpts%nkpt
        IF ( parent(ikpt) .eq. ikpt ) THEN
          ic = ic + 1
          DO iop=1,nsymop
            isym = psym(iop)
            rotkpt = matmul( rrot(:,:,isym), kpts%bk(:,ikpt) )

            !transfer rotkpt into BZ
            rotkpt = modulo1(rotkpt,kpts%nkpt3)

            !check if rotkpt is identical to bk(:,ikpt)
            IF( maxval( abs( rotkpt - kpts%bk(:,ikpt) ) ) .le. 1E-06) THEN
              n_q(ic) = n_q(ic) + 1 
            END IF
          END DO
        END IF
      END DO
      IF( ic .ne. nkpt_EIBZ ) STOP 'symm: failure EIBZ'


      ! calculate degeneracy:
      ! degenerat(i) = 1 state i  is not degenerat,
      ! degenerat(i) = j state i has j-1 degenerat states at {i+1,...,i+j-1}
      ! degenerat(i) = 0 state i is degenerat

      tolerance = 1E-07 !0.00001

      degenerat = 1
      IF ( irank2 == 0 ) THEN
        WRITE(6,'(A,f10.8)') ' Tolerance for determining degenerate states=', tolerance
     END IF

      DO i=1,nbands
        DO j=i+1,nbands
          IF( abs(eig_irr(i,nk)-eig_irr(j,nk)) .le. tolerance) THEN
            degenerat(i) = degenerat(i) + 1
          END IF
        END DO
      END DO

      DO i=1,ne_eig
        IF ( degenerat(i) .ne. 1 .or. degenerat(i) .ne. 0 ) THEN
          degenerat(i+1:i+degenerat(i)-1) = 0
        END IF
      END DO


      ! maximal number of degenerate bands -> maxndb
      maxndb = maxval(degenerat)

      ! number of different degenerate bands/states
      nddb   = count( degenerat .ge. 1)


      IF ( irank2 == 0 ) THEN
        WRITE(6,*) ' Degenerate states:'
        DO iband = 1,nbands/5+1
          WRITE(6,'(5i5)')degenerat(iband*5-4:min(iband*5,nbands))
        END DO
      END IF

      IF( irreps ) THEN
        ! calculate representation, i.e. the action of an element of
        ! the little group of k on \phi_n,k:
        ! P(R,T)\phi_n,k = \sum_{n'\in degenerat(n)} rep_v(n',n) *\phi_n',k
      
        ! read in cmt and z at current k-point (nk)
        irecl_cmt = dimension%neigd*hybrid%maxlmindx*atoms%nat*16
        OPEN(unit=777,file='cmt',form='unformatted',access='direct',&
     &       recl=irecl_cmt)
        READ(777,rec=nk) cmt
        CLOSE(777)

#ifdef CPP_INVERSION
        irecl_z   =  dimension%nbasfcn*dimension%neigd*8
#else
        irecl_z   =  dimension%nbasfcn*dimension%neigd*16
#endif
        OPEN(unit=778,file='z',form='unformatted',access='direct',&
     &      recl=irecl_z)
        READ(778,rec=nk) z
        CLOSE(778)

        ALLOCATE(rep_d(maxndb,nddb,nsymop),stat=ok )
        IF( ok .ne. 0) STOP 'symm: failure allocation rep_v'


        ALLOCATE( olappw(lapw%nv(jsp),lapw%nv(jsp)),&
     &            olapmt(hybrid%maxindx,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype),stat=ok)
        IF( ok .ne. 0) STOP 'symm: failure allocation olappw/olapmt'

        olappw = 0
        olapmt = 0
        CALL wfolap_init (olappw,olapmt,gpt,&
     &                   atoms,hybrid,&
     &                    cell,&
     &                    bas1,bas2)


#ifdef CPP_DEBUG1
        ! check orthogonality of wavefunctions
        OPEN(677,file='ortho')
        WRITE(677,*) 'iteration',it,'k-point',nk
        DO iband=1,nbands
          IF(degenerat(iband) .ge. 1) THEN
            ndb = degenerat(iband)
            DO i= 1,nbands!iband,iband+ndb-1
              WRITE(677,*) iband,i,&
     &          wfolap1( cmt(iband,:,:),z(:lapw%nv(jsp),iband),&
     &           cmt(i,:,:),z(:lapw%nv(jsp),i),&
     &          lapw%nv(jsp),lapw%nv(jsp),olappw,olapmt,&
     &          atoms%ntype,atoms%neq,atoms%nat,atoms%lmax,atoms%lmaxd,hybrid%nindx,hybrid%maxindx,hybrid%maxlmindx)


            END DO
          
          END IF
        END DO
        CLOSE(677)
#endif

        ALLOCATE( cmthlp(hybrid%maxlmindx,atoms%nat,maxndb), cpwhlp(lapw%nv(jsp),maxndb),&
     &            stat= ok )
        IF( ok .ne. 0 ) STOP 'symm: failure allocation cmthlp/cpwhlp' 

        DO isym=1,nsymop 
          iop= psym(isym) 

#ifdef CPP_DEBUG
          rotkpt   = matmul(rrot(:,:,iop),kpts%bk(:,nk))
          rotkpt   = modulo1(rotkpt,kpts%nkpt3)
          IF( maxval( abs( rotkpt - kpts%bk(:,nk) ) ) .gt. 1e-08) THEN
            STOP 'symm: error in psym'
          END IF
#endif
          ic = 0
          DO i=1,nbands
            ndb = degenerat(i)
            IF ( ndb .ge. 1 ) THEN
              ic     = ic + 1
              cmthlp = 0
              cpwhlp = 0
              STOP "WAVETRAVO_SYM real/complex"
              CALL waveftrafo_symm( &
     &                      cmthlp(:,:,:ndb),cpwhlp(:,:ndb),cmt,.false.,(/0.1/),z,&
     &                      i,ndb,nk,iop,atoms,&
     &                      hybrid,kpts,&
     &                      sym,&
     &                      jsp,dimension,&
     &                      cell,gpt,lapw )

            
              DO iband = 1,ndb
                carr1 = cmt(iband+i-1,:,:)
                carr2 = z(:lapw%nv(jsp),iband+i-1)
                rep_d(iband,ic,isym) &
     &        = wfolap( carr1,carr2,cmthlp(:,:,iband),cpwhlp(:,iband),&
     &                  lapw%nv(jsp),lapw%nv(jsp),olappw,olapmt,atoms,hybrid)
              END DO
            
            END IF
          END DO

        END DO

        DEALLOCATE( cmthlp,cpwhlp )


        ! calculate trace of irrecudible representation
        ALLOCATE( trace(sym%nsym,nddb), stat=ok )
        IF( ok .ne. 0) STOP 'symm: failure allocation trace'

        ic    = 0
        trace = 0
        DO iband = 1,nbands
          ndb = degenerat(iband)
          IF( ndb .ge. 1 ) THEN
            ic = ic + 1
            !calculate trace
            DO iop = 1,nsymop
              isym = psym(iop)
              DO i = 1,ndb
                trace(isym,ic) = trace(isym,ic) + rep_d(i,ic,iop)
              END DO
            END DO
          END IF
        END DO


        ! determine symmetry equivalent bands/irreducible representations by comparing the trace

        ALLOCATE( symequivalent(nddb,nddb), stat=ok )
        IF( ok .ne. 0) STOP 'symm: failure allocation symequivalent'

        ic1 = 0
        symequivalent = .false.
        DO iband1 = 1,nbands
          ndb1 = degenerat(iband1)
          IF( ndb1 .ge. 1 ) THEN
            ic1 = ic1 + 1
            ic2 = 0
            DO iband2 = 1,nbands
              ndb2 = degenerat(iband2)
              IF( ndb2 .ge. 1 ) THEN
                ic2 = ic2 + 1
                IF( ndb2 .eq. ndb1 ) THEN
                  ! note that this criterium is only valid for pure spatial rotations
                  ! if one combines spatial rotations with time reversal symmetry there
                  ! is no unique criteria to identify symequivalent state
                  ! however, also in the latter case the trace of the spatial rotations 
                  ! for two symmetry equivalent states must be equivalent
                  IF( all(abs(trace(:sym%nop,ic1)-trace(:sym%nop,ic2)) <= 1E-8))&
     &            THEN
                    symequivalent(ic2,ic1) = .true.
                  END IF
                END IF
              END IF
            END DO
          END IF
        END DO
      
        
      ELSE
        ! read in cmt and z at current k-point (nk)
        irecl_cmt = dimension%neigd*hybrid%maxlmindx*atoms%nat*16
        OPEN(unit=777,file='cmt',form='unformatted',access='direct',&
     &       recl=irecl_cmt)
        READ(777,rec=nk) cmt
        CLOSE(777)
        
        CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,gridf)
        
        IF( allocated(olapmt)) deallocate(olapmt)
        ALLOCATE( olapmt(hybrid%maxindx,hybrid%maxindx,0:atoms%lmaxd,atoms%ntype),stat=ok)
        IF( ok .ne. 0) STOP 'symm: failure allocation olapmt'
        olapmt = 0

        DO itype = 1,atoms%ntype
          DO l = 0,atoms%lmax(itype)
            nn = hybrid%nindx(l,itype)
            DO n2 = 1,nn
              DO n1 = 1,nn
                olapmt(n1,n2,l,itype) = intgrf ( &
     &                        bas1(:,n1,l,itype)*bas1(:,n2,l,itype)&
     &                       +bas2(:,n1,l,itype)*bas2(:,n2,l,itype),&
     &                        atoms%jri,atoms%jmtd,atoms%rmsh,atoms%dx,atoms%ntype,itype,gridf)
              END DO
            END DO
          END DO
        END DO
      
        ALLOCATE( wavefolap(nbands,nbands),carr(hybrid%maxindx),stat=ok )
        IF( ok .ne. 0) STOP 'symm: failure allocation wfolap/maxindx'
        wavefolap = 0
        
        iatom  = 0
        DO itype = 1,atoms%ntype
          DO ieq = 1,atoms%neq(itype)
            iatom = iatom + 1
            lm = 0
            DO l = 0,atoms%lmax(itype)
              DO M = -l,l
                nn     = hybrid%nindx(l,itype)
                DO iband1 = 1,nbands
                  carr(:nn) = matmul( olapmt(:nn,:nn,l,itype),&
     &                                cmt(iband1,lm+1:lm+nn,iatom)) 
                  DO iband2 = 1,iband1 
                    wavefolap(iband2,iband1)&
     &            = wavefolap(iband2,iband1)&
     &            + dot_product(cmt(iband2,lm+1:lm+nn,iatom),carr(:nn))
                  END DO
                END DO
                lm     = lm + nn
              END DO
            END DO
          END DO
        END DO
        
        DO iband1 = 1,nbands
          DO iband2 = 1,iband1
            wavefolap(iband1,iband2) = conjg(wavefolap(iband2,iband1))
          END DO
        END DO
        
        ALLOCATE( symequivalent(nddb,nddb), stat=ok )
        IF( ok .ne. 0) STOP 'symm: failure allocation symequivalent'
        symequivalent = .false.
        ic1 = 0
        DO iband1 = 1,nbands
          ndb1 = degenerat(iband1)
          IF( ndb1 .eq. 0 ) CYCLE
          ic1 = ic1 + 1
          ic2 = 0
          DO iband2 = 1,nbands
            ndb2 = degenerat(iband2)
            IF( ndb2 .eq. 0 ) CYCLE
            ic2 = ic2 + 1
            IF( any( abs(wavefolap(iband1:iband1+ndb1-1,&
     &                             iband2:iband2+ndb2-1)) > 1E-9 ) )THEN
!     &          .and. ndb1 .eq. ndb2 ) THEN
              symequivalent(ic2,ic1) = .true.
            END IF
          END DO
        END DO
      END IF
      
      !
      ! generate index field which contain the band combinations (n1,n2),
      ! which are non zero 
      !   
      
      ic1       = 0
      indx_sest = 0
      nsest     = 0
      DO iband1 = 1,nbands
        ndb1 = degenerat(iband1)
        IF( ndb1 .ge. 1 ) ic1 = ic1 + 1
        i = 0
        DO WHILE ( degenerat(iband1-i) == 0 )
          i = i+1
        END DO
        ndb1 = degenerat(iband1-i)
        ic2 = 0
        DO iband2 = 1,nbands
          ndb2 = degenerat(iband2)
          IF( ndb2 .ge. 1 ) ic2 = ic2 + 1
          i = 0
          DO WHILE ( degenerat(iband2-i) == 0 )
            i = i+1
          END DO
          ndb2 = degenerat(iband2-i)
          ! only upper triangular part
          IF( symequivalent(ic2,ic1) .and. iband2 .le. iband1 ) THEN
!            IF( ndb1 .ne. ndb2 ) STOP 'symm_hf: failure symequivalent'
            nsest(iband1)                   = nsest(iband1) + 1
            indx_sest(nsest(iband1),iband1) = iband2
          END IF
        END DO
      END DO

      !
      ! calculate representations for core states
      ! (note that for a core state, these are proportional to the Wigner D matrices)
      !
      ! Definition of the Wigner rotation matrices
      !
      !                     -1                l
      ! P(R) Y  (r)  = Y  (R  r)  =  sum     D    (R)  Y   (r)
      !       lm        lm              m'    m m'      lm'
      !

      pi = pimach()

      IF( lmaxcd .gt. atoms%lmaxd ) STOP &
     & 'symm_hf: The very impropable case that lmaxcd > atoms%lmaxd occurs'

      ALLOCATE(rep_c(-lmaxcd:lmaxcd,-lmaxcd:lmaxcd,0:lmaxcd,nsymop,atoms%nat)&
     &        , stat=ok )
      IF( ok .ne. 0) STOP 'symm_hf: failure allocation rep_c'

      iatom  = 0
      iatom0 = 0
      DO itype = 1,atoms%ntype
        DO ieq = 1,atoms%neq(itype)
          iatom = iatom + 1
          DO iop = 1,nsymop
            isym = psym(iop)
            IF( isym .le. sym%nop) THEN
              iisym = isym
            ELSE
              iisym = isym - sym%nop
            END IF

            ratom  = hybrid%map(iatom,isym)
            rotkpt = matmul( rrot(:,:,isym), kpts%bk(:,nk) )
            g      = nint(rotkpt - kpts%bk(:,nk))

            cdum   = exp(-2*pi*img*dot_product(rotkpt,sym%tau(:,iisym)))* &
     &               exp( 2*pi*img*dot_product(g,atoms%taual(:,ratom)))

            rep_c(:,:,:,iop,iatom) = &
     &         sym%d_wgn(-lmaxcd:lmaxcd,-lmaxcd:lmaxcd,0:lmaxcd,isym) * cdum
          END DO
        END DO
        iatom0 = iatom0 + atoms%neq(itype)
      END DO

      END SUBROUTINE symm_hf

      INTEGER FUNCTION symm_hf_nkpt_EIBZ(kpts,nk,sym)

      USE m_util, ONLY: modulo1
    USE m_types
      IMPLICIT NONE
      TYPE(t_sym),INTENT(IN)   :: sym
      TYPE(t_kpts),INTENT(IN)   :: kpts

!     - scalar input -
      INTEGER, INTENT(IN)   :: nk
!     - array input -

!     - local scalars -
      INTEGER               ::  isym,ic,iop,ikpt,ikpt1
      INTEGER               ::  nsymop,nrkpt
!     - local arrays -
      INTEGER               ::  rrot(3,3,sym%nsym)
      INTEGER               ::  neqvkpt(kpts%nkpt),list(kpts%nkpt),parent(kpts%nkpt),&
     &                          symop(kpts%nkpt)
      INTEGER, ALLOCATABLE  ::  psym(:)!,help(:)
      REAL                  ::  rotkpt(3)

      ! calculate rotations in reciprocal space
      DO isym = 1,sym%nsym
        IF( isym .le. sym%nop ) THEN
          rrot(:,:,isym) = transpose(sym%mrot(:,:,sym%invtab(isym)))
        ELSE
          rrot(:,:,isym) = -rrot(:,:,isym-sym%nop)
        END IF
      END DO

      ! determine little group of k., i.e. those symmetry operations
      ! which keep bk(:,nk,nw) invariant
      ! nsymop :: number of such symmetry-operations
      ! psym   :: points to the symmetry-operation

      ic = 0
      ALLOCATE(psym(sym%nsym))

      DO iop=1,sym%nsym
        rotkpt = matmul( rrot(:,:,iop), kpts%bk(:,nk) )

        !transfer rotkpt into BZ
        rotkpt = modulo1(rotkpt,kpts%nkpt3)

        !check if rotkpt is identical to bk(:,nk)
        IF( maxval( abs( rotkpt - kpts%bk(:,nk) ) ) .le. 1E-07) THEN
          ic = ic + 1
          psym(ic) = iop
        END IF
      END DO
      nsymop = ic

      ! reallocate psym
!       ALLOCATE(help(ic))
!       help = psym(1:ic)
!       DEALLOCATE(psym)
!       ALLOCATE(psym(ic))
!       psym = help
!       DEALLOCATE(help)

      ! determine extented irreducible BZ of k ( EIBZ(k) ), i.e.
      ! those k-points, which can generate the whole BZ by
      ! applying the symmetry operations of the little group of k

      neqvkpt = 0

!       list = (/ (ikpt-1, ikpt=1,nkpt) /)
      DO ikpt=1,kpts%nkpt
        list(ikpt) = ikpt-1
      END DO

      DO ikpt=2,kpts%nkpt
        DO iop=1,nsymop

          rotkpt = matmul( rrot(:,:,psym(iop)), kpts%bk(:,ikpt) )

          !transfer rotkpt into BZ
          rotkpt = modulo1(rotkpt,kpts%nkpt3)

          !determine number of rotkpt
          nrkpt = 0
          DO ikpt1=1,kpts%nkpt
            IF ( maxval( abs( rotkpt - kpts%bk(:,ikpt1) ) ) .le. 1E-06 ) THEN
              nrkpt = ikpt1
              EXIT
            END IF
          END DO
          IF( nrkpt .eq.0 ) STOP 'symm: Difference vector not found !'

          IF( list(nrkpt) .ne. 0) THEN
            list(nrkpt)   = 0
            neqvkpt(ikpt) = neqvkpt(ikpt) + 1
            parent(nrkpt) = ikpt
            symop(nrkpt)  = psym(iop)
          END IF
          IF ( all(list .eq. 0) ) EXIT

        END DO
      END DO

      ! for the Gamma-point holds:
      parent(1)  = 1
      neqvkpt(1) = 1

      ! determine number of members in the EIBZ(k)
      ic = 0
      DO ikpt=1,kpts%nkpt
        IF(parent(ikpt) .eq. ikpt) ic = ic + 1
      END DO
      symm_hf_nkpt_EIBZ = ic

      END FUNCTION symm_hf_nkpt_EIBZ

      END MODULE m_symm_hf
