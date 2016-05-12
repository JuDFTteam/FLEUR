      MODULE m_atomsym
      use m_juDFT
      CONTAINS
      SUBROUTINE atom_sym(
     >                    dispfh,outfh,errfh,dispfn,natmax,
     X                    natin,atomid,atompos,
     X                    ngen,mmrot,ttr,
     >                    cartesian,symor,as,bs,nop48,
     <                    ntype,nat,nops,mrot,tau,
     <                    neq,ntyrep,zatom,natype,natrep,natmap,pos)

!***********************************************************************
!        reads in atomic positions and space group information
!        for the case when natin < 0 or when cal_symm = .false.
!
!--->    atomic positions:
!--->
!--->    atomic positions input can be either in scaled cartesian
!--->    or lattice vector units, as determined by logical cartesian.
!--->    (for supercells, sometimes more natural to input positions
!--->    in scaled cartesian.) only represenative atoms need be
!--->    given, but will handle all (generates all of them and then
!--->    throws out duplicates).
!--->
!--->    for each atom, give identification (atomic) number and
!--->    position, e.g.,
!                 26    0.5  0.5  0.25
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!--->    space group symmetry information:
!--->
!--->    to input symmetry, give the generators (except identity)
!--->    of the group as a matrix plus a non-primitive translation,
!--->    all in lattice coordinates.
!--->
!--->    for each generator, there will be 4 lines of input
!--->    corresponding to the matrix as one would write it:
!--->
!--->        m(1,1)  m(1,2) m(1,3)    ! integer defining
!--->        m(2,1)  m(2,2) m(2,3)    ! rotation matrix for
!--->        m(3,1)  m(3,2) m(3,3)    ! column vectors
!--->        tau(1)  tau(2) tau(3)    ! non-primitive translation
!                                                               mw 12-99
!***********************************************************************

      USE m_closure, ONLY : closure
      IMPLICIT NONE

!==> Arguments

      INTEGER, INTENT (IN)    :: dispfh,outfh,errfh
      LOGICAL, INTENT (IN)    :: cartesian         ! how atomic positions are given
      LOGICAL, INTENT (IN)    :: symor             ! whether to reduce to symmorphic subgroup
      INTEGER, INTENT (INOUT) :: natin             ! formerly 'ntype0'
      INTEGER, INTENT (IN)    :: ngen              ! Number of generators
      INTEGER, INTENT (IN)    :: nop48, natmax     ! dimensioning
      INTEGER, INTENT (INOUT) :: mmrot(3,3,nop48)
      REAL,    INTENT (INOUT) :: ttr(3,nop48)
      REAL,    INTENT (IN)    :: as(3,3),bs(3,3)
      REAL,    INTENT (INOUT) :: atomid(natmax)    ! formerly 'zat0'
      REAL,    INTENT (INOUT) :: atompos(3,natmax) ! formerly 'postype0'
      INTEGER, INTENT (OUT)   :: ntype,nat,nops
      CHARACTER(len=4), INTENT (IN) :: dispfn
!--> actually, intent out:
      INTEGER, ALLOCATABLE :: neq(:), ntyrep(:)    ! here, these variables are allocated with
      REAL,    ALLOCATABLE :: zatom(:)             ! the dimension 'ntype'
      INTEGER, ALLOCATABLE :: natype(:),natrep(:),natmap(:)  ! or  'nat'
      REAL,    ALLOCATABLE :: pos(:,:)                       ! or  '3,nat'
      INTEGER, ALLOCATABLE :: mrot(:,:,:)                    ! or  '3,3,nop'
      REAL,    ALLOCATABLE :: tau(:,:)                       ! or  '3,nop'

!---->Automatic arrays
      INTEGER nneq(natmax),icount(nop48*natmax),imap(natmax)
      INTEGER ity(nop48*natmax),index_op(nop48)
      REAL    tpos(3,nop48*natmax),atomid2(natmax)
      REAL    tr(3),tt(3),disp(3,natmax)

      INTEGER mp(3,3),mtmp(3,3)
      REAL    ttau(3),orth(3,3)

      INTEGER mmrot2(3,3,ngen)
      REAL    ttr2(3,ngen)

      LOGICAL l_exist,lclose,l_inipos
      INTEGER n,na,ng,ncyl,nc,no,nop0,nn,nt,i,j,mops
      INTEGER ios,istep0
      REAL    eps7
      CHARACTER(len=30) :: filen
      REAL,    ALLOCATABLE :: inipos(:,:)

      eps7= 1.0e-7 ; istep0 = 0
!
!---> take atomic positions and shift to (-1/2,1/2] in lattice coords.
!
      natin = abs(natin)
      DO n=1,natin
         IF (cartesian) THEN  ! convert to lattice coords. if necessary
            atompos(:,n) = matmul( bs, atompos(:,n) )
         ENDIF
         atompos(:,n) = atompos(:,n) - anint( atompos(:,n) - eps7 )
      ENDDO

!--->  store the positions (in lattice coord.s) given in the input file
!      ALLOCATE ( inipos(3,natin) )
!      DO n = 1, natin
!         inipos(:,n) = atompos(:,n)
!      ENDDO
!      l_inipos = .false.

!---> check whether a displacement file exists (dispfh)
!---> these displacement should be in the same units as input file,
!---> i.e., lattice  or scaled cartesian

      l_exist = .false.
      INQUIRE (file=dispfn,exist=l_exist)
      IF (l_exist) THEN
        OPEN (dispfh,file=dispfn,status='old',form='formatted')

        DO n = 1, natin
          READ (dispfh,*,iostat=ios) tr(1:3)
          IF (ios .NE. 0) EXIT
          disp(1:3,n) = tr(1:3)
          IF (cartesian) THEN    ! convert to lattice coords. if necessary
            tr = matmul( bs, tr )
          ENDIF
          atompos(:,n) = atompos(:,n) + tr(:)
          atompos(:,n) = atompos(:,n) - anint( atompos(:,n)- eps7 )
        ENDDO
        CLOSE (dispfh)
        IF ( ios==0 ) THEN
          WRITE(filen,'(a,"_",i2.2)') trim(adjustl(dispfn)),istep0+1
          OPEN(dispfh,file=filen,status='unknown',form='formatted')
          DO n=1,natin
            WRITE (dispfh,'(3f16.8)') disp(1:3,n)
          END DO
          CLOSE (dispfh)
        ENDIF
      ENDIF ! l_exist

      IF ( mmrot(1,1,1)==0 ) THEN  ! &gen was used

!--->   save generators
        IF (cartesian) THEN    ! convert to lattice coords. if necessary
          DO ng = 2, ngen+1
            mmrot2(:,:,1) = matmul( bs, mmrot(:,:,ng) )
            mmrot(:,:,ng) = matmul( mmrot2(:,:,1), as )
            ttr2(:,1) = matmul( bs, ttr(:,ng) )
            ttr(:,ng) = ttr2(:,1)
          ENDDO
        ENDIF
        mmrot2(:,:,1:ngen) = mmrot(:,:,2:ngen+1) 
        ttr2(:,1:ngen)     = ttr(:,2:ngen+1) 

!--->   set-up identity
        nops = 1
        mmrot(:,:,1) = reshape( (/ 1,0,0, 0,1,0, 0,0,1 /), (/ 3,3 /) )
        ttr(:,1) = 0.00 

        DO ng = 1, ngen

          n = nops + 1

!--->    next generator
          mmrot(:,:,n) = mmrot2(:,:,ng) 
          ttr(:,n)     = ttr2(:,ng) 

!--->    checking for order of cyclic (sub)group
          ncyl=1
          mtmp = mmrot(:,:,n)
          ttau = ttr(:,n)

          DO                    ! multiply on left until identity
            mp = matmul( mmrot(:,:,n) , mtmp )
            
            IF (mp(1,1)==1 .and. mp(2,2)==1 .and. mp(3,3)==1 .and.
     &          mp(1,2)==0 .and. mp(2,1)==0 .and. mp(3,2)==0 .and.
     &          mp(3,2)==0 .and. mp(3,1)==0 .and. mp(1,3)==0)  EXIT

!--->    new cyclic operation
            ncyl = ncyl+1
            mmrot(:,:,nops+ncyl) = mp
            ttr(:,nops+ncyl) = matmul( mmrot(:,:,n) , ttau ) + ttr(:,n)

            mtmp = mp
            ttau = ttr(:,nops+ncyl)
          ENDDO

!--->    now multiply these new operations into previous ones
!--->    (not the identity since that already included above)
          nop0=nops
          nops=nops+ncyl

          DO nc=1,ncyl
            DO nn=2,nop0 ! don't multiply by identity, nn=1
              nops = nops+1
              IF (nops.gt.48) then
                WRITE (6,'(" Error in generators: nops>48")')
                 CALL juDFT_error("atom_sym: nops>48",
     &                    calledby="atom_sym")
              ENDIF
              mmrot(:,:,nops) =
     &                matmul( mmrot(:,:,nop0+nc) , mmrot(:,:,nn) )
              ttr(:,nops) = ttr(:,nop0+nc) +
     &                matmul( mmrot(:,:,nop0+nc),ttr(:,nn) )
            ENDDO
          ENDDO

        ENDDO ! end loop on generators

      ELSE ! &sym was used

        nops = ngen + 1
        
      ENDIF

!---> rewrite all the non-primitive translations so in (-1/2,1/2]
      ttr(:,1:nops) = ttr(:,1:nops) - anint( ttr(:,1:nops)- eps7 )

!--->  allocate arrays for space group information (mod_spgsym)
!      if( nopd < nops )then
!        nopd = nops
!      endif
      ALLOCATE ( mrot(3,3,nops), tau(3,nops) )

      DO n = 1, nops
         mrot(:,:,n) = mmrot(:,:,n)
         tau(:,n) = ttr(:,n)
         index_op(n) = n
      ENDDO

!--->    check that the group is closed, etc.
      WRITE(6,*) 'nops=',nops

      mops = nops   ! different values appear in the call in spg_gen
      CALL closure(
     >             mops,mrot,tau,nops,index_op,
     <             lclose)

      IF ( .not. lclose ) then
         WRITE (6,'(/," Group did not close; check input")')
          CALL juDFT_error("atom_sym: not closed",calledby="atom_sym")
      ENDIF

      IF (symor) THEN
!--->   reduce symmetry to the largest symmorphic subgroup
        j = 1
        DO i = 1, nops
          IF ( all ( abs( tau(:,i) ) < eps7 ) ) THEN
            IF ( j<i ) then
              mrot(:,:,j) = mrot(:,:,i)
            ENDIF
            j = j + 1
          ENDIF
        END DO
        tau = 0.00
        IF ( nops > j-1 ) THEN
          WRITE (outfh,*) 'System has non-symmorphic symmetry with',
     &                   nops, 'operations.'
          nops = j - 1
          WRITE (outfh,*)'Symmetry reduced to symmorphic subgroup with',
     &   nops, 'operations.'
        ENDIF

      ENDIF ! symor

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!---> generate other symmetry equivalent positions from representatives
      ntype = 0
      nat = 0
      icount(:) = 0
      imap(:) = 0

      repres_atoms: DO nt = 1, natin  ! loop over input atoms

!--->    check whether this atom exists already
         DO n = 1, nat
            tt = ( atompos(:,nt) - tpos(:,n) )
     &           - anint( atompos(:,nt) - tpos(:,n) )
            IF ( all( abs(tt) < eps7 ) ) THEN
               icount(n) = icount(n) + 1
               imap(nt) = n
               IF ( abs( atomid(nt)-atomid(ity(n)) ) < eps7 )
     &             CYCLE repres_atoms
               CALL juDFT_error("ERROR! mismatch between atoms."
     +              ,calledby ="atom_sym")
            ENDIF
         ENDDO

!--->    add in the representative
         ntype = ntype+1
         nneq(ntype) = 1
         tpos(:,nat+1) = atompos(:,nt)
         ity(nat+1) = nt
         icount(nat+1) = icount(nat+1) + 1
         imap(nt) = -(nat+1)

!--->    loop over operations
         opts: DO no = 2, nops
            tr = matmul( mrot(:,:,no) , atompos(:,nt) ) + tau(:,no)
            tr = tr - anint( tr - eps7 )
!--->       check whether this is a new atom
            DO n = 1, nneq(ntype)
               tt = ( tr-tpos(:,nat+n) ) - anint( tr-tpos(:,nat+n) )
               IF ( all( abs(tt) < eps7 ) ) THEN

                  nn = ity(nat+n)
                  IF ( abs( atomid(nt)-atomid(nn) ) < eps7 ) CYCLE opts
                  WRITE (6,'(" Mismatch between atoms and",
     &                       " symmetry input")')
                  CALL juDFT_error("atom_sym: mismatch rotated",calledby
     +                 ="atom_sym")
               ENDIF
            ENDDO
!--->       new position
            nneq(ntype) = nneq(ntype) + 1
            tpos(:,nat+nneq(ntype)) = tr
            ity(nat+nneq(ntype)) = nt

         ENDDO opts

         nat = nat + nneq(ntype)
      ENDDO repres_atoms

!--->    allocate arrays in mod_crystal and fill as needed

!      if( ntypd < ntype )then
!        ntypd = ntype
!      endif
!      if ( natd < nat ) THEN
!        natd = nat
!      endif
      ALLOCATE ( neq(ntype), ntyrep(ntype), zatom(ntype) )
      ALLOCATE ( natype(nat), natrep(nat), natmap(nat), pos(3,nat) )

      na=0
      DO n = 1, ntype
         neq(n) = nneq(n)
         zatom(n) = atomid( ity(na+1) )
         ntyrep(n) = na + 1
         DO nn=1,neq(n)
            natype(na+nn) = n
            natrep(na+nn) = na + 1
            natmap(na+nn) = na + nn
            pos(:,na+nn)  = tpos(:,na+nn)
            atomid2(na+nn) = atomid(n)
         ENDDO
         na = na + neq(n)
      ENDDO
      atomid = atomid2


!---> check, that we have  
!--->    either a list of representatives
!--->    or a list of exactly all the atoms in the unit cell.
!---> stop if input is inconsistent 

      IF ( natin > ntype ) THEN
        IF ( .not. all( icount(1:nat)==1 ) ) THEN

          WRITE (outfh,'("ERROR!: atom_sym",/,5x,
     &       "  atom   -   type   -   rep   - times in input")')
          WRITE (outfh,'(4i10)') 
     &       (na, natype(na),natrep(na),icount(na),na=1,nat)

          WRITE (outfh,'("input atom --> atom, is rep, duplicate")') 
          DO nt = 1, natin
            WRITE (errfh,'(2i8,7x,l1,8x,l1)') nt,abs(imap(nt)),
     &         (imap(nt)<0),(icount(abs(imap(nt)))>1)
          ENDDO

          CALL juDFT_error
     +         ("atom_sym:too many representative atoms given in input."
     +         ,calledby ="atom_sym")

        ENDIF
!       l_inipos = .true.
      ENDIF
!     DEALLOCATE ( inipos )

      END SUBROUTINE atom_sym
      END MODULE m_atomsym
