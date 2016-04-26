      MODULE m_spggen
      use m_juDFT
!********************************************************************
!     calculates the space group operations given the lattice vectors
!     and the atomic positions.                              mw 12-99
!********************************************************************
      CONTAINS
      SUBROUTINE spg_gen(
     >                   dispfh,outfh,errfh,dispfn,
     >                   cartesian,symor,as,bs,scale,
     >                   atomid,atompos,natin,nop48,neig12,
     <                   ntype,nat,nops,mrot,tau,
     <                   neq,ntyrep,zatom,natype,natrep,natmap,pos)

      USE m_bravaissymm
      USE m_closure, ONLY : close_pt, closure
      USE m_supercheck

      IMPLICIT NONE

      INTEGER, INTENT (IN)    :: dispfh,outfh,errfh
      INTEGER, INTENT (IN)    :: natin
      INTEGER, INTENT (IN)    :: nop48,neig12
      LOGICAL, INTENT (IN)    :: cartesian,symor
      REAL,    INTENT (IN)    :: as(3,3),bs(3,3),scale(3)
      REAL,    INTENT (INOUT) :: atompos(3,natin)
      REAL,    INTENT (IN)    :: atomid(natin)
      INTEGER, INTENT (OUT)   :: ntype,nat,nops
      CHARACTER(len=4), INTENT (IN) :: dispfn
!--> actually, intent out:
      INTEGER, ALLOCATABLE :: neq(:), ntyrep(:)    ! here, these variables are allocated with
      REAL,    ALLOCATABLE :: zatom(:)             ! the dimension 'ntype'
      INTEGER, ALLOCATABLE :: natype(:),natrep(:),natmap(:)  ! or  'nat'
      REAL,    ALLOCATABLE :: pos(:,:)                       ! or  '3,nat'
      INTEGER, ALLOCATABLE :: mrot(:,:,:)                    ! or  '3,3,nop'
      REAL,    ALLOCATABLE :: tau(:,:)                       ! or  '3,nop'


!   point group operations from bravais lattice
      INTEGER mops,mmrot(3,3,nop48),m_inv(nop48)
      INTEGER ncyl(nop48)
      INTEGER index_op(nop48),num_tr(nop48)

      INTEGER mtmp(3,3),mp(3,3)
      INTEGER i,j,k,n,mop,nc,new,ns,nt,ntypm,ind(1)
      INTEGER ity(natin),npaint,ipaint(natin)
      INTEGER ios,istep0
      CHARACTER(len=30) :: filen

      REAL    posr(3,natin),rtau(3),tr(3),disp(3,natin)
      REAL    ttau(3,nop48),trs(3,natin)
      REAL    eps7

      LOGICAL lnew,lclose
      LOGICAL l_exist

      INTEGER, ALLOCATABLE :: mtable(:,:)
      REAL,    ALLOCATABLE :: inipos(:,:)

      eps7= 1.0e-7 ; istep0 = 0

!---> determine the point group of the Bravais lattice
      CALL bravais_symm(
     >                  as,bs,scale,nop48,neig12,
     <                  mops,mmrot)

      ALLOCATE ( mtable(mops,mops) )
      CALL close_pt(                        ! check closure and get
     >              mops,mmrot,             ! multiplication table
     <              mtable)
      DO mop = 1, mops                      ! inverses
         DO i = 1, mops
            IF ( mtable(i,mop) == 1 ) THEN
               m_inv(mop) = i
               EXIT
            ENDIF
         ENDDO
      ENDDO

!---> read in atomic positions and shift to (-1/2,1/2] in lattice
!---> coords. also read in identification (atomic) number (atomid)
!---> to distinguish different atom types (need not be atomic number)
      DO n=1,natin
         IF (cartesian) THEN  ! convert to lattice coords. if necessary
            atompos(:,n) = matmul( bs, atompos(:,n) )
         ENDIF
         pos(:,n) = atompos(:,n) - anint( atompos(:,n) - eps7 )
      ENDDO

!---> store the positions (in lattice coord.s) given in the input file
!     ALLOCATE ( inipos(3,natin) )
!     DO n=1,natin
!        inipos(:,n) = pos(:,n)
!     enddo
!     l_inipos = .true.

!--->  check whether a displacement file exists (dispfh)
!---> these displacements should be in the same units as input file,
!---> i.e., lattice  or scaled cartesian
!
      l_exist = .false.
      INQUIRE (file=dispfn,exist=l_exist)
      IF (l_exist) then
        OPEN(dispfh,file=dispfn,status='old',form='formatted')
        DO n = 1, natin
          READ( dispfh,*,iostat=ios) tr(1:3)
          IF (ios .ne. 0) EXIT
          disp(1:3,n) = tr(1:3)
          IF (cartesian) THEN    ! convert to lattice coords. if necessary
            tr = matmul( bs, tr )
          ENDIF
          pos(:,n) = pos(:,n) + tr(:)
          pos(:,n) = pos(:,n) - anint( pos(:,n)- eps7 )
        ENDDO
        CLOSE (dispfh)
        IF ( ios==0 ) THEN
          WRITE (filen,'(a,"_",i2.2)') trim(adjustl(dispfn)),istep0+1
          OPEN (dispfh,file=filen,status='unknown',form='formatted')
          DO n = 1, natin
            WRITE (dispfh,'(3f16.8)') disp(1:3,n)
          ENDDO
          CLOSE (dispfh)
        ENDIF
      ENDIF ! l_exist

!---> sanity check: atoms must be distinct
      DO n = 1, natin
         DO i = n+1, natin
            IF ( all( abs( pos(:,n) - pos(:,i) ) < eps7 ) ) then
               WRITE(6,'(/," Error in atomic positions: atoms",i3,
     &             " and",i3," (at least) are the same")') n,i
                CALL juDFT_error("atom positions",calledby="spg_gen")
            ENDIF
         ENDDO
      ENDDO

!---> determine the number of distinct atoms based on atomic number,
!---> etc. (not necessarily symmetry inequivalent)

      ntypm = 1
      ity(1) = 1
      DO n=2, natin
         lnew = .true.
         DO i=1,n-1
            IF ( abs( atomid(i)-atomid(n) ) < eps7 ) THEN
               ity(n) = ity(i)
               lnew = .false.
               EXIT
            ENDIF
         ENDDO
         IF (lnew) then
            ntypm = ntypm + 1
            ity(n) = ntypm
         ENDIF
      ENDDO

!---> determine the order of cyclic (sub)group for each operation

      ncyl(1) = 1
      DO mop=1,mops
         ncyl(mop) = 2
         i = mop
         DO              ! multiply on left until identity
           i = mtable(mop, i)
           IF ( i == 1 ) EXIT
           ncyl(mop) = ncyl(mop)+1
         ENDDO
      ENDDO

!---> check if this is a supercell

      nat = natin
      CALL super_check(
     >                 nat,pos,ity,ntypm,
     <                 ns,trs)

!---> determine the independent atoms in the supercell and
!---> mark off the equivalent atoms

      IF ( ns > 1 ) then
         ipaint(1:nat) = 0
         DO n = 1, nat
            IF ( ipaint(n) < 0 ) CYCLE
            ipaint(n) = 0
            DO i = 2, ns        ! want to mark the equivalent atoms
               tr(1:3) = pos(1:3,n) + trs(1:3,i)
               tr(:) = tr(:) - anint( tr(:) - eps7 )
               DO j = n+1, nat
                  IF ( all( abs( tr(:) - pos(:,j) ) < eps7 ) ) THEN
                     ipaint(j) = -1
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF


!---> loop over operations (first operation is identity) and
!---> determine the possible non-primitive translations
!---> for supercells, it is possible to have more than 1 translation
!---> that also satisfies the cyclic condition; to generate subgroup
!---> consistent with the overall cell (and other conventions), we
!---> break the supercell symmetry if necessary; in this case, only
!---> a single non-primitive tranlation per rotation is allowed.
!---> (if one wants the full symmetry of the supercell, can easily
!---> modify the code here to calculate all translations.)

      npaint = 0
  100 CONTINUE

      num_tr(1) = 1
      ttau(1:3,1) = (/ 0.00 , 0.00 , 0.00 /)

      num_tr(2:mops) = 0
      DO mop = 2,  nop48            ! DO mop = 2, mops seems to be equivalent,
       IF (mop <= mops ) THEN       ! but not for ifc (7.0) with -O3

         IF ( num_tr(mop) < 0 ) cycle  ! this operation removed previously

!--->    rotate all atoms:
         DO n=1,nat
            posr(:,n) = matmul( mmrot(:,:,mop) , pos(:,n) )
         ENDDO

!--->    generate possible non-symmorphic pieces
         DO j=1,nat
            DO i=1,nat
               IF ( ity(i) .ne. ity(j) ) CYCLE

               tr(1:3) = pos(1:3,j) - posr(1:3,i)
               tr(1:3) = tr(1:3) - anint( tr(1:3) - eps7 )
               rtau(1:3) = tr(1:3)

!--->          check that this is consistent with cyclic part of the group
               DO nc = 1,ncyl(mop)-1
                  rtau(:) = matmul( mmrot(:,:,mop) , rtau(:) ) + tr(:)
               ENDDO
               rtau(1:3) = rtau(1:3) - anint( rtau(1:3) - eps7 )
               IF ( any( abs(rtau(:)) > eps7 ) ) CYCLE  ! skip if not 0

!--->          test if this vector brings the system into registry
               IF ( l_rotmatch(tr) ) THEN   ! new translation
                  num_tr(mop) = num_tr(mop) + 1
                  ttau(1:3,mop) = tr(1:3)
                  EXIT                  ! have one, go to next operation
               ENDIF

            ENDDO
            IF ( num_tr(mop) > 0 ) EXIT  ! have one, go to next operation
         ENDDO
         IF ( num_tr(mop) < 1 ) num_tr(m_inv(mop)) = -1 ! remove inverse also
       ENDIF
      ENDDO

      nops = 0
      DO mop = 1, mops
         IF ( num_tr(mop) < 1 ) CYCLE
         nops = nops + 1
         index_op(nops) = mop
      ENDDO


!---> check closure of group
      CALL closure(
     >             mops,mmrot,ttau,nops,index_op,
     <             lclose)

      IF ( ( ns==1 ) .AND. ( .not. lclose ) ) THEN
         WRITE (6,'(/," Congratulations, you have found a system (not",
     &    " a supercell) that breaks the algorithms. Sorry...")')
          CALL juDFT_error("Program failed :(",calledby="spg_gen")
      ENDIF

!---> supercells: if we have determined a (sub)group directly, great;
!---> however, may not have found a subgroup consistent with the
!---> supercell, so need to search harder. to do this, we will
!---> "paint" the different atoms in the primitive cell and see
!---> which choice gives the largest number of operations. because
!---> of this requirement, there is a "go to" construct. there is
!---> extra work done for the case of supercells, but seems to
!---> be necessary

      IF ( ns > 1 ) THEN

         IF ( ( .not. lclose ) .or. ( npaint > 0 ) ) THEN
            IF ( npaint == 0 ) THEN
               WRITE (6,'(/," Generating subgroups of supercell...")')
            ELSE
               IF (lclose) THEN
                  ipaint(npaint) = nops
                  WRITE (6,'(5x,"painting atom",i5,":",i3,
     &                  " operations")') npaint,nops
               ELSE
                  ipaint(npaint) = -1
                  WRITE (6,'(5x,"painting atom",i5,": not a group")')
     &                 npaint
               ENDIF
            ENDIF
            ity(1:nat) = abs( ity(1:nat) )
            npaint = -1
            DO i = 1, nat
               IF ( ipaint(i) .ne. 0 ) CYCLE
               npaint = i
               ity(i) = -ity(i)
               exit
            ENDDO
            IF ( npaint > 0 ) go to 100
         ENDIF

         IF ( npaint == -1 ) THEN  ! (re)calculate group best one
            ind = maxloc( ipaint(1:3) )  ! ind MUST be an array
            IF ( ipaint(ind(1)) < 1 ) THEN
               WRITE (6,'(/," Algorithm didnot work. Sorry...")')
               CALL juDFT_error("Supercell failure ;(",
     &                   calledby ="spg_gen")
            ENDIF
            ity(ind(1)) = -abs(ity(ind(1)))
            npaint = -2
            GO TO 100
         ENDIF

      ENDIF

      ity(1:nat) = abs( ity(1:nat) )

!---> set up various mapping functions, etc., in mod_crystal
!---> allocate arrays for space group information (mod_spgsym)
      DEALLOCATE ( mtable)          ! deallocate to avoid fragmenting heap
!      if( nopd < nops )then
!        nopd = nops
!      endif
      ALLOCATE ( mrot(3,3,nops),tau(3,nops) )

      DO n=1,nops
         mrot(:,:,n) = mmrot(:,:,index_op(n))
         tau(:,n) = ttau(:,index_op(n))
      ENDDO

      IF ( symor ) THEN
!--->   reduce symmetry to the largest symmorphic subgroup
        j = 1
        DO i = 1, nops
          IF ( all ( abs( tau(:,i) ) < eps7 ) ) THEN
            IF ( j<i ) THEN
              mrot(:,:,j) = mrot(:,:,i)
            ENDIF
            j = j + 1
          ENDIF
        ENDDO
        tau = 0.00
        IF ( nops > j-1 ) THEN
          WRITE (outfh,*) 'System has non-symmorphic symmetry with',
     &                   nops, 'operations.'
          nops = j - 1
          WRITE(outfh,*) 'Symmetry reduced to symmorphic subgroup with',
     &   nops, 'operations.'
        ENDIF

      ENDIF ! symor

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!--->    at this point, the symmetry is correct (assumed here)
!--->    now fill in arrays in mod_crystal

      natype(1:nat) = 0
      ntype = 0
      DO i =1,nat
         IF ( natype(i) .ne. 0 ) cycle
         ntype = ntype + 1   ! new atom type
         natype(i) = ntype   ! atom type
         natrep(i) = i       ! this is the representative
!--->    rotate representative and get symmetry equavalent atoms
         DO n=1,nops
            tr(:) = matmul( mrot(:,:,n) , pos(:,i) ) + tau(:,n)
            tr(:) = tr(:) - anint( tr(:) -eps7 )
!--->       this (rotated) atom done already? (correct symmetry assumed)
            DO j=i+1,nat
               IF ( natype(j) .ne. 0 ) CYCLE
               IF ( ity(j) .ne. ity(i) ) CYCLE
               IF ( any( abs( tr(:) - pos(:,j) ) > eps7 ) ) CYCLE
               natrep(j) = i      ! representative atom
               natype(j) = ntype  ! atom type
               EXIT
            ENDDO
         ENDDO
      ENDDO

!      if( ntypd < ntype )then
!        ntypd = ntype
!      endif
      ALLOCATE( neq(ntype),ntyrep(ntype),zatom(ntype) )

      neq(1:ntype) = 0
      ntyrep(1:ntype) = 0
      DO n=1,nat
         neq( natype(n) ) = neq( natype(n) ) + 1
         zatom( natype(n) ) = atomid(n)
         IF ( ntyrep( natype(n) ) == 0 ) ntyrep( natype(n) ) = n
      ENDDO

      natmap(1:nat) = 0
      j = 0
      DO nt = 1,ntype
         DO n=1,nat
            IF ( natype(n) == nt ) THEN
               j = j+ 1
               natmap(j) = n
            ENDIF
         ENDDO
      ENDDO

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      CONTAINS ! INTERNAL SUBROUTINES

      LOGICAL FUNCTION l_rotmatch(tr)
!********************************************************************
!     determines whether the vector tr will bring the rotated
!     positions posr into regestry with the original
!********************************************************************
      IMPLICIT NONE

      REAL, INTENT (IN) :: tr(3)

      INTEGER i,j,in
      REAL    rp(3)

      l_rotmatch = .false.

      DO i=1,nat
!--->       rotated and shifted atom, reduced to (-1/2,1/2]
         rp(:) = posr(:,i) + tr(:) - anint( posr(:,i) + tr(:) - eps7 )
!--->       find which atom, if any, this matches
         in = 0
         DO j=1,nat
            IF ( ity(i).ne.ity(j) ) CYCLE
            IF ( all( abs(pos(:,j)-rp(:)) < eps7 ) ) THEN
               in = j
               EXIT
            ENDIF
         ENDDO
         IF (in == 0 ) RETURN
      ENDDO

!--->   here only if everything matches
      l_rotmatch = .true.

      END FUNCTION l_rotmatch

      END SUBROUTINE spg_gen
      END MODULE m_spggen

