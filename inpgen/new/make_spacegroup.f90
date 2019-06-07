

MODULE m_make_spacegroup
  USE m_juDFT
  !********************************************************************
  !  Generate the spacegroup given the cell and the atomic positions/ids
  !  Takes into account the symmetry reductions needed for films,soc&SSDW
  !  Based on code by
  !         mw 12-99
  !         Modified by GM (2016)
  !********************************************************************
CONTAINS
  SUBROUTINE make_spacegroup(film,noco,cell,atompos,atomid,sym)

    USE m_bravaissymm
    USE m_closure, ONLY : close_pt, closure
    USE m_supercheck

    IMPLICIT NONE
    LOGICAL,INTENT(in)         :: film
    TYPE(t_noco),INTENT(in)    :: noco
    TYPE(t_cell),intent(in)    :: cell
    REAL,INTENT(IN)            :: pos(:,:),atomid(:)
    TYPE(t_sym),INTENT(inout)  :: sym

    INTEGER, ALLOCATABLE :: ntyrep(:)    ! here, these variables are allocated with
    INTEGER, ALLOCATABLE :: natype(:),natrep(:)  ! or  'nat'
    

  
    INTEGER,PARAMETER::nop48=48
    !   point group operations from bravais lattice
    INTEGER mops,mmrot(3,3,nop48),m_inv(nop48)
    INTEGER ncyl(nop48)
    INTEGER index_op(nop48),num_tr(nop48)

    INTEGER mtmp(3,3),mp(3,3),ntype,nat
    INTEGER u,v,w
    INTEGER i,j,k,n,mop,nc,NEW,ns,nt,ntypm,ind(1),iTr,maxTrVecs
    INTEGER ity(size(atomid)),npaint,ipaint(size(atomid))
    INTEGER ios,istep0, binSize, maxBinSize
    INTEGER locBinDim(3), secondAtom(size(atomid))
    INTEGER binDim(3), iBin(3)
    INTEGER trIndices(size(atomid))
    CHARACTER(len=30) :: filen

    REAL    posr(3,size(atomid)),rtau(3),tr(3)
    REAL    ttau(3,nop48),trs(3,size(atomid))
    REAL    eps7
    REAL    trVecs(3,size(atomid))

    LOGICAL lnew,lclose, foundOne, boundary(3)
    LOGICAL l_exist

    INTEGER, ALLOCATABLE :: mtable(:,:), binSizes(:,:,:)
    INTEGER, ALLOCATABLE :: atomIndexBins(:,:,:,:)
    REAL,    ALLOCATABLE :: inipos(:,:)

    eps7= 1.0e-7 ; istep0 = 0

    !---> determine the point group of the Bravais lattice
    CALL bravais_symm(cell, mops,mmrot)
    !reduce symmetry in special setups
    ALLOCATE ( error(mpos) ); error=.FALSE.
    ! reduce symmetry if SSDW calculation
    IF (noco%l_ss) CALL ss_sym(mops,mmrot,noco%qss,error)
    IF (noco%l_soc) CALL soc_sym(mops,mmrot,noco%theta,noco%phi,amat, error)
    IF (film) CALL film_sym(mops,mmrot,error)
    n=0 !Keep only operations without error
    DO i=1,mops
       IF (.NOT.error(i)) THEN
          n=n+1
          mmrot(:,:,n)=mmrot(:,:,i)
       ENDIF
    ENDIF
    mops=n


    
    ALLOCATE ( mtable(mops,mops) )
    ! check closure and get  multiplication table
    CALL close_pt(mops,mmrot, mtable)
    DO mop = 1, mops                      ! inverses
       DO i = 1, mops
          IF ( mtable(i,mop) == 1 ) THEN
             m_inv(mop) = i
             EXIT
          ENDIF
       ENDDO
    ENDDO

    !---> sanity check: atoms must be distinct
    DO n = 1, size(atomid)
       DO i = n+1, size(atomid)
          IF ( ALL( ABS( pos(:,n) - pos(:,i) ) < eps7 ) ) THEN
             WRITE(6,'(/," Error in atomic positions: atoms",i5," and",i5," (at least) are the same")') n,i
             CALL juDFT_error("atom positions",calledby="spg_gen")
          ENDIF
       ENDDO
    ENDDO

    !---> determine the number of distinct atoms based on atomic number,
    !---> etc. (not necessarily symmetry inequivalent)

    ntypm = 1
    ity(1) = 1
    DO n=2, size(atomid)
       lnew = .TRUE.
       DO i=1,n-1
          IF ( ABS( atomid(i)-atomid(n) ) < eps7 ) THEN
             ity(n) = ity(i)
             lnew = .FALSE.
             EXIT
          ENDIF
       ENDDO
       IF (lnew) THEN
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

    nat = size(atomid)
    CALL super_check(nat,pos,ity,ntypm,ns,trs)

    !---> determine the independent atoms in the supercell and
    !---> mark off the equivalent atoms

    IF ( ns > 1 ) THEN
       ipaint(1:nat) = 0
       DO n = 1, nat
          IF ( ipaint(n) < 0 ) CYCLE
          ipaint(n) = 0
          DO i = 2, ns        ! want to mark the equivalent atoms
             tr(1:3) = pos(1:3,n) + trs(1:3,i)
             tr(:) = tr(:) - ANINT( tr(:) - eps7 )
             DO j = n+1, nat
                IF ( ALL( ABS( tr(:) - pos(:,j) ) < eps7 ) ) THEN
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

          IF ( num_tr(mop) < 0 ) CYCLE  ! this operation removed previously

          !--->    rotate all atoms:
          DO n=1,nat
             posr(:,n) = MATMUL( mmrot(:,:,mop) , pos(:,n) )
             posr(:,n) = posr(:,n) - ANINT(posr(:,n) - eps7)
          ENDDO

          !        Start code section A (replacing the commented part following the section)

          !        Determine possible translation vectors. Each translation vector has
          !        to work for all atoms.

          trVecs = 0.0
          maxTrVecs = 0
          secondAtom = 0

          !!       1. Determine all possible translation vectors for the first atom

          trIndices = -1
          j = 1
          DO i = 1, nat
             IF (ity(i).NE.ity(j)) CYCLE
             tr(1:3) = pos(1:3,i) - posr(1:3,j)
             tr(1:3) = tr(1:3) - ANINT(tr(1:3) - eps7)
             maxTrVecs = maxTrVecs + 1
             trVecs(:,maxTrVecs) = tr(:)
             trIndices(maxTrVecs) = maxTrVecs
             secondAtom(maxTrVecs) = i
          END DO

          !!       2. Sort atoms into "location bins"
          !!          (position vectors are in the region -0.5 to 0.5)

          binDim(:) = CEILING(size(atomid)**(1.0/3.0)*0.5)
          !TODO: The dimension binDim should better be adjusted to the unit cell shape.
          !      This simple setting might be bad for very elongated unit cells.

          ALLOCATE(binSizes(-binDim(1)-1:binDim(1)+1, -binDim(2)-1:binDim(2)+1, -binDim(3)-1:binDim(3)+1))
          binSizes = 0
          DO n = 1, nat
             iBin(:) = NINT(binDim(:) * pos(:,n) / 0.501)

             boundary(:) = (ABS(pos(:,n))-0.5).LE.eps7
             DO i = -1, 1, 2
                IF((i.EQ.-1).AND.(.NOT.boundary(1))) CYCLE
                DO j = -1, 1, 2
                   IF((j.EQ.-1).AND.(.NOT.boundary(2))) CYCLE
                   DO k = -1, 1, 2
                      IF((k.EQ.-1).AND.(.NOT.boundary(3))) CYCLE
                      binSize = binSizes(i*iBin(1),j*iBin(2),k*iBin(3))
                      binSize = binSize + 1
                      binSizes(i*iBin(1),j*iBin(2),k*iBin(3)) = binSize
                   END DO
                END DO
             END DO

          END DO

          maxBinSize = 0
          DO i = -binDim(1)-1, binDim(1)+1
             DO j = -binDim(2)-1, binDim(2)+1
                DO k = -binDim(3)-1, binDim(3)+1
                   IF (binSizes(i,j,k).GT.maxBinSize) THEN
                      maxBinSize = binSizes(i,j,k)
                   END IF
                END DO
             END DO
          END DO

          ALLOCATE(atomIndexBins(maxBinSize,-binDim(1)-1:binDim(1)+1, -binDim(2)-1:binDim(2)+1, -binDim(3)-1:binDim(3)+1))

          binSizes = 0
          DO n = 1, nat
             iBin(:) = NINT(binDim(:) * pos(:,n) / 0.501)

             boundary(:) = (ABS(pos(:,n))-0.5).LE.eps7
             DO i = -1, 1, 2
                IF((i.EQ.-1).AND.(.NOT.boundary(1))) CYCLE
                DO j = -1, 1, 2
                   IF((j.EQ.-1).AND.(.NOT.boundary(2))) CYCLE
                   DO k = -1, 1, 2
                      IF((k.EQ.-1).AND.(.NOT.boundary(3))) CYCLE
                      binSize = binSizes(i*iBin(1),j*iBin(2),k*iBin(3))
                      binSize = binSize + 1
                      binSizes(i*iBin(1),j*iBin(2),k*iBin(3)) = binSize
                      atomIndexBins(binSize,i*iBin(1),j*iBin(2), k*iBin(3)) = n
                   END DO
                END DO
             END DO
          END DO

          !!       3. Check for every other atom which of the first atom's translation
          !!          vectors are applicable.

          DO j = 2, nat
             iTr = 1
             DO WHILE (iTr.LE.maxTrVecs)
                tr(1:3) = posr(1:3,j) + trVecs(1:3,trIndices(iTr))
                tr(1:3) = tr(1:3) - ANINT(tr(1:3) - eps7)

                iBin(:) = NINT(binDim(:) * tr(:) / 0.501)

                foundOne = .FALSE.

                !Search for atoms in the nearest bin
                DO k = 1,binSizes(iBin(1),iBin(2),iBin(3))
                   i = atomIndexBins(k,iBin(1),iBin(2),iBin(3))
                   rtau(:) = tr(:)-pos(:,i)
                   rtau(:) = rtau(:) - ANINT(rtau(:) - eps7)
                   IF(ALL(ABS(rtau(:)).LE.eps7)) THEN
                      IF (ity(i).NE.ity(j)) EXIT
                      foundOne = .TRUE.
                      EXIT
                   END IF
                END DO

                IF(.NOT.foundOne) THEN
                   !Search for atoms in the surrounding bins
                   binLoop: DO u = -1, 1
                      DO v = -1, 1
                         DO w = -1, 1
                            IF((u.EQ.0).AND.(v.EQ.0).AND.(w.EQ.0)) CYCLE
                            DO k = 1,binSizes(iBin(1)+u,iBin(2)+v,iBin(3)+w)
                               i = atomIndexBins(k,iBin(1)+u,iBin(2)+v, iBin(3)+w)
                               rtau(:) = tr(:)-pos(:,i)
                               rtau(:) = rtau(:) - ANINT(rtau(:) - eps7)
                               IF(ALL(ABS(rtau(:)).LE.eps7)) THEN
                                  IF (ity(i).NE.ity(j)) EXIT binLoop
                                  foundOne = .TRUE.
                                  EXIT binLoop
                               END IF
                            END DO
                         END DO
                      END DO
                   END DO binLoop
                END IF

                IF (foundOne) THEN
                   iTr = iTr + 1
                ELSE
                   trIndices(iTr) = trIndices(maxTrVecs)
                   maxTrVecs = maxTrVecs - 1
                END IF
             END DO
          END DO

          !        Check which translation vectors are consistent with the cyclic
          !        part of the group

          DO iTr = 1, maxTrVecs
             j = 1
             i = secondAtom(trIndices(iTr))
             tr(:) = trVecs(:,trIndices(iTr))
             rtau(1:3) = tr(1:3)

             !--->       check that this is consistent with cyclic part of the group
             DO nc = 1,ncyl(mop)-1
                rtau(:) = MATMUL( mmrot(:,:,mop) , rtau(:) ) + tr(:)
             END DO
             rtau(1:3) = rtau(1:3) - ANINT( rtau(1:3) - eps7 )
             IF ( ANY( ABS(rtau(:)) > eps7 ) ) CYCLE  ! skip if not 0

             num_tr(mop) = num_tr(mop) + 1
             ttau(1:3,mop) = tr(1:3)
             EXIT                  ! have one, go to next operation
          END DO

          DEALLOCATE(atomIndexBins,binSizes)

          !        End code section A (replacing the commented part following the section)

          !--->    generate possible non-symmorphic pieces
          !        DO j=1,nat
          !           DO i=1,nat
          !              IF ( ity(i) .ne. ity(j) ) CYCLE
          !
          !              tr(1:3) = pos(1:3,j) - posr(1:3,i)
          !              tr(1:3) = tr(1:3) - anint( tr(1:3) - eps7 )
          !              rtau(1:3) = tr(1:3)

          !--->          check that this is consistent with cyclic part of the group
          !              DO nc = 1,ncyl(mop)-1
          !                 rtau(:) = matmul( mmrot(:,:,mop) , rtau(:) ) + tr(:)
          !              ENDDO
          !              rtau(1:3) = rtau(1:3) - anint( rtau(1:3) - eps7 )
          !              IF ( any( abs(rtau(:)) > eps7 ) ) CYCLE  ! skip if not 0

          !--->          test if this vector brings the system into registry
          !              IF ( l_rotmatch(tr) ) THEN   ! new translation
          !                 num_tr(mop) = num_tr(mop) + 1
          !                 ttau(1:3,mop) = tr(1:3)
          !                 EXIT                  ! have one, go to next operation
          !              ENDIF
          !
          !           ENDDO
          !           IF ( num_tr(mop) > 0 ) EXIT  ! have one, go to next operation
          !        ENDDO

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
    CALL closure(mops,mmrot,ttau,nops,index_op, lclose)

    IF ( ( ns==1 ) .AND. ( .NOT. lclose ) ) THEN
       WRITE (6,'(/," Congratulations, you have found a system (not"," a supercell) that breaks the algorithms. Sorry...")')
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

       IF ( ( .NOT. lclose ) .OR. ( npaint > 0 ) ) THEN
          IF ( npaint == 0 ) THEN
             WRITE (6,'(/," Generating subgroups of supercell...")')
          ELSE
             IF (lclose) THEN
                ipaint(npaint) = nops
                WRITE (6,'(5x,"painting atom",i5,":",i3," operations")') npaint,nops
             ELSE
                ipaint(npaint) = -1
                WRITE (6,'(5x,"painting atom",i5,": not a group")') npaint
             ENDIF
          ENDIF
          ity(1:nat) = ABS( ity(1:nat) )
          npaint = -1
          DO i = 1, nat
             IF ( ipaint(i) .NE. 0 ) CYCLE
             npaint = i
             ity(i) = -ity(i)
             EXIT
          ENDDO
          IF ( npaint > 0 ) go to 100
       ENDIF

       IF ( npaint == -1 ) THEN  ! (re)calculate group best one
          ind = MAXLOC( ipaint(1:3) )  ! ind MUST be an array
          IF ( ipaint(ind(1)) < 1 ) THEN
             WRITE (6,'(/," Algorithm didnot work. Sorry...")')
             CALL juDFT_error("Supercell failure ;(", calledby ="spg_gen")
          ENDIF
          ity(ind(1)) = -ABS(ity(ind(1)))
          npaint = -2
          GO TO 100
       ENDIF

    ENDIF

    ity(1:nat) = ABS( ity(1:nat) )

    !---> set up various mapping functions, etc., in mod_crystal
    !---> allocate arrays for space group information (mod_spgsym)
    DEALLOCATE ( mtable)          ! deallocate to avoid fragmenting heap
    !      if( nopd < nops )then
    !        nopd = nops
    !      endif


    
    ALLOCATE ( sym%mrot(3,3,nops),sym%tau(3,nops) )
    sym%nops=nops
    DO n=1,nops
       sym%mrot(:,:,n) = mmrot(:,:,index_op(n))
       sym%tau(:,n) = ttau(:,index_op(n))
    ENDDO
    IF ( sym%symor ) THEN
       !--->   reduce symmetry to the largest symmorphic subgroup
       j = 1
       DO i = 1, sym%nops
          IF ( ALL ( ABS( sym%tau(:,i) ) < eps7 ) ) THEN
             IF ( j<i ) THEN
                sym%mrot(:,:,j) = sym%mrot(:,:,i)
             ENDIF
             j = j + 1
          ENDIF
       ENDDO
       sym%tau = 0.00
       IF ( sym%nops > j-1 ) THEN
          WRITE (6,*) 'System has non-symmorphic symmetry with', nops, 'operations.'
          sym%nops = j - 1
          WRITE(6,*) 'Symmetry reduced to symmorphic subgroup with', nops, 'operations.'
       ENDIF

    ENDIF ! sym%symor
    WHERE ( ABS( sym%tau ) < eps7 ) sym%tau = 0.00
    

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

      l_rotmatch = .FALSE.

      DO i=1,nat
         !--->       rotated and shifted atom, reduced to (-1/2,1/2]
         rp(:) = posr(:,i) + tr(:) - ANINT( posr(:,i) + tr(:) - eps7 )
         !--->       find which atom, if any, this matches
         in = 0
         DO j=1,nat
            IF ( ity(i).NE.ity(j) ) CYCLE
            IF ( ALL( ABS(pos(:,j)-rp(:)) < eps7 ) ) THEN
               in = j
               EXIT
            ENDIF
         ENDDO
         IF (in == 0 ) RETURN
      ENDDO

      !--->   here only if everything matches
      l_rotmatch = .TRUE.

    END FUNCTION l_rotmatch

  END SUBROUTINE make_spacegroup
END MODULE m_make_spacegroup
