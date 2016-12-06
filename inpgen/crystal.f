      MODULE m_crystal
      use m_juDFT
!********************************************************************
!      generate space group operations from lattice information
!********************************************************************
      CONTAINS
      SUBROUTINE crystal(
     >                   dbgfh,errfh,outfh,dispfh,dispfn,
     >                   cal_symm,cartesian,symor,oldfleur,
     >                   natin,natmax,nop48,
     >                   atomid,atompos,a1,a2,a3,aa,scale,i_c,
     <                   invs,zrfs,invs2,nop,nop2,
     <                   ngen,mmrot,ttr,ntype,nat,nops,
     <                   neq,ntyrep,zatom,natype,natrep,natmap,
     <                   mrot,tau,pos,amat,bmat,omtil)

      USE m_setab
      USE m_lattice, ONLY : angles
      USE m_atomsym
      USE m_spggen
      USE m_closure, ONLY : check_close
      USE m_symproperties
      USE m_generator

      IMPLICIT NONE

!===> Arguments
      LOGICAL, INTENT(IN)    :: cal_symm,cartesian,oldfleur
      INTEGER, INTENT(IN)    :: ngen,natmax,nop48,i_c
      INTEGER, INTENT(IN)    :: dbgfh,errfh,outfh,dispfh ! file handles, mainly 6
      REAL,    INTENT(IN)    :: aa
      LOGICAL, INTENT(INOUT) :: symor                    ! on input: if true, reduce symmetry if oldfleur
                                                         ! on output : if its symmorphic
      INTEGER, INTENT(INOUT) :: natin                    ! might change if atoms are  to be completed
                                                         ! by symmetry operations (if natin < 0 )
      REAL,    INTENT(INOUT) :: atomid(natmax)
      REAL,    INTENT(INOUT) :: atompos(3,natmax)
      REAL,    INTENT(IN)    :: a1(3),a2(3),a3(3)
      REAL,    INTENT(INOUT) :: scale(3)
      INTEGER, INTENT(INOUT) :: mmrot(3,3,nop48)         ! calculated here, if cal_symm is true
      REAL,    INTENT(INOUT) :: ttr(3,nop48)             ! or completed, or if only generators 
      INTEGER, INTENT (OUT)  :: ntype,nat,nops,nop,nop2
      LOGICAL, INTENT (OUT)  :: invs,zrfs,invs2
      REAL,    INTENT (OUT)  :: amat(3,3),bmat(3,3),omtil
      CHARACTER(len=4), INTENT (IN) :: dispfn
!--> actually, intent out:
      INTEGER, ALLOCATABLE :: neq(:), ntyrep(:)              ! these variables are allocated with
      REAL,    ALLOCATABLE :: zatom(:)                       ! dim 'ntype'
      INTEGER, ALLOCATABLE :: natype(:),natrep(:),natmap(:)  ! or  'nat'
      REAL,    ALLOCATABLE :: pos(:,:)                       ! or  '3,nat'
      INTEGER, ALLOCATABLE :: mrot(:,:,:)                    ! or  '3,3,nop'
      REAL,    ALLOCATABLE :: tau(:,:)                       ! or  '3,nop' here, or in atom_sym

                                                         ! are given (if mmrot(1,1,1) = 0 )
! additional arguments are all variables in mod_lattice

!===> Parameters
!  dimensions for group operations, etc.; passed down to other routines
!                                         to force automatic storage
      INTEGER, PARAMETER  :: neig12=12

!===> Local Variables
      LOGICAL lerr,err_setup,invsym
      INTEGER i,j,k,n,m,na,nt,mdet,mtr,nop0,fh,inversionOp
      REAL    t,volume,eps7,eps12

      INTEGER optype(nop48)
      REAL    orth(3,3),ttau(3),rdummy(3,3),as(3,3),bs(3,3)
      REAL    amatinv(3,3),aamat(3,3),bbmat(3,3)
      REAL,   ALLOCATABLE :: poscc(:,:)
      INTEGER, ALLOCATABLE :: multtab(:,:),inv_op(:)

      eps7 = 1.0e-7 ; eps12 = 1.0e-12 
!
!---> set up the as and bs matrices needed to go between
!---> (scaled) cartesian and lattice coordinates
!
      CALL setab_scaled(
     >                  a1,a2,a3,
     <                  as,bs,volume)

      IF (volume < 0.00 ) THEN
        WRITE(6,'(/," Input coordinate system islefthanded;",
     &          " interchange a1 and a2 and try again.")')
         CALL juDFT_error("lefthanded system",calledby="crystal")
      ENDIF
!
!---> modify scale as necessary; note that scale(:) will
!---> be needed to convert to scaled cartesian coordinates
!
      DO i=1,3
         IF ( scale(i) .LT. 0.00 ) THEN
            scale(i) = sqrt( abs(scale(i)) )
         ELSEIF ( abs(scale(i)) .LT. 1.0e-10 ) THEN
            scale(i) = 1.00
         ENDIF
      ENDDO
!
!--->    generate lattice matrices (everything in mod_lattice)
!
      CALL setab(
     >           a1,a2,a3,aa,scale,
     <           amat,bmat,aamat,bbmat,amatinv,omtil)


!---> output: lattice

      WRITE (6,'(//," Lattice information:",/,1x,20("-"))')
      WRITE (6,'(/," overall lattice constant a0     =",
     &                f15.6," bohr")') aa
      WRITE (6,'(/," real-space primitive lattice vectors in units",
     &            " of a_{x,y,z}")')
      WRITE (6,'("      a_1:",3f12.6)') a1
      WRITE (6,'("      a_2:",3f12.6)') a2
      WRITE (6,'("      a_3:",3f12.6)') a3
      WRITE (6,'(/," lattice constants a_x, a_y, a_z =   ",3f12.6)')
     &     aa*scale(:)
      WRITE (6,'(" volume of unit cell (a.u.^3)    =",f15.6)')
     &      omtil

!---> lattice matrices

      WRITE (dbgfh,'(/,"dbg: lattice matrices")')
      WRITE (dbgfh,*) volume
      WRITE (dbgfh,'("dbg:   as      :",3(/,10x,3f12.6) )') as
      WRITE (dbgfh,'("dbg:   bs      :",3(/,10x,3f12.6) )') bs
      WRITE (dbgfh,'("dbg:   amat    :",3(/,10x,3f12.6) )') amat
      WRITE (dbgfh,'("dbg:   bmat    :",3(/,10x,3f12.6) )') bmat
      WRITE (dbgfh,'("dbg:   amatinv :",3(/,10x,3f12.6) )') amatinv
      WRITE (dbgfh,'("dbg:   aamat   :",3(/,10x,3f12.6) )') aamat
      WRITE (dbgfh,'("dbg:   bbmat   :",3(/,10x,3f12.6) )') bbmat
      WRITE (dbgfh,'(/,"dbg:   lattice vectors :")')
      CALL angles( amat )
      WRITE (dbgfh,'(/,"dbg:   reciprocal lattice vectors :")')
      rdummy = transpose( bmat )
      CALL angles( rdummy )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!---> atomic positions:
!--->
!---> atomic positions input can be either in scaled cartesian
!---> or lattice vector units, as determined by logical cartesian.
!---> (for supercells, sometimes more natural to input positions
!---> in scaled cartesian.)
!
!---> if natin < 0, then the representative atoms only are given;
!---> this requires that the space group symmetry be given as input.
!
!---> read in number of atoms or types

      IF ( .not.cal_symm ) THEN  ! generate atoms and 
                                 ! read in symmetry information
         CALL atom_sym(
     >                 dispfh,outfh,errfh,dispfn,natmax,
     X                 natin,atomid,atompos,
     X                 ngen,mmrot,ttr,
     >                 cartesian,i_c,symor,as,bs,nop48,
     <                 ntype,nat,nops,mrot,tau,
     <                 neq,ntyrep,zatom,natype,natrep,natmap,pos)

      ELSE
!----->  allocate arrays in mod_crystal
         nat = natin
         ALLOCATE( natype(nat),natrep(nat),natmap(nat),pos(3,nat) )

!--->    calculate space group symmetry
         CALL spg_gen(
     >                dispfh,outfh,errfh,dispfn,
     >                cartesian,symor,as,bs,scale,
     >                atomid,atompos,natin,nop48,neig12,
     <                ntype,nat,nops,mrot,tau,
     <                neq,ntyrep,zatom,natype,natrep,natmap,pos)
         ! Check whether there is an inversion center that is not at the
         ! origin and if one is found shift the crystal such that the
         ! inversion is with respect to the origin. Then recalculate
         ! symmetry operations.
         inversionOp = -1
         symOpLoop: DO k = 1, nops
            DO i = 1, 3
               DO j = 1, 3
                  IF (i.EQ.j) THEN
                     IF (mrot(i,j,k).NE.-1) CYCLE symOpLoop
                  ELSE
                     IF (mrot(i,j,k).NE.0) CYCLE symOpLoop
                  END IF
                  IF ((i.EQ.3).AND.(j.EQ.3)) THEN
                     inversionOp = k
                     EXIT symOpLoop
                  END IF
               END DO
            END DO
         END DO symOpLoop
         IF (inversionOp.GT.0) THEN
            IF(ANY(ABS(tau(:,inversionOp)).GT.eps7)) THEN
               WRITE(*,*) 'Found inversion center at finite position.'
               WRITE(*,*) 'Shifting crystal by:'
               WRITE(*,'(3f15.10)') 0.5*tau(:,inversionOp)
               WRITE(*,*) ''
               DO k = 1, ABS(natin)
                  atompos(:,k) = atompos(:,k) + 0.5*tau(:,inversionOp)
               END DO
               DEALLOCATE(neq,ntyrep,zatom,mrot,tau)
               CALL spg_gen(
     >                      dispfh,outfh,errfh,dispfn,
     >                      .FALSE.,symor,as,bs,scale,
     >                      atomid,atompos,natin,nop48,neig12,
     <                      ntype,nat,nops,mrot,tau,
     <                      neq,ntyrep,zatom,natype,natrep,natmap,pos)
            END IF
         END IF
      ENDIF

      WHERE ( abs( tau ) < eps7 ) tau = 0.00

!---> atom positions in cartesian coordinates

      ALLOCATE ( poscc(3,nat) )
      WHERE ( abs( pos ) < eps12 ) pos = 0.00
      DO n=1,nat
        poscc(:,n) = matmul( amat , pos(:,n) )
      ENDDO
      WHERE ( abs( poscc ) < eps12 ) poscc = 0.00

!---> check order of atoms

      lerr = .false.
      na = 0
      DO nt = 1, ntype
        DO n = 1, neq(nt)
          na = na + 1
          IF ( natmap(na).ne.na ) lerr = .true.
        ENDDO
      ENDDO
      IF ( lerr ) THEN
        WRITE (errfh,*)
        WRITE (errfh,*) '_err: crystal: ERROR. ',
     &  'Order of atoms is incompatible with fleur21 code.'
        WRITE (errfh,*) '_err: Change order of atoms in input.'
        err_setup = .true.

        WRITE (outfh,1030) '!===> suggested order of atoms'
        WRITE (outfh,1020) nat, ' ! number of atoms'
        WRITE (outfh,1010) '! atomid','x','y','z','type','oldidx'
        na = 0
        DO nt = 1, ntype
          DO n = 1, neq(nt)
            na = na + 1
            i = natmap(na)
            IF ( cartesian ) THEN
              WRITE (outfh,1000) atomid(i),matmul(as,pos(:,i)),nt,i
            ELSE
              WRITE (outfh,1000) atomid(i),pos(:,i),nt,i
            ENDIF
          ENDDO
        ENDDO
 1000   FORMAT (f9.2,3f19.12,' !',2i5)
 1010   FORMAT (a9,a12,a19,a19,a13,a7)
 1020   FORMAT (i10,5x,a)
 1030   FORMAT (a)
      ENDIF

!---> closure, multiplication table and some mapping functions

      ALLOCATE ( inv_op(nops),multtab(nops,nops) )
      CALL check_close( 
     >                 nops,mrot,tau,
     <                 multtab,inv_op,optype)

!---> determine properties of symmmetry operations, 
!---> rearrange mrot,tau as needed

      CALL symproperties(
     >                   nop48,optype,oldfleur,nops,multtab,amat,
     X                   symor,mrot,tau,
     <                   invsym,invs,zrfs,invs2,nop,nop2)


!---> redo to ensure proper mult. table and mapping functions
      CALL check_close(
     >                 nops,mrot,tau,
     <                 multtab,inv_op,optype)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!---> output: space group information
      WRITE (6,'(/," Space group information:",/,1x,24("-"))')
      WRITE (6,'(i5," operations")') nops
      IF ( symor ) THEN
        WRITE (6,'(6x,"space group is symmorphic")')
      ELSE
        WRITE (6,'(6x,"space group is nonsymmorphic")')
      ENDIF
      IF ( invs ) THEN
        WRITE (6,'(18x,"has inversion symmetry")')
        IF ( .not.invsym ) THEN
          WRITE (6,
     &          '(18x,"but inversion is NOT a symmetry operation")')
        ENDIF
      ELSE
        WRITE (6,'(18x,"does NOT have inversion symmetry")')
      ENDIF
      IF ( invs2 ) WRITE (6,'(18x,"has 2d inversion")')
      IF ( zrfs )  WRITE (6,'(18x,"has z-reflection")')

      WRITE (6,'(//," Operations: (in International notation)",/,
     &             " ---------------------------------------",/,3x,
     &    "lattice coordinates",17x,"(scaled) Cartesian coordinates")')

      DO n=1,nops

         IF ( optype(n) .lt. 0 ) THEN
             WRITE (6,'(16x,"_")')
         ELSE
             WRITE (6,*)
         ENDIF
         WRITE (6,
     &            '(" operation",i3,":  ",i1,"  (inverse =",i3,")")')
     &         n,abs(optype(n)),inv_op(n)

         orth = matmul( amat, matmul( mrot(:,:,n) , amatinv ) )
         ttau = matmul( amat, tau(:,n) )
         where( abs( ttau ) < 1.0e-13 ) ttau = 0.00
         WRITE (6,'("  (",3i3," )  (",f6.3," )", 7x,
     &             "  (",3f9.5," )  (",f6.3," )")')
     &    ((mrot(j,i,n),i=1,3),tau(j,n),
     &            (orth(j,i),i=1,3),ttau(j),j=1,3)
      ENDDO

      WRITE(outfh,'(/,"   Multiplcation table: {R_j|t_j}{R_i|t_i}")')
      DO j=1,nops
         WRITE(outfh,'(6x,"operation j=",i2," :",12i4,:/,(22x,12i4))')
     &         j,multtab(j,1:nops)
      ENDDO

!--->    determine a set of generators for this group
      CALL generator(nops,mrot,tau,outfh,errfh)

!---> output: the atomic positions, etc.

      WRITE (6,'(//," Atomic positions:",/,1x,17("-"))')
      WRITE (6,'(" atom types =",i5/,"      total =",i5)') ntype,nat
      WRITE (6,'(/,7x,"lattice coordinates",15x,
     &          "(scaled) Cartesian coordinates   atom")')

      na = 0
      DO nt=1,ntype
         WRITE (6,'(/," atom type",i4,":",2x,
     &            "atomic identification number =",f5.1,
     &            "     representative =",i4)')
     &             nt,zatom(nt),ntyrep(nt)
         DO n=1,neq(nt)
            WRITE (6,'(3f10.6,10x,3f10.6,i7)')
     &           pos(:,natmap(na+n)),poscc(:,natmap(na+n)),natmap(na+n)
         ENDDO
         na = na + neq(nt)
      ENDDO

      DEALLOCATE ( poscc,inv_op,multtab )
      RETURN

      CONTAINS   ! INTERNAL subroutines

      SUBROUTINE setab_scaled(
     >                        a1,a2,a3,
     <                        as,bs,volume)

!*****************************************************************
!     set up matrices needed to convert rotation matrices between
!     (scaled) cartesian and lattice coordinates
!*****************************************************************
      IMPLICIT NONE

      REAL, INTENT (IN)  :: a1(3),a2(3),a3(3)
      REAL, INTENT (OUT) :: as(3,3),bs(3,3),volume

      as(1,1) = a1(1)
      as(2,1) = a1(2)
      as(3,1) = a1(3)
      as(1,2) = a2(1)
      as(2,2) = a2(2)
      as(3,2) = a2(3)
      as(1,3) = a3(1)
      as(2,3) = a3(2)
      as(3,3) = a3(3)

      volume  = a1(1)*a2(2)*a3(3) + a2(1)*a3(2)*a1(3) +
     &          a3(1)*a1(2)*a2(3) - a1(3)*a2(2)*a3(1) -
     &          a2(3)*a3(2)*a1(1) - a3(3)*a1(2)*a2(1)

      bs(1,1) = (a2(2)*a3(3) - a2(3)*a3(2))/volume ! b1(1)
      bs(1,2) = (a2(3)*a3(1) - a2(1)*a3(3))/volume ! b1(2)
      bs(1,3) = (a2(1)*a3(2) - a2(2)*a3(1))/volume ! b1(3)
      bs(2,1) = (a3(2)*a1(3) - a3(3)*a1(2))/volume ! b2(1)
      bs(2,2) = (a3(3)*a1(1) - a3(1)*a1(3))/volume ! b2(2)
      bs(2,3) = (a3(1)*a1(2) - a3(2)*a1(1))/volume ! b2(3)
      bs(3,1) = (a1(2)*a2(3) - a1(3)*a2(2))/volume ! b3(1)
      bs(3,2) = (a1(3)*a2(1) - a1(1)*a2(3))/volume ! b3(2)
      bs(3,3) = (a1(1)*a2(2) - a1(2)*a2(1))/volume ! b3(3)

      RETURN
      END SUBROUTINE setab_scaled

      END SUBROUTINE crystal
      END MODULE m_crystal
