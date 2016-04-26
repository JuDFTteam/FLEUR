      MODULE m_localsym
!*********************************************************************
!     generate the lattice harmonics appropriate to the local symmetry
!     of the atoms, given the space group  operations, bravias lattice,
!     and the atomic positions.
!
!     the coordinate system for vector and operations must be given.
!     here it is assumed that vectors and operations are given in
!     terms of INTERNAL coordinates, i.e., in terms of a1,a2,a3. this
!     is the more complicated case; for everything in cartesian
!     coordinates, the needed changes are rather obvious.
!
!     input:
!        lmax     max. l to calculate for each atom type
!        lmaxd    max. l (needed to set array sizes)
!        a1,a2,a3 primitive translation vectors; CARTESIAN coord.
!
!        nops     number of operations in space group
!        mrot     rotation matrices in INTERNAL coordinates:
!                      ( 1,1 1,2 1,3 )( t_1 )   ( t_1' )
!                      ( 2,1 2,2 2,3 )( t_2 ) = ( t_2' )
!                      ( 3,1 3,2 3,3 )( t_3 )   ( t_3' )
!        tau      non-primitive translations, in INTERNAL coord.
!
!        ntype    number of atom types
!        neq      number of equivalent atoms of each type
!        ntyrep   representative atom for each type
!        pos      atomic positions in INTERNAL coord. (in fleur: taual)
!
!     the results for the lattice harmonics and local symmetry
!     information is put into the module mod_harmonics which can then
!     be use'd as needed.
!                                              m. weinert 12-99
!*********************************************************************
      CONTAINS
      SUBROUTINE local_sym(
     >                     lmaxd,lmax,nops,mrot,tau,
     >                     natd,ntype,neq,amat,bmat,pos,
     X                     nlhd,memd,ntypsd,l_dim,
     <                     nlhtyp,ntypsy,nlh,llh,nmem,mlh,clnu)

      USE m_ptsym
      USE m_lhcal
      USE m_constants, ONLY : pimach
      IMPLICIT NONE

!---> Arguments
      INTEGER, INTENT (IN) :: lmaxd,nops,ntype,natd
      INTEGER, INTENT (IN) :: neq(ntype),lmax(ntype),mrot(3,3,nops)
      REAL,    INTENT (IN) :: tau(3,nops),pos(3,natd)
      REAL,    INTENT (IN) :: amat(3,3),bmat(3,3)
      LOGICAL, INTENT (IN) :: l_dim
      INTEGER              :: nlhd,memd,ntypsd
      INTEGER              :: nlhtyp(ntype)
      INTEGER, INTENT(OUT) :: ntypsy(natd)
      INTEGER, INTENT(OUT) :: llh(0:nlhd,ntypsd),nmem(0:nlhd,ntypsd)
      INTEGER, INTENT(OUT) ::  mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      COMPLEX, INTENT(OUT) :: clnu(memd,0:nlhd,ntypsd)

!---> Locals
      INTEGER :: lmax0,mem_maxd,nlhd_max
      INTEGER :: lh,lm0,m,n,nsym,na,nsymt,nn
      REAL    :: orth(3,3,nops),amatinv(3,3)
      INTEGER :: nlhs(natd),locops(nops,natd),nrot(natd)
      INTEGER :: lnu((lmaxd+1)**2,natd)
      INTEGER :: mem((lmaxd+1)**2,natd)
      INTEGER :: lmnu(2*lmaxd+1,(lmaxd+1)**2,natd)
      COMPLEX :: c(2*lmaxd+1,(lmaxd+1)**2,natd)

      INTEGER, ALLOCATABLE :: typsym(:)

      amatinv = bmat / ( 2 * pimach() )
      mem_maxd = 2*lmaxd+1
      nlhd_max = (lmaxd+1)**2
      ALLOCATE ( typsym(natd) )

      WRITE (6,'(//," Local symmetries:",/,1x,17("-"))')
!
!===> determine the point group symmetries for each atom given
!===> the space group operations and atomic positions
!===> operations and positions are in internal (lattice) coordinates
!
      CALL ptsym(
     >           ntype,natd,neq,pos,nops,mrot,tau,lmax,
     <           nsymt,typsym,nrot,locops)

      WRITE (6,'("   symmetry kinds =",i4)') nsymt
      DO nsym = 1, nsymt
         WRITE (6,'(/,"   symmetry",i3,":",i4," operations in",
     &       " local point group",/,8x,"atoms:")') nsym,nrot(nsym)
         na = 0
         DO n=1,ntype
           DO nn = 1, neq(n)
             na = na + 1 
             IF ( typsym(na) == nsym ) WRITE (6,'(i14)') na
           ENDDO
         ENDDO
      ENDDO
!
!===>  generate the lattice harmonics for each local symmetry
!
      DO nsym = 1, nsymt

!--->    need to generate transformation matrices in cartesian
!--->    coordinates (rotations in real space)
         DO n = 1, nrot(nsym)
            orth(:,:,n) = matmul( amat,
     &         matmul( real( mrot(:,:,locops(n,nsym))),amatinv ) )
         ENDDO

!--->    get max. l for this symmetry type
         lmax0 = 0
         na = 0
         DO n=1,ntype
           DO nn = 1, neq(n)
             na = na + 1
             IF (typsym(na).EQ.nsym) lmax0 = max(lmax0,lmax(n))
           ENDDO
         ENDDO

!--->     generate the lattice harmonics
         CALL lhcal(
     >              mem_maxd,nlhd_max,lmax0,nrot(nsym),orth,
     <              nlhs(nsym),lnu(1,nsym),mem(1,nsym),
     <              lmnu(1,1,nsym),c(1,1,nsym))

      ENDDO
!
!====>  allocate arrays in module mod_harmonics and store for later use
!====>  this part can be changed depending on program to interface to;
!====>  this version is consistent with fleur.
!
      nlhd = 0
      memd = 0
      DO nsym = 1, nsymt
         nlhd = max(nlhd,nlhs(nsym))
         DO lh=1,nlhs(nsym)
            memd = max(memd,mem(lh,nsym))
         ENDDO
      ENDDO
      nlhd = nlhd - 1

      IF ( nsymt > ntypsd ) ntypsd = nsymt
      IF ( l_dim ) THEN
         DEALLOCATE ( typsym )
         RETURN
      ENDIF
      clnu = cmplx( 0.0,0.0 )
      mlh = 0
      DO nsym = 1,nsymt
         nlh(nsym) = nlhs(nsym)-1
         DO lh = 1, nlhs(nsym)
            llh(lh-1,nsym)  = lnu(lh,nsym)
            nmem(lh-1,nsym) = mem(lh,nsym)
            lm0 = lnu(lh,nsym)*(lnu(lh,nsym)+1) + 1
            DO m = 1, mem(lh,nsym)
               mlh(m,lh-1,nsym) = lmnu(m,lh,nsym) - lm0
               clnu(m,lh-1,nsym) = c(m,lh,nsym)
            ENDDO
         ENDDO
      ENDDO

      WHERE ( abs(aimag(clnu)) < 1.e-13 ) clnu = cmplx( real(clnu),0.0)
      WHERE ( abs( real(clnu)) < 1.e-13 ) clnu = cmplx(0.0,aimag(clnu))
!
!--->    different atom types may have the same symmetry, but different
!--->    lmax. to deal with this possibility, define nlhtyp(ntype) to
!--->    give the number of harmonics for each atom type.
!
      na = 0
      DO n = 1, ntype
         nlhtyp(n) = 0
         DO nn = 1,neq(n)
            na = na + 1
            DO lh = 1, nlh( typsym(na) )
               IF ( llh(lh,typsym(na)) .GT. lmax(n) ) EXIT
               nlhtyp(n) = nlhtyp(n) + 1
            ENDDO
         ENDDO
      ENDDO

      na = 0
      DO n = 1, ntype
         DO nn = 1,neq(n)
            na = na + 1
            ntypsy(na) = typsym(na)
!            ntypsy(na) = typsym(na-nn+1)
         ENDDO
      ENDDO

!---> output results
      DO n = 1, nsymt
        WRITE (6,'(/," --- Local symmetry",i3,":",i4,
     &       " lattice harmonics ",30("-"))') n,nlh(n)+1
        DO lh = 0,nlh(n)
          WRITE (6,'(/,5x,"lattice harmonic",i4,":  l=",i2,
     &         ",",i3," members:")') lh+1,llh(lh,n),nmem(lh,n)
          IF ( mod(nmem(lh,n),2)==1 ) THEN
            WRITE (6,'(5x,i5,2f14.8,5x,i5,2f14.8)')
     &                     mlh(1,lh,n),clnu(1,lh,n)
            IF ( nmem(lh,n) > 1 ) THEN
              WRITE (6,'(5x,i5,2f14.8,5x,i5,2f14.8)')
     &             (mlh(m,lh,n),clnu(m,lh,n),m=2,nmem(lh,n))
            ENDIF
          ELSE
            WRITE (6,'(5x,i5,2f14.8,5x,i5,2f14.8)')
     &            (mlh(m,lh,n),clnu(m,lh,n),m=1,nmem(lh,n))
          ENDIF
        ENDDO
      ENDDO

      DEALLOCATE ( typsym )
      RETURN
      END SUBROUTINE local_sym
      END MODULE m_localsym
