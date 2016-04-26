      MODULE m_socorssdw
      use m_juDFT
      CONTAINS
      SUBROUTINE soc_or_ssdw(
     >                       l_soc,l_ss,theta,phi,qss,amat,
     >                       mrot,tau,nop,nop2,nat,atomid,atompos,
     <                       mmrot,ttr,no3,no2,ntype,neq,natmap,
     <                       ntyrep,natype,natrep,zatom,pos)

      USE m_sssym
      USE m_socsym
      IMPLICIT NONE

      LOGICAL, INTENT (IN) :: l_soc,l_ss
      INTEGER, INTENT (IN) :: nop,nop2,nat
      REAL,    INTENT (IN) :: theta,phi
      INTEGER, INTENT (OUT):: no3,no2,ntype

      INTEGER, INTENT (IN) :: mrot(3,3,nop)
      REAL,    INTENT (IN) :: tau(3,nop),qss(3),amat(3,3)
      REAL,    INTENT (IN) :: atomid(nat),atompos(3,nat)
      INTEGER, INTENT (OUT):: mmrot(3,3,nop)
      REAL,    INTENT (OUT):: ttr(3,nop)
!--> actually, intent out:
      INTEGER, ALLOCATABLE :: neq(:), ntyrep(:)              ! these variables are allocated with
      REAL,    ALLOCATABLE :: zatom(:)                       ! dim 'ntype'
      INTEGER, ALLOCATABLE :: natype(:),natrep(:),natmap(:)  ! or  'nat'
      REAL,    ALLOCATABLE :: pos(:,:)                       ! or  '3,nat'

      INTEGER n,nt,i,j,nops,ntypm,ity(nat)
      REAL    tr(3),eps7
      LOGICAL lnew
      LOGICAL, ALLOCATABLE :: error(:)

      eps7 = 1.0e-7 

      ALLOCATE ( error(nop) )
      error(:) = .false.
      IF (l_ss) THEN                     ! reduce symmetry if SSDW calculation
        CALL ss_sym(
     >              nop,mrot,qss,
     <              error)
      ENDIF
      IF (l_soc) THEN                    ! reduce symmetry if SOC calculation
        CALL soc_sym(
     >               nop,mrot,theta,phi,amat,
     <               error)
      ENDIF
      IF (l_ss.AND.l_soc)  CALL juDFT_error("no spin-spirals with SOC!"
     +     ,calledby ="soc_or_ssdw")
      no2 = 0                      ! No. of 2D sym.op's allowed by SOC or SS
      DO n = 1, nop2
        IF ( .not.error(n) ) THEN
           no2 = no2 + 1
           mmrot(:,:,no2) = mrot(:,:,n)
           ttr(:,no2) = tau(:,n)
        ENDIF
      ENDDO
      no3 = no2                    ! same for 3D sym.op's
      DO n = nop2+1,nop
        IF ( .not.error(n) ) THEN
           no3 = no3 + 1
           mmrot(:,:,no3) = mrot(:,:,n)
           ttr(:,no3) = tau(:,n)
        ENDIF
      ENDDO
      DEALLOCATE (error)

!---> determine the number of distinct atoms based on atomic number,
!---> etc. (not necessarily symmetry inequivalent)

      ntypm = 1
      ity(1) = 1
      DO n=2, nat
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

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!--->  at this point, the symmetry is correct (assumed here)

      nops = no3
      natype(1:nat) = 0
      ntype = 0
      DO i =1,nat
         IF ( natype(i) .ne. 0 ) cycle
         ntype = ntype + 1   ! new atom type
         natype(i) = ntype   ! atom type
         natrep(i) = i       ! this is the representative
!--->    rotate representative and get symmetry equavalent atoms
         DO n=1,nops
            tr(:) = matmul( mmrot(:,:,n) , pos(:,i) ) + ttr(:,n)
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
         zatom( natype(n) ) = NINT(atomid(n))
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


      END SUBROUTINE soc_or_ssdw
      END MODULE m_socorssdw
