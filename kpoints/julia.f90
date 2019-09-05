!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_julia

USE m_juDFT

CONTAINS

SUBROUTINE julia(sym,cell,input,noco,banddos,kpts,l_q,l_fillArrays)
!----------------------------------------------------------------------+
! Generate a k-point file with approx. nkpt k-pts or a Monkhorst-Pack  |
! set with nmod(i) divisions in i=x,y,z direction. Interface to kptmop |
! and kpttet routines of the MD-programm.                              |
!                                                          G.B. 07/01  |
!----------------------------------------------------------------------+

   USE m_constants
   USE m_bravais
   USE m_divi
   USE m_brzone
   USE m_brzone2
   USE m_kptmop
   USE m_kpttet
   USE m_bandstr1
   USE m_types

   IMPLICIT NONE

   TYPE(t_sym),     INTENT(IN)    :: sym
   TYPE(t_cell),    INTENT(IN)    :: cell
   TYPE(t_input),   INTENT(IN)    :: input
   TYPE(t_noco),    INTENT(IN)    :: noco
   TYPE(t_banddos), INTENT(IN)    :: banddos
   TYPE(t_kpts),    INTENT(INOUT) :: kpts

   LOGICAL, INTENT (IN)           :: l_q, l_fillArrays
   
   INTEGER, PARAMETER :: nop48  = 48
   INTEGER, PARAMETER :: mface  = 51
   INTEGER, PARAMETER :: mdir   = 10
   INTEGER, PARAMETER :: nbsz   =  3
   INTEGER, PARAMETER :: ibfile =  6
   INTEGER, PARAMETER :: nv48   = (2*nbsz+1)**3+48

   INTEGER ndiv3              ! max. number of tetrahedrons (< 6*(kpts%nkpt+1)
   INTEGER ntet               ! actual number of tetrahedrons

   REAL, ALLOCATABLE    :: vkxyz(:,:)  ! vector of kpoint generated; in cartesian representation
   REAL, ALLOCATABLE    :: wghtkp(:)   !   associated with k-points for BZ integration
   INTEGER, ALLOCATABLE :: ntetra(:,:) ! corners of the tetrahedrons
   REAL, ALLOCATABLE    :: voltet(:)   ! voulmes of the tetrahedrons
   REAL, ALLOCATABLE    :: vktet(:,:)

   REAL    divis(4)           ! Used to find more accurate representation of k-points
                              ! vklmn(i,kpt)/divis(i) and weights as wght(kpt)/divis(4)
   INTEGER nkstar             ! number of stars for k-points generated in full stars
   REAL    bltv(3,3)          ! cartesian Bravais lattice basis (a.u.)
   REAL    rltv(3,3)          ! reciprocal lattice basis (2\pi/a.u.)
   REAL    ccr(3,3,nop48)     ! rotation matrices in cartesian repr.
   REAL    rlsymr(3,3,nop48)  ! rotation matrices in reciprocal lattice basis representation
   REAL    talfa(3,nop48)     ! translation vector associated with (non-symmorphic)
                              ! symmetry elements in Bravais lattice representation
   INTEGER ncorn,nedge,nface  ! number of corners, faces and edges of the IBZ
   REAL    fnorm(3,mface)     ! normal vector of the planes bordering the IBZ
   REAL    fdist(mface)       ! distance vector of the planes bordering t IBZ
   REAL    cpoint(3,mface)    ! cartesian coordinates of corner points of IBZ
   REAL    xvec(3)            ! arbitrary vector lying in the IBZ

   INTEGER idsyst   ! crystal system identification in MDDFT programs
   INTEGER idtype   ! lattice type identification in MDDFT programs

   INTEGER idimens  ! number of dimensions for k-point set (2 or 3)
   INTEGER nreg     ! 1 kpoints in full BZ; 0 kpoints in irrBZ
   INTEGER nfulst   ! 1 kpoints ordered in full stars
                    !    (meaningful only for nreg =1; full BZ)
   INTEGER nbound   ! 0 no primary points on BZ boundary;
                    ! 1 with boundary points (not for BZ integration!!!)
   INTEGER ikzero   ! 0 no shift of k-points;
                    ! 1 shift of k-points for better use of sym in irrBZ
   REAL    kzero(3) ! shifting vector to bring one k-point to or 
                    ! away from (0,0,0) (for even/odd nkpt3)

   INTEGER i,j,k,l,idiv,mkpt,addSym,nsym
   INTEGER iofile,iokpt,kpri,ktest,kmidtet
   INTEGER idivis(3)
   LOGICAL random,trias
   REAL help(3),binv(3,3),rlsymr1(3,3),ccr1(3,3)

   random = .false.  ! do not use random tetra-points

   !------------------------------------------------------------
   !
   !        idsyst         idtype 
   !
   !   1  cubic          primitive
   !   2  tetragonal     body centered
   !   3  orthorhombic   face centered
   !   4  hexagonal      A-face centered
   !   5  trigonal       B-face centered
   !   6  monoclinic     C-face centered
   !   7  triclinic 
   !
   ! --->   for 2 dimensions only the following Bravais lattices exist:
   !
   !    TYPE                    EQUIVALENT 3-DIM        idsyst/idtype
   !   square               = p-tetragonal ( 1+2 axis )      2/1
   !   rectangular          = p-orthorhomb ( 1+2 axis )      3/1
   !   centered rectangular = c-face-orthorhomb( 1+2 axis)   3/6
   !   hexagonal            = p-hexagonal  ( 1+2 axis )      4/1
   !   oblique              = p-monoclinic ( 1+2 axis )      6/1
   !
   !------------------------------------------------------------

   IF(l_q) THEN
      trias=input%tria
      if (input%tria) call judft_error("tria=T not implemented for q-point generator",calledby='julia')
      !input%tria=.false.
   ENDIF
       
   IF (cell%latnam.EQ.'squ') THEN
      idsyst = 2
      idtype = 1
         IF (.not.input%film) THEN
            IF (abs(cell%amat(1,1)-cell%amat(3,3)) < 0.0000001) THEN
               idsyst = 1
               idtype = 1
            END IF
         END IF
      END IF
   IF (cell%latnam.EQ.'p-r') THEN
      idsyst = 3
      idtype = 1
   END IF
   IF ((cell%latnam.EQ.'c-b').OR.(cell%latnam.EQ.'c-r')) THEN
      idsyst = 3
      idtype = 6
   END IF
   IF ((cell%latnam.EQ.'hex').OR.(cell%latnam.EQ.'hx3')) THEN
      idsyst = 4
      idtype = 1
   END IF
   IF (cell%latnam.EQ.'obl') THEN
      idsyst = 6
      idtype = 1
   END IF
   IF (cell%latnam.EQ.'any') THEN
      CALL bravais(cell%amat,idsyst,idtype) 
   END IF
   nsym = sym%nop
   IF (input%film) nsym = sym%nop2        

   ! Want to make a Bandstructure?
   IF (banddos%ndir == -4) THEN
      CALL bandstr1(idsyst,idtype,cell%bmat,kpts,input,l_fillArrays,banddos)
      RETURN
   END IF

   ! Some variables we do not use

   iofile = 6
   iokpt  = 6
   kpri   = 0 ! 3
   ktest  = 0 ! 5
   kmidtet = 0
   nreg    = 0
   nfulst  = 0
   ikzero  = 0
   kzero(1) = 0.0
   kzero(2) = 0.0
   kzero(3) = 0.0 
   nbound  = 0
   IF (input%tria) THEN
      IF (input%film) nbound  = 1
      ! IF ((idsyst==1).AND.(idtype==1)) nbound  = 1
      ! IF ((idsyst==2).AND.(idtype==1)) nbound  = 1
      ! IF ((idsyst==3).AND.(idtype==1)) nbound  = 1
      ! IF ((idsyst==3).AND.(idtype==6)) nbound  = 1
      ! IF ((idsyst==4).AND.(idtype==1)) nbound  = 1
      IF (nbound == 0) random = .true.
   END IF
   idimens = 3
   IF (input%film) idimens = 2

   ! Lattice information

   DO j = 1, 3
      DO k = 1, 3
         bltv(j,k) = cell%amat(k,j)
         binv(j,k) = cell%bmat(k,j) / tpi_const
         rltv(j,k) = cell%bmat(k,j)
         DO i = 1,nsym
            rlsymr(k,j,i) = real(sym%mrot(j,k,i))
         END DO
      END DO
   END DO

   ccr = 0.0
   DO i = 1, nsym
      DO j = 1, 3
         talfa(j,i) = 0.0
         DO k = 1, 3
            talfa(j,i) = bltv(j,k) * sym%tau(k,i)
            help(k) = 0.0
            DO l = 1, 3
               help(k) =  help(k) + rlsymr(l,k,i) * binv(j,l)
            END DO
         END DO
         DO k = 1, 3
            ccr(j,k,i) = 0.0
            DO l = 1, 3
               ccr(j,k,i) = ccr(j,k,i) + bltv(l,k) * help(l)
            END DO
         END DO
      END DO
   END DO
   DO i = 1, nsym
      rlsymr1(:,:) = rlsymr(:,:,i)
      ccr1(:,:)    = ccr(:,:,i)
      DO j = 1, 3
         DO k = 1, 3
            rlsymr(k,j,i) = rlsymr1(j,k)
            ccr(k,j,i)    = ccr1(j,k)
         END DO
      END DO
   END DO

   IF ((.not.noco%l_ss).AND.(.not.noco%l_soc).AND.(2*nsym<nop48)) THEN
      IF ((input%film.AND.(.not.sym%invs2)).OR.((.not.input%film).AND.(.not.sym%invs))) THEN
         addSym = 0
         ! Note: We have to add the negative of each symmetry operation
         !       to exploit time reversal symmetry. However, if the new
         !       symmetry operation is the identity matrix it is excluded.
         !       This is the case iff it is (-Id) + a translation vector.
         DO i = 1, nsym
            ! This test assumes that ccr(:,:,1) is the identity matrix.
            IF(.NOT.ALL(ABS(ccr(:,:,1)+ccr(:,:,i)).LT.10e-10) ) THEN
               ccr(:,:,nsym+addSym+1 ) = -ccr(:,:,i)
               rlsymr(:,:,nsym+addSym+1 ) = -rlsymr(:,:,i)
               addSym = addSym + 1
            END IF
         END DO
         nsym = nsym + addSym
      END IF
   END IF

   ! brzone and brzone2 find the corner-points, the edges, and the
   ! faces of the irreducible wedge of the brillouin zone (IBZ).
   ! In these subroutines many special cases can occur. Due to this the very 
   ! sophisticated old routine brzone had a few bugs. The new routine
   ! brzone2 was written with a different algorithm that is slightly slower
   ! but should be more stable. To make comparisons possible the old
   ! routine is only commented out. Both routines are directly 
   ! interchangable. GM, 2016.

!   CALL brzone(rltv,nsym,ccr,mface,nbsz,nv48,cpoint,xvec,ncorn,nedge,nface,fnorm,fdist)

   CALL brzone2(rltv,nsym,ccr,mface,nbsz,nv48,cpoint,xvec,ncorn,nedge,nface,fnorm,fdist)

   IF (input%tria.AND.random) THEN
      ! Calculate the points for tetrahedron method
      mkpt = kpts%nkpt
      ndiv3 = 6*(mkpt+1)
      ALLOCATE (vkxyz(3,mkpt),wghtkp(mkpt))
      ALLOCATE (voltet(ndiv3),vktet(3,mkpt),ntetra(4,ndiv3))
      vkxyz = 0.0
      CALL kpttet(iofile,ibfile,iokpt,kpri,ktest,kmidtet,mkpt,ndiv3,&
                  nreg,nfulst,rltv,cell%omtil,nsym,ccr,mdir,mface,&
                  ncorn,nface,fdist,fnorm,cpoint,voltet,ntetra,ntet,vktet,&
                  kpts%nkpt,divis,vkxyz,wghtkp)
   ELSE
      ! If just the total number of k-points is given, determine 
      ! the divisions in each direction (nkpt3):

      ! IF (tria) THEN
      !    nkpt = nkpt/4
      !    nkpt3(:) = nkpt3(:) / 2
      ! END IF
      IF (sum(kpts%nkpt3).EQ.0) THEN
         CALL divi(kpts%nkpt,cell%bmat,input%film,sym%nop,sym%nop2,kpts%nkpt3)
      END IF

      ! Now calculate Monkhorst-Pack k-points:
      IF (kpts%nkpt3(2).EQ.0) kpts%nkpt3(2) = kpts%nkpt3(1)
      IF ((.not.input%film).AND.(kpts%nkpt3(3).EQ.0)) kpts%nkpt3(3) = kpts%nkpt3(2)
      IF (nbound.EQ.1) THEN
         mkpt = (2*kpts%nkpt3(1)+1)*(2*kpts%nkpt3(2)+1)
         IF (.not.input%film) mkpt = mkpt*(2*kpts%nkpt3(3)+1)
      ELSE
         mkpt = kpts%nkpt3(1)*kpts%nkpt3(2)
         IF (.not.input%film) mkpt = mkpt*kpts%nkpt3(3)
      END IF
      ALLOCATE (vkxyz(3,mkpt),wghtkp(mkpt) )
      vkxyz = 0.0

      CALL kptmop(iofile,iokpt,kpri,ktest,idsyst,idtype,kpts%nkpt3,ikzero,kzero,&
                  rltv,bltv,nreg,nfulst,nbound,idimens,xvec,fnorm,fdist,ncorn,nface,&
                  nedge,cpoint,nsym,ccr,rlsymr,talfa,mkpt,mface,mdir,&
                  kpts%nkpt,divis,vkxyz,nkstar,wghtkp)
   END IF

   idivis(1) = int(divis(1)) 
   idivis(2) = int(divis(2)) 
   idivis(3) = int(divis(3)) 
   idiv = lcm(3,idivis)
   IF (idiv.GE.200) idiv = 1
   DO j=1,kpts%nkpt
      wghtkp(j) = wghtkp(j) * divis(4)
      DO k = 1,3
         help(k) = 0.0
         DO l = 1,3
            help(k) = help(k) + cell%amat(l,k) * vkxyz(l,j)
         END DO
      END DO
      DO i=1,3
         vkxyz(i,j) = help(i) * idiv / tpi_const
      END DO
   END DO

   ! if (l_q) write qpts file:
   IF(l_q)THEN
      IF(input%film) THEN
         CALL juDFT_error("For the case of input%film q-points generator not implemented!", calledby = "julia")
      END IF
    
      OPEN(113,file='qpts',form='formatted',status='new')
      WRITE(113,'(i5)') kpts%nkpt+1
      WRITE(113,8050) 0.,0.,0.
      DO j = 1, kpts%nkpt
         WRITE (113,FMT=8050) (vkxyz(i,j)/real(idiv),i=1,3)
      END DO
      CLOSE(113)
         !input%tria=trias
         RETURN
   END IF
   8050 FORMAT (2(f14.10,1x),f14.10)

   ! write k-points file or write data into arrays
   IF (l_fillArrays) THEN
      IF (ALLOCATED(kpts%bk)) THEN
         DEALLOCATE(kpts%bk)
      END IF
      IF (ALLOCATED(kpts%wtkpt)) THEN
         DEALLOCATE(kpts%wtkpt)
      END IF
      ALLOCATE(kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt))
      IF (idiv.NE.0) kpts%posScale = REAL(idiv)
      DO j = 1, kpts%nkpt
         kpts%bk(1,j) = vkxyz(1,j)
         kpts%bk(2,j) = vkxyz(2,j)
         kpts%bk(3,j) = vkxyz(3,j)
         kpts%wtkpt(j) = wghtkp(j)
      END DO
      IF (input%tria.AND.random) THEN
         kpts%ntet = ntet
         IF (ALLOCATED(kpts%ntetra)) THEN
            DEALLOCATE(kpts%ntetra)
         END IF
         IF (ALLOCATED(kpts%voltet)) THEN
            DEALLOCATE(kpts%voltet)
         END IF
         ALLOCATE(kpts%ntetra(4,kpts%ntet))
         ALLOCATE(kpts%voltet(kpts%ntet))
         DO j = 1, ntet
            DO i = 1, 4
               kpts%ntetra(i,j) = ntetra(i,j)
            END DO
            kpts%voltet(j) = ABS(voltet(j))
         END DO
      END IF
   ELSE
      OPEN (41,file='kpts',form='formatted',status='new')
      IF (input%film) THEN
         WRITE (41,FMT=8110) kpts%nkpt,real(idiv),.false.
         DO j = kpts%nkpt, 1, -1
            WRITE (41,FMT=8040) (vkxyz(i,j),i=1,2),wghtkp(j)
         END DO
      ELSE
         WRITE (41,FMT=8100) kpts%nkpt,real(idiv)
         DO j = 1, kpts%nkpt
            WRITE (41,FMT=8040) (vkxyz(i,j),i=1,3),wghtkp(j)
         END DO
         IF (input%tria.AND.random) THEN
            WRITE (41,'(i5)') ntet
            WRITE (41,'(4(4i6,4x))') ((ntetra(i,j),i=1,4),j=1,ntet)
            WRITE (41,'(4f20.13)') (ABS(voltet(j)),j=1,ntet)
         END IF
      END IF
      8100 FORMAT (i5,f20.10)
      8110 FORMAT (i5,f20.10,3x,l1)
      8040 FORMAT (4f10.5)
      CLOSE (41)
   END IF

   DEALLOCATE (vkxyz,wghtkp)
   IF (input%tria.AND..not.input%film) DEALLOCATE (voltet,vktet,ntetra)
   RETURN

   CONTAINS

   INTEGER FUNCTION lcm( n, ints )
      ! Compute least common multiple (lcm) of n positive integers.
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: ints(n)

      INTEGER :: i,j,m

      IF (any(ints(1:n) <= 0)) THEN
         m = 0
      ELSE
         m = maxval( ints(1:n) )
         DO i = 1, n
            DO j = 1, ints(i) / 2
               IF (mod(m*j,ints(i)) == 0) EXIT
            END DO
            m = m*j
         END DO
      END IF
      lcm = m
      RETURN
   END FUNCTION lcm

END SUBROUTINE julia

END MODULE m_julia
