!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_hsdata
!-----------------------------------------------
! DESC:Module to store, handle and communicate the FLEUR Hamiltonian
!
!      In the port to current FLEUR the Hamiltonian and overlap are
!      kept as FULL Hermitian complex matrices (the modern hsmt/hs_int
!      fill the upper triangle data(j,i)=<G_j|M|G_i>, j<=i; the lower
!      triangle is completed in gf_storeHS). The old packed storage
!      is gone.
!                 Daniel Wortmann, (05-08-24)
!-----------------------------------------------
      use m_juDFT
      USE m_gf_types, ONLY: t_lapw, t_lapw_gf
      USE m_gf_energies, ONLY: gf_z
      IMPLICIT NONE
      PRIVATE
      COMPLEX,ALLOCATABLE,TARGET,SAVE :: H(:,:),S(:,:)
      INTEGER,SAVE                    :: spin,k,stored_layer
      PUBLIC gf_storeHS,gf_getlargeH_eS,gf_writeHS,gf_readHS            &
     &     ,gf_getHSaddr,gf_hs_map
      CONTAINS
      !<-- S: gf_getHSaddr(H,S)
      SUBROUTINE gf_getHSaddr(layer,HP,SP)
!-----------------------------------------------
!     get address of H and S
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN)     :: layer
      COMPLEX,POINTER:: hp(:,:),sp(:,:)

      IF ((layer /= stored_layer).OR..NOT.allocated(h)) THEN
         IF (.NOT.gf_readHS(layer,spin,k)) THEN
             write(*,*) "Stored layer:",stored_layer
             write(*,*) "Required:",layer
             CALL juDFT_error(              &
     &        "Did not find the Data of the FLAPW-Hamiltonian")
         endif
      ENDIF
      IF(.NOT.allocated(h)) CALL juDFT_error("gf_getHSaddr",calledby="gf_hsdata")
      IF(.NOT.allocated(s)) CALL juDFT_error("gf_getHSaddr",calledby="gf_hsdata")
      hp=>h
      sp=>s
      END SUBROUTINE

      !<-- S: gf_storeHS(jspin,nk,layer,hmat,smat)
      SUBROUTINE gf_storeHS(jspin,nk,layer,hmat,smat,perm)
!-----------------------------------------------
!     stores the assembled H and S of one (layer,spin,k). The input
!     matrices have only their upper triangle filled (modern hsmt
!     convention); the Hermitian completion happens here. perm maps
!     the input basis index onto the stored one (noco reordering).
!-----------------------------------------------
      IMPLICIT NONE
      !<-- Arguments
      INTEGER,INTENT(IN)  :: jspin,nk,layer
      COMPLEX,INTENT(IN)  :: hmat(:,:),smat(:,:)
      INTEGER,INTENT(IN),OPTIONAL :: perm(:)
      !>
      INTEGER             :: n1,n2,nbas,err,i1,i2
      INTEGER,ALLOCATABLE :: p(:)
      nbas = SIZE(hmat,1)
      IF (allocated(H)) THEN
         IF (SIZE(H,1) /= nbas)  DEALLOCATE(H,s)
      ENDIF
      IF (.NOT.allocated(H)) THEN
         ALLOCATE(H(nbas,nbas),S(nbas,nbas),STAT=err)
         IF (err /= 0) CALL juDFT_error                                   &
     &        ("Not enough memory for Hamiltonian")
      ENDIF
      spin            = jspin
      k               = nk
      stored_layer    = layer
      IF (PRESENT(perm)) THEN
         p = perm
      ELSE
         p = [(n1,n1=1,nbas)]
      ENDIF
      DO n1 = 1,nbas
         i1 = p(n1)
         DO n2 = 1,n1
            i2 = p(n2)
            H(i2,i1) = hmat(n2,n1)
            S(i2,i1) = smat(n2,n1)
            H(i1,i2) = CONJG(hmat(n2,n1))
            S(i1,i2) = CONJG(smat(n2,n1))
         ENDDO
      ENDDO
      END SUBROUTINE

      !>
      !<-- S: gf_getlargeH_eS(layer,jspin,nk,en,lapw,lapw_gf,l_sph,Hinv)
      SUBROUTINE gf_getlargeH_eS(layer,jspin,nk,en,lapw,lapw_gf,l_sph,  &
     &     Hinv)
!-----------------------------------------------
!     Constructs the large 2D G^-1 = (e*S-H) matrix, in the full basis
!     or (l_sph) restricted to the rkmax sphere plus the LOs
!           (last modified: 05-08-24) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<-- Arguments
      INTEGER,INTENT(IN)     :: layer,jspin,nk,en
      TYPE(t_lapw),INTENT(IN)    :: lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      LOGICAL,INTENT(IN)     :: l_sph
      COMPLEX,INTENT(OUT)    :: Hinv(:,:)
      !>
      !<-- Locals
      INTEGER             :: i,j
      INTEGER             :: map(lapw%nmat)
      COMPLEX             :: z
      !>

      !<-- Check if correct h,s in memory, otherwise try to read it

      IF (.NOT.allocated(h).OR..NOT.(layer == stored_layer.AND.spin == &
     &     jspin.AND.k == nk)) THEN
         IF (.NOT.gf_readHS(layer,jspin,nk)) CALL juDFT_error             &
     &        ("Did not find the Data of the FLAPW-Hamiltonian")
      ENDIF
      !>
      !<--Check if enough data is present
      IF (lapw%nmat>SIZE(h,1)) CALL                                     &
     &     juDFT_error("Not enough data in FLAPW-Hamiltonian")
      !>
                       !The energy
      z=gf_z(en,layer)
      map = gf_hs_map(jspin,lapw,lapw_gf,l_sph)
      IF (SIZE(Hinv,1) /= MAXVAL(map)) CALL                             &
     &     juDFT_error("Wrong size of G in gf_getlargeH_eS")
      DO i = 1,lapw%nmat
         IF (map(i)==0) CYCLE
         DO j = 1,lapw%nmat
            IF (map(j)==0) CYCLE
            Hinv(map(j),map(i)) = z*S(j,i) - H(j,i)
         ENDDO
      ENDDO
      END SUBROUTINE
      !>
      !<-- F: gf_hs_map(jspin,lapw,lapw_gf,l_sph)
      FUNCTION gf_hs_map(jspin,lapw,lapw_gf,l_sph) RESULT(map)
!-----------------------------------------------
!     index map from the stored H/S onto the working matrix: identity,
!     or (l_sph) the cylinder basis functions outside the rkmax sphere
!     are compressed out (map=0) and the LOs follow the sphere
!-----------------------------------------------
      IMPLICIT NONE
      INTEGER,INTENT(IN)         :: jspin
      TYPE(t_lapw),INTENT(IN)    :: lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
      LOGICAL,INTENT(IN)         :: l_sph
      INTEGER                    :: map(lapw%nmat)

      INTEGER :: i,nv_lo,nv_skip,n,nn

      IF (.NOT.l_sph) THEN
         map = (/(i,i=1,lapw%nmat)/)
         RETURN
      ENDIF
      nv_skip=lapw%nv(jspin)-lapw_gf%nv_sphere(jspin)  !no of lapw to be skipped
      if (lapw%nmat>lapw%nv_tot) THEN
         !There are some LO
         nv_lo=lapw%nmat-lapw%nv_tot
      else
         nv_lo=0
      endif
      if (lapw%nv_tot>lapw%nv(jspin)) nv_lo=nv_lo/2  !noco calculation
      n=0
      nn=0

      map(n+1:n+lapw_gf%nv_sphere(jspin))=(/(i,i=nn+1,nn+lapw_gf%nv_sphere(jspin))/)
      n=n+lapw_gf%nv_sphere(jspin)
      nn=nn+lapw_gf%nv_sphere(jspin)
      !skip some lapw
      if (nv_skip>0) THEN
         map(n+1:n+nv_skip)=0
         n=n+nv_skip
      endif
      if (lapw%nv_tot>lapw%nv(1)) THEN
         !noco
         map(n+1:n+lapw_gf%nv_sphere(jspin))=(/(i,i=nn+1,nn+lapw_gf%nv_sphere(jspin))/)
         n=n+lapw_gf%nv_sphere(jspin)
         nn=nn+lapw_gf%nv_sphere(jspin)
         !skip some lapw
         if (nv_skip>0) THEN
            map(n+1:n+nv_skip)=0
            n=n+nv_skip
         endif
      endif
      !add the LO
      if (nv_lo>0) THEN
         map(n+1:n+nv_lo)=(/(i,i=nn+1,nn+nv_lo)/)
         n=n+nv_lo
         nn=nn+nv_lo
         if (lapw%nv_tot>lapw%nv(1)) THEN
            map(n+1:n+nv_lo)=(/(i,i=nn+1,nn+nv_lo)/)
            n=n+nv_lo
            nn=nn+nv_lo
         endif
      endif
      END FUNCTION
      !>

      !<-- S: gf_writeHS(savemem)
      SUBROUTINE gf_writeHS(savemem)
!-----------------------------------------------
!      writes the hamiltonian and Overlapp matrix to a unformated file
!      The filename will contain the kpt&spin information
!      so that all the different H and S used will be kept in the
!      current directory.
!           (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      LOGICAL,INTENT(IN)     :: savemem
      !<-- Locals
      !>
      OPEN(111,FILE = priv_makefilename(stored_layer,spin,k),FORM='unformatted',STATUS  ='REPLACE')
      WRITE(111) SIZE(h,1),stored_layer,spin,k
      WRITE(111) h
      WRITE(111) s
      CLOSE(111)
      IF (savemem) THEN
         DEALLOCATE(h,s)
      ENDIF
      END SUBROUTINE
      !>
      !<-- F: gf_readHS(layer,jspin,nk)
      FUNCTION gf_readHS(layer,jspin,nk)
!-----------------------------------------------
!    reads the H and S matrices from the file
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)     :: layer,jspin,nk
      LOGICAL                :: gf_readHS
      !>
      !<-- Locals
      INTEGER             :: n
      !>
      OPEN(111,FILE = priv_makefilename(layer,jspin,nk),FORM='unformatted',STATUS  ='OLD',ERR = 100)
      READ(111,ERR = 100) n,stored_layer,spin,k
      IF (spin/=jspin.OR.k/=nk.OR.layer/=stored_layer) CALL       &
     &     juDFT_error("Wrong data in FLEUR-Hamiltonian file")
      !<-- The storage might have the wrong size or might not be allocated
      IF (.NOT.allocated(H)) THEN
         ALLOCATE(H(n,n),S(n,n))
      ELSEIF (n /= SIZE(h,1)) THEN
         DEALLOCATE(H,S)
         ALLOCATE(H(n,n),S(n,n))
      ENDIF
      !>
      READ(111,ERR = 100) h
      READ(111,ERR = 100) s
      CLOSE(111,ERR = 100)
      gf_readHS=.TRUE.
      RETURN
                        !something went wrong
  100 gf_readHS=.FALSE.

      END FUNCTION
      !>
      !<-- F: priv_makefilename(layer,jspin,nk)
      FUNCTION priv_makefilename(layer,jspin,nk)
!-----------------------------------------------
!   generates a filename for the current jspin,nk
!             (last modified: 2004-00-00) D. Wortmann
!-----------------------------------------------
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)     :: layer,nk,jspin
      CHARACTER(LEN = 15)       :: priv_makefilename
      !>
#ifdef CPP_MPI
      WRITE(priv_makefilename,"(a,i4.0,i1,i4)") "HSmat",layer,jspin,nk
#else
      WRITE(priv_makefilename,"(a,i4.0,i1,i4)") "HSmat",layer
#endif
      END FUNCTION
      !>
      END
