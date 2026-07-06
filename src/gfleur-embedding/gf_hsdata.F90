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
      IMPLICIT NONE
      PRIVATE
      COMPLEX,ALLOCATABLE,TARGET,SAVE :: H(:,:),S(:,:)
      INTEGER,SAVE                    :: spin,k,stored_layer
      PUBLIC gf_storeHS,gf_getlargeH_eS,gf_writeHS,gf_readHS            &
     &     ,gf_getHSaddr
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
      SUBROUTINE gf_storeHS(jspin,nk,layer,hmat,smat)
!-----------------------------------------------
!     stores the assembled H and S of one (layer,spin,k). The input
!     matrices have only their upper triangle filled (modern hsmt
!     convention); the Hermitian completion happens here.
!-----------------------------------------------
      IMPLICIT NONE
      !<-- Arguments
      INTEGER,INTENT(IN)  :: jspin,nk,layer
      COMPLEX,INTENT(IN)  :: hmat(:,:),smat(:,:)
      !>
      INTEGER             :: n1,n2,nbas,err
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
      DO n1 = 1,nbas
         DO n2 = 1,n1
            H(n2,n1) = hmat(n2,n1)
            S(n2,n1) = smat(n2,n1)
            H(n1,n2) = CONJG(hmat(n2,n1))
            S(n1,n2) = CONJG(smat(n2,n1))
         ENDDO
      ENDDO
      END SUBROUTINE

      !>
      !<-- S: gf_getlargeH_eS(jspin,nk,H,S)
      SUBROUTINE gf_getlargeH_eS(layer,jspin,nk,en,Hinv)
!-----------------------------------------------
!     Constructs the large 2D G^-1 = (e*S-H) matrix
!           (last modified: 05-08-24) D. Wortmann
!-----------------------------------------------
      USE m_gf_energies,ONLY:gf_z
      IMPLICIT NONE
      !<-- Arguments
      INTEGER,INTENT(IN)     :: layer,jspin,nk,en
      COMPLEX,INTENT(OUT)    :: Hinv(:,:)
      !>
      !<-- Locals
      INTEGER             :: n
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
      IF (SIZE(hinv,1)>SIZE(h,1)) CALL                                  &
     &     juDFT_error("Not enough data in FLAPW-Hamiltonian")
      !>
                       !The energy
      z=gf_z(en,layer)
      n = SIZE(Hinv,1)
      Hinv(:,:) = z*S(:n,:n) - H(:n,:n)
      END SUBROUTINE
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
