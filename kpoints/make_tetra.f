      MODULE m_maketetra
      CONTAINS
      SUBROUTINE make_tetra(
     >     nkpt,bk,ntria,itria,atr,
     <     ntetra,itetra,voltet)
C----------------------------------------------------------------
c     
c     Make tetrahedrons out of triangles; assume that layers of
c     k-points exist, that have been grouped to triangles in
c     subroutine fertri. Build prisms from them and cut each into
c     three tetrahedrons.
c     
c     bk(1-3,nk)    ... coordinates of k-point nk , nk = 1,nkpt
c     ntria         ... number of triangles per layer
c     itria(1-3,nt) ... index of k-points forming triangle nt
c     atr(nt)       ... area of triangle nt
c     
c     ntetra)       ... number of tetrahedrons
c     itetra(1-4,nt)... index of k-points forming tetrahedron nt
c     voltet(nt)    ... volume of tetrahedron nt
c     omega_bz      ... volume of irreducible part of BZ
c     
c     Note that the following assumes layers of equally distributed k-points !
c     
C----------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nkpt,ntria
      INTEGER, INTENT (OUT) :: ntetra

      INTEGER, INTENT (IN)  :: itria(:,:) !(3,2*nkptd)
      REAL,    INTENT (IN)  :: atr(:) !(2*nkptd)
      REAL,    INTENT (IN)  :: bk(:,:) !(3,nkptd)
      INTEGER, INTENT (OUT) :: itetra(:,:) !(4,6*nkptd)
      REAL,    INTENT (OUT) :: voltet(:) !(6*nkpt)


      INTEGER ikpt,nkpp,itri,itet,ip1,ip2,ip3,ip4,ip5,ip6,i,ilay
      REAL h,h_thrd,tol,omega_bz
c     
c     determine distance between planes (h) and number of k-points per plane (nkpp)
c     
      tol = 1.0e-15
      DO ikpt = 2,nkpt
         h = abs(bk(3,ikpt)-bk(3,1))
         IF (h.GT.tol) EXIT
      ENDDO
      nkpp = ikpt - 1
      h_thrd = h / 3.0
c     
c     make tetrahedrons
c     
      ntetra = 0
      DO ilay = 0, (nkpt/nkpp)-2
         DO itri = 1,ntria
            ip1 = itria(1,itri) + nkpp*ilay ; ip4 = ip1 + nkpp
            ip2 = itria(2,itri) + nkpp*ilay ; ip5 = ip2 + nkpp
            ip3 = itria(3,itri) + nkpp*ilay ; ip6 = ip3 + nkpp
c     
            ntetra = ntetra + 1
            itetra(1,ntetra) = ip1 ; itetra(2,ntetra) = ip2
            itetra(3,ntetra) = ip3 ; itetra(4,ntetra) = ip4
            voltet(ntetra) = h_thrd * atr(itri)
c     
            ntetra = ntetra + 1
            itetra(1,ntetra) = ip4 ; itetra(2,ntetra) = ip5
            itetra(3,ntetra) = ip6 ; itetra(4,ntetra) = ip2
            voltet(ntetra) = h_thrd * atr(itri)
c     
            ntetra = ntetra + 1
            itetra(1,ntetra) = ip2 ; itetra(2,ntetra) = ip3
            itetra(3,ntetra) = ip4 ; itetra(4,ntetra) = ip6
            voltet(ntetra) = h_thrd * atr(itri)
c     
         ENDDO
      ENDDO

      omega_bz = 0.0
      DO itet = 1,ntetra
         omega_bz = omega_bz + voltet(itet)
      ENDDO
      DO itet = 1,ntetra
         voltet(itet) =  voltet(itet) /omega_bz
      ENDDO

      RETURN
      END SUBROUTINE make_tetra
      END MODULE m_maketetra
