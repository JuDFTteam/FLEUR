MODULE m_BfieldtoVmat
   USE m_types
   USE m_constants

   IMPLICIT NONE

CONTAINS

   SUBROUTINE BfieldtoVmat(sym, stars, atoms, sphhar, vacuum, &
                          vScal, bx, by, bz, vMat)
      USE m_fft3d

      TYPE(t_sym),      INTENT(IN)  :: sym
      TYPE(t_stars),    INTENT(IN)  :: stars
      TYPE(t_atoms),    INTENT(IN)  :: atoms
      TYPE(t_sphhar),   INTENT(IN)  :: sphhar
      TYPE(t_vacuum),   INTENT(IN)  :: vacuum
      TYPE(t_potden),   INTENT(IN)  :: vScal, bx, by, bz
      TYPE(t_potden),   INTENT(OUT) :: vMat

      ! Local scalars: iteration indices, matrix elements etc.
      INTEGER iden, ifft3, ityp, iri, ilh, imesh
      REAL zero, rho_11, rho_22, rerho_21, imrho_21, cdn11, cdn22, recdn21, imcdn21
      COMPLEX czero

      ! Local arrays: densities in real space and off-diagonal elements.
      REAL,    ALLOCATABLE        :: ris(:,:), ris2(:,:), fftwork(:)
      REAL,    ALLOCATABLE        :: rho(:,:,:,:)
      COMPLEX, ALLOCATABLE        :: qpw(:,:), qpww(:,:)

      zero  = 0.0; czero = CMPLX(0.0,0.0)
      ifft3 = 27*stars%mx1*stars%mx2*stars%mx3

      ! Allocation of arrays and initialization of those that make up the real
      ! space density matrix.
      ALLOCATE (rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,4), qpw(stars%ng3,4), &
                qpww(stars%ng3,4), fftwork(0:27*stars%mx1*stars%mx2*stars%mx3-1), &
                ris(0:27*stars%mx1*stars%mx2*stars%mx3-1,4), ris2(0:27*stars%mx1*stars%mx2*stars%mx3-1,4))

      rho(:,:,:,:) = zero; qpw(:,:) = czero

      rho(:,0:,:,1) = vScal%mt(:,0:,:,1)
      rho(:,0:,:,2) = bx%mt(:,0:,:,1)
      rho(:,0:,:,3) = by%mt(:,0:,:,1)
      rho(:,0:,:,4) = bz%mt(:,0:,:,1)
      qpw(1:,1)     = vScal%pw(1:,1)
      qpw(1:,2)     = bx%pw(1:,1)
      qpw(1:,3)     = by%pw(1:,1)
      qpw(1:,4)     = bz%pw(1:,1)
      qpww(1:,1)    = vScal%pw_w(1:,1)
      qpww(1:,2)    = bx%pw_w(1:,1)
      qpww(1:,3)    = by%pw_w(1:,1)
      qpww(1:,4)    = bz%pw_w(1:,1)

      ! Calculate the charge and magnetization densities in the muffin tins.

      DO ityp = 1,atoms%ntype
         DO ilh = 0,sphhar%nlh(sym%ntypsy(atoms%firstAtom(ityp)))
            DO iri = 1,atoms%jri(ityp)
               cdn11 = rho(iri,ilh,ityp,1)+rho(iri,ilh,ityp,4)
               cdn22 = rho(iri,ilh,ityp,1)-rho(iri,ilh,ityp,4)
               recdn21 = rho(iri,ilh,ityp,2)
               imcdn21 = rho(iri,ilh,ityp,3)

               rho(iri,ilh,ityp,1) = cdn11
               rho(iri,ilh,ityp,2) = cdn22
               rho(iri,ilh,ityp,3) = recdn21
               rho(iri,ilh,ityp,4) = imcdn21
            END DO
         END DO
      END DO

      ! Fourier transform the diagonal part of the density matrix in the
      ! interstitial (qpw) to real space (ris).
      DO iden = 1,4
         CALL fft3d(ris(0:,iden),fftwork,qpw(1,iden),stars,1)
         CALL fft3d(ris2(0:,iden),fftwork,qpww(1,iden),stars,1)
      END DO

      DO imesh = 0,ifft3-1
         rho_11  = ris(imesh,1)+ris(imesh,4)
         rho_22  = ris(imesh,1)-ris(imesh,4)
         rerho_21  = ris(imesh,2)
         imrho_21  = ris(imesh,3)

         ris(imesh,1) = rho_11
         ris(imesh,2) = rho_22
         ris(imesh,3) = rerho_21
         ris(imesh,4) = imrho_21

         rho_11  = ris2(imesh,1)+ris2(imesh,4)
         rho_22  = ris2(imesh,1)-ris2(imesh,4)
         rerho_21  = ris2(imesh,2)
         imrho_21  = ris2(imesh,3)

         ris2(imesh,1) = rho_11
         ris2(imesh,2) = rho_22
         ris2(imesh,3) = rerho_21
         ris2(imesh,4) = imrho_21
      END DO

      DO iden = 1,2
         fftwork=zero
         CALL fft3d(ris(0:,iden),fftwork,qpw(1,iden),stars,-1)
         fftwork=zero
         CALL fft3d(ris2(0:,iden),fftwork,qpww(1,iden),stars,-1)
      END DO

      CALL fft3d(ris(0:,3),ris(0:,4),qpw(1,3),stars,-1)
      CALL fft3d(ris2(0:,3),ris2(0:,4),qpww(1,3),stars,-1)

      CALL vMat%init_potden_simple(stars%ng3,atoms%jmtd,atoms%msh,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,atoms%n_vPairs,2,.TRUE.,.TRUE.,&
                                      POTDEN_TYPE_POTTOT,vacuum%nmzd,&
                                      vacuum%nmzxyd,stars%ng2)
                                     
      ALLOCATE (vMat%pw_w, mold=vMat%pw)
      vMat%mt(:,0:,1:,1:4) = rho(:,0:,1:,1:4)
      vMat%pw(1:,1:3) = qpw(1:,1:3)
      vMat%pw_w(1:,1:3) = qpww(1:,1:3)

      DEALLOCATE (rho, qpw, qpww, fftwork, ris, ris2)

   END SUBROUTINE BfieldtoVmat

END MODULE m_BfieldtoVmat
