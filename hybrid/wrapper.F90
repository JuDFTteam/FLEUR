# define test
module m_wrapper

   interface blockmat
      module procedure blockmat_d, blockmat_z
   end interface

   interface packmat
      module procedure packmat_d, packmat_z
   end interface

   interface packmatcoul
      module procedure packmatcoul_d, packmatcoul_z
   end interface

   interface unpackmat
      module procedure unpackmat_d, unpackmat_z
   end interface

   interface dotprod
      module procedure dotprod_dd, dotprod_dz, dotprod_zd, dotprod_zz
   end interface

   interface matvec
      module procedure matvec_dpd, matvec_dpz, matvec_zpd, matvec_zpz
   end interface

   interface matmat
      module procedure &
         matmat_dpdp, matmat_dpzp, matmat_zpdp, matmat_zpzp, &
         matmat_dmzp, matmat_dpzm, matmat_dmzm, &
         matmat_zmdp, matmat_zpdm, matmat_zmdm, &
         matmat_zmzp, matmat_zpzm, matmat_zmzm
   end interface

   interface matmatmatd
      module procedure matmatmatd_ddd, matmatmatd_dzd
   end interface

   interface diagonalize
!       module procedure  diagonalize_de,  diagonalize_dv,  diagonalize_dpe,  diagonalize_dpv,
!      &                  diagonalize_ze,  diagonalize_zv,  diagonalize_zpe,  diagonalize_zpv,
!      &                  diagonalize_deo, diagonalize_dvo, diagonalize_dpeo, diagonalize_dpvo,
!      &                  diagonalize_zeo, diagonalize_zvo, diagonalize_zpeo, diagonalize_zpvo
      module procedure &
         diagonalize_de, diagonalize_dv, diagonalize_dpe, &
         diagonalize_dpv, diagonalize_ze, diagonalize_zv, &
         diagonalize_zpe, diagonalize_zpv, diagonalize_deo, &
         diagonalize_dvo, diagonalize_dpeo, diagonalize_dpvo, &
         diagonalize_zeo, diagonalize_zvo, diagonalize_zpeo, &
         diagonalize_zpvo, diagonalize_dvs, diagonalize_dvos, &
         diagonalize_dpvs, diagonalize_dpvos, diagonalize_zvs, &
         diagonalize_zvos, diagonalize_zpvs, diagonalize_zpvos, &
         diagonalize_dvx, diagonalize_dvox, diagonalize_dpvx, &
         diagonalize_dpvox, diagonalize_zvx, diagonalize_zvox, &
         diagonalize_zpvx, diagonalize_zpvox
   end interface

   interface geteigen
      module procedure geteigen_zpvo
   end interface

   interface inverse
      module procedure inverse_d, inverse_dp, inverse_z, inverse_zp, &
         inverse_d1, inverse_dp1, inverse_z1, inverse_zp1
   end interface

   interface sqrtmat
      module procedure sqrtmat_d, sqrtmat_dp, sqrtmat_z, sqrtmat_zp, &
         sqrtmat_d1, sqrtmat_dp1, sqrtmat_z1, sqrtmat_zp1
   end interface

contains

!     --------

   function identity(n)
      implicit none
      integer, intent(in)  :: n
      integer ::              :: identity(n, n)
      integer ::              :: i, j
      identity = 0
      do i = 1, n
         identity(i, i) = 1
      enddo
   end function identity

!     --------

   function blockmat_d(a, b)
      implicit none
      real, intent(in) :: a(:, :), b(:, :)
      real ::             :: blockmat_d(size(a, 1) + size(b, 1), &
                                        size(a, 1) + size(b, 1))
      integer ::             :: na, nb
      na = size(a, 1); nb = size(b, 1)
      if (size(a, 2) /= na) &
         stop 'blockmat_d: dimensions of first array differ.'
      if (size(b, 2) /= nb) &
         stop 'blockmat_d: dimensions of second array differ.'
      blockmat_d = 0.0
      blockmat_d(:na, :na) = a
      blockmat_d(na + 1:, na + 1:) = b
   end function blockmat_d

   function blockmat_z(a, b)
      implicit none
      complex, intent(in) :: a(:, :), b(:, :)
      complex             :: blockmat_z(size(a, 1) + size(b, 1), &
                                        size(a, 1) + size(b, 1))
      integer ::                :: na, nb
      na = size(a, 1); nb = size(b, 1)
      if (size(a, 2) /= na) &
         stop 'blockmat_z: dimensions of first array differ.'
      if (size(b, 2) /= nb) &
         stop 'blockmat_z: dimensions of second array differ.'
      blockmat_z = 0.0
      blockmat_z(:na, :na) = a
      blockmat_z(na + 1:, na + 1:) = b
   end function blockmat_z

!     --------

   function packmat_d(mat)
      implicit none
      real, intent(in)  :: mat(:, :)
      real ::              :: packmat_d(size(mat, 1)*(size(mat, 1) + 1)/2)
      integer ::              :: n, nn, i, j, k
      n = size(mat, 1); nn = n*(n + 1)/2
      if (size(mat, 2) /= n) stop 'packmat_d: array dimensions differ.'
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            packmat_d(k) = mat(i, j)
# ifdef test
            if (abs(mat(j, i) - mat(i, j)) > 10.0**-8) &
               STOP 'packmat_d: input matrix not symmetric'
# endif
         enddo
      enddo
   end function packmat_d

   function packmatcoul_d(mat)
      implicit none
      real, intent(in)  :: mat(:, :)
      real ::              :: packmatcoul_d( &
                                      size(mat, 1)*(size(mat, 1) + 1)/2)
      integer ::              :: n, nn, i, j, k
      n = size(mat, 1); nn = n*(n + 1)/2
      if (size(mat, 2) /= n) stop 'packmat_d: array dimensions differ.'
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1

            packmatcoul_d(k) = (mat(i, j) + mat(j, i))/2.

            !           if(abs(mat(j,i)-mat(i,j)).gt.10.0**-6) then
            !             write(*,*) 'packmatcoul_d: input matrix not symmetric; deviation .gt. 1E-06'
            !           endif
         enddo
      enddo
   end function packmatcoul_d

   function unpackmat_d(mat)
      implicit none
      real, intent(in)  :: mat(:)
      real ::              :: unpackmat_d( &
                                      nint(sqrt(0.25 + 2*size(mat)) - 0.5), &
                                      nint(sqrt(0.25 + 2*size(mat)) - 0.5))
      integer ::              :: n, nn, i, j, k
      nn = size(mat); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            unpackmat_d(i, j) = mat(k)
            unpackmat_d(j, i) = mat(k)
         enddo
      enddo
   end function unpackmat_d

   function packmat_z(mat)
      implicit none
      complex, intent(in) :: mat(:, :)
      complex             :: packmat_z(size(mat, 1)*(size(mat, 1) + 1)/2)
      integer ::                :: n, nn, i, j, k
      n = size(mat, 1); nn = n*(n + 1)/2
      if (size(mat, 2) /= n) stop 'packmat_z: array dimensions differ.'
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            packmat_z(k) = mat(i, j)
# ifdef test
            if (abs(conjg(mat(j, i)) - mat(i, j)) > 10.0**-8) &
               stop 'packmat_z: input matrix not Hermitian.'
# endif
         enddo
      enddo
   end function packmat_z

   function packmatcoul_z(mat)
      implicit none
      complex, intent(in)  :: mat(:, :)
      complex              :: packmatcoul_z( &
                              size(mat, 1)*(size(mat, 1) + 1)/2)
      integer ::                 :: n, nn, i, j, k
      n = size(mat, 1); nn = n*(n + 1)/2
      if (size(mat, 2) /= n) stop 'packmat_z: array dimensions differ.'
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            packmatcoul_z(k) = (mat(i, j) + conjg(mat(j, i)))/2.

            if (abs(conjg(mat(j, i)) - mat(i, j)) > 10.0**-4) then
               stop 'packmatcoul_z: input matrix not Hermitian; & &
                  deviation > 1E-04.'
            endif
         enddo
      enddo
   end function packmatcoul_z

   function unpackmat_z(mat)
      implicit none
      complex, intent(in)  :: mat(:)
      complex              :: unpackmat_z( &
                              nint(sqrt(0.25 + 2*size(mat)) - 0.5), &
                              nint(sqrt(0.25 + 2*size(mat)) - 0.5))
      integer ::                 :: n, nn, i, j, k
      nn = size(mat); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            unpackmat_z(i, j) = mat(k)
            unpackmat_z(j, i) = conjg(mat(k))
         enddo
      enddo
   end function unpackmat_z

!     --------

   function dotprod_dd(vec1, vec2)
      implicit none
      real, intent(in) :: vec1(:), vec2(:)
      real ::             :: dotprod_dd
      integer ::             :: n
      real ::             :: ddot
      n = size(vec1)
      if (size(vec2) /= n) &
         stop 'dotprod_dd: sizes of input vectors differ.'
      dotprod_dd = ddot(n, vec1, 1, vec2, 1)
   end function dotprod_dd

   function dotprod_dz(vec1, vec2)
      implicit none
      real, intent(in) :: vec1(:)
      complex, intent(in) :: vec2(:)
      complex             :: dotprod_dz
      integer ::                :: n
      real ::                :: ddot
      n = size(vec1)
      if (size(vec2) /= n) &
         stop 'dotprod_dz: sizes of input vectors differ.'
      dotprod_dz = ddot(n, vec1, 1, real(vec2), 1) &
                   + (0.0, 1.0)*ddot(n, vec1, 1, aimag(vec2), 1)
   end function dotprod_dz

   function dotprod_zd(vec1, vec2)
      implicit none
      complex, intent(in) :: vec1(:)
      real, intent(in) :: vec2(:)
      complex             :: dotprod_zd
      integer ::                :: n
      real ::                :: ddot
      n = size(vec1)
      if (size(vec2) /= n) &
         stop 'dotprod_zd: sizes of input vectors differ.'
      dotprod_zd = ddot(n, real(vec1), 1, vec2, 1) &
                   - (0.0, 1.0)*ddot(n, aimag(vec1), 1, vec2, 1)
   end function dotprod_zd

   function dotprod_zz(vec1, vec2)
      implicit none
      complex, intent(in) :: vec1(:), vec2(:)
      complex             :: dotprod_zz
      integer ::                :: n
      complex             :: zdotc
      n = size(vec1)
      if (size(vec2) /= n) &
         stop 'dotprod_zz: sizes of input vectors differ.'
      dotprod_zz = zdotc(n, vec1, 1, vec2, 1)
   end function dotprod_zz

!     --------

   function matvec_dpd(mat, vec)
      implicit none
      real, intent(in)  :: mat(:), vec(:)
      real ::              :: matvec_dpd(size(vec))
      integer ::              :: nn, n
      n = size(vec)
      nn = n*(n + 1)/2
      if (size(mat) /= nn) stop 'matvec_dpd: input array has wrong size.'
      call dspmv('U', n, 1.0, mat, vec, 1, 0.0, matvec_dpd, 1)
   end function matvec_dpd

   function matvec_dpz(mat, vec)
      implicit none
      real, intent(in) :: mat(:)
      complex, intent(in) :: vec(:)
      complex             :: matvec_dpz(size(vec))
      real, allocatable :: vecr(:), veci(:)
      integer ::                :: nn, n
      n = size(vec); allocate (vecr(n), veci(n))
      nn = n*(n + 1)/2
      if (size(mat) /= nn) stop 'matvec_dpz: input array has wrong size.'
      call dspmv('U', n, 1.0, mat, real(vec), 1, 0.0, vecr, 1)
      call dspmv('U', n, 1.0, mat, aimag(vec), 1, 0.0, veci, 1)
      matvec_dpz = vecr + (0.0, 1.0)*veci
      deallocate (vecr, veci)
   end function matvec_dpz

   function matvec_zpd(mat, vec)
      implicit none
      complex, intent(in) :: mat(:)
      real, intent(in) :: vec(:)
      complex             :: matvec_zpd(size(vec))
      real, allocatable :: vecr(:), veci(:)
      integer ::                :: nn, n
      n = size(vec); allocate (vecr(n), veci(n))
      nn = n*(n + 1)/2
      if (size(mat) /= nn) stop 'matvec_zpd: input array has wrong size.'
      call dspmv('U', n, 1.0, real(mat), vec, 1, 0.0, vecr, 1)
      call dspmv('U', n, 1.0, aimag(mat), vec, 1, 0.0, veci, 1)
      matvec_zpd = vecr + (0.0, 1.0)*veci
      deallocate (vecr, veci)
   end function matvec_zpd

   function matvec_zpz(mat, vec)
      implicit none
      complex, intent(in)  :: mat(:), vec(:)
      complex              :: matvec_zpz(size(vec))
      integer ::                 :: nn, n
      n = size(vec)
      nn = n*(n + 1)/2
      if (size(mat) /= nn) stop 'matvec_zpz: input array has wrong size.'
      call zhpmv('U', n, (1.0, 0.0), mat, vec, 1, (0.0, 0.0), matvec_zpz, 1)
   end function matvec_zpz

!     --------

   function matmat_dpdp(mat1, mat2)
      implicit none
      real, intent(in)  :: mat1(:), mat2(:)
      real ::              :: matmat_dpdp( &
                                      nint(sqrt(0.25 + 2*size(mat1)) - 0.5), &
                                      nint(sqrt(0.25 + 2*size(mat1)) - 0.5))
      real, allocatable :: vec(:), vec2(:)
      integer ::              :: nn, n, k1, i, j, k
      nn = size(mat1)
      n = nint(sqrt(0.25 + 2*nn) - 0.5); allocate (vec(n), vec2(n))
      if (size(mat2) /= nn) &
         stop 'matmat_dpdp: second input array has wrong size.'
      k = 0
      do i = 1, n
         vec2(:i) = mat2(k + 1:k + i)
         k1 = k + 2*i
         do j = i + 1, n
            vec2(j) = mat2(k1)
            k1 = k1 + j
         enddo
         call dspmv('U', n, 1.0, mat1, vec2, 1, 0.0, vec, 1)
         matmat_dpdp(:, i) = vec
         k = k + i
      enddo
      deallocate (vec, vec2)
   end function matmat_dpdp

   function matmat_dpzp(mat1, mat2)
      implicit none
      real, intent(in)  :: mat1(:)
      complex, intent(in)  :: mat2(:)
      complex              :: matmat_dpzp( &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5), &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5))
      real, allocatable :: vecr(:), veci(:)
      complex, allocatable :: vec2(:)
      integer ::                 :: nn, n, k1, i, j, k
      nn = size(mat1)
      n = nint(sqrt(0.25 + 2*nn) - 0.5)
      allocate (vecr(n), veci(n), vec2(n))
      if (size(mat2) /= nn) &
         stop 'matmat_dpzp: second input array has wrong size.'
      k = 0
      do i = 1, n
         vec2(:i) = mat2(k + 1:k + i)
         k1 = k + 2*i
         do j = i + 1, n
            vec2(j) = conjg(mat2(k1))
            k1 = k1 + j
         enddo
         call dspmv('U', n, 1.0, mat1, real(vec2), 1, 0.0, vecr, 1)
         call dspmv('U', n, 1.0, mat1, aimag(vec2), 1, 0.0, veci, 1)
         matmat_dpzp(:, i) = vecr + (0.0, 1.0)*veci
         k = k + i
      enddo
      deallocate (vecr, veci, vec2)
   end function matmat_dpzp

   function matmat_zpdp(mat1, mat2)
      implicit none
      complex, intent(in)  :: mat1(:)
      real, intent(in)  :: mat2(:)
      complex              :: matmat_zpdp( &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5), &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5))
      real, allocatable :: vecr(:), veci(:)
      complex, allocatable :: vec1(:)
      integer ::                 :: nn, n, k1, i, j, k
      nn = size(mat1)
      n = nint(sqrt(0.25 + 2*nn) - 0.5)
      allocate (vecr(n), veci(n), vec1(n))
      if (size(mat2) /= nn) &
         stop 'matmat_zpdp: second input array has wrong size.'
      k = 0
      do i = 1, n
         vec1(:i) = conjg(mat1(k + 1:k + i))
         k1 = k + 2*i
         do j = i + 1, n
            vec1(j) = mat1(k1)
            k1 = k1 + j
         enddo
         call dspmv('U', n, 1.0, mat2, real(vec1), 1, 0.0, vecr, 1)
         call dspmv('U', n, 1.0, mat2, aimag(vec1), 1, 0.0, veci, 1)
         matmat_zpdp(i, :) = vecr + (0.0, 1.0)*veci
         k = k + i
      enddo
      deallocate (vecr, veci, vec1)
   end function matmat_zpdp

   function matmat_zpzp(mat1, mat2)
      implicit none
      complex, intent(in)  :: mat1(:), mat2(:)
      complex              :: matmat_zpzp( &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5), &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5))
      complex, allocatable :: vec(:), vec2(:)
      integer ::                 :: nn, n, k1, i, j, k
      nn = size(mat1)
      n = nint(sqrt(0.25 + 2*nn) - 0.5); allocate (vec(n), vec2(n))
      if (size(mat2) /= nn) &
         stop 'matmat_zpzp: second input array has wrong size.'
      k = 0
      do i = 1, n
         vec2(:i) = mat2(k + 1:k + i)
         k1 = k + 2*i
         do j = i + 1, n
            vec2(j) = conjg(mat2(k1))
            k1 = k1 + j
         enddo
         call zhpmv('U', n, (1.0, 0.0), mat1, vec2, 1, (0.0, 0.0), vec, 1)
         matmat_zpzp(:, i) = vec
         k = k + i
      enddo
      deallocate (vec, vec2)
   end function matmat_zpzp

   function matmat_dpdm(mat1, mat2)
      implicit none
      real, intent(in)  :: mat1(:), mat2(:, :)
      real ::              :: matmat_dpdm(size(mat2, 1), size(mat2, 1))
      real, allocatable :: vec(:), vec2(:)
      integer ::              :: nn, n, k1, i
      n = size(mat2, 1); nn = n*(n + 1)/2; allocate (vec(n), vec2(n))
      if (size(mat2, 2) /= n) &
         stop 'matmat_dpdm: dimensions of second input array differ.'
      if (size(mat1) /= nn) &
         stop 'matmat_dpdm: first input array has wrong size.'
      do i = 1, n
         vec2 = mat2(:, i)
         call dspmv('U', n, 1.0, mat1, vec2, 1, 0.0, vec, 1)
         matmat_dpdm(:, i) = vec
      enddo
      deallocate (vec, vec2)
   end function matmat_dpdm

   function matmat_dmdp(mat1, mat2)
      implicit none
      real, intent(in)  :: mat1(:, :), mat2(:)
      real ::              :: matmat_dmdp(size(mat1, 1), size(mat1, 1))
      real, allocatable :: vec(:), vec2(:)
      integer ::              :: nn, n, k1, i
      n = size(mat1, 1); nn = n*(n + 1)/2; allocate (vec(n), vec2(n))
      if (size(mat1, 2) /= n) &
         stop 'matmat_dmdp: dimensions of first input array differ.'
      if (size(mat2) /= nn) &
         stop 'matmat_dmdp: second input array has wrong size.'
      do i = 1, n
         vec2 = mat1(i, :)
         call dspmv('U', n, 1.0, mat2, vec2, 1, 0.0, vec, 1)
         matmat_dmdp(i, :) = vec
      enddo
      deallocate (vec, vec2)
   end function matmat_dmdp

   function matmat_dmdm(mat1, mat2)
      implicit none
      real, intent(in) :: mat1(:, :), mat2(:, :)
      real ::             :: matmat_dmdm(size(mat1, 1), size(mat1, 1))
      integer ::             :: n
      n = size(mat1, 1)
      if (size(mat1, 2) /= n) &
         stop 'matmat_dmdm: dimensions of first input array differ.'
      if (size(mat2, 1) /= n) &
         stop 'matmat_dmdm: second input array has wrong dimensions.'
      if (size(mat2, 2) /= n) &
         stop 'matmat_dmdm: dimensions of second input array differ.'
      call dgemm('N', 'N', n, n, n, 1.0, mat1, n, mat2, n, 0.0, matmat_dmdm, n)
   end function matmat_dmdm

   function matmat_dpzm(mat1, mat2)
      implicit none
      real, intent(in)  :: mat1(:)
      complex, intent(in)  :: mat2(:, :)
      complex              :: matmat_dpzm(size(mat2, 1), size(mat2, 1))
      real, allocatable :: vecr(:), veci(:)
      complex, allocatable :: vec2(:)
      integer ::                 :: nn, n, k1, i
      n = size(mat2, 1)
      nn = n*(n + 1)/2; allocate (vecr(n), veci(n), vec2(n))
      if (size(mat2, 2) /= n) &
         stop 'matmat_dpzm: dimensions of second input array differ.'
      if (size(mat1) /= nn) &
         stop 'matmat_dpzm: first input array has wrong size.'
      do i = 1, n
         vec2 = mat2(:, i)
         call dspmv('U', n, 1.0, mat1, real(vec2), 1, 0.0, vecr, 1)
         call dspmv('U', n, 1.0, mat1, aimag(vec2), 1, 0.0, veci, 1)
         matmat_dpzm(:, i) = vecr + (0.0, 1.0)*veci
      enddo
      deallocate (vecr, veci, vec2)
   end function matmat_dpzm

   function matmat_dmzp(mat1, mat2)
      implicit none
      real, intent(in)  :: mat1(:, :)
      complex, intent(in)  :: mat2(:)
      complex              :: matmat_dmzp(size(mat1, 1), size(mat1, 1))
      complex, allocatable :: vec1(:), vec(:)
      integer ::                 :: nn, n, k1, i
      n = size(mat1, 1); nn = n*(n + 1)/2; allocate (vec(n), vec1(n))
      if (size(mat1, 2) /= n) &
         stop 'matmat_dmzp: dimensions of first input array differ.'
      if (size(mat2) /= nn) &
         stop 'matmat_dmzp: second input array has wrong size.'
      do i = 1, n
         vec1 = mat1(i, :)
         call zhpmv('U', n, (1.0, 0.0), mat2, vec1, 1, (0.0, 0.0), vec, 1)
         matmat_dmzp(i, :) = conjg(vec)
      enddo
      deallocate (vec, vec1)
   end function matmat_dmzp

   function matmat_dmzm(mat1, mat2)
      implicit none
      real, intent(in) :: mat1(:, :)
      complex, intent(in) :: mat2(:, :)
      complex             :: matmat_dmzm(size(mat1, 1), size(mat2, 2))
      real ::                :: matr(size(mat1, 1), size(mat2, 2)), &
                                        mati(size(mat1, 1), size(mat2, 2))
      integer ::                :: n, n1, n2
      n1 = size(mat1, 1)
      n = size(mat1, 2)
      n2 = size(mat2, 2)
      if (size(mat2, 1) /= n) &
         stop 'matmat_dmzm: dimensions of matrices are inconsistent.'
      call dgemm('N', 'N', n1, n2, n, 1.0, mat1, n1, real(mat2), n, 0.0, matr, n1)
      call dgemm('N', 'N', n1, n2, n, 1.0, mat1, n1, aimag(mat2), n, 0.0, mati, n1)
      matmat_dmzm = matr + (0.0, 1.0)*mati
   end function matmat_dmzm

   function matmat_zpdm(mat1, mat2)
      implicit none
      complex, intent(in)  :: mat1(:)
      real, intent(in)  :: mat2(:, :)
      complex              :: matmat_zpdm(size(mat2, 1), size(mat2, 1))
      complex, allocatable :: vec(:), vec2(:)
      integer ::                 :: nn, n, k1, i
      n = size(mat2, 1); nn = n*(n + 1)/2; allocate (vec(n), vec2(n))
      if (size(mat2, 2) /= n) &
         stop 'matmat_zpdm: dimensions of second input array differ.'
      if (size(mat1) /= nn) &
         stop 'matmat_zpdm: first input array has wrong size.'
      do i = 1, n
         vec2 = mat2(:, i)
         call zhpmv('U', n, (1.0, 0.0), mat1, vec2, 1, (0.0, 0.0), vec, 1)
         matmat_zpdm(:, i) = vec
      enddo
      deallocate (vec, vec2)
   end function matmat_zpdm

   function matmat_zmdp(mat1, mat2)
      implicit none
      complex, intent(in)  :: mat1(:, :)
      real, intent(in)  :: mat2(:)
      complex              :: matmat_zmdp(size(mat1, 1), size(mat1, 1))
      complex, allocatable :: vec1(:)
      real, allocatable :: vecr(:), veci(:)
      integer ::                 :: nn, n, k1, i
      n = size(mat1, 1); nn = n*(n + 1)/2
      allocate (vecr(n), veci(n), vec1(n))
      if (size(mat1, 2) /= n) &
         stop 'matmat_zmdp: dimensions of first input array differ.'
      if (size(mat2) /= nn) &
         stop 'matmat_zmdp: second input array has wrong size.'
      do i = 1, n
         vec1 = conjg(mat1(i, :))
         call dspmv('U', n, 1.0, mat2, real(vec1), 1, 0.0, vecr, 1)
         call dspmv('U', n, 1.0, mat2, aimag(vec1), 1, 0.0, veci, 1)
         matmat_zmdp(i, :) = vecr - (0.0, 1.0)*veci
      enddo
      deallocate (vecr, veci, vec1)
   end function matmat_zmdp

   function matmat_zmdm(mat1, mat2)
      implicit none
      complex, intent(in) :: mat1(:, :)
      real, intent(in) :: mat2(:, :)
      complex             :: matmat_zmdm(size(mat1, 1), size(mat2, 2))
      real ::                :: matr(size(mat1, 1), size(mat2, 2)), &
                                        mati(size(mat1, 1), size(mat2, 2))
      integer ::                :: n, n1, n2
      n1 = size(mat1, 1)
      n = size(mat1, 2)
      n2 = size(mat2, 2)
      if (size(mat2, 1) /= n) &
         stop 'matmat_zmdm: dimensions of matrices are inconsistent.'
      call dgemm('N', 'N', n1, n2, n, 1.0, real(mat1), n1, mat2, n, 0.0, matr, n1)
      call dgemm('N', 'N', n1, n2, n, 1.0, aimag(mat1), n1, mat2, n, 0.0, mati, n1)
      matmat_zmdm = matr + (0.0, 1.0)*mati
   end function matmat_zmdm

   function matmat_zpzm(mat1, mat2)
      implicit none
      complex, intent(in)  :: mat1(:), mat2(:, :)
      complex              :: matmat_zpzm(size(mat2, 1), size(mat2, 2))
      complex, allocatable :: vec(:), vec2(:)
      integer ::                 :: nn, n, k1, i, n2
      n = size(mat2, 1); nn = n*(n + 1)/2; allocate (vec(n), vec2(n))
      n2 = size(mat2, 2)
      if (size(mat1) /= nn) &
         stop 'matmat_zpzm: first input array has wrong size.'
      do i = 1, n2
         vec2 = mat2(:, i)
         call zhpmv('U', n, (1.0, 0.0), mat1, vec2, 1, (0.0, 0.0), vec, 1)
         matmat_zpzm(:, i) = vec
      enddo
      deallocate (vec, vec2)
   end function matmat_zpzm

   function matmat_zmzp(mat1, mat2)
      implicit none
      complex, intent(in)  :: mat1(:, :), mat2(:)
      complex              :: matmat_zmzp(size(mat1, 1), size(mat1, 1))
      complex, allocatable :: vec(:), vec2(:)
      integer ::                 :: nn, n, k1, i
      n = size(mat1, 1); nn = n*(n + 1)/2; allocate (vec(n), vec2(n))
      if (size(mat1, 2) /= n) &
         stop 'matmat_zmzp: dimensions of first input array differ.'
      if (size(mat2) /= nn) &
         stop 'matmat_zmzp: second input array has wrong size.'
      do i = 1, n
         vec2 = conjg(mat1(i, :))
         call zhpmv('U', n, (1.0, 0.0), mat2, vec2, 1, (0.0, 0.0), vec, 1)
         matmat_zmzp(i, :) = conjg(vec)
      enddo
      deallocate (vec, vec2)
   end function matmat_zmzp

   function matmat_zmzm(mat1, mat2)
      implicit none
      complex, intent(in) :: mat1(:, :), mat2(:, :)
      complex             :: matmat_zmzm(size(mat1, 1), size(mat2, 2))
      integer ::                :: n1, n, n2
      complex, parameter     :: one = (1, 0), zero = 0
      n1 = size(mat1, 1)
      n = size(mat1, 2)
      n2 = size(mat2, 2)
      if (size(mat2, 1) /= n) &
         stop 'matmat_zmzm: dimensions of matrices are inconsistent.'
      call zgemm('N', 'N', n1, n2, n, one, mat1, n1, mat2, n, zero, matmat_zmzm, n1)
   end function matmat_zmzm

!     --------

   function matmatmatd_ddd(diag1, mat, diag2)
      implicit none
      real, intent(in) :: diag1(:), mat(:, :), diag2(:)
      real ::             :: matmatmatd_ddd(size(diag1), size(diag1))
      integer ::             :: n, i
      n = size(diag1)
      if (size(mat, 1) /= n) &
         stop 'matmatmatd_ddd: input matrix has wrong size.'
      if (size(mat, 2) /= n) &
         stop 'matmatmatd_ddd: dimensions of input matrix differ.'
      if (size(diag2) /= n) &
         stop 'matmatmatd_ddd: second diagonal matrix has wrong size.'
      do i = 1, n
         matmatmatd_ddd(:, i) = diag2(i)*mat(:, i)
      enddo
      do i = 1, n
         matmatmatd_ddd(i, :) = diag1(i)*matmatmatd_ddd(i, :)
      enddo
   end function matmatmatd_ddd

   function matmatmatd_dzd(diag1, mat, diag2)
      implicit none
      real, intent(in) :: diag1(:), diag2(:)
      complex, intent(in) :: mat(:, :)
      complex             :: matmatmatd_dzd(size(diag1), size(diag1))
      integer ::                :: n, i
      n = size(diag1)
      if (size(mat, 1) /= n) &
         stop 'matmatmatd_ddd: input matrix has wrong size.'
      if (size(mat, 2) /= n) &
         stop 'matmatmatd_ddd: dimensions of input matrix differ.'
      if (size(diag2) /= n) &
         stop 'matmatmatd_ddd: second diagonal matrix has wrong size.'
      do i = 1, n
         matmatmatd_dzd(:, i) = diag2(i)*mat(:, i)
      enddo
      do i = 1, n
         matmatmatd_dzd(i, :) = diag1(i)*matmatmatd_dzd(i, :)
      enddo
   end function matmatmatd_dzd

!     --------

   subroutine diagonalize_de(eval, mat)
      implicit none
      real, intent(out) :: eval(:)
      real, intent(in)  :: mat(:, :)
      real, allocatable :: mat1(:, :), work(:)
      integer ::              :: n, info
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_de: zero dimension in eigenvalue problem.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_de: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_de: eigenvalue array has wrong size.'
      allocate (mat1(n, n), work(3*n)); mat1 = mat
      call dsyev('N', 'U', n, mat1, n, eval, work, 3*n, info)
      if (info /= 0) stop 'diagonalize_de: dsyev failed.'
      deallocate (mat1, work)
   end subroutine diagonalize_de

   subroutine diagonalize_dv(evec, eval, mat)
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:, :)
      real, allocatable :: work(:)
      integer ::              :: n, info
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_dv: zero dimension in eigenvalue problem.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_dv: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_dv: eigenvalue array has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_dv: eigenvector array has wrong dimensions.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_dv: dimensions of eigenvector array differ.'
      allocate (work(3*n)); evec = mat
      call dsyev('V', 'U', n, evec, n, eval, work, 3*n, info)
      if (info /= 0) stop 'diagonalize_dv: dsyev failed.'
      deallocate (work)
   end subroutine diagonalize_dv

   subroutine diagonalize_dpe(eval, mat)
      implicit none
      real, intent(out) :: eval(:)
      real, intent(in)  :: mat(:)
      real, allocatable :: mat1(:), work(:)
      integer ::              :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_dpe: zero dimension in eigenvalue problem.'
      if (size(mat) /= nn) &
         stop 'diagonalize_dpe: input matrix has wrong size.'
      allocate (mat1(nn), work(3*n)); mat1 = mat
      call dspev('N', 'U', n, mat1, eval, work, n, work, info)
      if (info /= 0) stop 'diagonalize_dpe: dspev failed.'
      deallocate (mat1, work)
   end subroutine diagonalize_dpe

   subroutine diagonalize_dpv(evec, eval, mat)
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:)
      real, allocatable :: mat1(:), work(:)
      integer ::              :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_dpv: zero dimension in eigenvalue problem.'
      if (size(mat) /= nn) &
         stop 'diagonalize_dpv: input matrix has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_dpv: eigenvector array has wrong dimensions.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_dpv: dimensions of eigenvector array differ.'
      allocate (mat1(nn), work(3*n)); mat1 = mat
      call dspev('V', 'U', n, mat1, eval, evec, n, work, info)
      if (info /= 0) stop 'diagonalize_dpv: dspev failed.'
      deallocate (mat1, work)
   end subroutine diagonalize_dpv

   subroutine diagonalize_ze(eval, mat)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(in)  :: mat(:, :)
      complex, allocatable :: mat1(:, :), work(:)
      real, allocatable :: rwork(:)
      integer ::                 :: n, info
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_ze: zero dimension in eigenvalue problem.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_ze: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_ze: eigenvalue array has wrong size.'
      allocate (mat1(n, n), work(3*n), rwork(3*n)); mat1 = mat
      call zheev('N', 'U', n, mat1, n, eval, work, 3*n, rwork, info)
      if (info /= 0) stop 'diagonalize_ze: zheev failed.'
      deallocate (mat1, work, rwork)
   end subroutine diagonalize_ze

   subroutine diagonalize_zv(evec, eval, mat)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:, :)
      complex, allocatable :: work(:)
      real, allocatable :: rwork(:)
      integer ::                 :: n, info
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_zv: zero dimension in eigenvalue problem.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_zv: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_zv: eigenvalue array has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_zv: eigenvector array has wrong dimensions.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_zv: dimensions of eigenvector array differ.'
      allocate (work(3*n), rwork(3*n)); evec = mat
      call zheev('V', 'U', n, evec, n, eval, work, 3*n, rwork, info)
      if (info /= 0) stop 'diagonalize_zv: zheev failed.'
      deallocate (work, rwork)
   end subroutine diagonalize_zv

   subroutine diagonalize_zpe(eval, mat)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(in)  :: mat(:)
      complex, allocatable :: mat1(:), work(:)
      real, allocatable :: rwork(:)
      integer ::                 :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_zpe: zero dimension in eigenvalue problem.'
      if (size(mat) /= nn) &
         stop 'diagonalize_zpe: input matrix has wrong size.'
      allocate (mat1(nn), work(3*n), rwork(3*n)); mat1 = mat
      call zhpev('N', 'U', n, mat1, eval, work, n, work, rwork, info)
      if (info /= 0) stop 'diagonalize_zpe: zhpev failed.'
      deallocate (mat1, work, rwork)
   end subroutine diagonalize_zpe

   subroutine diagonalize_zpv(evec, eval, mat)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:)
      complex, allocatable :: mat1(:), work(:)
      real, allocatable :: rwork(:)
      integer ::                 :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_zpv: zero dimension in eigenvalue problem.'
      if (size(mat) /= nn) &
         stop 'diagonalize_zpv: input matrix has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_zpv: eigenvector array has wrong dimensions.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_zpv: dimensions of eigenvector array differ.'
      allocate (mat1(nn), work(3*n), rwork(3*n)); mat1 = mat
      call zhpev('V', 'U', n, mat1, eval, evec, n, work, rwork, info)
      if (info /= 0) stop 'diagonalize_zpv: zhpev failed.'
      deallocate (mat1, work, rwork)
   end subroutine diagonalize_zpv

   subroutine diagonalize_deo(eval, mat, olap)
      implicit none
      real, intent(out) :: eval(:)
      real, intent(in)  :: mat(:, :), olap(:, :)
      real, allocatable :: mat1(:, :), olap1(:, :), work(:)
      integer ::              :: n, info
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_deo: zero dimension in eigenvalue problem.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_deo: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_deo: eigenvalue array has wrong size.'
      if (size(olap, 1) /= n) &
         stop 'diagonalize_deo: overlap matrix has wrong size.'
      if (size(olap, 2) /= n) &
         stop 'diagonalize_deo: dimensions of overlap matrix differ.'
      allocate (mat1(n, n), olap1(n, n), work(3*n))
      mat1 = mat; olap1 = olap
      call dsygv(1, 'N', 'U', n, mat1, n, olap1, n, eval, work, 3*n, info)
      if (info /= 0) stop 'diagonalize_deo: dsygv failed.'
      deallocate (mat1, olap1, work)
   end subroutine diagonalize_deo

   subroutine diagonalize_dvo(evec, eval, mat, olap)
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:, :), olap(:, :)
      real, allocatable :: olap1(:, :), work(:)
      integer ::              :: n, info
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_dvo: zero dimension in eigenvalue problem.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_dvo: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_dvo: eigenvalue array has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_dvo: eigenvector array has wrong dimensions.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_dvo: dimensions of eigenvector array differ.'
      if (size(olap, 1) /= n) &
         stop 'diagonalize_dvo: overlap matrix has wrong dimensions.'
      if (size(olap, 2) /= n) &
         stop 'diagonalize_dvo: dimensions of overlap matrix differ.'
      allocate (olap1(n, n), work(3*n)); evec = mat; olap1 = olap
      call dsygv(1, 'V', 'U', n, evec, n, olap1, n, eval, work, 3*n, info)
      if (info /= 0) stop 'diagonalize_dvo: dsygv failed.'
      deallocate (olap1, work)
   end subroutine diagonalize_dvo

   subroutine diagonalize_dpeo(eval, mat, olap)
      implicit none
      real, intent(out) :: eval(:)
      real, intent(in)  :: mat(:), olap(:)
      real, allocatable :: mat1(:), olap1(:), work(:)
      integer ::              :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_dpeo: zero dimension in eigenvalue problem.'
      if (size(mat) /= nn) &
         stop 'diagonalize_dpeo: input matrix has wrong size.'
      if (size(olap) /= nn) &
         stop 'diagonalize_dpeo: overlap matrix has wrong size.'
      allocate (mat1(nn), olap1(nn), work(3*n))
      mat1 = mat; olap1 = olap
      call dspgv(1, 'N', 'U', n, mat1, olap1, eval, work, n, work, info)
      if (info /= 0) stop 'diagonalize_dpeo: dspgv failed.'
      deallocate (mat1, olap1, work)
   end subroutine diagonalize_dpeo

   subroutine diagonalize_dpvo(evec, eval, mat, olap)
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:), olap(:)
      real, allocatable :: mat1(:), olap1(:), work(:)
      integer ::              :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_dpvo: zero dimension in eigenvalue problem.'
      if (size(mat) /= nn) &
         stop 'diagonalize_dpvo: input matrix has wrong size.'
      if (size(olap) /= nn) &
         stop 'diagonalize_dpvo: overlap matrix has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_dpvo: eigenvector array has wrong dimensions.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_dpvo: dimensions of eigenvector array differ.'
      allocate (mat1(nn), olap1(nn), work(3*n))
      mat1 = mat; olap1 = olap
      call dspgv(1, 'V', 'U', n, mat1, olap1, eval, evec, n, work, info)
      if (info /= 0) stop 'diagonalize_dpvo: dspgv failed.'
      deallocate (mat1, olap1, work)
   end subroutine diagonalize_dpvo

   subroutine diagonalize_zeo(eval, mat, olap)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(in)  :: mat(:, :), olap(:, :)
      complex, allocatable :: mat1(:, :), olap1(:, :), work(:)
      real, allocatable :: rwork(:)
      integer ::                 :: n, info
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_zeo: zero dimension in eigenvalue problem.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_zeo: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_zeo: eigenvalue array has wrong size.'
      if (size(olap, 1) /= n) &
         stop 'diagonalize_zeo: overlap matrix has wrong size.'
      if (size(olap, 2) /= n) &
         stop 'diagonalize_zeo: dimensions of overlap matrix differ.'
      allocate (mat1(n, n), olap1(n, n), work(3*n), rwork(3*n))
      mat1 = mat; olap1 = olap
      call zhegv(1, 'N', 'U', n, mat1, n, olap1, n, eval, work, 3*n, rwork, info)
      if (info /= 0) stop 'diagonalize_zeo: zhegv failed.'
      deallocate (mat1, olap1, work, rwork)
   end subroutine diagonalize_zeo

   subroutine diagonalize_zvo(evec, eval, mat, olap)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:, :), olap(:, :)
      complex, allocatable :: olap1(:, :), work(:)
      real, allocatable :: rwork(:)
      integer ::                 :: n, info
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_zvo: zero dimension in eigenvalue problem.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_zvo: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_zvo: eigenvalue array has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_zvo: eigenvector array has wrong dimensions.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_zvo: dimensions of eigenvector array differ.'
      if (size(olap, 1) /= n) &
         stop 'diagonalize_zvo: overlap matrix has wrong dimensions.'
      if (size(olap, 2) /= n) &
         stop 'diagonalize_zvo: dimensions of overlap matrix differ.'
      allocate (olap1(n, n), work(3*n), rwork(3*n))
      evec = mat; olap1 = olap
      call zhegv(1, 'V', 'U', n, evec, n, olap1, n, eval, work, 3*n, rwork, info)
      if (info /= 0) stop 'diagonalize_zvo: zhegv failed.'
      deallocate (olap1, work, rwork)
   end subroutine diagonalize_zvo

   subroutine diagonalize_zpeo(eval, mat, olap)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(in)  :: mat(:), olap(:)
      complex, allocatable :: mat1(:), olap1(:), work(:)
      real, allocatable :: rwork(:)
      integer ::                 :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_zpeo: zero dimension in eigenvalue problem.'
      if (size(mat) /= nn) &
         stop 'diagonalize_zpeo: input matrix has wrong size.'
      if (size(olap) /= nn) &
         stop 'diagonalize_zpeo: overlap matrix has wrong size.'
      allocate (mat1(nn), olap1(nn), work(3*n), rwork(3*n))
      mat1 = mat; olap1 = olap
      call zhpgv(1, 'N', 'U', n, mat1, olap1, eval, work, n, work, rwork, info)
      if (info /= 0) stop 'diagonalize_zpeo: zhpev failed.'
      deallocate (mat1, olap1, work, rwork)
   end subroutine diagonalize_zpeo

   subroutine diagonalize_zpvo(evec, eval, mat, olap)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:), olap(:)
      complex, allocatable :: mat1(:), olap1(:), work(:)
      real, allocatable :: rwork(:)
      integer ::                 :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_zpvo: zero dimension in eigenvalue problem.'
      if (size(mat) /= nn) &
         stop 'diagonalize_zpvo: input matrix has wrong size.'
      if (size(olap) /= nn) &
         stop 'diagonalize_zpvo: overlap matrix has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_zpvo: eigenvector array has wrong dimensions.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_zpvo: dimensions of eigenvector array differ.'
      allocate (mat1(nn), olap1(nn), work(3*n), rwork(3*n))
      mat1 = mat; olap1 = olap
      call zhpgv(1, 'V', 'U', n, mat1, olap1, eval, evec, n, work, rwork, info)
      if (info /= 0) stop 'diagonalize_zpvo: zhpgv failed.'
      deallocate (mat1, olap1, work, rwork)
   end subroutine diagonalize_zpvo

   subroutine diagonalize_dvs(evec, eval, mat, m)
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:, :)
      real, allocatable :: work(:), mat1(:, :)
      real ::              :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer ::              :: n, ma, idum, info
      ma = abs(m)
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_dvs: zero dimension in eigenvalue problem.'
      if (ma > n) stop 'diagonalize_dvs: number of selected eigenvalues& &
         exceeds maximal number.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_dvs: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_dvs: eigenvalue array has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_dvs: first dimension of eigenvector is wrong.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_dvs: second dimension of eigenvector is wrong.'
      allocate (work(8*n), iwork(5*n), mat1(n, n), ifail(n)); mat1 = mat
      abstol = 2*dlamch('S')
      if (m > 0) then
         call dsyevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, n - ma + 1, n, abstol, idum, &
                     eval, evec, n, work, 8*n, iwork, ifail, info)
      else
         call dsyevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, 1, ma, abstol, idum, &
                     eval, evec, n, work, 8*n, iwork, ifail, info)
      endif
      if (info /= 0) stop 'diagonalize_dvs: dsyevx failed.'
      deallocate (work, iwork, mat1, ifail)
   end subroutine diagonalize_dvs

   subroutine diagonalize_dvos(evec, eval, mat, olap, m)
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:, :), olap(:, :)
      real, allocatable :: work(:), mat1(:, :), olap1(:, :)
      real ::              :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer ::              :: n, ma, idum, info
      ma = abs(m)
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_dvos: zero dimension in eigenvalue problem.'
      if (ma > n) stop 'diagonalize_dvos: number of selected& &
         eigenvalues exceeds maximal number.'
      if (size(mat, 2) /= n) stop 'diagonalize_dvos: dimensions of& &
         input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_dvos: eigenvalue array has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_dvos: first dimension of eigenvector is wrong.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_dvos: second dimension of eigenvector is wrong.'
      if (size(olap, 1) /= n) stop 'diagonalize_dvos: first dimension & &
         of overlap matrix is wrong.'
      if (size(olap, 2) /= n) stop 'diagonalize_dvos: second dimension of& &
         overlap matrix is wrong.'
      allocate (work(8*n), iwork(5*n), mat1(n, n), olap1(n, n), ifail(n))
      mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (m > 0) then
         call dsygvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, n - ma + 1, n, &
                     abstol, idum, eval, evec, n, work, 8*n, iwork, ifail, info)
      else
         call dsygvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, 1, ma, &
                     abstol, idum, eval, evec, n, work, 8*n, iwork, ifail, info)
      endif
      if (info /= 0) stop 'diagonalize_dvos: dsygvx failed.'
      deallocate (work, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_dvos

   subroutine diagonalize_dpvs(evec, eval, mat, m)
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:)
      real, allocatable :: work(:), mat1(:)
      real ::              :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer ::              :: n, nn, ma, idum, info
      ma = abs(m)
      n = size(eval)
      nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_dpvs: zero dimension in eigenvalue problem.'
      if (ma > n) stop 'diagonalize_dpvs: number of selected & &
         eigenvalues exceeds maximal number.'
      if (size(mat) /= nn) &
         stop 'diagonalize_dpvs: input matrix has wrong size.'
      if (size(evec, 1) /= n) stop 'diagonalize_dpvs: first dimension & &
         of eigenvector is wrong.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_dpvs: second dimension of eigenvector is wrong.'
      allocate (work(8*n), iwork(5*n), mat1(nn), ifail(n)); mat1 = mat
      abstol = 2*dlamch('S')
      if (m > 0) then
         call dspevx('V', 'I', 'U', n, mat1, 0.0, 0.0, n - ma + 1, n, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      else
         call dspevx('V', 'I', 'U', n, mat1, 0.0, 0.0, 1, ma, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      endif
      if (info /= 0) stop 'diagonalize_dpvs: dspevx failed.'
      deallocate (work, iwork, mat1, ifail)
   end subroutine diagonalize_dpvs

   subroutine diagonalize_dpvos(evec, eval, mat, olap, m)
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:), olap(:)
      real, allocatable :: work(:), mat1(:), olap1(:)
      real ::              :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer ::              :: n, nn, ma, idum, info
      ma = abs(m)
      n = size(eval)
      nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_dpvos: zero dimension in eigenvalue problem.'
      if (ma > n) stop 'diagonalize_dpvos: number of selected & &
         eigenvalues exceeds maximal number.'
      if (size(mat) /= nn) &
         stop 'diagonalize_dpvos: input matrix has wrong size.'
      if (size(olap) /= nn) &
         stop 'diagonalize_dpvos: overlap matrix has wrong size.'
      if (size(evec, 1) /= n) stop 'diagonalize_dpvos: first dimension& &
         of eigenvector is wrong.'
      if (size(evec, 2) /= n) stop 'diagonalize_dpvos: second dimension& &
         of eigenvector is wrong.'
      allocate (work(8*n), iwork(5*n), mat1(nn), olap1(nn), ifail(n))
      mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (m > 0) then
         call dspgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, n - ma + 1, n, abstol, &
                     idum, eval, evec, n, work, iwork, ifail, info)
      else
         call dspgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, 1, ma, abstol, &
                     idum, eval, evec, n, work, iwork, ifail, info)
      endif
      if (info /= 0) stop 'diagonalize_dpvos: dspgvx failed.'
      deallocate (work, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_dpvos

   subroutine diagonalize_zvs(evec, eval, mat, m)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:, :)
      complex, allocatable :: work(:), mat1(:, :)
      real, allocatable :: rwork(:)
      real ::                 :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer ::                 :: n, ma, idum, info
      ma = abs(m)
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_zvs: zero dimension in eigenvalue problem.'
      if (ma > n) stop 'diagonalize_zvs: number of selected eigenvalues& &
         exceeds maximal number.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_zvs: dimensions of input matrix differ.'
      if (size(eval) /= n) &
         stop 'diagonalize_zvs: eigenvalue array has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_zvs: first dimension of eigenvector is wrong.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_zvs: second dimension of eigenvector is wrong.'
      allocate (work(2*n), rwork(7*n), iwork(5*n), mat1(n, n), ifail(n))
      mat1 = mat
      abstol = 2*dlamch('S')
      if (m > 0) then
         call zheevx('V', 'I', 'U', n, mat1, n, 00., 0.0, n - ma + 1, n, abstol, idum, &
                     eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      else
         call zheevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, 1, ma, abstol, idum, &
                     eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      endif
      if (info /= 0) stop 'diagonalize_zvs: zheevx failed.'
      deallocate (work, rwork, iwork, mat1, ifail)
   end subroutine diagonalize_zvs

   subroutine diagonalize_zvos(evec, eval, mat, olap, m)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:, :), olap(:, :)
      complex, allocatable :: work(:), mat1(:, :), olap1(:, :)
      real, allocatable :: rwork(:)
      real ::                 :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer ::                 :: n, ma, idum, info
      ma = abs(m)
      n = size(mat, 1)
      if (n == 0) &
         stop 'diagonalize_zvos: zero dimension in eigenvalue problem.'
      if (ma > n) stop 'diagonalize_zvos: number of selected & &
         eigenvalues exceeds maximal number.'
      if (size(mat, 2) /= n) stop 'diagonalize_zvos: dimensions of input& &
         matrix differ.'
      if (size(eval) /= n) stop 'diagonalize_zvos: eigenvalue array& &
         has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_zvos: first dimension of eigenvector is wrong.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_zvos: second dimension of eigenvector is wrong.'
      if (size(olap, 1) /= n) stop 'diagonalize_zvos: first dimension of& &
         overlap matrix is wrong.'
      if (size(olap, 2) /= n) stop 'diagonalize_zvos: second dimension of& &
         overlap matrix is wrong.'
      allocate (work(2*n), rwork(7*n), iwork(5*n), mat1(n, n), olap1(n, n), &
                ifail(n)); mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (m > 0) then
         call zhegvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, n - ma + 1, n, &
                     abstol, idum, eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      else
         call zhegvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, 1, ma, &
                     abstol, idum, eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      endif
      if (info /= 0) stop 'diagonalize_zvos: zhegvx failed.'
      deallocate (work, rwork, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_zvos

   subroutine diagonalize_zpvs(evec, eval, mat, m)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:)
      complex, allocatable :: work(:), mat1(:)
      real, allocatable :: rwork(:)
      real ::                 :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer ::                 :: n, nn, ma, idum, info
      ma = abs(m)
      n = size(eval)
      nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_zpvs: zero dimension in eigenvalue problem.'
      if (ma > n) stop 'diagonalize_zpvs: number of selected& &
         eigenvalues exceeds maximal number.'
      if (size(mat) /= nn) &
         stop 'diagonalize_zpvs: input matrix has wrong size.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_zpvs: first dimension of eigenvector is wrong.'
      if (size(evec, 2) /= n) &
         stop 'diagonalize_zpvs: second dimension of eigenvector is wrong.'
      allocate (work(2*n), rwork(7*n), iwork(5*n), mat1(nn), ifail(n))
      mat1 = mat
      abstol = 2*dlamch('S')
      if (m > 0) then
         call zhpevx('V', 'I', 'U', n, mat1, 0.0, 0.0, n - ma + 1, n, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      else
         call zhpevx('V', 'I', 'U', n, mat1, 0.0, 0.0, 1, ma, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      endif
      if (info /= 0) stop 'diagonalize_zpvs: zhpevx failed.'
      deallocate (work, rwork, iwork, mat1, ifail)
   end subroutine diagonalize_zpvs

   subroutine diagonalize_zpvos(evec, eval, mat, olap, m)
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:), olap(:)
      complex, allocatable :: work(:), mat1(:), olap1(:)
      real, allocatable :: rwork(:)
      real ::                 :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer ::                 :: n, nn, ma, idum, info
      ma = abs(m)
      n = size(eval)
      nn = n*(n + 1)/2
      if (n == 0) &
         stop 'diagonalize_zpvos: zero dimension in eigenvalue problem.'
      if (ma > n) stop 'diagonalize_zpvos: number of& &
         selected eigenvalues exceeds maximal number.'
      if (size(mat) /= nn) &
         stop 'diagonalize_zpvos: input matrix has wrong size.'
      if (size(olap) /= nn) stop 'diagonalize_zpvos: overlap matrix& &
         has wrong size.'
      if (size(evec, 1) /= n) stop 'diagonalize_zpvos: first dimension& &
         of eigenvector is wrong.'
      if (size(evec, 2) /= n) stop 'diagonalize_zpvos: second dimension& &
         of eigenvector is wrong.'
      allocate (work(2*n), rwork(7*n), iwork(5*n), mat1(nn), olap1(nn), &
                ifail(n)); mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (m > 0) then
         call zhpgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, n - ma + 1, n, abstol, &
                     idum, eval, evec, n, work, rwork, iwork, ifail, info)
      else
         call zhpgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, 1, ma, abstol, &
                     idum, eval, evec, n, work, rwork, iwork, ifail, info)
      endif
      if (info /= 0) stop 'diagonalize_zpvos: zhpgvx failed.'
      deallocate (work, rwork, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_zpvos

! routines for diagonalization: eigenvalue range [r1,r2) or index range [ir1,ir2].
! the number of actually found eigenvectors is returned in ir2.

   subroutine diagonalize_dvx(evec, eval, mat, ir1, ir2, r1, r2)
      implicit none
      real, intent(out)   :: eval(:), evec(:, :)
      real, intent(in)    :: mat(:, :)
      real, allocatable   :: work(:), mat1(:, :)
      real ::                :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer ::                :: n, m, idum, info
      n = size(mat, 1)
      m = ir2 - ir1 + 1
      if (n == 0) &
         stop 'diagonalize_dvx: zero dimension in eigenvalue problem.'
      if (m < 0) &
         stop 'diagonalize_dvx: negative index range.'
      if (m > n) stop 'diagonalize_dvx: number of selected& &
         eigenvalues exceeds maximal number.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_dvx: dimensions of input matrix differ.'
      if (size(eval) < m) &
         stop 'diagonalize_dvx: eigenvalue array too small.'
      if (size(evec, 1) /= n) stop 'diagonalize_dvx: first dimension of& &
         eigenvector is wrong.'
      if (size(evec, 2) < m) stop 'diagonalize_dvx: second dimension of& &
         eigenvector too small.'
      allocate (work(8*n), iwork(5*n), mat1(n, n), ifail(n)); mat1 = mat
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call dsyevx('V', 'V', 'U', n, mat1, n, r1, r2, 0, 0, abstol, idum, &
                     eval, evec, n, work, 8*n, iwork, ifail, info)
      else
         call dsyevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, ir1, ir2, abstol, idum, &
                     eval, evec, n, work, 8*n, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) stop 'diagonalize_dvx: dsyevx failed.'
      deallocate (work, iwork, mat1, ifail)
   end subroutine diagonalize_dvx

   subroutine diagonalize_dvox(evec, eval, mat, olap, ir1, ir2, r1, r2)
      implicit none
      real, intent(out)   :: eval(:), evec(:, :)
      real, intent(in)    :: mat(:, :), olap(:, :)
      real, allocatable   :: work(:), mat1(:, :), olap1(:, :)
      real ::                :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer ::                :: n, m, ma, idum, info
      n = size(mat, 1)
      m = ir2 - ir1 + 1
      if (n == 0) &
         stop 'diagonalize_dvox: zero dimension in eigenvalue problem.'
      if (m < 0) &
         stop 'diagonalize_dvox: negative index range.'
      if (m > n) stop 'diagonalize_dvox: number of selected eigenvalues& &
         exceeds maximal number.'
      if (size(mat, 2) /= n) stop 'diagonalize_dvox: dimensions of input& &
         matrix differ.'
      if (size(eval) < m) &
         stop 'diagonalize_dvox: eigenvalue array too small.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_dvox: first dimension of eigenvector is wrong.'
      if (size(evec, 2) < m) stop 'diagonalize_dvox: second dimension of& &
         eigenvector too small.'
      if (size(olap, 1) /= n) stop 'diagonalize_dvox: first dimension of& &
         overlap matrix is wrong.'
      if (size(olap, 2) /= n) stop 'diagonalize_dvox: second dimension of& &
         overlap matrix is wrong.'
      allocate (work(8*n), iwork(5*n), mat1(n, n), olap1(n, n), ifail(n))
      mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call dsygvx(1, 'V', 'V', 'U', n, mat1, n, olap1, n, r1, r2, 0, 0, &
                     abstol, idum, eval, evec, n, work, 8*n, iwork, ifail, info)
      else
         call dsygvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, ir1, ir2, &
                     abstol, idum, eval, evec, n, work, 8*n, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) stop 'diagonalize_dvos: dsygvx failed.'
      deallocate (work, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_dvox

   subroutine diagonalize_dpvx(evec, eval, mat, ir1, ir2, r1, r2)
      implicit none
      real, intent(out)   :: eval(:), evec(:, :)
      real, intent(in)    :: mat(:)
      real, allocatable   :: work(:), mat1(:)
      real ::                :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer ::                :: n, m, nn, idum, info
      n = size(evec, 1)
      nn = n*(n + 1)/2
      m = ir2 - ir1 + 1
      if (n == 0) &
         stop 'diagonalize_dpvx: zero dimension in eigenvalue problem.'
      if (m < 0) &
         stop 'diagonalize_dpvx: negative index range.'
      if (m > n) stop 'diagonalize_dpvx: number of selected& &
         eigenvalues exceeds maximal number.'
      if (size(mat) /= nn) &
         stop 'diagonalize_dpvx: input matrix has wrong size.'
      if (size(eval) < m) &
         stop 'diagonalize_dpvx: eigenvalue array too small.'
      if (size(evec, 2) < m) stop 'diagonalize_dpvx: second dimension of& &
         eigenvector too small.'
      allocate (work(8*n), iwork(5*n), mat1(nn), ifail(n)); mat1 = mat
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call dspevx('V', 'V', 'U', n, mat1, r1, r2, 0, 0, abstol, idum, eval, &
                     evec, n, work, iwork, ifail, info)
      else
         call dspevx('V', 'I', 'U', n, mat1, 0.0, 0.0, ir1, ir2, abstol, idum, eval, &
                     evec, n, work, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) stop 'diagonalize_dpvx: dspevx failed.'
      deallocate (work, iwork, mat1, ifail)
   end subroutine diagonalize_dpvx

   subroutine diagonalize_dpvox(evec, eval, mat, olap, ir1, ir2, r1, r2)
      implicit none
      real, intent(out)   :: eval(:), evec(:, :)
      real, intent(in)    :: mat(:), olap(:)
      real, allocatable   :: work(:), mat1(:), olap1(:)
      real ::                :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer ::                :: n, nn, m, idum, info
      n = size(evec, 1)
      nn = n*(n + 1)/2
      m = ir2 - ir1 + 1
      if (n == 0) &
         stop 'diagonalize_dpvox: zero dimension in eigenvalue problem.'
      if (m < 0) &
         stop 'diagonalize_dpvox: negative index range.'
      if (m > n) stop 'diagonalize_dpvox: number of selected& &
         eigenvalues exceeds maximal number.'
      if (size(mat) /= nn) &
         stop 'diagonalize_dpvox: input matrix has wrong size.'
      if (size(olap) /= nn) &
         stop 'diagonalize_dpvox: overlap matrix has wrong size.'
      if (size(eval) < m) &
         stop 'diagonalize_dpvox: eigenvalue array too small.'
      if (size(evec, 2) < m) stop 'diagonalize_dpvox: second dimension& &
         of eigenvector too small.'
      allocate (work(8*n), iwork(5*n), mat1(nn), olap1(nn), ifail(n))
      mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call dspgvx(1, 'V', 'V', 'U', n, mat1, olap1, r1, r2, 0, 0, abstol, &
                     idum, eval, evec, n, work, iwork, ifail, info)
      else
         call dspgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, ir1, ir2, abstol, &
                     idum, eval, evec, n, work, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) stop 'diagonalize_dpvox: dspgvx failed.'
      deallocate (work, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_dpvox

   subroutine diagonalize_zvx(evec, eval, mat, ir1, ir2, r1, r2)
      implicit none
      real, intent(out)   :: eval(:)
      complex, intent(out)   :: evec(:, :)
      complex, intent(in)    :: mat(:, :)
      complex, allocatable   :: work(:), mat1(:, :)
      real, allocatable   :: rwork(:)
      real ::                   :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer ::                   :: n, m, idum, info
      n = size(mat, 1)
      m = ir2 - ir1 + 1
      if (n == 0) &
         stop 'diagonalize_zvx: zero dimension in eigenvalue problem.'
      if (m < 0) stop 'diagonalize_zvx: negative index range.'
      if (m > n) stop 'diagonalize_zvx: number of selected& &
         eigenvalues exceeds maximal number.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_zvx: dimensions of input matrix differ.'
      if (size(eval) < m) &
         stop 'diagonalize_zvx: eigenvalue array too small.'
      if (size(evec, 1) /= n) &
         stop 'diagonalize_zvx: first dimension of eigenvector is wrong.'
      if (size(evec, 2) < m) &
         stop 'diagonalize_zvx: second dimension of eigenvector too small.'
      allocate (work(2*n), rwork(7*n), iwork(5*n), mat1(n, n), ifail(n))
      mat1 = mat
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call zheevx('V', 'V', 'U', n, mat1, n, r1, r2, 0, 0, abstol, idum, &
                     eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      else
         call zheevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, ir1, ir2, abstol, idum, &
                     eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) stop 'diagonalize_zvx: zheevx failed.'
      deallocate (work, rwork, iwork, mat1, ifail)
   end subroutine diagonalize_zvx

   subroutine diagonalize_zvox(evec, eval, mat, olap, ir1, ir2, r1, r2)
      implicit none
      real, intent(out)   :: eval(:)
      complex, intent(out)   :: evec(:, :)
      complex, intent(in)    :: mat(:, :), olap(:, :)
      complex, allocatable   :: work(:), mat1(:, :), olap1(:, :)
      real, allocatable   :: rwork(:)
      real ::                   :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer ::                   :: n, m, idum, info
      n = size(mat, 1)
      m = ir2 - ir1 + 1
      if (n == 0) &
         stop 'diagonalize_zvox: zero dimension in eigenvalue problem.'
      if (m < 0) &
         stop 'diagonalize_zvox: negative index range.'
      if (m > n) stop 'diagonalize_zvox: number of selected& &
         eigenvalues exceeds maximal number.'
      if (size(mat, 2) /= n) &
         stop 'diagonalize_zvox: dimensions of input matrix differ.'
      if (size(eval) < m) stop 'diagonalize_zvox: eigenvalue array& &
         too small.'
      if (size(evec, 1) /= n) stop 'diagonalize_zvox: first dimension of& &
         eigenvector is wrong.'
      if (size(evec, 2) < m) stop 'diagonalize_zvox: second dimension of& &
         eigenvector too small.'
      if (size(olap, 1) /= n) stop 'diagonalize_zvox: first dimension of& &
         overlap matrix is wrong.'
      if (size(olap, 2) /= n) stop 'diagonalize_zvox: second dimension of& &
         overlap matrix is wrong.'
      allocate (work(2*n), rwork(7*n), iwork(5*n), mat1(n, n), olap1(n, n), &
                ifail(n)); mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call zhegvx(1, 'V', 'V', 'U', n, mat1, n, olap1, n, r1, r2, 0, 0, &
                     abstol, idum, eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      else
         call zhegvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, ir1, ir2, &
                     abstol, idum, eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) stop 'diagonalize_zvox: zhegvx failed.'
      deallocate (work, rwork, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_zvox

   subroutine diagonalize_zpvx(evec, eval, mat, ir1, ir2, r1, r2)
      implicit none
      real, intent(out)   :: eval(:)
      complex, intent(out)   :: evec(:, :)
      complex, intent(in)    :: mat(:)
      complex, allocatable   :: work(:), mat1(:)
      real, allocatable   :: rwork(:)
      real ::                   :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer ::                   :: n, nn, m, idum, info
      n = size(evec, 1)
      nn = n*(n + 1)/2
      m = ir2 - ir1 + 1
      if (n == 0) &
         stop 'diagonalize_zpvx: zero dimension in eigenvalue problem.'
      if (m < 0) &
         stop 'diagonalize_zpvx: negative index range.'
      if (m > n) stop 'diagonalize_zpvx: number of selected& &
         eigenvalues exceeds maximal number.'
      if (size(mat) /= nn) &
         stop 'diagonalize_zpvx: input matrix has wrong size.'
      if (size(eval) < m) &
         stop 'diagonalize_zpvx: eigenvalue array too small.'
      if (size(evec, 2) < m) stop 'diagonalize_zpvx: second dimension& &
         of eigenvector too small.'
      allocate (work(2*n), rwork(7*n), iwork(5*n), mat1(nn), ifail(n))
      mat1 = mat
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call zhpevx('V', 'V', 'U', n, mat1, r1, r2, 0, 0, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      else
         call zhpevx('V', 'I', 'U', n, mat1, 0.0, 0.0, ir1, ir2, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) stop 'diagonalize_zpvx: zhpevx failed.'
      deallocate (work, rwork, iwork, mat1, ifail)
   end subroutine diagonalize_zpvx

   subroutine diagonalize_zpvox(evec, eval, mat, olap, ir1, ir2, r1, r2)
      implicit none
      real, intent(out)   :: eval(:)
      complex, intent(out)   :: evec(:, :)
      complex, intent(in)    :: mat(:), olap(:)
      complex, allocatable   :: work(:), mat1(:), olap1(:)
      real, allocatable   :: rwork(:)
      real ::                   :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer ::                   :: n, nn, m, idum, info
      n = size(evec, 1)
      nn = n*(n + 1)/2
      m = ir2 - ir1 + 1
      if (n == 0) &
         stop 'diagonalize_zpvox: zero dimension in eigenvalue problem.'
      if (m > n) stop 'diagonalize_zpvox: number of selected& &
         eigenvalues exceeds maximal number.'
      if (m < 0) &
         stop 'diagonalize_zpvox: negative index range.'
      if (size(mat) /= nn) &
         stop 'diagonalize_zpvox: input matrix has wrong size.'
      if (size(olap) /= nn) &
         stop 'diagonalize_zpvox: overlap matrix has wrong size.'
      if (size(eval) < m) stop 'diagonalize_zpvox: first dimension& &
         of eigenvector too small.'
      if (size(evec, 2) < m) stop 'diagonalize_zpvox: second dimension& &
         of eigenvector too small.'
      allocate (work(2*n), rwork(7*n), iwork(5*n), mat1(nn), olap1(nn), &
                ifail(n)); mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call zhpgvx(1, 'V', 'V', 'U', n, mat1, olap1, r1, r2, 0, 0, abstol, &
                     idum, eval, evec, n, work, rwork, iwork, ifail, info)
      else
         call zhpgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, ir1, ir2, abstol, &
                     idum, eval, evec, n, work, rwork, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) stop 'diagonalize_zpvox: zhpgvx failed.'
      deallocate (work, rwork, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_zpvox

!     --------

   subroutine geteigen_zpvo(evec, eval, mat, olap)
      implicit none
      real, intent(out) :: eval
      complex, intent(out) :: evec(:)
      complex, intent(in)  :: mat(:), olap(:)
      complex, allocatable :: mat1(:), olap1(:), work(:), evec1(:, :)
      real, allocatable :: eval1(:), rwork(:)
      integer, allocatable :: iwork(:), ifail(:)
      integer ::                 :: n, nn, info
      real ::                 :: dlamch
      stop 'geteigen: disabled!'
      n = size(evec); nn = n*(n + 1)/2
      if (size(mat) /= nn) &
         stop 'diagonalize_zpvo: input matrix has wrong size.'
      if (size(olap) /= nn) &
         stop 'diagonalize_zpvo: overlap matrix has wrong size.'
      allocate (mat1(nn), olap1(nn), eval1(n), evec1(n, n), work(2*n), &
                rwork(7*n), iwork(5*n), ifail(n))
      mat1 = mat; olap1 = olap
!      call zhpgvx(1,'V','I','U',n,mat1,olap1,0.0,0.0,n,n,2*dlamch('S'),1,eval1,evec1,n,work,rwork,iwork,ifail,info)
      if (info /= 0) stop 'geteigen_zpvo: zhpgvx failed.'
      evec = evec1(:, 1)
      eval = eval1(1)
      deallocate (mat1, olap1, eval1, evec1, work, rwork, iwork, ifail)
   end subroutine geteigen_zpvo

!     --------

   subroutine inverse_d(mati, mat)
      implicit none
      real, intent(out) :: mati(:, :)
      real, intent(in)  :: mat(:, :)
      integer ::              :: n, info, i, j
      n = size(mat, 1)
      if (size(mat, 2) /= n) &
         stop 'inverse_d: dimensions of input array differ.'
      if (size(mati, 1) /= n) &
         stop 'inverse_d: output array has wrong dimensions.'
      if (size(mati, 2) /= n) &
         stop 'inverse_d: dimensions of output array differ.'
      mati = mat
      call dpotrf('U', n, mati, n, info)
      if (info /= 0) stop 'inverse_d: dpotrf failed.'
      call dpotri('U', n, mati, n, info)
      if (info /= 0) stop 'inverse_d: dpotri failed.'
      do i = 1, n; do j = 1, i; mati(i, j) = mati(j, i); enddo; enddo
   end subroutine inverse_d

   subroutine inverse_dp(mati, mat)
      implicit none
      real, intent(out) :: mati(:)
      real, intent(in)  :: mat(:)
      integer ::              :: n, nn, info
      nn = size(mat, 1); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      if (size(mati) /= nn) &
         stop 'inverse_dp: output array has wrong size.'
      mati = mat
      call dpptrf('U', n, mati, info)
      if (info /= 0) stop 'inverse_dp: dpptrf failed.'
      call dpptri('U', n, mati, info)
      if (info /= 0) stop 'inverse_dp: dpptri failed.'
   end subroutine inverse_dp

   subroutine inverse_z(mati, mat)
      implicit none
      complex, intent(out) :: mati(:, :)
      complex, intent(in)  :: mat(:, :)
      integer ::                 :: n, info, i, j
      n = size(mat, 1)
      if (size(mat, 2) /= n) &
         stop 'inverse_z: dimensions of input array differ.'
      if (size(mati, 1) /= n) &
         stop 'inverse_z: output array has wrong dimensions.'
      if (size(mati, 2) /= n) &
         stop 'inverse_z: dimensions of output array differ.'
      mati = mat
      call zpotrf('U', n, mati, n, info)
      if (info /= 0) then
         WRITE (*, *) 'info', info
         stop 'inverse_z: zpotrf failed.'
      endif
      call zpotri('U', n, mati, n, info)
      if (info /= 0) stop 'inverse_z: zpotri failed.'
      do i = 1, n
         do j = 1, i
            mati(i, j) = conjg(mati(j, i))
         enddo
      enddo
   end subroutine inverse_z

   subroutine inverse_zp(mati, mat)
      implicit none
      complex, intent(out) :: mati(:)
      complex, intent(in)  :: mat(:)
      integer ::                 :: n, nn, info
      nn = size(mat, 1); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      if (size(mati) /= nn) &
         stop 'inverse_zp: output array has wrong size.'
      mati = mat
      call zpptrf('U', n, mati, n, info)
      if (info /= 0) stop 'inverse_zp: zpptrf failed.'
      call zpptri('U', n, mati, n, info)
      if (info /= 0) stop 'inverse_zp: zpptri failed.'
   end subroutine inverse_zp

   subroutine inverse_d1(mat)
      implicit none
      real, intent(inout) :: mat(:, :)
      real, allocatable   :: work(:)
      integer, allocatable   :: ipiv(:)
      integer ::                :: n, info, i, j
      n = size(mat, 1)
      if (size(mat, 2) /= n) stop 'inverse_d1: array dimensions differ.'
      allocate (ipiv(n), work(n))
      call dgetrf(n, n, mat, n, ipiv, info)
      if (info /= 0) stop 'inverse_d1: dpotrf failed.'
      call dgetri(n, mat, n, ipiv, work, n, info)
      if (info /= 0) stop 'inverse_d1: dpotri failed.'
   end subroutine inverse_d1

   subroutine inverse_dp1(mat)
      implicit none
      real, intent(inout) :: mat(:)
      integer ::                :: n, nn, info
      nn = size(mat, 1); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      call dpptrf('U', n, mat, info)
      if (info /= 0) stop 'inverse_dp1: dpptrf failed.'
      call dpptri('U', n, mat, info)
      if (info /= 0) stop 'inverse_dp1: dpptri failed.'
   end subroutine inverse_dp1

   subroutine inverse_z1(mat)
      implicit none
      complex, intent(inout) :: mat(:, :)
      complex, allocatable   :: work(:)
      integer, allocatable   :: ipiv(:)
      integer ::                   :: n, info
      n = size(mat, 1)
      if (size(mat, 2) /= n) stop 'inverse_z1: array dimensions differ.'
      allocate (ipiv(n), work(n))
      call zgetrf(n, n, mat, n, ipiv, info)
      if (info /= 0) stop 'inverse_z1: zgetrf failed.'
      call zgetri(n, mat, n, ipiv, work, n, info)
      if (info /= 0) stop 'inverse_z1: zgetri failed.'
   end subroutine inverse_z1

   subroutine inverse_zp1(mat)
      implicit none
      complex, intent(inout) :: mat(:)
      complex, allocatable   :: work(:)
      integer, allocatable   :: ipiv(:)
      integer ::                   :: n, nn, info
      nn = size(mat, 1); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      allocate (ipiv(n), work(n))
      call zsptrf('U', n, mat, ipiv, info)
      if (info /= 0) stop 'inverse_zp1: zpptrf failed.'
      call zsptri('U', n, mat, ipiv, work, info)
      if (info /= 0) stop 'inverse_zp1: zpptri failed.'
   end subroutine inverse_zp1

   !     --------

   subroutine sqrtmat_d(matout, matin)
      implicit none
      real, intent(out) :: matout(:, :)
      real, intent(in)  :: matin(:, :)
      real, allocatable :: eval(:), evec(:, :)
      integer ::              :: n, i
      n = size(matin, 1)
      if (size(matin, 2) /= n) &
         stop 'sqrtmat_d: dimensions of input array differ.'
      allocate (evec(n, n), eval(n))
      call diagonalize(evec, eval, matin)
      if (any(eval < 0.0)) stop 'sqrtmat_d: negative eigenvalue.'
      do i = 1, n
         evec(:, i) = sqrt(sqrt(eval(i)))*evec(:, i)
      enddo
      matout = matmul(evec, transpose(evec))
      deallocate (evec, eval)
   end subroutine sqrtmat_d

   subroutine sqrtmat_dp(matout, matin)
      implicit none
      real, intent(out) :: matout(:)
      real, intent(in)  :: matin(:)
      real, allocatable :: eval(:), evec(:, :)
      integer ::              :: nn, n, i, j
      nn = size(matin, 1); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      if (size(matout) /= nn) &
         stop 'sqrtmat_dp: output array has wrong size.'
      allocate (evec(n, n), eval(n))
      call diagonalize(evec, eval, matin)
      if (any(eval < 0.0)) stop 'sqrtmat_dp: negative eigenvalue.'
      do i = 1, n
         evec(:, i) = sqrt(sqrt(eval(i)))*evec(:, i)
      enddo
      matout = packmat(matmul(evec, transpose(evec)))
      deallocate (evec, eval)
   end subroutine sqrtmat_dp

   subroutine sqrtmat_z(matout, matin)
      implicit none
      complex, intent(out) :: matout(:, :)
      complex, intent(in)  :: matin(:, :)
      complex, allocatable :: evec(:, :)
      real, allocatable :: eval(:)
      integer ::                 :: n, i
      n = size(matin, 1)
      if (size(matin, 2) /= n) &
         stop 'sqrtmat_z: dimensions of input array differ.'
      allocate (evec(n, n), eval(n))
      call diagonalize(evec, eval, matin)
      if (any(eval < 0.0)) stop 'sqrtmat_z: negative eigenvalue.'
      do i = 1, n
         evec(:, i) = sqrt(sqrt(eval(i)))*evec(:, i)
      enddo
      matout = matmul(evec, conjg(transpose(evec)))
      deallocate (evec, eval)
   end subroutine sqrtmat_z

   subroutine sqrtmat_zp(matout, matin)
      implicit none
      complex, intent(out) :: matout(:)
      complex, intent(in)  :: matin(:)
      complex, allocatable :: evec(:, :)
      real, allocatable :: eval(:)
      integer ::                 :: nn, n, i, j
      nn = size(matin, 1); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      if (size(matout) /= nn) &
         stop 'sqrtmat_zp: output array has wrong size.'
      allocate (eval(n), evec(n, n))
      call diagonalize(evec, eval, matin)
      if (any(eval < 0.0)) stop 'sqrtmat_zp: negative eigenvalue.'
      do i = 1, n
         evec(:, i) = sqrt(sqrt(eval(i)))*evec(:, i) !sqrt(sqrt(eval(i))) * evec(:,i)
      enddo
      matout = packmat(matmul(evec, conjg(transpose(evec))))
      deallocate (evec, eval)
   end subroutine sqrtmat_zp

   subroutine sqrtmat_d1(mat)
      implicit none
      real, intent(inout) :: mat(:, :)
      real, allocatable   :: eval(:), evec(:, :)
      integer ::                :: n, i
      n = size(mat, 1)
      if (size(mat, 2) /= n) stop 'sqrtmat_d1: array dimensions differ.'
      allocate (evec(n, n), eval(n))
      call diagonalize(evec, eval, mat)
      if (any(eval < 0.0)) stop 'sqrtmat_d1: negative eigenvalue.'
      do i = 1, n
         evec(:, i) = sqrt(sqrt(eval(i)))*evec(:, i)
      enddo
      mat = matmul(evec, transpose(evec))
      deallocate (evec, eval)
   end subroutine sqrtmat_d1

   subroutine sqrtmat_dp1(mat)
      implicit none
      real, intent(inout) :: mat(:)
      real, allocatable   :: eval(:), evec(:, :)
      integer ::                :: nn, n, i, j
      nn = size(mat, 1); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      allocate (evec(n, n), eval(n))
      call diagonalize(evec, eval, mat)
      if (any(eval < 0.0)) stop 'sqrtmat_dp1: negative eigenvalue.'
      do i = 1, n
         evec(:, i) = sqrt(sqrt(eval(i)))*evec(:, i)
      enddo
      mat = packmat(matmul(evec, transpose(evec)))
      deallocate (eval, evec)
   end subroutine sqrtmat_dp1

   subroutine sqrtmat_z1(mat)
      implicit none
      complex, intent(inout) :: mat(:, :)
      complex, allocatable   :: evec(:, :)
      real, allocatable   :: eval(:)
      integer ::                   :: n, i
      n = size(mat, 1)
      if (size(mat, 2) /= n) stop 'sqrtmat_z1: array dimensions differ.'
      allocate (evec(n, n), eval(n))
      call diagonalize(evec, eval, mat)
      if (any(eval < 0.0)) stop 'sqrtmat_z1: negative eigenvalue.'
      do i = 1, n
         evec(:, i) = sqrt(sqrt(eval(i)))*evec(:, i)
      enddo
      mat = matmul(evec, conjg(transpose(evec)))
      deallocate (eval, evec)
   end subroutine sqrtmat_z1

   subroutine sqrtmat_zp1(mat)
      implicit none
      complex, intent(inout) :: mat(:)
      complex, allocatable   :: evec(:, :)
      real, allocatable   :: eval(:)
      integer ::                   :: nn, n, i, j
      nn = size(mat, 1); n = nint(sqrt(0.25 + 2*nn) - 0.5)
      allocate (eval(n), evec(n, n))
      call diagonalize(evec, eval, mat)
      if (any(eval < 0.0)) stop 'sqrtmat_zp1: negative eigenvalue.'
      do i = 1, n
         evec(:, i) = sqrt(sqrt(eval(i)))*evec(:, i)
      enddo
      mat = packmat(matmat(evec, conjg(transpose(evec))))
      deallocate (evec, eval)
   end subroutine sqrtmat_zp1

   !     --------

end module m_wrapper

