module m_wrapper
   interface packmat
      module procedure packmat_d, packmat_z
   end interface

   interface packmatcoul
      module procedure packmatcoul_d, packmatcoul_z
   end interface

   interface unpackmat
      module procedure unpackmat_d, unpackmat_z
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
contains

!     --------

   function identity(n)
      implicit none
      integer, intent(in)  :: n
      integer        :: identity(n, n)
      integer        :: i
      identity = 0
      do i = 1, n
         identity(i, i) = 1
      enddo
   end function identity

!     --------

   function packmat_d(mat)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat(:, :)
      real        :: packmat_d(size(mat, 1)*(size(mat, 1) + 1)/2)
      integer        :: n, nn, i, j, k
      n = size(mat, 1); nn = n*(n + 1)/2
      if (size(mat, 2) /= n) call juDFT_error( 'packmat_d: array dimensions differ.')
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            packmat_d(k) = mat(i, j)
         enddo
      enddo
   end function packmat_d

   function packmatcoul_d(mat)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat(:, :)
      real        :: packmatcoul_d( &
                                      size(mat, 1)*(size(mat, 1) + 1)/2)
      integer        :: n, nn, i, j, k
      n = size(mat, 1); nn = n*(n + 1)/2
      if (size(mat, 2) /= n) call juDFT_error( 'packmat_d: array dimensions differ.')
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1

            packmatcoul_d(k) = (mat(i, j) + mat(j, i))/2.

            !           if(abs(mat(j,i)-mat(i,j)).gt.1e-6) then
            !             write(*,*) 'packmatcoul_d: input matrix not symmetric; deviation .gt. 1E-06'
            !           endif
         enddo
      enddo
   end function packmatcoul_d

   function unpackmat_d(mat)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat(:)
      real        :: unpackmat_d( &
                                      nint(sqrt(0.25 + 2*size(mat)) - 0.5), &
                                      nint(sqrt(0.25 + 2*size(mat)) - 0.5))
      integer        :: n, nn, i, j, k
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
      use m_juDFT
      implicit none
      complex, intent(in) :: mat(:, :)
      complex             :: packmat_z(size(mat, 1)*(size(mat, 1) + 1)/2)
      integer        :: n, nn, i, j, k
      n = size(mat, 1); nn = n*(n + 1)/2
      if (size(mat, 2) /= n) call juDFT_error( 'packmat_z: array dimensions differ.')
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            packmat_z(k) = mat(i, j)
         enddo
      enddo
   end function packmat_z

   function packmatcoul_z(mat)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat(:, :)
      complex              :: packmatcoul_z( &
                              size(mat, 1)*(size(mat, 1) + 1)/2)
      integer        :: n, nn, i, j, k
      n = size(mat, 1); nn = n*(n + 1)/2
      if (size(mat, 2) /= n) call juDFT_error( 'packmat_z: array dimensions differ.')
      k = 0
      do j = 1, n
         do i = 1, j
            k = k + 1
            packmatcoul_z(k) = (mat(i, j) + conjg(mat(j, i)))/2.

            if (abs(conjg(mat(j, i)) - mat(i, j)) > 1e-4) then
               call juDFT_error( 'packmatcoul_z: input matrix not Hermitian; '&
                  // 'deviation > 1E-04.')
            endif
         enddo
      enddo
   end function packmatcoul_z

   function unpackmat_z(mat)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat(:)
      complex              :: unpackmat_z( &
                              nint(sqrt(0.25 + 2*size(mat)) - 0.5), &
                              nint(sqrt(0.25 + 2*size(mat)) - 0.5))
      integer        :: n, nn, i, j, k
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

   function matvec_dpd(mat, vec)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat(:), vec(:)
      real        :: matvec_dpd(size(vec))
      integer        :: nn, n
      n = size(vec)
      nn = n*(n + 1)/2
      if (size(mat) /= nn) call juDFT_error( 'matvec_dpd: input array has wrong size.')
      call dspmv('U', n, 1.0, mat, vec, 1, 0.0, matvec_dpd, 1)
   end function matvec_dpd

   function matvec_dpz(mat, vec)
      use m_juDFT
      implicit none
      real, intent(in) :: mat(:)
      complex, intent(in) :: vec(:)
      complex             :: matvec_dpz(size(vec))
      real, allocatable :: vecr(:), veci(:)
      integer        :: nn, n
      n = size(vec); allocate(vecr(n), veci(n))
      nn = n*(n + 1)/2
      if (size(mat) /= nn) call juDFT_error( 'matvec_dpz: input array has wrong size.')
      call dspmv('U', n, 1.0, mat, real(vec), 1, 0.0, vecr, 1)
      call dspmv('U', n, 1.0, mat, aimag(vec), 1, 0.0, veci, 1)
      matvec_dpz = vecr + (0.0, 1.0)*veci
      deallocate(vecr, veci)
   end function matvec_dpz

   function matvec_zpd(mat, vec)
      use m_juDFT
      implicit none
      complex, intent(in) :: mat(:)
      real, intent(in) :: vec(:)
      complex             :: matvec_zpd(size(vec))
      real, allocatable :: vecr(:), veci(:)
      integer        :: nn, n
      n = size(vec); allocate(vecr(n), veci(n))
      nn = n*(n + 1)/2
      if (size(mat) /= nn) call juDFT_error( 'matvec_zpd: input array has wrong size.')
      call dspmv('U', n, 1.0, real(mat), vec, 1, 0.0, vecr, 1)
      call dspmv('U', n, 1.0, aimag(mat), vec, 1, 0.0, veci, 1)
      matvec_zpd = vecr + (0.0, 1.0)*veci
      deallocate(vecr, veci)
   end function matvec_zpd

   function matvec_zpz(mat, vec)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat(:), vec(:)
      complex              :: matvec_zpz(size(vec))
      integer        :: nn, n
      n = size(vec)
      nn = n*(n + 1)/2
      if (size(mat) /= nn) call juDFT_error( 'matvec_zpz: input array has wrong size.')
      call zhpmv('U', n, (1.0, 0.0), mat, vec, 1, (0.0, 0.0), matvec_zpz, 1)
   end function matvec_zpz

!     --------

   function matmat_dpdp(mat1, mat2)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat1(:), mat2(:)
      real        :: matmat_dpdp( &
                                      nint(sqrt(0.25 + 2*size(mat1)) - 0.5), &
                                      nint(sqrt(0.25 + 2*size(mat1)) - 0.5))
      real, allocatable :: vec(:), vec2(:)
      integer        :: nn, n, k1, i, j, k
      nn = size(mat1)
      n = nint(sqrt(0.25 + 2*nn) - 0.5); allocate(vec(n), vec2(n))
      if (size(mat2) /= nn) &
         call juDFT_error( 'matmat_dpdp: second input array has wrong size.')
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
      deallocate(vec, vec2)
   end function matmat_dpdp

   function matmat_dpzp(mat1, mat2)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat1(:)
      complex, intent(in)  :: mat2(:)
      complex              :: matmat_dpzp( &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5), &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5))
      real, allocatable :: vecr(:), veci(:)
      complex, allocatable :: vec2(:)
      integer        :: nn, n, k1, i, j, k
      nn = size(mat1)
      n = nint(sqrt(0.25 + 2*nn) - 0.5)
      allocate(vecr(n), veci(n), vec2(n))
      if (size(mat2) /= nn) &
         call juDFT_error( 'matmat_dpzp: second input array has wrong size.')
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
      deallocate(vecr, veci, vec2)
   end function matmat_dpzp

   function matmat_zpdp(mat1, mat2)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat1(:)
      real, intent(in)  :: mat2(:)
      complex              :: matmat_zpdp( &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5), &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5))
      real, allocatable :: vecr(:), veci(:)
      complex, allocatable :: vec1(:)
      integer        :: nn, n, k1, i, j, k
      nn = size(mat1)
      n = nint(sqrt(0.25 + 2*nn) - 0.5)
      allocate(vecr(n), veci(n), vec1(n))
      if (size(mat2) /= nn) &
         call juDFT_error( 'matmat_zpdp: second input array has wrong size.')
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
      deallocate(vecr, veci, vec1)
   end function matmat_zpdp

   function matmat_zpzp(mat1, mat2)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat1(:), mat2(:)
      complex              :: matmat_zpzp( &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5), &
                              nint(sqrt(0.25 + 2*size(mat1)) - 0.5))
      complex, allocatable :: vec(:), vec2(:)
      integer        :: nn, n, k1, i, j, k
      nn = size(mat1)
      n = nint(sqrt(0.25 + 2*nn) - 0.5); allocate(vec(n), vec2(n))
      if (size(mat2) /= nn) &
         call juDFT_error( 'matmat_zpzp: second input array has wrong size.')
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
      deallocate(vec, vec2)
   end function matmat_zpzp

   function matmat_dpdm(mat1, mat2)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat1(:), mat2(:, :)
      real        :: matmat_dpdm(size(mat2, 1), size(mat2, 1))
      real, allocatable :: vec(:), vec2(:)
      integer        :: nn, n, i
      n = size(mat2, 1); nn = n*(n + 1)/2; allocate(vec(n), vec2(n))
      if (size(mat2, 2) /= n) &
         call juDFT_error( 'matmat_dpdm: dimensions of second input array differ.')
      if (size(mat1) /= nn) &
         call juDFT_error( 'matmat_dpdm: first input array has wrong size.')
      do i = 1, n
         vec2 = mat2(:, i)
         call dspmv('U', n, 1.0, mat1, vec2, 1, 0.0, vec, 1)
         matmat_dpdm(:, i) = vec
      enddo
      deallocate(vec, vec2)
   end function matmat_dpdm

   function matmat_dmdp(mat1, mat2)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat1(:, :), mat2(:)
      real        :: matmat_dmdp(size(mat1, 1), size(mat1, 1))
      real, allocatable :: vec(:), vec2(:)
      integer        :: nn, n, i
      n = size(mat1, 1); nn = n*(n + 1)/2; allocate(vec(n), vec2(n))
      if (size(mat1, 2) /= n) &
         call juDFT_error( 'matmat_dmdp: dimensions of first input array differ.')
      if (size(mat2) /= nn) &
         call juDFT_error( 'matmat_dmdp: second input array has wrong size.')
      do i = 1, n
         vec2 = mat1(i, :)
         call dspmv('U', n, 1.0, mat2, vec2, 1, 0.0, vec, 1)
         matmat_dmdp(i, :) = vec
      enddo
      deallocate(vec, vec2)
   end function matmat_dmdp

   function matmat_dmdm(mat1, mat2)
      use m_juDFT
      implicit none
      real, intent(in) :: mat1(:, :), mat2(:, :)
      real        :: matmat_dmdm(size(mat1, 1), size(mat1, 1))
      integer        :: n
      n = size(mat1, 1)
      if (size(mat1, 2) /= n) &
         call juDFT_error( 'matmat_dmdm: dimensions of first input array differ.')
      if (size(mat2, 1) /= n) &
         call juDFT_error( 'matmat_dmdm: second input array has wrong dimensions.')
      if (size(mat2, 2) /= n) &
         call juDFT_error( 'matmat_dmdm: dimensions of second input array differ.')
      call dgemm('N', 'N', n, n, n, 1.0, mat1, n, mat2, n, 0.0, matmat_dmdm, n)
   end function matmat_dmdm

   function matmat_dpzm(mat1, mat2)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat1(:)
      complex, intent(in)  :: mat2(:, :)
      complex              :: matmat_dpzm(size(mat2, 1), size(mat2, 1))
      real, allocatable :: vecr(:), veci(:)
      complex, allocatable :: vec2(:)
      integer        :: nn, n, i
      n = size(mat2, 1)
      nn = n*(n + 1)/2; allocate(vecr(n), veci(n), vec2(n))
      if (size(mat2, 2) /= n) &
         call juDFT_error( 'matmat_dpzm: dimensions of second input array differ.')
      if (size(mat1) /= nn) &
         call juDFT_error( 'matmat_dpzm: first input array has wrong size.')
      do i = 1, n
         vec2 = mat2(:, i)
         call dspmv('U', n, 1.0, mat1, real(vec2), 1, 0.0, vecr, 1)
         call dspmv('U', n, 1.0, mat1, aimag(vec2), 1, 0.0, veci, 1)
         matmat_dpzm(:, i) = vecr + (0.0, 1.0)*veci
      enddo
      deallocate(vecr, veci, vec2)
   end function matmat_dpzm

   function matmat_dmzp(mat1, mat2)
      use m_juDFT
      implicit none
      real, intent(in)  :: mat1(:, :)
      complex, intent(in)  :: mat2(:)
      complex              :: matmat_dmzp(size(mat1, 1), size(mat1, 1))
      complex, allocatable :: vec1(:), vec(:)
      integer        :: nn, n, i
      n = size(mat1, 1); nn = n*(n + 1)/2; allocate(vec(n), vec1(n))
      if (size(mat1, 2) /= n) &
         call juDFT_error( 'matmat_dmzp: dimensions of first input array differ.')
      if (size(mat2) /= nn) &
         call juDFT_error( 'matmat_dmzp: second input array has wrong size.')
      do i = 1, n
         vec1 = mat1(i, :)
         call zhpmv('U', n, (1.0, 0.0), mat2, vec1, 1, (0.0, 0.0), vec, 1)
         matmat_dmzp(i, :) = conjg(vec)
      enddo
      deallocate(vec, vec1)
   end function matmat_dmzp

   function matmat_dmzm(mat1, mat2)
      use m_juDFT
      implicit none
      real, intent(in) :: mat1(:, :)
      complex, intent(in) :: mat2(:, :)
      complex             :: matmat_dmzm(size(mat1, 1), size(mat2, 2))
      real        :: matr(size(mat1, 1), size(mat2, 2)), &
                                        mati(size(mat1, 1), size(mat2, 2))
      integer        :: n, n1, n2
      n1 = size(mat1, 1)
      n = size(mat1, 2)
      n2 = size(mat2, 2)
      if (size(mat2, 1) /= n) &
         call juDFT_error( 'matmat_dmzm: dimensions of matrices are inconsistent.')
      call dgemm('N', 'N', n1, n2, n, 1.0, mat1, n1, real(mat2), n, 0.0, matr, n1)
      call dgemm('N', 'N', n1, n2, n, 1.0, mat1, n1, aimag(mat2), n, 0.0, mati, n1)
      matmat_dmzm = matr + (0.0, 1.0)*mati
   end function matmat_dmzm

   function matmat_zpdm(mat1, mat2)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat1(:)
      real, intent(in)  :: mat2(:, :)
      complex              :: matmat_zpdm(size(mat2, 1), size(mat2, 1))
      complex, allocatable :: vec(:), vec2(:)
      integer        :: nn, n, i
      n = size(mat2, 1); nn = n*(n + 1)/2; allocate(vec(n), vec2(n))
      if (size(mat2, 2) /= n) &
         call juDFT_error( 'matmat_zpdm: dimensions of second input array differ.')
      if (size(mat1) /= nn) &
         call juDFT_error( 'matmat_zpdm: first input array has wrong size.')
      do i = 1, n
         vec2 = mat2(:, i)
         call zhpmv('U', n, (1.0, 0.0), mat1, vec2, 1, (0.0, 0.0), vec, 1)
         matmat_zpdm(:, i) = vec
      enddo
      deallocate(vec, vec2)
   end function matmat_zpdm

   function matmat_zmdp(mat1, mat2)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat1(:, :)
      real, intent(in)  :: mat2(:)
      complex              :: matmat_zmdp(size(mat1, 1), size(mat1, 1))
      complex, allocatable :: vec1(:)
      real, allocatable :: vecr(:), veci(:)
      integer        :: nn, n, i
      n = size(mat1, 1); nn = n*(n + 1)/2
      allocate(vecr(n), veci(n), vec1(n))
      if (size(mat1, 2) /= n) &
         call juDFT_error( 'matmat_zmdp: dimensions of first input array differ.')
      if (size(mat2) /= nn) &
         call juDFT_error( 'matmat_zmdp: second input array has wrong size.')
      do i = 1, n
         vec1 = conjg(mat1(i, :))
         call dspmv('U', n, 1.0, mat2, real(vec1), 1, 0.0, vecr, 1)
         call dspmv('U', n, 1.0, mat2, aimag(vec1), 1, 0.0, veci, 1)
         matmat_zmdp(i, :) = vecr - (0.0, 1.0)*veci
      enddo
      deallocate(vecr, veci, vec1)
   end function matmat_zmdp

   function matmat_zmdm(mat1, mat2)
      use m_juDFT
      implicit none
      complex, intent(in) :: mat1(:, :)
      real, intent(in) :: mat2(:, :)
      complex             :: matmat_zmdm(size(mat1, 1), size(mat2, 2))
      real        :: matr(size(mat1, 1), size(mat2, 2)), &
                                        mati(size(mat1, 1), size(mat2, 2))
      integer        :: n, n1, n2
      n1 = size(mat1, 1)
      n = size(mat1, 2)
      n2 = size(mat2, 2)
      if (size(mat2, 1) /= n) &
         call juDFT_error( 'matmat_zmdm: dimensions of matrices are inconsistent.')
      call dgemm('N', 'N', n1, n2, n, 1.0, real(mat1), n1, mat2, n, 0.0, matr, n1)
      call dgemm('N', 'N', n1, n2, n, 1.0, aimag(mat1), n1, mat2, n, 0.0, mati, n1)
      matmat_zmdm = matr + (0.0, 1.0)*mati
   end function matmat_zmdm

   function matmat_zpzm(mat1, mat2)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat1(:), mat2(:, :)
      complex              :: matmat_zpzm(size(mat2, 1), size(mat2, 2))
      complex, allocatable :: vec(:), vec2(:)
      integer        :: nn, n, i, n2
      n = size(mat2, 1); nn = n*(n + 1)/2; allocate(vec(n), vec2(n))
      n2 = size(mat2, 2)
      if (size(mat1) /= nn) &
         call juDFT_error( 'matmat_zpzm: first input array has wrong size.')
      do i = 1, n2
         vec2 = mat2(:, i)
         call zhpmv('U', n, (1.0, 0.0), mat1, vec2, 1, (0.0, 0.0), vec, 1)
         matmat_zpzm(:, i) = vec
      enddo
      deallocate(vec, vec2)
   end function matmat_zpzm

   function matmat_zmzp(mat1, mat2)
      use m_juDFT
      implicit none
      complex, intent(in)  :: mat1(:, :), mat2(:)
      complex              :: matmat_zmzp(size(mat1, 1), size(mat1, 1))
      complex, allocatable :: vec(:), vec2(:)
      integer        :: nn, n, i
      n = size(mat1, 1); nn = n*(n + 1)/2; allocate(vec(n), vec2(n))
      if (size(mat1, 2) /= n) &
         call juDFT_error( 'matmat_zmzp: dimensions of first input array differ.')
      if (size(mat2) /= nn) &
         call juDFT_error( 'matmat_zmzp: second input array has wrong size.')
      do i = 1, n
         vec2 = conjg(mat1(i, :))
         call zhpmv('U', n, (1.0, 0.0), mat2, vec2, 1, (0.0, 0.0), vec, 1)
         matmat_zmzp(i, :) = conjg(vec)
      enddo
      deallocate(vec, vec2)
   end function matmat_zmzp

   function matmat_zmzm(mat1, mat2)
      use m_juDFT
      implicit none
      complex, intent(in) :: mat1(:, :), mat2(:, :)
      complex             :: matmat_zmzm(size(mat1, 1), size(mat2, 2))
      integer        :: n1, n, n2
      complex, parameter     :: one = (1, 0), zero = 0
      n1 = size(mat1, 1)
      n = size(mat1, 2)
      n2 = size(mat2, 2)
      if (size(mat2, 1) /= n) &
         call juDFT_error( 'matmat_zmzm: dimensions of matrices are inconsistent.')
      call zgemm('N', 'N', n1, n2, n, one, mat1, n1, mat2, n, zero, matmat_zmzm, n1)
   end function matmat_zmzm
end module m_wrapper
