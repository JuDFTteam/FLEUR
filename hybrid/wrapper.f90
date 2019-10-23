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
contains

!     --------

   function identity(n)
      implicit none
      integer, intent(in)  :: n
      integer        :: identity(n, n)
      integer        :: i, j
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
      integer        :: nn, n, k1, i
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
      integer        :: nn, n, k1, i
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
      integer        :: nn, n, k1, i
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
      integer        :: nn, n, k1, i
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
      integer        :: nn, n, k1, i
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
      integer        :: nn, n, k1, i
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
      integer        :: nn, n, k1, i, n2
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
      integer        :: nn, n, k1, i
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

!     --------

   subroutine diagonalize_de(eval, mat)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      real, intent(in)  :: mat(:, :)
      real, allocatable :: mat1(:, :), work(:)
      integer        :: n, info
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_de: zero dimension in eigenvalue problem.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_de: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_de: eigenvalue array has wrong size.')
      allocate(mat1(n, n), work(3*n)); mat1 = mat
      call dsyev('N', 'U', n, mat1, n, eval, work, 3*n, info)
      if (info /= 0) call juDFT_error( 'diagonalize_de: dsyev failed.')
      deallocate(mat1, work)
   end subroutine diagonalize_de

   subroutine diagonalize_dv(evec, eval, mat)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:, :)
      real, allocatable :: work(:)
      integer        :: n, info
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_dv: zero dimension in eigenvalue problem.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_dv: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_dv: eigenvalue array has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_dv: eigenvector array has wrong dimensions.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_dv: dimensions of eigenvector array differ.')
      allocate(work(3*n)); evec = mat
      call dsyev('V', 'U', n, evec, n, eval, work, 3*n, info)
      if (info /= 0) call juDFT_error( 'diagonalize_dv: dsyev failed.')
      deallocate(work)
   end subroutine diagonalize_dv

   subroutine diagonalize_dpe(eval, mat)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      real, intent(in)  :: mat(:)
      real, allocatable :: mat1(:), work(:)
      integer        :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_dpe: zero dimension in eigenvalue problem.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_dpe: input matrix has wrong size.')
      allocate(mat1(nn), work(3*n)); mat1 = mat
      call dspev('N', 'U', n, mat1, eval, work, n, work, info)
      if (info /= 0) call juDFT_error( 'diagonalize_dpe: dspev failed.')
      deallocate(mat1, work)
   end subroutine diagonalize_dpe

   subroutine diagonalize_dpv(evec, eval, mat)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:)
      real, allocatable :: mat1(:), work(:)
      integer        :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_dpv: zero dimension in eigenvalue problem.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_dpv: input matrix has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_dpv: eigenvector array has wrong dimensions.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_dpv: dimensions of eigenvector array differ.')
      allocate(mat1(nn), work(3*n)); mat1 = mat
      call dspev('V', 'U', n, mat1, eval, evec, n, work, info)
      if (info /= 0) call juDFT_error( 'diagonalize_dpv: dspev failed.')
      deallocate(mat1, work)
   end subroutine diagonalize_dpv

   subroutine diagonalize_ze(eval, mat)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(in)  :: mat(:, :)
      complex, allocatable :: mat1(:, :), work(:)
      real, allocatable :: rwork(:)
      integer        :: n, info
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_ze: zero dimension in eigenvalue problem.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_ze: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_ze: eigenvalue array has wrong size.')
      allocate(mat1(n, n), work(3*n), rwork(3*n)); mat1 = mat
      call zheev('N', 'U', n, mat1, n, eval, work, 3*n, rwork, info)
      if (info /= 0) call juDFT_error( 'diagonalize_ze: zheev failed.')
      deallocate(mat1, work, rwork)
   end subroutine diagonalize_ze

   subroutine diagonalize_zv(evec, eval, mat)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:, :)
      complex, allocatable :: work(:)
      real, allocatable :: rwork(:)
      integer        :: n, info
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_zv: zero dimension in eigenvalue problem.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_zv: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_zv: eigenvalue array has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_zv: eigenvector array has wrong dimensions.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_zv: dimensions of eigenvector array differ.')
      allocate(work(3*n), rwork(3*n)); evec = mat
      call zheev('V', 'U', n, evec, n, eval, work, 3*n, rwork, info)
      if (info /= 0) call juDFT_error( 'diagonalize_zv: zheev failed.')
      deallocate(work, rwork)
   end subroutine diagonalize_zv

   subroutine diagonalize_zpe(eval, mat)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(in)  :: mat(:)
      complex, allocatable :: mat1(:), work(:)
      real, allocatable :: rwork(:)
      integer        :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_zpe: zero dimension in eigenvalue problem.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_zpe: input matrix has wrong size.')
      allocate(mat1(nn), work(3*n), rwork(3*n)); mat1 = mat
      call zhpev('N', 'U', n, mat1, eval, work, n, work, rwork, info)
      if (info /= 0) call juDFT_error( 'diagonalize_zpe: zhpev failed.')
      deallocate(mat1, work, rwork)
   end subroutine diagonalize_zpe

   subroutine diagonalize_zpv(evec, eval, mat)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:)
      complex, allocatable :: mat1(:), work(:)
      real, allocatable :: rwork(:)
      integer        :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_zpv: zero dimension in eigenvalue problem.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_zpv: input matrix has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_zpv: eigenvector array has wrong dimensions.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_zpv: dimensions of eigenvector array differ.')
      allocate(mat1(nn), work(3*n), rwork(3*n)); mat1 = mat
      call zhpev('V', 'U', n, mat1, eval, evec, n, work, rwork, info)
      if (info /= 0) call juDFT_error( 'diagonalize_zpv: zhpev failed.')
      deallocate(mat1, work, rwork)
   end subroutine diagonalize_zpv

   subroutine diagonalize_deo(eval, mat, olap)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      real, intent(in)  :: mat(:, :), olap(:, :)
      real, allocatable :: mat1(:, :), olap1(:, :), work(:)
      integer        :: n, info
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_deo: zero dimension in eigenvalue problem.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_deo: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_deo: eigenvalue array has wrong size.')
      if (size(olap, 1) /= n) &
         call juDFT_error( 'diagonalize_deo: overlap matrix has wrong size.')
      if (size(olap, 2) /= n) &
         call juDFT_error( 'diagonalize_deo: dimensions of overlap matrix differ.')
      allocate(mat1(n, n), olap1(n, n), work(3*n))
      mat1 = mat; olap1 = olap
      call dsygv(1, 'N', 'U', n, mat1, n, olap1, n, eval, work, 3*n, info)
      if (info /= 0) call juDFT_error( 'diagonalize_deo: dsygv failed.')
      deallocate(mat1, olap1, work)
   end subroutine diagonalize_deo

   subroutine diagonalize_dvo(evec, eval, mat, olap)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:, :), olap(:, :)
      real, allocatable :: olap1(:, :), work(:)
      integer        :: n, info
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_dvo: zero dimension in eigenvalue problem.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_dvo: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_dvo: eigenvalue array has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_dvo: eigenvector array has wrong dimensions.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_dvo: dimensions of eigenvector array differ.')
      if (size(olap, 1) /= n) &
         call juDFT_error( 'diagonalize_dvo: overlap matrix has wrong dimensions.')
      if (size(olap, 2) /= n) &
         call juDFT_error( 'diagonalize_dvo: dimensions of overlap matrix differ.')
      allocate(olap1(n, n), work(3*n)); evec = mat; olap1 = olap
      call dsygv(1, 'V', 'U', n, evec, n, olap1, n, eval, work, 3*n, info)
      if (info /= 0) call juDFT_error( 'diagonalize_dvo: dsygv failed.')
      deallocate(olap1, work)
   end subroutine diagonalize_dvo

   subroutine diagonalize_dpeo(eval, mat, olap)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      real, intent(in)  :: mat(:), olap(:)
      real, allocatable :: mat1(:), olap1(:), work(:)
      integer        :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_dpeo: zero dimension in eigenvalue problem.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_dpeo: input matrix has wrong size.')
      if (size(olap) /= nn) &
         call juDFT_error( 'diagonalize_dpeo: overlap matrix has wrong size.')
      allocate(mat1(nn), olap1(nn), work(3*n))
      mat1 = mat; olap1 = olap
      call dspgv(1, 'N', 'U', n, mat1, olap1, eval, work, n, work, info)
      if (info /= 0) call juDFT_error( 'diagonalize_dpeo: dspgv failed.')
      deallocate(mat1, olap1, work)
   end subroutine diagonalize_dpeo

   subroutine diagonalize_dpvo(evec, eval, mat, olap)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:), olap(:)
      real, allocatable :: mat1(:), olap1(:), work(:)
      integer        :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_dpvo: zero dimension in eigenvalue problem.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_dpvo: input matrix has wrong size.')
      if (size(olap) /= nn) &
         call juDFT_error( 'diagonalize_dpvo: overlap matrix has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_dpvo: eigenvector array has wrong dimensions.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_dpvo: dimensions of eigenvector array differ.')
      allocate(mat1(nn), olap1(nn), work(3*n))
      mat1 = mat; olap1 = olap
      call dspgv(1, 'V', 'U', n, mat1, olap1, eval, evec, n, work, info)
      if (info /= 0) call juDFT_error( 'diagonalize_dpvo: dspgv failed.')
      deallocate(mat1, olap1, work)
   end subroutine diagonalize_dpvo

   subroutine diagonalize_zeo(eval, mat, olap)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(in)  :: mat(:, :), olap(:, :)
      complex, allocatable :: mat1(:, :), olap1(:, :), work(:)
      real, allocatable :: rwork(:)
      integer        :: n, info
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_zeo: zero dimension in eigenvalue problem.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_zeo: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_zeo: eigenvalue array has wrong size.')
      if (size(olap, 1) /= n) &
         call juDFT_error( 'diagonalize_zeo: overlap matrix has wrong size.')
      if (size(olap, 2) /= n) &
         call juDFT_error( 'diagonalize_zeo: dimensions of overlap matrix differ.')
      allocate(mat1(n, n), olap1(n, n), work(3*n), rwork(3*n))
      mat1 = mat; olap1 = olap
      call zhegv(1, 'N', 'U', n, mat1, n, olap1, n, eval, work, 3*n, rwork, info)
      if (info /= 0) call juDFT_error( 'diagonalize_zeo: zhegv failed.')
      deallocate(mat1, olap1, work, rwork)
   end subroutine diagonalize_zeo

   subroutine diagonalize_zvo(evec, eval, mat, olap)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:, :), olap(:, :)
      complex, allocatable :: olap1(:, :), work(:)
      real, allocatable :: rwork(:)
      integer        :: n, info
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_zvo: zero dimension in eigenvalue problem.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_zvo: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_zvo: eigenvalue array has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_zvo: eigenvector array has wrong dimensions.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_zvo: dimensions of eigenvector array differ.')
      if (size(olap, 1) /= n) &
         call juDFT_error( 'diagonalize_zvo: overlap matrix has wrong dimensions.')
      if (size(olap, 2) /= n) &
         call juDFT_error( 'diagonalize_zvo: dimensions of overlap matrix differ.')
      allocate(olap1(n, n), work(3*n), rwork(3*n))
      evec = mat; olap1 = olap
      call zhegv(1, 'V', 'U', n, evec, n, olap1, n, eval, work, 3*n, rwork, info)
      if (info /= 0) call juDFT_error( 'diagonalize_zvo: zhegv failed.')
      deallocate(olap1, work, rwork)
   end subroutine diagonalize_zvo

   subroutine diagonalize_zpeo(eval, mat, olap)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(in)  :: mat(:), olap(:)
      complex, allocatable :: mat1(:), olap1(:), work(:)
      real, allocatable :: rwork(:)
      integer        :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_zpeo: zero dimension in eigenvalue problem.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_zpeo: input matrix has wrong size.')
      if (size(olap) /= nn) &
         call juDFT_error( 'diagonalize_zpeo: overlap matrix has wrong size.')
      allocate(mat1(nn), olap1(nn), work(3*n), rwork(3*n))
      mat1 = mat; olap1 = olap
      call zhpgv(1, 'N', 'U', n, mat1, olap1, eval, work, n, work, rwork, info)
      if (info /= 0) call juDFT_error( 'diagonalize_zpeo: zhpev failed.')
      deallocate(mat1, olap1, work, rwork)
   end subroutine diagonalize_zpeo

   subroutine diagonalize_zpvo(evec, eval, mat, olap)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:), olap(:)
      complex, allocatable :: mat1(:), olap1(:), work(:)
      real, allocatable :: rwork(:)
      integer        :: n, nn, info
      n = size(eval, 1); nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_zpvo: zero dimension in eigenvalue problem.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_zpvo: input matrix has wrong size.')
      if (size(olap) /= nn) &
         call juDFT_error( 'diagonalize_zpvo: overlap matrix has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_zpvo: eigenvector array has wrong dimensions.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_zpvo: dimensions of eigenvector array differ.')
      allocate(mat1(nn), olap1(nn), work(3*n), rwork(3*n))
      mat1 = mat; olap1 = olap
      call zhpgv(1, 'V', 'U', n, mat1, olap1, eval, evec, n, work, rwork, info)
      if (info /= 0) call juDFT_error( 'diagonalize_zpvo: zhpgv failed.')
      deallocate(mat1, olap1, work, rwork)
   end subroutine diagonalize_zpvo

   subroutine diagonalize_dvs(evec, eval, mat, m)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:, :)
      real, allocatable :: work(:), mat1(:, :)
      real        :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer        :: n, ma, idum, info
      ma = abs(m)
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_dvs: zero dimension in eigenvalue problem.')
      if (ma > n) call juDFT_error( 'diagonalize_dvs: number of selected eigenvalues' &
         // ' exceeds maximal number.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_dvs: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_dvs: eigenvalue array has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_dvs: first dimension of eigenvector is wrong.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_dvs: second dimension of eigenvector is wrong.')
      allocate(work(8*n), iwork(5*n), mat1(n, n), ifail(n)); mat1 = mat
      abstol = 2*dlamch('S')
      if (m > 0) then
         call dsyevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, n - ma + 1, n, abstol, idum, &
                     eval, evec, n, work, 8*n, iwork, ifail, info)
      else
         call dsyevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, 1, ma, abstol, idum, &
                     eval, evec, n, work, 8*n, iwork, ifail, info)
      endif
      if (info /= 0) call juDFT_error( 'diagonalize_dvs: dsyevx failed.')
      deallocate(work, iwork, mat1, ifail)
   end subroutine diagonalize_dvs

   subroutine diagonalize_dvos(evec, eval, mat, olap, m)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:, :), olap(:, :)
      real, allocatable :: work(:), mat1(:, :), olap1(:, :)
      real        :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer        :: n, ma, idum, info
      ma = abs(m)
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_dvos: zero dimension in eigenvalue problem.')
      if (ma > n) call juDFT_error( 'diagonalize_dvos: number of selected' &
         // ' eigenvalues exceeds maximal number.')
      if (size(mat, 2) /= n) call juDFT_error( 'diagonalize_dvos: dimensions of ' &
         //'input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_dvos: eigenvalue array has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_dvos: first dimension of eigenvector is wrong.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_dvos: second dimension of eigenvector is wrong.')
      if (size(olap, 1) /= n) call juDFT_error( 'diagonalize_dvos: first dimension '&
         // 'of overlap matrix is wrong.')
      if (size(olap, 2) /= n) call juDFT_error( 'diagonalize_dvos: second dimension of' &
         // 'overlap matrix is wrong.')
      allocate(work(8*n), iwork(5*n), mat1(n, n), olap1(n, n), ifail(n))
      mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (m > 0) then
         call dsygvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, n - ma + 1, n, &
                     abstol, idum, eval, evec, n, work, 8*n, iwork, ifail, info)
      else
         call dsygvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, 1, ma, &
                     abstol, idum, eval, evec, n, work, 8*n, iwork, ifail, info)
      endif
      if (info /= 0) call juDFT_error( 'diagonalize_dvos: dsygvx failed.')
      deallocate(work, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_dvos

   subroutine diagonalize_dpvs(evec, eval, mat, m)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:)
      real, allocatable :: work(:), mat1(:)
      real        :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer        :: n, nn, ma, idum, info
      ma = abs(m)
      n = size(eval)
      nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_dpvs: zero dimension in eigenvalue problem.')
      if (ma > n) call juDFT_error( 'diagonalize_dpvs: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_dpvs: input matrix has wrong size.')
      if (size(evec, 1) /= n) call juDFT_error( 'diagonalize_dpvs: first dimension ' &
         // 'of eigenvector is wrong.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_dpvs: second dimension of eigenvector is wrong.')
      allocate(work(8*n), iwork(5*n), mat1(nn), ifail(n)); mat1 = mat
      abstol = 2*dlamch('S')
      if (m > 0) then
         call dspevx('V', 'I', 'U', n, mat1, 0.0, 0.0, n - ma + 1, n, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      else
         call dspevx('V', 'I', 'U', n, mat1, 0.0, 0.0, 1, ma, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      endif
      if (info /= 0) call juDFT_error( 'diagonalize_dpvs: dspevx failed.')
      deallocate(work, iwork, mat1, ifail)
   end subroutine diagonalize_dpvs

   subroutine diagonalize_dpvos(evec, eval, mat, olap, m)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:), evec(:, :)
      real, intent(in)  :: mat(:), olap(:)
      real, allocatable :: work(:), mat1(:), olap1(:)
      real        :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer        :: n, nn, ma, idum, info
      ma = abs(m)
      n = size(eval)
      nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_dpvos: zero dimension in eigenvalue problem.')
      if (ma > n) call juDFT_error( 'diagonalize_dpvos: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_dpvos: input matrix has wrong size.')
      if (size(olap) /= nn) &
         call juDFT_error( 'diagonalize_dpvos: overlap matrix has wrong size.')
      if (size(evec, 1) /= n) call juDFT_error( 'diagonalize_dpvos: first dimension ' &
         // 'of eigenvector is wrong.')
      if (size(evec, 2) /= n) call juDFT_error( 'diagonalize_dpvos: second dimension ' &
         //'of eigenvector is wrong.')
      allocate(work(8*n), iwork(5*n), mat1(nn), olap1(nn), ifail(n))
      mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (m > 0) then
         call dspgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, n - ma + 1, n, abstol, &
                     idum, eval, evec, n, work, iwork, ifail, info)
      else
         call dspgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, 1, ma, abstol, &
                     idum, eval, evec, n, work, iwork, ifail, info)
      endif
      if (info /= 0) call juDFT_error( 'diagonalize_dpvos: dspgvx failed.')
      deallocate(work, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_dpvos

   subroutine diagonalize_zvs(evec, eval, mat, m)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:, :)
      complex, allocatable :: work(:), mat1(:, :)
      real, allocatable :: rwork(:)
      real        :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer        :: n, ma, idum, info
      ma = abs(m)
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_zvs: zero dimension in eigenvalue problem.')
      if (ma > n) call juDFT_error( 'diagonalize_zvs: number of selected eigenvalues ' &
         // 'exceeds maximal number.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_zvs: dimensions of input matrix differ.')
      if (size(eval) /= n) &
         call juDFT_error( 'diagonalize_zvs: eigenvalue array has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_zvs: first dimension of eigenvector is wrong.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_zvs: second dimension of eigenvector is wrong.')
      allocate(work(2*n), rwork(7*n), iwork(5*n), mat1(n, n), ifail(n))
      mat1 = mat
      abstol = 2*dlamch('S')
      if (m > 0) then
         call zheevx('V', 'I', 'U', n, mat1, n, 00., 0.0, n - ma + 1, n, abstol, idum, &
                     eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      else
         call zheevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, 1, ma, abstol, idum, &
                     eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      endif
      if (info /= 0) call juDFT_error( 'diagonalize_zvs: zheevx failed.')
      deallocate(work, rwork, iwork, mat1, ifail)
   end subroutine diagonalize_zvs

   subroutine diagonalize_zvos(evec, eval, mat, olap, m)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:, :), olap(:, :)
      complex, allocatable :: work(:), mat1(:, :), olap1(:, :)
      real, allocatable :: rwork(:)
      real        :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer        :: n, ma, idum, info
      ma = abs(m)
      n = size(mat, 1)
      if (n == 0) &
         call juDFT_error( 'diagonalize_zvos: zero dimension in eigenvalue problem.')
      if (ma > n) call juDFT_error( 'diagonalize_zvos: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (size(mat, 2) /= n) call juDFT_error( 'diagonalize_zvos: dimensions of input ' &
         // 'matrix differ.')
      if (size(eval) /= n) call juDFT_error( 'diagonalize_zvos: eigenvalue array '&
         // 'has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_zvos: first dimension of eigenvector is wrong.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_zvos: second dimension of eigenvector is wrong.')
      if (size(olap, 1) /= n) call juDFT_error( 'diagonalize_zvos: first dimension of ' &
         // 'overlap matrix is wrong.')
      if (size(olap, 2) /= n) call juDFT_error( 'diagonalize_zvos: second dimension of ' &
         // 'overlap matrix is wrong.')
      allocate(work(2*n), rwork(7*n), iwork(5*n), mat1(n, n), olap1(n, n), &
                ifail(n)); mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (m > 0) then
         call zhegvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, n - ma + 1, n, &
                     abstol, idum, eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      else
         call zhegvx(1, 'V', 'I', 'U', n, mat1, n, olap1, n, 0.0, 0.0, 1, ma, &
                     abstol, idum, eval, evec, n, work, 2*n, rwork, iwork, ifail, info)
      endif
      if (info /= 0) call juDFT_error( 'diagonalize_zvos: zhegvx failed.')
      deallocate(work, rwork, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_zvos

   subroutine diagonalize_zpvs(evec, eval, mat, m)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:)
      complex, allocatable :: work(:), mat1(:)
      real, allocatable :: rwork(:)
      real        :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer        :: n, nn, ma, idum, info
      ma = abs(m)
      n = size(eval)
      nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_zpvs: zero dimension in eigenvalue problem.')
      if (ma > n) call juDFT_error( 'diagonalize_zpvs: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_zpvs: input matrix has wrong size.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_zpvs: first dimension of eigenvector is wrong.')
      if (size(evec, 2) /= n) &
         call juDFT_error( 'diagonalize_zpvs: second dimension of eigenvector is wrong.')
      allocate(work(2*n), rwork(7*n), iwork(5*n), mat1(nn), ifail(n))
      mat1 = mat
      abstol = 2*dlamch('S')
      if (m > 0) then
         call zhpevx('V', 'I', 'U', n, mat1, 0.0, 0.0, n - ma + 1, n, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      else
         call zhpevx('V', 'I', 'U', n, mat1, 0.0, 0.0, 1, ma, abstol, idum, &
                     eval, evec, n, work, iwork, ifail, info)
      endif
      if (info /= 0) call juDFT_error( 'diagonalize_zpvs: zhpevx failed.')
      deallocate(work, rwork, iwork, mat1, ifail)
   end subroutine diagonalize_zpvs

   subroutine diagonalize_zpvos(evec, eval, mat, olap, m)
      use m_juDFT
      implicit none
      real, intent(out) :: eval(:)
      complex, intent(out) :: evec(:, :)
      complex, intent(in)  :: mat(:), olap(:)
      complex, allocatable :: work(:), mat1(:), olap1(:)
      real, allocatable :: rwork(:)
      real        :: abstol, dlamch
      integer, intent(in)  :: m
      integer, allocatable :: iwork(:), ifail(:)
      integer        :: n, nn, ma, idum, info
      ma = abs(m)
      n = size(eval)
      nn = n*(n + 1)/2
      if (n == 0) &
         call juDFT_error( 'diagonalize_zpvos: zero dimension in eigenvalue problem.')
      if (ma > n) call juDFT_error( 'diagonalize_zpvos: number of' &
         // 'selected eigenvalues exceeds maximal number.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_zpvos: input matrix has wrong size.')
      if (size(olap) /= nn) call juDFT_error( 'diagonalize_zpvos: overlap matrix' &
         // 'has wrong size.')
      if (size(evec, 1) /= n) call juDFT_error( 'diagonalize_zpvos: first dimension& &
         of eigenvector is wrong.')
      if (size(evec, 2) /= n) call juDFT_error( 'diagonalize_zpvos: second dimension& &
         of eigenvector is wrong.')
      allocate(work(2*n), rwork(7*n), iwork(5*n), mat1(nn), olap1(nn), &
                ifail(n)); mat1 = mat; olap1 = olap
      abstol = 2*dlamch('S')
      if (m > 0) then
         call zhpgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, n - ma + 1, n, abstol, &
                     idum, eval, evec, n, work, rwork, iwork, ifail, info)
      else
         call zhpgvx(1, 'V', 'I', 'U', n, mat1, olap1, 0.0, 0.0, 1, ma, abstol, &
                     idum, eval, evec, n, work, rwork, iwork, ifail, info)
      endif
      if (info /= 0) call juDFT_error( 'diagonalize_zpvos: zhpgvx failed.')
      deallocate(work, rwork, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_zpvos

! routines for diagonalization: eigenvalue range [r1,r2) or index range [ir1,ir2].
! the number of actually found eigenvectors is returned in ir2.

   subroutine diagonalize_dvx(evec, eval, mat, ir1, ir2, r1, r2)
      use m_juDFT
      implicit none
      real, intent(out)   :: eval(:), evec(:, :)
      real, intent(in)    :: mat(:, :)
      real, allocatable   :: work(:), mat1(:, :)
      real        :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer        :: n, m, idum, info
      n = size(mat, 1)
      m = ir2 - ir1 + 1
      if (n == 0) &
         call juDFT_error( 'diagonalize_dvx: zero dimension in eigenvalue problem.')
      if (m < 0) &
         call juDFT_error( 'diagonalize_dvx: negative index range.')
      if (m > n) call juDFT_error( 'diagonalize_dvx: number of selected ' &
         //'eigenvalues exceeds maximal number.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_dvx: dimensions of input matrix differ.')
      if (size(eval) < m) &
         call juDFT_error( 'diagonalize_dvx: eigenvalue array too small.')
      if (size(evec, 1) /= n) call juDFT_error( 'diagonalize_dvx: first dimension of& &
         eigenvector is wrong.')
      if (size(evec, 2) < m) call juDFT_error( 'diagonalize_dvx: second dimension of& &
         eigenvector too small.')
      allocate(work(8*n), iwork(5*n), mat1(n, n), ifail(n)); mat1 = mat
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call dsyevx('V', 'V', 'U', n, mat1, n, r1, r2, 0, 0, abstol, idum, &
                     eval, evec, n, work, 8*n, iwork, ifail, info)
      else
         call dsyevx('V', 'I', 'U', n, mat1, n, 0.0, 0.0, ir1, ir2, abstol, idum, &
                     eval, evec, n, work, 8*n, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) call juDFT_error( 'diagonalize_dvx: dsyevx failed.')
      deallocate(work, iwork, mat1, ifail)
   end subroutine diagonalize_dvx

   subroutine diagonalize_dvox(evec, eval, mat, olap, ir1, ir2, r1, r2)
      use m_juDFT
      implicit none
      real, intent(out)   :: eval(:), evec(:, :)
      real, intent(in)    :: mat(:, :), olap(:, :)
      real, allocatable   :: work(:), mat1(:, :), olap1(:, :)
      real        :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer        :: n, m, ma, idum, info
      n = size(mat, 1)
      m = ir2 - ir1 + 1
      if (n == 0) &
         call juDFT_error( 'diagonalize_dvox: zero dimension in eigenvalue problem.')
      if (m < 0) &
         call juDFT_error( 'diagonalize_dvox: negative index range.')
      if (m > n) call juDFT_error( 'diagonalize_dvox: number of selected eigenvalues& &
         exceeds maximal number.')
      if (size(mat, 2) /= n) call juDFT_error( 'diagonalize_dvox: dimensions of input& &
         matrix differ.')
      if (size(eval) < m) &
         call juDFT_error( 'diagonalize_dvox: eigenvalue array too small.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_dvox: first dimension of eigenvector is wrong.')
      if (size(evec, 2) < m) call juDFT_error( 'diagonalize_dvox: second dimension of ' &
         // 'eigenvector too small.')
      if (size(olap, 1) /= n) call juDFT_error( 'diagonalize_dvox: first dimension of ' &
         // 'overlap matrix is wrong.')
      if (size(olap, 2) /= n) call juDFT_error( 'diagonalize_dvox: second dimension of' &
         // 'overlap matrix is wrong.')
      allocate(work(8*n), iwork(5*n), mat1(n, n), olap1(n, n), ifail(n))
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
      if (info /= 0) call juDFT_error( 'diagonalize_dvos: dsygvx failed.')
      deallocate(work, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_dvox

   subroutine diagonalize_dpvx(evec, eval, mat, ir1, ir2, r1, r2)
      use m_juDFT
      implicit none
      real, intent(out)   :: eval(:), evec(:, :)
      real, intent(in)    :: mat(:)
      real, allocatable   :: work(:), mat1(:)
      real        :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer        :: n, m, nn, idum, info
      n = size(evec, 1)
      nn = n*(n + 1)/2
      m = ir2 - ir1 + 1
      if (n == 0) &
         call juDFT_error( 'diagonalize_dpvx: zero dimension in eigenvalue problem.')
      if (m < 0) &
         call juDFT_error( 'diagonalize_dpvx: negative index range.')
      if (m > n) call juDFT_error( 'diagonalize_dpvx: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_dpvx: input matrix has wrong size.')
      if (size(eval) < m) &
         call juDFT_error( 'diagonalize_dpvx: eigenvalue array too small.')
      if (size(evec, 2) < m) call juDFT_error( 'diagonalize_dpvx: second dimension of ' &
         // 'eigenvector too small.')
      allocate(work(8*n), iwork(5*n), mat1(nn), ifail(n)); mat1 = mat
      abstol = 2*dlamch('S')
      if (r1 < r2) then
         call dspevx('V', 'V', 'U', n, mat1, r1, r2, 0, 0, abstol, idum, eval, &
                     evec, n, work, iwork, ifail, info)
      else
         call dspevx('V', 'I', 'U', n, mat1, 0.0, 0.0, ir1, ir2, abstol, idum, eval, &
                     evec, n, work, iwork, ifail, info)
      endif
      ir2 = idum
      if (info /= 0) call juDFT_error( 'diagonalize_dpvx: dspevx failed.')
      deallocate(work, iwork, mat1, ifail)
   end subroutine diagonalize_dpvx

   subroutine diagonalize_dpvox(evec, eval, mat, olap, ir1, ir2, r1, r2)
      use m_juDFT
      implicit none
      real, intent(out)   :: eval(:), evec(:, :)
      real, intent(in)    :: mat(:), olap(:)
      real, allocatable   :: work(:), mat1(:), olap1(:)
      real        :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer        :: n, nn, m, idum, info
      n = size(evec, 1)
      nn = n*(n + 1)/2
      m = ir2 - ir1 + 1
      if (n == 0) &
         call juDFT_error( 'diagonalize_dpvox: zero dimension in eigenvalue problem.')
      if (m < 0) &
         call juDFT_error( 'diagonalize_dpvox: negative index range.')
      if (m > n) call juDFT_error( 'diagonalize_dpvox: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_dpvox: input matrix has wrong size.')
      if (size(olap) /= nn) &
         call juDFT_error( 'diagonalize_dpvox: overlap matrix has wrong size.')
      if (size(eval) < m) &
         call juDFT_error( 'diagonalize_dpvox: eigenvalue array too small.')
      if (size(evec, 2) < m) call juDFT_error( 'diagonalize_dpvox: second dimension ' &
         // 'of eigenvector too small.')
      allocate(work(8*n), iwork(5*n), mat1(nn), olap1(nn), ifail(n))
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
      if (info /= 0) call juDFT_error( 'diagonalize_dpvox: dspgvx failed.')
      deallocate(work, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_dpvox

   subroutine diagonalize_zvx(evec, eval, mat, ir1, ir2, r1, r2)
      use m_juDFT
      implicit none
      real, intent(out)   :: eval(:)
      complex, intent(out)   :: evec(:, :)
      complex, intent(in)    :: mat(:, :)
      complex, allocatable   :: work(:), mat1(:, :)
      real, allocatable   :: rwork(:)
      real        :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer        :: n, m, idum, info
      n = size(mat, 1)
      m = ir2 - ir1 + 1
      if (n == 0) &
         call juDFT_error( 'diagonalize_zvx: zero dimension in eigenvalue problem.')
      if (m < 0) call juDFT_error( 'diagonalize_zvx: negative index range.')
      if (m > n) call juDFT_error( 'diagonalize_zvx: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_zvx: dimensions of input matrix differ.')
      if (size(eval) < m) &
         call juDFT_error( 'diagonalize_zvx: eigenvalue array too small.')
      if (size(evec, 1) /= n) &
         call juDFT_error( 'diagonalize_zvx: first dimension of eigenvector is wrong.')
      if (size(evec, 2) < m) &
         call juDFT_error( 'diagonalize_zvx: second dimension of eigenvector too small.')
      allocate(work(2*n), rwork(7*n), iwork(5*n), mat1(n, n), ifail(n))
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
      if (info /= 0) call juDFT_error( 'diagonalize_zvx: zheevx failed.')
      deallocate(work, rwork, iwork, mat1, ifail)
   end subroutine diagonalize_zvx

   subroutine diagonalize_zvox(evec, eval, mat, olap, ir1, ir2, r1, r2)
      use m_juDFT
      implicit none
      real, intent(out)   :: eval(:)
      complex, intent(out)   :: evec(:, :)
      complex, intent(in)    :: mat(:, :), olap(:, :)
      complex, allocatable   :: work(:), mat1(:, :), olap1(:, :)
      real, allocatable   :: rwork(:)
      real        :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer        :: n, m, idum, info
      n = size(mat, 1)
      m = ir2 - ir1 + 1
      if (n == 0) &
         call juDFT_error( 'diagonalize_zvox: zero dimension in eigenvalue problem.')
      if (m < 0) &
         call juDFT_error( 'diagonalize_zvox: negative index range.')
      if (m > n) call juDFT_error( 'diagonalize_zvox: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (size(mat, 2) /= n) &
         call juDFT_error( 'diagonalize_zvox: dimensions of input matrix differ.')
      if (size(eval) < m) call juDFT_error( 'diagonalize_zvox: eigenvalue array& &
         too small.')
      if (size(evec, 1) /= n) call juDFT_error( 'diagonalize_zvox: first dimension of ' &
         // 'eigenvector is wrong.')
      if (size(evec, 2) < m) call juDFT_error( 'diagonalize_zvox: second dimension of ' &
         // 'eigenvector too small.')
      if (size(olap, 1) /= n) call juDFT_error( 'diagonalize_zvox: first dimension of ' &
         // 'overlap matrix is wrong.')
      if (size(olap, 2) /= n) call juDFT_error( 'diagonalize_zvox: second dimension of ' &
         // 'overlap matrix is wrong.')
      allocate(work(2*n), rwork(7*n), iwork(5*n), mat1(n, n), olap1(n, n), &
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
      if (info /= 0) call juDFT_error( 'diagonalize_zvox: zhegvx failed.')
      deallocate(work, rwork, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_zvox

   subroutine diagonalize_zpvx(evec, eval, mat, ir1, ir2, r1, r2)
      use m_juDFT
      implicit none
      real, intent(out)   :: eval(:)
      complex, intent(out)   :: evec(:, :)
      complex, intent(in)    :: mat(:)
      complex, allocatable   :: work(:), mat1(:)
      real, allocatable   :: rwork(:)
      real        :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer        :: n, nn, m, idum, info
      n = size(evec, 1)
      nn = n*(n + 1)/2
      m = ir2 - ir1 + 1
      if (n == 0) &
         call juDFT_error( 'diagonalize_zpvx: zero dimension in eigenvalue problem.')
      if (m < 0) &
         call juDFT_error( 'diagonalize_zpvx: negative index range.')
      if (m > n) call juDFT_error( 'diagonalize_zpvx: number of selected& &
         eigenvalues exceeds maximal number.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_zpvx: input matrix has wrong size.')
      if (size(eval) < m) &
         call juDFT_error( 'diagonalize_zpvx: eigenvalue array too small.')
      if (size(evec, 2) < m) call juDFT_error( 'diagonalize_zpvx: second dimension& &
         of eigenvector too small.')
      allocate(work(2*n), rwork(7*n), iwork(5*n), mat1(nn), ifail(n))
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
      if (info /= 0) call juDFT_error( 'diagonalize_zpvx: zhpevx failed.')
      deallocate(work, rwork, iwork, mat1, ifail)
   end subroutine diagonalize_zpvx

   subroutine diagonalize_zpvox(evec, eval, mat, olap, ir1, ir2, r1, r2)
      use m_juDFT
      implicit none
      real, intent(out)   :: eval(:)
      complex, intent(out)   :: evec(:, :)
      complex, intent(in)    :: mat(:), olap(:)
      complex, allocatable   :: work(:), mat1(:), olap1(:)
      real, allocatable   :: rwork(:)
      real        :: abstol, dlamch
      real, intent(in)    :: r1, r2
      integer, intent(in)    :: ir1
      integer, intent(inout) :: ir2
      integer, allocatable   :: iwork(:), ifail(:)
      integer        :: n, nn, m, idum, info
      n = size(evec, 1)
      nn = n*(n + 1)/2
      m = ir2 - ir1 + 1
      if (n == 0) &
         call juDFT_error( 'diagonalize_zpvox: zero dimension in eigenvalue problem.')
      if (m > n) call juDFT_error( 'diagonalize_zpvox: number of selected ' &
         // 'eigenvalues exceeds maximal number.')
      if (m < 0) &
         call juDFT_error( 'diagonalize_zpvox: negative index range.')
      if (size(mat) /= nn) &
         call juDFT_error( 'diagonalize_zpvox: input matrix has wrong size.')
      if (size(olap) /= nn) &
         call juDFT_error( 'diagonalize_zpvox: overlap matrix has wrong size.')
      if (size(eval) < m) call juDFT_error( 'diagonalize_zpvox: first dimension ' &
         // 'of eigenvector too small.')
      if (size(evec, 2) < m) call juDFT_error( 'diagonalize_zpvox: second dimension ' &
         // 'of eigenvector too small.')
      allocate(work(2*n), rwork(7*n), iwork(5*n), mat1(nn), olap1(nn), &
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
      if (info /= 0) call juDFT_error( 'diagonalize_zpvox: zhpgvx failed.')
      deallocate(work, rwork, iwork, mat1, olap1, ifail)
   end subroutine diagonalize_zpvox
end module m_wrapper
