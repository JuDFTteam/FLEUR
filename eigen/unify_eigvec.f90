module m_unique_eigvec
contains
   subroutine unique_eigvec(eigval, z, ierr)
      use omp_lib
      implicit none
      real(kind=8), intent(in) :: eigval(:)
      real(kind=8)             :: z(:,:)
      integer, intent(out)     :: ierr


      integer :: beg_group, end_group, n_g
      integer, allocatable :: groups(:)

      groups = make_groups(eigval)
   
      do n_g = 1, size(groups)
         beg_group = sum(groups(1:n_g-1)) + 1
         end_group = sum(groups(1:n_g))

         if(end_group <= size(z,2))then 
            call unify_group_operator(beg_group, end_group, z, ierr)
         endif
      enddo
   end subroutine unique_eigvec

   subroutine unify_group_operator(beg_group, end_group, z, ierr)
      implicit none
      integer, intent(in)         :: beg_group, end_group
      real(kind=8), intent(inout) :: z(:,:) 
      integer, intent(out)        :: ierr
      real(kind=8), allocatable   :: mtx(:,:), eigval(:), work(:), rq_mat(:,:), tau(:), Q(:,:), tmp(:,:)
      integer :: dim, lwork, info, i, j
      real(kind=8) :: work_size(1)
      real(kind=8), parameter :: zero = 0.0, one = 1.0

      dim = end_group - beg_group + 1 
      allocate(rq_mat(dim, dim), tau(dim),Q(dim,dim))
      call make_rq_mat(z(:,beg_group:end_group), rq_mat, ierr)
      if(ierr /= 0) return 
      
      rq_mat = transpose(rq_mat)

      call dgeqrfp(dim, dim, rq_mat, dim, tau, work_size, -1, info)
      if(info /= 0) then
         write (*,*) "Problem setting up dgeqrfp"
         ierr = 2
         return
      endif
      lwork = int(work_size(1))
      allocate(work(lwork))
      call dgeqrfp(dim, dim, rq_mat, dim, tau, work, lwork, info)
      if(info /= 0) then
         write (*,*) "Problem executing dgeqrfp"
         ierr = 4
         return
      endif
      deallocate(work)

      q = rq_mat
 
      do i =1,3 
         do j=i,dim
            q(i,j) = 0.0
         enddo
      enddo


      call dorgqr(dim, dim, dim, Q, dim, tau, work_size, -1, info)
      if(info /= 0) then
         write (*,*) "Problem setting up dorgqr"
         ierr = 8
         return
      endif
      lwork = int(work_size(1))
      allocate(work(lwork))
      call dorgqr(dim, dim, dim, Q, dim, tau, work, lwork, info)
      if(info /= 0) then
         write (*,*) "Problem executing dorgqr"
         ierr = 16
         return
      endif
      tmp = z(:,beg_group:end_group)
      ! z(:,beg_group:end_group) = matmul(z(:,beg_group:end_group), q)
      call dgemm("N", "N", size(z,1), dim, dim, one, tmp, size(z,1), q, dim, zero, z(:,beg_group:end_group), size(z,1))
   end subroutine unify_group_operator

   function make_groups(eigval) result(groups)
      implicit none
      real(8), intent(in)  :: eigval(:)
      integer, allocatable :: groups(:)
      integer              :: color(size(eigval))

      integer :: i, beg, g_cnt, n_groups

      g_cnt = 1
      beg  = 1 

      do while (beg <= size(eigval)) 
         i = beg
         do while (abs(eigval(i) - eigval(beg)) < 1e-8  ) 
            color(i) = g_cnt
            i = i + 1 
            if(i > size(eigval)) exit
         enddo 
         g_cnt = g_cnt + 1 
         beg = i 
      enddo


      n_groups = color(size(eigval))
      allocate(groups(n_groups))

      do i = 1,n_groups 
         groups(i) = count(i == color)
      enddo 
   end function make_groups

   subroutine make_rq_mat(eigvecs, rq_mat, ierr)
      implicit none
      real(kind=8), intent(in)                 :: eigvecs(:,:)
      real(kind=8), intent(inout), allocatable :: rq_mat(:,:)
      integer, intent(out) :: ierr
      integer :: n, i, j, best_j
      logical :: l_full_search
      real(kind=8) :: cutoff, lindep, best_lindep
      real(kind=8), allocatable :: tmp(:,:)

      n = size(eigvecs,2)
      cutoff = 1e-6 * sqrt(1.0/size(eigvecs,1))
      if(allocated(rq_mat)) deallocate(rq_mat)
      allocate(rq_mat(n,n), tmp(n,n))

      rq_mat = 0.0
      lindep = 0.0  
      j = 0
      do i = 1,n 
         if(allocated(tmp)) deallocate(tmp)
         allocate(tmp(i,n))
         tmp = 0.0
         lindep = 0.0
         l_full_search = .False.

         best_j = 0
         best_lindep = 0.0
         do while(lindep < cutoff)
            j = j + 1
            if(j > size(eigvecs,1)) then
               rq_mat(i,:) = eigvecs(best_j,:)
               if(l_full_search) then
                  exit
               else 
                  j = 1 
                  l_full_search = .True.
               endif
            endif
            
            rq_mat(i,:) = eigvecs(j,:)
            tmp = rq_mat(1:i,:)
            lindep = linear_independency(tmp)
            if(lindep > best_lindep) then
               best_j = j 
               best_lindep = lindep 
            endif
         enddo
      enddo
      if(lindep < 1e-9) then
         write (*,*) "RQ matrix seems linear dependent", lindep
         ierr = 1
      else
         write (*,*) "lindep", lindep
         ierr = 0
      endif
   end subroutine make_rq_mat

   function linear_independency(mat)
      implicit none
      real(kind=8), intent(inout) :: mat(:,:)
      real(kind=8) :: linear_independency
      
      integer :: info, lwork, iwork(8*minval(shape(mat))), ldmat, m, n
      integer, parameter :: ldu = 1, ldvt = 1

      real(kind=8) :: s(minval(shape(mat))), u(1,1), vt(1,1), work_size(1)
      real(kind=8), allocatable :: work(:)

      ldmat = size(mat, 1)
      m = size(mat,1)
      n = size(mat,2)

      !call dgesdd(jobz,m,n, a,   lda, s, u, ldu, vt, ldvt, work,lwork, iwork, info)
      call dgesdd("N", m, n, mat, ldmat, s, u, ldu, vt, ldvt, work_size, -1, iwork, info)
      
      lwork = int(work_size(1))
      allocate(work(lwork))

      call dgesdd("N", m, n, mat, ldmat, s, u, ldu, vt, ldvt, work, lwork, iwork, info)

      linear_independency = s(size(s))
   end function linear_independency


   subroutine print_mtx(mtx)
      implicit none 
      real(kind=8), intent(in) :: mtx(:,:)
      integer ::i

      do i = 1, size(mtx,1) 
         write (*,*) mtx(i,:)
      enddo
      write (*,*) "####"
   end subroutine print_mtx 
end module m_unique_eigvec