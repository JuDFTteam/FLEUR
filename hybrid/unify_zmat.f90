module m_unify_zmat
   use m_types
contains

   subroutine unify_zmat(eigval, z)
      implicit none
      real, intent(in) :: eigval(:)
      type(t_mat)      :: z 

      integer :: beg_group, end_group, n_g
      integer, allocatable :: groups(:)

      groups = make_groups(eigval)
   
      do n_g = 1, size(groups)
         beg_group = sum(groups(1:n_g-1)) + 1
         end_group = sum(groups(1:n_g))
         if(end_group <= z%matsize2) call unify_group(beg_group, end_group, z)
      enddo
   end subroutine unify_zmat

   subroutine unify_group(beg_group, end_group, z)
      use m_legendre_poly
      implicit none
      integer, intent(in)        :: beg_group, end_group
      type(t_mat), intent(inout) :: z 
      type(t_mat)   :: targ, lhs, new_basis
      real, allocatable :: x(:)
      integer :: nrhs

      nrhs = (end_group - beg_group) + 1

      ! leg poly can be symmetric and we don't want that
      x = linspace(-1.1, 1.0, z%matsize1)

      call targ%init(z%l_real, z%matsize1, nrhs)
      call set_sin_targ(targ)

      call lhs%init(z%l_real, z%matsize1, nrhs )
      if(lhs%l_real) then 
         lhs%data_r = z%data_r(:,beg_group:end_group)
      else 
         lhs%data_c = z%data_c(:,beg_group:end_group)
      endif

      call lhs%leastsq(targ)

      if(lhs%l_real) then 
         lhs%data_r = z%data_r(:,beg_group:end_group)
      else 
         lhs%data_c = z%data_c(:,beg_group:end_group)
      endif

      call new_basis%init(lhs) 
      call lhs%multiply(targ, res=new_basis)
      call mod_gram_schmidt(new_basis)

      if(lhs%l_real) then 
         z%data_r(:,beg_group:end_group) = new_basis%data_r
      else 
         z%data_c(:,beg_group:end_group) = new_basis%data_c
      endif
   end subroutine

   subroutine set_sin_targ(targ)
      use m_constants
      implicit none 
      type(t_mat), intent(inout) :: targ 
      integer :: i 
      real, allocatable :: x(:) 

      x = linspace(0.0,4.0, targ%matsize1)
      do i = 1, targ%matsize2 
         if(targ%l_real) then 
            targ%data_r(:,i) = sin(i*x)
         else 
            targ%data_c(:,i) = sin(i*x) + imagunit * cos(i*x)
         endif 
      enddo

   end subroutine set_sin_targ

   function proj_r(u,v) result(p)
      implicit none 
      real, intent(in) :: u(:), v(:)
      real :: p(size(u)) 

      p = dot_product(u,v) / dot_product(u,u) * u 
   end function

   function proj_c(u,v) result(p)
      implicit none 
      complex, intent(in) :: u(:), v(:)
      complex :: p(size(u)) 

      p = dot_product(u,v) / dot_product(u,u) * u 
   end function

   subroutine mod_gram_schmidt(M)
      implicit none 
      type(t_mat), intent(inout) :: M 
      
      integer :: k, i 

      do k = 1,M%matsize2 
         do i = 1,k-1 
            if(M%l_real) then 
               M%data_r(:,k) = M%data_r(:,k) - proj_r(M%data_r(:,i), M%data_r(:,k)) 
            else
               M%data_c(:,k) = M%data_c(:,k) - proj_c(M%data_c(:,i), M%data_c(:,k)) 
            endif 
         enddo 


         if(M%l_real) then
            M%data_r(:,k) = M%data_r(:,k) / norm2(M%data_r(:,k))
         else 
            M%data_c(:,k) = M%data_c(:,k) / norm2(abs(M%data_c(:,k)))
         endif
      enddo 
   end subroutine mod_gram_schmidt

   function linspace(beg, fin, n_steps) result(res)
      implicit none 
      real, intent(in)    :: beg, fin 
      integer, intent(in) :: n_steps
      real                :: res(n_steps), step

      integer :: i

      step = (fin - beg)  / (n_steps - 1.0)

      do i = 0,n_steps-1 
         res(i+1) = beg + step * i 
      enddo 
   end function linspace

   function make_groups(eigval) result(groups)
      implicit none
      real, intent(in)     :: eigval(:)
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
end module m_unify_zmat