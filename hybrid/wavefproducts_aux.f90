module m_wavefproducts_aux

CONTAINS
   subroutine prep_list_of_gvec(lapw, lapw_nkqpt, hybrid,g_bounds, g_t, iq,jsp, pointer,gpt0, ngpt0)
      use m_types
      use m_juDFT
      implicit none
      type(t_lapw),         intent(in)    :: lapw, lapw_nkqpt
      type(t_hybrid),       intent(in)    :: hybrid
      integer,              intent(in)    :: g_bounds(3), g_t(3), iq, jsp
      integer, allocatable, intent(inout) :: pointer(:,:,:), gpt0(:,:)
      integer,              intent(inout) :: ngpt0

      integer :: ic, ig1, igptm, iigptm, ok, g(3)


      ALLOCATE (pointer(-g_bounds(1):g_bounds(1), &
                        -g_bounds(2):g_bounds(2),&
                        -g_bounds(3):g_bounds(3)), stat=ok)
      IF (ok /= 0) call juDFT_error('wavefproducts_noinv2: error allocation pointer')
      ALLOCATE (gpt0(3, size(pointer)), stat=ok)
      IF (ok /= 0) call juDFT_error('wavefproducts_noinv2: error allocation gpt0')

      call timestart("prep list of Gvec")
      pointer = 0
      ic = 0
      DO ig1 = 1, lapw%nv(jsp)
         DO igptm = 1, hybrid%ngptm(iq)
            iigptm = hybrid%pgptm(igptm, iq)
            g = lapw%gvec(:,ig1,jsp) + hybrid%gptm(:, iigptm) - g_t
            IF (pointer(g(1), g(2), g(3)) == 0) THEN
               ic = ic + 1
               gpt0(:, ic) = g
               pointer(g(1), g(2), g(3)) = ic
            END IF
         END DO
      END DO
      ngpt0 = ic
      call timestop("prep list of Gvec")
   end subroutine prep_list_of_gvec

   function calc_number_of_basis_functions(lapw, atoms, noco) result(nbasfcn)
      use m_types
      implicit NONE
      type(t_lapw), intent(in)  :: lapw
      type(t_atoms), intent(in) :: atoms
      type(t_noco), intent(in)  :: noco
      integer                   :: nbasfcn

      if(noco%l_noco) then
         nbasfcn = lapw%nv(1) + lapw%nv(2) + 2*atoms%nlotot
      else
         nbasfcn = lapw%nv(1) + atoms%nlotot
      endif
   end function calc_number_of_basis_functions
end module m_wavefproducts_aux
