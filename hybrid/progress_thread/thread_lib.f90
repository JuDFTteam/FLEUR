module m_thread_lib
   use iso_c_binding
   implicit none

   interface
      subroutine start_prog_thread (threadId)  bind ( C, name="start_prog_thread" )
         import
         type(c_ptr) :: threadId
      end subroutine start_prog_thread

      subroutine stop_prog_thread (threadId)  bind ( C, name="stop_prog_thread" )
         import
         type(c_ptr) :: threadId
      end subroutine stop_prog_thread
   end interface
end module m_thread_lib
