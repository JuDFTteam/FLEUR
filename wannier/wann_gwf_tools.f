c************************************c
c  routines to                       c
c      a) find index of (k,q) point  c
c      b) write plot template file   c
c************************************c
      module m_wann_gwf_tools
      implicit none
      contains

      subroutine get_dimension(l_dim,dim)
      implicit none
      logical,intent(in) :: l_dim(3)
      integer,intent(inout) :: dim
      integer :: i

      dim=3
      do i=1,3
         if(l_dim(i))dim=dim+1
      enddo

      end subroutine get_dimension


      subroutine get_shift(l_dim,shift)
      implicit none
      logical,intent(in) :: l_dim(3)
      integer,intent(inout) :: shift(3)

      shift=0
      if(l_dim(1))then
         shift(2)=shift(2)+1
         shift(3)=shift(3)+1
      endif

      if(l_dim(2))then
         shift(3)=shift(3)+1
      endif

      end subroutine get_shift


      integer function get_index_kq(ikpt,iqpt,nkpts)
      implicit none
      integer,intent(in) :: ikpt,iqpt,nkpts
      get_index_kq = ikpt+(iqpt-1)*nkpts
      end function get_index_kq


      integer function get_index_k(ikqpt,nkpts)
      implicit none
      integer,intent(in) :: ikqpt,nkpts
      get_index_k = ikqpt -(get_index_q(ikqpt,nkpts)-1)*nkpts
      end function get_index_k


      integer function get_index_q(ikqpt,nkpts)
      implicit none
      integer,intent(in) :: ikqpt,nkpts
      get_index_q = (ikqpt-1)/nkpts + 1
      end function get_index_q


      integer function get_index_nn_k(bpt,nntot,ikpt_b,gb_kq,gb,dim)
      implicit none
      integer,intent(in) :: nntot,ikpt_b,dim
      integer,intent(in) :: bpt(nntot)
      integer,intent(in) :: gb(3,nntot),gb_kq(dim)
      integer :: nn

c      if(ANY(gb_kq(4:).ne.0)) write(*,*)'problem get_index_nn_k'

      do nn=1,nntot
         if((bpt(nn).eq.ikpt_b) .and. (gb(1,nn).eq.gb_kq(1))
     >            .and.(gb(2,nn).eq.gb_kq(2))
     >            .and.(gb(3,nn).eq.gb_kq(3))) exit
      enddo
      if((nn.eq.(nntot+1)).or.(ANY(gb_kq(4:).ne.0))) then
c       write(*,*)'nn not found!'
       nn=-1
      endif

      get_index_nn_k = nn
      end function get_index_nn_k


      integer function get_index_nn_q(bpt,nntot,ikpt_b,gb_kq,gb,
     >                                dim,shift,l_dim)
      implicit none
      integer,intent(in) :: nntot,ikpt_b,dim,shift(3)
      integer,intent(in) :: bpt(nntot)
      integer,intent(in) :: gb(3,nntot),gb_kq(dim)
      logical,intent(in) :: l_dim(3)
      integer :: nn,ind(3),g(3)

      g = 0
      do nn=1,3
         ind(nn)=4+shift(nn)
         if(l_dim(nn))g(nn)=gb_kq(ind(nn))
      enddo

c      if(any(gb_kq(1:3).ne.0)) write(*,*)'problem get_index_nn_q'

      do nn=1,nntot
         if((bpt(nn).eq.ikpt_b) .and. (gb(1,nn).eq.g(1))
     >              .and.(gb(2,nn).eq.g(2)).and.(gb(3,nn).eq.g(3))) exit
      enddo
      if((nn.eq.(nntot+1)).or.(ANY(gb_kq(1:3).ne.0))) then
c       write(*,*)'nn not found!'
       nn=-1
      endif

      get_index_nn_q = nn
      end function get_index_nn_q


!      integer function get_index_nn_kq(bpt,nntot,kqb,gb_kq,gb,gb_q)
!      implicit none
!      integer,intent(in) :: nntot,kqb
!      integer,intent(in) :: bpt(nntot)
!      integer,intent(in) :: gb_kq(4,nntot),gb(3),gb_q(3)
!      integer :: nn
!
!      do nn=1,nntot
!         if((bpt(nn).eq.kqb) .and. (gb_kq(4,nn).eq.gb_q(3))
!     >                       .and. (gb_kq(3,nn).eq.gb(3))
!     >                       .and. (gb_kq(2,nn).eq.gb(2))
!     >                       .and. (gb_kq(1,nn).eq.gb(1))) exit
!      enddo
!      if(nn==(nntot+1)) stop 'nn not found!'
!
!      get_index_nn_kq = nn
!      end function get_index_nn_kq


      subroutine gwf_plottemplate()
      use m_fleurenv
      implicit none
      integer :: i,nwfs,numbands
      logical :: l_exist

      inquire(file='proj',exist=l_exist)
      if(.not.l_exist) then
         call fleur_err('Where is proj?',
     >                  calledby='gwf_plottemplate')
      endif

      open(8888,file='proj',status='old')
      read(8888,*)nwfs,numbands
      close(8888)

      open(8888,file='printhdwf',status='unknown')
      write(8888,'(i4,3x,a1,3x,a4,i3)')nwfs,'F','nga=',5
      do i=1,nwfs
         write(8888,'(i4,3x,i1,3x,i4)')i,4,0
      enddo

      close(8888)

      write(*,*)'******************************'
      write(*,*)'* created printhdwf template *'
      write(*,*)'******************************'

      end subroutine gwf_plottemplate


      end module m_wann_gwf_tools
