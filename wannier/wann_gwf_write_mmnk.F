c************************c
c  routine to write the  c
c  composite overlaps    c
c************************c
      module m_wann_gwf_write_mmnk
      contains
      subroutine wann_gwf_write_mmnk(mmnk_filename,
     >               nkpts,nntot,
     >               nbnd,bpt,gb,mmnk,l_unformat)

      implicit none
      character(len=*),intent(in) :: mmnk_filename
      integer, intent(in)    :: nkpts,nntot,nbnd
      integer, intent(in)    :: bpt(nntot,nkpts),gb(4,nntot,nkpts)
      complex, intent(inout) :: mmnk(nbnd,nbnd,nntot,nkpts)
      logical, intent(in)    :: l_unformat

      integer :: ikpt,i,j
      integer :: ikpt_b

      if(.not.l_unformat)then

       open (305,file=trim(mmnk_filename)//'.mmn')
       write (305,*) 'Overlaps at neighboring (k,q)-points'
       write (305,'(i5,i7,i5)') nbnd,nkpts,nntot
       do ikpt = 1,nkpts
        do ikpt_b = 1,nntot
         write (305,'(i7,i7,3x,4i4)') ikpt,bpt(ikpt_b,ikpt),
     &                                  gb(1:4,ikpt_b,ikpt)
         do i = 1,nbnd
          do j = 1,nbnd
           write (305,'(2f24.18)')
     &         real(mmnk(j,i,ikpt_b,ikpt)),-aimag(mmnk(j,i,ikpt_b,ikpt))
          enddo
         enddo
        enddo
       enddo !ikpt
       close (305)

      else

       open(305,file=trim(mmnk_filename)//'.mmn',form='unformatted')
       write(305) nbnd,nkpts,nntot
       do ikpt = 1,nkpts
        do ikpt_b = 1,nntot
         write (305) ikpt,bpt(ikpt_b,ikpt),
     &               gb(1:4,ikpt_b,ikpt)
         write (305) ((conjg(mmnk(j,i,ikpt_b,ikpt)),
     &                 j=1,nbnd),i=1,nbnd)
        enddo
       enddo !ikpt
       close (305)

      endif !l_unformat


      end subroutine wann_gwf_write_mmnk

      end module m_wann_gwf_write_mmnk
