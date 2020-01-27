c*********************************************c
c    Check whether uHu matrix is Hermitian    c
c*********************************************c
c                    J.-P. Hanke, Dec. 2015   c
c*********************************************c
      module m_wann_uHu_symcheck
      contains
      subroutine wann_uHu_symcheck(uHu,nbnd,nntot,nntot2,fullnkpts)

      implicit none 
      integer,intent(in) :: nbnd,nntot,nntot2,fullnkpts
      complex,intent(in) :: uHu(nbnd,nbnd,nntot2,nntot,fullnkpts)

      logical :: l_check,l_fdiag,l_foffd
      integer :: i,j,ikpt,ikpt_b,ikpt_b2
      integer :: maxi,maxj,maxb1,maxb2,maxk
      real :: eps,diff,maxdiff

      eps = 1.e-10
      l_check = .true.
      l_fdiag = .false.
      l_foffd = .false.
      maxdiff = 0.0

      open(999,file='out_symcheck')

      do ikpt=1,fullnkpts
         do ikpt_b=1,nntot
            do ikpt_b2=1,nntot2
               do i=1,nbnd
                  do j=1,i
                     diff= abs(uHu(j,i,ikpt_b2,ikpt_b,ikpt)
     >                   -conjg(uHu(i,j,ikpt_b,ikpt_b2,ikpt)))
                     if(diff.gt.eps) then
                      write(999,'(5i7)')j,i,ikpt_b2,ikpt_b,ikpt
                      write(999,'(4f16.10)')
     >                  uHu(j,i,ikpt_b2,ikpt_b,ikpt),diff,
     >                  diff/abs(uHu(j,i,ikpt_b2,ikpt_b,ikpt))*100.0
                      if(diff.gt.maxdiff) then
                         maxj = j
                         maxi = i
                         maxb2= ikpt_b2
                         maxb1= ikpt_b
                         maxk = ikpt
                         maxdiff=diff
                      endif
                      l_check = .false.
                      if(i.eq.j .and. ikpt_b.eq.ikpt_b2) then
                       l_fdiag=.true.
                      else
                       l_foffd=.true.
                      endif
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      close(999)

      if(l_check)write(*,*)'*** uHu_symcheck: fine ***'
      if(.not. l_check) then
         write(*,*)'*** uHu_symcheck: problem ***'
         if(l_fdiag) write(*,*)'              --> fail on diagonal'
         if(l_foffd) write(*,*)'              --> fail on off-diag'
         write(*,*)'              maxdiff:',maxdiff
         write(*,'(a,5(i8,1x))')' index: ',maxj,maxi,maxb2,maxb1,maxk
      endif

      end subroutine wann_uHu_symcheck
      end module m_wann_uHu_symcheck
