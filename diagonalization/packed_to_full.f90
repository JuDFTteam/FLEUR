!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Gruenberg Institut, Forschungszentrum Juelich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_packed_to_full
! Contains a service routine to expand a packed storage matrix to full storage
! A real and complex version is provided with a common interface
  use m_juDFT
  implicit none
  private
  interface packed_to_full
     module procedure packed_to_full_r,packed_to_full_c
  end interface packed_to_full
  public packed_to_full
contains
  subroutine packed_to_full_r(n,packed,full)
    integer,intent(in)           :: n
    real,intent(in)              :: packed(:)
    real,allocatable,intent(out) :: full(:,:)

    integer:: i,err,i1,i2

     ALLOCATE ( full(n,n), stat=err )
     if (err/=0) call judft_error("Allocation of full matrix failed",calledby="packed_to_full")
     i=0
     DO i1=1,n
        DO i2=1,i1
           i=i+1
           full(i2,i1)=packed(i)
           full(i1,i2)=packed(i)
        ENDDO
     ENDDO
   end subroutine packed_to_full_r

  subroutine packed_to_full_c(n,packed,full)
    integer,intent(in)           :: n
    complex,intent(in)              :: packed(:)
    complex,allocatable,intent(out) :: full(:,:)

    integer:: i,err,i1,i2

     ALLOCATE ( full(n,n), stat=err )
     if (err/=0) call judft_error("Allocation of full matrix failed",calledby="packed_to_full")
     i=0
     DO i1=1,n
        DO i2=1,i1
           i=i+1
           full(i2,i1)=packed(i)
           full(i1,i2)=conjg(packed(i))
        ENDDO
     ENDDO
   end subroutine packed_to_full_c
 end MODULE m_packed_to_full
