!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_elemental
   PRIVATE
   !Module to call the elemental library for parallel diagonalization
   !complex only version at present!
#ifdef CPP_ELEMENTAL      
! This is the interface defined in fleur_elemental.cpp
interface 
 subroutine fl_el_initialize(n,hmat,smat,mpi_com) bind (c)
   use, intrinsic :: iso_c_binding
   integer(kind=c_int),value                              :: n,mpi_com
   complex(kind=c_double_complex),dimension(*),intent(in) :: hmat,smat
 end subroutine
end interface

interface
 subroutine fl_el_diagonalize(neig,direct,nex,deg,tol,mode,opt,z) bind (c)
   use, intrinsic :: iso_c_binding
   integer(kind=c_int),value :: neig,direct,nex,deg,mode,opt
   real(kind=c_double),value :: tol
   complex(kind=c_double_complex),dimension(*),intent(in):: z
 end subroutine
end interface

interface
 subroutine fl_el_eigenvalues(neig,eig) bind (c)
   use, intrinsic :: iso_c_binding
   integer(kind=c_int),value                   :: neig
   real(kind=c_double),dimension(*),intent(out):: eig
 end subroutine
end interface

interface
 subroutine fl_el_eigenvectors(neig,eig,z) bind (c)
   use, intrinsic :: iso_c_binding
   integer(kind=c_int),value                             :: neig
   real(kind=c_double),dimension(*),intent(out)           :: eig
   complex(kind=c_double_complex),dimension(*),intent(out):: z
 end subroutine
end interface
#endif
 PUBLIC elemental

 CONTAINS

 SUBROUTINE elemental(m,n,SUB_COMM,a,b,z,eig,num,direct)
! 
!----------------------------------------------------
!- Parallel eigensystem solver - driver routine based on chani; dw'12
!
! m ........ actual (=leading) dimension of full a & b matrices
!            must be problem size, as input a, b  are one-dimensional
!            and shall be redistributed to two-dimensional matrices
!            actual (=leading) dimension of eigenvector z(,)
! n ........ number of columns of full (sub)matrix ( about n/np)
! SUB_COMM.. communicator for MPI
! a,b   .... packed (sub)matrices, here expanded to non-packed
! z,eig .... eigenvectors and values, output
! num ...... number of ev's searched (and found) on this node
!            On input, overall number of ev's searched,
!            On output, local number of ev's found
!
!----------------------------------------------------
!
   use m_juDFT
   IMPLICIT NONE
   INTEGER, INTENT (IN)                  :: m,n,direct
   INTEGER, INTENT (IN)                  :: SUB_COMM
   INTEGER, INTENT (INOUT)               :: num
   REAL,    INTENT   (OUT)               :: eig(:)
   COMPLEX, ALLOCATABLE, INTENT (INOUT)  :: a(:),b(:)
   COMPLEX,              INTENT (INOUT)  :: z(:,:)

   INTEGER ::neig,nex,deg,mode,opt
   REAL    :: tol
   INTEGER :: isize,ierr
#ifdef CPP_ELEMENTAL
   !Initialize the matrices in elemental
   CALL fl_el_initialize(m,a,b,SUB_COMM)
   DEALLOCATE(a,b)

   !call diagonalization
   neig=num
   nex=40
   deg=15
   tol=1E-8
   mode=0 !random vectors!?
   opt=1 !opt multiple

   print*, "Elemental: seeking",neig," eigenvalues with direct:",direct
   call fl_el_diagonalize(neig,direct,nex,deg,tol,mode,opt,z(:m,:))

   CALL MPI_COMM_SIZE(SUB_COMM,isize,ierr)
   num=min(size(z,2),neig/isize)
   print*, "Elemental provides ",num," local eigenvalues"

   call fl_el_eigenvectors(num,eig,z(:m,:))

#endif
  END subroutine
  END module
