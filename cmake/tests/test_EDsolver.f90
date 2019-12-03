! ==========================================================================
! Example calls of the exact-diagonalization DMFT solver. The files
! selfen_from_bathGF_diag.dat and selfen_from_cfg.dat produced by this
! program should be identical. See module EDsolver for details.
! --------------------------------------------------------------------------
! Copyright (c) J. Kolorenc <kolorenc@fzu.cz> 2013, 2014
! All rights reserved.
! ==========================================================================

program EDdemo
  use EDsolver, only: EDsolver_from_bathGF_diag, &
       EDsolver_from_bathGF_offdiag, EDsolver_from_cfg
  implicit none

  integer, parameter :: dp=kind(1.0d0)
  integer, parameter :: dpc=kind((1.0d0,1.0d0))

  integer, parameter :: Ndos=1001
  integer, parameter :: maxNbpar=3
  integer, parameter :: bathMethod=1     ! large-z assymptotics
  integer, parameter :: rediagonalize=0  ! no diagonalization in the second call
  real(dp), parameter :: Emin =-15.0_dp
  real(dp), parameter :: Emax =  5.0_dp
  real(dp), parameter :: eps  =  0.1_dp
  real(dp), parameter :: eDC  = 58.0_dp

  complex(dpc), dimension(:,:,:), allocatable :: G0, G0save, selfen
  complex(dpc), dimension(:), allocatable :: z, wt
  real(dp), dimension(:,:), allocatable :: re, im
  real(dp), dimension(:,:), allocatable :: bpar
  integer, dimension(:), allocatable :: Nbpar
  real(dp) :: a, b, c, d
  integer :: Ncirc, dim
  integer :: i
  real(dp) :: time0, time1
  real(dp), dimension(3) :: time

  ! precomputed bath Green's function on a semicircle in the complex plane
  open(unit=11,file="G0circ32.dat",action="read")
  read(unit=11,fmt=*) Ncirc, dim

  allocate( G0(dim,dim,Ncirc+Ndos), G0save(dim,dim,Ncirc), &
       selfen(dim,dim,Ncirc+Ndos), z(Ncirc+Ndos), wt(Ncirc), &
       Nbpar(dim), bpar(maxNbpar,dim), re(dim,dim), im(dim,dim) )
  Nbpar=3

  do i=1, Ncirc
     read(unit=11,fmt=*) a, b, c, d
     z(i)=cmplx(a,b,dpc)
     wt(i)=cmplx(c,d,dpc)
     read(unit=11,fmt=*) re
     read(unit=11,fmt=*) im
     G0(:,:,i)=cmplx(re,im,dpc)
  end do

  close(unit=11)

  ! the energies where we want the selfenergy to be computed
  do i=1, Ndos
     z(Ncirc+i)=cmplx( Emin + (Emax-Emin)/(Ndos-1)*(i-1), eps, dpc )
  end do

  ! EDsolver_from_bathGF_diag() first constructs the discrete impurity model
  ! from the supplied bath Green's function and then calculates the selfenergy
  G0save=G0(:,:,1:Ncirc)    ! G0 is destroyed by the subroutine
  call cpu_time(time0)
  call EDsolver_from_bathGF_diag(dim,Ncirc+Ndos,z,G0,Ncirc,wt,1,1,bathMethod, &
       0.0_dp,Nbpar,bpar,eDC,selfen)
  call cpu_time(time1)
  time(1)=time1-time0
  call write_selfen(z(Ncirc+1:Ncirc+Ndos),selfen(:,:,Ncirc+1:Ncirc+Ndos), &
       "selfen_from_bathGF_diag.dat")

  ! EDsolver_from_cfg() takes an existing impurity model from configuration
  ! files stored on the disk and then calculates the selfenergy
  call cpu_time(time0)
  call EDsolver_from_cfg( dim, Ndos, z(Ncirc+1:Ncirc+Ndos), &
       selfen(:,:,Ncirc+1:Ncirc+Ndos), rediagonalize )
  call cpu_time(time1)
  time(2)=time1-time0
  call write_selfen(z(Ncirc+1:Ncirc+Ndos),selfen(:,:,Ncirc+1:Ncirc+Ndos), &
       "selfen_from_cfg.dat")

  ! EDsolver_from_bathGF_offdiag() first constructs the discrete impurity model
  ! from the supplied bath Green's function and then calculates the selfenergy
  G0(:,:,1:Ncirc)=G0save    ! restore the bath Green's function
  call cpu_time(time0)
  call EDsolver_from_bathGF_offdiag(dim,Ncirc+Ndos,z,G0,Ncirc,wt, &
       Nbpar(1),bpar,eDC,0,selfen)
  call cpu_time(time1)
  time(3)=time1-time0
  call write_selfen(z(Ncirc+1:Ncirc+Ndos),selfen(:,:,Ncirc+1:Ncirc+Ndos), &
       "selfen_from_bathGF_offdiag.dat")
  call write_selfen_full_mtrx(z(Ncirc+1:Ncirc+Ndos), &
       selfen(:,:,Ncirc+1:Ncirc+Ndos), "selfen_full_mtrx.dat")

  ! from_cfg is faster because the diagonalization was not performed again,
  ! from_bathGF_offdiag is a lot slower because it does not explore any
  ! symmetry
  print *
  print '(1x,a,f9.4)', "time in EDsolver_from_bathGF_diag()   :", time(1)
  print '(1x,a,f9.4)', "time in EDsolver_from_cfg()           :", time(2)
  print '(1x,a,f9.4)', "time in EDsolver_from_bathGF_offdiag():", time(3)

contains

  subroutine write_selfen(z,selfen,filename)
    ! {{{
    complex(dpc), dimension(:), intent(in) :: z
    complex(dpc), dimension(:,:,:), intent(in) :: selfen
    character(len=*), intent(in) :: filename
    integer :: i
    open(unit=11,file=filename,action="write")
    do i=1, size(z)
       write(unit=11,fmt="(5f15.7)") real(z(i),dp), &
            real(selfen(3,3,i),dp), aimag(selfen(3,3,i)), &
            real(selfen(2,2,i),dp), aimag(selfen(2,2,i))
    end do
    close(unit=11)
    ! }}}
  end subroutine write_selfen

  subroutine write_selfen_full_mtrx(z,selfen,filename)
    ! {{{
    complex(dpc), dimension(:), intent(in) :: z
    complex(dpc), dimension(:,:,:), intent(in) :: selfen
    character(len=*), intent(in) :: filename
    integer :: i
    open(unit=11,file=filename,action="write")
    do i=1, size(z)
       write(unit=11,fmt="(2f15.7)") real(z(i),dp), aimag(z(i))
       write(unit=11,fmt="(10f15.7)") real(selfen(:,:,i),dp)
       write(unit=11,fmt="(10f15.7)") aimag(selfen(:,:,i))
    end do
    close(unit=11)
    ! }}}
  end subroutine write_selfen_full_mtrx

end program EDdemo

! Local variables:
! folded-file: t
! End:
