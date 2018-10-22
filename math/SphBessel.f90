!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_SphBessel

  !-------------------------------------------------------------------------
  ! SphBessel calculates spherical Bessel functions of the first, 
  ! second and third kind (Bessel, Neumann and Hankel functions).
  ! ModSphBessel calculates modified spherical Bessel functions
  ! of the first and second kind.
  !
  ! jl  :          spherical Bessel function of the first  kind (Bessel)
  ! nl  :          spherical Bessel function of the second kind (Neumann)
  ! hl  :          spherical Bessel function of the third  kind (Hankel)
  ! il  : modified spherical Bessel function of the first  kind
  ! kl  : modified spherical Bessel function of the second kind
  !
  ! z   : Bessel functions are calculated for this value
  ! lmax: Bessel functions are calculated for all the indices l
  !       from 0 to lmax
  !
  ! intent(in):
  ! z   : complex or real scalar
  ! lmax: integer
  ! 
  ! intent(out):
  ! * SphBessel( jl, nl, hl, z, lmax )
  ! jl:   complex or real, dimension(0:lmax)
  ! nl:   complex or real, dimension(0:lmax)
  ! hl:   complex,         dimension(0:lmax)
  ! * ModSphBessel( il, kl, z, lmax )
  ! il:   complex or real, dimension(0:lmax)
  ! kl:   complex or real, dimension(0:lmax)
  ! 
  ! All subroutines are pure and therefore can be called for a range of 
  ! z-values concurrently, f.e. this way:
  ! allocate( il(0:lmax, size(z)), kl(0:lmax, size(z)) )
  ! do concurrent (i = 1: size(z))
  !   call ModSphBessel( il(:,i), kl(:,i), z(i), lmax )
  ! end do
  !
  ! details on implementation:
  ! For |z| <= 1 the taylor expansions of jl and nl are used.
  ! For |z| >  1 the explicit expressions for hl(+), hl(-) are used.
  ! For modified spherical Bessel functions il and kl the relations
  !                   il(z) =  I^{-l} * jl(I*z)
  !                   kl(z) = -I^{l}  * hl(I*z)
  ! are used.
  !
  ! authors:
  ! originally written by             R. Zeller (1990)
  ! modernised and extended by        M. Hinzen (2016)
  !-------------------------------------------------------------------------

  use m_constants, only: ImagUnit
  implicit none


  interface SphBessel
    module procedure  SphBesselComplex, SphBesselReal
  end interface

  interface ModSphBessel
    ! variant Complex2 takes workspace as an argument.
    ! this is not possible for the subroutine working on reals.
    module procedure  ModSphBesselComplex, ModSphBesselReal, ModSphBesselComplex2
  end interface


contains


  pure subroutine SphBesselComplex ( jl, nl, hl, z, lmax )

    complex,                    intent(in)  :: z
    integer,                    intent(in)  :: lmax
    complex, dimension(0:lmax), intent(out) :: jl, nl, hl

    complex                                 :: termj, termn, z2, zj, zn
    real                                    :: rl, rn
    real,    dimension(0:lmax)              :: rnm
    integer                                 :: l, m, n


    zj = 1.0
    zn = 1.0 / z
    z2 = z * z
    jl(:) = 1.0
    nl(:) = 1.0
    if ( abs( z ) < lmax + 1.0 ) then
      SERIAL_L_LOOP: do l = 0, lmax
        rl = l + l
        termj = 1.0
        termn = 1.0
        EXPANSION: do n = 1, 25
          rn = n + n
          termj = -termj / ( rl + rn + 1.0 ) / rn * z2
          termn =  termn / ( rl - rn + 1.0 ) / rn * z2
          jl(l) = jl(l) + termj
          nl(l) = nl(l) + termn
        end do EXPANSION
        jl(l) =  jl(l) * zj
        nl(l) = -nl(l) * zn
        hl(l) =  jl(l) + nl(l) * ImagUnit
        zj = zj * z / ( rl + 3.0 )
        zn = zn / z * ( rl + 1.0 )
      end do SERIAL_L_LOOP
    end if

    rnm(:) = 1.0
    !PARALLEL_L_LOOP: do concurrent (l = 0: lmax)
    PARALLEL_L_LOOP: do l = 0,lmax
      if ( abs( z ) >= l + 1.0 ) then
        hl(l) = 0.0
        nl(l) = 0.0
        SERIAL_M_LOOP: do m = 0, l
          hl(l) = hl(l) + (-1) ** m * rnm(l)
          nl(l) = nl(l) + rnm(l)
          rnm(l) = rnm(l) / ( m + 1.0 ) * ( l * ( l + 1 ) - m * ( m + 1 ) ) / ( ImagUnit * ( z + z ) )
        end do SERIAL_M_LOOP
        hl(l) = hl(l) * (-ImagUnit) ** l * exp(  ImagUnit * z ) / (  ImagUnit * z )
        nl(l) = nl(l) *   ImagUnit  ** l * exp( -ImagUnit * z ) / ( -ImagUnit * z )
        jl(l) = ( hl(l) + nl(l) ) / 2.0
        nl(l) = ( hl(l) - jl(l) ) * (-ImagUnit)
      end if
    end do PARALLEL_L_LOOP

  end subroutine SphBesselComplex



  pure subroutine SphBesselReal ( jl, nl, hl, x, lmax )

    real,                       intent(in)  :: x
    integer,                    intent(in)  :: lmax
    real,    dimension(0:lmax), intent(out) :: jl, nl
    complex, dimension(0:lmax), intent(out) :: hl

    complex, dimension(0:lmax)              :: jl_complex, nl_complex
    complex                                 :: z


    z = x           ! internal conversion from real to complex

    call SphBesselComplex( jl_complex, nl_complex, hl, z, lmax )

    jl = jl_complex ! internal conversion from complex to real
    nl = nl_complex ! internal conversion from complex to real

  end subroutine SphBesselReal



  pure subroutine ModSphBesselComplex ( il, kl, z, lmax )

    complex,                    intent(in)  :: z
    integer,                    intent(in)  :: lmax
    complex, dimension(0:lmax), intent(out) :: il, kl

    complex, dimension(0:lmax)              :: nl
    integer                                 :: l 


    call SphBesselComplex( il, nl, kl, ImagUnit * z, lmax )

    do l = 0, lmax
      il(l) = (-ImagUnit) ** l * il(l)
      kl(l) = - ImagUnit  ** l * kl(l)
    end do     

  end subroutine ModSphBesselComplex


  !another implementation of ModSphBesselComplex, where nl is allocated outside for performance reasons
  pure subroutine ModSphBesselComplex2 ( il, kl, nl, z, lmax )

    complex,                    intent(in)  :: z
    integer,                    intent(in)  :: lmax
    complex, dimension(0:lmax), intent(out) :: il, kl, nl

    integer                                 :: l 


    call SphBesselComplex( il, nl, kl, ImagUnit * z, lmax )

    do l = 0, lmax
      il(l) = (-ImagUnit) ** l * il(l)
      kl(l) = - ImagUnit  ** l * kl(l)
    end do     

  end subroutine ModSphBesselComplex2



  pure subroutine ModSphBesselReal ( il, kl, x, lmax )

    real,                       intent(in)  :: x
    integer,                    intent(in)  :: lmax
    real,    dimension(0:lmax), intent(out) :: il, kl

    complex, dimension(0:lmax)              :: jl, nl, hl
    integer                                 :: l 
    complex                                 :: z


    z = ImagUnit * x
    call SphBesselComplex( jl, nl, hl, z, lmax )

    do l = 0, lmax
      il(l) = (-ImagUnit) ** l * jl(l)
      kl(l) = - ImagUnit  ** l * hl(l)
    end do     

  end subroutine ModSphBesselReal


end module m_SphBessel
