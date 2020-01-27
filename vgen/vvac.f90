module m_vvac
  ! ****************************************************************
  ! calculates the g(2-dim)=0 part of the vacuum coulomb potential *
  ! for general symmetry.          c.l.fu, r.podloucky             *
  ! ****************************************************************
  contains

  subroutine vvac( vacuum, stars, cell, sym, input, field, psq, rht, vz, rhobar, sig1dh, vz1dh )

    use m_constants
    use m_qsf
    use m_types
    implicit none

    type(t_vacuum), intent(in)    :: vacuum
    type(t_stars),  intent(in)    :: stars
    type(t_cell),   intent(in)    :: cell
    type(t_sym),    intent(in)    :: sym

    complex,        intent(out)   :: rhobar
    real,           intent(out)   :: sig1dh, vz1dh
    type(t_input),  intent(in)    :: input
    type(t_field),  intent(inout) :: field ! efield is modified here

    real,           intent(in)    :: rht(vacuum%nmzd,2) 
    complex,        intent(in)    :: psq(stars%ng3)
    real,           intent(out)   :: vz(vacuum%nmzd,2)

    complex                       :: sumq, vcons
    real                          :: bj0, bj1, dh, qzh, sigmaa(2)
    integer                       :: ig3, imz, ivac, ncsh
    real                          :: f(vacuum%nmzd), sig(vacuum%nmzd), vtemp(vacuum%nmzd)
    intrinsic cos, sin


    vz(:,1:vacuum%nvac) = 0.0 ! initialize potential

    ! obtain mesh point (ncsh) of charge sheet for external electric field
    ncsh = field%efield%zsigma / vacuum%delz + 1.01
    sigmaa(1) = ( field%efield%sigma + field%efield%sig_b(1) ) / cell%area
    sigmaa(2) = ( field%efield%sigma + field%efield%sig_b(2) ) / cell%area

    ! g=0 vacuum potential due to neutral charge density
    ! inside slab and zero charge density outside

    vcons = - 2. * fpi_const * ImagUnit
    dh = cell%z1
    rhobar = - psq(1)
    sumq = cmplx( 0.0, 0.0 )

    do ig3 = 2, stars%ng3
      if (stars%ig2(ig3) == 1) then           ! select g_|| = 0
        qzh = stars%kv3(3,ig3) * cell%bmat(3,3) * dh
        bj0 = sin(qzh) / qzh
        rhobar = rhobar - psq(ig3) * bj0 * stars%nstr(ig3)
        if ( .not.( sym%zrfs .or. sym%invs ) ) then
          bj1 = ( sin(qzh) - qzh * cos(qzh) ) / ( qzh * qzh )
          sumq = sumq + bj1 * psq(ig3) * dh * dh
        endif
      endif
    enddo

    ivac = 2                        ! lower (ivac=2) vacuum
    if ( vacuum%nvac == 2 ) then
      vz(1:vacuum%nmz,ivac) = vcons * sumq
    endif

    ! g=0 vacuum potential due to
    ! negative of rhobar + vacuum (g=0) charge ----> v2(z)

    ivac = 1 ! upper vacuum

    ! =================== dirichlet ==============================
    if ( field%efield%dirichlet ) then
      vz(ncsh+1:vacuum%nmz,ivac) = field%efield%sig_b(1)
      call qsf( vacuum%delz, rht(1,ivac), sig, ncsh, 1 )
      sig(1:ncsh) = sig(ncsh) - sig(1:ncsh)
      call qsf( vacuum%delz, sig, vtemp, ncsh, 1 )
      do imz = 1, ncsh
        vz(imz,ivac) = - fpi_const * ( vtemp(ncsh) - vtemp(imz) ) + field%efield%sig_b(1)
      enddo
      sig1dh = sig(1)
      vz1dh = vz(1,ivac)   ! potential on vacuum boundary
      if ( vacuum%nvac == 1 ) return

      ivac = 2     ! lower vacuum
      call qsf( vacuum%delz, rht(1,ivac), sig, ncsh, 1 )
      f(1:ncsh) = sig(1:ncsh) - rhobar*vacuum%dvac + sig1dh
      call qsf( vacuum%delz, f, vtemp, ncsh, 1 )
      do imz = 1,ncsh
        vz(imz,ivac) = - fpi_const * ( vtemp(imz) + sig1dh * vacuum%dvac - rhobar * vacuum%dvac * vacuum%dvac / 2. ) + vz(imz,ivac) + vz1dh
      end do

      ! force matching on the other side
      field%efield%vslope = ( field%efield%sig_b(2) - vz(ncsh,1) ) / ( 2 * vacuum%delz * ( ncsh + 1 ) + vacuum%dvac )
      ivac = 1
      do imz = 1, ncsh
        vz(imz,ivac) = vz(imz,ivac) + field%efield%vslope * vacuum%delz * ( ncsh - imz + 1 )
      end do
      ivac = 2
      do imz = 1, ncsh
        vz(imz,ivac) = vz(imz,ivac) + field%efield%vslope * ( vacuum%dvac + vacuum%delz * imz + vacuum%delz * ncsh )
      end do
      vz(ncsh+1:vacuum%nmz,ivac) = field%efield%sig_b(2)

      ! =================== neumann ==============================
    else ! neumann

      call qsf( vacuum%delz, rht(1,ivac), sig, vacuum%nmz, 1 )
      sig1dh = sig(vacuum%nmz) - sigmaa(1)  ! need to include contribution from electric field
      sig(1:vacuum%nmz) = sig(vacuum%nmz) - sig(1:vacuum%nmz)
      call qsf( vacuum%delz, sig, vtemp, vacuum%nmz, 1 )
      ! external electric field contribution 
      do imz = 1, ncsh
        vz(imz,ivac) = - fpi_const * ( vtemp(vacuum%nmz) - vtemp(imz) ) + vz(imz,ivac) - fpi_const * ( imz - ncsh ) * vacuum%delz * sigmaa(1)
      enddo
      do imz = ncsh + 1, vacuum%nmz
        vz(imz,ivac) = - fpi_const * ( vtemp(vacuum%nmz) - vtemp(imz) ) + vz(imz,ivac)
      enddo
      vz1dh = vz(1,ivac)   ! potential on vacuum boundary
      if ( vacuum%nvac == 1 ) return

      ivac = 2 ! lower vacuum
      call qsf( vacuum%delz, rht(1,ivac), sig, vacuum%nmz, 1 )
      f(1:vacuum%nmz) = sig(1:vacuum%nmz) - rhobar * vacuum%dvac + sig1dh
      call qsf( vacuum%delz, f, vtemp, vacuum%nmz, 1 )

      ! external electric field contribution
      do imz = 1, ncsh
        vz(imz,ivac) = - fpi_const * ( vtemp(imz) + sig1dh * vacuum%dvac - rhobar * vacuum%dvac * vacuum%dvac / 2. ) + vz1dh + vz(imz,ivac)
      enddo
      do imz = ncsh + 1, vacuum%nmz
        vz(imz,ivac) = - fpi_const * ( vtemp(imz) + sig1dh * vacuum%dvac - rhobar * vacuum%dvac * vacuum%dvac / 2. ) + vz1dh + vz(imz,ivac) + fpi_const * ( imz - ncsh ) * vacuum%delz * sigmaa(2)
      enddo

    end if ! dirichlet (vs. neumann)

  end subroutine vvac

end module m_vvac
