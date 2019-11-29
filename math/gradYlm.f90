MODULE m_gradYlm

CONTAINS

  ! Not from me, taken from Spex, ask Christoph! Calculates derivative of a function with scalar argument lying on a muffin-tin mesh
  subroutine Derivative(f, itype, atoms, df)

    use m_types

    implicit none

     type(t_atoms), intent(in)  :: atoms

     integer,       intent(in)  :: itype
     real,          intent(in)  :: f(atoms%jri(itype))
     real,          intent(out) :: df(atoms%jri(itype))
     real                       :: h, r, d21, d32, d43, d31, d42, d41, df1, df2, s
     real                       :: y0, y1, y2
     integer                    :: i, n

     n = atoms%jri(itype)
     h = atoms%dx(itype)
     r = atoms%rmsh(1, itype)

     ! use Lagrange interpolation of 3rd order (and averaging) for points 3 to n
     d21 = r * (exp(h)-1) ; d32 = d21 * exp(h) ; d43 = d32 * exp(h)
     d31 = d21 + d32      ; d42 = d32 + d43
     d41 = d31 + d43
     df(1) =   d31*d41 / (d21*d32*d42) * f(2) + ( -1.0/d21 - 1.0/d31 - 1.0/d41) * f(1)&
    &        - d21*d41 / (d31*d32*d43) * f(3) + d21*d31 / (d41*d42*d43) * f(4)
     df(2) = - d32*d42 / (d21*d31*d41) * f(1) + (  1.0/d21 - 1.0/d32 - 1.0/d42) * f(2)&
    &        + d21*d42 / (d31*d32*d43) * f(3) - d21*d32 / (d41*d42*d43) * f(4)
     df1   =   d32*d43 / (d21*d31*d41) * f(1) - d31*d43 / (d21*d32*d42) * f(2) +&
    &  ( 1.0/d31 + 1.0/d32 - 1.0/d43 ) * f(3) + d31*d32 / (d41*d42*d43) * f(4)
     do i = 3, n - 2
        d21 = d32 ; d32 = d43 ; d43 = d43 * exp(h)
        d31 = d42 ; d42 = d42 * exp(h)
        d41 = d41 * exp(h)
        df2   = - d32*d42 / (d21*d31*d41) * f(i-1) + ( 1.0/d21 - 1.0/d32 - 1.0/d42) * f(i) + &
    &             d21*d42 / (d31*d32*d43) * f(i+1) - d21*d32 / (d41*d42*d43) * f(i+2)
        df(i) = ( df1 + df2 ) / 2
        df1   = d32*d43 / (d21*d31*d41) * f(i-1) - d31*d43 / (d21*d32*d42) * f(i) +&
    &    ( 1.0/d31 + 1.0/d32 - 1.0/d43 ) * f(i+1) + d31*d32 / (d41*d42*d43) * f(i+2)
     enddo
     df(n-1) = df1
     df(n)   = - d42*d43 / (d21*d31*d41) * f(n-3) + d41*d43 / (d21*d32*d42) * f(n-2) -&
    &            d41*d42 / (d31*d32*d43) * f(n-1) + ( 1.0/d41 + 1.0/d42 + 1.0/d43 ) * f(n)
     ! for first two points use Lagrange interpolation of second order for log(f(i))
     ! or, as a fall-back, Lagrange interpolation with the conditions f(1), f(2), f(3), f'(3).
     s = sign(1.0,f(1))
     if(sign(1.0,f(2)) /= s .or. sign(1.0,f(3))  /= s .or. any(f(:3) == 0)) then
        d21   = r * (exp(h)-1)
        d32   = d21 * exp(h)
        d31   = d21 + d32
        s     = df(3) / (d31*d32) - f(1) / (d21*d31**2) + f(2) / (d21*d32**2) - f(3) / (d31**2*d32) - f(3) / (d31*d32**2)
        df(1) = - (d21+d31) / (d21*d31) * f(1) + d31 / (d21*d32) * f(2) - d21 / (d31*d32) * f(3) + d21*d31 * s

        df(2) = - (d21-d32) / (d21*d32) * f(2) - d32 / (d21*d31) * f(1) + d21 / (d31*d32) * f(3) - d21*d32 * s
     else
        y0    = log(abs(f(1)))
        y1    = log(abs(f(2)))
        y2    = log(abs(f(3)))
        df(1) = ( - 3*y0/2 + 2*y1 - y2/2 ) * f(1) / (h*r)
        df(2) = (y2-y0)/2                  * f(2) / (h*r*exp(h))
     endif
  end subroutine Derivative

   SUBROUTINE gradYlm(atoms, r2FshMt, r2GrFshMt)
   ! Based on work for juphon by C. Gerhorst.
  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the spherical harmonic expansion coefficients of the muffin-tin gradient applied to an arbitrary function multiplied
  !> by $r^2$. The resulting gradient expansion coefficients are multiplied by a factor $r^2$.
  !>
  !> @note
  !> The ingoing function is assumed to be multiplied with $r^2$. The outgoing resulting function is also multiplied by $r^2$.
  !>
  !> @param[in]  atoms     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in]  r2FshMt   : Lattice harmonic coefficients of muffin-tin quantity multiplied by a factor of r**2.
  !> @param[out] r2GrFshMt : Spherical harmonic coefficients of muffin-tin quantity's gradient multiplied by a factor of r**2
  !---------------------------------------------------------------------------------------------------------------------------------

    use m_constants, only : fpi_const, ImagUnit
    use m_gaunt, only : gaunt1
    use m_types
    USE m_grdchlh

    implicit none

    ! Type parameter
    ! ***************
    type(t_atoms),                     intent(in)  :: atoms

    ! Array parameters
    ! ****************
    complex,       allocatable,        intent(in)  :: r2FshMt(:, :, :)
    complex,       allocatable,        intent(out) :: r2GrFshMt(:, :, :, :)


    ! Local Scalar Variables
    ! **********************
    ! pfac     : Prefactor
    ! tGaunt   : Gaunt coefficient
    ! itype    : Loop index for atom types
    ! ieqat    : Loop index for equivalent atoms
    ! iatom    : Loop index for all atoms
    ! imesh    : Loop index for radial mesh point
    ! mqn_m    : Magnetic quantum number m
    ! oqn_l    : Orbital quantum number l
    ! mqn_mpp  : Magnetic quantum number double primed to index the natural coordinates
    ! lmInput  : Collective index for orbital and magnetic quantum number used for input function
    ! lmOutput : Collective index for orbital and magnetic quantum number used for output function
    real                                           :: pfac
    real                                           :: tGaunt
    integer                                        :: itype
    integer                                        :: ieqat
    integer                                        :: iatom
    integer                                        :: imesh
    integer                                        :: mqn_m
    integer                                        :: oqn_l
    integer                                        :: mqn_mpp
    integer                                        :: lmOutput
    integer                                        :: lmInput

    ! Local Array Variables
    ! *********************
    ! rDerFlhMtre  : Real part of the radial derrivative applied to the incoming fuction coefficients
    ! rDerFlhMtim  : Imaginary part of the radial derrivative applied to the incoming fuction coefficients
    ! rDerFlhMt    : Radial derrivative of the incoming fuction
    ! r2GrFshMtNat : Expansion coefficients of the muffin-tin gradient applied to the incoming function. The coefficients are given
    !                in natural coordinates and multiplied by $r^2$
    real,          allocatable                     :: rDerFshMtre(:)
    real,          allocatable                     :: rDerFshMtim(:)
    real,          allocatable                     :: rDerJunk(:)
    complex,       allocatable                     :: rDerFshMt(:)
    complex,       allocatable                     :: r2GrFshMtNat(:, :, :, :)
    !Matrix syntax idea from http://stackoverflow.com/questions/3708307/how-to-initialize-two-dimensional-arrays-in-fortran
    complex, parameter, dimension(3, 3)  :: Tmatrix = transpose(reshape([                                                          &
                                                                   & cmplx(1 / sqrt(2.), 0), cmplx(0, 0), cmplx(-1 / sqrt(2.), 0), &
                                                                   & -ImagUnit / sqrt(2.), cmplx(0, 0), -ImagUnit / sqrt(2.), &
                                                                   & cmplx(0, 0), cmplx(1, 0), cmplx(0, 0) &
                                                                   & ], [3, 3] )) !< see my notes


    ! Initialization of additionaly required arrays.
    allocate( r2GrFshMt(atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    allocate( r2GrFshMtNat(atoms%jmtd, ( atoms%lmaxd + 2 )**2, atoms%nat, 3) )
    allocate( rDerFshMtre(atoms%jmtd), rDerFshMtim(atoms%jmtd), rDerJunk(atoms%jmtd), rDerFshMt(atoms%jmtd) )

    rDerFshMtre(:) = 0.
    rDerFshMtim(:) = 0.
    rDerFshMt(:) = 0.
    r2GrFshMt = cmplx(0., 0.)
    r2GrFshMtNat = cmplx(0., 0.)

    pfac = sqrt( fpi_const / 3 )
    do mqn_mpp = -1, 1
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l

              ! l + 1 block
              if ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) then
                lmOutput = ( oqn_l + 1 ) * ( oqn_l + 2 ) + 1 + mqn_m - mqn_mpp
                lmInput = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
                rDerFshMtre(:) = 0.
                rDerFshMtim(:) = 0.
                ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe store it?
                !CALL grdchlh(1,1,atoms%jri(itype),atoms%dx(itype),atoms%rmsh(1, itype),REAL(r2FshMt(:, lmInput, itype)),6,rDerFshMtre,rDerJunk)
                !CALL grdchlh(1,1,atoms%jri(itype),atoms%dx(itype),atoms%rmsh(1, itype),AIMAG(r2FshMt(:, lmInput, itype)),6,rDerFshMtim,rDerJunk)
                call Derivative( real(r2FshMt(:, lmInput, itype)), itype, atoms, rDerFshMtre )
                call Derivative( aimag(r2FshMt(:, lmInput, itype)), itype, atoms, rDerFshMtim )
                do imesh = 1, atoms%jri(itype)
                  rDerFshMt(imesh) = cmplx(rDerFshMtre(imesh), rDerFshMtim(imesh))
                end do ! imesh
                tGaunt = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd + 1 )
                do imesh = 1, atoms%jri(itype)
                  r2GrFshMtNat(imesh, lmOutput, iatom, mqn_mpp + 2) = r2GrFshMtNat(imesh, lmOutput, iatom, mqn_mpp + 2) + pfac *   &
                    & (-1)**mqn_mpp * tGaunt &
                    & * (rDerFshMt(imesh) - ((oqn_l + 2)* r2FshMt(imesh, lmInput, iatom) / atoms%rmsh(imesh, itype)))
                end do ! imesh
              end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 )

              ! l - 1 block
              if ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) then
                lmInput = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
                lmOutput = (oqn_l - 1) * oqn_l + 1 + mqn_m - mqn_mpp
                rDerFshMtre(:) = 0.
                rDerFshMtim(:) = 0.
                ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe store it?
                !CALL grdchlh(1,1,atoms%jri(itype),atoms%dx(itype),atoms%rmsh(1, itype),REAL(r2FshMt(:, lmInput, itype)),6,rDerFshMtre,rDerJunk)
                !CALL grdchlh(1,1,atoms%jri(itype),atoms%dx(itype),atoms%rmsh(1, itype),AIMAG(r2FshMt(:, lmInput, itype)),6,rDerFshMtim,rDerJunk)
                call Derivative( real(r2FshMt(:, lmInput, itype)), itype, atoms, rDerFshMtre )
                call Derivative( aimag(r2FshMt(:, lmInput, itype)), itype, atoms, rDerFshMtim )
                do imesh = 1, atoms%jri(itype)
                  rDerFshMt(imesh) = cmplx(rDerFshMtre(imesh), rDerFshMtim(imesh))
                end do ! imesh
                tGaunt = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd + 1 )
                do imesh = 1, atoms%jri(itype)
                  r2GrFshMtNat(imesh, lmOutput, iatom, mqn_mpp + 2) = r2GrFshMtNat(imesh, lmOutput, iatom, mqn_mpp + 2) + pfac *   &
                    & (-1)**mqn_mpp * tGaunt * ( rDerFshMt(imesh) + ( (oqn_l - 1) * r2FshMt(imesh, lmInput, iatom)&
                    & / atoms%rmsh(imesh, itype) ) )
                enddo ! imesh
              end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 )
            end do ! mqn_m
          end do ! oqn_l
        end do ! ieqat
      end do ! itype
    end do ! mqn_mpp

    ! Conversion from natural to cartesian coordinates
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype) + 1
          do mqn_m = -oqn_l, oqn_l
            lmOutput = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              r2GrFshMt(imesh, lmOutput, iatom, 1:3) = matmul( Tmatrix(1:3, 1:3), r2GrFshMtNat(imesh, lmOutput, iatom, 1:3) )
            end do ! imesh
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype

   END SUBROUTINE gradYlm

   SUBROUTINE divYlm(gradMtx, gradMty, gradMtz, divMt)
      COMPLEX, INTENT(IN) :: gradMtx(:,:,:,:), gradMty(:,:,:,:), gradMtz(:,:,:,:) !r,lm,n,x
      COMPLEX, INTENT(OUT) :: divMt(:,:,:)

      divMT(:,:,:)=gradMtx(:,:,:,1)+gradMty(:,:,:,2)+gradMtz(:,:,:,3)

   END SUBROUTINE divYlm
END MODULE m_gradYlm
