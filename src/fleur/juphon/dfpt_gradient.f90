MODULE m_dfpt_gradient
   USE m_juDFT
   USE m_constants
   USE m_types

   IMPLICIT NONE
CONTAINS
    subroutine mt_gradient_new(atoms, sphhar, sym, r2FlhMt, GrFshMt)

      use m_gaunt, only : gaunt1

      type(t_atoms),               intent(in)  :: atoms
      type(t_sphhar),              intent(in)  :: sphhar
      type(t_sym),                 intent(in)  :: sym

      real,                        intent(in)  :: r2FlhMt(:, 0:, :)
      complex,                     intent(out) :: GrFshMt(:, :, :, :)

      real                                     :: pfac
      real                                     :: tGaunt
      integer                                  :: itype
      integer                                  :: imesh
      integer                                  :: mqn_m
      integer                                  :: oqn_l
      integer                                  :: mqn_mpp
      integer                                  :: lm
      integer                                  :: symType
      integer                                  :: ilh
      integer                                  :: imem

      real,           allocatable              :: rDerFlhMt(:)
      complex,        allocatable              :: r2GrFshMtNat(:, :, :, :)

      allocate( r2GrFshMtNat(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3) )
      allocate( rDerFlhMt(atoms%jmtd) )
      GrFshMt = cmplx(0., 0.)
      r2GrFshMtNat = cmplx(0., 0.)
      rDerFlhMt = 0.

      pfac = sqrt( fpi_const / 3. )
      do mqn_mpp = -1, 1
        do itype = 1, atoms%ntype
            symType = sym%ntypsy(itype)
            do ilh = 0, sphhar%nlh(symType)
              oqn_l = sphhar%llh(ilh, symType)
              do imem = 1, sphhar%nmem(ilh,symType)
                mqn_m = sphhar%mlh(imem,ilh,symType)

                ! l + 1 block
                ! oqn_l - 1 to l, so oqn_l should be < lmax not <= lmax
                if ( ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < atoms%lmax(itype)) ) then
                  lm = ( oqn_l + 1 ) * ( oqn_l + 2 ) + 1 + mqn_m - mqn_mpp
                  call derivative_loc( r2FlhMt(:, ilh, itype), itype, atoms, rDerFlhMt )
                  tGaunt = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                  do imesh = 1, atoms%jri(itype)
                    r2GrFshMtNat(imesh, lm, itype, mqn_mpp + 2) = r2GrFshMtNat(imesh, lm, itype, mqn_mpp + 2) + pfac * (-1)**mqn_mpp &
                      &* tGaunt * (rDerFlhMt(imesh) * sphhar%clnu(imem,ilh,symType) &
                      &- ((oqn_l + 2) * r2FlhMt(imesh, ilh, itype) * sphhar%clnu(imem,ilh,symType) / atoms%rmsh(imesh, itype)))
                  end do ! imesh
                end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l )

                ! l - 1 block
                if ( ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
                  if ( oqn_l - 1 == -1 ) then
                    write (*, *) 'oqn_l too low'
                  end if
                  lm = (oqn_l - 1) * oqn_l + 1 + mqn_m - mqn_mpp
                  ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe store it?
                  call derivative_loc( r2FlhMt(:, ilh, itype), itype, atoms, rDerFlhMt )
                  tGaunt = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                  do imesh = 1, atoms%jri(itype)
                    r2GrFshMtNat(imesh, lm, itype, mqn_mpp + 2) = r2GrFshMtNat(imesh, lm, itype, mqn_mpp + 2) + pfac * (-1)**mqn_mpp &
                      & * tGaunt * (rDerFlhMt(imesh)  * sphhar%clnu(imem,ilh,symType) &
                      & + ((oqn_l - 1) * r2FlhMt(imesh, ilh, itype) * sphhar%clnu(imem,ilh,symType) / atoms%rmsh(imesh, itype)))
                  end do ! imesh
                end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l )
              end do ! imem
            end do ! ilh
        end do ! itype
      end do ! mqn_mpp

      ! Conversion from natural to cartesian coordinates
      do itype = 1, atoms%ntype
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
              do imesh = 1, atoms%jri(itype)
                grFshMt(imesh, lm, itype, 1:3) = matmul( Tmatrix0(1:3, 1:3), r2GrFshMtNat(imesh, lm, itype, 1:3) ) / atoms%rmsh(imesh, itype)**2
              end do
            end do ! mqn_m
          end do ! oqn_l
      end do ! itype

   end subroutine mt_gradient_new

   subroutine derivative_loc(f, itype, atoms, df)

      integer,       intent(in)  :: itype
      type(t_atoms), intent(in)  :: atoms
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
      df(1) =   d31*d41 / (d21*d32*d42) * f(2) + ( -1d0/d21 - 1d0/d31 - 1d0/d41) * f(1)&
     &        - d21*d41 / (d31*d32*d43) * f(3) + d21*d31 / (d41*d42*d43) * f(4)
      df(2) = - d32*d42 / (d21*d31*d41) * f(1) + (  1d0/d21 - 1d0/d32 - 1d0/d42) * f(2)&
     &        + d21*d42 / (d31*d32*d43) * f(3) - d21*d32 / (d41*d42*d43) * f(4)
      df1   =   d32*d43 / (d21*d31*d41) * f(1) - d31*d43 / (d21*d32*d42) * f(2) +&
     &  ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(3) + d31*d32 / (d41*d42*d43) * f(4)
      do i = 3, n - 2
         d21 = d32 ; d32 = d43 ; d43 = d43 * exp(h)
         d31 = d42 ; d42 = d42 * exp(h)
         d41 = d41 * exp(h)
         df2   = - d32*d42 / (d21*d31*d41) * f(i-1) + ( 1d0/d21 - 1d0/d32 - 1d0/d42) * f(i) + &
     &             d21*d42 / (d31*d32*d43) * f(i+1) - d21*d32 / (d41*d42*d43) * f(i+2)
         df(i) = ( df1 + df2 ) / 2
         df1   = d32*d43 / (d21*d31*d41) * f(i-1) - d31*d43 / (d21*d32*d42) * f(i) +&
     &    ( 1d0/d31 + 1d0/d32 - 1d0/d43 ) * f(i+1) + d31*d32 / (d41*d42*d43) * f(i+2)
      enddo
      df(n-1) = df1
      df(n)   = - d42*d43 / (d21*d31*d41) * f(n-3) + d41*d43 / (d21*d32*d42) * f(n-2) -&
     &            d41*d42 / (d31*d32*d43) * f(n-1) + ( 1d0/d41 + 1d0/d42 + 1d0/d43 ) * f(n)
      ! for first two points use Lagrange interpolation of second order for log(f(i))
      ! or, as a fall-back, Lagrange interpolation with the conditions f(1), f(2), f(3), f'(3).
      s = sign(1d0,f(1))
      if(sign(1d0,f(2)) /= s .or. sign(1d0,f(3))  /= s .or. any(abs(f(:3)) < 1e0)) then
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
   end subroutine derivative_loc

    SUBROUTINE sh_to_lh(sym, atoms, sphhar, jspins, radfact, rhosh, rholhreal, rholhimag)

        ! WARNING: This routine will not fold back correctly for activated sym-
        !          metry and gradients (rho in l=0 and lattice harmonics do not
        !          allow l=1 --> gradient in l=1 is lost)

        TYPE(t_sym),    INTENT(IN)  :: sym
        TYPE(t_atoms),  INTENT(IN)  :: atoms
        TYPE(t_sphhar), INTENT(IN)  :: sphhar
        INTEGER,        INTENT(IN)  :: jspins, radfact
        COMPLEX,        INTENT(IN)  :: rhosh(:, :, :, :)
        REAL,           INTENT(OUT) :: rholhreal(:, 0:, :, :), rholhimag(:, 0:, :, :)

        INTEGER :: iSpin, iType, iAtom, ilh, iMem, ilm, iR
        INTEGER :: ptsym, l, m
        REAL    :: factor

        rholhreal = 0.0
        rholhimag = 0.0
        print*,"Hallo"
        DO iSpin = 1, jspins
            DO iType = 1, atoms%ntype
                iAtom = atoms%firstAtom(iType)
                ptsym = sym%ntypsy(iAtom)
                print*, "ptsym", ptsym
                DO ilh = 0, sphhar%nlh(ptsym)
                    l = sphhar%llh(iLH, ptsym)
                    print*,"l", l
                    print*,"sphhar%nmem(ilh, ptsym)",sphhar%nmem(ilh, ptsym)
                    DO iMem = 1, sphhar%nmem(ilh, ptsym)
                        print*,"iMem-loop"
                        m = sphhar%mlh(iMem, ilh, ptsym)
                      
                        print*, 'm',m
                        ilm = l * (l+1) + m + 1
                        print*, 'ilm', ilm
                        print*,"atoms%jri(iType)",atoms%jri(iType)
                        print*, "sphhar%clnu(1, ilh, ptsym)",sphhar%clnu(1, ilh, ptsym)
                        DO iR = 1, atoms%jri(iType)
                          !print*, "iR", iR
                          !print*,"rholhreal(iR, ilh, iType, iSpin)",rholhreal(iR, ilh, iType, iSpin)
                           IF ((radfact.EQ.0).AND.(l.EQ.0)) THEN
                               factor = atoms%rmsh(iR, iType) / sfp_const
                           ELSE IF (radfact.EQ.2) THEN
                               factor = atoms%rmsh(iR, iType)**2
                           ELSE
                               factor = 1.0
                           END IF
                            rholhreal(iR, ilh, iType, iSpin) = &
                          & rholhreal(iR, ilh, iType, iSpin) + &
                          &  real(rhosh(iR, ilm, iatom, iSpin) * conjg(sphhar%clnu(iMem, ilh, ptsym))) * factor
                            rholhimag(iR, ilh, iType, iSpin) = &
                          & rholhimag(iR, ilh, iType, iSpin) + &
                          & aimag(rhosh(iR, ilm, iatom, iSpin) * conjg(sphhar%clnu(iMem, ilh, ptsym))) * factor
                        END DO
                    END DO
                END DO
            END DO
        END DO

    END SUBROUTINE sh_to_lh
END MODULE m_dfpt_gradient