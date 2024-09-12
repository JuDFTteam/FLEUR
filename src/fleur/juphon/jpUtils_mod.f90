!----------------------------------------------------------------------------------------------------------------------------------------
! Forschungszentrum Jülich, juPhon Plugin for the FLEUR program
!----------------------------------------------------------------------------------------------------------------------------------------
!
! MODULE: Calculates Required Input Quantities for juPhon
!
!> @author
!> Christian-Roman Gerhorst
!
!> @brief
!> This module represents a toolbox of util routines which are needed within the program.
!>
!> @details checkisinteger
!>
!> @note
!> Additional information and formulas pointing out the routines of this module can be found within this
!> <a href='jpUtils.pdf'>document</a>.
!----------------------------------------------------------------------------------------------------------------------------------------
module mod_juPhonUtils

  implicit none

  contains

  function checkisinteger(string)
      implicit none

      logical                  :: checkisinteger
      character(*), intent(in) :: string
      character                :: c
      integer                  :: i
      logical                  :: number

      checkisinteger = .true.
      number         = .false.
      do i = 1, len_trim(string)
        c = string(i:i)
        if(all(c /= (/'0','1','2','3','4','5','6','7','8','9',' ','-'/))) checkisinteger = .false.
        if(any(c == (/'0','1','2','3','4','5','6','7','8','9'/)))         number         = .true.
      enddo
      if(.not.number) checkisinteger = .false.
    end function checkisinteger

    subroutine fopen(iunit, name, status, action, form, access, position, recl)
      implicit none

      !TODO comment on meaning of variables
      integer(4),                 intent(in) :: iunit
      character(len=*), optional             :: name
      character(256)                         :: name_temp
      character(len=*), optional             :: status, action, form, access, position
      character(11)                          :: status_temp, action_temp, form_temp, access_temp, position_temp
      integer,          optional             :: recl
      integer                                :: recl_temp
      logical                                :: ldum, l_opened
      integer(4)                             :: i

      ! set default values for local variables
      status_temp   = 'old'
      action_temp   = 'readwrite'
      form_temp     = 'formatted'
      access_temp   = 'sequential'
      position_temp = 'asis'
      recl_temp     = 0

      if(present(status))    status_temp   = status
      if(present(action))    action_temp   = action
      if(present(form))      form_temp     = form
      if(present(access))    access_temp   = access
      if(present(position))  position_temp = position

      if(present(recl)) then
        recl_temp = recl
      else if(access_temp ==  'direct') then
        NOstopNO'file: No record length given for direct-access file.' ! TODO replace NOstopNOby juDFT error? valid for whole file
      end if

      if(present(name)) then
        if(status_temp == 'scratch') then
          NOstopNO'file: Filename given for scratch file.'
        end if
        name_temp = name
      else
        write(name_temp, '(A,I3.3)') 'fort.', iunit
      end if

      inquire(iunit, opened=l_opened)
      if(l_opened) then
        write(6,'(/2A)')  'fopen: Could not open ', trim(name_temp)
        write(6,'(A)')    '       Unit already connected'
        NOstopNO             'fopen: Error while trying to open file'
      end if

      ! Some compilers don't like recl if access is sequential.
      if(access_temp ==  'direct') then
        if(present(name)) then
          open(iunit, file=name_temp, status=status_temp, action=action_temp, form=form_temp, access=access_temp, recl=recl_temp, iostat=i)
        else
          open(iunit,                 status=status_temp, action=action_temp, form=form_temp, access=access_temp, recl=recl_temp, iostat=i)
        end if
      else
        if(present(name)) then
          open(iunit, file=name_temp, status=status_temp, action=action_temp, form=form_temp, access=access_temp, iostat=i, position=position_temp)
        else
          open(iunit,                 status=status_temp, action=action_temp, form=form_temp, access=access_temp, iostat=i, position=position_temp)
        end if
      end if

      if(i /= 0) then
        write(6,'(/2A)')  'fopen: Could not open ', trim(name_temp)
        write(6,'(A,I4)') '       Error code is', i
        NOstopNO             'fopen: Error while trying to open a file.'
      endif
    end subroutine fopen


    ! TODO comment on this routine
    subroutine fclose(iunit, status)

      implicit none

      ! TODO comment on variables
      integer(4),       intent(in)  :: iunit
      character(len=*), optional    :: status
      character(11)                 :: status_temp
      integer(4)                    :: i

      ! set default values
      status_temp = 'keep'

      if(present(status)) status_temp = status

      close(iunit, status=status_temp, iostat=i)
      if( i /=  0) then
        write(6,'(/A,I3)') 'fclose: Could not close unit', iunit
        write(6,'(A,I3)')  '        Error code is', i
        NOstopNO              'fclose: Error while trying to close a file.'
      end if
    end subroutine fclose

    function intgrf(f, jri, jmtd, rmsh, dx, ntype, itype, gridf)

      implicit none

      real                 :: intgrf
      integer, intent(in)  :: itype, ntype, jmtd
      integer, intent(in)  :: jri(ntype)
      real   , intent(in)  :: dx(ntype), rmsh(jmtd,ntype)
      real   , intent(in)  :: gridf(jmtd, ntype)
      real   , intent(in)  :: f(*)
!     - local -
      integer              :: n
      real                 :: r1, h, a, x

      n  = jri(itype)
      r1 = rmsh(1, itype)
      h  = dx(itype)

      ! integral from 0 to r1 approximated by leading term in power series expansion
      ! f(r) = a+c*r**x
      if (f(1)* f(2) > 1d-10 .and. h > 0) then
        if(f(2) == f(1)) then
          intgrf = r1 * f(1)
        else
          x      = (f(3)-f(2)) / (f(2)-f(1))
          a      = (f(2)-x*f(1)) / (1-x)
          if( x > 0 ) then
            x      = log(x) / h
            intgrf = r1 * (f(1) + x * a) / (x + 1)
          else
            if(x > -1) write(6, '(A,ES9.1)') 'intgrf: Warning! Negative exponent x in extrapolation a+c*r**x:', x
            if(x <= -1) write(6, '(A,ES9.1)') 'intgrf: Negative exponent x in extrapolation a+c*r**x:', x
            WRITE(6, '(A)') 'intgrf: Switch to second integration approach for the integral from 0 to r1'

            x      = f(2) / f(1)
            if(x < 0) then
              if(x > -1) write(6,'(A,ES9.1)') 'intgrf: Warning! Negative exponent x in extrapolation c*r**x:',x
              if(x <= -1) write(6,'(A,ES9.1)') 'intgrf: Negative exponent x in extrapolation c*r**x:',x
              if(x <= -1) NOstopNO                'intgrf: Negative exponent x in extrapolation c*r**x as well as a+c*r**x'
            end if

            x      = log(x) / h
            intgrf = (r1 * f(1)) / (x + 1)
          end if

        end if
      else
        intgrf = 0
      end if
!     IF (f(1)*f(2).gt.1d-20) THEN
!        x = log ( f(2)/f(1) ) / h
!        IF(x.lt.-1) THEN
!          WRITE(0,'(A,ES9.1)') 'intgrf: Exponent x in extrapolation c*r**x below -1:',x
!          NOstopNO
!        END IF
!        IF(x.lt.-0.5) THEN
!          WRITE(0,'(A,ES9.1)') 'intgrf: Warning! Large negative exponent x in extrapolation c*r**x:',x
!        END IF
!        intgrf = r1*f(1) / (x+1)
!      ELSE
!        intgrf = 0
!      END IF

      ! integrate from r(1) to r(n) by multiplying with gridf
      intgrf = intgrf + dot_product(gridf(:,itype), f(:n))

      end function intgrf

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Initializes fast numerical integration intgrf (see below).

      subroutine intgrf_init(ntype, jmtd, jri, dx, rmsh, gridf)

      implicit none

      integer,              intent(in)  :: ntype,jmtd
      integer,              intent(in)  :: jri(ntype)
      real,                 intent(in)  :: dx(ntype),rmsh(jmtd,ntype)
      real,    allocatable, intent(out) :: gridf(:,:)

!     - local -
      integer                           :: i, j, itype, ierror
      integer                           :: n, nstep, n0 = 6
      integer, parameter                :: simpson(7) = (/41, 216, 27, 272, 27, 216, 41/)
      real                              :: r1, h, dr
      real                              :: r(7)
      real, parameter                   :: lagrange(7, 6) = reshape((/19087.,65112.,-46461., 37504.,-20211., 6312., -863., &
                                                      &-863.,25128., 46989.,-16256.,  7299.,-2088., 271., &
                                                      &271.,-2760., 30819., 37504., -6771., 1608., -191., &
                                                      &-191., 1608., -6771., 37504., 30819.,-2760., 271., &
                                                      &271.,-2088.,  7299.,-16256., 46989.,25128., -863., &
                                                      &-863., 6312.,-20211., 37504.,-46461.,65112.,19087./), &! The last row is actually never used.
                                                      &(/7,6/))

      n = jmtd
      allocate (gridf(n,ntype), stat=ierror)
      if( ierror /= 0 ) NOstopNO'util: error allocation gridf'

      gridf = 0

      do itype = 1,ntype

        n  = jri(itype)
        r1 = rmsh(1, itype)
        h  = dx(itype)

        nstep = (n - 1) / 6
        n0    = n - 6*nstep
        dr    = exp(h)

        ! Calculate Lagrange-integration coefficients from r(1) to r(n0)
        r(1)=r1
        if(n0 > 1) then
          do i=2,7
            r(i) = r(i - 1) * dr
          end do
          do i=1,7
            gridf(i, itype) = h / 60480 * r(i) * sum(lagrange(i, 1:n0-1))
          end do
          r(1)  = r(n0)
        end if

        ! Calculate Simpson-integration coefficients from r(n0) to r(n)
        do i = 1, nstep
          do j = 2, 7
            r(j) = r(j - 1) * dr
          end do
          do j = n0, n0 + 6
            gridf(j, itype) = gridf(j, itype) + h / 140 * r(j - n0 + 1) * simpson(j - n0 + 1)
          end do
          n0   = n0 + 6
          r(1) = r(7)
        end do

      end do

      end subroutine intgrf_init

      subroutine fitchk(f1,f2,n,av,rms,dmx)
!     ************************************************
!     compare functions f1 and f2
!     ************************************************
!     .. Scalar Arguments ..

      implicit none
      real av,dmx,rms
      integer n
!     ..
!     .. Array Arguments ..
      real f1(n),f2(n)
!     ..
!     .. Local Scalars ..
      real d
      integer i
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC abs,max,sqrt
!     ..
      av = 0.
      rms = 0.
      dmx = 0.
      DO 10 i = 1,n
         av = av + f1(i)
         d = (f1(i)-f2(i))**2
         dmx = max(d,dmx)
         rms = rms + d
   10 CONTINUE
      av = av/n
      IF (abs(av).LT.1.e-30) THEN
         rms = 0.
         dmx = 0.
         RETURN
      END IF
      rms = sqrt(rms/n)/av*100.
      dmx = sqrt(dmx)/av*100.
      RETURN
      end subroutine fitchk

!       Returns least common multiple of the integers iarr(1:n).
      function kgv(iarr,n)
        implicit none
        integer              :: kgv
        integer, intent(in)  :: n,iarr(n)
        logical              :: lprim(2:maxval(iarr))
        integer, allocatable :: prim(:),expo(:)
        integer              :: nprim,marr
        integer              :: i,j,ia,k
        ! Determine prime numbers
        marr  = maxval(iarr)
        lprim = .true.
        do i = 2,marr
          j = 2
          do while (i*j.le.marr)
            lprim(i*j) = .false.
            j          = j + 1
          enddo
        enddo
        nprim = count(lprim)
        allocate ( prim(nprim),expo(nprim) )
        j = 0
        do i = 2,marr
          if(lprim(i)) then
            j       = j + 1
            prim(j) = i
          endif
        enddo
        ! Determine least common multiple
        expo = 0
        do i = 1,n
          ia = iarr(i)
          if(ia.eq.0) cycle
          do j = 1,nprim
            k = 0
            do while(ia/prim(j)*prim(j).eq.ia)
              k  = k + 1
              ia = ia / prim(j)
            enddo
            expo(j) = max(expo(j),k)
          enddo
        enddo
        kgv = 1
        do j = 1,nprim
          kgv = kgv * prim(j)**expo(j)
        enddo
        deallocate ( prim,expo )
      end function kgv

!      Returns derivative of f in df.
! todo rename to Derivative1D
      subroutine derivative(f, itype, atoms, df)

      use m_types

      implicit none

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
      end subroutine derivative


!maps into first brillouin zone !TODO move into utils
  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst and Markus Betzinger
  !>
  !> @brief
  !> Maps into first brillouin zone.
  !>
  !> @details
  !>
  !>
  !> @param[in] input_type
  !--------------------------------------------------------------------------------------------------------------------------------------
function modulo1r(kpoint)
  implicit none
  real(8), intent(in) :: kpoint(3)
  real(8)             :: modulo1r(3)
  integer             :: i

  modulo1r = modulo (kpoint , 1d0)

  do i = 1,3
    if(abs(1-abs(modulo1r(i))).lt.1d-13) modulo1r(i) = 0d0
  enddo
end function modulo1r

    !--------------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculates the spherical-harmonics expansion-coefficients coordinates of a general MT function's gradient resolved in its
    !> scattering channels within the natural coordinate system.
    !>
    !> @details
    !> The MT gradient of a function expanded spherical harmonics scatters into two neighbouring l channels and respective outgoing
    !> m-channels. This routine gives the resulting output l and m quantum numbers and the respective gradient channels. Furthermore,
    !> this routine remains within the natural coordinate system.
    !>
    !> @note
    !> The input function should be dependent on the orbital quantum number l and the radial function index p and can contain a large
    !> and a small relativistic component.
    !>
    !> @todo
    !> review documentation
    !> move whole routine two spaces to the left
    !--------------------------------------------------------------------------------------------------------------------------------------
    subroutine CalcChannelsGrFlpNat( atoms, itype, nRadFun, flp1, flp2, delrFlp1, delrFlp2, grFlpChLout, grFlpChMout, grFlpCh1,    &
        & grFlpCh2 )

      use m_JPConstants, only : fpi
      use m_types, only : t_atoms
      use m_gaunt, only : gaunt1

      implicit none

      ! Type parameter
      type(t_atoms),               intent(in)  :: atoms

      ! Scalar parameter
      integer,                     intent(in)  :: itype

      ! Array parameter
      ! nRadFun     : number of radial functions per l and type
      ! flp1        : 1st(large) component of function with l and p(p loops over the radial function types)
      ! flp2        : 2nd(small) component of function with l and p(p loops over the radial function types)
      ! delrFlp1    : Radial derivative of the 1st component of function with l and p(p loops over the radial function types)
      ! delrFlp2    : Radial derivative of the 2nd component of function with l and p(p loops over the radial function types)
      ! grFlpChLout : Contains the resulting outgoing l of a respective scattering channel of the gradient
      ! grFlpChMout : Contains the resulting outgoing m of a respective scattering channel of the gradient
      ! grFlpCh1    : Contains the 1st (large) component of the gradient scattering channels of the input function
      ! grFlpCh2    : Contains the 2nd (small) component of the gradient scattering channels of the input function
      integer,                     intent(in)  :: nRadFun(0:, :)
      real,                        intent(in)  :: flp1(:, :, 0:)
      real,                        intent(in)  :: flp2(:, :, 0:)
      real,                        intent(in)  :: delrFlp1(:, :, 0:)
      real,                        intent(in)  :: delrFlp2(:, :, 0:)
      integer,                     intent(out) :: grFlpChLout(:, 0:)
      integer,                     intent(out) :: grFlpChMout(-atoms%lmaxd:, -1:)
      real,                        intent(out) :: grFlpCh1(:, :, :, -1:)
      real,                        intent(out) :: grFlpCh2(:, :, :, -1:)

      ! Local Variables:
      !
      ! pfac     : contains sqrt(4 pi / 3)
      ! lm          : encodes oqn_l and mqn_m
      ! idirec      : runs over 3 directions the atom can be displaced to
      ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
      ! imesh       : runs over mesh points of current grid
      ! mqn_m       : magnetic quantum number m
      ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
      ! oqn_l       : orbital quantum number l

      ! Local Scalar Variables
      real                                     :: pfac
      integer                                  :: oqn_l
      integer                                  :: mqn_m
      integer                                  :: iradf
      integer                                  :: lmp
      real                                     :: gauntCoeff
      integer                                  :: imesh
      integer                                  :: mqn_m2Pr


      ! Local Array Variables

      grFlpCh1(:, :, :, :) = 0.
      grFlpCh2(:, :, :, :) = 0.


      pfac = sqrt( fpi / 3. )

      ! Precalculate L-output channels
      do oqn_l = 0, atoms%lmax(itype)
        grFlpChLout(1, oqn_l) = oqn_l + 1
        grFlpChLout(2, oqn_l) = oqn_l - 1
      end do ! oqn_l
      ! Set this lout = 0 - 1 to abnormal value as not allowed
      grFlpChLout(2, 0) = -9999

      ! Precalculate M-output channels
      do mqn_m2Pr = -1, 1
        do mqn_m = -atoms%lmax(itype), atoms%lmax(itype)
          grFlpChMout(mqn_m, mqn_m2Pr) = mqn_m - mqn_m2Pr
        end do ! mqn_m
      end do ! mqn_m2Pr

      do mqn_m2Pr = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1

              ! scattering channel (l + 1)
              if ( (abs(mqn_m - mqn_m2Pr) <= oqn_l + 1) .and. (oqn_l < atoms%lmax(itype))) then
                gauntCoeff = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmax(itype) )
                do imesh = 1, atoms%jri(itype)
                  ! Consider large and small relativistic components of radial solution
                  grFlpCh1(imesh, 1, lmp, mqn_m2Pr) = grFlpCh1(imesh, 1, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                                 & * (delrFlp1(imesh, iradf, oqn_l) -  oqn_l      * flp1(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
                  grFlpCh2(imesh, 1, lmp, mqn_m2Pr) = grFlpCh2(imesh, 1, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                                 & * (delrFlp2(imesh, iradf, oqn_l) -  oqn_l      * flp2(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
                enddo ! imesh
              endif ! scattering channel (l + 1)

              ! scattering channel (l - 1)
              ! This condition ensures that oqn_l = 0 is not accepted due to the emerging false condition 0 or 1 <= -1
              if ( ( abs(mqn_m - mqn_m2Pr) <= oqn_l - 1 ) ) then
                gauntCoeff = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmax(itype))
                do imesh = 1, atoms%jri(itype)
                  ! Consider large and small relativistic components of radial solution
                  grFlpCh1(imesh, 2, lmp, mqn_m2Pr) = grFlpCh1(imesh, 2, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                                 & * (delrFlp1(imesh, iradf, oqn_l) + (oqn_l + 1) * flp1(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
                  grFlpCh2(imesh, 2, lmp, mqn_m2Pr) = grFlpCh2(imesh, 2, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                                 & * (delrFlp2(imesh, iradf, oqn_l) + (oqn_l + 1) * flp2(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
                enddo ! imesh
              endif ! scattering channel (l - 1)

            end do ! p
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2Pr

    end subroutine CalcChannelsGrFlpNat

    !--------------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculates the spherical-harmonics expansion-coefficients coordinates of a general MT function's gradient resolved in its
    !> scattering channels within the natural coordinate system.
    !>
    !> @details
    !> The MT gradient of a function expanded spherical harmonics scatters into two neighbouring l channels and respective outgoing
    !> m-channels. This routine gives the resulting output l and m quantum numbers and the respective gradient channels. Furthermore,
    !> this routine remains within the natural coordinate system.
    !>
    !> @note
    !> The input function should be dependent on the orbital quantum number l and the radial function index p and can contain a large
    !> and a small relativistic component.
    !>
    !> @todo
    !> review documentation
    !> move whole routine two spaces to the left
    !--------------------------------------------------------------------------------------------------------------------------------------
    subroutine CalcChannelsGrR2FlpNat( atoms, itype, nRadFun, flp1, flp2, delrFlp1, delrFlp2, grFlpChLout, grFlpChMout, grFlpCh1,    &
        & grFlpCh2 )

      use m_JPConstants, only : fpi
      use m_types, only : t_atoms
      use m_gaunt, only : gaunt1

      implicit none

      ! Type parameter
      type(t_atoms),               intent(in)  :: atoms

      ! Scalar parameter
      integer,                     intent(in)  :: itype

      ! Array parameter
      ! nRadFun     : number of radial functions per l and type
      ! flp1        : 1st(large) component of function with l and p(p loops over the radial function types)
      ! flp2        : 2nd(small) component of function with l and p(p loops over the radial function types)
      ! delrFlp1    : Radial derivative of the 1st component of function with l and p(p loops over the radial function types)
      ! delrFlp2    : Radial derivative of the 2nd component of function with l and p(p loops over the radial function types)
      ! grFlpChLout : Contains the resulting outgoing l of a respective scattering channel of the gradient
      ! grFlpChMout : Contains the resulting outgoing m of a respective scattering channel of the gradient
      ! grFlpCh1    : Contains the 1st (large) component of the gradient scattering channels of the input function
      ! grFlpCh2    : Contains the 2nd (small) component of the gradient scattering channels of the input function
      integer,                     intent(in)  :: nRadFun(0:, :)
      real,                        intent(in)  :: flp1(:, :, 0:)
      real,                        intent(in)  :: flp2(:, :, 0:)
      real,                        intent(in)  :: delrFlp1(:, :, 0:)
      real,                        intent(in)  :: delrFlp2(:, :, 0:)
      integer,                     intent(out) :: grFlpChLout(:, 0:)
      integer,                     intent(out) :: grFlpChMout(-atoms%lmaxd:, -1:)
      real,                        intent(out) :: grFlpCh1(:, :, :, -1:)
      real,                        intent(out) :: grFlpCh2(:, :, :, -1:)

      ! Local Variables:
      !
      ! pfac     : contains sqrt(4 pi / 3)
      ! lm          : encodes oqn_l and mqn_m
      ! idirec      : runs over 3 directions the atom can be displaced to
      ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
      ! imesh       : runs over mesh points of current grid
      ! mqn_m       : magnetic quantum number m
      ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
      ! oqn_l       : orbital quantum number l

      ! Local Scalar Variables
      real                                     :: pfac
      integer                                  :: oqn_l
      integer                                  :: mqn_m
      integer                                  :: iradf
      integer                                  :: lmp
      real                                     :: gauntCoeff
      integer                                  :: imesh
      integer                                  :: mqn_m2Pr


      ! Local Array Variables

      grFlpCh1(:, :, :, :) = 0.
      grFlpCh2(:, :, :, :) = 0.


      pfac = sqrt( fpi / 3. )

      ! Precalculate L-output channels
      do oqn_l = 0, atoms%lmax(itype)
        grFlpChLout(1, oqn_l) = oqn_l + 1
        grFlpChLout(2, oqn_l) = oqn_l - 1
      end do ! oqn_l
      ! Set this lout = 0 - 1 to abnormal value as not allowed
      grFlpChLout(2, 0) = -9999

      ! Precalculate M-output channels
      do mqn_m2Pr = -1, 1
        do mqn_m = -atoms%lmax(itype), atoms%lmax(itype)
          grFlpChMout(mqn_m, mqn_m2Pr) = mqn_m - mqn_m2Pr
        end do ! mqn_m
      end do ! mqn_m2Pr

      do mqn_m2Pr = -1, 1
        lmp = 0
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            do iradf = 1, nRadFun(oqn_l, itype)
              lmp = lmp + 1

              ! scattering channel (l + 1)
              if ( (abs(mqn_m - mqn_m2Pr) <= oqn_l + 1) .and. (oqn_l < atoms%lmax(itype)) ) then
                gauntCoeff = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmax(itype) )
                do imesh = 1, atoms%jri(itype)
                  ! Consider large and small relativistic components of radial solution
                  grFlpCh1(imesh, 1, lmp, mqn_m2Pr) = grFlpCh1(imesh, 1, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                                 & * (delrFlp1(imesh, iradf, oqn_l) -  ((oqn_l + 2)     * flp1(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype)))
                  grFlpCh2(imesh, 1, lmp, mqn_m2Pr) = grFlpCh2(imesh, 1, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                                 & * (delrFlp2(imesh, iradf, oqn_l) -  ((oqn_l + 2)     * flp2(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype)))
                enddo ! imesh
              endif ! scattering channel (l + 1)

              ! scattering channel (l - 1)
              ! This condition ensures that oqn_l = 0 is not accepted due to the emerging false condition 0 or 1 <= -1
              if ( ( abs(mqn_m - mqn_m2Pr) <= oqn_l - 1 ) ) then
                gauntCoeff = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmax(itype))
                do imesh = 1, atoms%jri(itype)
                  ! Consider large and small relativistic components of radial solution
                  grFlpCh1(imesh, 2, lmp, mqn_m2Pr) = grFlpCh1(imesh, 2, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                                 & * (delrFlp1(imesh, iradf, oqn_l) + (oqn_l - 1) * flp1(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
                  grFlpCh2(imesh, 2, lmp, mqn_m2Pr) = grFlpCh2(imesh, 2, lmp, mqn_m2Pr) + pfac * (-1)**mqn_m2Pr * gauntCoeff &
                                 & * (delrFlp2(imesh, iradf, oqn_l) + (oqn_l - 1) * flp2(imesh, iradf, oqn_l) / atoms%rmsh(imesh, itype))
                enddo ! imesh
              endif ! scattering channel (l - 1)

            end do ! p
          end do ! mqn_m
        end do ! oqn_l
      end do ! mqn_m2Pr

    end subroutine CalcChannelsGrR2FlpNat

    !--------------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculates the spherical-harmonics expansion-coefficients coordinates of MT gradient channels determined with the routine
    !> CalcChannelsGrFlpNat and remains within the natural coordinate system.
    !>
    !> @details
    !> The MT gradient of a function gradient expanded in spherical harmonics scatters into four neighbouring l channels and
    !> respective outgoing m-channels. This routine gives the resulting output l and m quantum numbers and the respective double
    !> gradient channels. Furthermore, this routine remains within the natural coordinate system.
    !>
    !> @note
    !> The input function should have been generated with the calcChannelsGrFlpNat routine can contain a large and a small
    !> relativistic component.
    !>
    !> @todo
    !> review documentation
    !> move whole routine two spaces to the left
    !--------------------------------------------------------------------------------------------------------------------------------------
    subroutine CalcChannelsGrGrtFlpNat( atoms, itype, lmpMax, nRadFun, grFlpChLout, grFlpChMout, grFlpCh1, grFlpCh2, delrGrFlpCh1, &
        & delrGrFlpCh2, grGrtFlpChLout, grGrtFlpChMout, grGrtFlpCh1, grGrtFlpCh2 )

      use m_jpConstants, only : fpi
      use m_types, only : t_atoms
      use m_gaunt, only : gaunt1

      implicit none

      ! Type parameter
      type(t_atoms),               intent(in)  :: atoms

      ! Scalar parameter
      integer,                     intent(in)  :: itype
      integer,                     intent(in)  :: lmpMax

      ! Array parameter
      ! nRadFun        : number of radial functions per l and type
      ! grFlpChLout    : Contains the resulting outgoing l of a respective scattering channel of the input gradient
      ! grFlpChMout    : Contains the resulting outgoing m of a respective scattering channel of the input gradient
      ! grFlpCh1       : Contains the 1st (large) component of the gradient scattering channels of the input function
      ! grFlpCh2       : Contains the 2nd (small) component of the gradient scattering channels of the input function
      ! delrgrFlpCh1   : Radial derivative of the 1st (large) component of the gradient scattering channels of the input function
      ! delrgrFlpCh2   : Radial derivative of the 2nd (small) component of the gradient scattering channels of the input function
      ! grGrtFlpChLout : Resulting output l channels to the scattering channels given in grGrtFlpCh1/2
      ! grGrtFlpChMout : Resulting output m channels to the scattering channels given in grGrtFlpCh1/2
      ! grGrtFlpCh1    : 1st(large) component of the double gradient of a function given in spherical harmonics
      ! grGrtFlpCh2    : 2nd(small) component of the double gradient of a function given in spherical harmonics
      ! Array parameter
      integer,                     intent(in)  :: nRadFun(0:, :)
      integer,                     intent(in)  :: grFlpChLout(:, 0:)
      integer,                     intent(in)  :: grFlpChMout(-atoms%lmaxd:, -1:)
      real,                        intent(in)  :: grFlpCh1(:, :, :, -1:)
      real,                        intent(in)  :: grFlpCh2(:, :, :, -1:)
      real,                        intent(in)  :: delrGrFlpCh1(:, :, :, -1:)
      real,                        intent(in)  :: delrGrFlpCh2(:, :, :, -1:)
      integer,                     intent(out) :: grGrtFlpChLout(:, 0:)
      integer,                     intent(out) :: grGrtFlpChMout(-atoms%lmaxd:, -1:, -1:)
      real,                        intent(out) :: grGrtFlpCh1(:, :, :, -1:, -1:)
      real,                        intent(out) :: grGrtFlpCh2(:, :, :, -1:, -1:)

      ! Scalar Variables
      real                                     :: pfac
      integer                                  :: mqn_m2PrC
      integer                                  :: mqn_m2PrR
      integer                                  :: lmp
      integer                                  :: oqn_l
      integer                                  :: mqn_m
      integer                                  :: iradf
      integer                                  :: iChGrGrt
      integer                                  :: iChGr
      integer                                  :: oqn_l3Pr
      integer                                  :: mqn_m3Pr
      real                                     :: gauntCoeff
      integer                                  :: imesh

      ! Array Variables

      pfac = sqrt( fpi / 3 )

      do mqn_m2PrC = -1, 1
        do mqn_m2PrR = -1, 1
          ! Precalculate output magnetic quantum number m
          do mqn_m = -atoms%lmax(itype), atoms%lmax(itype)
            ! We calculate grad grad^T. Hence, the incoming gradient is a column-vector
            mqn_m3Pr = grFlpChMout(mqn_m, mqn_m2PrC)
            ! The resulting m of the double gradient is calculated from the output m of the simple gradient.
            grGrtFlpChMout(mqn_m, mqn_m2PrR, mqn_m2PrC) = mqn_m3Pr - mqn_m2PrR
          end do ! mqn_m
          lmp = 0
          ! If we loop over the l, m quantum numbers of the original function and the scattering channels of its gradient, we
          ! account for all contributions.
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              mqn_m3Pr = grFlpChMout(mqn_m, mqn_m2PrC)
              do iradf = 1, nRadFun(oqn_l, itype)
                lmp = lmp + 1
                iChGrGrt = 0
                do iChGr = 1, 2
                  oqn_l3Pr = grFlpChLout(iChGr, oqn_l)
                  if ((oqn_l3Pr < 0 ) .or. (oqn_l3Pr > atoms%lmax(itype) + 1)) cycle
                  iChGrGrt = iChGrGrt + 1
                  ! Still for the Gaunt coefficient conditions we need the resulting m and l quantum numbers of the simple gradient
                  ! scattering channels
                  if ( ( abs(mqn_m3Pr - mqn_m2PrR) <= oqn_l3Pr + 1 ) .and. ( abs(mqn_m3Pr) <= oqn_l3Pr ) ) then
                    grGrtFlpChLout(iChGrGrt, oqn_l) = oqn_l3Pr + 1
                    gauntCoeff = Gaunt1( oqn_l3Pr + 1, oqn_l3Pr, 1, mqn_m3Pr - mqn_m2PrR, mqn_m3Pr, -mqn_m2PrR, atoms%lmax(itype) + 2)
                    do imesh = 1, atoms%jri(itype)
                      ! We do not need two quantities here for large and small component but to make it consistent we use two!
                      grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                  & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) &
                                                  & - oqn_l3Pr * grFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                      grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                  & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) &
                                                  & - oqn_l3Pr * grFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                    end do ! imesh
                  end if ! l + 1-scattering channel

                  iChGrGrt = iChGrGrt + 1
                  ! Still for the Gaunt coefficient conditions we need the resulting m and l quantum numbers of the simple gradient
                  ! scattering channels
                  if ( (abs(mqn_m3Pr - mqn_m2PrR) <= oqn_l3Pr - 1) .and. ( abs(mqn_m3Pr) <= oqn_l3Pr ) ) then
                    grGrtFlpChLout(iChGrGrt, oqn_l) = oqn_l3Pr - 1
                    gauntCoeff = Gaunt1( oqn_l3Pr - 1, oqn_l3Pr, 1, mqn_m3Pr - mqn_m2PrR, mqn_m3Pr, -mqn_m2PrR, atoms%lmax(itype) + 1)
                    do imesh = 1, atoms%jri(itype)
                      ! We do not need two quantities here for large and small component but to make it consistent we use two!
                      grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                  & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) &
                                                  & + (oqn_l3Pr + 1) * grFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                      grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                  & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) &
                                                  & + (oqn_l3Pr + 1) * grFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                    end do ! imesh
                  end if ! l - 1-scattering channel
                end do ! iChGr
              end do ! iradf
            end do ! mqn_m
          end do ! oqn_l
        end do ! mqn_m2PrR
      end do ! mqn_m2PrC

    end subroutine CalcChannelsGrGrtFlpNat

    !--------------------------------------------------------------------------------------------------------------------------------------
    !> @author
    !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
    !>
    !> @brief
    !> Calculates the spherical-harmonics expansion-coefficients coordinates of MT gradient channels determined with the routine
    !> CalcChannelsGrFlpNat and remains within the natural coordinate system.
    !>
    !> @details
    !> The MT gradient of a function gradient expanded in spherical harmonics scatters into four neighbouring l channels and
    !> respective outgoing m-channels. This routine gives the resulting output l and m quantum numbers and the respective double
    !> gradient channels. Furthermore, this routine remains within the natural coordinate system.
    !>
    !> @note
    !> The input function should have been generated with the calcChannelsGrFlpNat routine can contain a large and a small
    !> relativistic component.
    !>
    !> @todo
    !> review documentation
    !> move whole routine two spaces to the left
    !--------------------------------------------------------------------------------------------------------------------------------------
    subroutine CalcChannelsGrGrtR2FlpNat( atoms, itype, lmpMax, nRadFun, grFlpChLout, grFlpChMout, grFlpCh1, grFlpCh2, delrGrFlpCh1, &
        & delrGrFlpCh2, grGrtFlpChLout, grGrtFlpChMout, grGrtFlpCh1, grGrtFlpCh2 )

      use m_jpConstants, only : fpi
      use m_types, only : t_atoms
      use m_gaunt, only : gaunt1

      implicit none

      ! Type parameter
      type(t_atoms),               intent(in)  :: atoms

      ! Scalar parameter
      integer,                     intent(in)  :: itype
      integer,                     intent(in)  :: lmpMax

      ! Array parameter
      ! nRadFun        : number of radial functions per l and type
      ! grFlpChLout    : Contains the resulting outgoing l of a respective scattering channel of the input gradient
      ! grFlpChMout    : Contains the resulting outgoing m of a respective scattering channel of the input gradient
      ! grFlpCh1       : Contains the 1st (large) component of the gradient scattering channels of the input function
      ! grFlpCh2       : Contains the 2nd (small) component of the gradient scattering channels of the input function
      ! delrgrFlpCh1   : Radial derivative of the 1st (large) component of the gradient scattering channels of the input function
      ! delrgrFlpCh2   : Radial derivative of the 2nd (small) component of the gradient scattering channels of the input function
      ! grGrtFlpChLout : Resulting output l channels to the scattering channels given in grGrtFlpCh1/2
      ! grGrtFlpChMout : Resulting output m channels to the scattering channels given in grGrtFlpCh1/2
      ! grGrtFlpCh1    : 1st(large) component of the double gradient of a function given in spherical harmonics
      ! grGrtFlpCh2    : 2nd(small) component of the double gradient of a function given in spherical harmonics
      ! Array parameter
      integer,                     intent(in)  :: nRadFun(0:, :)
      integer,                     intent(in)  :: grFlpChLout(:, 0:)
      integer,                     intent(in)  :: grFlpChMout(-atoms%lmaxd:, -1:)
      real,                        intent(in)  :: grFlpCh1(:, :, :, -1:)
      real,                        intent(in)  :: grFlpCh2(:, :, :, -1:)
      real,                        intent(in)  :: delrGrFlpCh1(:, :, :, -1:)
      real,                        intent(in)  :: delrGrFlpCh2(:, :, :, -1:)
      integer,                     intent(out) :: grGrtFlpChLout(:, 0:)
      integer,                     intent(out) :: grGrtFlpChMout(-atoms%lmaxd:, -1:, -1:)
      real,                        intent(out) :: grGrtFlpCh1(:, :, :, -1:, -1:)
      real,                        intent(out) :: grGrtFlpCh2(:, :, :, -1:, -1:)

      ! Scalar Variables
      real                                     :: pfac
      integer                                  :: mqn_m2PrC
      integer                                  :: mqn_m2PrR
      integer                                  :: lmp
      integer                                  :: oqn_l
      integer                                  :: mqn_m
      integer                                  :: iradf
      integer                                  :: iChGrGrt
      integer                                  :: iChGr
      integer                                  :: oqn_l3Pr
      integer                                  :: mqn_m3Pr
      real                                     :: gauntCoeff
      integer                                  :: imesh

      ! Array Variables

      pfac = sqrt( fpi / 3 )

      do mqn_m2PrC = -1, 1
        do mqn_m2PrR = -1, 1
          ! Precalculate output magnetic quantum number m
          do mqn_m = -atoms%lmax(itype), atoms%lmax(itype)
            ! We calculate grad grad^T. Hence, the incoming gradient is a column-vector
            mqn_m3Pr = grFlpChMout(mqn_m, mqn_m2PrC)
            ! The resulting m of the double gradient is calculated from the output m of the simple gradient.
            grGrtFlpChMout(mqn_m, mqn_m2PrR, mqn_m2PrC) = mqn_m3Pr - mqn_m2PrR
          end do ! mqn_m
          lmp = 0
          ! If we loop over the l, m quantum numbers of the original function and the scattering channels of its gradient, we
          ! account for all contributions.
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
              mqn_m3Pr = grFlpChMout(mqn_m, mqn_m2PrC)
              do iradf = 1, nRadFun(oqn_l, itype)
                lmp = lmp + 1
                iChGrGrt = 0
                do iChGr = 1, 2
                  oqn_l3Pr = grFlpChLout(iChGr, oqn_l)
                  if ((oqn_l3Pr < 0 ) .or. (oqn_l3Pr > atoms%lmax(itype) + 1)) cycle
                  iChGrGrt = iChGrGrt + 1
                  ! Still for the Gaunt coefficient conditions we need the resulting m and l quantum numbers of the simple gradient
                  ! scattering channels
                  if ( ( abs(mqn_m3Pr - mqn_m2PrR) <= oqn_l3Pr + 1 ) .and. ( abs(mqn_m3Pr) <= oqn_l3Pr ) ) then
                    grGrtFlpChLout(iChGrGrt, oqn_l) = oqn_l3Pr + 1
                    gauntCoeff = Gaunt1( oqn_l3Pr + 1, oqn_l3Pr, 1, mqn_m3Pr - mqn_m2PrR, mqn_m3Pr, -mqn_m2PrR, atoms%lmax(itype) + 2)
                    do imesh = 1, atoms%jri(itype)
                      ! We do not need two quantities here for large and small component but to make it consistent we use two!
                      grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                  & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) &
                                                  & - ((oqn_l3Pr + 2) * grFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype)))
                      grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                  & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) &
                                                  & - ((oqn_l3Pr + 2) * grFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype)))
                    end do ! imesh
                  end if ! l + 1-scattering channel

                  iChGrGrt = iChGrGrt + 1
                  ! Still for the Gaunt coefficient conditions we need the resulting m and l quantum numbers of the simple gradient
                  ! scattering channels
                  if ( (abs(mqn_m3Pr - mqn_m2PrR) <= oqn_l3Pr - 1) .and. ( abs(mqn_m3Pr) <= oqn_l3Pr ) ) then
                    grGrtFlpChLout(iChGrGrt, oqn_l) = oqn_l3Pr - 1
                    gauntCoeff = Gaunt1( oqn_l3Pr - 1, oqn_l3Pr, 1, mqn_m3Pr - mqn_m2PrR, mqn_m3Pr, -mqn_m2PrR, atoms%lmax(itype) + 1)
                    do imesh = 1, atoms%jri(itype)
                      ! We do not need two quantities here for large and small component but to make it consistent we use two!
                      grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh1(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                  & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) &
                                                  & + (oqn_l3Pr - 1) * grFlpCh1(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                      grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) = grGrtFlpCh2(imesh, iChGrGrt, lmp, mqn_m2PrR, mqn_m2PrC) &
                                                  & + pfac * (-1)**mqn_m2PrR * gauntCoeff * (delrGrFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) &
                                                  & + (oqn_l3Pr - 1) * grFlpCh2(imesh, iChGr, lmp, mqn_m2PrC) / atoms%rmsh(imesh, itype))
                    end do ! imesh
                  end if ! l - 1-scattering channel
                end do ! iChGr
              end do ! iradf
            end do ! mqn_m
          end do ! oqn_l
        end do ! mqn_m2PrR
      end do ! mqn_m2PrC

    end subroutine CalcChannelsGrGrtR2FlpNat

    ! Returns the inverse of a matrix calculated by finding the LU decomposition. Depends on LAPACK. Taken from Fortran Wiki 25/08/19
    function inv(A) result(Ainv)

      implicit none

      ! Array parameters
      complex, intent(in) :: A(:, :)

      ! Array variables
      complex             :: Ainv(size(A,1), size(A, 2))
      complex             :: work(size(A, 1)) ! work array for LAPACK
      integer          :: ipiv(size(A, 1)) ! pivot indices

      ! Scalar variables
      integer          :: n
      integer          :: info

      ! External porcedures defined in LAPACK
      external ZGETRF
      external ZGETRI

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
      Ainv = A
      n = size(A, 1)

      ! DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
      call Zgetrf(n, n, Ainv, n, ipiv, info)

      if (info /= 0) then
        NOstopNO'Matrix is numerically singular'
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
      call Zgetri(n, Ainv, n, ipiv, work, n, info)

      if (info /= 0) then
        NOstopNO'Matrix inversion failed!'
      end if

    end function inv

    ! taken from old fleur, maybe a better version is in the new fleur
    SUBROUTINE sphpts(p,n,r,pos)
!     *******************************************************
!     generates points on sphere at pos with radius r
!     e. wimmer     feb. 1980
!     modified to give a better distribution of points
!     m. weinert    jan. 1982
!     *******************************************************

      USE m_constants, ONLY : pimach
      IMPLICIT NONE
!     .. Scalar Arguments ..
      REAL r
      INTEGER n
!     ..
!     .. Array Arguments ..
      REAL p(3,n),pos(3)
!     ..
!     .. Local Scalars ..
      REAL phi,pi2,t,tc,x,xr,y,yr,z
      INTEGER i,j
!     ..
!     .. External Functions ..
!      REAL qranf
!      EXTERNAL qranf
!     ..
!     .. Intrinsic Functions ..
!      INTRINSIC cos,sin,sqrt
!     ..
      pi2 = 2.0*pimach()
      j = 0
      xr = sqrt(13.e0)
      yr = sqrt(7.e0)
      DO i = 1,n
         tc = 2.e0*qranf(xr,j) - 1.e0
         phi = pi2*qranf(yr,j)
         t = sqrt(1.e0-tc*tc)
         x = t*cos(phi)
         y = t*sin(phi)
         z = tc
         p(1,i) = r*x + pos(1)
         p(2,i) = r*y + pos(2)
         p(3,i) = r*z + pos(3)
      end do
      RETURN
    END

      REAL FUNCTION qranf(x,j)
!     **********************************************************
!     quasi random generator in the interval (0.,1.)
!     **********************************************************

      IMPLICIT NONE
!     .. Scalar Arguments ..
      REAL x
      INTEGER j
!     ..
!     .. Intrinsic Functions ..
!     INTRINSIC aint
!     ..
      j = j + 1
      qranf = j*x
      qranf = qranf - aint(qranf)
      RETURN
      END

  ! calculates a 2argument volume integral in the IR
  subroutine Calc2ArgIntIR(cell, ngdp, f, w_g, integral)

    use m_types

    implicit none

    ! Type parameter
    type(t_cell), intent(in)  :: cell

    ! Scalar parameter
    integer,      intent(in)  :: ngdp
    complex,      intent(out) :: integral

    ! Array parameter
    complex,      intent(in)  :: f(:)
    complex,      intent(in)  :: w_g(:)


    ! The complex conjugation is done implicetly for the function passed as first argument of dot_product
    integral = cell%omtil * dot_product( f(:ngdp), w_g(:ngdp) )

  end subroutine Calc2ArgIntIR

  ! Calculates a 2 argument voluyme integral in the MT
  subroutine Calc2ArgCmplxIntMT( atoms, itype, f, g, integral )

    use m_types
    !use m_intgr, only : intgr3
    use m_intgr, only : intgr3LinIntp

    implicit none

    ! Type parameters
    type(t_atoms),             intent(in)  :: atoms

    ! Scalar parameters
    integer,                   intent(in)  :: itype
    complex,                   intent(out) :: integral

    ! Array parameters
    complex,                   intent(in)  :: f(:)
    complex,                   intent(in)  :: g(:)

    ! Scalar variables
    integer                                :: imesh
    real                                   :: integralReal
    real                                   :: integralImag

    ! Array variables
    real,          allocatable             :: intgrdR(:)
    real,          allocatable             :: intgrdI(:)


    allocate( intgrdR(atoms%jmtd), intgrdI(atoms%jmtd) )

    intgrdR(:) = 0.
    intgrdI(:) = 0.
    integral = cmplx(0., 0.)

    ! The intgr3LinIntp subroutine only accepts real quantities, so we split the integral into real and imaginary part.
    do imesh = 1, atoms%jri(itype)
      intgrdR(imesh) = real( atoms%rmsh(imesh, itype)**2 * conjg(f(imesh)) * g(imesh) )
      intgrdI(imesh) = aimag( atoms%rmsh(imesh, itype)**2 * conjg(f(imesh)) * g(imesh) )
    end do ! imesh
    !call intgr3(intgrdR, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), integralReal)
    !call intgr3(intgrdI, atoms%rmsh(1,itype), atoms%dx(itype), atoms%jri(itype), integralImag)
    call intgr3LinIntp( intgrdR, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralReal, 1 )
    call intgr3LinIntp( intgrdI, atoms%rmsh(1, itype), atoms%dx(itype), atoms%jri(itype), integralImag, 1 )

    integral = integral + cmplx( integralReal, integralImag )

  end subroutine Calc2ArgCmplxIntMT

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine calculates the muffin-tin derivative according to formula in PhD thesis Aaron.
  !> @details
  !>
  !> @param[in] atomsT     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in] latharT
  !> @param[in] clnu_atom
  !> @param[in] nmem_atom
  !> @param[in] mlh_atom
  !> @param[in] rho0MT
  !> @param[in] gradrho0MT
  !--------------------------------------------------------------------------------------------------------------------------------------
  !TODO put this into utils module
  !TODO rename to CalcGradFuncMT
  subroutine calcGrFinLH(atoms, sphhar, clnu_atom, nmem_atom, mlh_atom, rho0MT, grRho0MT)

    use m_JPConstants, only : Tmatrix, fpi, compPhon
    use m_gaunt, only : gaunt1
    use m_types

    implicit none

    type(t_atoms),                     intent(in)  :: atoms
    type(t_sphhar),                    intent(in)  :: sphhar
    complex,                           intent(in)  :: clnu_atom(:, 0:, :)
    integer,                           intent(in)  :: nmem_atom(0:, :)
    integer,                           intent(in)  :: mlh_atom(:, 0:, :)
    real,                              intent(in)  :: rho0MT(:, 0:, :)
    complex,        allocatable,       intent(out) :: grRho0MT(:, :, :, :)

    ! Local Variables:
    !
    ! sqr4pi3     : contains sqrt(4 pi / 3)
    ! iatom       : runs over all atoms
    ! lm          : encodes oqn_l and mqn_m
    ! itype       : runs over all atom types
    ! idirec      : runs over 3 directions the atom can be displaced to
    ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
    ! tempGaunt2  : auxillary variable to store a Gaunt coefficient
    ! imesh       : runs over mesh points of current grid
    ! rl2         : stores R^(l + 2)
    ! fint        : stores integral
    ! ll1         : auxillary variable to calculate lm
    ! mqn_m       : magnetic quantum number m
    ! sk3r        : stores argument of spherical Bessel function
    ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
    ! oqn_l       : orbital quantum number l
    ! iGvec       : indexes current G-vector
    ! grRho0MT  : stores the gradient of the unperturbed density in the MTs (equation 7.59, PhD thesis Aaron Klüppelberg)
    ! gradrho0PWi : stores the the second part of equation 7.58 using equation 7.60 from PhD thesis Aaron Klüppelberg
    ! pylm        : contains 4π i i^l G / |G| exp(i G τ)  Y*_lm(G / |G|)
    ! sbes        : stores the Bessel function
    ! cil         : stores everything except for pylm of the result of this routine
    ! f           : stores the integrand of the first term in (7.58, PhD thesis Aaron Klüppelberg)

    ! Local Scalar Variables
    real                                           :: pfac
    real                                           :: tGaunt1, tGaunt2
    integer                                        :: itype
    integer                                        :: ieqat
    integer                                        :: iatom
    integer                                        :: imesh
    integer                                        :: mqn_m
    integer                                        :: oqn_l
    integer                                        :: mqn_mpp
    integer                                        :: lm
    integer                                        :: lm_temp
    integer                                        :: lmMinus
    integer                                        :: lmPlus
    integer                                        :: symType
    integer                                        :: lh
    integer                                        :: mem

    ! Local Array Variables
    real                                           :: rDerRho0MT(atoms%jmtd)
    complex,        allocatable                    :: grRho0MTNat(:, :, :, :)


    !todo is r^2 rho numerically more stable?
    allocate( grRho0MT(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    grRho0MT = cmplx(0., 0.)
    allocate( grRho0MTNat(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    grRho0MTNat = cmplx(0., 0.)

    pfac = sqrt( fpi / 3. )
    do mqn_mpp = -1, 1
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          symType = atoms%ntypsy(iatom)
          do lh = 0, sphhar%nlh(symType)
            oqn_l = sphhar%llh(lh, symType)
            do mem = 1, nmem_atom(lh, iatom)
              mqn_m = mlh_atom(mem, lh, iatom)
              ! oqn_l - 1 to l, so oqn_l should be < lmax not <=lmax
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < atoms%lmax(itype)) ) then
                lm = ( oqn_l + 1 ) * ( oqn_l + 2 ) + 1 + mqn_m - mqn_mpp
                call Derivative( rho0MT(:, lh, itype), itype, atoms, rDerRho0MT )
                tGaunt1 = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) = grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp * tGaunt1 * (rDerRho0MT(imesh) * clnu_atom(mem, lh, iatom) &
                                                               &- (oqn_l * rho0MT(imesh, lh, itype) * clnu_atom(mem, lh, iatom) / atoms%rmsh(imesh, itype)))
                enddo ! imesh
              endif
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
                if ( oqn_l - 1 == -1 ) then
                  NOstopNO'oqn_l too low'
                end if
                lm = (oqn_l - 1) * oqn_l + 1 + mqn_m - mqn_mpp
                call Derivative( rho0MT(:, lh, itype), itype, atoms, rDerRho0MT ) ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe also more stable if it is stored todo analyze it later what is the better way
                tGaunt2 = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) = grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp * tGaunt2 &
                    &* (rDerRho0MT(imesh)  * clnu_atom(mem, lh, iatom) + ((oqn_l + 1) * rho0MT(imesh, lh, itype) * clnu_atom(mem, lh, iatom)&
                    &/ atoms%rmsh(imesh, itype))) !indices of rho
                enddo ! imesh
              endif
            end do ! mem
          end do ! lh
        end do ! ieqat
      end do ! itype
    end do ! mqn_mpp

   iatom = 0
   do itype = 1, atoms%ntype
     do ieqat = 1, atoms%neq(itype)
       iatom = iatom + 1
       do oqn_l = 0, atoms%lmax(itype)
         do mqn_m = -oqn_l, oqn_l
           lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
           do imesh = 1, atoms%jri(itype)
             grRho0MT(imesh, lm, iatom, 1:3) = matmul( Tmatrix(1:3, 1:3), grRho0MTNat(imesh, lm, iatom, 1:3) )

             if (compPhon) then
               write(109,*) imesh, lm, iatom, 1
               write(109,*) real(grRho0MT(imesh,lm,iatom,1)), aimag(grRho0MT(imesh,lm,iatom,1))
               write(109,*) imesh, lm, iatom, 2
               write(109,*) real(grRho0MT(imesh,lm,iatom,2)), aimag(grRho0MT(imesh,lm,iatom,2))
               write(109,*) imesh, lm, iatom, 3
               write(109,*) real(grRho0MT(imesh,lm,iatom,3)), aimag(grRho0MT(imesh,lm,iatom,3))
             end if

           enddo
         enddo ! mqn_m
       enddo ! oqn_l
     enddo ! ieqat
   enddo ! itype

  end subroutine calcGrFinLH

  !---------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> Calculates the spherical harmonic expansion coefficients of the muffin-tin gradient applied to an arbitrary function multiplied
  !> by $r^2$. The resulting gradient expansion coefficients are multiplied by a factor $r^2$.
  !>
  !> @note
  !> The ingoing function is assumed to be multiplied with $r$. The outgoing resulting function is also multiplied by $r$.
  !>
  !> @param[in]  atoms     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in]  lathar    : Contains entities concerning the lattice harmonics; more precise definition in type.F90 file.
  !> @param[in]  clnu_atom : Coefficients to transform from lattice harmonics to spherical harmonics without any symmetry.
  !> @param[in]  nmem_atom : Number of member per lattice harmonic in a system without any symmetry.
  !> @param[in]  mlh_atom  : Magnetic quantum numbers of the members of any lattice harmonic in a system without any symmetry.
  !> @param[in]  r2FlhMt   : Lattice harmonic coefficients of muffin-tin quantity multiplied by a factor of r**2.
  !> @param[out] r2GrFshMt : Spherical harmonic coefficients of muffin-tin quantity's gradient multiplied by a factor of r**2
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine calcGrR2FinLH(atoms, lathar, clnu_atom, nmem_atom, mlh_atom, r2FlhMt, r2GrFshMt)

    use m_JPConstants, only : Tmatrix, fpi, compPhon
    use m_gaunt, only : gaunt1
    use m_types

    implicit none

    ! Type parameters
    ! ***************
    type(t_atoms),               intent(in)  :: atoms
    type(t_sphhar),              intent(in)  :: lathar

    ! Array parameters
    ! ****************
    complex,                     intent(in)  :: clnu_atom(:, 0:, :)
    integer,                     intent(in)  :: nmem_atom(0:, :)
    integer,                     intent(in)  :: mlh_atom(:, 0:, :)
    real,                        intent(in)  :: r2FlhMt(:, 0:, :)
    complex,        allocatable, intent(out) :: r2GrFshMt(:, :, :, :)


    ! Local Scalar Variables
    ! **********************
    ! pfac    : Prefactor
    ! tGaunt  :  Gaunt coefficient
    ! itype   : Loop index for atom types
    ! ieqat   : Loop index for equivalent atoms
    ! iatom   : Loop index for all atoms
    ! imesh   : Loop index for radial mesh point
    ! mqn_m   : Magnetic quantum number m
    ! oqn_l   : Orbital quantum number l
    ! mqn_mpp : Magnetic quantum number double primed to index the natural coordinates
    ! lm      : Collective index for orbital and magnetic quantum number
    ! symType : Index of the symmetry
    ! ilh     : Loop index for different lattice harmonics (not their members!)
    ! imem    : Loop index for members of a lattice harmonics
    real                                     :: pfac
    real                                     :: tGaunt
    integer                                  :: itype
    integer                                  :: ieqat
    integer                                  :: iatom
    integer                                  :: imesh
    integer                                  :: mqn_m
    integer                                  :: oqn_l
    integer                                  :: mqn_mpp
    integer                                  :: lm
    integer                                  :: symType
    integer                                  :: ilh
    integer                                  :: imem

    ! Local Array Variables
    ! *********************
    ! rDerFlhMt    : Radial derrivative of the incoming fuction
    ! r2GrFshMtNat : Expansion coefficients of the muffin-tin gradient applied to the incoming function. The coefficients are given
    !                in natural coordinates and multiplied by $r^2$
    real,           allocatable              :: rDerFlhMt(:)
    complex,        allocatable              :: r2GrFshMtNat(:, :, :, :)


    ! Initialization of additionaly required arrays.
    allocate( r2GrFshMt(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3), &
            & r2GrFshMtNat(atoms%jmtd, ( atoms%lmaxd + 1)**2, atoms%nat, 3) )
    allocate( rDerFlhMt(atoms%jmtd) )
    r2GrFshMt = cmplx(0., 0.)
    r2GrFshMtNat = cmplx(0., 0.)
    rDerFlhMt = 0.

    pfac = sqrt( fpi / 3. )
    do mqn_mpp = -1, 1
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          symType = atoms%ntypsy(iatom)
          do ilh = 0, lathar%nlh(symType)
            oqn_l = lathar%llh(ilh, symType)
            do imem = 1, nmem_atom(ilh, iatom)
              mqn_m = mlh_atom(imem, ilh, iatom)

              ! l + 1 block
              ! oqn_l - 1 to l, so oqn_l should be < lmax not <= lmax
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < atoms%lmax(itype)) ) then
                lm = ( oqn_l + 1 ) * ( oqn_l + 2 ) + 1 + mqn_m - mqn_mpp
                call Derivative( r2FlhMt(:, ilh, itype), itype, atoms, rDerFlhMt )
                tGaunt = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  r2GrFshMtNat(imesh, lm, iatom, mqn_mpp + 2) = r2GrFshMtNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp &
                    &* tGaunt * (rDerFlhMt(imesh) * clnu_atom(imem, ilh, iatom) &
                    &- ((oqn_l + 2) * r2FlhMt(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom) / atoms%rmsh(imesh, itype)))
                end do ! imesh
              end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l )

              ! l - 1 block
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
                if ( oqn_l - 1 == -1 ) then
                  write (*, *) 'oqn_l too low'
                end if
                lm = (oqn_l - 1) * oqn_l + 1 + mqn_m - mqn_mpp
                ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe store it?
                call Derivative( r2FlhMt(:, ilh, itype), itype, atoms, rDerFlhMt )
                tGaunt = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  r2GrFshMtNat(imesh, lm, iatom, mqn_mpp + 2) = r2GrFshMtNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp &
                    & * tGaunt * (rDerFlhMt(imesh)  * clnu_atom(imem, ilh, iatom) &
                    & + ((oqn_l - 1) * r2FlhMt(imesh, ilh, itype) * clnu_atom(imem, ilh, iatom) / atoms%rmsh(imesh, itype)))
                end do ! imesh
              end if ! ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l )
            end do ! imem
          end do ! ilh
        end do ! ieqat
      end do ! itype
    end do ! mqn_mpp

    ! Conversion from natural to cartesian coordinates
    iatom = 0
    do itype = 1, atoms%ntype
      do ieqat = 1, atoms%neq(itype)
        iatom = iatom + 1
        do oqn_l = 0, atoms%lmax(itype)
          do mqn_m = -oqn_l, oqn_l
            lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
            do imesh = 1, atoms%jri(itype)
              r2GrFshMt(imesh, lm, iatom, 1:3) = matmul( Tmatrix(1:3, 1:3), r2GrFshMtNat(imesh, lm, iatom, 1:3) )
              if (compPhon) then
                write(109,*) imesh, lm, iatom, 1
                write(109,*) real(r2GrFshMt(imesh,lm,iatom,1))/atoms%rmsh(imesh, itype)**2, &
                           aimag(r2GrFshMt(imesh,lm,iatom,1))/atoms%rmsh(imesh, itype)**2
                write(109,*) imesh, lm, iatom, 2
                write(109,*) real(r2GrFshMt(imesh,lm,iatom,2))/atoms%rmsh(imesh, itype)**2, &
                           aimag(r2GrFshMt(imesh,lm,iatom,2))/atoms%rmsh(imesh, itype)**2
                write(109,*) imesh, lm, iatom, 3
                write(109,*) real(r2GrFshMt(imesh,lm,iatom,3))/atoms%rmsh(imesh, itype)**2, &
                            aimag(r2GrFshMt(imesh,lm,iatom,3))/atoms%rmsh(imesh, itype)**2
              end if
            end do
          end do ! mqn_m
        end do ! oqn_l
      end do ! ieqat
    end do ! itype

  end subroutine calcGrR2FinLH


  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine calculates the muffin-tin derivative according to formula in PhD thesis Aaron.
  !> @details
  !>
  !> @param[in] atomsT     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in] latharT
  !> @param[in] clnu_atom
  !> @param[in] nmem_atom
  !> @param[in] mlh_atom
  !> @param[in] rho0MT
  !> @param[in] gradrho0MT
  !--------------------------------------------------------------------------------------------------------------------------------------
  !TODO put this into utils module
  subroutine calcGrFLhNatAt(atoms, lathar, itype, iatom, clnu_atom, nmem_atom, mlh_atom, fLh, grFLhNat)

    use m_JPConstants, only : Tmatrix, fpi
    use m_gaunt, only : gaunt1
    use m_types

    implicit none

    ! Type parameters
    type(t_atoms),                     intent(in)  :: atoms
    type(t_sphhar),                    intent(in)  :: lathar

    ! Scalar parameters
    integer,                           intent(in)  :: itype
    integer,                           intent(in)  :: iatom

    ! Array parameters
    complex,                           intent(in)  :: clnu_atom(:, 0:, :)
    integer,                           intent(in)  :: nmem_atom(0:, :)
    integer,                           intent(in)  :: mlh_atom(:, 0:, :)
    real,                              intent(in)  :: fLh(:, 0:)
    complex,                           intent(out) :: grFLhNat(:, :, :)

    ! Local Variables:
    !
    ! pfac        : contains sqrt(4 pi / 3)
    ! lm          : encodes oqn_l and mqn_m
    ! gauntFactor : auxillary variable to store a Gaunt coefficient
    ! imesh       : runs over mesh points of current grid
    ! mqn_m       : magnetic quantum number m
    ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
    ! oqn_l       : orbital quantum number l

    ! Local Scalar Variables
    real                                           :: pfac
    real                                           :: gauntFactor
    integer                                        :: imesh
    integer                                        :: mqn_m
    integer                                        :: oqn_l
    integer                                        :: mqn_m2Pr
    integer                                        :: lm
    integer                                        :: lm_temp
    integer                                        :: symType
    integer                                        :: lh
    integer                                        :: mem

    ! Local Array Variables
    real,                  allocatable             :: delrFLh(:)

    allocate( delrFLh(atoms%jmtd))
    delrFLh = 0.

    pfac = sqrt( fpi / 3 )
    do mqn_m2Pr = -1, 1
      symType = atoms%ntypsy(iatom)
      do lh = 0, lathar%nlh(symType)
        oqn_l = lathar%llh(lh, symType)
        call Derivative( fLh(:, lh), itype, atoms, delrFLh )
        do mem = 1, nmem_atom(lh, iatom)
          mqn_m = mlh_atom(mem, lh, iatom)

          ! l + 1 scattering channel
          if ( ( abs(mqn_m - mqn_m2Pr) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < atoms%lmax(itype))) then
            lm = ( oqn_l + 1 ) * ( oqn_l + 2 ) + mqn_m - mqn_m2Pr
            gauntFactor = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmaxd )
            do imesh = 1, atoms%jri(itype)
              grFLhNat(imesh, lm, mqn_m2Pr + 2) = grFLhNat(imesh, lm, mqn_m2Pr + 2) + pfac * (-1)**mqn_m2Pr * gauntFactor &
                                           & * ( delrFLh(imesh) * clnu_atom(mem, lh, iatom) &
                                           & - ( oqn_l * fLh(imesh, lh) * clnu_atom(mem, lh, iatom) / atoms%rmsh(imesh, itype) ) )
            enddo ! imesh
          endif ! l + 1 scattering channel

          ! l - 1 scattering channel
          if ( ( abs(mqn_m - mqn_m2Pr) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
            lm = (oqn_l - 1) * oqn_l + mqn_m - mqn_m2Pr
            gauntFactor = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_m2Pr, mqn_m, -mqn_m2Pr, atoms%lmaxd )
            do imesh = 1, atoms%jri(itype)
              grFLhNat(imesh, lm, mqn_m2Pr + 2) = grFLhNat(imesh, lm, mqn_m2Pr + 2) + pfac * (-1)**mqn_m2Pr * gauntFactor &
                                     & * ( delrFLh(imesh)  * clnu_atom(mem, lh, iatom) &
                                     & + ( (oqn_l + 1) * fLh(imesh, lh) * clnu_atom(mem, lh, iatom) / atoms%rmsh(imesh, itype) ) )
            enddo ! imesh
          endif ! l - 1 scattering channel

        end do ! mem
      end do ! lh
    end do ! mqn_m2Pr

  end subroutine calcGrFLhNatAt

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine calculates the muffin-tin derivative according to formula in PhD thesis Aaron.
  !> @details
  !>
  !> @param[in] atomsT     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in] latharT
  !> @param[in] clnu_atom
  !> @param[in] nmem_atom
  !> @param[in] mlh_atom
  !> @param[in] rho0MT
  !> @param[in] gradrho0MT
  !--------------------------------------------------------------------------------------------------------------------------------------
  !TODO put this into utils module
  !TODO rename to CalcGradFuncMT
  subroutine calcGrFinSH(atoms, rho0MT, grRho0MT)

    use m_JPConstants, only : Tmatrix, fpi
    use m_gaunt, only : gaunt1
    use m_types

    implicit none

    type(t_atoms),                     intent(in)  :: atoms
    complex,                           intent(in)  :: rho0MT(:, :, :)
    complex,        allocatable,       intent(out) :: grRho0MT(:, :, :, :)

    ! Local Variables:
    !
    ! sqr4pi3     : contains sqrt(4 pi / 3)
    ! iatom       : runs over all atoms
    ! lm          : encodes oqn_l and mqn_m
    ! itype       : runs over all atom types
    ! idirec      : runs over 3 directions the atom can be displaced to
    ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
    ! tempGaunt2  : auxillary variable to store a Gaunt coefficient
    ! imesh       : runs over mesh points of current grid
    ! rl2         : stores R^(l + 2)
    ! fint        : stores integral
    ! ll1         : auxillary variable to calculate lm
    ! mqn_m       : magnetic quantum number m
    ! sk3r        : stores argument of spherical Bessel function
    ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
    ! oqn_l       : orbital quantum number l
    ! iGvec       : indexes current G-vector
    ! grRho0MT  : stores the gradient of the unperturbed density in the MTs (equation 7.59, PhD thesis Aaron Klüppelberg)
    ! gradrho0PWi : stores the the second part of equation 7.58 using equation 7.60 from PhD thesis Aaron Klüppelberg
    ! pylm        : contains 4π i i^l G / |G| exp(i G τ)  Y*_lm(G / |G|)
    ! sbes        : stores the Bessel function
    ! cil         : stores everything except for pylm of the result of this routine
    ! f           : stores the integrand of the first term in (7.58, PhD thesis Aaron Klüppelberg)

    ! Local Scalar Variables
    real                                           :: pfac
    real                                           :: tGaunt1, tGaunt2
    integer                                        :: itype
    integer                                        :: ieqat
    integer                                        :: iatom
    integer                                        :: imesh
    integer                                        :: mqn_m
    integer                                        :: oqn_l
    integer                                        :: mqn_mpp
    integer                                        :: lm
    integer                                        :: lmRho
    integer                                        :: lm_temp
    integer                                        :: lmMinus
    integer                                        :: lmPlus

    ! Local Array Variables
    real,       allocatable                         :: rDerRho0MTre(:)
    real,       allocatable                         :: rDerRho0MTim(:)
    complex,    allocatable                         :: rDerRho0MT(:)
    complex,    allocatable                         :: grRho0MTNat(:, :, :, :)


    !todo is r^2 rho numerically more stable?
    allocate( grRho0MT(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    allocate( grRho0MTNat(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    allocate( rDerRho0MTre(atoms%jmtd), rDerRho0MTim(atoms%jmtd), rDerRho0MT(atoms%jmtd) )
    rDerRho0MTre(:) = 0.
    rDerRho0MTim(:) = 0.
    rDerRho0MT(:) = 0.

    grRho0MTNat = cmplx(0., 0.)
    grRho0MT = cmplx(0., 0.)

    pfac = sqrt( fpi / 3 )
    do mqn_mpp = -1, 1
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
            ! todo second if clause not required!
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < atoms%lmax(itype) ) ) then
                lm = ( oqn_l + 1 ) * ( oqn_l + 2 ) + 1 + mqn_m - mqn_mpp
                lmRho = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
                rDerRho0MTre(:) = 0.
                rDerRho0MTim(:) = 0.
                rDerRho0MT(:) = 0.
                call Derivative( real(rho0MT(:, lmRho, itype)), itype, atoms, rDerRho0MTre )
                call Derivative( aimag(rho0MT(:, lmRho, itype)), itype, atoms, rDerRho0MTim )
                do imesh = 1, atoms%jri(itype)
                  rDerRho0MT(imesh) = cmplx(rDerRho0MTre(imesh), rDerRho0MTim(imesh))
                end do ! imesh
                tGaunt1 = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) = grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp * tGaunt1 * (rDerRho0MT(imesh) &
                                                               &- (oqn_l * rho0MT(imesh, lmRho, iatom) / atoms%rmsh(imesh, itype)))
                enddo ! imesh
              endif
              !todo second clause might be for security reasons if ms of lattice harmoincs are wrongly stored
            ! todo second if clause not required!
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
                if ( oqn_l - 1 == -1 ) then
                  write (*, *) 'oqn_l too low'
                end if
                lmRho = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
                lm = (oqn_l - 1) * oqn_l + 1 + mqn_m - mqn_mpp
                rDerRho0MTre(:) = 0.
                rDerRho0MTim(:) = 0.
                rDerRho0MT(:) = 0.
                call Derivative( real(rho0MT(:, lmRho, itype)), itype, atoms, rDerRho0MTre ) ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe also more stable if it is stored todo analyze it later what is the better way
                call Derivative( aimag(rho0MT(:, lmRho, itype)), itype, atoms, rDerRho0MTim ) ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe also more stable if it is stored todo analyze it later what is the better way
                do imesh = 1, atoms%jri(1)
                  rDerRho0MT(imesh) = cmplx(rDerRho0MTre(imesh), rDerRho0MTim(imesh))
                end do ! imesh
                tGaunt2 = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) = grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp * tGaunt2 &
                    &* (rDerRho0MT(imesh) + ((oqn_l + 1) * rho0MT(imesh, lmRho, iatom)&
                    &/ atoms%rmsh(imesh, itype))) !indices of rho
                enddo ! imesh
              endif
            end do ! mem
          end do ! lh
        end do ! ieqat
      end do ! itype
    end do ! mqn_mpp

   iatom = 0
   do itype = 1, atoms%ntype
     do ieqat = 1, atoms%neq(itype)
       iatom = iatom + 1
       do oqn_l = 0, atoms%lmax(itype)
         do mqn_m = -oqn_l, oqn_l
           lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
           do imesh = 1, atoms%jri(itype)
             grRho0MT(imesh, lm, iatom, :) = matmul( Tmatrix, grRho0MTNat(imesh, lm, iatom, :) )
           enddo
         enddo ! mqn_m
       enddo ! oqn_l
     enddo ! ieqat
   enddo ! itype

  end subroutine calcGrFinSH

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine calculates the muffin-tin derivative according to formula in PhD thesis Aaron.
  !> @details
  !>
  !> @param[in] atomsT     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in] latharT
  !> @param[in] clnu_atom
  !> @param[in] nmem_atom
  !> @param[in] mlh_atom
  !> @param[in] rho0MT
  !> @param[in] gradrho0MT
  !--------------------------------------------------------------------------------------------------------------------------------------
  !TODO put this into utils module
  !TODO rename to CalcGradFuncMT
  subroutine calcGrR2FinSH(atoms, rho0MT, grRho0MT)

    use m_JPConstants, only : Tmatrix, fpi
    use m_gaunt, only : gaunt1
    use m_types

    implicit none

    type(t_atoms),                     intent(in)  :: atoms
    complex,                           intent(in)  :: rho0MT(:, :, :)
    complex,                           intent(out) :: grRho0MT(:, :, :, :)

    ! Local Variables:
    !
    ! sqr4pi3     : contains sqrt(4 pi / 3)
    ! iatom       : runs over all atoms
    ! lm          : encodes oqn_l and mqn_m
    ! itype       : runs over all atom types
    ! idirec      : runs over 3 directions the atom can be displaced to
    ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
    ! tempGaunt2  : auxillary variable to store a Gaunt coefficient
    ! imesh       : runs over mesh points of current grid
    ! rl2         : stores R^(l + 2)
    ! fint        : stores integral
    ! ll1         : auxillary variable to calculate lm
    ! mqn_m       : magnetic quantum number m
    ! sk3r        : stores argument of spherical Bessel function
    ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
    ! oqn_l       : orbital quantum number l
    ! iGvec       : indexes current G-vector
    ! grRho0MT  : stores the gradient of the unperturbed density in the MTs (equation 7.59, PhD thesis Aaron Klüppelberg)
    ! gradrho0PWi : stores the the second part of equation 7.58 using equation 7.60 from PhD thesis Aaron Klüppelberg
    ! pylm        : contains 4π i i^l G / |G| exp(i G τ)  Y*_lm(G / |G|)
    ! sbes        : stores the Bessel function
    ! cil         : stores everything except for pylm of the result of this routine
    ! f           : stores the integrand of the first term in (7.58, PhD thesis Aaron Klüppelberg)

    ! Local Scalar Variables
    real                                           :: pfac
    real                                           :: tGaunt1, tGaunt2
    integer                                        :: itype
    integer                                        :: ieqat
    integer                                        :: iatom
    integer                                        :: imesh
    integer                                        :: mqn_m
    integer                                        :: oqn_l
    integer                                        :: mqn_mpp
    integer                                        :: lm
    integer                                        :: lmRho
    integer                                        :: lm_temp
    integer                                        :: lmMinus
    integer                                        :: lmPlus

    ! Local Array Variables
    real,       allocatable                         :: rDerRho0MTre(:)
    real,       allocatable                         :: rDerRho0MTim(:)
    complex,    allocatable                         :: rDerRho0MT(:)
    complex,    allocatable                         :: grRho0MTNat(:, :, :, :)


    !todo is r^2 rho numerically more stable?
    allocate( grRho0MTNat(atoms%jmtd, ( atoms%lmaxd + 1 )**2, atoms%nat, 3) )
    allocate( rDerRho0MTre(atoms%jmtd), rDerRho0MTim(atoms%jmtd), rDerRho0MT(atoms%jmtd) )
    rDerRho0MTre(:) = 0.
    rDerRho0MTim(:) = 0.
    rDerRho0MT(:) = 0.
    grRho0MTNat(:, :, :, :) = cmplx(0., 0.)

    grRho0MT = cmplx(0., 0.)

    pfac = sqrt( fpi / 3 )
    do mqn_mpp = -1, 1
      iatom = 0
      do itype = 1, atoms%ntype
        do ieqat = 1, atoms%neq(itype)
          iatom = iatom + 1
          do oqn_l = 0, atoms%lmax(itype)
            do mqn_m = -oqn_l, oqn_l
            ! todo second if clause not required!
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < atoms%lmax(itype)) ) then
                lm = ( oqn_l + 1 ) * ( oqn_l + 2 ) + 1 + mqn_m - mqn_mpp
                lmRho = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
                rDerRho0MTre(:) = 0.
                rDerRho0MTim(:) = 0.
                call Derivative( real(rho0MT(:, lmRho, itype)), itype, atoms, rDerRho0MTre )
                call Derivative( aimag(rho0MT(:, lmRho, itype)), itype, atoms, rDerRho0MTim )
                do imesh = 1, atoms%jri(itype)
                  rDerRho0MT(imesh) = cmplx(rDerRho0MTre(imesh), rDerRho0MTim(imesh))
                end do ! imesh
                tGaunt1 = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) = grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp * tGaunt1 * (rDerRho0MT(imesh) &
                                                               &- ((oqn_l + 2)* rho0MT(imesh, lmRho, iatom) / atoms%rmsh(imesh, itype)))
                enddo ! imesh
              endif
              !todo second clause might be for security reasons if ms of lattice harmoincs are wrongly stored
            ! todo second if clause not required!
              if ( ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
                if ( oqn_l - 1 == -1 ) then
                  write (*, *) 'oqn_l too low'
                end if
                lmRho = oqn_l * ( oqn_l + 1 ) + 1 + mqn_m
                lm = (oqn_l - 1) * oqn_l + 1 + mqn_m - mqn_mpp
                rDerRho0MTre(:) = 0.
                rDerRho0MTim(:) = 0.
                call Derivative( real(rho0MT(:, lmRho, itype)), itype, atoms, rDerRho0MTre ) ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe also more stable if it is stored todo analyze it later what is the better way
                call Derivative( aimag(rho0MT(:, lmRho, itype)), itype, atoms, rDerRho0MTim ) ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe also more stable if it is stored todo analyze it later what is the better way
                do imesh = 1, atoms%jri(1)
                  rDerRho0MT(imesh) = cmplx(rDerRho0MTre(imesh), rDerRho0MTim(imesh))
                end do ! imesh
                tGaunt2 = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, atoms%lmaxd )
                do imesh = 1, atoms%jri(itype)
                  grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) = grRho0MTNat(imesh, lm, iatom, mqn_mpp + 2) + pfac * (-1)**mqn_mpp * tGaunt2 &
                    &* (rDerRho0MT(imesh) + ((oqn_l - 1) * rho0MT(imesh, lmRho, iatom)&
                    &/ atoms%rmsh(imesh, itype))) !indices of rho
                enddo ! imesh
              endif
            end do ! mem
          end do ! lh
        end do ! ieqat
      end do ! itype
    end do ! mqn_mpp

   iatom = 0
   do itype = 1, atoms%ntype
     do ieqat = 1, atoms%neq(itype)
       iatom = iatom + 1
       do oqn_l = 0, atoms%lmax(itype)! + 1
         do mqn_m = -oqn_l, oqn_l
           lm = oqn_l * (oqn_l + 1) + 1 + mqn_m
           do imesh = 1, atoms%jri(itype)
             grRho0MT(imesh, lm, iatom, 1:3) = matmul( Tmatrix(1:3, 1:3), grRho0MTNat(imesh, lm, iatom, 1:3) )
           enddo
         enddo ! mqn_m
       enddo ! oqn_l
     enddo ! ieqat
   enddo ! itype
   write(*, *) 'was once here'

  end subroutine calcGrr2FinSH

  !--------------------------------------------------------------------------------------------------------------------------------------
  !> @author
  !> Christian-Roman Gerhorst, Forschungszentrum Jülich: IAS1 / PGI1
  !>
  !> @brief
  !> This subroutine calculates the muffin-tin derivative according to formula in PhD thesis Aaron.
  !> @details
  !>
  !> @param[in] atomsT     : Contains atoms-related quantities; definition of its members in types.F90 file.
  !> @param[in] latharT
  !> @param[in] clnu_atom
  !> @param[in] nmem_atom
  !> @param[in] mlh_atom
  !> @param[in] rho0MT
  !> @param[in] gradrho0MT
  !--------------------------------------------------------------------------------------------------------------------------------------
  !deprecated
  subroutine grFlmpYlmPerType( atoms, itype, lmaxIn, nRadFun, fLm, grFlm )

    use m_JPConstants, only : Tmatrix, fpi
    use m_gaunt, only : gaunt1
    use m_types

    implicit none

    ! Type parameter
    type(t_atoms),                     intent(in)  :: atoms

    ! Scalar parameter
    integer,                           intent(in)  :: itype
    integer,                           intent(in)  :: lmaxIn

    ! Array parameter
    integer,                           intent(in)  :: nRadFun(0:, :)
    real,                              intent(in)  :: fLm(:, :, 0:)
    complex,                           intent(out) :: grFlm(:, :, :, 0:)

    ! Local Variables:
    !
    ! sqr4pi3     : contains sqrt(4 pi / 3)
    ! lm          : encodes oqn_l and mqn_m
    ! itype       : runs over all atom types
    ! idirec      : runs over 3 directions the atom can be displaced to
    ! tempGaunt1  : auxillary variable to store a Gaunt coefficient
    ! tempGaunt2  : auxillary variable to store a Gaunt coefficient
    ! imesh       : runs over mesh points of current grid
    ! rl2         : stores R^(l + 2)
    ! fint        : stores integral
    ! ll1         : auxillary variable to calculate lm
    ! mqn_m       : magnetic quantum number m
    ! sk3r        : stores argument of spherical Bessel function
    ! mqn_mpp     : magnetic quantum number m", also used for indexing 3 directions the atom can be displaced to
    ! oqn_l       : orbital quantum number l
    ! iGvec       : indexes current G-vector
    ! grFlm  : stores the gradient of the unperturbed density in the MTs (equation 7.59, PhD thesis Aaron Klüppelberg)
    ! gradrho0PWi : stores the the second part of equation 7.58 using equation 7.60 from PhD thesis Aaron Klüppelberg
    ! pylm        : contains 4π i i^l G / |G| exp(i G τ)  Y*_lm(G / |G|)
    ! sbes        : stores the Bessel function
    ! cil         : stores everything except for pylm of the result of this routine
    ! f           : stores the integrand of the first term in (7.58, PhD thesis Aaron Klüppelberg)

    ! Local Scalar Variables
    real                                           :: pfac
    real                                           :: tGaunt1, tGaunt2
    integer                                        :: imesh
    integer                                        :: mqn_m
    integer                                        :: oqn_l
    integer                                        :: mqn_mpp
    integer                                        :: lmIn
    integer                                        :: lmOut
    integer                                        :: lm_temp
    integer                                        :: lmMinus
    integer                                        :: lmPlus
    integer                                        :: p

    ! Local Array Variables
    real                                           :: rDerRho0MT(atoms%jri(itype))

    grFlm = cmplx(0., 0.)

    pfac = sqrt( fpi / 3 )
    do oqn_l = 0, lmaxIn
      do mqn_m = - oqn_l, oqn_l
        lmIn = oqn_l * (oqn_l + 1) + mqn_m
        do mqn_mpp = -1, 1
          if ( ( abs(mqn_m - mqn_mpp) <= oqn_l + 1 ) .and. ( abs(mqn_m) <= oqn_l ) .and. (oqn_l < lmaxIn) ) then
            lmOut = ( oqn_l + 1 ) * ( oqn_l + 2 ) + mqn_m - mqn_mpp
            do p = 1, nRadFun(oqn_l, itype)
              call Derivative( fLm(:, p, lmIn), itype, atoms, rDerRho0MT )
              tGaunt1 = Gaunt1( oqn_l + 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, lmaxIn )
              do imesh = 1, atoms%jri(itype)
                grFlm(imesh, p, mqn_mpp + 2, lmOut) = grFlm(imesh, p, mqn_mpp + 2, lmOut) + pfac * (-1)**mqn_mpp * &
                  & tGaunt1 * (rDerRho0MT(imesh) - (oqn_l * fLm(imesh, p, lmIn) / atoms%rmsh(imesh, itype)))
              end do ! imesh
            end do ! p
          endif
          if ( ( abs(mqn_m - mqn_mpp) <= oqn_l - 1 ) .and. ( abs(mqn_m) <= oqn_l ) ) then
            if ( oqn_l - 1 == -1 ) then
              write (*, *) 'oqn_l too low'
            end if
            lmOut = (oqn_l - 1) * oqn_l + mqn_m - mqn_mpp
            do p = 1, nRadFun(oqn_l, itype)
              ! This is also a trade of between storage and performance, because derivative is called redundantly, maybe also more
              ! stable if it is stored todo analyze it later what is the better way
              call Derivative( fLm(:, p, lmIn), itype, atoms, rDerRho0MT )
              tGaunt2 = Gaunt1( oqn_l - 1, oqn_l, 1, mqn_m - mqn_mpp, mqn_m, -mqn_mpp, lmaxIn )
              do imesh = 1, atoms%jri(itype)
                ! indices of rho
                grFlm(imesh, p, mqn_mpp + 2, lmOut) = grFlm(imesh, p, mqn_mpp + 2, lmOut) + pfac * (-1)**mqn_mpp * &
                  & tGaunt2 * (rDerRho0MT(imesh) + ((oqn_l + 1) * fLm(imesh, p, lmIn) / atoms%rmsh(imesh, itype)))
              end do ! imesh
            end do ! p
          end if
        end do ! mqn_mpp
      end do ! mqn_m
    end do ! oqn_l

    do oqn_l = 0, lmaxIn! + 1
      do mqn_m = -oqn_l, oqn_l
        lmOut = oqn_l * (oqn_l + 1) + mqn_m
        do p = 1, nRadFun(oqn_l, itype)
          do imesh = 1, atoms%jri(itype)
            grFlm(imesh, p, :, lmOut) = matmul( Tmatrix, grFlm(imesh, p, :, lmOut) )
          end do
        end do
      end do ! mqn_m
    end do ! oqn_l

  end subroutine grFlmpYlmPerType

  ! Idea taken from Numerical recipes in FORTRAN 90: the art of parallel scientific computing, Cambridge University Press, 1996, s.970, chapter 22 between equations 22.1.5 and 22.1.6
  ! This function takes 2 real vectors. None of the parameters should be transposed. This is taken care of within the routine.
  ! this routine has been tested by comparing results to analytically calculated results
  function outerProduct(a, b)

    implicit none

    real, intent(in) :: a(:)
    real, intent(in) :: b(:)

    real             :: outerProduct(size(a), size(b))

    outerProduct(:, :) = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))

  end function outerProduct

  ! This functions calculates the complex matrix element ij of the outer product of a and b.
  !DIR$ ATTRIBUTES INLINE :: outerProductME
  function outerProductME(a, b, i, j)

    use m_juDFT_NOstopNO, only : juDFT_error

    implicit none

    complex,    intent(in)  :: a(:)
    complex,    intent(in)  :: b(:)
    integer,    intent(in)  :: i
    integer,    intent(in)  :: j

    complex                 :: outerProductME

    if (i > size(a) .or. j > size(b)) then
      call juDFT_error( 'Wished matrix element is out of scope of outer product matrix.', calledby='outerProductME', &
        & hint='Choose smaller indices.')
    end if

    outerProductME = a(i) * b(j)

  end function outerProductME

  !Note: if we perform an FFT and FFT back we get cut-off errors. Therefore we read in the real space xc stuff and do only one fourier transformation
  !this routine at the moment only works for real input, which is sufficient for the xc potentials.
  ! This routine seems not be used any more, but we still leave it here, because it might be somewhere called from dead in the test
  ! routines that are not reviewed yet
  subroutine calcFFTreal2reciprocal(ngpqdp, nfft, gfft, gpqdp, w_gfft)

    use m_cfft
    implicit none


    ! Scalar parameters
    integer,                                      intent(in)  :: ngpqdp

    integer,                                      intent(in)  :: nfft(:)
    real,                            intent(in)              :: gfft(0:, :)
    integer,                                      intent(in)  :: gpqdp(:, :)
    complex,                                      intent(out) :: w_gfft(:)

    ! Scalar variables
    integer                                                   :: iG
    integer                                                   :: idir
    integer                                                   :: ifftd
    real                                                      :: scaling

    ! Array variables
    integer,                         allocatable              :: igfft(:)
    integer                                                   :: gabs(3)
    real,                            allocatable              :: heregfft(:, :)

    ifftd = product(nfft)
    allocate(heregfft(0:ifftd-1, 2))

    heregfft = 0.
    heregfft(0:ifftd-1, 1) = gfft(0:ifftd-1, 1)


    allocate( igfft(ngpqdp) )
    gabs = 0
    igfft = 0
    do iG = 1, ngpqdp
      do idir = 1, 3
        if ( gpqdp(idir, iG) >= 0 ) then
          gabs(idir) = gpqdp(idir, iG)
        else
          gabs(idir) = gpqdp(idir, iG) + nfft(idir)
        end if
      end do
      igfft(iG) = gabs(1) + gabs(2) * nfft(1) + gabs(3) * nfft(1) * nfft(2)
    end do


    call cfft(heregfft(0, 1), heregfft(0, 2), ifftd, nfft(1), nfft(1), -1)
    call cfft(heregfft(0, 1), heregfft(0, 2), ifftd, nfft(2), nfft(1) * nfft(2), -1)
    call cfft(heregfft(0, 1), heregfft(0, 2), ifftd, nfft(3), ifftd, -1)

    w_gfft=0
    scaling = 1. / real(ifftd)
    do iG=1, ngpqdp
      w_gfft(iG) = w_gfft(iG) +  cmplx( heregfft(igfft(iG), 1), heregfft(igfft(iG), 2) ) * scaling
    end do


  end subroutine calcFFTreal2reciprocal

end module mod_juPhonUtils
