      MODULE m_efield
      USE m_juDFT
      USE m_constants
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: e_field
      CONTAINS
      SUBROUTINE e_field(atoms,  stars, sym, vacuum, cell, input,efield)
!
!*********************************************************************
!     sets the values of the sheets of charge for external electric
!     fields by requiring charge neutrality.
!     the position of the sheets of charge relative to the vacuum
!     boundary is fixed (10 a.u.), but can be set to a different
!     value in the file apwefl.
!
!     modified and fixed 10-99
!*********************************************************************

      USE m_types
      !USE m_setcor, ONLY: setcor
      IMPLICIT NONE
!     ..
!     .. Scalar Arguments ..
      TYPE(t_atoms), INTENT (IN)    :: atoms
      Type(t_stars),INTENT(IN)      :: stars
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_vacuum),INTENT(IN)     :: vacuum
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_input),INTENT(IN)      :: input
      TYPE(t_efield),INTENT(INOUT)  :: efield
!     ..
!     ..
!     .. Local Scalars ..
      REAL, parameter :: eps = 1.0e-7
      REAL  qn,qe,bmu
      INTEGER n,iwd,nst,nc
!     ..
!     ..
!     .. Local Parameters ..
      INTEGER, PARAMETER :: pTOP = 1, pBOT = 2, pTOPBOT = 3
      INTEGER, PARAMETER :: pADD = 1, pREPLACE = 2
      INTEGER, PARAMETER :: pALL = 1, pZERO = 2, pNONZERO = 3
!
!--->    obtain total nuclear charge
!
      qn=0.0
      DO n=1,atoms%ntype
         qn=qn+atoms%neq(n)*atoms%zatom(n)
      ENDDO
!
!--->    obtain total electronic charge (in electron units)
!
      qe=0.0
!---> core electrons
      DO n = 1,atoms%ntype
         IF (atoms%zatom(n).GE.1.0) THEN
            qe=qe+atoms%neq(n)*atoms%econf(n)%core_electrons
            WRITE (6,*) 'neq= ',atoms%neq(n),'  ncore= ',qe
         ENDIF
      ENDDO
!---> semi-core and valence electrons
      qe=qe+input%zelec
      WRITE (6,"(A, F12.8)") 'zelec=  ',input%zelec

      WRITE (6, '(/,/,a)') ' parameters for external electric field:'
      WRITE (6, '(3x,a,f12.5)') 'total electronic charge   =', qe
      WRITE (6, '(3x,a,f12.5)') 'total nuclear charge      =', qn

      CALL read_efield (efield, stars%mx1, stars%mx2, vacuum%nvac, cell%area)

      ! Sign convention of electric field: E>0 repels electrons,
      ! consistent with conventional choice that current is
      ! opposite to electron flow. With this choice, qe < qn
      ! corresponds to E>0
      ! In case of Dirichlet boundary conditions, we ignore the
      ! value of "sigma" and take the surplus charge into account
      ! in vvac.
      if (efield%autocomp .or. efield%dirichlet) efield%sigma = 0.5*(qe-qn)

      CALL print_efield (6, efield, cell%area, vacuum%nvac, cell%amat,vacuum%dvac)

      IF (.NOT. efield%dirichlet&
     &    .AND. ABS (SUM (efield%sig_b) + 2*efield%sigma - (qe-qn)) > eps) THEN
        IF (ABS (SUM (efield%sig_b) - (qe-qn)) < eps) THEN
          CALL juDFT_error&
     &          ("E field: top+bottom+film charge does not add up to "&
     &           //"zero.",&
     &           hint="Consider setting 'autocomp = false' in apwefl. "&
     &           //"(By default, the number of electrons in the energy "&
     &           //"window is automatically compensated via the charge "&
     &           //"sheets.)",&
     &           calledby ="efield")
        END IF
        CALL juDFT_error&
     &        ("E field: top+bottom+film charge does not add up to zero"&
     &        ,calledby ="efield")
      ENDIF

      IF (ABS (efield%sigma) > 0.49 .OR. ANY (ABS (efield%sig_b) > 0.49)) THEN
        WRITE ( 6,*) 'If you really want to calculate an e-field this'
        WRITE ( 6,*) 'big, uncomment STOP in efield.f !'
        CALL juDFT_error("E-field too big or No. of e- not correct"&
     &       ,calledby ="efield")
      ENDIF

      IF (efield%l_segmented) THEN
        CALL V_seg_EF(&
     &                efield,&
     &                vacuum, stars)

        IF (efield%plot_rho)&
     &    CALL print_rhoEF(&
     &                     efield, stars%mx1, stars%mx2, vacuum%nvac, stars%ng2, sym%nop, sym%nop2,&
     &                     stars%ng2, stars%kv2, sym%mrot, sym%symor, sym%tau, sym%invtab,&
     &                     cell%area, stars%nstr2, cell%amat)
      END IF

      IF (ALLOCATED (efield%sigEF)) DEALLOCATE (efield%sigEF)

      IF (efield%dirichlet .AND. ALLOCATED (efield%rhoEF))&
     &  call set_dirchlet_coeff (efield, vacuum%dvac, stars%ng2, stars%sk2, vacuum%nvac)

      CONTAINS

      SUBROUTINE set_dirchlet_coeff (E, dvac, nq2, sk2, nvac)
        TYPE(t_efield), INTENT(INOUT), TARGET :: E
        REAL, INTENT(IN) :: dvac
        INTEGER, INTENT(IN) :: nq2, nvac
        REAL, INTENT(IN) :: sk2(:)

        COMPLEX, POINTER :: V(:,:)
        REAL :: x, y, g
        INTEGER :: ig, ivac1, ivac2

        ALLOCATE (E%C1(nq2-1), E%C2(nq2-1))
        E%l_dirichlet_coeff = .true.

        V => E%rhoEF
        ivac1 = 1
        ivac2 = nvac

        DO ig = 1, nq2-1
          g = sk2(ig+1)
          x = EXP ( g*(E%zsigma+Dvac/2.0))
          y = EXP (-g*(E%zsigma+Dvac/2.0))
          E%C1(ig) = (V(ig, ivac1)*x - V(ig, ivac2)*y) / (x**2-y**2)
          E%C2(ig) = (V(ig, ivac2)*x - V(ig, ivac1)*y) / (x**2-y**2)
        END DO
      END SUBROUTINE set_dirchlet_coeff

! Read the electric field parameters from the 'apwefl' file
! (except for sigma, which is determined by the charge in the film).
!
! The position of the charged sheet can be given (in a.u.) via
! zsigma, an additional, in x-y homogeneous field can be given via
! sig_b(1) and sig_b(2). NOTE: The position zsigma counts from the
! interstitial-vacuum boundary ("z1") and not from the center (z=0).
!
! There are two input formats. The old format is
!    zsigma
!    sig_b(nvac)
! where sig_b is optional. zsigma is 10.0 by default while b_sig is zero.
!
! In the new format, lines starting with "!" or "#" are treated as comments
! and all keywords are case insensitive.
!
! Namelist (optionally, defaults to the shown values):
!   &efield zsigma=10.0, sig_b=0.0,0.0, plot_charge=f, plot_rho=f,
!           autocomp=t, dirichlet=f, eV=f /
!
! By default, Neumann boundary conditions are assumed and additional charge
! is placed as surface charge density at the charge sheets located at zsigma.
! If Dirichlet is true, a (segmented) metal plate at zsigma is assumed,
! which has a charge of zero - or the value given by sig_b or via adding
! different shapes. In case of Dirichlet boundary conditions, the value is
! regarded as potential in Hartree (or [electron] Volt, if "eV" is set to
! true.)
!
! A fancy potential can now be created by placing charge in a rectangular
! or circular region via the "rect" and "circ" directive. Both directives
! are followed by "top, "bot" or "TopBot"/"BotTop" and then a point (x,y),
! which denotes for "rect" the bottom-left coordinate and for "circ" the
! centre. Next parameter is the width/height for "rect" and the diameter
! radius for "circ". Followed by a parameter for the charge (or the
! potential in case of Dirichlet boundary conditions). Optional
! arguments: The charge can be added ("add", default) to the previous
! charge in this sheet - or it replaces ("replace") the charge. The
! charge can be placed in the whole area ("all", default) or only where
! previously no charge ("zero") or a "nonzero" charge was.
! The sig_b charge is added after all the shapes and is thus treated like
! a "rect 0,0 top/bot  1.0,1.0 b_sig(1)/b_sig(2)" would be.
!
! Note: All coordinates are relative coordinates. The regions can exceed the
! x-y extend of the film; e.g. using
!   circ top 0,0  0.25  0.5
! one places 1/2 an electron in a quarter circle with origin (0,0).
!
! WARNING: "circ" creates a perfect circle on the grid, however, it only
! matches a circle and not an ellipse if the k1d/k2d ratio matches the
! crystal a/b ratio.
! TODO: Support input in absolute values instead of relative ones.
! (Could be an option to the namelist or another (optional) tag.)
!
! rectSinX: Create a sinusoidal potential in x direction (constant in y
! direction for any x value), i.e.
!   V(x,y) = A * sin(2*pi*n*x + 2*pi*delta)
! A is the amplitude; however, the argument in apwefl is not A directly
! but   A' = A * Lz, where Lz is the number of points in z direction.
! Contrary to "circ" and "rect", charges are mask out without being
! redistributed to non-masked positions. Note:
!  Int_0^(2pi) A*|sin(x)| dx = 4*A)
! n is the order and delta the offset. "x" is relative to the rectangle
! and thus between 0 and 1. The syntax is:
!   rectSinX top/bot  x,y  w,h,  A', n, delta  [options]
!
!   datafile top/bot filename [zero_avg] [option]
! Read data points from "filename"; if zero_avg is present, average the
! read data to zero; replace/add is supported, but zero/not_zero is not.
!
! polygon: Create a polygon-shaped charge distribution; note that the
! currently used algorithm does not always give the perfactly shaped
! polygon - and the edge points are not always included in the polygon.
! Syntax:
!   polygon top/bot n_points, x1,y1, ..., xn,yn charge [options]
!
!
! Example 1: To have two top plates (segments):
!   rect top 0,  0  0.5,1.0  0.2
!   rect top 0.5,0  0.5,1.0 -0.2
!
! Example 2: To have a charged ring with 0.2e and -0.2e of charge
! evenly distributed outside this ring:
!   circ topBot 0.5,0.50.2  1           ! Create temporary an inner ring
!   circ topBot 0.5,0.50.3  0.2 zero    ! Create outer ring
!   circ topBot 0.5,0.50.2  0   replace ! Delete inner ring
!   rect topBot 0,0     1,1 -0.2 zero    ! Place smeared opposite charge
!
!
! If autocomp is .true., the extra charge in the film (Window) is distributed
! over both external sheets; if it is .false. one needs to compensate it
! manually.
!
! Note: The namelist has to be placed before the shapes when Dirichlet
! boundary conditions are used as otherwise the potential is divided by
! the number of grid points...

      SUBROUTINE read_efield (E, k1d, k2d, nvac, area)
        TYPE(t_efield), INTENT(INOUT) :: E
        INTEGER, INTENT(IN) :: k1d, k2d, nvac
        REAL, INTENT(IN) :: area

        REAL ::   tmp
        INTEGER :: i

        ! New format
        ALLOCATE(E%sigEF(3*k1d, 3*k2d, nvac))
        E%sigEF = 0.0
        if (allocated(e%shapes)) then
        DO i=1,SIZE(e%shapes)
           CALL read_shape (E, e%shapes(i), nvac)
        END DO
        endif
        IF (e%l_eV) THEN
          E%sig_b(:) = E%sig_b/hartree_to_ev_const
          E%sigEF(:,:,:) = E%sigEF/hartree_to_ev_const
        END IF
        ! Save average sigEF potential in sig_b and remove it from
        ! sigEF to avoid double counting; i.e. make g_|| = 0 of sigEF == 0.
        IF (E%dirichlet) THEN
          tmp = sum(E%sigEF(:,:,1))/SIZE(E%sigEF)*nvac
          E%sig_b(1) = E%sig_b(1) + tmp
          E%sigEF(:,:,1) = E%sigEF(:,:,1) - tmp
        ELSE
          tmp = sum(E%sigEF(:,:,1))
          E%sig_b(1) = E%sig_b(1) + tmp
          E%sigEF(:,:,1) = E%sigEF(:,:,1) - tmp/SIZE(E%sigEF)*nvac
        END IF
        IF (nvac > 1) THEN
          IF (E%dirichlet) THEN
            tmp = sum(E%sigEF(:,:,2))/SIZE(E%sigEF)*nvac
            E%sig_b(2) = E%sig_b(2) + tmp
            E%sigEF(:,:,2) = E%sigEF(:,:,2) - tmp
          ELSE
            tmp = sum(E%sigEF(:,:,2))
            E%sig_b(2) = E%sig_b(2) + tmp
            E%sigEF(:,:,2) = E%sigEF(:,:,2) - tmp/SIZE(E%sigEF)*nvac
          END IF
        ELSE IF (E%sig_b(2) /= 0.0) THEN
           CALL juDFT_error&
     &          ("z-mirror/inversion symmetry but second e-field"//&
     &          "sheet specified",calledby ="efield")
        END IF

        IF (ALL (E%sigEF == 0.0)) THEN
          DEALLOCATE (E%sigEF)
        ELSE
          E%l_segmented = .TRUE.
        END IF
      END SUBROUTINE read_efield

      SUBROUTINE read_shape (E, orig_str, nvac)
        USE m_constants, ONLY : pimach
        TYPE(t_efield), INTENT(INOUT) :: E
        CHARACTER(*), INTENT(IN) :: orig_str
        INTEGER, INTENT(IN) :: nvac

        REAL :: TWO_PI
        INTEGER :: ipos, action, iopt, ivac, ix, iy, nx, ny, cnt, i, j
        INTEGER :: np
        CHARACTER(10) :: tag, pos
        CHARACTER(200) :: str
        REAL :: x, y, h, w, radius, charge, order, shift, tmp
        REAL, ALLOCATABLE :: poly(:,:)
        REAL :: data(UBOUND(E%sigEF, dim=1), UBOUND(E%sigEF, dim=2))
        LOGICAL :: mask(UBOUND(E%sigEF, dim=1), UBOUND(E%sigEF, dim=2))

        TWO_PI = 2.0 * pimach()
        str = lower_case (orig_str)

        IF (str(1:8) == 'rectsinx') THEN
          READ (str, *) tag, pos, x, y, w, h, charge, order, shift
          IF (tag /= 'rectsinx')&
     &        CALL juDFT_error("Internal error in read_shape (rectSinX)"&
     &        ,calledby ="efield")
        ELSE IF (str(1:4) == 'circ') THEN
          READ (str, *) tag, pos, x, y, radius, charge
          IF (tag /= 'circ')  CALL juDFT_error&
     &         ("Internal error in read_shape (circ)",calledby ="efield"&
     &         )
        ELSE IF (str(1:4) == 'rect') THEN
          READ (str, *) tag, pos, x, y, w, h, charge
          IF (tag /= 'rect')  CALL juDFT_error&
     &         ("Internal error in read_shape (rect)",calledby ="efield"&
     &         )
        ELSE IF (str(1:7) == 'polygon') THEN
          READ (str, *) tag, pos, np
          ALLOCATE (poly(2,np))
          READ (str, *) tag, pos, np, poly, charge
          IF (tag /= 'polygon')&
     &         CALL juDFT_error("Internal error in read_shape (polygon)"&
     &         ,calledby ="efield")
        ELSE IF (str(1:8) == "datafile") THEN
          READ (str, *) tag, pos
          IF (tag /= 'datafile')  CALL juDFT_error&
     &         ("Internal error in read_shape",calledby ="efield")
        ELSE
           CALL juDFT_error("ERROR reading ",calledby ="efield")
        END IF

        IF (pos == 'top') THEN
          ipos = pTOP
        ELSEIF (pos == 'bot') THEN
          ipos = pBOT
        ELSEIF (pos == 'topbot' .OR. pos == 'bottop') THEN
          ipos = pTOPBOT
        ELSE
           CALL juDFT_error("ERROR reading ",calledby ="efield")
        END IF

        IF (nvac == 1 .AND. ipos /= pTOP)&
     &       CALL juDFT_error("ERROR reading ",calledby ="efield")

        action = pADD
        IF (INDEX (str, 'replace') > 0) action = pREPLACE
        iopt = pALL
        IF (INDEX (str, ' zero') > 0) iopt = pZERO
        IF (INDEX (str, 'nonzero') > 0) iopt = pNONZERO

        mask = .false.
        nx = UBOUND (mask, dim=1)
        ny = UBOUND (mask, dim=2)

        IF (tag == 'rect') THEN
          mask(MAX (FLOOR(x*nx-0.5)+2,1) : MIN (FLOOR((x+w)*nx+0.5),nx),&
     &         MAX (FLOOR(y*ny-0.5)+2,1) : MIN (FLOOR((y+h)*ny+0.5),ny))&
     &      = .true.
          WRITE (6, *) tag, pos2str (ipos), x, y, h, w, charge,&
     &                 action2str (action), opt2str (iopt)
        ELSE IF (tag == 'rectsinx') THEN
          mask(MAX (FLOOR(x*nx-0.5)+2,1) : MIN (FLOOR((x+w)*nx+0.5),nx),&
     &         MAX (FLOOR(y*ny-0.5)+2,1) : MIN (FLOOR((y+h)*ny+0.5),ny))&
     &      = .true.
          WRITE (6, *) tag, pos2str (ipos), x, y, h, w, charge, order,&
     &                 shift, action2str (action), opt2str (iopt)
        ELSE IF (tag == 'circ') THEN
          DO iy = 1, ny
            DO ix = 1, nx
              IF (SQRT ((REAL (ix)/nx - x)**2 + (REAL(iy)/ny - y)**2)&
     &            <= radius)&
     &          mask(ix, iy) = .true.
            END DO
          END DO
          WRITE (6, *) tag, pos2str (ipos), x, y, radius, charge,&
     &                 action2str (action), opt2str (iopt)
        ELSE IF (tag == 'datafile') THEN
          WRITE (6, '(1x,a,1x,a)', advance="no") tag, pos2str (ipos)
          call readDataFile (str, data)
          WRITE (6, *) action2str (action), opt2str (iopt)
          mask = .true.
        ELSE IF (tag == 'polygon') THEN
          DO iy = 1, ny
            DO ix = 1, nx
              mask(ix, iy) = in_polygon (np, poly, (/real(ix-1)/nx,&
     &                                               real(iy-1)/ny/))
            END DO
          END DO
        ELSE
           CALL juDFT_error("Internal error in read_shape",calledby&
     &          ="efield")
        END IF

        IF (action == pADD) THEN
          ! All relevant vacua
          DO ivac = 2 - MOD (ipos, 2), MIN (ipos, 2)
            IF (iopt == pZERO) THEN
              mask = mask .and. (E%sigEF(:,:,ivac) == 0.0)
            ELSE IF (iopt == pNONZERO) THEN
              mask = mask .and. (E%sigEF(:,:,ivac) /= 0.0)
            END IF

            IF (tag == 'rectsinx') THEN
              ! number of grid points in y direction
              IF (E%dirichlet) THEN
                cnt = 1
              ELSE
                cnt = MAXVAL (COUNT (mask, dim=2))
              END IF
              i = MAX (FLOOR(x*nx-0.5)+2,1)
              j = MIN (FLOOR((x+w)*nx+0.5),nx)
              DO ix = i, j
                tmp = REAL (ix-1)/REAL (j) ! Range: [0,1)
                DO iy = 1, ny
                  IF (.NOT. mask (ix, iy)) CYCLE
                  E%sigEF(ix,iy,ivac) = E%sigEF(ix,iy,ivac)&
     &                   + charge/cnt * SIN (TWO_PI*(order*tmp+shift))
                END DO
              END DO
            ELSE IF (tag == 'datafile') THEN
              WHERE (mask) E%sigEF(:,:,ivac) = E%sigEF(:,:,ivac) + data
            ELSE ! circ, rect
              IF (E%dirichlet) THEN
                cnt = 1
              ELSE
                cnt = COUNT (mask)
              END IF
              WHERE (mask) E%sigEF(:,:,ivac) = E%sigEF(:,:,ivac)&
     &                                         + charge/cnt
            END IF
         END DO ! ivac
        ELSE ! pREPLACE
          ! All relevant vacua
          DO ivac = 2 - MOD (ipos, 2), MIN (ipos, 2)
            IF (iopt == pZERO) THEN
              mask = mask .and. (E%sigEF(:,:,ivac) == 0.0)
            ELSE IF (iopt == pNONZERO) THEN
              mask = mask .and. (E%sigEF(:,:,ivac) /= 0.0)
            END IF

            IF (tag == 'rectsinx') THEN
              ! number of grid points in y direction
              IF (E%dirichlet) THEN
                cnt = 1
              ELSE
                cnt = MAXVAL (COUNT (mask, dim=2))
              END IF
              i = MAX (FLOOR(x*nx-0.5)+2,1)
              j = MIN (FLOOR((x+w)*nx+0.5),nx)
              DO ix = i, j
                tmp = REAL (ix-1)/REAL (j) ! Range: [0,1)
                DO iy = 1, ny
                  IF (.NOT. mask (ix, iy)) CYCLE
                    E%sigEF(ix,iy,ivac) =&
     &                     + charge/cnt * SIN (TWO_PI*(order*tmp+shift))
                END DO
              END DO
            ELSE IF (tag == 'datafile') THEN
              WHERE (mask) E%sigEF(:,:,ivac) = data
            ELSE ! circ, rect
              IF (E%dirichlet) THEN
                cnt = 1
              ELSE
                cnt = COUNT (mask)
              END IF
             WHERE (mask) E%sigEF(:,:,ivac) = charge/cnt
            END IF
          END DO ! ivac
        END IF
      END SUBROUTINE read_shape

      ! Read electric-field data points from "file"
      ! Format: First line, number of x and y points
      !         Second line: charge(x=1,y=1)
      !         Third line:  charge(x=1,y=2)
      !         etc.
      SUBROUTINE readDataFile (str, data)
        CHARACTER(*), intent(in) :: str
        REAL, intent(out) :: data(:,:)

        integer :: i, j, nx, ny
        logical :: l_exist
        character(len (str)) :: file, dummy

        read(str, *) dummy, dummy, file

        nx = ubound(data, dim=1)
        ny = ubound(data, dim=2)

        INQUIRE (file=file, exist=l_exist)
        IF (.not. l_exist) GOTO 110

        OPEN (1243, file=file, status='old')
        READ (1243, *, end=110) i, j  ! nx ny
        IF (i /= nx .or. j /= ny) THEN
          WRITE (*,'(3a)') "ERROR: Invalid number of points in '",&
     &           trim(file), "' during electrical field dataFile"
          WRITE (*,"(a,i0,a,i0,a)") "Expected: ", nx, ",", ny
          WRITE (*,"(a,i0,a,i0,a)") "Read    : ", i, ", ", j
          CALL juDFT_error("ERROR: dataFile of electric fiel&
     &d",calledby ="efield")
        END IF

        DO j = 1, ny
          DO i = 1, nx
            READ (1243, *) data(i,j)
          END DO
        END DO
        CLOSE (1243)

        WRITE (6, '(1x,3a,1x)', advance="no") '"', trim (file), '"'
        IF (INDEX (lower_case (str), 'zero_avg') > 0) THEN
          data(:,:) = data - sum(data)/(nx*ny)
          WRITE (6, '(a,1x)', advance="no") "zero_avg"
        END IF

        RETURN ! No error

        ! ERROR
110     WRITE (*,'(3a)') "ERROR: While reading '",trim(file),&
     &                   "' during electrical field dataFile"
        WRITE (*,"(a,i0,a,i0,a)") "Expecting  ", nx, ",", ny, " points"
        CALL juDFT_error("ERROR: dataFile of electric field",calledby&
     &       ="efield")
      END SUBROUTINE readDataFile


! Give a polygon and a trial point. Now a line is taken from the edge
! of the grid and it is counted how many times an edges of the polygon
! is crossed. If it was an odd time until one reaches the point, it is
! in the polygon - otherwise not.
! Note: The technique does not always find the optimal polygon.
!
! Cf. http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
!
      PURE FUNCTION in_polygon (n, poly, point)
        INTEGER, intent(in) :: n
        REAL, intent(in) :: poly(2,n), point(2)
        LOGICAL :: in_polygon
        INTEGER :: i, j

        i = 1
        j = n
        in_polygon = .FALSE.
        DO
          IF (i > n) EXIT
          IF ((poly(2,i) > point(2)) .NEQV. (poly(2,j) > point(2))) THEN
            IF (point(1) <  (poly(1,j) - poly(1,i))&
     &                    * (point(2)-poly(2,i)) / (poly(2,j)-poly(2,i))&
     &                    + poly(1,i))&
     &      in_polygon = .NOT. in_polygon
          END IF
          j = i
          i = i+1
        END DO
      END FUNCTION in_polygon

      PURE FUNCTION pos2str (ipos)
        INTEGER, INTENT(IN) :: ipos
        CHARACTER(6) :: pos2str
        IF (ipos == pTOP) THEN
          pos2str = 'Top'
        ELSE IF (ipos == pBOT) THEN
          pos2str = 'Bot'
        ELSE
          pos2str = 'TopBot'
        END IF
      END FUNCTION pos2str


      PURE FUNCTION action2str (action)
        INTEGER, INTENT(IN) :: action
        CHARACTER(7) :: action2str
        IF (action == pADD) THEN
          action2str = 'Add'
        ELSE
          action2str = 'Replace'
        END IF
      END FUNCTION action2str


      PURE FUNCTION opt2str (iopt)
        INTEGER, INTENT(IN) :: iopt
        CHARACTER(7) :: opt2str
        IF (iopt == pALL) THEN
          opt2str = 'All'
        ELSE IF (iopt == pZERO) THEN
          opt2str = 'Zero'
        ELSE
          opt2str = 'NonZero'
        END IF
      END FUNCTION opt2str


      SUBROUTINE print_efield (unit, E, area, nvac, amat,dvac)
          INTEGER, INTENT(IN)        :: unit, nvac
          TYPE(t_efield), INTENT(IN) :: E
          REAL, INTENT(IN)           :: area, amat(3,3),dvac

          ! electric constant
          REAL, PARAMETER :: epsilon0 = 8.854187817620e-12 ! F/m
          ! Bohr radius
          REAL, PARAMETER :: a0       = 0.52917720859e-10  ! m
          ! elemental charge
          REAL, PARAMETER :: ec       = 1.602176487e-19    ! C

          ! htr -> electron Volt
          REAL, PARAMETER :: htr_eV   = 27.21138386 ! eV

          ! Conversion of surface charge sigma to an electric field
          ! near a plate: E = -sigma[1/m^2]*e/eps0
          !                 = -sigma[1/a0^2]*e/(eps0 * a0^2)
          ! Assuming field in V/cm and charge in atomic units
          ! The 10^-2 is due to centimetres
          !
          ! NOTE: This equation assumes that one has two plates, one
          ! charged "sigma" and one "-sigma"; then E = sigma/eps0.
          ! If one plate were 0 and the other sigma, one had to use
          !   E = sigma/(2*eps0)
          ! Cf. http://hyperphysics.phy-astr.gsu.edu/hbase/electric/elesht.html
          !
          ! BE CAREFUL TO AVOID DOUBLE COUNTING
          REAL, PARAMETER :: sig_to_E = -ec/(epsilon0*a0**2)*1e-2

          REAL :: tmp, pt_rel(3), pt_abs(3)
          INTEGER :: ivac, i, j, nx, ny

          IF (SUM (ABS (E%sig_b)) < 1e-15&
             &.AND. ABS (E%sigma) < 1e-15&
             &.AND. .NOT. ALLOCATED (E%sigEF)&
     &        .AND. .NOT. E%dirichlet) RETURN

          WRITE (unit,'(3x,a,2(f12.5,a))')'z-z1 of external sheet    =',&
     &          E%zsigma, ' a0 = ', E%zsigma*a0*1e10, ' A'


          IF (E%dirichlet) THEN

          IF (E%sigma > 0)&
     &      WRITE (unit,'(3x,a,f12.5)') 'Surplus charge: ', 2.0*E%sigma

          IF (ALLOCATED(E%sigEF))&
     &      WRITE (unit,'(3x,a)') 'Average potential:'
          WRITE (unit, '(3x,a,f12.5,a, f12.5,a)') 'on sheet 1: ',&
     &       E%sig_b(1),' htr = ',E%sig_b(1)*hartree_to_ev_const,' V'
          IF (nvac > 1) THEN
            WRITE (unit, '(3x,a,f12.5,a, f12.5,a)') 'on sheet 2: ',&
     &         E%sig_b(2),' htr = ',E%sig_b(2)*hartree_to_ev_const,' V'

            WRITE (unit,'(3x,a,f14.5,a)')&
     &            'Average field (plate to plate):',&
     &            (E%sig_b(2)-E%sig_b(1))/((2*E%zsigma+Dvac)*a0*100),&
     &            ' V/cm'
          END IF ! nvac > 1

          ELSE ! Neumann

          IF (ALLOCATED(E%sigEF))&
     &      WRITE (unit,'(3x,a)') 'Average charges:'
          tmp = E%sigma+E%sig_b(1)
          WRITE (unit, '((3x,a,f12.5,5x,a,f12.5,a))')&
     &          'charge on external sheet 1 =', tmp,&
     &          '(surface density=', tmp/area, '  e/a.u.**2)'
          tmp = tmp/area*sig_to_E
          WRITE (unit, '(3x,a,e13.5,5x,a)')&
     &          'external field on sheet 1 =', tmp, 'V/cm'

          IF (nvac > 1) THEN
            tmp = E%sigma+E%sig_b(2)
            WRITE (unit, '((3x,a,f12.5,5x,a,f12.5,a))')&
     &            'charge on external sheet 2 =', tmp,&
     &            '(surface density=', tmp/area, '  e/a.u.**2)'
            tmp = tmp/area*sig_to_E
            WRITE (unit, '(3x,a,e13.5,5x,a)')&
     &            'external field on sheet 2 =', tmp, 'V/cm'

            ! Cf. comment above, where "sig_to_E" is defined
            WRITE (unit, '(a,/,a)') 'NOTE: The equation for the E field'&
     &         // ' assumes two oppositely charged plates.',&
     &            '      You may need to divide by two before '&
     &         //'summing the external fields to avoid double counting'
          ELSE
            WRITE (unit, '(a)') 'NOTE: The equation for the E field '&
     &            // 'already assumes two oppositely charged plates - '&
     &            // 'be careful to avoid double counting'
          END IF ! nvac > 1

          END IF ! Neumann (vs. Dirichlet)

          IF (ALLOCATED (E%sigEF) .AND. E%plot_charge) THEN
            nx = UBOUND (E%sigEF, dim=1) - 1
            ny = UBOUND (E%sigEF, dim=2) - 1
            DO ivac = 1, nvac
              IF (ivac == 1) THEN
                OPEN (748, file='efield-1.dat')
              ELSE
                OPEN (748, file='efield-2.dat')
              END IF
              IF (E%dirichlet) THEN
                WRITE (748, '(3a)') '# X[a0]     Y[a0]   X[rel]   ',&
     &                'Y[rel]    potential[htr]  potential[V]'
              ELSE ! Neumann
                WRITE (748, '(3a)') '# X[a0]     Y[a0]   X[rel]   ',&
     &                'Y[rel]    charge per grid point',&
     &                '    E_field (V/cm) in vicinity of the grid point'
              END IF
              pt_rel(3) = 0.0
              DO i = 1, nx+1
                pt_rel(1) = REAL(i-1)/nx
                DO j = 1, ny+1
                  pt_rel(2) = REAL(j-1)/ny
                  !CALL cotra0 (pt_rel, pt_abs, amat)
                  pt_abs=matmul(amat,pt_rel)
                  IF (E%dirichlet) THEN
                    WRITE (748, '(4f12.5,2g16.5)')&
     &                pt_abs(1:2), pt_rel(1:2),&
     &                E%sigEF(i,j,ivac)&
     &                + E%sig_b(ivac),&
     &                (E%sigEF(i,j,ivac)&
     &                + E%sig_b(ivac))*hartree_to_ev_const
                  ELSE ! Neumann
                    WRITE (748, '(4f12.5,2g16.5)')&
     &                pt_abs(1:2), pt_rel(1:2),&
     &                E%sigEF(i,j,ivac)&
     &                + (E%sigma + E%sig_b(ivac))/((nx+1)*(ny+1)),&
     &                (E%sigEF(i,j,ivac)*((nx+1)*(ny+1))&
     &                + E%sigma + E%sig_b(ivac))/area*sig_to_E
                  END IF
                END DO
                WRITE (748, *)
              END DO
              CLOSE (748)
            END DO
          END IF
        END SUBROUTINE print_efield

        SUBROUTINE V_seg_EF(&
     &                      efield,&
     &                      vacuum, stars)
          USE m_fft2d
          use m_types
          ! Dummy variables:
          TYPE(t_efield), INTENT(INOUT) :: efield
          TYPE(t_vacuum), INTENT(IN) :: vacuum
          TYPE(t_stars), INTENT(IN) :: stars

          ! Local variables:
          INTEGER :: i, ivac
          REAL :: fg, fgi
          REAL :: rhoRS(3*stars%mx1,3*stars%mx2), rhoRSimag(3*stars%mx1,3*stars%mx2) ! Real space density

          ALLOCATE (efield%rhoEF (stars%ng2-1, vacuum%nvac))

          DO ivac = 1, vacuum%nvac
            rhoRSimag = 0.0 ! Required in the loop. fft2d overrides this
            ! The fft2d algorithm puts the normalization to the  isn=-1
            ! (r->k) transformation, thus we need to multiply by
            ! ifft2 = 3*k1d*3*k2d to compensate.
            IF (efield%dirichlet) THEN
              rhoRS(:,:) = efield%sigEF(:,:,ivac)
            ELSE
              rhoRS(:,:) = efield%sigEF(:,:,ivac)*(3*stars%mx1*3*stars%mx2)
            END IF

            ! FFT rhoRS(r_2d) -> rhoEF(g_2d)
            CALL fft2d (&
     &                  stars,&
     &                  rhoRS, rhoRSimag,&
     &                  fg, fgi,&
     &                  efield%rhoEF(:,ivac), 1, -1)
!           FFT gives the the average charge per grid point
!           while sig_b stores the (total) charge per sheet
            IF (efield%dirichlet .and. ABS (fg) > 1.0e-15) THEN
              PRINT *, 'INFO: Difference of average potential: fg=',&
     &                 fg,', sig_b=', efield%sig_b(ivac),&
     &                 ", ivac=", ivac
            ELSE IF (ABS (fg/(3*stars%mx1*3*stars%mx2)) > 1.0e-15) THEN
              PRINT *, 'INFO: Difference of average potential: fg=',&
     &                 fg/(3*stars%mx1*3*stars%mx2),', sig_b=', efield%sig_b(ivac),&
     &                 ", ivac=", ivac
            END IF
          END DO ! ivac
        END SUBROUTINE V_seg_EF


        SUBROUTINE print_rhoEF(&
     &                         efield, k1d, k2d, nvac, n2d, nop, nop2,&
     &                         nq2, kv2, mrot, symor, tau, invtab, area,&
     &                         nstr2, amat)
          USE m_starf, ONLY: starf2
          ! Arguments
          TYPE(t_efield), INTENT(IN) :: efield
          INTEGER, INTENT(IN) :: k1d, k2d, nvac, n2d, nq2, nop, nop2
          REAL,    INTENT (IN) :: area
          LOGICAL, INTENT (IN) :: symor
          INTEGER, INTENT (IN) :: nstr2(n2d), kv2(2,n2d), mrot(3,3,nop),&
     &                            invtab(nop)
          REAL,    INTENT (IN) :: tau(3,nop), amat(3,3)

          ! Local variables
          real :: pt_rel(3), pt_abs(3), rcc(3)
          integer :: ix, iy, ivac, k, nix, niy
          real :: rhoOut, rhoOutIm
          complex :: sf2(n2d)

          nix = 3*k1d-1
          niy = 3*k2d-1

          DO ivac = 1,nvac
            IF (ivac == 1) THEN
              OPEN (754,file='efield-fft-1.dat')
            ELSE
              OPEN (754,file='efield-fft-2.dat')
            END IF

            IF (efield%dirichlet) THEN
              WRITE (754, '(2a)') '# X[a0]   Y[a0]  X[rel]  Y[rel]   ',&
     &           'potential per grid point (Re, Im, Abs), grid (ix, iy)'
            ELSE ! Neumann
              WRITE (754, '(2a)') '# X[a0]   Y[a0]  X[rel]  Y[rel]   ',&
     &           'charge per grid point (Re, Im, Abs), grid (ix, iy)'
            END IF
            pt_rel(3) = 0.0
            DO ix = 0, nix
              pt_rel(1) = REAL (ix)/REAL (nix)
              DO iy = 0, niy
                pt_rel(2) = REAL (iy)/REAL (niy)
                CALL starf2 (&
     &                       nop2,nq2,kv2,mrot,symor,tau,pt_rel,invtab,&
     &                       sf2)
                IF (efield%dirichlet) THEN
                  rhoOut = efield%sig_b(ivac)
                ELSE
                  rhoOut = (efield%sigma + efield%sig_b(ivac))/area
                END IF
                rhoOutIm = 0.0
                ! As we moved the normalization to g-space while the
                ! fft2d routine works in r space, we need to compensate
                ! by dividing by ifft2 = 3*k1d*3*k2d
                DO k = 2, nq2
                  IF (efield%dirichlet) THEN
                    rhoOut = rhoOut&
     &                       + REAL (efield%rhoEF(k-1,ivac)*sf2(k))&
     &                         * nstr2(k)
                    rhoOutIm = rhoOutIm&
     &                        + AIMAG (efield%rhoEF(k-1,ivac)*sf2(k))&
     &                          * nstr2(k)
                  ELSE
                    rhoOut = rhoOut&
     &                       + REAL (efield%rhoEF(k-1,ivac)*sf2(k))&
     &                         * nstr2(k)/(3*k1d*3*k2d)
                    rhoOutIm = rhoOutIm&
     &                        + AIMAG (efield%rhoEF(k-1,ivac)*sf2(k))&
     &                          * nstr2(k)/(3*k1d*3*k2d)
                  END IF
                END DO ! k
                !CALL cotra0 (pt_rel, pt_abs, amat)
                pt_abs=matmul(amat,pt_rel)
                WRITE (754,'(4f14.7,3g16.7,2i6)')&
     &            pt_abs(1:2),pt_rel(1:2),&
     &            rhoOut, rhoOutIm, SQRT (rhoOut**2+rhoOutIm**2), ix, iy
              END DO ! iy
              WRITE (754,*)
            END DO ! ix
            CLOSE (754)
          END DO ! ivac
        END SUBROUTINE print_rhoEF
      END SUBROUTINE e_field

      SUBROUTINE read_namelist (iou, E, eV)
        USE m_types, only: t_efield
        INTEGER, INTENT(IN) :: iou
        TYPE(t_efield), INTENT(INOUT) :: E
        LOGICAL, INTENT(INOUT) :: eV

        REAL :: zsigma, sig_b(2)
        LOGICAL :: plot_charge
        LOGICAL :: plot_rho
        LOGICAL :: autocomp
        LOGICAL :: dirichlet
        NAMELIST /efield/ zsigma, sig_b, plot_charge, plot_rho,&
     &                    autocomp, dirichlet, eV

        zsigma = E%zsigma
        plot_charge = E%plot_charge
        plot_rho = E%plot_rho
        autocomp = E%autocomp
        dirichlet = E%dirichlet
        sig_b = E%sig_b

        READ (iou,nml=efield)

        E%zsigma = zsigma
        E%plot_charge = plot_charge
        E%plot_rho = plot_rho
        E%autocomp = autocomp
        E%dirichlet = dirichlet
        E%sig_b = sig_b
      END SUBROUTINE read_namelist


      ELEMENTAL FUNCTION lower_case(string)
        CHARACTER(len=*), INTENT(IN) :: string
        CHARACTER(len=len(string))   :: lower_case

        INTEGER :: i

        DO i = 1, LEN (string)
          IF (      IACHAR ('A') <= IACHAR (string(i:i)) &
              .and. IACHAR ('Z') >= IACHAR (string(i:i))) THEN
            lower_case(i:i) = ACHAR (IACHAR (string(i:i))&
                              + IACHAR ('a') - IACHAR ('A'))
          ELSE
            lower_case(i:i) = string(i:i)
          END IF
        END do
      END FUNCTION lower_case
    END MODULE m_efield
