!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


module m_dfpt_gradient
    use m_juDFT
    use m_constants
    use m_types

    implicit none

contains

    subroutine dfpt_generate_gradient(fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVc3,grVext3,grgrVext3x3)
        
        use m_dfpt_vgen
        use m_vgen_coulomb
        use m_dfpt_potdenLocal
        use m_grdchlh

        type(t_fleurinput),intent(in) :: fi 
        type(t_mpi), intent(in)       :: fmpi
        type(t_sphhar),intent(in)     :: sphhar
        type(t_hybdat),intent(in)     :: hybdat
        class(t_xcpot), intent(in)     :: xcpot
        type(t_nococonv), intent(in)  :: nococonv
        type(t_stars), intent(in)     :: stars
        type(t_potden), intent(in)    :: rho
        type(t_potden), intent(in)    :: vTot
        type(t_potden), intent(inout)    :: grRho3(3)
        type(t_potden), intent(inout)    :: grVtot3(3)
        type(t_potden), intent(inout)    :: grVc3(3)
        type(t_potden), intent(inout)    :: grVext3(3)
        type(t_potden),intent(inout)     :: grgrVext3x3(3,3)



        type(t_stars) :: starsLocal
        type(t_potden) :: potdummy,potdummyLocal
        type(t_potden) :: rho_tmp ! a copy of the starting density. This is done to not mess with the starting density 
        type(t_atoms)  :: atomsLocal
        integer :: iDir, iDir2, iSpin, xInd, yInd, zInd, iStar , iVac, zlim
        complex :: sigma_loc(2)

        complex, allocatable :: grrhodummy(:, :, :, :, :) 
        real    :: dr_re(fi%vacuum%nmzd), dr_im(fi%vacuum%nmzd), drr_dummy(fi%vacuum%nmzd)

        ! remove this asap should not be needed any longer
        sigma_loc = cmplx(0.0,0.0)


        call rho_tmp%copyPotDen(rho)

        ALLOCATE(grrhodummy(fi%atoms%jmtd, (fi%atoms%lmaxd+1)**2, fi%atoms%nat, SIZE(rho%mt,4), 3))

        ! For the gradient of the external Potential we need a higher cutoff in the expansion, in order to confine the multipole moments in the MTs
        ! This is crucial for the Film-Mode and will be visible in the z-Eigenmodes.  
        call create_typesLocal(fi,fmpi,fi%sym,fi%cell,fi%input,sphhar, fi%vacuum , fi%noco ,starsLocal, grVext3(1),atomsLocal)
        call grVext3(2)%copyPotDen(grVext3(1))
        call grVext3(3)%copyPotDen(grVext3(1))

        ! initialize all potentials
        ! currently we do this by copy, is this faster than just init?

        call potdummy%copyPotDen(rho)
        call potdummy%resetpotden()
        call potdummyLocal%copyPotDen(grVext3(1))
        call potdummyLocal%resetpotden()


        do iDir = 1 , 3 
            call grRho3(iDir)%copyPotDen(rho)
            call grRho3(iDir)%resetpotden()
            do iDir2 = 1 , 3 
                call grgrVext3x3(iDir2,iDir)%copyPotDen(grVext3(1))
                call grgrVext3x3(iDir2,iDir)%resetpotden()
            end do 
            call grVtot3(iDir)%copyPotDen(vTot)
            call grVtot3(iDir)%resetpotden()
            call grVc3(iDir)%copyPotDen(vTot)
            call grVc3(iDir)%resetpotden() 
        end do 

        ! compute numerical gradient of the density 

        do iSpin = 1, size(rho%mt,4)
            call mt_gradient_new(fi%atoms, sphhar, fi%sym, rho%mt(:, :, :, iSpin), grrhodummy(:, :, :, iSpin, :))
        end do 
        do  zInd = -stars%mx3, stars%mx3
            do yInd = -stars%mx2, stars%mx2
                do xInd = -stars%mx1, stars%mx1
                    iStar = stars%ig(xInd, yInd, zInd)
                    if (iStar.eq.0) cycle
                    grRho3(1)%pw(iStar,:) = rho%pw(iStar,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,zInd]),fi%cell%bmat)))
                    grRho3(2)%pw(iStar,:) = rho%pw(iStar,:) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,zInd]),fi%cell%bmat)))
                    grRho3(3)%pw(iStar,:) = rho%pw(iStar,:) * cmplx(0.0,dot_product([0.0,0.0,1.0],matmul(real([xInd,yInd,zInd]),fi%cell%bmat)))
                end do 
            end do 
        end do 

        if (fi%input%film) then
            do yInd = -stars%mx2, stars%mx2
                do xInd = -stars%mx1, stars%mx1
                    iStar = stars%ig(xInd, yInd, 0)
                    if (iStar.EQ.0) cycle
                    iStar = stars%i2g(xInd, yInd)
                    grRho3(1)%vac(:,iStar,:,:) = rho%vac(:,iStar,:,:) * cmplx(0.0,dot_product([1.0,0.0,0.0],matmul(real([xInd,yInd,0]),fi%cell%bmat)))
                    grRho3(2)%vac(:,iStar,:,:) = rho%vac(:,iStar,:,:) * cmplx(0.0,dot_product([0.0,1.0,0.0],matmul(real([xInd,yInd,0]),fi%cell%bmat)))
                    do iVac = 1, fi%vacuum%nvac
                        do iSpin = 1, SIZE(rho%vac,4)
                            zlim = merge(fi%vacuum%nmz,fi%vacuum%nmzxy,iStar==1)
                            call grdchlh(fi%vacuum%delz, REAL(rho%vac(:zlim,iStar,iVac,iSpin)),dr_re(:zlim),drr_dummy(:zlim))
                            call grdchlh(fi%vacuum%delz,AIMAG(rho%vac(:zlim,iStar,iVac,iSpin)),dr_im(:zlim),drr_dummy(:zlim))
                            grRho3(3)%vac(:,iStar,iVac,iSpin) = (3-2*iVac)*(dr_re + ImagUnit * dr_im)
                        end do 
                    end do 
                end do 
            end do 
        end if 

        do iDir = 1 ,3 
            call sh_to_lh(fi%sym, fi%atoms, sphhar, SIZE(rho%mt,4), 2, grrhodummy(:, :, :, :, iDir), grRho3(iDir)%mt, potdummy%mt)
        end do 


        !  compute gradient external potential
        do iDir = 1 , 3 
            call vgen_coulomb(1, fmpi, fi%input, fi%field, fi%vacuum, fi%sym, fi%juphon, starsLocal, fi%cell, &
                         & sphhar, atomsLocal, .FALSE., potdummyLocal, grVext3(iDir), sigma_loc, &
                         & dfptdenimag=potdummyLocal, dfptvCoulimag=potdummyLocal,dfptden0=potdummyLocal,stars2=starsLocal,iDtype=0,iDir=iDir)
        end do 

        ! Coulomb/Effective potential gradients
        do iDir = 1, 3
            call potdummy%resetPotDen()
            sigma_loc  = cmplx(0.0,0.0)
            call dfpt_vgen(hybdat, fi%field, fi%input, xcpot, fi%atoms, sphhar, stars, fi%vacuum, fi%sym, &
                            fi%juphon, fi%cell, fmpi, fi%noco, nococonv, rho_tmp, vTot, &
                            stars, potdummy, grVtot3(iDir), .TRUE., potdummy, grRho3(iDir), 0, iDir, [0,0], sigma_loc)
            sigma_loc  = cmplx(0.0,0.0)
            call dfpt_vgen(hybdat, fi%field, fi%input, xcpot, fi%atoms, sphhar, stars, fi%vacuum, fi%sym, &
                            fi%juphon, fi%cell, fmpi, fi%noco, nococonv, rho_tmp, vTot, &
                            stars, potdummy, grVc3(iDir), .FALSE., potdummy, grRho3(iDir), 0, iDir, [0,0], sigma_loc)
        end do 


        ! Hessian of Potential for external Potential 

        do iDir2 = 1, 3
            do iDir = 1, 3
                sigma_loc = cmplx(0.0,0.0)
                call potdummyLocal%resetPotDen()
                call vgen_coulomb(1, fmpi, fi%input, fi%field, fi%vacuum, fi%sym, fi%juphon, starsLocal, fi%cell, &
                            & sphhar, atomsLocal, .TRUE., potdummyLocal, grgrVext3x3(iDir2,iDir), sigma_loc, &
                            & dfptdenimag=potdummyLocal, dfptvCoulimag=potdummyLocal,dfptden0=potdummyLocal,stars2=starsLocal,iDtype=0,iDir=iDir,iDir2=iDir2, &
                            & sigma_disc2=merge(sigma_loc,[cmplx(0.0,0.0),cmplx(0.0,0.0)],iDir2==3.AND.iDir==3.AND..FALSE.))
            end do 
        end do 




              
      call grRho3(1)%distribute(fmpi%mpi_comm)
      call grRho3(2)%distribute(fmpi%mpi_comm)
      call grRho3(3)%distribute(fmpi%mpi_comm)
      call grVext3(1)%distribute(fmpi%mpi_comm)
      call grVext3(2)%distribute(fmpi%mpi_comm)
      call grVext3(3)%distribute(fmpi%mpi_comm)

    end subroutine dfpt_generate_gradient 

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

    subroutine  sh_to_lh(sym, atoms, sphhar, jspins, radfact, rhosh, rholhreal, rholhimag)

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

      DO iSpin = 1, jspins
          DO iType = 1, atoms%ntype
              iAtom = atoms%firstAtom(iType)
              ptsym = sym%ntypsy(iAtom)
              DO ilh = 0, sphhar%nlh(ptsym)
                  l = sphhar%llh(iLH, ptsym)
                  DO iMem = 1, sphhar%nmem(ilh, ptsym)
                      m = sphhar%mlh(iMem, ilh, ptsym)
                      ilm = l * (l+1) + m + 1
                      DO iR = 1, atoms%jri(iType)
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

    end subroutine sh_to_lh

end module m_dfpt_gradient