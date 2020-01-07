!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_plot
   USE m_types
   USE m_juDFT
   USE m_constants

   IMPLICIT NONE

   !-----------------------------------------------------------------------------
   ! A general purpose plotting routine for FLEUR.
   !
   ! Based on the older plotting routines pldngen.f90 and plotdop.f90 that were
   ! originally called by optional.F90 and are now used in the scf-loop instead.
   ! At the cost of reduced postprocess functionality, this allowed us to re-
   ! move I/O (using plot.hdf/text files) from the plotting routine completely.
   !
   ! TODO:
   ! - plot_inp files are still in use and should be replaced by a plotting type.
   !
   ! A. Neukirchen & R. Hilgers, September 2019
   !------------------------------------------------

   PUBLIC :: checkplotinp, vectorsplit, matrixsplit, savxsf, vectorplot, &
             matrixplot, makeplots, procplot, getMTSphere

CONTAINS

   SUBROUTINE checkplotinp()

      !--------------------------------------------------------------------------
      ! Checks for existing plot input. If an ancient plotin file is used, an
      ! error is called. If no usable plot_inp exists, a new one is generated.
      !--------------------------------------------------------------------------

      LOGICAL :: oldform,newform

      oldform=.FALSE.
      INQUIRE(file="plotin",exist = oldform)

      IF (oldform) THEN
         CALL juDFT_error("Use of plotin file no longer supported", &
                          calledby="plot")
      END IF
      INQUIRE(file="plot_inp",exist = newform)
      IF (newform) THEN
        CALL juDFT_error("Use of plot_inp file no longer supported", &
        calledby="plot")
      END IF


   END SUBROUTINE checkplotinp

   SUBROUTINE vectorsplit(stars, atoms, sphhar, vacuum, input, factor, denmat, &
                          cden,mden)

      !--------------------------------------------------------------------------
      ! Takes a spin-polarized density vector and rearranges it into two plottable
      ! seperate ones, i.e. (rho_up, rho_down) ---> n, m.
      !
      ! Can also be applied to a potential with an additional factor, i.e. while
      ! the densities n/m are defined as (rho_up+-rho_down), V_eff/B_eff are de-
      ! fined as (V_up+-V_down)/2.
      !--------------------------------------------------------------------------

      TYPE(t_stars),  INTENT(IN)  :: stars
      TYPE(t_atoms),  INTENT(IN)  :: atoms
      TYPE(t_sphhar), INTENT(IN)  :: sphhar
      TYPE(t_vacuum), INTENT(IN)  :: vacuum
      TYPE(t_input),  INTENT(IN)  :: input
      REAL,           INTENT(IN)  :: factor
      TYPE(t_potden), INTENT(IN)  :: denmat
      TYPE(t_potden), INTENT(OUT) :: cden, mden

      IF (factor==1.0) THEN
         CALL cden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,&
                                      stars%ng2)
         CALL mden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,&
                                      stars%ng2)
      ELSE
         CALL cden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_POTTOT,vacuum%nmzd,&
                                      vacuum%nmzxyd,stars%ng2)
         CALL mden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_POTTOT,vacuum%nmzd,&
                                      vacuum%nmzxyd,stars%ng2)
      END IF

      cden%mt(:,0:,1:,1) = (denmat%mt(:,0:,1:,1)+denmat%mt(:,0:,1:,2))/factor
      cden%pw(1:,1) = (denmat%pw(1:,1)+denmat%pw(1:,2))/factor
      cden%vacz(1:,1:,1) = (denmat%vacz(1:,1:,1)+denmat%vacz(1:,1:,2))/factor
      cden%vacxy(1:,1:,1:,1) = (denmat%vacxy(1:,1:,1:,1)+denmat%vacxy(1:,1:,1:,2))/factor

      mden%mt(:,0:,1:,1) = (denmat%mt(:,0:,1:,1)-denmat%mt(:,0:,1:,2))/factor
      mden%pw(1:,1) = (denmat%pw(1:,1)-denmat%pw(1:,2))/factor
      mden%vacz(1:,1:,1) = (denmat%vacz(1:,1:,1)-denmat%vacz(1:,1:,2))/factor
      mden%vacxy(1:,1:,1:,1) = (denmat%vacxy(1:,1:,1:,1)-denmat%vacxy(1:,1:,1:,2))/factor

   END SUBROUTINE vectorsplit

   SUBROUTINE matrixsplit(sym,stars, atoms, sphhar, vacuum, input, noco, factor, &
                          denmat, cden, mxden, myden, mzden)
      USE m_fft2d
      USE m_fft3d
      USE m_rotdenmat

      !--------------------------------------------------------------------------
      ! Takes a 2x2 density matrix and rearranges it into four plottable seperate
      ! ones, i.e. ((rho_11, rho_12),(rho_21, rho_22)) ---> n, mx, my, mz.
      !
      ! This is an adaptation of the old pldngen.f90 by P. Kurz.
      !
      ! Can also be applied to a potential with additional factors, i.e. while
      ! the densities n/m are defined as (rho_up+-rho_down), V_eff/B_eff are de-
      ! fined as (V_up+-V_down)/2 and while rho_21 is of the form a-i*b, for the
      ! potential we have a+i*b.
      !
      ! TODO: Find out, whether the latter form of rho is still valid. Having the
      ! additive rho_21=m_x+i*m_y would be nicer and more convenient/consistent.
      !--------------------------------------------------------------------------

      IMPLICIT NONE
      TYPE(t_sym),    INTENT(IN)  :: sym
      TYPE(t_stars),  INTENT(IN)  :: stars
      TYPE(t_atoms),  INTENT(IN)  :: atoms
      TYPE(t_sphhar), INTENT(IN)  :: sphhar
      TYPE(t_vacuum), INTENT(IN)  :: vacuum
      TYPE(t_input),  INTENT(IN)  :: input
      TYPE(t_noco),   INTENT(IN)  :: noco
      REAL,           INTENT(IN)  :: factor
      TYPE(t_potden), INTENT(IN)  :: denmat
      TYPE(t_potden), INTENT(OUT) :: cden, mxden, myden, mzden

      ! Local scalars: iteration indices, matrix elements etc.
      INTEGER iden, ivac, ifft2, ifft3, imz, ityp, iri, ilh, imesh
      REAL cdnup, cdndown, chden, mgden, theta, phi, rhotot, mx, my, mz
      REAL zero, rziw, vz_r, vz_i, rho_11, rho_22, rho_21r, rho_21i, cdn11, cdn22
      COMPLEX czero, cdn21

      ! Local arrays: densities in real space and off-diagonal elements.
      REAL,    ALLOCATABLE        :: rvacxy(:,:,:,:), ris(:,:), fftwork(:)
      REAL,    ALLOCATABLE        :: rho(:,:,:,:), rht(:,:,:)
      COMPLEX, ALLOCATABLE        :: qpw(:,:), rhtxy(:,:,:,:)
      COMPLEX, ALLOCATABLE        :: cdom(:), cdomvz(:,:), cdomvxy(:,:,:)

      zero  = 0.0; czero = CMPLX(0.0,0.0)
      ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
      ifft2 = 9*stars%mx1*stars%mx2

      ! Allocation of arrays and initialization of those that make up the real
      ! space density matrix.
      ALLOCATE (rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,4), qpw(stars%ng3,4), &
                rht(vacuum%nmzd,2,4), rhtxy(vacuum%nmzxyd,stars%ng2-1,2,4), &
                cdomvz(vacuum%nmzd,2), cdomvxy(vacuum%nmzxyd,stars%ng2-1,2), &
                cdom(stars%ng3), fftwork(0:27*stars%mx1*stars%mx2*stars%mx3-1), &
                ris(0:27*stars%mx1*stars%mx2*stars%mx3-1,4), &
                rvacxy(0:9*stars%mx1*stars%mx2-1,vacuum%nmzxyd,2,4))

      rho(:,:,:,:) = zero; qpw(:,:) = czero; cdom(:) = czero

      IF (input%film) THEN
         cdomvz(:,:) = czero;    rhtxy(:,:,:,:) = czero
         cdomvxy(:,:,:) = czero; rht(:,:,:) = zero
      END IF

      rho(:,0:,1:,:input%jspins)    = denmat%mt(:,0:,1:,:input%jspins)
      qpw(1:,:input%jspins)         = denmat%pw(1:,:input%jspins)
      rht(1:,1:,:input%jspins)      = denmat%vacz(1:,1:,:input%jspins)
      rhtxy(1:,1:,1:,:input%jspins) = denmat%vacxy(1:,1:,1:,:input%jspins)

      IF(noco%l_noco) THEN
         cdom = denmat%pw(:,3)
         cdomvz(:,:) = CMPLX(denmat%vacz(:,:,3),denmat%vacz(:,:,4))
         cdomvxy = denmat%vacxy(:,:,:,3)
      END IF

      ! Calculate the charge and magnetization densities in the muffin tins.
      DO ityp = 1,atoms%ntype
         theta   = noco%beta(ityp)
         phi     = noco%alph(ityp)
         DO ilh = 0,sphhar%nlh(sym%ntypsy(ityp))
            DO iri = 1,atoms%jri(ityp)
               IF (SIZE(denmat%mt,4).LE.2) THEN
                  cdnup   = rho(iri,ilh,ityp,1)
                  cdndown = rho(iri,ilh,ityp,2)

                  chden   = cdnup + cdndown
                  mgden   = cdnup - cdndown
                  rho(iri,ilh,ityp,1) = chden
                  rho(iri,ilh,ityp,2) = mgden*COS(phi)*SIN(theta)
                  rho(iri,ilh,ityp,3) = mgden*SIN(phi)*SIN(theta)
                  rho(iri,ilh,ityp,4) = mgden*COS(theta)
               ELSE
                  cdn11 = rho(iri,ilh,ityp,1)
                  cdn22 = rho(iri,ilh,ityp,2)
                  cdn21 = CMPLX(denmat%mt(iri,ilh,ityp,3), &
                                denmat%mt(iri,ilh,ityp,4))

                  CALL rot_den_mat(noco%alph(ityp), noco%beta(ityp), &
                                   cdn11,cdn22,cdn21)

                  rho(iri,ilh,ityp,1) = cdn11 + cdn22
                  rho(iri,ilh,ityp,2) = 2.0*REAL(cdn21)
                  !! Note: The minus sign in the following line is temporary to
                  !!       adjust for differences in the offdiagonal part of the
                  !!       density between this FLEUR version and ancient v0.26 .
                  !!
                  !! TODO: Should that still be here? It effectively amounts to a
                  !!       conjugation of the density matrix.
                  rho(iri,ilh,ityp,3) = -2.0*AIMAG(cdn21)
                  rho(iri,ilh,ityp,4) = cdn11 - cdn22
               END IF
            END DO
         END DO
      END DO

      ! Fourier transform the diagonal part of the density matrix in the
      ! interstitial (qpw) to real space (ris).
      DO iden = 1,2
         CALL fft3d(ris(0,iden),fftwork,qpw(1,iden),stars,1)
      END DO

      ! Also do that for the off-diagonal part. Real part goes into index
      ! 3 and imaginary part into index 4.
      CALL fft3d(ris(0,3),ris(0,4),cdom(1),stars,1)

      ! Calculate the charge and magnetization densities in the interstitial.
      DO imesh = 0,ifft3-1
         rho_11  = ris(imesh,1)
         rho_22  = ris(imesh,2)
         rho_21r = ris(imesh,3)
         rho_21i = ris(imesh,4)
         rhotot  = rho_11 + rho_22
         mx      =  2*rho_21r
         my      = -2*rho_21i
         mz      = (rho_11-rho_22)

         ris(imesh,1) = rhotot
         ris(imesh,2) = mx
         ris(imesh,3) = my
         ris(imesh,4) = mz
      END DO

      ! Invert the transformation to put the four densities back into
      ! reciprocal space.
      DO iden = 1,4
         fftwork=zero
         CALL fft3d(ris(0,iden),fftwork,qpw(1,iden),stars,-1)
      END DO

      ! As above, but for the vacuum components in xy- and z-direction.
      IF (input%film) THEN
         DO iden = 1,2
            DO ivac = 1,vacuum%nvac
               DO imz = 1,vacuum%nmzxyd
                  rziw = 0.0
                  CALL fft2d(stars, rvacxy(0,imz,ivac,iden), fftwork, &
                             rht(imz,ivac,iden), rziw, rhtxy(imz,1,ivac,iden), &
                             vacuum%nmzxyd, 1)
               END DO
            END DO
         END DO

         DO ivac = 1,vacuum%nvac
            DO imz = 1,vacuum%nmzxyd
               rziw = 0.0
               vz_r = REAL(cdomvz(imz,ivac))
               vz_i = AIMAG(cdomvz(imz,ivac))
               CALL fft2d(stars, rvacxy(0,imz,ivac,3), rvacxy(0,imz,ivac,4), &
                          vz_r, vz_i, cdomvxy(imz,1,ivac), vacuum%nmzxyd, 1)
            END DO
         END DO

         DO ivac = 1,vacuum%nvac
            DO imz = 1,vacuum%nmzxyd
               DO imesh = 0,ifft2-1
                  rho_11  = rvacxy(imesh,imz,ivac,1)
                  rho_22  = rvacxy(imesh,imz,ivac,2)
                  rho_21r = rvacxy(imesh,imz,ivac,3)
                  rho_21i = rvacxy(imesh,imz,ivac,4)
                  rhotot  = rho_11 + rho_22
                  mx      =  2*rho_21r
                  my      = -2*rho_21i
                  mz      = (rho_11-rho_22)

                  rvacxy(imesh,imz,ivac,1) = rhotot
                  rvacxy(imesh,imz,ivac,2) = mx
                  rvacxy(imesh,imz,ivac,3) = my
                  rvacxy(imesh,imz,ivac,4) = mz
               END DO
            END DO

            DO imz = vacuum%nmzxyd+1,vacuum%nmzd
               rho_11  = rht(imz,ivac,1)
               rho_22  = rht(imz,ivac,2)
               rho_21r = REAL(cdomvz(imz,ivac))
               rho_21i = AIMAG(cdomvz(imz,ivac))
               rhotot  = rho_11 + rho_22
               mx      =  2*rho_21r
               my      = -2*rho_21i
               mz      = (rho_11-rho_22)

               rht(imz,ivac,1) = rhotot
               rht(imz,ivac,2) = mx
               rht(imz,ivac,3) = my
               rht(imz,ivac,4) = mz
            END DO
         END DO

         DO iden = 1,4
            DO ivac = 1,vacuum%nvac
               DO imz = 1,vacuum%nmzxyd
                  fftwork=zero
                  CALL fft2d(stars, rvacxy(0,imz,ivac,iden), fftwork, &
                             rht(imz,ivac,iden), rziw, rhtxy(imz,1,ivac,iden), &
                             vacuum%nmzxyd, -1)
               END DO
            END DO
         END DO
      END IF



      ! Initialize and save the four output densities.
      IF (factor==1.0) THEN
         CALL cden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,&
                                      stars%ng2)
         CALL mxden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,&
                                      stars%ng2)
         CALL myden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,&
                                      stars%ng2)
         CALL mzden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_DEN,vacuum%nmzd,vacuum%nmzxyd,&
                                      stars%ng2)

         cden%mt(:,0:,1:,1) = rho(:,0:,1:,1)
         cden%pw(1:,1) = qpw(1:,1)
         cden%vacz(1:,1:,1) = rht(1:,1:,1)
         cden%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,1)

         mxden%mt(:,0:,1:,1) = rho(:,0:,1:,2)
         mxden%pw(1:,1) = qpw(1:,2)
         mxden%vacz(1:,1:,1) = rht(1:,1:,2)
         mxden%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,2)

         myden%mt(:,0:,1:,1) = rho(:,0:,1:,3)
         myden%pw(1:,1) = qpw(1:,3)
         myden%vacz(1:,1:,1) = rht(1:,1:,3)
         myden%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,3)

         mzden%mt(:,0:,1:,1) = rho(:,0:,1:,4)
         mzden%pw(1:,1) = qpw(1:,4)
         mzden%vacz(1:,1:,1) = rht(1:,1:,4)
         mzden%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,4)
      ELSE
         CALL cden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_POTTOT,vacuum%nmzd,&
                                      vacuum%nmzxyd,stars%ng2)
         CALL mxden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_POTTOT,vacuum%nmzd,&
                                      vacuum%nmzxyd,stars%ng2)
         CALL myden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_POTTOT,vacuum%nmzd,&
                                      vacuum%nmzxyd,stars%ng2)
         CALL mzden%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd,&
                                      atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE.,&
                                      POTDEN_TYPE_POTTOT,vacuum%nmzd,&
                                      vacuum%nmzxyd,stars%ng2)

         ! Correction for the case of plotting the total potential.
         ! Needed due to the different definitons of density/potential matrices in
         ! FLEUR:
         !
         ! rhoMat = 0.5*((n    +m_z,m_x+i*m_y),(m_x-i*m_y,n    -m_z))
         !   vMat =     ((V_eff+B_z,B_x-i*B_y),(B_x+i*m_y,V_eff-B_z))

         rho(:,0:,1:,:) = rho(:,0:,1:,:)/2.0
         qpw(1:,:) = qpw(1:,:)/2.0
         rht(1:,1:,:) = rht(1:,1:,:)/2.0
         rhtxy(1:,1:,1:,:) = rhtxy(1:,1:,1:,:)/2.0

         rho(:,0:,1:,3) = -rho(:,0:,1:,3)
         qpw(1:,3) = -qpw(1:,3)
         rht(1:,1:,3) = -rht(1:,1:,3)
         rhtxy(1:,1:,1:,3) = -rhtxy(1:,1:,1:,3)

         cden%mt(:,0:,1:,1) = rho(:,0:,1:,1)
         cden%pw(1:,1) = qpw(1:,1)
         cden%vacz(1:,1:,1) = rht(1:,1:,1)
         cden%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,1)

         mxden%mt(:,0:,1:,1) = rho(:,0:,1:,2)
         mxden%pw(1:,1) = qpw(1:,2)
         mxden%vacz(1:,1:,1) = rht(1:,1:,2)
         mxden%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,2)

         myden%mt(:,0:,1:,1) = rho(:,0:,1:,3)
         myden%pw(1:,1) = qpw(1:,3)
         myden%vacz(1:,1:,1) = rht(1:,1:,3)
         myden%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,3)

         mzden%mt(:,0:,1:,1) = rho(:,0:,1:,4)
         mzden%pw(1:,1) = qpw(1:,4)
         mzden%vacz(1:,1:,1) = rht(1:,1:,4)
         mzden%vacxy(1:,1:,1:,1) = rhtxy(1:,1:,1:,4)

      END IF

      DEALLOCATE (rho, qpw, rht, rhtxy, cdomvz, cdomvxy, &
                  cdom, fftwork, ris, rvacxy)

   END SUBROUTINE matrixsplit

   SUBROUTINE savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                     noco, score, potnorm, denName, denf, denA1, denA2, denA3)
      USE m_outcdn
      USE m_xsf_io

      ! Takes one/several t_potden variable(s), i.e. scalar fields in MT-sphere/
      ! plane wave representation and makes it/them into plottable .xsf file(s)
      ! according to a scheme given in plot_inp.
      !
      ! This is an adaptation of the old plotdop.f90 by P. Kurz.
      !
      ! Naming convention:
      !  The output .xsf files will have names composed of the plotted density
      !  (c.f. identifier) and an optional tag. For a scalar density, only said
      !  density is plotted. For a spin-polarized density with an up and down
      !  component, f denotes the sum of both densities and g denotes their
      !  difference u-d possibly with additional modifying factors as explained
      !  above in vector-/matrixsplit. f and g are to be understood as scalar
      !  fields f(r_vec) and g(r_vec). For a density matrix, f is still the
      !  scalar part. Aditionally, three vector components A_[1-3] arise.
      !
      !   E.g.: scalar density rho (n(r_vec))
      !          ---> rho.xsf (n(r_vec))
      !         spin-polarized density rho (n_up(r_vec),n_down(r_vec))
      !          ---> rho_f.xsf, rho_g.xsf [n(r_vec),m(r_vec)]
      !         matrix density rho ((n_11(r_vec),n_12(r_vec)),
      !                             (n_21(r_vec),n_22(r_vec)))
      !          ---> rho_f.xsf, rho_A1.xsf, rho_A2.xsf, rho_A3.xsf
      !                                [n(r_vec),m_vec(r_vec)]
      type(t_sliceplot),           INTENT(IN) :: sliceplot
      TYPE(t_stars),               INTENT(IN) :: stars
      TYPE(t_atoms),               INTENT(IN) :: atoms
      TYPE(t_sphhar),              INTENT(IN) :: sphhar
      TYPE(t_vacuum),              INTENT(IN) :: vacuum
      TYPE(t_input),               INTENT(IN) :: input
      TYPE(t_oneD),                INTENT(IN) :: oneD
      TYPE(t_sym),                 INTENT(IN) :: sym
      TYPE(t_cell),                INTENT(IN) :: cell
      TYPE(t_noco),                INTENT(IN) :: noco
      LOGICAL,                     INTENT(IN) :: score, potnorm
      CHARACTER(len=20),           INTENT(IN) :: denName
      TYPE(t_potden),              INTENT(IN) :: denf
      TYPE(t_potden),    OPTIONAL, INTENT(IN) :: denA1
      TYPE(t_potden),    OPTIONAL, INTENT(IN) :: denA2
      TYPE(t_potden),    OPTIONAL, INTENT(IN) :: denA3

      REAL    :: tec, qint, phi0, angss
      INTEGER :: i, j, ix, iy, iz, na, nplo, iv, iflag, nfile
      INTEGER :: nplot, nt, jm, jspin, numInDen, numOutFiles
      LOGICAL :: twodim, cartesian, xsf, unwind, polar

      TYPE(t_potden),     ALLOCATABLE :: den(:)
      REAL,               ALLOCATABLE :: xdnout(:)
      REAL                            :: pt(3), vec1(3), vec2(3), vec3(3), &
                                         zero(3), help(3), qssc(3), point(3)
      INTEGER                         :: grid(3)
      REAL                            :: rhocc(atoms%jmtd)
      CHARACTER (len=20), ALLOCATABLE :: outFilenames(:)
      CHARACTER (len=30)              :: filename
      CHARACTER (len=7)               :: textline

      REAL,               PARAMETER   :: eps = 1.0e-15

      NAMELIST /plot/twodim,cartesian,unwind,vec1,vec2,vec3,grid,zero,phi0,filename

      nfile = 120

      IF (PRESENT(denA2)) THEN
         ALLOCATE(den(4))
         den(1)          = denf
         den(2)          = denA1
         den(3)          = denA2
         den(4)          = denA3
         numInDen        = 4
         numOutFiles     = 4
      ELSE IF (PRESENT(denA1)) THEN
         ALLOCATE(den(2))
         den(1)          = denf
         den(2)          = denA1
         numInDen        = 2
         numOutFiles     = 2
      ELSE
         ALLOCATE(den(1))
         den(1)          = denf
         numInDen        = 1
         numOutFiles     = 1
      END IF

      DO i = 1, numInDen

         ! TODO: Does this work as intended?
         IF ((numInDen.NE.4).AND.(score)) THEN
            OPEN (17,file='cdnc',form='unformatted',status='old')
            REWIND 17
            DO jspin = 1, input%jspins
               DO nt = 1, atoms%ntype
                  jm = atoms%jri(nt)
                  READ (17) (rhocc(j),j=1,jm)
                  DO j = 1, jm
                     den(i)%mt(j,0,nt,jspin) = den(i)%mt(j,0,nt,jspin) - rhocc(j)/2.0/SQRT(pi_const)
                  END DO
                  READ (17) tec
               END DO
               READ (17) qint
               den(i)%pw(1,jspin) = den(i)%pw(1,jspin) - qint/cell%volint
            END DO
            CLOSE (17)
         ELSE IF (score) THEN
            CALL juDFT_error('Subtracting core charge in noco calculations not supported', calledby = 'plot')
         END IF
      END DO

      IF (noco%l_ss) THEN
         qssc = MATMUL(TRANSPOSE(cell%bmat),noco%qss)
      END IF

      ! Open the plot_inp file for input

      !OPEN (18,file='plot_inp')
      !READ(18,'(i2,5x,l1,1x,a)') nplot,xsf,textline
      polar = .FALSE.
      !IF ((noco%l_noco).AND.(numInDen.EQ.4)) THEN
      !   polar = (textline(1:7)=='polar=T').OR.(textline(1:7)=='polar=t')
      !     IF (polar) THEN
      !      numOutFiles = 7
      !   END IF
      !END IF
      xsf=sliceplot%format==PLOT_XSF_FORMAT
      ALLOCATE(outFilenames(numOutFiles))
      ALLOCATE(xdnout(numOutFiles))

      IF (numOutFiles.EQ.1) THEN
         outFilenames(1) = TRIM(denName)
      ELSE IF (numOutFiles.EQ.2) THEN
         outFilenames(1) = TRIM(denName) // '_f'
         outFilenames(2) = TRIM(denName) // '_g'
      ELSE
         outFilenames(1) = TRIM(denName) // '_f'
         outFilenames(2) = TRIM(denName) // '_A1'
         outFilenames(3) = TRIM(denName) // '_A2'
         outFilenames(4) = TRIM(denName) // '_A3'
         IF (polar) THEN
            outFilenames(5) = TRIM(denName) // '_Aabs'
            outFilenames(6) = TRIM(denName) // '_Atheta'
            outFilenames(7) = TRIM(denName) // '_Aphi'
         END IF
      END IF

      ! If xsf is specified we create input files for xcrysden
      IF (xsf) THEN
         DO i = 1, numOutFiles
            OPEN(nfile+i,file=TRIM(ADJUSTL(outFilenames(i)))//'.xsf',form='formatted')
            CALL xsf_WRITE_atoms(nfile+i,atoms,input%film,oneD%odi%d1,cell%amat)
         END DO
      END IF

      ! Loop over all plots
      DO nplo = 1, size(sliceplot%plot)
        twodim=sliceplot%plot(nplo)%twodim
        cartesian=sliceplot%plot(nplo)%cartesian
        grid=sliceplot%plot(nplo)%grid
        vec1=sliceplot%plot(nplo)%vec1
        vec2=sliceplot%plot(nplo)%vec2
        vec3=sliceplot%plot(nplo)%vec3
        zero=sliceplot%plot(nplo)%zero
        filename=sliceplot%plot(nplo)%filename
         IF (twodim.AND.ANY(grid(1:2)<1)) &
            CALL juDFT_error("Illegal grid size in plot",calledby="plot")
         IF (.NOT.twodim.AND.ANY(grid<1)) &
            CALL juDFT_error("Illegal grid size in plot",calledby="plot")
         IF (twodim) grid(3) = 1

         !calculate cartesian coordinates if needed
         IF (.NOT.cartesian) THEN
            vec1=matmul(cell%amat,vec1)
            vec2=matmul(cell%amat,vec2)
            vec3=matmul(cell%amat,vec3)
            zero=matmul(cell%amat,zero)
         END IF

         !Open the file
         IF (filename =="default") WRITE(filename,'(a,i2)') "plot",nplo

         IF (xsf) THEN
            DO i = 1, numOutFiles
               CALL xsf_WRITE_header(nfile+i,twodim,filename,vec1,vec2,vec3,zero,grid)
            END DO
         ELSE
            OPEN (nfile,file = TRIM(ADJUSTL(denName))//'_'//filename,form='formatted')
         END IF

         IF (twodim) THEN
            IF (numOutFiles.EQ.1) THEN
               WRITE(nfile,'(3a15)') 'x','y','f'
            ELSE IF (numOutFiles.EQ.2) THEN
               WRITE(nfile,'(4a15)') 'x','y','f','g'
            ELSE IF (numOutFiles.EQ.4) THEN
               WRITE(nfile,'(6a15)') 'x','y','f','A1','A2','A3'
            ELSE
               WRITE(nfile,'(9a15)') 'x','y','f','A1','A2','A3','|A|','theta','phi'
            END IF
         ELSE
            IF (numOutFiles.EQ.1) THEN
               WRITE(nfile,'(4a15)') 'x','y','z','f'
            ELSE IF (numOutFiles.EQ.2) THEN
               WRITE(nfile,'(5a15)') 'x','y','z','f','g'
            ELSE IF (numOutFiles.EQ.4) THEN
               WRITE(nfile,'(7a15)') 'x','y','z','f','A1','A2','A3'
            ELSE
               WRITE(nfile,'(10a15)') 'x','y','z','f','A1','A2','A3','|A|','theta','phi'
            END IF
         END IF

         !loop over all points
         DO iz = 0, grid(3)-1
            DO iy = 0, grid(2)-1
               DO ix = 0, grid(1)-1

                  point = zero + vec1*REAL(ix)/(grid(1)-1) +&
                                 vec2*REAL(iy)/(grid(2)-1)
                  IF (.NOT.twodim) point = point + vec3*REAL(iz)/(grid(3)-1)

                  ! Set region specific parameters for point

                  ! Get MT sphere for point if point is in MT sphere
                  CALL getMTSphere(input,cell,atoms,oneD,point,nt,na,pt)
                  IF (na.NE.0) THEN
                  ! In MT sphere
                     iv = 0
                     iflag = 1
                  ELSE IF (input%film.AND..NOT.oneD%odi%d1.AND.ABS(point(3))>=cell%z1) THEN
                     ! In vacuum in 2D system
                     iv = 1
                     iflag = 0
                     pt(:) = point(:)
                  ELSE IF ((oneD%odi%d1).AND.(SQRT((point(1))**2 + (point(2))**2)>=cell%z1)) THEN
                     ! In vacuum in 1D system
                     iv = 1
                     iflag = 0
                     pt(:) = point(:)
                  ELSE
                     ! In interstitial region
                     iv = 0
                     iflag = 2
                     pt(:) = point(:)
                  END IF

                  DO i = 1, numInDen
                     CALL outcdn(pt,nt,na,iv,iflag,1,potnorm,stars,&
                                 vacuum,sphhar,atoms,sym,cell,oneD,&
                                 den(i),xdnout(i))
                  END DO

                  IF (na.NE.0) THEN
                     IF (noco%l_ss) THEN
                        ! rotate magnetization "backward"
                        angss = DOT_PRODUCT(qssc,pt-atoms%pos(:,na))
                        help(1) = xdnout(2)
                        help(2) = xdnout(3)
                        xdnout(2) = +help(1)*COS(angss)+help(2)*SIN(angss)
                        xdnout(3) = -help(1)*SIN(angss)+help(2)*COS(angss)
                        ! xdnout(2)=0. ; xdnout(3)=0. ; xdnout(4)=0.
                     END IF
                  END IF

                  IF (noco%l_ss .AND. (.NOT. unwind)) THEN
                     ! rotate magnetization
                     angss = DOT_PRODUCT(qssc,point)
                     help(1) = xdnout(2)
                     help(2) = xdnout(3)
                     xdnout(2) = +help(1)*COS(angss) -help(2)*SIN(angss)
                     xdnout(3) = +help(1)*SIN(angss) +help(2)*COS(angss)
                  END IF

                  IF (polar) THEN
                     xdnout(5) = SQRT(ABS(xdnout(2)**2+xdnout(3)**2+xdnout(4)**2))
                     IF (xdnout(5)<eps) THEN
                        xdnout(5)= 0.0
                        xdnout(6)= -tpi_const
                        xdnout(7)= -tpi_const
                     ELSE
                        DO j = 1, 3
                           help(j) = xdnout(1+j)/xdnout(5)
                        END DO
                        IF (help(3)<0.5) THEN
                           xdnout(6)= ACOS(help(3))
                        ELSE
                           xdnout(6)= pi_const/2.0-ASIN(help(3))
                        END IF
                        IF (SQRT(ABS(help(1)**2+help(2)**2)) < eps) THEN
                           xdnout(7)= -tpi_const
                        ELSE
                           IF ( ABS(help(1)) > ABS(help(2)) ) THEN
                              xdnout(7)= ABS(ATAN(help(2)/help(1)))
                           ELSE
                              xdnout(7)= pi_const/2.0-ABS(ATAN(help(1)/help(2)))
                           END IF
                           IF (help(2)<0.0) THEN
                              xdnout(7)= -xdnout(7)
                           END IF
                           IF (help(1)<0.0) THEN
                              xdnout(7)= pi_const-xdnout(7)
                           END IF
                           phi0=0
                           DO WHILE (xdnout(7)-pi_const*phi0 > +pi_const)
                              xdnout(7)= xdnout(7)-tpi_const
                           END DO
                           DO WHILE (xdnout(7)-pi_const*phi0 < -pi_const)
                              xdnout(7)= xdnout(7)+tpi_const
                           END DO
                           IF (ABS(xdnout(2)-xdnout(3))<eps) THEN
                              IF (xdnout(2)>0) THEN
                                 xdnout(7)=pi_const/4.0
                              ELSE
                                 xdnout(7)=-3*pi_const/4.0
                              END IF
                           END IF
                        END IF
                     END IF
                     xdnout(6)= xdnout(6)/pi_const
                     xdnout(7)= xdnout(7)/pi_const
                  END IF ! (polar)
                  IF (xsf) THEN
                     DO i = 1, numOutFiles
                        WRITE(nfile+i,*) xdnout(i)
                     END DO
                  ELSE
                     WRITE(nfile,'(10e15.7)') point ,xdnout
                  END IF
               END DO !x-loop
            END DO !y-loop
         END DO !z-loop


         IF (xsf) THEN
            DO i = 1, numOutFiles
               CALL xsf_WRITE_endblock(nfile+i,twodim)
            END DO
         ELSE
            CLOSE(nfile)
         END IF

      END DO !nplot

      IF (xsf) THEN
         DO i = 1, numOutFiles
            CLOSE(nfile+i)
         END DO
      END IF

      DEALLOCATE(xdnout, outFilenames)

   END SUBROUTINE savxsf

   SUBROUTINE vectorplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                         noco, factor, score, potnorm, denmat, denName)

      ! Takes a spin-polarized t_potden variable, i.e. a 2D vector in MT-sphere/
      ! plane wave representation, splits it into two spinless ones, which are
      ! then passed on to the savxsf routine to get 2 .xsf files out.
      TYPE(t_sliceplot), INTENT(IN) :: sliceplot
      TYPE(t_stars),     INTENT(IN) :: stars
      TYPE(t_atoms),     INTENT(IN) :: atoms
      TYPE(t_sphhar),    INTENT(IN) :: sphhar
      TYPE(t_vacuum),    INTENT(IN) :: vacuum
      TYPE(t_input),     INTENT(IN) :: input
      TYPE(t_oned),      INTENT(IN) :: oneD
      TYPE(t_sym),       INTENT(IN) :: sym
      TYPE(t_cell),      INTENT(IN) :: cell
      TYPE(t_noco),      INTENT(IN) :: noco
      REAL,              INTENT(IN) :: factor
      LOGICAL,           INTENT(IN) :: score, potnorm
      TYPE(t_potden),    INTENT(IN) :: denmat
      CHARACTER(len=20), INTENT(IN) :: denName

      TYPE(t_potden)                :: cden, mden

      CALL vectorsplit(stars, atoms, sphhar, vacuum, input, factor, denmat, &
                       cden, mden)

      CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, &
                  score, potnorm, denName, cden, mden)

   END SUBROUTINE vectorplot

   SUBROUTINE matrixplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                         noco, factor, score, potnorm, denmat, denName)

      ! Takes a 2x2 t_potden variable, i.e. a sum of Pauli matrices in MT-
      ! sphere/ plane wave representation and splits it into four spinless ones,
      ! which are then passed on to the savxsf routine to get 4 .xsf files out.
      TYPE(t_sliceplot), INTENT(IN) :: sliceplot
      TYPE(t_stars),     INTENT(IN) :: stars
      TYPE(t_atoms),     INTENT(IN) :: atoms
      TYPE(t_sphhar),    INTENT(IN) :: sphhar
      TYPE(t_vacuum),    INTENT(IN) :: vacuum
      TYPE(t_input),     INTENT(IN) :: input
      TYPE(t_oneD),      INTENT(IN) :: oneD
      TYPE(t_sym),       INTENT(IN) :: sym
      TYPE(t_cell),      INTENT(IN) :: cell
      TYPE(t_noco),      INTENT(IN) :: noco
      REAL,              INTENT(IN) :: factor
      LOGICAL,           INTENT(IN) :: score, potnorm
      TYPE(t_potden),    INTENT(IN) :: denmat
      CHARACTER(len=20), INTENT(IN) :: denName

      TYPE(t_potden)                   :: cden, mxden, myden, mzden

      CALL matrixsplit(sym,stars, atoms, sphhar, vacuum, input, noco, factor, &
                       denmat, cden, mxden, myden, mzden)

      CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, &
                  score, potnorm, denName, cden, mxden, myden, mzden)

   END SUBROUTINE matrixplot

   SUBROUTINE procplot(stars, atoms, sphhar, sliceplot,vacuum, input, oneD, sym, cell, &
                       noco, denmat, plot_const)

      ! According to iplot, we process which exact plots we make after we assured
      ! that we do any. n-th digit (from the back) of iplot ==1 --> plot with
      ! identifier n is done.
      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_sphhar),    INTENT(IN)    :: sphhar
      TYPE(t_sliceplot), INTENT(IN)    :: sliceplot
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_oneD),      INTENT(IN)    :: oneD
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_potden),    INTENT(IN)    :: denmat
      INTEGER,           INTENT(IN)    :: plot_const

      INTEGER            :: i
      REAL               :: factor
      CHARACTER (len=20) :: denName
      LOGICAL            :: score, potnorm

      ! Plotting the input density matrix as n / n, m / n, mx, my, mz.
      ! Plot identifier: PLOT_INPDEN = 1
      ! Additive term for iplot: 2
      IF (plot_const.EQ.1) THEN
         factor = 1.0
         IF(sliceplot%slice) denName='slice'
         IF(.NOT.sliceplot%slice)denName = 'denIn'
         score = .FALSE.
         potnorm = .FALSE.
         IF (input%jspins.EQ.2) THEN
            IF (noco%l_noco) THEN
               CALL matrixplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            ELSE
               CALL vectorplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            END IF
         ELSE
            CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                        noco, score, potnorm, denName, denmat)
         END IF
      END IF

      ! Plotting the output density matrix as n / n, m / n, mx, my, mz.
      ! Plot identifier: PLOT_OUTDEN_Y_CORE = 2
      ! No core subtraction done!
      ! Additive term for iplot: 4
      IF (plot_const.EQ.2) THEN
         factor = 1.0
         denName = 'denOutWithCore'
         score = .FALSE.
         potnorm = .FALSE.
         IF (input%jspins.EQ.2) THEN
            IF (noco%l_noco) THEN
               CALL matrixplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            ELSE
               CALL vectorplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            END IF
         ELSE
            CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                        noco, score, potnorm, denName, denmat)
         END IF
      END IF

      ! Plotting the output density matrix as n / n, m / n, mx, my, mz.
      ! Plot identifier: PLOT_OUTDEN_N_CORE = 3
      ! Core subtraction done!
      ! Additive term for iplot: 8
      IF (plot_const.EQ.3) THEN
         factor = 1.0
         denName = 'denOutNOCore'
         score = .TRUE.
         potnorm = .FALSE.
         IF (input%jspins.EQ.2) THEN
            IF (noco%l_noco) THEN
               CALL matrixplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            ELSE
               CALL vectorplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            END IF
         ELSE
            CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                        noco, score, potnorm, denName, denmat)
         END IF
      END IF

      ! Plotting the mixed density matrix as n / n, m / n, mx, my, mz.
      ! Plot identifier: PLOT_MIXDEN_Y_CORE = 4
      ! No core subtraction done!
      ! Additive term for iplot: 16
      IF (plot_const.EQ.4) THEN
         factor = 1.0
         denName = 'denOutMixWithCore'
         score = .FALSE.
         potnorm = .FALSE.
         IF (input%jspins.EQ.2) THEN
            IF (noco%l_noco) THEN
               CALL matrixplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            ELSE
               CALL vectorplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            END IF
         ELSE
            CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                        noco, score, potnorm, denName, denmat)
         END IF
      END IF

      ! Plotting the mixed density matrix as n / n, m / n, mx, my, mz.
      ! Plot identifier: PLOT_MIXDEN_N_CORE = 5
      ! Core subtraction done!
      ! Additive term for iplot: 32
      IF (plot_const.EQ.5) THEN
         factor = 1.0
         denName = 'denOutMixNoCore'
         score = .TRUE.
         potnorm = .FALSE.
         IF (input%jspins.EQ.2) THEN
            IF (noco%l_noco) THEN
               CALL matrixplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            ELSE
               CALL vectorplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            END IF
         ELSE
            CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                        noco, score, potnorm, denName, denmat)
         END IF
      END IF


      ! Plotting the total potential as vTot / v_eff, B_eff / v_eff,
      ! B_xc_1, B_xc_2, B_xc_3.
      ! Plot identifier: PLOT_POT_TOT = 7
      ! No core subtraction done!
      ! Additive term for iplot: 128
      IF (plot_const.EQ.7) THEN
         factor = 2.0
         denName = 'vTot'
         score = .FALSE.
         potnorm = .TRUE.
         IF (input%jspins.EQ.2) THEN
            IF (noco%l_noco) THEN
               CALL matrixplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            ELSE
               CALL vectorplot(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, &
                               cell, noco, factor, score, potnorm, denmat, &
                               denName)
            END IF
         ELSE
            CALL savxsf(sliceplot,stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                        noco, score, potnorm, denName, denmat)
         END IF
      END IF

   END SUBROUTINE procplot

   SUBROUTINE makeplots(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                        noco, denmat, plot_const, sliceplot)

      ! Checks, based on the iplot switch that is given in the input, whether or
      ! not plots should be made. Before the plot command is processed, we check
      ! whether the plot_inp is there and no oldform is given. If that is not the
      ! case, we create a plot_inp.
      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_sphhar),    INTENT(IN)    :: sphhar
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_oneD),      INTENT(IN)    :: oneD
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_potden),    INTENT(IN)    :: denmat
      INTEGER,           INTENT(IN)    :: plot_const
      TYPE(t_sliceplot), INTENT(IN)    :: sliceplot

      LOGICAL :: allowplot

      ! The check is done via bitwise operations. If the i-th position of iplot
      ! in binary representation (2^n == n-th position) has a 1, the correspon-
      ! ding plot with number 2^n is plotted.
      !
      ! E.g.: If the plots with identifying constants 1,2 and 4 are to be plotted
      ! and none else, iplot would need to be 2^1 + 2^2 + 2^3 = 2 + 4 + 8 = 14.
      ! iplot=1 or any odd number will *always* plot all possible options.

      CALL timestart("Plotting iplot plots")

      allowplot=BTEST(sliceplot%iplot,plot_const).OR.(MODULO(sliceplot%iplot,2).EQ.1)
      IF (allowplot) THEN
         CALL checkplotinp()
         CALL  procplot(stars, atoms, sphhar,sliceplot, vacuum, input, oneD, sym, cell, &
                        noco, denmat, plot_const)
      END IF
      CALL timestop("Plotting iplot plots")
   END SUBROUTINE makeplots

   SUBROUTINE getMTSphere(input, cell, atoms, oneD, point, iType, iAtom, pt)

      ! Old subroutine originally from plotdop.f90, which is needed in savxsf.

      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_cell),  INTENT(IN)    :: cell
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_oneD),  INTENT(IN)    :: oneD

      INTEGER,       INTENT(OUT)   :: iType, iAtom
      REAL,          INTENT(OUT)   :: pt(3)
      REAL,          INTENT(IN)    :: point(3)

      INTEGER                      :: ii1, ii2, ii3, i1, i2, i3, nq
      REAL                         :: s

      ii1 = 3
      ii2 = 3
      ii3 = 3
      IF (input%film .AND. .NOT.oneD%odi%d1) ii3 = 0
      IF (oneD%odi%d1) THEN
         ii1 = 0
         ii2 = 0
      END IF

      DO i1 = -ii1, ii1
         DO i2 = -ii2, ii2
            DO i3 = -ii3, ii3
               pt = point+MATMUL(cell%amat,(/i1,i2,i3/))
               iAtom = 0
               DO iType = 1, atoms%ntype
                  DO nq = 1, atoms%neq(iType)
                     iAtom = iAtom + 1
                     s = SQRT(DOT_PRODUCT(atoms%pos(:,iAtom)-pt,atoms%pos(:,iAtom)-pt))
                     IF (s<atoms%rmsh(atoms%jri(iType),iType)) THEN
                        ! Return with the current iType, iAtom, pt
                        RETURN
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO !i1

      ! If no MT sphere encloses the point return 0 for iType, iAtom
      iType = 0
      iAtom = 0
      pt(:) = point(:)

   END SUBROUTINE getMTSphere

END MODULE m_plot
