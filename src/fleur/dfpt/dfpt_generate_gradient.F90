!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


module m_dfpt_generate_gradient 

    use m_juDFT
    use m_constants
    use m_types


    implicit none 

contains 

    subroutine dfpt_generate_gradient(sternheimerJob,fi,fmpi,sphhar,hybdat,xcpot,nococonv,stars,rho,vTot,grRho3,grVtot3,grVc3,grVext3,grgrVext3x3)
        
        use m_dfpt_vgen
        use m_vgen_coulomb
        use m_dfpt_potdenLocal
        use m_dfpt_gradient 
        use m_grdchlh
        

        type(t_sternheimerjob),intent(in) :: sternheimerJob
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

        complex, allocatable :: grrhodummy(:, :, :, :, :) 
        real    :: dr_re(fi%vacuum%nmzd), dr_im(fi%vacuum%nmzd), drr_dummy(fi%vacuum%nmzd)

        
        call rho_tmp%copyPotDen(rho)

        ALLOCATE(grrhodummy(fi%atoms%jmtd, (fi%atoms%lmaxd+1)**2, fi%atoms%nat, SIZE(rho%mt,4), 3))

        ! For the gradient of the external Potential we need a higher cutoff in the expansion, in order to confine the multipole moments in the MTs
        ! This is crucial for the Film-Mode and will be visible in the z-Eigenmodes.  
        call create_typesLocal(fi,fmpi,fi%sym,fi%cell,fi%input,sphhar, fi%vacuum , fi%noco ,starsLocal, grVext3(1),atomsLocal,l_dfpt=.TRUE.)
        call grVext3(2)%copyPotDen(grVext3(1))
        call grVext3(3)%copyPotDen(grVext3(1))

        ! initialize all potentials
        ! currently we do this by copy, is this faster than just init?

        call potdummy%copyPotDen(rho)
        call potdummy%resetpotden()
        call potdummyLocal%copyPotDen(grVext3(1))
        call potdummyLocal%resetpotden()


        do iDir = 1 , 3 
            call grRho3(iDir)%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,rho%potdenType,l_dfpt=.TRUE.)
            do iDir2 = 1 , 3 
                call grgrVext3x3(iDir2,iDir)%copyPotDen(grVext3(1))
                call grgrVext3x3(iDir2,iDir)%resetpotden()
            end do 
            call grVtot3(iDir)%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,vTot%potdenType,l_dfpt=.TRUE.)
            call grVc3(iDir)%init(stars,fi%atoms,sphhar,fi%vacuum,fi%noco,fi%input%jspins,vTot%potdenType,l_dfpt=.TRUE.)
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
            call vgen_coulomb(1, fmpi, fi%input, fi%field, fi%vacuum, fi%sym, starsLocal, fi%cell, &
                         & sphhar, atomsLocal, .FALSE., potdummyLocal, grVext3(iDir), sternheimerJob=sternheimerJob,dfpt=fi%dfpt, &
                         & dfptden0=potdummyLocal,stars2=starsLocal,iDtype=0,iDir=iDir)
        end do 

        ! Coulomb/Effective potential gradients
        do iDir = 1, 3
            call dfpt_vgen(sternheimerJob,hybdat, fi%field, fi%input, xcpot, fi%atoms, sphhar, stars, fi%vacuum, fi%sym, &
                            fi%dfpt, fi%cell, fmpi, fi%noco, nococonv, rho_tmp, vTot, &
                            stars, grVtot3(iDir), .TRUE., grRho3(iDir), 0, iDir, [0,0])
            call dfpt_vgen(sternheimerJob,hybdat, fi%field, fi%input, xcpot, fi%atoms, sphhar, stars, fi%vacuum, fi%sym, &
                            fi%dfpt, fi%cell, fmpi, fi%noco, nococonv, rho_tmp, vTot, &
                            stars, grVc3(iDir), .FALSE., grRho3(iDir), 0, iDir, [0,0])
        end do 


        ! Hessian of Potential for external Potential 

        do iDir2 = 1, 3
            do iDir = 1, 3
                call potdummyLocal%resetPotDen()
                call vgen_coulomb(1, fmpi, fi%input, fi%field, fi%vacuum, fi%sym, starsLocal, fi%cell, &
                            & sphhar, atomsLocal, .TRUE., potdummyLocal, grgrVext3x3(iDir2,iDir), sternheimerJob=sternheimerJob,dfpt=fi%dfpt, &
                            & dfptden0=potdummyLocal,stars2=starsLocal,iDtype=0,iDir=iDir,iDir2=iDir2)
            end do 
        end do 




              
      CALL grRho3(1)%distribute(fmpi%mpi_comm)
      CALL grRho3(2)%distribute(fmpi%mpi_comm)
      CALL grRho3(3)%distribute(fmpi%mpi_comm)
      CALL grVext3(1)%distribute(fmpi%mpi_comm)
      CALL grVext3(2)%distribute(fmpi%mpi_comm)
      CALL grVext3(3)%distribute(fmpi%mpi_comm)

    end subroutine dfpt_generate_gradient 

end module m_dfpt_generate_gradient