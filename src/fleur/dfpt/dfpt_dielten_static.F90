!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_dielten_static
    use m_juDFT
    use m_constants
    use m_xmlOutput
    use m_types_fleurinput
    use m_dfpt_NAC, only: read_bec, read_dielten
    use m_dfpt_dynmat_io, only: read_dynmat_q

    implicit none
    private
    public :: dfpt_dielten_static

contains

    subroutine dfpt_dielten_static(fi)
        !! Static dielectric tensor eps0 = eps_inf + 4pi/Omega sum_m p_m p_m^T / omega_m^2 from the Gamma
        !! eigenpairs in dynmat.hdf (phonon run with l_polar=F), born_eff_charge and diel_tensor.
        !! p_m(beta) = sum_{kappa,alpha} Z*_kappa(alpha,beta) e_m(kappa alpha) / sqrt(M_kappa)

        type(t_fleurinput), intent(in) :: fi

        real,    allocatable :: borneffcharge(:,:,:), eigenVals(:), weightTrans(:), mass(:)
        complex, allocatable :: dynMat(:,:), eigenVecs(:,:)
        logical, allocatable :: l_skip(:)

        real    :: dieltenInf(3,3), dieltenStatic(3,3)
        complex :: pMode(3), transProj(3)
        integer :: nModes, iQ, iMode, iAcoustic, iAtom, iDir, iDir1
        character(len=20) :: attributes(1)

        nModes = 3*fi%atoms%nat

        allocate(borneffcharge(fi%atoms%ntype,3,3))
        call read_bec(fi, borneffcharge)
        call read_dielten(fi, dieltenInf)

        ! Gamma is in the phonon q-point list, checked in precheck_dfpt
        iQ = findloc(norm2(fi%dfpt%qvec%bk(:,:fi%dfpt%qvec%nkpt),dim=1) < 1e-8, .true., dim=1)
        allocate(dynMat(nModes,nModes), eigenVals(nModes), eigenVecs(nModes,nModes))
        call read_dynmat_q("dynmat.hdf", iQ, dynMat, eigVals=eigenVals, eigVecs=eigenVecs)

        allocate(mass(fi%atoms%nat))
        do iAtom = 1, fi%atoms%nat
            mass(iAtom) = atomicMasses_const(fi%atoms%nz(fi%atoms%itype(iAtom)))*massInElectronMasses
        end do

        ! Acoustic modes: the three modes with the largest overlap with a rigid translation
        allocate(weightTrans(nModes), l_skip(nModes))
        do iMode = 1, nModes
            transProj = cmplx(0.0,0.0)
            do iAtom = 1, fi%atoms%nat
                transProj(:) = transProj(:) + sqrt(mass(iAtom))*eigenVecs(3*(iAtom-1)+1:3*iAtom,iMode)
            end do
            weightTrans(iMode) = sum(abs(transProj)**2)/sum(mass)
        end do
        l_skip = .false.
        do iAcoustic = 1, 3
            l_skip(maxloc(weightTrans, dim=1, mask=.not.l_skip)) = .true.
        end do

        dieltenStatic = dieltenInf
        do iMode = 1, nModes
            if (l_skip(iMode)) cycle
            if (eigenVals(iMode) <= 0.0) then
                call juDFT_warn("Unstable optical mode at Gamma is left out of the static dielectric tensor.", calledby="dfpt_dielten_static")
                cycle
            end if
            pMode = cmplx(0.0,0.0)
            do iAtom = 1, fi%atoms%nat
                do iDir = 1, 3
                    pMode(:) = pMode(:) + borneffcharge(fi%atoms%itype(iAtom),iDir,:)*eigenVecs(3*(iAtom-1)+iDir,iMode)/sqrt(mass(iAtom))
                end do
            end do
            do iDir1 = 1, 3
                do iDir = 1, 3
                    dieltenStatic(iDir,iDir1) = dieltenStatic(iDir,iDir1) &
                                              + fpi_const/fi%cell%omtil*real(pMode(iDir)*conjg(pMode(iDir1)))/eigenVals(iMode)
                end do
            end do
        end do

        write(*,*) '-------------------------'
        write(*,*) "Static Dielectric tensor"
        open(110, file="diel_tensor_static", status='replace', action='write', form='formatted')
        do iDir = 1, 3
            write(110,'(2(es16.8,1x),es16.8)') dieltenStatic(iDir,:)
            write(*,'(2(es16.8,1x),es16.8)') dieltenStatic(iDir,:)
        end do
        close(110)
        write(*,*) '-------------------------'

        CALL openXMLElementNoAttributes('DFPTEfield')
        do iDir = 1, 3
            attributes = ''
            WRITE(attributes(1),'(i0)') iDir
            CALL writeXMLElementPoly('dieltensorStatic',(/'iDir'/), attributes, dieltenStatic(iDir,:))
        end do
        CALL closeXMLElement('DFPTEfield')

    end subroutine dfpt_dielten_static

end module m_dfpt_dielten_static
