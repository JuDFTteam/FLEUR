!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_denmat_to_full
    use m_judft
    implicit none
    
contains

    subroutine denmatrix_to_full_density(itype,input,sphhar,atoms,sym,radfun,denmat,rho,rhoIm,moments)
    !! Current situation:
    !!
    !! This routine calculates density contributions
    !! $$\rho_{L}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha}(r)=
    !! \sum_{l',l,\lambda',\lambda,s}d_{l',l,L,\lambda',\lambda}^{\sigma_{\alpha}',\sigma_{\alpha},\alpha}
    !! u_{l',\lambda',s}^{\sigma_{\alpha}',\alpha}(r)u_{l,\lambda,s}^{\sigma_{\alpha},\alpha}(r)$$
    !! \(s\) is the index for the big/small components yielded by the
    !! scalar-relativistic Schrödinger equation.
    
          use m_types
          use m_types_denmatrix
          use m_types_radfun
          
          implicit none
          integer,intent(in)             :: itype
          type(t_input), intent(IN)      :: input
          type(t_sphhar), intent(IN)     :: sphhar
          type(t_atoms), intent(IN)      :: atoms
          type(t_sym), intent(IN)        :: sym
          type(t_denmatrix), intent(IN)  :: denmat(:,:)
          type(t_radfun),intent(in)      :: radfun
          real,intent(inout)             :: rho(:,0:,:,:)
    
          
          real,intent(inout),optional    :: rhoIm(:,0:)
          type(t_moments), optional, intent(INOUT) :: moments
    
          
          integer :: lh, l, lp,llp, j, i, ii,ispinpr,ispin,spin
          complex :: cs
          integer,parameter:: lcf=3
    
        DO ispinpr=lbound(denmat,1),ubound(denmat,1)
            DO ispin=lbound(denmat,2),ubound(denmat,2)
                spin=merge(ispin,3,ispin==ispinpr)
                    do lh = 0, sphhar%nlh(sym%ntypsy(atoms%firstAtom(itype)))
                        do l = 0, atoms%lmax(itype)
                            do lp = 0, merge(l, atoms%lmax(itype), present(moments))
                            llp = (l*(l + 1))/2 + lp
                            if (.not. present(moments)) llp = lp*(atoms%lmax(itype) + 1) + l
                            if (lh > 0 .and. atoms%l_outputCFpot(itype) .and. atoms%l_outputCFremove4f(itype) &
                                .and. (l == lcf .and. lp == lcf)) cycle !Exclude non-spherical contributions for CF
                
                            do j = 1, atoms%jri(itype)
                                cs = 0.0
                                do i = 1, radfun%n_r(lp) !Loop over radial functions
                                    do ii = 1, radfun%n_r(l)
                                        cs = cs + denmat(ispinpr, ispin)%mat(i, ii, lp, l, lh) &   !density matrix
                                            *(radfun%r(i, j, 1, lp, ispin)*radfun%R(ii, j, 1, l, ispinpr) &  !large components
                                            + radfun%R(i, j, 2, lp, ispin)*radfun%R(ii,j, 2, l, ispinpr))    !small components
                                    end do
                                end do
                                rho(j, lh,itype,spin) = rho(j, lh,itype,spin) + real(cs)/atoms%neq(itype)
                                if (spin==3) rho(j, lh,itype,4) = rho(j, lh,itype,4) + aimag(cs)/atoms%neq(itype) !Store imaginary part as 4th spin
                                if (present(rhoIm)) rhoIm(j, lh) = rhoIm(j, lh) + aimag(cs)/atoms%neq(itype)
                                if ((l <= input%lResMax) .and. (lp <= input%lResMax) .and. ispin == ispinpr .and. present(moments)) &
                                moments%rhoLRes(j, lh, llp, itype, ispin) = moments%rhoLRes(j, lh, llp, itype, ispin) + real(cs)/atoms%neq(itype)
                
                            end do
                            end do
                        end do
                
                ENDDO
            enddo
        enddo
       end subroutine denmatrix_to_full_density
    end module