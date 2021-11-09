module m_multipoleExpansion

    use m_types
    use m_uj2f
    use m_wigner
    use m_factorials
    use m_xmlOutput

    implicit none

    contains

    !Multipole expansion of the DFT+U potential see Phys. Rev. B 80, 035121 (2009).

    subroutine dftuMultipoleExp(atoms, input, den)

        type(t_atoms), intent(in)  :: atoms
        type(t_input), intent(in)  :: input
        type(t_potden), intent(in) :: den

        real :: e_u
        real :: f(0:6)
        integer :: i_u,l,k,p,r,t,atomtype
        real :: e_hartree(0:2*lmaxU_const,0:1,0:2*lmaxU_const+1)
        real :: e_xc(0:2*lmaxU_const,0:1,0:2*lmaxU_const+1)
        real :: tensor_square(0:2*lmaxU_const,0:1,0:2*lmaxU_const+1)
        complex :: w(0:2*lmaxU_const,0:1,0:2*lmaxU_const+1,-(lmaxU_const+1)**2:(lmaxU_const+1)**2)
        character(len=20) :: attributes(7)

        call timestart('DFT+U: Multipole Expansion')
        call openXMLElementNoAttributes('dftUMultipoleExpansion')

        do i_u = 1, atoms%n_u
            atomtype = atoms%lda_u(i_u)%atomType
            l = atoms%lda_u(i_u)%l

            call uj2f(input%jspins, atoms%lda_u(i_u), f)

            e_hartree = 0.0
            e_xc = 0.0
            tensor_square = 0.0
            w = 0.0
            e_u = 0.0
            write(oUnit,*)
            write(oUnit,'("Multipole expansion atom type: ", I4," Orbital: ", I2)') atomtype, l
            !Loop over tensor indices
            do k = 0, 2*l
                do p = 0,1
                    do r = abs(k-p), k+p
                        !1. Calculate the 3-index tensor
                        w(k,p,r,:) = w_tensor(l,input%jspins,k,p,r,den%mmpMat(:,:,i_u,:))

                        !2. Calculate the sqare of the tensor
                        tensor_square(k,p,r) = cmplx_0
                        do t = -r, r
                            tensor_square(k,p,r) = tensor_square(k,p,r) + (-1)**t * (w(k,p,r,-t)*w(k,p,r,t))
                        enddo
                        !3. Calculate the contributions to the hartree and exchange energy
                        e_hartree(k,p,r) = hartree_coeff(l,k,p,r,f) * tensor_square(k,p,r)
                        e_xc(k,p,r) = xc_coeff(l,k,p,r,f) * tensor_square(k,p,r)
                        e_u = e_u + e_hartree(k,p,r) + e_xc(k,p,r)

                    enddo
                enddo
            enddo

            !4. Output
            attributes = ''
            write (attributes(1), '(i0)') atomType
            write (attributes(2), '(i0)') l
            write (attributes(3), '(f14.8)') e_u
            write (attributes(4), '(A)') 'Htr'
            call openXMLElementForm('dftUMultipoleExpansionfor', (/'atomType  ', 'l         ', 'dftUEnergy', 'units     '/), &
                                    attributes(:4), reshape((/8, 1, 10,5, 6, 1, 14,3/), (/4, 2/)))
            do k = 0, 2*l
                do p = 0,1
                    do r = abs(k-p), k+p
                        !out
                        write(oUnit,'(tr5,"k = ",I2," p = ",I2," r = ",I2)') k,p,r
                        write(oUnit,'(tr5,"w_krp.w_krp = ",f14.8," Hartree = ",f14.8, " Exchange = ",f14.8)') tensor_square(k,p,r), e_hartree(k,p,r), e_xc(k,p,r)
                        do t = -r, r
                            write(oUnit,'(tr9,"t =",I2,"; ",2f14.8)') t, w(k,p,r,t)
                        enddo
                        write(oUnit,*)

                        !out.xml
                        attributes = ''
                        write (attributes(1), '(i0)') k
                        write (attributes(2), '(i0)') p
                        write (attributes(3), '(i0)') r
                        write (attributes(4), '(f14.8)') tensor_square(k,p,r)
                        write (attributes(5), '(f14.8)') e_hartree(k,p,r)
                        write (attributes(6), '(f14.8)') e_xc(k,p,r)
                        write (attributes(7), '(A)') 'Htr'
                        call openXMLElementForm('multipoleTensor', (/'k        ', 'p        ', 'r        ', 'w_square ', 'e_hartree', 'e_xc     ', 'units    '/), &
                                                attributes, reshape((/1,1,1,8,9,4,5,2, 1, 2, 14, 14, 14,3/), (/7, 2/)))
                        do t = -r, r
                            attributes = ''
                            write (attributes(1), '(i0)') t
                            call writeXMLElementPoly('component', (/'t'/), &
                                                        attributes(:1), contentList=[real(w(k,p,r,t)), aimag(w(k,p,r,t))])
                        enddo

                        call closeXMLElement('multipoleTensor')
                    
                    enddo
                enddo
            enddo
            call closeXMLElement('dftUMultipoleExpansionfor')

            write(oUnit,'("Norm of tensor components = ",f14.8)') sqrt(sum(tensor_square(:,:,:)**2))
            write(oUnit,'("DFT+U energy correction = ",f14.8,/)') e_u
        enddo

        call closeXMLElement('dftUMultipoleExpansion')
        call timestop('DFT+U: Multipole Expansion')

    end subroutine

    function w_tensor(l,jspins,k,p,r,denmat)

        !Calculate the 3-index tensor of the DFT+U density matrix
        ! Phys. Rev. B 80, 035121 (2009) Eqs. 20-22, 24,25,27,28

        integer, intent(in) :: l,jspins
        integer, intent(in) :: k,p,r
        complex, intent(in) :: denmat(-lmaxU_const:,-lmaxU_const:,:)

        complex w_tensor(-(lmaxU_const+1)**2:(lmaxU_const+1)**2)

        integer :: x,y,t,ispin,jspin,spin_index,m,mp
        real :: nlk,nsp,sa,sb,tensor
        complex :: nkpr, prefactor

        nlk = factorial(2*l)/sqrt(factorial(2*l-k)*factorial(2*l+k+1))
        nsp = 1/sqrt(factorial(p+2))
        nkpr = three_tensor_norm(k,p,r)

        prefactor = (-1)**(k+p)/(nlk*nsp*nkpr)
        w_tensor = cmplx_0
        do t = -r, r
            do x = -k,k
                do y = -p,p
                    do ispin = 1,jspins
                        do jspin = merge(1,ispin,size(denmat,3)==3),merge(jspins,ispin,size(denmat,3)==3)
                            sa = ispin-1.5 
                            sb = jspin-1.5 
                            if(ispin==jspin) then
                                spin_index = ispin
                            else if(ispin>jspin) then
                                spin_index = 3
                            else
                                spin_index = 4
                            endif

                            tensor = (-1)**(x+y) * wigner3j(k,r,p, -x, t, -y)
                            tensor = tensor * (-1)**(jspin-1) * wigner3j(0.5, real(p), 0.5, -sb, real(y), sa)
                            if(abs(tensor)<1e-12) cycle

                            do m= -l, l
                                do mp = -l,l
                                    tensor = tensor * (-1)**(l-mp) * wigner3j(l, k, l, -mp, x, m)
                                    if (spin_index<4) then 
                                        w_tensor(t) = w_tensor(t) + prefactor * tensor * denmat(m,mp,spin_index)
                                    else
                                        w_tensor(t) = w_tensor(t) + prefactor * tensor * conjg(denmat(mp,m,spin_index))
                                    endif
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end function

    real function hartree_coeff(l,k,p,r,f)

        !Contributions to the hartree energy Eq. 23

        integer, intent(in) :: l,k,p,r
        real,    intent(in) :: f(0:) !Slater parameter
    
        real :: nlk

        hartree_coeff = 0.0
        if (mod(k,2) == 1) return
        if (p /= 0) return
        nlk = factorial(2*l)/sqrt(factorial(2*l-k)*factorial(2*l+k+1))

        hartree_coeff = (2*l+1)**2/2.0 * nlk**2 &
                       * wigner3j(l, k, l, 0, 0, 0)**2 * f(k)
        
    end function

    real function xc_coeff(l,k,p,r,f)

        !Contributions to the XC energy Eq. 30-32

        integer, intent(in) :: l,k,p,r
        real,    intent(in) :: f(0:) !Slater parameter

        integer kp
        real :: nlk,prefactor

        nlk = factorial(2*l)/sqrt(factorial(2*l-k)*factorial(2*l+k+1))
        prefactor = (-1)**(k+1) * abs(three_tensor_norm(k,p,r))**2 * &
                    nlk**2 * ((2*l+1)**2*(2*k+1)*(2*r+1))/4.0
        xc_coeff = 0.0
        do kp = 0, 4*l, 2
            xc_coeff = xc_coeff + prefactor * f(kp/2) &
                                 * wigner3j(l, kp/2, l, 0, 0, 0)**2 &
                                 * wigner6j(l, l, k, l, l, kp/2)
        enddo

    end function

    complex function three_tensor_norm(k,p,r)

        !Normalization factor for the three index tensor Eq. 28

        integer, intent(in) :: k,p,r

        integer :: g

        g = k + p + r

        three_tensor_norm = ImagUnit**g * sqrt(factorial(g-2*k)) &
                            * sqrt(factorial(g-2*p)) &
                            * sqrt(factorial(g-2*r)) & 
                            / sqrt(factorial(g+1)) &
                        * doubleFactorial(g) &
                        / ( &
                              doubleFactorial(g-2*k) &
                            * doubleFactorial(g-2*p) &
                            * doubleFactorial(g-2*r) &
                            )

    end function

end module m_multipoleExpansion