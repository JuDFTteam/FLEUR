!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_gf_current_single
    implicit none
    contains
    subroutine gf_current_single(layers,lapw,cell,sym,mpi,bkpts,nk,en,jspin,sigmaL)
        use m_gf_types
        USE m_gf_writetrans,ONLY:writetrans
        USE m_gf_embedding
        implicit none
        type(t_layers),intent(in):: layers
        type(t_lapw),intent(in) :: lapw
        type(t_cell),intent(in) :: cell
        type(t_sym),intent(in)  :: sym
        type(t_mpi),intent(in)  :: mpi
        real,intent(in)         :: bkpts(:,:)
        integer,intent(in)      :: nk,en,jspin
        complex,intent(in)      :: sigmaL(:,:)

        COMPLEX :: sigmaR(size(sigmaL,1),size(sigmaL,2))
        real    :: j(2)
        logical :: l_noco

        l_noco=.not.(lapw_gf%nv2_tot==lapw_gf%nv2(1))

        CALL gf_getemb2(sigmaR,2,layers%num_layers,en,nk,jspin,lapw)
        CALL gf_landauer1plane(l_noco,                                 &
     &        sigmaL,                                                       &
     &        sigmaR,                                                       &
     &        j,0.0)
       CALL writetrans(en,nk,jspin,bkpts,sym,cell,2,j,mpi)


    end subroutine gf_current_single


    SUBROUTINE gf_Landauer1Plane(l_noco,                              &
     &                             g1,g2,                               &
     &                             j,imaginary)
!c*********************************************************************
!c     subroutine to calculate the current from two embedding
!c     potentials on the same plane
!c
!c                                      Daniel Wortmann
!c*********************************************************************
      USE m_gf_math
      IMPLICIT NONE
      LOGICAL,INTENT(IN)  :: l_noco
      COMPLEX,INTENT(IN)  :: g1(:,:),g2(:,:)
      REAL,INTENT(OUT)    :: j(2)
      real,intent(in)     :: imaginary

      INTEGER             :: n
      COMPLEX             :: G(SIZE(g1,1),SIZE(g1,1))
      COMPLEX             :: A(SIZE(g1,1),SIZE(g1,1)),B(SIZE(g1,1)      &
     &     ,SIZE(g1,1))


      IF(.FALSE.)THEN
         G = mat_inverse(G1+G2)
         j(1) = 2.0*REAL(trace(MATMUL(MATMUL(imag2d(G1),G),             &
     &     matmul(imag2d(g2)                                            &
     &     ,TRANSPOSE(CONJG(G))))))

      ELSE
        ! longer version
        ! Re(G1 G12 G2 G21^* - G1 G12 G2^* G21^*)
        G = mat_inverse(G1+G2+cmplx(0.0,imaginary))
        A=matmul(matmul(g1,g),g2)
        B=matmul(matmul(g1,g),transpose(conjg(g2)))
        if (abs(imaginary)>epsilon(1.0)) G = mat_inverse(G1+G2-cmplx(0.0,imaginary))
        A=matmul(A,transpose(conjg(g)))
        B=matmul(B,transpose(conjg(g)))

        j(1) =-2.*REAL(trace(A-B))
      ENDIF

      IF (l_noco) THEN
         n = SIZE(g1,1)/2
         j(2) = 2.0*REAL(trace(MATMUL(MATMUL(imag2d(G1(:n,:n)),G(:n,:n))&
     &        ,MATMUL(imag2d(g2(:n,:n)),TRANSPOSE(CONJG(G(:n,:n)))))))
      ELSE
         j(2) = 0.0
      ENDIF

      END SUBROUTINE gf_landauer1plane
end
