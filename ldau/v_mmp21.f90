MODULE m_vmmp21
!     ************************************************************
!     This subroutine calculates the potential matrix v^{s}_{m,m'}
!     for a given atom 'iType' and l-quantum number 'l'. The l,u,j's for
!     all atoms are stored in lda_u, the density matrix is ns_mmp,
!     and the e-e- interaction potential is u(m1,m2,m3,m4,iType).
!     For details see Eq.(16) of Shick et al. PRB 60, 10765 (1999)
!
!     Additionally, the total energy contribution of LDA+U (Eq.24)
!     is calculated (e_ldau).
!     Part of the LDA+U package                   G.B., Oct. 2000
!     ************************************************************
      USE m_types
      USE m_constants

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE v_mmp_21(atoms,n_mmp21,u,v_mmp21,e_ldau)

      TYPE(t_atoms),    INTENT(IN)     :: atoms
      REAL,             INTENT(IN)     :: u(-lmaxU_const:,-lmaxU_const:,-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX,          INTENT(IN)     :: n_mmp21(-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX,          INTENT(INOUT)  :: v_mmp21(-lmaxU_const:,-lmaxU_const:,:)
      REAL,             INTENT(INOUT)  :: e_ldau

      INTEGER :: l,m,mp,p,q,iType,i_u
      REAL    :: e_dc, e_off
      !
      ! Loop over atoms
      !
      e_off = 0.0
      DO i_u = 1, atoms%n_u+atoms%n_hia
         iType = atoms%lda_u(i_u)%atomType
         l     = atoms%lda_u(i_u)%l
         !
         !-------------------------------------+
         !  offd    --                    offd |
         ! V     =  >  (  <m,p|V|q,m'> ) n     |
         !  m,m'    --                    p,q  |
         !        p,q                          |
         !-------------------------------------+
         v_mmp21(:,:,i_u) = cmplx_0

         DO m = -l,l
            DO mp =-l,l
               DO p = -l,l
                  DO q = -l,l
                     v_mmp21(m,mp,i_u) = v_mmp21(m,mp,i_u) + &
                                          u(m,p,q,mp,i_u) * n_mmp21(p,q,i_u)
                  ENDDO
               ENDDO
            ENDDO ! m' loop
         ENDDO ! m  loop

         !----------------------------------------------------------------------+
         !              s                                                       !
         !  ee      1  ---   s        s                     1        s  1       !
         ! E  (n) = -  >    n      ( V     + d     ( U (n - -) - J (n - -) ))   !
         !          2  ---   m,m'     m,m'    m,m'          2           2       !
         !             m,m'                                                     !
         !----------------------------------------------------------------------+

         e_dc = 0.0
         DO m = -l,l
            DO mp =-l,l
               DO p = -l,l
                  DO q = -l,l
                     e_dc = e_dc + u(m,p,q,mp,i_u) * (&
                                    n_mmp21(m,mp,i_u)*conjg(n_mmp21(q,p,i_u)) +&
                                    conjg(n_mmp21(mp,m,i_u))*n_mmp21(p,q,i_u) )
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         e_off = e_off + e_dc

      ENDDO ! loop over atoms

      e_ldau = e_ldau + e_off / 2

   END SUBROUTINE v_mmp_21
END MODULE m_vmmp21
