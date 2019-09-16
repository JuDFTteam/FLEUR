MODULE m_vmmp21
!     ************************************************************
!     This subroutine calculates the potential matrix v^{s}_{m,m'}
!     for a given atom 'n' and l-quantum number 'l'. The l,u,j's for
!     all atoms are stored in lda_u, the density matrix is ns_mmp,
!     and the e-e- interaction potential is u(m1,m2,m3,m4,n).
!     For details see Eq.(16) of Shick et al. PRB 60, 10765 (1999)
!
!     Additionally, the total energy contribution of LDA+U (Eq.24)
!     is calculated (e_ldau).
!     Part of the LDA+U package                   G.B., Oct. 2000
!     ************************************************************
      CONTAINS
      SUBROUTINE v_mmp_21(u_in,n_u,irank,n_mmp21,u,f0,f2,v_mmp21,e_off)

      USE m_types
      USE m_constants

      IMPLICIT NONE
!
! ..  Arguments ..
      INTEGER, INTENT (IN) :: n_u,irank
      REAL,    INTENT (IN) :: u(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,&
                                -lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_u)
      REAL,    INTENT (IN) :: f0(n_u),f2(n_u)
      REAL,    INTENT (OUT):: e_off 
      TYPE (t_utype), INTENT (IN) :: u_in(n_u)

      COMPLEX           :: n_mmp21(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_u)
      COMPLEX,INTENT(OUT)::v_mmp21(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_u)
!
! ..  Local Variables ..
      INTEGER n,ispin,jspin,l,m,mp,p,q,iType,i_u
      REAL rho_tot,u_htr,j_htr,e_ee,ns_sum,spin_deg,e_dc,e_dcc
      REAL a1,a2,alpha
      LOGICAL l_mix
!
! Use around-mean-field limit (true) of atomic limit (false)
!
!
! Loop over atoms
!
      e_off = 0.0
      DO i_u = 1, n_u
        iType = u_in(i_u)%atomType
        l = u_in(i_u)%l
        l_mix = .false. ! lda_u(itype)%l_amf ! switched off
        u_htr = f0(n)/hartree_to_ev_const
        IF (l.EQ.1) THEN
          j_htr = f2(i_u)/(5*hartree_to_ev_const)
        ELSEIF (l.EQ.2) THEN
          j_htr = 1.625*f2(i_u)/(14*hartree_to_ev_const)
        ELSEIF (l.EQ.3) THEN
          j_htr = (286.+195*451/675+250*1001/2025)*f2(i_u)/(6435*hartree_to_ev_const)
        ENDIF
!
!------------------------------------------------------------------------+
!  offd    --                    offd                  1        s   1    |
! V     =  >  (  <m,p|V|q,m'> ) n     + d    ( -U (n - -) + J (n  - -) ) |
!  m,m'    --                    p,q     m,m'          2            2    |
!        p,q                                                             |
!------------------------------------------------------------------------+     
! initialise v_mmp
!
        v_mmp21(:,:,n) = cmplx(0.0,0.0)
!
! outer spin loop - set up v_mmp
!
        DO m = -l,l
          DO mp =-l,l 
             DO p = -l,l
               DO q = -l,l
                 v_mmp21(m,mp,i_u) = v_mmp21(m,mp,i_u) - &
                                   u(m,p,q,mp,i_u) * n_mmp21(p,q,i_u) 
               ENDDO
             ENDDO
          ENDDO ! m' loop
        ENDDO   ! m  loop

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

        IF (irank==0) THEN
          WRITE (6,'(a6,i4)') 'atom: ',itype
          WRITE (6,'(a6,f12.6)') 'e_off:',e_off
          WRITE (6,'(14f10.5)') v_mmp21(:,:,i_u)
        ENDIF
     
        ENDIF
      ENDDO ! loop over atoms

      e_off = e_off / 2
      
      END SUBROUTINE v_mmp_21
END MODULE m_vmmp21
