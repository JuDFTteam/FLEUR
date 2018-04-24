!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_vmtxc
CONTAINS
  SUBROUTINE vmtxc( DIMENSION,sphhar,atoms,den,xcpot,input,sym, vxc, exc,vx)
    !     ********************************************************************
    !     Calculate the LDA xc-potential in the MT-spheres
    !     *********************************************************

    USE m_lhglpts
    USE m_gaussp
    USE m_types
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN)      :: xcpot
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms
    TYPE(t_potden),INTENT(IN)      :: den
    TYPE(t_potden),INTENT(INOUT)   :: vxc,exc
    TYPE(t_potden),INTENT(INOUT),optional     :: vx
    
    !     ..
    !     .. Local Scalars ..
    REAL elh,rr2,r2rho,vlh
    INTEGER jr,js,k,lh,n,nd,i,nsp,nat
    !     ..
    !     .. Local Arrays ..
    REAL v_x(DIMENSION%nspd,DIMENSION%jspd),v_xc(DIMENSION%nspd,DIMENSION%jspd),pos0(3)
    REAL rhoxc(DIMENSION%nspd,DIMENSION%jspd),e_ex(DIMENSION%nspd),wt(DIMENSION%nspd),rx(3,DIMENSION%nspd)
    REAL ylh(DIMENSION%nspd,0:sphhar%nlhd,sphhar%ntypsd)
    !     ..
    !     generates nspd points on a sherical shell with radius 1.0
    !     angular mesh equidistant in phi, 
    !     theta are zeros of the legendre polynomials
    !     
    CALL gaussp(atoms%lmaxd,rx,wt)

    nsp = DIMENSION%nspd
    !
    !     generates the lattice harmonics on the angular mesh 
    !     
    CALL lhglpts( sphhar,atoms, rx,nsp, sym, ylh)
    !
    !
    !     loop over topologically non-equivalent atoms
    !
    !$OMP PARALLEL DO DEFAULT(NONE)&
    !$OMP SHARED(atoms,input,sphhar,xcpot,den,ylh,nsp,wt,vxc,vx,exc)&
    !$OMP PRIVATE(elh,rr2,r2rho,vlh,jr,js,k,lh,n,nd,i,nat)&
    !$OMP PRIVATE(v_x,v_xc,pos0,rhoxc,e_ex)
    DO n = 1,atoms%ntype
       nat=SUM(atoms%neq(:n-1))+1
       nd = atoms%ntypsy(nat)
       !
       !     loop over radial mesh
       !
       DO  jr = 1,atoms%jri(n)
          rr2 = 1.e0/ (atoms%rmsh(jr,n)*atoms%rmsh(jr,n))
          DO  js = 1,input%jspins
             rhoxc(:,js) = 0.
             DO  lh = 0,sphhar%nlh(nd)
                r2rho = rr2*den%mt(jr,lh,n,js)
                !
                !     generate the densities on an angular mesh
                !
                DO  k = 1,nsp
                   rhoxc(k,js) = rhoxc(k,js) + ylh(k,lh,nd)*r2rho
                ENDDO
             ENDDO
          ENDDO
          !     calculate the ex.-cor. potential
          !
          CALL xcpot%get_vxc(input%jspins,rhoxc(:nsp,:), v_xc,v_x) 
          !     ----> now determine the corresponding potential number 
          DO js = 1,input%jspins
             !
             ! ---> multiplikate vxc with the weights of the k-points    
             !
             DO  k = 1,nsp
                v_xc(k,js) = v_xc(k,js)*wt(k)

                v_x (k,js) = v_x (k,js)*wt(k)
             ENDDO
             DO  lh = 0,sphhar%nlh(nd)

                !
                ! ---> determine the corresponding potential number through gauss integration
                !
                vlh=DOT_PRODUCT( v_xc(:nsp,js),ylh(:nsp,lh,nd))

                !
                ! ---> add to the given potential

                vxc%mt(jr,lh,n,js) = vxc%mt(jr,lh,n,js) + vlh

                IF (PRESENT(vx) ) THEN
                   vlh = 0
                   DO  k = 1,nsp
                      vlh = vlh + v_x(k,js)*ylh(k,lh,nd)
                   END DO

                   vx%mt(jr,lh,n,js) = vx%mt(jr,lh,n,js) + vlh
                ENDIF
             ENDDO
          ENDDO
          IF (ALLOCATED(exc%mt)) THEN 
             !
             !     calculate the ex.-cor energy density
             !
             CALL xcpot%get_exc(input%jspins,rhoxc(:nsp,:), e_ex) 
             !     ----> now determine the corresponding energy density number 
             !
             ! ---> multiplikate exc with the weights of the k-points    
             !
             DO  k = 1,nsp
                e_ex(k) = e_ex(k)*wt(k)
             ENDDO
             DO  lh = 0,sphhar%nlh(nd)
                
                !
                ! ---> determine the corresponding potential number through gauss integration
                !
                elh=DOT_PRODUCT(e_ex(:nsp),ylh(:nsp,lh,nd))
                
                !
                ! ---> add to the given potential
                
                exc%mt(jr,lh,n,1) = elh
             ENDDO
          END IF
          
       ENDDO

    ENDDO
    !$OMP END PARALLEL DO
    !

    RETURN
  END SUBROUTINE vmtxc
END MODULE m_vmtxc
