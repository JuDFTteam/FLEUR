MODULE m_vmtxc
CONTAINS
  SUBROUTINE vmtxc( DIMENSION,sphhar,atoms, rho,xcpot,input,sym, vr, excr,vxr)
    !     ********************************************************************
    !     fit spherical-harmonics expansion of exchange-correlation potential*
    !     inside muffint-tin spheres and add it to coulomb potential         *
    !                                       c.l.fu and r.podloucky           *
    !     ********************************************************************
    !     instead of vmtxcor.f: the different exchange-correlation 
    !     potentials defined through the key icorr are called through 
    !     the driver subroutine vxcall.f, subroutines vectorized
    !     ** r.pentcheva 22.01.96
    !     *********************************************************
    !     angular mesh calculated on speacial gauss-legendre points
    !     in order to use orthogonality of lattice harmonics and
    !     avoid a least sqare fit
    !     ** r.pentcheva 04.03.96
    !     *********************************************************
    !     In the previous version the parameter statement nspd has
    !     been used synonymous to nsp. I have introduced the variable
    !     nsp in order to change this. It has cause problems with
    !     subroutine: sphglpts.f
    !     ** s.bluegel, IFF, 30.July.97
    !     *********************************************************

    USE m_lhglpts
    USE m_gaussp
    USE m_xcall, ONLY : vxcall,excall
    USE m_types
    IMPLICIT NONE
    TYPE(t_xcpot),INTENT(IN)       :: xcpot
    TYPE(t_dimension),INTENT(IN)   :: DIMENSION
    TYPE(t_input),INTENT(IN)       :: input
    TYPE(t_sym),INTENT(IN)         :: sym
    TYPE(t_sphhar),INTENT(IN)      :: sphhar
    TYPE(t_atoms),INTENT(IN)       :: atoms
    !     ..
    !     .. Scalar Arguments ..
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: rho(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,DIMENSION%jspd)
    REAL,    INTENT (INOUT):: vr (atoms%jmtd,0:sphhar%nlhd,atoms%ntype,DIMENSION%jspd)
    REAL,    INTENT (OUT)  :: excr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype)
    REAL,    INTENT (INOUT),OPTIONAL:: vxr(atoms%jmtd,0:sphhar%nlhd,atoms%ntype,DIMENSION%jspd)

    !     ..
    !     .. Local Scalars ..
    REAL elh,rr2,r2rho,vlh
    INTEGER jr,js,k,lh,n,nd,i,nsp,nat
    !     ..
    !     .. Local Arrays ..
    REAL vx(DIMENSION%nspd,DIMENSION%jspd),vxc(DIMENSION%nspd,DIMENSION%jspd),pos0(3)
    REAL rhoxc(DIMENSION%nspd,DIMENSION%jspd),exc(DIMENSION%nspd),wt(DIMENSION%nspd),rx(3,DIMENSION%nspd)
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
    !$OMP SHARED(atoms,input,sphhar,xcpot,rho,ylh,nsp,wt,vr,vxr,excr)&
    !$OMP PRIVATE(elh,rr2,r2rho,vlh,jr,js,k,lh,n,nd,i,nat)&
    !$OMP PRIVATE(vx,vxc,pos0,rhoxc,exc)
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
                r2rho = rr2*rho(jr,lh,n,js)
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
          CALL vxcall (6,xcpot%icorr,xcpot%krla,input%jspins, nsp,nsp,rhoxc, vx,vxc) 
          !     ----> now determine the corresponding potential number 
          DO js = 1,input%jspins
             !
             ! ---> multiplikate vxc with the weights of the k-points    
             !
             DO  k = 1,nsp
                vxc(k,js) = vxc(k,js)*wt(k)

                vx (k,js) = vx (k,js)*wt(k)
             ENDDO
             DO  lh = 0,sphhar%nlh(nd)

                !
                ! ---> determine the corresponding potential number through gauss integration
                !
                vlh=DOT_PRODUCT( vxc(:nsp,js),ylh(:nsp,lh,nd))

                !
                ! ---> add to the given potential

                vr(jr,lh,n,js) = vr(jr,lh,n,js) + vlh

                IF (PRESENT(vxr) ) THEN
                   vlh = 0
                   DO  k = 1,nsp
                      vlh = vlh + vx(k,js)*ylh(k,lh,nd)
                   END DO

                   vxr(jr,lh,n,js) = vxr(jr,lh,n,js) + vlh
                ENDIF
             ENDDO
          ENDDO
          IF (input%total) THEN 
             !
             !     calculate the ex.-cor energy density
             !
             CALL excall(6,xcpot%icorr,xcpot%krla,input%jspins, nsp,nsp,rhoxc, exc) 
          END IF
          !     ----> now determine the corresponding energy density number 
          !
          ! ---> multiplikate exc with the weights of the k-points    
          !
          DO  k = 1,nsp
             exc(k) = exc(k)*wt(k)
          ENDDO
          DO  lh = 0,sphhar%nlh(nd)

             !
             ! ---> determine the corresponding potential number through gauss integration
             !
             elh=DOT_PRODUCT(exc(:nsp),ylh(:nsp,lh,nd))

             !
             ! ---> add to the given potential

             excr(jr,lh,n) = elh
          ENDDO

       ENDDO

    ENDDO
    !$OMP END PARALLEL DO
    !

    RETURN
  END SUBROUTINE vmtxc
END MODULE m_vmtxc
