MODULE m_cdnovlp
  USE m_fleurenv
CONTAINS
  SUBROUTINE cdnovlp()

    !*****************************************************************

    !     Added calculation of force contribution from coretails
    !     outside of their native muffin-tin spheres, i.e. in the
    !     interstitial region and other muffin-tins; only for bulk.
    !     refer to Klüppelberg et al., PRB 91 035105 (2015)
    !     Aaron Klueppelberg, Oct. 2015
    !*****************************************************************
    !
    USE m_intgr, ONLY : intgr3,intgz0
    USE m_constants
    USE m_spgrot
    USE m_cylbes
    USE m_dcylbs
    USE m_rcerf
    USE m_force_a4_add
    USE m_od_cylbes
    USE m_od_chirot
    USE m_od_types, ONLY : od_inp, od_sym
#ifdef CPP_MPI
    USE m_mpi_bc_st
#endif
    USE m_sphbes
    USE m_phasy1
    USE m_loddop

    IMPLICIT NONE


    !     ..
    !     .. Local Scalars ..
    COMPLEX sf,carg,value,slope, phase
    REAL    rdum,radint
    REAL    ai,ar,a4,dif,dxx,f11,f12,g,gr,gz,qfin,qfout,dtildh, rkappa,sign,signz,tol,x,z,zvac,alpha3, g2,phi,gamma,qq
    INTEGER ig3,imz,ir,ivac,j,j1,k,kz,k1,k2,l_cutoff,m0, n,nz,nrz,nat1,nat2,nzvac,n_out_p,nat,irec2,irec3,irec1,m,gzi
    LOGICAL tail
    !     ..
    !     .. Local Arrays ..
    REAL,    ALLOCATABLE :: qf(:)
    REAL    acoff(ntypd),alpha(ntypd),rho_out(2), rvec(3),kre(3)
    REAL    phas(nop),rhohelp(msh),rat(msh,ntypd)
    INTEGER kr(3,nop),mshc(ntypd)
    REAL    kro(3,ods%nop),fJ(-odi%M:odi%M),dfJ(-odi%M:odi%M)
    COMPLEX phaso(ods%nop), rhoc2(msh,ntypd)

    INTEGER :: nu,maxl,dir,iter,nd,lh,l,lm,jm,sym,ierr(3)
    REAL    :: s13,s23,factor,i2
    REAL   , ALLOCATABLE :: vz(:,:,:),vr(:,:,:,:),vr0(:,:,:)
    REAL   , ALLOCATABLE :: vr2(:,:,:),bsl(:,:),integrandr(:)
    REAL   , ALLOCATABLE :: vrigrid(:),integrand(:,:)
    COMPLEX :: ycomp1(3,-1:1),gv(3),sm,kcmplx(3)
    COMPLEX, ALLOCATABLE :: vpw(:,:),vzxy(:,:,:,:),ffonat(:,:)
    COMPLEX, ALLOCATABLE :: pylm2(:,:,:)
    CHARACTER(8) :: dop,iop,name(10)

    INTEGER :: minpar,left
    LOGICAL :: l_f2

    INTEGER :: nkpt_pT(0:isize-1),ioffset_pT(0:isize-1),kp
    INTEGER, ALLOCATABLE :: n1(:),n2(:)

    !     ..
    !     .. External Functions ..
    REAL     diflgr
    EXTERNAL diflgr
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC cabs,cmplx,int,log,exp,max,min,real,sin,sqrt
    !     ..
    DATA  tol /1.0e-10/
    !
   


    !       Allocate the force arrays in the routine force_a4_add.f
    CALL alloc_fa4_arrays(ntypd,jspd)
    force_a4_mt(:,:,jspin) =  0.0
    force_a4_is(:,:,jspin) =  0.0


    !       lattice/spherical harmonics related variables
    ycomp1=0.0
    s13 = sqrt(1.0/3.0)
    s23 = sqrt(2.0/3.0)
    ycomp1(3,0) = cmplx(2.0*s13,0.0)
    ycomp1(1,-1) = cmplx(s23,0.0)
    ycomp1(2,-1) = cmplx(0.0,-s23)
    ycomp1(1,1) = cmplx(-s23,0.0)
    ycomp1(2,1) = cmplx(0.0,-s23)

    !       read in effective potential
    ALLOCATE ( vpw(n3d,jspd),vzxy(nmzxyd,odi%n2d-1,2,jspd) )
    ALLOCATE ( vz(nmzd,2,4),vr(jmtd,0:nlhd,ntypd,jspd) )
    ALLOCATE ( vr0(jmtd,ntypd,jspd),vr2(jmtd,0:nlhd,ntypd) )
    nu=10
    OPEN (nu,file='pottot',form='unformatted',status='old')
    CALL loddop(jspd,n3d,odi%n2d,nmzxyd,nmzd,jmtd,nlhd,ntypd, jspin,nq3,odi%nq2,nvac,ntype,invs,invs2,film, nlh,jri,ntypsd,ntypsy,nu,natd,neq, iop,dop,iter,vr,vpw,vz,vzxy,name)
    CLOSE(nu)

    !       the l = 0 component of the potential is multiplied by r/sqrt(4 pi), 
    !       for simple use, this is corrected here
    DO n = 1,ntype
       vr2(:,0,n) = sfp*vr(:,0,n,jspin)/rmsh(:,n)
       vr2(:,1:,n) = vr(:,1:,n,jspin)
    END DO ! n

    ALLOCATE ( ffonat(3,n3d*nop),integrandr(jmtd) )
    ALLOCATE ( pylm2( (lmaxd+1)**2,3,nop ) )
    ALLOCATE ( bsl(jmtd,0:atoms%lmax),integrand(jmtd,0:atoms%lmax) )

    !       (f)orce(f)actor(on)(at)oms calculation, parallelization in k
    ffonat = 0.0
    nat = 1
    DO n = 1,ntypd
       nd = ntypsy(nat)
       g = -0.1 ! g is the norm of a star and can't be negative, this is to initialize a check if the norm between stars has changed
       !       on each processor, calculate a certain consecutive set of k
       kp = 0
       DO k = 2,stars%ng3 ! for k = 1 (G = 0), grad rho_core^alpha is zero
          IF (ABS(g-sk3(k)).GT.tol) THEN ! only calculate new spherical Bessel functions if the length of the star vector has changed
             g=sk3(k)
             !generate spherical Bessel functions up to maxl for the radial grid
             DO j = 1,jri(n)
                gr = g * rmsh(j,n)
                CALL sphbes(MAXVAL(llh(0:nlh(nd))),gr,bsl(j,:))
                bsl(j,:) = bsl(j,:) * rmsh(j,n)**2
             END DO
          END IF
          !         as phasy1, but with i\vec{G} in it, i.e. gradient of plane wave, only for atom n and star k
          CALL phasy2(nop,lmaxd,lmax(n),fpi_const,taual(:,nat), bmat,kv3(:,k),tau,mrot,symor,invtab, pylm2)
          !         construct and evaluate radial integral int_0^R_{beta} r^2 j_{l}(Gr) V_{eff,l}^{beta}(r) dr
          !         then, multiply by pylm2 times clnu
          DO lh = 0,nlh(nd)
             l = llh(lh,nd)
             integrandr(:) = bsl(:,l) * vr2(:,lh,n)
             CALL intgr3(integrandr,rmsh(1,n),dx(n),jri(n),factor)
             DO j = 1,nop
                sym = kp*nop + j
                DO dir = 1,3
                   sm = 0.0
                   DO jm = 1,nmem(lh,nd)
                      lm = l*(l+1) + mlh(jm,lh,nd) + 1
                      sm = sm + CONJG(clnu(jm,lh,nd)) * pylm2(lm,dir,j)
                   END DO ! jm
                   ffonat(dir,sym) = ffonat(dir,sym) + factor * sm
                END DO ! dir
             END DO ! sym
          END DO ! lh

          kp = kp+1
       END DO
       nat = nat+neq(n)
    END DO ! n atom type

    !
    !----> prepare local array to store pw-expansion of pseudo core charge
    ! 
    !----> (1) set up radial mesh beyond muffin-tin radius
    !      (2) cut_off core tails from noise 
    !
    ntype_loop:DO  n = 1 , ntype
       IF (ncst(n).GT.0) THEN
          rat(:jri(n),n) = rmsh(:jri(n),n)
          dxx = exp(dx(n))
          DO j = jri(n) + 1 , msh
             rat(j,n) = rat(j-1,n)*dxx
          ENDDO
          rh(:msh,n)=rh(:msh,n)/ (fpi_const*rat(:msh,n)*rat(:msh,n))
          DO j = msh , jri(n) , -1
             IF ( rh(j,n) .gt. tol ) then
                mshc(n) = j
                CYCLE ntype_loop
             END IF
          ENDDO
          mshc(n) = jri(n)
       ENDIF
    END DO ntype_loop
    !
    !-----> the core density inside the spheres is replaced by a
    !       gaussian-like pseudo density : n(r) = acoff*exp(-alpha*r*r)
    !       acoff and alpha determined to obtain a continous and 
    !       differentiable density at the sphere boundary.
    !       IF mshc = jri  either core tail too small or no core (i.e. H)
    !
    nat = 1
    DO  n = 1,ntype
       IF ((mshc(n).GT.jri(n)).AND.(ncst(n).GT.0)) THEN

          j1 = jri(n) - 1
          alpha(n) = log( rh(j1,n) / rh(jri(n),n) )
          alpha(n) = alpha(n) / ( rmt(n)*rmt(n)*( 1.0-exp( -2.0*dx(n) ) ) )
          acoff(n) = rh(jri(n),n) * exp( alpha(n)*rmt(n)*rmt(n) )
          WRITE (6,FMT=8010) alpha(n),acoff(n)
          DO j = 1,jri(n) - 1
             rh(j,n) = acoff(n) * exp( -alpha(n)*rat(j,n)**2 )
          ENDDO

          !           Subtract pseudo density contribution from own mt sphere from mt forces
          ALLOCATE ( integrand(jmtd,3) )
          integrandr = 0.0
          integrand  = 0.0
          DO j = 1,jri(n)
             integrandr(j) = -alpha(n) * acoff(n) * rmsh(j,n)**3 *sfp * exp(-alpha(n) * rmsh(j,n)**2) !*2 ! factor of two missing? grad e^{-alpha*r^2} = -2alpha\vec{r}e^{-alpha*r^2}
          END DO ! j radial mesh

          DO lh = 0,nlh(ntypsy(nat))
             IF (llh(lh,ntypsy(nat)).ne.1) CYCLE

             gv = 0.0
             DO jm = 1,nmem(lh,ntypsy(nat))
                m = mlh(jm,lh,ntypsy(nat))

                DO dir = 1,3
                   gv(dir) = gv(dir) + ycomp1(dir,m) * clnu(jm,lh,ntypsy(nat)) ! why not conjg?
                END DO ! dir ection

             END DO ! jm

             DO dir = 1,3
                DO j = 1,jri(n)
                   integrand(j,dir) = integrand(j,dir) - integrandr(j) * vr2(j,lh,n) * real(gv(dir))
                END DO ! j radial mesh
             END DO ! dir ection

          END DO ! lh lattice harmonics

          DO dir = 1,3
             CALL intgr3(integrand(:,dir),rmsh(1,n),dx(n),jri(n), force_a4_mt(dir,n,jspin))
          END DO ! dir ection
          DEALLOCATE ( integrand )

       ELSE
          alpha(n) = 0.0
       ENDIF
       nat = nat+neq(n)
60  ENDDO
    !
8000 FORMAT (/,10x,'core density and its first derivative ', 'at sph. bound. for atom type', i2,' is',3x,2e15.7)
8010 FORMAT (/,10x,'alpha=',f10.5,5x,'acoff=',f10.5)
    !
    !=====> calculate the fourier transform of the core-pseudocharge
    !
    !
    !*****> start loop over the atom type
    !
    nat1 = 1
    DO  n = 1,ntype

       IF ( ( mshc(n) .GT. jri(n) ).AND. ( alpha(n) .GT. tol ) )    THEN

          !
          n_out_p = mshc(n)-jri(n)+1
          !
          ! (1) Form factor for each atom type
          !
          f11 = tpi_const * rmt(n) * rh(jri(n),n) / alpha(n)
          f12 = acoff(n) * ( pi_const/alpha(n) ) * sqrt( pi_const/alpha(n) )
          ar  = sqrt( alpha(n) ) * rmt(n) 
          !
          DO  k = 1,nq3
             g = sk3(k)
             !    first G=0
             IF ( k.eq.1 ) THEN
                ai = 0.0
                !
                ! ---->     calculate form factor inside the mt-sphere
                !           (use analytic integration of gaussian)
                !
                qfin = - f11 + f12 * rcerf(ar,ai)
                !
                ! ---->     calculate form factor outside the mt-sphere
                !           (do numerical integration of tails)
                !

                DO ir = 1 , n_out_p
                   j  = jri(n)+ir-1
                   rhohelp(ir) = rat(j,n) * rat(j,n) * rh(j,n)
                END DO
                CALL intgr3(rhohelp,rat(jri(n),n),dx(n), n_out_p,qfout)
                !--->             have to remove the small r-correction from intgr3
                !                   qfout=qfout-rmt(n)*rhohelp(1)
                if (rhohelp(1)*rhohelp(2).GT.0.0) THEN
                   alpha3 = 1.0 + log(rhohelp(2)/rhohelp(1))/dx(n)
                   IF (alpha3.GT.0.0) qfout = qfout - rat(jri(n),n)*rhohelp(1)/alpha3
                end if


                !       WRITE (999,'(2(i5,1x),f14.10)') n,k,qfout

                qfout = fpi_const * qfout
                !
             ELSE 
                !    then  G>0
                ai = 0.5*g/sqrt(alpha(n))
                gr = g*rmt(n)
                a4 = 0.25/alpha(n)
                !
                ! ---->     calculate form factor inside the mt-sphere
                !           (use analytic integration of gaussian)
                !
                qfin = - f11 * sin(gr)/gr + f12 * rcerf(ar,ai) * exp(-a4*g*g) 
                !
                ! ---->     calculate form factor outside the mt-sphere
                !           (do numerical integration of tails)


                DO ir = 1 , n_out_p
                   j  = jri(n)+ir-1
                   rhohelp(ir) = rat(j,n) * rh(j,n) * sin(g*rat(j,n))
                END DO
                CALL intgr3(rhohelp,rat(jri(n),n),dx(n), n_out_p,qfout)
                !--->             have to remove the small r-correction from intgr3
                !roa...correction.from.intgr3.......................
                if (rhohelp(1)*rhohelp(2).GT.0.0) THEN
                   alpha3 = 1.0 + log(rhohelp(2)/rhohelp(1))/dx(n)
                   IF (alpha3.GT.0.0) qfout = qfout - rat(jri(n),n)*rhohelp(1)/alpha3
                endif
                !roa...end.correction...............................



                !       WRITE (999,'(2(i5,1x),f14.10)') n,k,qfout

                qfout = fpi_const * qfout / g
                !
             END IF
             !
             qf(k)    = (qfin + qfout)/omtil
          end do
          !
          ! (2) structure constant for each atom of atom type
          !
          !
          !    first G=0
          !
          k=1
          !                  IF (l_f2) THEN
          !                    kcmplx is zero for k=1 anyway, so no contribution in interstitial by default
          !                  END IF

          !
          !    then  G>0
          !
          DO  k = 2,nq3
             CALL spgrot(nop,symor,tpi_const,mrot,tau,invtab, kv3(1,k), kr,phas)
             !
             ! ----> start loop over equivalent atoms
             !

             !              generate phase factors for each G, not only for each star, to incorporate the atomic phase factors
             kcmplx = cmplx(0.0,0.0)
             DO j = 1,nop
                x = -tpi_const* ( kr(1,j) * taual(1,nat1) + kr(2,j) * taual(2,nat1) + kr(3,j) * taual(3,nat1) )
                phase = cmplx(cos(x),sin(x))
                !                generate muffin-tin part of core force component
                force_a4_mt(:,n,jspin) = force_a4_mt(:,n,jspin) + qf(k) * phase * nstr(k) * ffonat(:,(k-1)*nop+j)
                kcmplx(:) = kcmplx(:) + kr(:,j) * phase * phas(j) ! should be conjg(phas(j)), but in FLEUR, only real phas(j) are accepted
             END DO ! j
             kcmplx = matmul(kcmplx,bmat) * nstr(k) / nop
             !              generate interstitial part of core force component
             force_a4_is(:,n,jspin) = force_a4_is(:,n,jspin) + qf(k) * CONJG(vpw(k,jspin))*omtil*CMPLX(0.,kcmplx(:))
          END DO
       END IF
       nat1 = nat1 + neq(n)
100 END DO

  END SUBROUTINE cdnovlp
   END MODULE m_cdnovlp
