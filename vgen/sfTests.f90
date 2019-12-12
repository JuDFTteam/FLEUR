MODULE m_sfTests
   IMPLICIT NONE
CONTAINS
   SUBROUTINE plotBtest(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, &
                        noco, xcB, div, phi, cvec, corrB, div2)
      USE m_plot

      TYPE(t_stars),     INTENT(IN)    :: stars
      TYPE(t_atoms),     INTENT(IN)    :: atoms
      TYPE(t_sphhar),    INTENT(IN)    :: sphhar
      TYPE(t_vacuum),    INTENT(IN)    :: vacuum
      TYPE(t_input),     INTENT(IN)    :: input
      TYPE(t_oneD),      INTENT(IN)    :: oneD
      TYPE(t_sym),       INTENT(IN)    :: sym
      TYPE(t_cell),      INTENT(IN)    :: cell
      TYPE(t_noco),      INTENT(IN)    :: noco
      TYPE(t_potden),               INTENT(IN)    :: div, phi, div2
      TYPE(t_potden), DIMENSION(3), INTENT(IN)    :: xcB, cvec, corrB

      LOGICAL                          :: xsf

      CALL checkplotinp()

      CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, &
                  .FALSE., .FALSE., 'bInitial            ', xcB(1), xcB(1), xcB(2), xcB(3))

      CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, &
                  .FALSE., .FALSE., 'div                 ', div)

      CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, &
                  .FALSE., .TRUE., 'phiDiv              ', phi)

      CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, &
                  .FALSE., .FALSE., 'gradPhiDiv          ', cvec(1), cvec(1), cvec(2), cvec(3))

      CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, &
                  .FALSE., .FALSE., 'bCorrected          ', corrB(1), corrB(1), corrB(2), corrB(3))

      CALL savxsf(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, &
                  .FALSE., .FALSE., 'divCorrected        ', div2)

      INQUIRE(file="bInitial_f.xsf",exist=xsf)

      IF (xsf) THEN
         OPEN  (120, FILE='bInitial_f.xsf', STATUS='OLD')
         CLOSE (120, STATUS="DELETE")
         OPEN  (120, FILE='gradPhiDiv_f.xsf', STATUS='OLD')
         CLOSE (120, STATUS="DELETE")
         OPEN  (120, FILE='bCorrected_f.xsf', STATUS='OLD')
         CLOSE (120, STATUS="DELETE")
      END IF

   END SUBROUTINE plotBtest

   SUBROUTINE buildAtest(stars,atoms,sphhar,vacuum,input,noco,sym,cell,itest,Avec,denMat,factor)
      USE m_mt_tofrom_grid
      USE m_pw_tofrom_grid
      USE m_xcBfield

      TYPE(t_stars),                INTENT(IN)     :: stars
      TYPE(t_atoms),                INTENT(IN)     :: atoms
      TYPE(t_sphhar),               INTENT(IN)     :: sphhar
      TYPE(t_vacuum),               INTENT(IN)     :: vacuum
      TYPE(t_input),                INTENT(IN)     :: input
      TYPE(t_noco),                 INTENT(IN)     :: noco
      TYPE(t_sym),                  INTENT(IN)     :: sym
      TYPE(t_cell),                 INTENT(IN)     :: cell
      INTEGER,                      INTENT(IN)     :: itest
      TYPE(t_potden), DIMENSION(3), INTENT(OUT)    :: Avec
      TYPE(t_potden), OPTIONAL,     INTENT(IN)     :: denMat
      REAL,           OPTIONAL,     INTENT(IN)     :: factor

      INTEGER                                      :: nsp, n, kt, kt2, ir, i, j, k, ifftxc3, ind
      REAL                                         :: r, th, ph, x, y, z, dx, dy, dz
      REAL, ALLOCATABLE                            :: thet(:), phi(:), A_temp(:,:,:)!space grid, index
      TYPE(t_gradients)                            :: grad

      REAL                            :: vec1(3), vec2(3), vec3(3), zero(3), point(3)
      INTEGER                         :: grid(3)


      IF (itest.EQ.0) THEN
         RETURN
      END IF

      IF (PRESENT(denMat)) THEN
         CALL makeVectorField(sym,stars,atoms,sphhar,vacuum,input,noco,denMat,factor,Avec)
         RETURN
      END IF

      nsp = atoms%nsp()
      ifftxc3=stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft

      DO i=1,3
         CALL Avec(i)%init_potden_simple(stars%ng3,atoms%jmtd,sphhar%nlhd, &
                       atoms%ntype,atoms%n_u,1,.FALSE.,.FALSE., &
                       POTDEN_TYPE_POTTOT,vacuum%nmzd,vacuum%nmzxyd,stars%ng2)
         ALLOCATE(Avec(i)%pw_w,mold=Avec(i)%pw)
         Avec(i)%mt    = 0
         ! Temporary.
         Avec(i)%pw    = 0
         Avec(i)%pw_w  = CMPLX(0.0,0.0)
         Avec(i)%vacxy = 0
         Avec(i)%vacz  = 0
      END DO

      ALLOCATE (thet(atoms%nsp()),phi(atoms%nsp()))

      CALL init_mt_grid(1, atoms, sphhar, .TRUE., sym, thet, phi)
      CALL init_pw_grid(.TRUE.,stars,sym,cell)
      !--------------------------------------------------------------------------
      ! Test case 1.
      ! In MT: radial function f(r)=r*R_MT in direction e_r.
      ! In interstitial: f(r_vec)=r_vec.
      ! So by construction, the divergence for a set of atoms with R_MT=rmt(n)
      ! for every n should be an isotropic 3*R_MT.
      ! TODO: Add interstitial part and if-test (1---> this test).

      DO n=1,atoms%ntype
         ALLOCATE (A_temp(atoms%jri(n)*nsp,3,1))
         kt = 0
         DO ir = 1, atoms%jri(n)
            r = atoms%rmsh(ir, n)
            !print *, 'ir'
            !print *, ir
            DO k = 1, nsp
               th = thet(k)
               ph = phi(k)
               !print *, k
	            !print *, th/pi_const
	            !print *, ph/pi_const
	            !(r/atoms%rmt(n))
	            A_temp(kt+k,1,1)=(r**2)*r**2!*SIN(th)*COS(ph)
	            A_temp(kt+k,2,1)=(r**2)*r**2!*SIN(th)*SIN(ph)
	            A_temp(kt+k,3,1)=(r**2)*r**2!*COS(th)
            ENDDO ! k
            kt = kt + nsp
         END DO ! ir

         CALL mt_from_grid(atoms,sym, sphhar, n, 1, A_temp(:,1,:), Avec(1)%mt(:,0:,n,:))
         CALL mt_from_grid(atoms,sym, sphhar, n, 1, A_temp(:,2,:), Avec(2)%mt(:,0:,n,:))
         CALL mt_from_grid(atoms,sym, sphhar, n, 1, A_temp(:,3,:), Avec(3)%mt(:,0:,n,:))

         !print *, 'A_z*r^2 3rd entry before fromto grid'
         !print *, A_temp(3,3,1)

         !CALL mt_to_grid(.FALSE., 1, atoms, sphhar, Avec(3)%mt(:,0:,n,:), n, grad, A_temp(:,3,:))

         !kt=0
         !DO ir = 1, atoms%jri(n)
         !   r = atoms%rmsh(ir, n)
         !   DO k = 1, nsp
	      !      A_temp(kt+k,3,1)=A_temp(kt+k,3,1)*r**2
         !   ENDDO ! k
         !   kt = kt + nsp
         !END DO ! ir

         !print *, 'A_z*r^2 3rd entry after fromto grid'
         !print *, A_temp(3,3,1)
         DEALLOCATE (A_temp)
      END DO ! n

      ALLOCATE (A_temp(ifftxc3,3,1))

      grid = (/stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft/)
      vec1 = (/1.,0.,0./)
      vec2 = (/0.,1.,0./)
      vec3 = (/0.,0.,1./)
      zero = (/0.,0.,0./)
      vec1=matmul(cell%amat,vec1)
      vec2=matmul(cell%amat,vec2)
      vec3=matmul(cell%amat,vec3)
      zero=matmul(cell%amat,zero)

      DO k = 0, grid(3)-1
         DO j = 0, grid(2)-1
            DO i = 0, grid(1)-1

               point = zero + vec1*REAL(i)/(grid(1)-1) +&
                              vec2*REAL(j)/(grid(2)-1) +&
                              vec3*REAL(k)/(grid(3)-1)

               ind = k*grid(2)*grid(1) + j*grid(1) + i + 1

               A_temp(ind,1,1)=SIN(i*2*pi_const/grid(1))
               A_temp(ind,2,1)=SIN(j*2*pi_const/grid(2))
               A_temp(ind,3,1)=SIN(k*2*pi_const/grid(3))

            END DO
         END DO
      END DO !z-loop

      DO i=1,3
         CALL pw_from_grid(.TRUE.,stars,.TRUE.,A_temp(:,i,:),Avec(i)%pw,Avec(i)%pw_w)
      END DO

      ! End test case 1.
      !--------------------------------------------------------------------------
      CALL finish_mt_grid
      CALL finish_pw_grid()

   END SUBROUTINE buildAtest

   SUBROUTINE sftest(mpi,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,itest,denMat,factor)
      USE m_xcBfield
      USE m_divergence
      USE m_analysistests
      TYPE(t_field),                INTENT(INOUT)  :: field
      TYPE(t_mpi),                  INTENT(IN)     :: mpi
      TYPE(t_stars),                INTENT(IN)     :: stars
      TYPE(t_atoms),                INTENT(IN)     :: atoms
      TYPE(t_sphhar),               INTENT(IN)     :: sphhar
      TYPE(t_vacuum),               INTENT(IN)     :: vacuum
      TYPE(t_input),                INTENT(IN)     :: input
      TYPE(t_oneD),                 INTENT(IN)     :: oneD
      TYPE(t_sym),                  INTENT(IN)     :: sym
      TYPE(t_cell),                 INTENT(IN)     :: cell
      TYPE(t_noco),                 INTENT(IN)     :: noco
      INTEGER,                      INTENT(IN)     :: itest
      TYPE(t_potden), OPTIONAL,     INTENT(IN)     :: denMat
      REAL,           OPTIONAL,     INTENT(IN)     :: factor

      TYPE(t_potden), DIMENSION(3)                 :: aVec, cvec, corrB
      TYPE(t_potden)                               :: div, phi, checkdiv
      INTEGER                                      :: i, n, lh, l
      REAL                                         :: g(atoms%jmtd)

      REAL :: radii(atoms%jmtd,atoms%ntype), funcsr(atoms%jmtd,atoms%ntype,3), trueder(atoms%jmtd,atoms%ntype,3), trueint(atoms%jmtd,atoms%ntype,3), &
              testders(atoms%jmtd,atoms%ntype,3,5), testints(atoms%jmtd,atoms%ntype,3,5)

      CALL buildfunctions(atoms,radii,funcsr,trueder,trueint)
      CALL derivtest(atoms,radii,funcsr,testders)
      CALL integtest(atoms,radii,funcsr,testints)

      ! Test: Build a field and either compare with theoretical results or check,
      !       whether the sourcefree routine made it sourcefree.

      IF (PRESENT(denMat)) THEN
        CALL buildAtest(stars,atoms,sphhar,vacuum,input,noco,sym,cell,1,aVec,denMat,factor)
        CALL sourcefree(mpi,field,stars,atoms,sphhar,vacuum,input,oneD,sym,cell,noco,aVec,div,phi,cvec,corrB,checkdiv)
        CALL plotBtest(stars, atoms, sphhar, vacuum, input, oneD, sym, cell, noco, aVec, div, phi, cvec, corrB, checkdiv)
      ELSE
        CALL buildAtest(stars,atoms,sphhar,vacuum,input,noco,sym,cell,0,aVec)
        RETURN
      END IF

      !testDen(3)%mt(:,1,:,1)=testDen(3)%mt(:,0,:,1)*atoms%rmsh
      !testDen(3)%mt(:,1:,:,:)=0.0
      !testDen(3)%mt(:,2:,:,:)=0.0
      !testDen(3)%mt(:,0,:,:)=0.0
      !testDen(2)%mt(:,:,:,:)=0.0
      !testDen(1)%mt(:,:,:,:)=0.0
      !testDen(3)%mt(:,0,:,1)*atoms%rmsh
      !testDen(3)%mt(:,0,:,1)=0.0

      END SUBROUTINE sftest

      SUBROUTINE difftester(atoms,n,l,f,g)
         USE m_types
         USE m_grdchlh

         TYPE(t_atoms), INTENT(IN) :: atoms
         INTEGER, INTENT(IN) :: n,l
         REAL, INTENT(IN) :: f(atoms%jri(n))
         REAL, INTENT(OUT) :: g(atoms%jri(n))

         REAL :: dfr(atoms%jri(n)), d2fr2(atoms%jri(n))

         CALL grdchlh(1, 1, atoms%jri(n), atoms%dx(n), atoms%rmsh(1, n), f, 5, dfr, d2fr2)

         g=(dfr-l*f/atoms%rmsh(:, n))!/(atoms%rmsh(1, n)**l) !must NOT be divergent towards the core

      END SUBROUTINE difftester

!   SUBROUTINE lhlmtest()
!     ALLOCATE (flh(atoms%jri(1),0:sphhar%nlh(atoms%ntypsy(1))),flm(atoms%jri(1),sphhar%nlh(atoms%ntypsy(1))+1),flh2(atoms%jri(1),0:sphhar%nlh(atoms%ntypsy(1))))
!     flh=inDen%mt(:,:,1,1)
!     flh(:,1)=-flh(:,0)
!     flh(:,2)=0*flh(:,0)
!     flh(:,3)=flh(:,0)
!     flh(:,4)=flh(:,0)
!     flh(:,5)=2*flh(:,0)
!     flh(:,6)=3*flh(:,0)
!     flh(:,7)=4*flh(:,0)
!     flh(:,8)=5*flh(:,0)
!     CALL lh_to_lm(atoms, sphhar, 1, flh, flm)
!     CALL lh_from_lm(atoms, sphhar, 1, flm, flh2)

!   END SUBROUTINE lhlmtest

END MODULE m_sfTests
