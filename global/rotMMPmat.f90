MODULE m_rotMMPmat

   USE m_constants
   USE m_types_sym

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: rotMMPmat

   INTERFACE rotMMPmat
      PROCEDURE :: rotMMPmat_dwgn, rotMMPmat_angle, rotMMPmat_sym_op
      PROCEDURE :: rotMMPmat_angle_one_spin, rotMMPmat_dwgn_one_spin, rotMMPmat_sym_op_one_spin
   END INTERFACE

   CONTAINS

   PURE FUNCTION rotMMPmat_dwgn(mmpmat,dwgn,dwgnp,su, spin_rotation, real_space_rotation) Result(mmpmatOut)

      COMPLEX,           INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX, OPTIONAL, INTENT(IN)  :: dwgn(-lmaxU_const:,-lmaxU_const:)
      COMPLEX, OPTIONAL, INTENT(IN)  :: dwgnp(-lmaxU_const:,-lmaxU_const:)
      COMPLEX, OPTIONAL, INTENT(IN)  :: su(:,:)
      LOGICAL, OPTIONAL, INTENT(IN)  :: spin_rotation
      LOGICAL, OPTIONAL, INTENT(IN)  :: real_space_rotation

      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:,:)

      COMPLEX :: d(2,2)
      INTEGER :: ispin,m,mp
      LOGICAL :: spin_rotation_arg, real_space_rotation_arg

      spin_rotation_arg = .FALSE.
      IF(PRESENT(spin_rotation)) spin_rotation_arg = spin_rotation

      real_space_rotation_arg = .TRUE.
      IF(PRESENT(real_space_rotation)) real_space_rotation_arg = real_space_rotation

      IF(.NOT.ALLOCATED(mmpmatOut)) ALLOCATE(mmpmatOut,mold=mmpmat)
      mmpmatOut = mmpmat


      IF(spin_rotation_arg.AND.PRESENT(su)) THEN
         DO m = -lmaxU_const, lmaxU_const
            DO mp = -lmaxU_const, lmaxU_const

               d = cmplx_0
               d(1,1) = mmpmatOut(m,mp,1)
               d(2,2) = mmpmatOut(m,mp,2)
               IF(SIZE(mmpmat,3)>=3) THEN
                  d(2,1) = mmpmatOut(m,mp,3)
                  IF(SIZE(mmpmat,3)==3) THEN
                     d(1,2) = conjg(mmpmatOut(mp,m,3))
                  ELSE
                     d(1,2) = mmpmatOut(m,mp,4)
                  ENDIF
               ENDIF

               d = matmul(conjg(transpose(su)),d)
               d = matmul(d,su)

               mmpmatOut(m,mp,1) = d(1,1)
               mmpmatOut(m,mp,2) = d(2,2)
               IF(SIZE(mmpmat,3)>=3) THEN
                  mmpmatOut(m,mp,3) = d(2,1)
                  IF(SIZE(mmpmat,3)==4) THEN
                     mmpmatOut(mp,m,4) = d(1,2)
                  ENDIF
               ENDIF

            ENDDO
         ENDDO
      ENDIF

      IF(real_space_rotation_arg.AND.PRESENT(dwgn)) THEN
         DO ispin = 1, SIZE(mmpmat,3)
            mmpmatOut(:,:,ispin) = matmul(conjg(transpose(dwgn)),mmpmatOut(:,:,ispin))
            IF(PRESENT(dwgnp)) THEN
               mmpmatOut(:,:,ispin) = matmul(mmpmatOut(:,:,ispin),dwgnp)
            ELSE
               mmpmatOut(:,:,ispin) = matmul(mmpmatOut(:,:,ispin),dwgn)
            ENDIF
         ENDDO
      ENDIF

   END FUNCTION rotMMPmat_dwgn

   PURE FUNCTION rotMMPmat_angle(mmpmat,alpha,beta,gamma,l,lp,spin_rotation,real_space_rotation,inverse) Result(mmpmatOut)

      COMPLEX,           INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      REAL,              INTENT(IN)  :: alpha,beta,gamma !Euler angles
      INTEGER,           INTENT(IN)  :: l
      INTEGER, OPTIONAL, INTENT(IN)  :: lp
      LOGICAL, OPTIONAL, INTENT(IN)  :: spin_rotation
      LOGICAL, OPTIONAL, INTENT(IN)  :: real_space_rotation
      LOGICAL, OPTIONAL, INTENT(IN)  :: inverse

      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:,:)
      COMPLEX :: su(2,2), eia
      REAL    :: co_bh,si_bh, alphaArg, betaArg, gammaArg
      COMPLEX :: d(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX :: dp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      LOGICAL :: spin_rotation_arg, inverseArg

      IF(.NOT.ALLOCATED(mmpmatOut)) ALLOCATE(mmpmatOut,mold=mmpmat)
      mmpmatOut = mmpmat

      IF(ABS(alpha)<1e-10.AND.ABS(beta)<1e-10.AND.ABS(gamma)<1e-10) RETURN

      inverseArg = .FALSE.
      IF(PRESENT(inverse)) inverseArg = inverse

      IF(inverseArg) THEN
         alphaArg = -gamma; betaArg = -beta; gammaArg = -alpha
      ELSE
         alphaArg = alpha; betaArg = beta; gammaArg = gamma
      ENDIF

      d = d_wigner_mat(alphaArg,betaArg,gammaArg,l)
      IF(PRESENT(lp)) THEN
         dp = d_wigner_mat(alphaArg,betaArg,gammaArg,lp)
      ENDIF

      co_bh = cos(betaArg*0.5)
      si_bh = sin(betaArg*0.5)
      eia = exp( ImagUnit * alphaArg/2.0 )
      su(1,1) =  conjg(eia)*co_bh
      su(2,1) = -conjg(eia)*si_bh
      su(1,2) = eia*si_bh
      su(2,2) = eia*co_bh

      IF(PRESENT(lp)) THEN
         mmpmatOut = rotMMPmat_dwgn(mmpmat,d,dwgnp=dp,su=su,spin_rotation=spin_rotation,&
                                    real_space_rotation=real_space_rotation)
      ELSE
         mmpmatOut = rotMMPmat_dwgn(mmpmat,d,su=su,spin_rotation=spin_rotation,&
                                    real_space_rotation=real_space_rotation)
      ENDIF

   END FUNCTION rotMMPmat_angle

   PURE FUNCTION rotMMPmat_sym_op(mmpmat,sym,iop,l,lp,inverse,reciprocal,spin_rotation,real_space_rotation) Result(mmpmatOut)

      COMPLEX,           INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_sym),       INTENT(IN)  :: sym
      INTEGER,           INTENT(IN)  :: iop
      INTEGER,           INTENT(IN)  :: l
      INTEGER, OPTIONAL, INTENT(IN)  :: lp
      LOGICAL, OPTIONAL, INTENT(IN)  :: inverse, reciprocal, spin_rotation
      LOGICAL, OPTIONAL, INTENT(IN)  :: real_space_rotation

      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:,:)
      COMPLEX :: dwgn(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      COMPLEX :: dwgnp(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)
      INTEGER :: iopArg
      LOGICAL :: reciprocalArg, inverseArg

      IF(.NOT.ALLOCATED(mmpmatOut)) ALLOCATE(mmpmatOut,mold=mmpmat)
      mmpmatOut = mmpmat

      IF(iop==1) RETURN

      reciprocalArg = .FALSE.
      IF(PRESENT(reciprocal)) reciprocalArg = reciprocal

      inverseArg = .FALSE.
      IF(PRESENT(inverse)) inverseArg = inverse

      iopArg = iop
      IF(reciprocalArg) THEN
         IF(iop <= sym%nop) THEN
            iopArg = sym%invtab(iop)
         ELSE
            iopArg = sym%invtab(iop-sym%nop)
         ENDIF
      ELSE IF(inverseArg) THEN
         iopArg = sym%invtab(iop)
      ENDIF

      dwgn = sym%d_wgn(:,:,l,iopArg)
      IF(PRESENT(lp)) THEN
         dwgnp = sym%d_wgn(:,:,lp,iopArg)
      ELSE
         dwgnp = dwgn
      ENDIF

      IF(reciprocalArg) THEN
         IF(iop <= sym%nop) THEN
            dwgn = transpose(dwgn)
            dwgnp = transpose(dwgnp)
         ELSE
            dwgn = -transpose(dwgn)
            dwgnp = -transpose(dwgnp)
         ENDIF
      ENDIF

      mmpmatOut = rotMMPmat_dwgn(mmpmat,dwgn=dwgn,dwgnp=dwgnp,spin_rotation=spin_rotation,&
                                 real_space_rotation=real_space_rotation)

   END FUNCTION rotMMPmat_sym_op

   PURE FUNCTION rotMMPmat_sym_op_one_spin(mmpmat,sym,iop,l,lp,inverse,reciprocal) Result(mmpmatOut)

      COMPLEX,           INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:)
      TYPE(t_sym),       INTENT(IN)  :: sym
      INTEGER,           INTENT(IN)  :: iop
      INTEGER,           INTENT(IN)  :: l
      INTEGER, OPTIONAL, INTENT(IN)  :: lp
      LOGICAL, OPTIONAL, INTENT(IN)  :: inverse
      LOGICAL, OPTIONAL, INTENT(IN)  :: reciprocal

      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:), mmpmatOutsplit(:,:,:)
      COMPLEX :: mmpmatsplit(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,1)

      ALLOCATE(mmpmatOut,mold=mmpMat)
      mmpmatsplit(:,:,1) = mmpmat
      mmpmatOutsplit = rotMMPmat_sym_op(mmpmatsplit,sym,iop,l,lp=lp,inverse=inverse,reciprocal=reciprocal)
      mmpmatOut = mmpmatOutsplit(:,:,1)

   END FUNCTION rotMMPmat_sym_op_one_spin

   PURE FUNCTION rotMMPmat_angle_one_spin(mmpmat,alpha,beta,gamma,l,lp) Result(mmpmatOut)

      COMPLEX,           INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:)
      REAL,              INTENT(IN)  :: alpha,beta,gamma !Euler angles
      INTEGER,           INTENT(IN)  :: l
      INTEGER, OPTIONAL, INTENT(IN)  :: lp

      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:), mmpmatOutsplit(:,:,:)
      COMPLEX :: mmpmatsplit(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,1)

      ALLOCATE(mmpmatOut,mold=mmpMat)
      mmpmatsplit(:,:,1) = mmpmat
      mmpmatOutsplit = rotMMPmat_angle(mmpmatsplit,alpha,beta,gamma,l,lp=lp)
      mmpmatOut = mmpmatOutsplit(:,:,1)

   END FUNCTION rotMMPmat_angle_one_spin

   PURE FUNCTION rotMMPmat_dwgn_one_spin(mmpmat,dwgn,dwgnp) Result(mmpmatOut)

      COMPLEX,           INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:)
      COMPLEX,           INTENT(IN)  :: dwgn(-lmaxU_const:,-lmaxU_const:)
      COMPLEX, OPTIONAL, INTENT(IN)  :: dwgnp(-lmaxU_const:,-lmaxU_const:)

      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:), mmpmatOutsplit(:,:,:)
      COMPLEX :: mmpmatsplit(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,1)

      ALLOCATE(mmpmatOut,mold=mmpMat)
      mmpmatsplit(:,:,1) = mmpmat
      mmpmatOutsplit = rotMMPmat_dwgn(mmpmatsplit,dwgn=dwgn,dwgnp=dwgnp)
      mmpmatOut = mmpmatOutsplit(:,:,1)

   END FUNCTION rotMMPmat_dwgn_one_spin

   PURE FUNCTION d_wigner_mat(alpha, beta, gamma, l) RESULT(d)

      REAL,              INTENT(IN)  :: alpha,beta,gamma !Euler angles
      INTEGER,           INTENT(IN)  :: l

      COMPLEX :: d(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

      INTEGER :: m,mp,x_lo,x_up,x,e_c,e_s
      REAL    :: fac_l_m,fac_l_mp,fac_lmpx,fac_lmx,fac_x,fac_xmpm
      REAL    :: co_bh,si_bh,zaehler,nenner,cp,sp
      COMPLEX :: phase_g,phase_a,bas


      co_bh = cos(beta*0.5)
      si_bh = sin(beta*0.5)
      d = cmplx_0

      DO m = -l,l
         fac_l_m = fac(l+m) * fac(l-m)
         phase_g = exp( - ImagUnit * gamma * m )

         DO mp = -l,l
            fac_l_mp = fac(l+mp) * fac(l-mp)

            zaehler = sqrt( real(fac_l_m * fac_l_mp) )
            phase_a = exp( - ImagUnit * alpha * mp )
            x_lo = max(0, m-mp)
            x_up = min(l-mp, l+m)

            bas = zaehler * phase_a * phase_g
            DO x = x_lo,x_up
               fac_lmpx = fac(l-mp-x)
               fac_lmx  = fac(l+m-x)
               fac_x    = fac(x)
               fac_xmpm = fac(x+mp-m)
               nenner = fac_lmpx * fac_lmx * fac_x * fac_xmpm
               e_c = 2*l + m - mp - 2*x
               e_s = 2*x + mp - m
               IF (e_c.EQ.0) THEN
                  cp = 1.0
               ELSE
                  cp = co_bh ** e_c
               ENDIF
               IF (e_s.EQ.0) THEN
                  sp = 1.0
               ELSE
                  sp = si_bh ** e_s
               ENDIF
               d(m,mp) = d(m,mp) + bas * (-1)**x * cp * sp / nenner
            ENDDO

         ENDDO ! loop over mp
      ENDDO   ! loop over m
      DO m = -l,l
         DO mp = -l,l
            d( m,mp ) = d( m,mp ) * (-1)**(m-mp)
         ENDDO
      ENDDO

   END FUNCTION d_wigner_mat

   ELEMENTAL REAL FUNCTION  fac(n)

      INTEGER, INTENT (IN) :: n
      INTEGER :: i

      fac = 0
      IF (n.LT.0) RETURN
      fac = 1
      IF (n.EQ.0) RETURN
      DO i = 2,n
        fac = fac * i
      ENDDO

   END FUNCTION  fac


END MODULE m_rotMMPmat