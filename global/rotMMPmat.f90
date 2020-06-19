MODULE m_rotMMPmat

   USE m_types
   USE m_constants

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: rotMMPmat

   INTERFACE rotMMPmat
      PROCEDURE :: rotMMPmat_dwgn, rotMMPmat_angle, rotMMPmat_types
   END INTERFACE

   CONTAINS

   PURE FUNCTION rotMMPmat_dwgn(mmpmat,dwgn,su) Result(mmpmatOut)

      COMPLEX,           INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX, OPTIONAL, INTENT(IN)  :: dwgn(-lmaxU_const:,-lmaxU_const:)
      COMPLEX, OPTIONAL, INTENT(IN)  :: su(:,:)

      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:,:)

      COMPLEX :: d(2,2)
      INTEGER :: ispin,m,mp

      IF(.NOT.ALLOCATED(mmpmatOut)) ALLOCATE(mmpmatOut,mold=mmpmat)
      mmpmatOut = mmpmat

      IF(PRESENT(dwgn)) THEN
         DO ispin = 1, SIZE(mmpmat,3)
            mmpmatOut(:,:,ispin) = matmul(conjg(transpose(dwgn)),mmpmatOut(:,:,ispin))
            mmpmatOut(:,:,ispin) = matmul(mmpmatOut(:,:,ispin),dwgn)
         ENDDO
      ENDIF

      IF(SIZE(mmpmat,3)==3 .AND. PRESENT(su)) THEN
         DO m = -lmaxU_const, lmaxU_const
            DO mp = -lmaxU_const, lmaxU_const
               d(1,1) = mmpmatOut(m,mp,1)
               d(2,2) = mmpmatOut(m,mp,2)
               d(2,1) = mmpmatOut(m,mp,3)
               d(1,2) = conjg(mmpmatOut(mp,m,3))

               d = matmul(su,d)
               d = matmul(d,conjg(transpose(su)))

               mmpmatOut(m,mp,1) = d(1,1)
               mmpmatOut(m,mp,2) = d(2,2)
               mmpmatOut(m,mp,3) = d(2,1)
            ENDDO
         ENDDO
      ENDIF

      mmpmatOut = conjg(mmpmatOut)

   END FUNCTION rotMMPmat_dwgn

   PURE FUNCTION rotMMPmat_angle(mmpmat,alpha,beta,gamma,l) Result(mmpmatOut)

      COMPLEX,    INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:)
      REAL,       INTENT(IN)  :: alpha,beta,gamma !Euler angles
      INTEGER,    INTENT(IN)  :: l


      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:,:)
      COMPLEX :: su(2,2)
      INTEGER :: m,mp,x_lo,x_up,x,e_c,e_s
      REAL    :: fac_l_m,fac_l_mp,fac_lmpx,fac_lmx,fac_x,fac_xmpm
      REAL    :: co_bh,si_bh,zaehler,nenner,cp,sp
      COMPLEX :: phase_g,phase_a,bas,eia
      COMPLEX :: d(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const)

      IF(.NOT.ALLOCATED(mmpmatOut)) ALLOCATE(mmpmatOut,mold=mmpmat)
      mmpmatOut = mmpmat

      IF(ABS(alpha)<1e-10.AND.ABS(beta)<1e-10.AND.ABS(gamma)<1e-10) RETURN

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

     eia = exp( cmplx( 0.0 , alpha ) )
     su(1,1) =  conjg(eia)*co_bh
     su(2,1) = -conjg(eia)*si_bh
     su(1,2) = eia*si_bh
     su(2,2) = eia*co_bh

     mmpmatOut = rotMMPmat_dwgn(mmpmat,d,su)

   END FUNCTION rotMMPmat_angle

   PURE FUNCTION rotMMPmat_types(mmpmat,atoms) Result(mmpmatOut)

      COMPLEX,        INTENT(IN)  :: mmpmat(-lmaxU_const:,-lmaxU_const:,:,:)
      TYPE(t_atoms),  INTENT(IN)  :: atoms

      COMPLEX, ALLOCATABLE :: mmpmatOut(:,:,:,:)

      INTEGER :: i_u

      IF(.NOT.ALLOCATED(mmpmatOut)) ALLOCATE(mmpmatOut,mold=mmpmat)

      DO i_u = 1, atoms%n_u+atoms%n_hia
         mmpmatOut(:,:,i_u,:) = rotMMPmat_angle(mmpmat(:,:,i_u,:),atoms%lda_u(i_u)%phi,&
                                                atoms%lda_u(i_u)%theta,0.0,atoms%lda_u(i_u)%l)
      ENDDO

   END FUNCTION  rotMMPmat_types

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