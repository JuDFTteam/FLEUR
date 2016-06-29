!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_apwsdim
      CONTAINS
      SUBROUTINE apws_dim(&
     &                    bkpt,cell,input,noco,oneD,&
     &                    nv,nv2,kq1d,kq2d,kq3d)
!
!*********************************************************************
!     determines dimensions of the lapw basis set with |k+G|<rkmax.
!     bkpt is the k-point given in internal units
!
!     dimensions kq(i)d for charge density FFT added.
!        s. bluegel, JRCAT, Feb. 97
!*********************************************************************
      USE m_boxdim
      USE m_ifft,     ONLY : ifft235
      USE m_types

      IMPLICIT NONE
      REAL,INTENT(IN)              :: bkpt(3)
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_oneD),INTENT(IN)      :: oneD
      INTEGER,INTENT(OUT) :: nv,nv2,kq1d,kq2d,kq3d   
      

      INTEGER j1,j2,j3,mk1,mk2,mk3,iofile,ksfft
      INTEGER ispin,nvh(2),nv2h(2)
      
      REAL arltv1,arltv2,arltv3,rkm,rk2,r2,s(3),gmaxp
! ..
!
!------->          ABBREVIATIONS
!

!   iofile      : device number for in and output
!   gmax        : cut-off wavevector for charge density
!   rkmax       : cut-off for |g+k|
!   gmaxp       : gmaxp = gmax/rkmax, ideal: gmaxp=2
!   arltv(i)    : length of reciprical lattice vector along
!                 direction (i)
!
!---> Determine rkmax box of size mk1, mk2, mk3,
!     for which |G(mk1,mk2,mk3) + (k1,k2,k3)| < rkmax
!
      CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)

!     (add 1+1 due to integer rounding, strange k_vector in BZ)
      mk1 = int(input%rkmax/arltv1) + 2
      mk2 = int(input%rkmax/arltv2) + 2
      mk3 = int(input%rkmax/arltv3) + 2

      rkm = input%rkmax
      rk2 = rkm*rkm
!---> obtain vectors
!---> in a spin-spiral calculation different basis sets are used for
!---> the two spin directions, because the cutoff radius is defined
!---> by |G + k +/- qss/2| < rkmax.
      DO ispin = 1,2
         nv = 0
         nv2 = 0
         DO j1 = -mk1,mk1
            s(1) = bkpt(1) + j1 + (2*ispin - 3)/2.0*noco%qss(1)
            DO j2 = -mk2,mk2
               s(2) = bkpt(2) + j2 + (2*ispin - 3)/2.0*noco%qss(2)
!--->          nv2 for films
               s(3) = 0.0
               !r2 = dotirp(s,s,cell%bbmat)
               r2 = dot_product(matmul(s,cell%bbmat),s)
               IF (r2.LE.rk2) nv2 = nv2 + 1
               DO j3 = -mk3,mk3
                  s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*noco%qss(3)
                  !r2 = dotirp(s,s,cell%bbmat)
                  r2 = dot_product(matmul(s,cell%bbmat),s)
                  IF (r2.LE.rk2) THEN
                     nv = nv + 1
                  END IF
               END DO
            END DO
         END DO
!-odim
         IF (oneD%odd%d1) THEN
           nv2 = 0
           s(1) = 0.0
           s(2) = 0.0
           DO j3 = -mk3,mk3
              s(3) = bkpt(3) + j3 + (2*ispin - 3)/2.0*noco%qss(3)
              !r2 = dotirp(s,s,cell%bbmat)
              r2 = dot_product(matmul(s,cell%bbmat),s)
               
              IF (r2.LE.rk2) THEN
                 nv2 = nv2 + 1
              END IF
           END DO
         END IF
!+odim
         nvh(ispin)  = nv
         nv2h(ispin) = nv2
      END DO
      nv  = max(nvh(1),nvh(2))
      nv2 = max(nv2h(1),nv2h(2))

!---> Determine the dimensions kq1d, kq2d, kq3d
!     of the dimension of the charge density fft-box
!     needed for the fast calculation of pw density
!     (add 1 due to integer rounding,
!      factor 2 due to positive domain)
!
      gmaxp = 2.0
      CALL boxdim(cell%bmat,arltv1,arltv2,arltv3)
!
      mk1 = int(gmaxp*input%rkmax/arltv1) + 1
      mk2 = int(gmaxp*input%rkmax/arltv2) + 1
      mk3 = int(gmaxp*input%rkmax/arltv3) + 1

!---> add + 1 in spin spiral calculation, to make sure that all G's are 
!---> still within the FFT-box after being shifted by the spin spiral
!---> q-vector.
      IF (noco%l_ss) THEN
         mk1 = mk1 + 1
         mk2 = mk2 + 1
         mk3 = mk3 + 1
      ENDIF         
!
      kq1d = 2*mk1
      kq2d = 2*mk2
      kq3d = 2*mk3
!
!---> fft's are usually fastest for low primes
!     (restrict kqid to: kqid=  (2**P) * (3**Q) * (5**R)
!
      ksfft = 1
!     ksfft=(0,1) : KEY OF SELECTING FFT-PRDOGRAM AND RADIX-TYPE
!                      0  PROGRAM, RADIX-2 ONLY
!                      1  PROGRAM, RADIX-2, RADIX-3,RADIX-5

      kq1d = ifft235(6,ksfft,kq1d,gmaxp)
      kq2d = ifft235(6,ksfft,kq2d,gmaxp)
      kq3d = ifft235(6,ksfft,kq3d,gmaxp)

      RETURN
      END SUBROUTINE
      END
