MODULE m_kk_cutoff

   USE m_trapz
   USE m_types
   USE m_juDFT
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE kk_cutoff(im,noco,l_mperp,l,jspins,eMesh,cutoff,scalingFactor)

      !This Subroutine determines the cutoff energy for the kramers-kronig-integration
      !This cutoff energy is defined so that the integral over the projDOS up to this cutoff
      !is equal to 2*(2l+1) (the number of states in the correlated shell) or not to small

      COMPLEX,             INTENT(IN)     :: im(:,-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_noco),        INTENT(IN)     :: noco
      LOGICAL,             INTENT(IN)     :: l_mperp
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: jspins
      REAL,                INTENT(IN)     :: eMesh(:)
      INTEGER,             INTENT(INOUT)  :: cutoff(:,:)
      REAL,                INTENT(INOUT)  :: scalingFactor(:)

      INTEGER :: m,ispin,spins_cut,ne
      REAL    :: lowerBound,upperBound,integral,n_states
      REAL    :: ec,del,eb,et
      REAL, ALLOCATABLE :: projDOS(:,:)

      ne = SIZE(eMesh)
      del = eMesh(2)-eMesh(1)
      eb = eMesh(1)
      et = eMesh(ne)


      !Calculate the trace over m,mp of the Imaginary Part matrix to obtain the projected DOS
      !n_f(e) = -1/pi * TR[Im(G_f(e))]
      ALLOCATE(projDOS(ne,jspins),source=0.0)
      DO ispin = 1, jspins
         DO m = -l , l
            projDOS(:,ispin) = projDOS(:,ispin) - 1/pi_const * REAL(im(:,m,m,ispin))
         ENDDO
      ENDDO


      spins_cut = MERGE(1,jspins,noco%l_noco.AND.l_mperp)
      n_states = (2*l+1) * MERGE(2.0,2.0/jspins,noco%l_noco.AND.l_mperp)

      cutoff(:,1) = 1   !we don't modify the lower bound
      cutoff(:,2) = ne
      scalingFactor(:) = 1.0

      DO ispin = 1, spins_cut

         !----------------------------------------
         !Check the integral up to the hard cutoff
         !----------------------------------------
         IF(spins_cut.EQ.1 .AND.jspins.EQ.2) projDOS(:,1) = projDOS(:,1) + projDOS(:,2)

         !Initial complete integral
         integral =  trapz(projDOS(:,ispin),del,ne)

#ifdef CPP_DEBUG
         WRITE(*,*) "Integral over DOS: ", integral
#endif

         IF(integral.LT.n_states) THEN
            !If we are calculating the greens function for a d-band this is expected to happen
            IF(l.EQ.2) THEN
               scalingFactor(ispin) = n_states/integral

#ifdef CPP_DEBUG
               WRITE(*,9000) l,ispin,scalingFactor(ispin)
9000           FORMAT("Scaling the DOS for l=",I1," and spin ",I1,"   factor: ",f14.8)
#endif
               IF(scalingFactor(ispin).GT.1.25) CALL juDFT_warn("scaling factor >1.25 -> increase elup(<1htr) or numbands",calledby="kk_cutoff")
            ELSE IF(integral.LT.n_states-0.1) THEN
               ! If the integral is to small we terminate here to avoid problems
               CALL juDFT_warn("Integral over DOS too small for f -> increase elup(<1htr) or numbands", calledby="kk_cutoff")
            ENDIF
         ELSE

            !IF the integral is bigger than 2l+1, search for the cutoff using the bisection method
            lowerBound = eb
            upperBound = et

            DO WHILE(ABS(upperBound-lowerBound).GT.del/2.0)

               ec = (lowerBound+upperBound)/2.0
               cutoff(ispin,2) = INT((ec-eb)/del)+1

               !Integrate the DOS up to the cutoff
               integral =  trapz(projDOS(:,ispin),del,cutoff(ispin,2))

               IF(ABS(integral-n_states).LT.1e-12) THEN
                  EXIT
               ELSE IF(integral.LT.n_states) THEN
                  !integral to small -> choose the right interval
                  lowerBound = ec
               ELSE IF(integral.GT.n_states) THEN
                  !integral to big   -> choose the left interval
                  upperBound = ec
               END IF

            ENDDO

#ifdef CPP_DEBUG
            WRITE(*,*) "CALCULATED CUTOFF: ", cutoff(ispin,2)
            WRITE(*,*) "CORRESPONDING ENERGY", ec
            WRITE(*,*) "INTEGRAL OVER projDOS with cutoff: ", integral
#endif
            !Copy cutoff to second spin if only one was calculated
            IF(spins_cut.EQ.1 .AND. jspins.EQ.2) cutoff(2,2) = cutoff(1,2)
         ENDIF
      ENDDO

   END SUBROUTINE kk_cutoff

   SUBROUTINE kk_cutoffRadial(uu,ud,du,dd,noco,scalarGF,l_mperp,&
                              l,input,eMesh,cutoff,scalingFactor)

      COMPLEX,                   INTENT(IN)     :: uu(:,-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX,                   INTENT(IN)     :: ud(:,-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX,                   INTENT(IN)     :: du(:,-lmaxU_const:,-lmaxU_const:,:)
      COMPLEX,                   INTENT(IN)     :: dd(:,-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_noco),              INTENT(IN)     :: noco
      TYPE(t_scalarGF),          INTENT(IN)     :: scalarGF
      LOGICAL,                   INTENT(IN)     :: l_mperp
      INTEGER,                   INTENT(IN)     :: l
      TYPE(t_input),             INTENT(IN)     :: input
      REAL,                      INTENT(IN)     :: eMesh(:)
      INTEGER,                   INTENT(INOUT)  :: cutoff(:,:)
      REAL,                      INTENT(INOUT)  :: scalingFactor(:)

      COMPLEX, ALLOCATABLE :: im(:,:,:,:)

      INTEGER :: jspin,m,mp,spin1,spin2

      !calculate the spherical average from the original greens function
      ALLOCATE(im(SIZE(uu,1),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,SIZE(uu,4)),source=cmplx_0)
      DO jspin = 1, SIZE(im,4)
         IF(jspin < 3) THEN
            spin1 = jspin
            spin2 = jspin
         ELSE
            spin1 = 2
            spin2 = 1
         ENDIF
         !$OMP parallel do default(none) &
         !$OMP shared(scalarGF,jspin,spin1,spin2,l,im,uu,ud,du,dd) &
         !$OMP private(m,mp) collapse(2)
         DO m = -l,l
            DO mp = -l,l
               im(:,m,mp,jspin) =     uu(:,m,mp,jspin) * scalarGF%uun(spin1,spin2) &
                                    + ud(:,m,mp,jspin) * scalarGF%udn(spin1,spin2) &
                                    + du(:,m,mp,jspin) * scalarGF%dun(spin1,spin2) &
                                    + dd(:,m,mp,jspin) * scalarGF%ddn(spin1,spin2)
            ENDDO
         ENDDO
         !$OMP end parallel do
      ENDDO

      CALL kk_cutoff(im,noco,l_mperp,l,input%jspins,eMesh,cutoff,scalingFactor)


   END SUBROUTINE kk_cutoffRadial

END MODULE m_kk_cutoff
