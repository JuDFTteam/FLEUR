MODULE m_kk_cutoff

   USE m_trapz
   USE m_types
   USE m_juDFT
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE kk_cutoff(im,noco,l_mperp,l,jspins,eMesh,cutoff)

      !This Subroutine determines the cutoff energy for the kramers-kronig-integration
      !This cutoff energy is defined so that the integral over the projDOS up to this cutoff
      !is equal to 2*(2l+1) (the number of states in the correlated shell) or not to small

      REAL,                INTENT(INOUT)  :: im(:,-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_noco),        INTENT(IN)     :: noco
      LOGICAL,             INTENT(IN)     :: l_mperp
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: jspins
      REAL,                INTENT(IN)     :: eMesh(:)
      INTEGER,             INTENT(INOUT)  :: cutoff(:,:)

      INTEGER :: m,ispin,spins_cut,ne
      REAL    :: lowerBound,upperBound,integral,n_states,scale
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
            projDOS(:,ispin) = projDOS(:,ispin) - 1/pi_const * im(:,m,m,ispin)
         ENDDO
      ENDDO


      spins_cut = MERGE(1,jspins,noco%l_noco.AND.l_mperp)
      n_states = (2*l+1) * MERGE(2.0,2.0/jspins,noco%l_noco.AND.l_mperp)

      cutoff(:,1) = 1   !we don't modify the lower bound
      cutoff(:,2) = ne

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
               scale = (2*l+1)/integral

#ifdef CPP_DEBUG
               WRITE(*,9000) l,ispin,scale
9000           FORMAT("Scaling the DOS for l=",I1," and spin ",I1,"   factor: ",f14.8)
#endif

               IF(scale.GT.1.25) CALL juDFT_warn("scaling factor >1.25 -> increase elup(<1htr) or numbands",calledby="kk_cutoff")
               im(:,-l:l,-l:l,ispin) = scale * im(:,-l:l,-l:l,ispin)
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
END MODULE m_kk_cutoff
