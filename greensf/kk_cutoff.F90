MODULE m_kk_cutoff

   USE m_kkintgr
   USE m_types
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE kk_cutoff(im,noco,l,jspins,ne,del,e_bot,e_top,cutoff)

      !This Subroutine determines the cutoff energy for the kramers-kronig-integration
      !This cutoff energy is defined so that the integral over the fDOS up to this cutoff
      !is equal to 2*(2l+1) (the number of states in the correlated shell) or not to small

      REAL,                INTENT(INOUT)  :: im(:,-lmaxU_const:,-lmaxU_const:,:)
      TYPE(t_noco),        INTENT(IN)     :: noco
      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: jspins
      INTEGER,             INTENT(IN)     :: ne
      REAL,                INTENT(IN)     :: del
      REAL,                INTENT(IN)     :: e_bot
      REAL,                INTENT(IN)     :: e_top
      INTEGER,             INTENT(INOUT)  :: cutoff(:,:)

      CHARACTER(len=5) :: filename
      INTEGER :: i,m,n_c,ispin,spins_cut
      REAL    :: lowerBound,upperBound,integral,n_states,scale,e_cut
      REAL    :: fDOS(ne,jspins)

      fDOS = 0.0

      !Calculate the trace over m,mp of the Greens-function matrix to obtain the fDOS
      !n_f(e) = -1/pi * TR[Im(G_f(e))]
      DO ispin = 1, jspins
         DO m = -l , l
            DO i = 1, ne
               fDOS(i,ispin) = fDOS(i,ispin) + im(i,m,m,ispin)
            ENDDO
         ENDDO
      ENDDO

      fDOS = -1/pi_const*fDOS

!#ifdef CPP_DEBUG
      !DO ispin = 1, jspins
      !   WRITE(filename,9010) ispin
!9010     FORMAT("fDOS",I1)
      !   OPEN(unit=1337,file=filename,status="replace")
      !   DO i = 1, ne
      !      WRITE(1337,"(2f14.8)") (i-1)*del+e_bot,fDOS(i,ispin)
      !   ENDDO
      !   CLOSE(unit=1337)
      !ENDDO
!#endif

      spins_cut = MERGE(1,jspins,noco%l_noco.AND.noco%l_mperp)
      n_states = (2*l+1) * MERGE(2.0,2.0/jspins,noco%l_noco.AND.noco%l_mperp)

      cutoff(:,1) = 1   !we don't modify the lower bound
      cutoff(:,2) = ne

      DO ispin = 1, spins_cut

         !----------------------------------------
         !Check the integral up to the hard cutoff
         !----------------------------------------
         IF(spins_cut.EQ.1.AND.jspins.EQ.2) fDOS(:,1) = fDOS(:,1) + fDOS(:,2)
         integral =  trapz(fDOS(:,ispin),del,ne)

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
         ELSE IF((integral.GT.n_states).AND.((integral-n_states).GT.0.00001)) THEN
            !IF the integral is bigger than 2l+1, search for the cutoff using the bisection method

            lowerBound = e_bot
            upperBound = e_top

            DO WHILE(upperBound-lowerBound.GT.del)

               e_cut = (lowerBound+upperBound)/2.0
               cutoff(ispin,2) = INT((e_cut-e_bot)/del)+1
               integral =  trapz(fDOS(:,ispin),del,cutoff(ispin,2))

               IF(integral.LT.n_states) THEN
                  !integral to small -> choose the right interval
                  lowerBound = e_cut
               ELSE IF(integral.GT.n_states) THEN
                  !integral to big   -> choose the left interval
                  upperBound = e_cut
               END IF

            ENDDO

#ifdef CPP_DEBUG
            WRITE(*,*) "CALCULATED CUTOFF: ", cutoff(ispin,2)
            WRITE(*,*) "CORRESPONDING ENERGY", e_cut
            WRITE(*,*) "INTEGRAL OVER fDOS with cutoff: ", integral
#endif
            IF(spins_cut.EQ.1.AND.jspins.EQ.2) cutoff(2,2) = cutoff(1,2)
         ENDIF
      ENDDO

   END SUBROUTINE kk_cutoff
END MODULE m_kk_cutoff
