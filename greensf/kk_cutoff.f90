MODULE m_kk_cutoff


   LOGICAL, PARAMETER :: l_debug=.TRUE.

   CONTAINS

   SUBROUTINE kk_cutoff(im,atoms,l,jspins,ne,del,e_bot,e_top,cutoff)
      !This Subroutine determines the cutoff energy for the kramers-kronig-integration
      !This cutoff energy is defined so that the integral over the fDOS up to this cutoff 
      !is equal to 2*(2l+1) (the number of states in the correlated shell) or not to small

      USE m_kkintgr
      USE m_types

      IMPLICIT NONE

      REAL,                INTENT(INOUT)  :: im(ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
      TYPE(t_atoms),       INTENT(IN)     :: atoms

      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: jspins
      INTEGER,             INTENT(IN)     :: ne
      REAL,                INTENT(IN)     :: del
      REAL,                INTENT(IN)     :: e_bot
      REAL,                INTENT(IN)     :: e_top

      INTEGER,             INTENT(OUT)    :: cutoff(jspins,2)


      INTEGER i,m,n_c,ispin
      REAL integral
      REAL a,b, n_states, scale

      REAL :: fDOS(ne,jspins)

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
      DO ispin = 1, jspins
         fDOS(:,ispin) = -1/pi_const*fDOS(:,ispin)

         integral =  trapz(fDOS(1:ne,ispin), del, ne)

         n_states = (2*l+1)
      
         IF(l_debug) WRITE(*,*) "Integral over DOS: ", integral

         cutoff(ispin,1) = 1   !we don't modify the lower bound
         cutoff(ispin,2) = ne

         IF(integral.LT.n_states) THEN
            !If we are calculating the greens function for a d-band this is expected
            IF(l.EQ.2) THEN
               scale = (2*l+1)/integral
               IF(l_debug) WRITE(*,9000) l,ispin,scale 
               IF(scale.GT.1.25) CALL juDFT_warn("scaling factor >1.25 -> increase elup(<1htr) or numbands",calledby="kk_cutoff")
               im(:,-l:l,-l:l,ispin) = scale * im(:,-l:l,-l:l,ispin)
            ELSE IF(integral.LT.n_states-0.1) THEN
            ! If the integral is to small we stop here to avoid problems
               CALL juDFT_warn("Integral over DOS too small for f -> increase elup(<1htr) or numbands", calledby="kk_cutoff")
            ENDIF
         ELSE IF((integral.GT.n_states).AND.((integral-n_states).GT.0.001)) THEN
            !IF the integral is bigger than 2l+1, search for the cutoff using the bisection method   

            a = e_bot
            b = e_top

            DO

               cutoff(ispin,2) = INT(((a+b)/2.0-e_bot)/del)+1
               integral =  trapz(fDOS(1:cutoff(ispin,2),ispin),del,cutoff(ispin,2))

               IF((ABS(integral-n_states).LT.0.001).OR.(ABS(a-b)/2.0.LT.del)) THEN
                  !The integral is inside the desired accuracy
                  EXIT
               ELSE IF((integral-n_states).LT.0) THEN
                  !integral to small -> choose the right interval
                  a = (a+b)/2.0
               ELSE IF((integral-n_states).GT.0) THEN
                  !integral to big   -> choose the left interval
                  b = (a+b)/2.0
               END IF

            ENDDO

            IF(l_debug) THEN
               WRITE(*,*) "CALCULATED CUTOFF: ", cutoff(ispin,2)
               WRITE(*,*) "INTEGRAL OVER fDOS with cutoff: ", integral
            ENDIF
         ENDIF
      ENDDO

9000  FORMAT("Scaling the DOS for l=",I1," and spin ",I1,"   factor: ",f14.8)

   END SUBROUTINE kk_cutoff
END MODULE m_kk_cutoff