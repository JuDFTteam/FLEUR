MODULE m_kk_cutoff


   LOGICAL, PARAMETER :: l_debug=.TRUE.

   CONTAINS

   SUBROUTINE kk_cutoff(im,atoms,noco,l,jspins,ne,del,e_bot,e_top,cutoff)
      !This Subroutine determines the cutoff energy for the kramers-kronig-integration
      !This cutoff energy is defined so that the integral over the fDOS up to this cutoff 
      !is equal to 2*(2l+1) (the number of states in the correlated shell) or not to small

      USE m_kkintgr
      USE m_types

      IMPLICIT NONE

      REAL,                INTENT(INOUT)  :: im(ne,-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,jspins)
      TYPE(t_atoms),       INTENT(IN)     :: atoms
      TYPE(t_noco),        INTENT(IN)     :: noco

      INTEGER,             INTENT(IN)     :: l
      INTEGER,             INTENT(IN)     :: jspins
      INTEGER,             INTENT(IN)     :: ne
      REAL,                INTENT(IN)     :: del
      REAL,                INTENT(IN)     :: e_bot
      REAL,                INTENT(IN)     :: e_top

      INTEGER,             INTENT(OUT)    :: cutoff(jspins,2)


      INTEGER i,m,n_c,ispin,spins_cut
      REAL integral
      REAL a,b, n_states, scale
      CHARACTER(len=5) :: filename

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

      fDOS = -1/pi_const*fDOS

      IF(l_debug) THEN
         DO ispin = 1, jspins
            WRITE(filename,9010) ispin
            OPEN(unit=1337,file=filename,status="replace")
            DO i = 1, ne
               WRITE(1337,"(2f14.8)") (i-1)*del+e_bot,fDOS(i,ispin)
            ENDDO
            CLOSE(unit=1337)
         ENDDO
      ENDIF

      spins_cut = MERGE(1,jspins,noco%l_soc.OR.noco%l_noco)
      n_states = (2*l+1) * MERGE(2.0,2.0/jspins,noco%l_soc.OR.noco%l_noco)

      spins_cut=jspins
      n_states = (2*l+1) * 2.0/jspins

      cutoff(:,1) = 1   !we don't modify the lower bound
      cutoff(:,2) = ne

      DO ispin = 1, spins_cut
         IF(spins_cut.EQ.1.AND.jspins.EQ.2) fDOS(:,1) = fDOS(:,1) + fDOS(:,2) 
         integral =  trapz(fDOS(1:ne,ispin), del, ne)      
         IF(l_debug) WRITE(*,*) "Integral over DOS: ", integral
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
            IF(spins_cut.EQ.1.AND.jspins.EQ.2) cutoff(2,2) = cutoff(1,2)
         ENDIF
      ENDDO

9000  FORMAT("Scaling the DOS for l=",I1," and spin ",I1,"   factor: ",f14.8)
9010  FORMAT("fDOS",I1)

   END SUBROUTINE kk_cutoff
END MODULE m_kk_cutoff