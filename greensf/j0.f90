MODULE m_j0

   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_kkintgr
   USE m_ind_greensf

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE eff_excinteraction(g0,atoms,input,ef,g0Coeffs)

      TYPE(t_greensf),        INTENT(IN)  :: g0
      TYPE(t_atoms),          INTENT(IN)  :: atoms
      TYPE(t_greensfCoeffs),  INTENT(IN)  :: g0Coeffs !For determining the onsite exchange splitting from the difference in the COM of the bands
      TYPE(t_input),          INTENT(IN)  :: input
      REAL,                   INTENT(IN)  :: ef

      COMPLEX integrand, sumup, sumdwn, sumupdwn
      INTEGER i,iz,m,l,mp,ispin,n,i_gf,matsize,ipm,ie,n_cut
      REAL beta,j0(0:lmaxU_const),exc_split(0:lmaxU_const),tmp,avgexc
      INTEGER nType,l_min,l_max,i_j0
      CHARACTER(len=30) :: filename

      LOGICAL l_matinv
      TYPE(t_mat) :: calcup,calcdwn
      TYPE(t_mat) :: delta
      TYPE(t_mat) :: calc

      REAL :: int_com(g0Coeffs%ne,input%jspins),int_norm(g0Coeffs%ne,input%jspins)
      REAL, PARAMETER    :: boltzmannConst = 3.1668114e-6 ! value is given in Hartree/Kelvin

      l_matinv = .FALSE. !Determines how the onsite exchange splitting is calculated

      DO i_j0 = 1, atoms%n_j0

         nType = atoms%j0(i_j0)%atomType
         l_min = atoms%j0(i_j0)%l_min
         l_max = atoms%j0(i_j0)%l_max
         WRITE(6,9010) nType
         WRITE(6,9020) l_min,l_max,atoms%j0(i_j0)%l_avgexc
         j0 = 0.0
         !Calculate the onsite exchange splitting by determining the difference in the center of mass
         !of the bands under consideration
         !The cutoffs were determined so that both the integral over the DOS for spin up/down is equal to 2*l+1
         IF(.NOT.l_matinv) THEN
            exc_split = 0.0
            DO l = l_min, l_max
               i_gf = ind_greensf(atoms,l,nType)
               !-------------------------------------------------
               ! Evaluate center of mass of the bands in question
               ! and tkae the differrence between spin up/down
               ! We use the same cutoffs as for the kk-integration
               !-------------------------------------------------
               ! E_COM = 1/(int dE N_l(E)) * int dE E * N_l(E)
               !-------------------------------------------------
               DO ispin = 1, input%jspins
                  int_com = 0.0
                  int_norm = 0.0
                  n_cut = g0Coeffs%kkintgr_cutoff(i_gf,ispin,2)
                  DO ie = 1, n_cut
                     DO m = -l, l
                        int_com(ie,ispin) = int_com(ie,ispin) + ((ie-1)*g0Coeffs%del+g0Coeffs%e_bot)&
                                                                *g0Coeffs%projdos(ie,i_gf,m,m,ispin)
                        int_norm(ie,ispin) = int_norm(ie,ispin) + g0Coeffs%projdos(ie,i_gf,m,m,ispin)
                     ENDDO
                  ENDDO
                  exc_split(l) = exc_split(l) + (-1)**(ispin) * 1.0/(trapz(int_norm(:n_cut,ispin),g0Coeffs%del,n_cut)) * trapz(int_com(:n_cut,ispin),g0Coeffs%del,n_cut) 
               ENDDO
            ENDDO
            IF(atoms%j0(i_j0)%l_avgexc) THEN
               avgexc = SUM(exc_split(l_min:l_max))/(l_max-l_min+1)
               DO l = l_min,l_max
                  exc_split(l) = avgexc
               ENDDO
            ENDIF
         ENDIF
      DO l = l_min,l_max
         WRITE(filename,9060) i_j0, l
         IF(atoms%j0(i_j0)%l_eDependence) OPEN(unit=1337,file=filename,status="replace")
         WRITE(6,9030) l,exc_split(l)*hartree_to_ev_const

         matsize = (2*l+1)
         CALL delta%init(.false.,matsize,matsize)
         CALL calc%init(.false.,matsize,matsize)
         DO i = 1, matsize
            delta%data_c(i,i) = exc_split(l)
         ENDDO

         DO iz = 1, g0%nz
            !
            !calculate the onsite exchange matrix if we use matrix inversion
            !
            IF(l_matinv) THEN
               !---------------------------------------------
               !\Delta = (G_up)^-1-(G_down)^-1
               !---------------------------------------------
               !Symmetrize the green's function for up/down
               !spin with respect to the complex plane
               !Here we assume that the onsite Hamiltonian
               !is real
               !---------------------------------------------
               !G^(up/down)^-1 = 1/2 * (G+^(up/down)^-1 + G-^(up/down)^-1)
               !---------------------------------------------
               delta%data_c = 0.0
               DO ispin = 1, input%jspins
                  DO ipm = 1, 2
                     CALL g0%get_gf(calc,atoms,input,iz,l,nType,ipm.EQ.2,spin=ispin)
                     CALL calc%inverse()
                     delta%data_c = delta%data_c + 1/2.0 * (-1)**(ispin-1) * calc%data_c
                  ENDDO
               ENDDO
            ENDIF
            !
            !  Tr_L[\Delta (G_up-G_down) + \Delta G_up \Delta G_down]
            !
            ! calculated for G^+/- and then substract to obtain imaginary part
            integrand = 0.0
            !These three are only for visualization of the individual terms
            sumup = 0.0
            sumdwn = 0.0
            sumupdwn = 0.0
            DO ipm = 1, 2
               CALL g0%get_gf(calcup,atoms,input,iz,l,nType,ipm.EQ.2,spin=1)
               CALL g0%get_gf(calcdwn,atoms,input,iz,l,nType,ipm.EQ.2,spin=2)
               calcup%data_c  = matmul(delta%data_c,calcup%data_c)
               calcdwn%data_c = matmul(delta%data_c,calcdwn%data_c)
               calc%data_c    = matmul(calcup%data_c,calcdwn%data_c)
               !Calculate the trace
               DO i = 1,matsize

                  !Calculate all three terms explicitly (only for visualization)
                  sumup     = sumup     + (-1)**(ipm) * 1./(2.0*fpi_const) * hartree_to_ev_const * calcup%data_c(i,i) &
                                                      * MERGE(g0%de(iz),conjg(g0%de(iz)),ipm.EQ.1)
                  sumdwn    = sumdwn    + (-1)**(ipm) * 1./(2.0*fpi_const) * hartree_to_ev_const * calcdwn%data_c(i,i) &
                                                      * MERGE(g0%de(iz),conjg(g0%de(iz)),ipm.EQ.1)
                  sumupdwn  = sumupdwn  + (-1)**(ipm) * 1./(2.0*fpi_const) * hartree_to_ev_const * calc%data_c(i,i) &
                                                      * MERGE(g0%de(iz),conjg(g0%de(iz)),ipm.EQ.1)

                  integrand = integrand + (-1)**(ipm-1) * (calcup%data_c(i,i)-calcdwn%data_c(i,i)+calc%data_c(i,i)) &
                                                      * MERGE(g0%de(iz),conjg(g0%de(iz)),ipm.EQ.1)
               ENDDO
               CALL calcup%free()
               CALL calcdwn%free()
            ENDDO

            IF(ABS(REAL(integrand)) > 1e-5) THEN
               CALL juDFT_error("integrand still has a Real part", calledby="eff_excinteraction")
            ENDIF
            j0(l) = j0(l) + AIMAG(integrand)

            IF(atoms%j0(i_j0)%l_eDependence) THEN
            WRITE(1337,"(5f14.8)") REAL(g0%e(iz)-ef)*hartree_to_ev_const, -1/(2.0*fpi_const)*hartree_to_ev_const *j0(l),&
                                   AIMAG(sumup),AIMAG(sumdwn),AIMAG(sumupdwn)
            ENDIF


            IF(l_matinv) CALL calc%free()

         ENDDO
         j0(l) = -1/(2.0*fpi_const)*hartree_to_ev_const * j0(l)
         WRITE(6,9040) l,j0(l),ABS(j0(l))*2/3*1/(boltzmannConst*hartree_to_ev_const)


         IF(atoms%j0(i_j0)%l_eDependence) CLOSE(unit=1337)
         CALL delta%free()
         CALL calc%free()
         ENDDO
         WRITE(6,9050) SUM(j0(l_min:l_max)), ABS(SUM(j0(l_min:l_max)))*2/3 * 1/(boltzmannConst*hartree_to_ev_const)
      ENDDO

9000  FORMAT("Effective Magnetic Exchange Interaction J0 (Compare Condens. Matter 26, 476003 (2014) EQ.1)")
9010  FORMAT("J0 calculation for atom ",I3)
9020  FORMAT("lmin: ", I3, "  lmax: ",I3," Averaged Exchange Splitting: ", L1)
9030  FORMAT("Onsite Exchange Splitting for l=",I1" : ", f14.8," eV")
9040  FORMAT("J0 (l=", I3, "): " f14.8 " eV", "  Tc: " f14.8," K")
9050  FORMAT("J0 : " f14.8 " eV", "  Tc: " f14.8," K")
9060  FORMAT("j0_",I1,".",I1)

   END SUBROUTINE eff_excinteraction
END MODULE m_j0