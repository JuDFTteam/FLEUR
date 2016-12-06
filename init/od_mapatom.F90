!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_od_mapatom
      use m_juDFT
      CONTAINS 
      SUBROUTINE od_mapatom(oneD,atoms,sym,cell)

!      written by Y.Mokrousov in order to generate the arrays connected
!      to the operations, transforming atoms into each other,
!      for more details look in mapatom.F.    year 2004 

      USE m_types
      IMPLICIT NONE
      TYPE(t_oneD),INTENT(IN)    :: oneD
      TYPE(t_atoms),INTENT(INOUT):: atoms
      TYPE(t_sym),INTENT(INOUT)  :: sym
      TYPE(t_cell),INTENT(IN)    :: cell


      REAL ij,pps(3),norm,aamat(3,3) 
      INTEGER i,j,n1,k,n,n2,np1,na,ix,iy,iz,nat1,nat2,na2
      REAL mt(3,3),sum_tau_lat(3),sum_taual(3)
      REAL, PARAMETER :: del = 1.0e-4

      aamat=matmul(transpose(cell%amat),cell%amat)
     
      n1 = 1
      DO n = 1,atoms%ntype
        n2 = n1 + atoms%neq(n) - 1
        IF (atoms%neq(n).EQ.1) THEN
           atoms%ngopr(n2) = 1
           WRITE (6,FMT=8010) n2,n2,atoms%ngopr(n2)
           n1 = n1 + atoms%neq(n)
           CYCLE
        END IF
        DO  na = n1,n2
           DO np1 = 1,oneD%odd%nop
                  pps =matmul(sym%mrot(:,:,np1),atoms%taual(:,n1))
                  pps(3) = pps(3)+sym%tau(3,np1)/cell%amat(3,3)
                  IF (all(abs(atoms%taual(:,na)-pps(:)).LE.1.e-4 )) THEN
                      atoms%ngopr(na) = np1
                      WRITE (6,FMT=8010) na,n1,atoms%ngopr(na)
 8010                 FORMAT (5x,'atom',i3,' can be mapped into atom', i3,' through group  operation',i4)
                      CYCLE
                   END IF
           END DO
      ENDDO
        n1 = n1 + atoms%neq(n)
    ENDDO

!---> defining inverse operations for the Hamiltonian and forces
!     where we do not need to consider the translational part

      DO n1 = 1,oneD%odd%nop
         n2loop:DO n2 = 1,oneD%odd%nop
            mt=matmul(sym%mrot(:,:,n1),sym%mrot(:,:,n2))
            DO n = 1,oneD%odd%nop
               if (all(abs(mt(:,:) - sym%mrot(:,:,n)).LE.1.e-06)) THEN
                     sym%multab(n1,n2) = n
                     IF (n.EQ.1) sym%invtab(n1) = n2
                     cycle n2loop
               endif
            enddo
            WRITE (6,FMT=8050) n1,n2
 8050       FORMAT (' error - n1,n2=',2i3)
            CALL juDFT_error("mult",calledby ="od_mapatom")
        ENDDO n2loop
    ENDDO

 8060 FORMAT (1x,24i3)

      WRITE (6,FMT='(//," inverse operations",//)')

      DO n1 = 1,oneD%odd%nop
         WRITE (6,FMT=8060) n1,sym%invtab(n1)
      END DO

      DO na = 1,atoms%natd
         atoms%invsat(na) = 0
         sym%invsatnr(na) = 0
      END DO

      IF (oneD%odd%invs) THEN
         WRITE (6,FMT=*)
         nat1 = 1
         DO n = 1,atoms%ntype
            nat2 = nat1 + atoms%neq(n) - 1
            DO na = nat1,nat2 - 1
               IF (atoms%invsat(na).EQ.0) THEN
                  naloop:DO na2 = na + 1,nat2
                     DO i = 1,3
                        sum_taual(i) = atoms%taual(i,na) + atoms%taual(i,na2)
                     END DO
                     DO ix = -2,2
                       sum_tau_lat(1) = sum_taual(1) + real(ix)
                       DO iy = -2,2
                         sum_tau_lat(2) = sum_taual(2) + real(iy)
                         DO iz = -2,2
                           sum_tau_lat(3) = sum_taual(3) + real(iz)
                           norm = sqrt(dot_product(matmul(sum_tau_lat,aamat),sum_tau_lat))
                           IF (norm.LT.del) THEN
                              atoms%invsat(na) = 1
                              atoms%invsat(na2) = 2
                              sym%invsatnr(na)  = na2
                              sym%invsatnr(na2) = na
                              WRITE (6,FMT=9000) n,na,na2
                              CYCLE naloop
                           END IF
                        END DO
                      END DO
                    END DO
                  ENDDO naloop
               END IF
            END DO
            nat1 = nat1 + atoms%neq(n)
         END DO
      END IF
      WRITE (6,FMT=*) atoms%invsat
 9000 FORMAT ('atom type',i3,': atom',i3,' can be mapped into atom',i3, ' via 3d inversion')

      END SUBROUTINE od_mapatom
      END MODULE m_od_mapatom
