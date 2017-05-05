      MODULE m_mapatom
      use m_juDFT
!*******************************************************************
!     determines the group operation which maps the representive
!     atom into its equivalent atoms     c.l.fu
!*******************************************************************
      CONTAINS
      SUBROUTINE mapatom(sym,atoms,cell,input,noco)
!
!     if (l_f) setup multab,invtab,invarop,invarind for force_a12 & 21
!***********************************************************************
! the contribution to the hamiltonian and to the overlap matrix in a
! system with inversion symmetry from one muffin tin is the complex
! conjugate of the contribution from the "invers" muffin tin. this fact
! can be exploited to save cpu-time in hssphn. Thus, it is nessessary to
! know whether an atom can be mapped onto an equivalent atom via 3d
! inversion. where both atoms have to belong to the same unit cell, i.e.
! the are not related to each other by a lattice translation. therefore,
! an array invsatom is set up.
! invsatom(natom) =
! 0 if the atom cannot be mapped onto an eqivalent atom via inversion
! 1 if the atom can be mapped onto an eqivalent atom via inversion, and
!   has a smaller atom index than the related atom
! 2 if the atom can be mapped onto an eqivalent atom via inversion, and
!   has a bigger atom index than the related atom
! p.kurz aug. 1996
!***********************************************************************
!
      USE m_socsym
      USE m_types
      IMPLICIT NONE
      TYPE(t_sym),INTENT(INOUT)   :: sym
      TYPE(t_atoms),INTENT(INOUT) :: atoms
      TYPE(t_cell),INTENT(IN)  :: cell
      TYPE(t_input),INTENT(IN) :: input
      TYPE(t_noco),INTENT(IN)  :: noco

!     .. Local Scalars ..
      REAL s3,norm
      INTEGER i,icount,j,j1,j2,j3,jop,n,na,nat1,nat2,nb,na_r
      INTEGER k,ij,n1,n2,ix,iy,iz,na2
      REAL, PARAMETER :: del = 1.0e-4
!     ..
!     .. Local Arrays ..
      INTEGER mt(3,3),mp(3,3)
      REAL aamat(3,3),sum_tau_lat(3),sum_taual(3)
      REAL gam(3),gaminv(3),gamr(3),sr(3),ttau(3)
      LOGICAL error(sym%nop)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC real,sqrt
!     ..
      
    !  CALL dotset(&
    ! &            cell,&
    ! &            aamat, )
      aamat=matmul(transpose(cell%amat),cell%amat)
    
      IF (noco%l_soc) THEN  ! check once more here...
        CALL soc_sym(&
     &               sym%nop,sym%mrot,noco%theta,noco%phi,cell%amat,&
     &               error)
      ELSE
        error(:) = .false.
      ENDIF
                               
      WRITE (6,FMT=8000)
 8000 FORMAT (/,/,5x,'group operations on equivalent atoms:')
      nat1 = 1
      DO n = 1,atoms%ntype
         nat2 = nat1 + atoms%neq(n) - 1
         atoms%ngopr(nat1) = 1
!+gu
         na_r = nat1
         DO na = nat1,nat2
            IF (atoms%ntypsy(na).NE.atoms%ntypsy(na_r)) na_r = na
!-gu
            DO i = 1,3
               gam(i) = atoms%taual(i,na)
            END DO
            sym%invarind(na) = 0
            icount = 0
            DO  jop = 1,sym%nop
               DO i = 1,3
                  gamr(i) = 0.
                  DO j = 1,3
                     gamr(i) = gamr(i) + sym%mrot(i,j,jop)*gam(j)
                  END DO
                  gamr(i) = gamr(i) + sym%tau(i,jop)
               END DO
               DO i = 1,3
                  gaminv(i) = gamr(i) - atoms%taual(i,na)
                  gamr(i)   = gamr(i) - atoms%taual(i,nat1) ! cf local_sym
               END DO
               IF (icount.EQ.0) THEN
                  DO j3 = -2,2
                     sr(3) = gamr(3) + real(j3)
                     DO j2 = -2,2
                        sr(2) = gamr(2) + real(j2)
                        DO j1 = -2,2
                           sr(1) = gamr(1) + real(j1)
                           s3 = sqrt(dot_product(matmul(sr,aamat),sr))
                           IF ((s3.LT.del).AND.(.not.error(jop))) THEN
                              icount = icount + 1
                              atoms%ngopr(na) = jop
                           END IF
                        END DO
                     END DO
                  END DO
               END IF
!
! search for operations which leave taual invariant
!
               IF (input%l_f.OR.(atoms%n_u.GT.0)) THEN 
                  DO j3 = -2,2
                     sr(3) = gaminv(3) + real(j3)
                     DO j2 = -2,2
                        sr(2) = gaminv(2) + real(j2)
                        DO j1 = -2,2
                           sr(1) = gaminv(1) + real(j1)
                           s3 = sqrt(dot_product(matmul(sr,aamat),sr))
                           IF (s3.LT.del) THEN
                              sym%invarind(na) = sym%invarind(na) + 1
                              sym%invarop(na,sym%invarind(na)) = jop
                           END IF
                        END DO
                     END DO
                  END DO
               ENDIF
!
! end of operations
          ENDDO
            IF (icount.LE.0) THEN
             write(6,*) "Mapping failed for atom:",nat1
             write(6,*) "No of symmetries tested:",sym%nop
             CALL juDFT_error("mapatom",calledby="mapatom")
           ENDIF
            WRITE (6,FMT=8010) nat1,na,atoms%ngopr(na)
 8010       FORMAT (5x,'atom',i3,' can be mapped into atom',i3,&
     &             ' through group  operation',i4)
!
! end of equivalent atoms
       ENDDO
!
         nat1 = nat1 + atoms%neq(n)
!
! end of different types of atoms
    ENDDO

!------------------------- FORCE PART -------------------------------
!+gu this is the remainder of spgset necessary for force calculations
!
      IF (input%l_f.OR.(atoms%n_u.GT.0)) THEN

      WRITE (6,FMT=&
     &  '(//,"list of operations which leave taual invariant",/)')
      DO na = 1,nat2
         WRITE (6,FMT='("atom nr.",i3,3x,(t14,"ops are:",24i3))') na,&
     &     (sym%invarop(na,nb),nb=1,sym%invarind(na))
      END DO

      ENDIF
!------------------------- FORCE PART ENDS --------------------------
!
!     check closure  ; note that:  {R|t} tau = R^{-1} tau -  R^{-1} t
!
!--->    loop over all operations
!
      WRITE (6,FMT=8040)
 8040 FORMAT (/,/,' multiplication table',/,/)
      sym%multab = 0
      DO j=1,sym%nop

!--->    multiply {R_j|t_j}{R_i|t_i}
         DO i=1,sym%nop
            mp = matmul( sym%mrot(:,:,j) , sym%mrot(:,:,i) )
            ttau = sym%tau(:,j) + matmul( sym%mrot(:,:,j) , sym%tau(:,i) )
            ttau = ttau - anint( ttau - 1.e-7 )

!--->    determine which operation this is
            DO k=1,sym%nop
              IF( all( mp(:,:) == sym%mrot(:,:,k) ) .AND.&
     &            ALL( abs( ttau(:)-sym%tau(:,k) ) < 1.e-7 ) ) THEN
                 IF (sym%multab(j,i) .EQ. 0 ) THEN
                    sym%multab(j,i) = k
                    IF (k .EQ. 1) sym%invtab(j)=i
                 ELSE
                    WRITE(6,'(" Symmetry error: multiple ops")')
                     CALL juDFT_error("Multiple ops",calledby="mapatom")
                 ENDIF
              ENDIF
            ENDDO

            IF (sym%multab(j,i).EQ.0) THEN
               WRITE (6,'(" Group not closed")')
               WRITE (6,'("  j , i =",2i4)') j,i
               CALL juDFT_error("mapatom: group not closed",calledby&
     &              ="mapatom")
            ENDIF
         ENDDO
      ENDDO

      DO n1 = 1,sym%nop
         WRITE (6,FMT=8060) (sym%multab(n1,n2),n2=1,sym%nop)
      END DO
 8060 FORMAT (1x,48i3)
      WRITE (6,FMT='(//," inverse operations",//)')
      DO n1 = 1,sym%nop
         WRITE (6,FMT=8060) n1,sym%invtab(n1)
      END DO

      DO na = 1,atoms%nat
         atoms%invsat(na) = 0
         sym%invsatnr(na) = 0
      END DO

      IF (sym%invs) THEN
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
                              cycle naloop
                           END IF
                        END DO
                      END DO
                    END DO
               END DO naloop
               END IF
            END DO
            nat1 = nat1 + atoms%neq(n)
         END DO
      WRITE (6,FMT=*) atoms%invsat
 9000 FORMAT ('atom type',i3,': atom',i3,' can be mapped into atom',i3,&
     &       ' via 3d inversion')
      END IF
 
      END  SUBROUTINE mapatom
      END  MODULE m_mapatom
