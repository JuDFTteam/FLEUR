    !
    ! plot matrix
    !
    IF (l_plot) THEN
       IF (mpi%isize /= 1) STOP 'coulombmatrix: l_plot only works with one process'
       DO ikpt = 1, kpts%nkpt

          ic = 0
          DO i = 1, nbasm1(ikpt)
             DO j = 1, i
                ic = ic + 1
                IF (ABS(coulomb(ic, ikpt)) > 1E-8) THEN
                   WRITE (600 + ikpt, '(2i6)') i, j
                   WRITE (600 + ikpt, '(2i6)') j, i
                END IF
             END DO
          END DO

          ALLOCATE (coulhlp(nbasm1(ikpt), nbasm1(ikpt)))
          coulhlp = unpackmat(coulomb(:nbasm1(ikpt)*(nbasm1(ikpt) + 1)/2, ikpt))

          ic = 0
          DO itype = 1, atoms%ntype
             DO ineq = 1, atoms%neq(itype)
                DO l = 0, hybrid%lcutm1(itype)
                   DO M = -l, l
                      WRITE (700 + ikpt, *) l, M
                      DO n = 1, hybrid%nindxm1(l, itype)
                         WRITE (700 + ikpt, '(16f8.4)') coulhlp(ic + n, ic + 1:ic + hybrid%nindxm1(l, itype))
                      END DO
                      ic = ic + hybrid%nindxm1(l, itype)
                   ENDDO
                END DO
             END DO
          END DO
          ALLOCATE (coulhlp1(nbasm1(ikpt), nbasm1(ikpt)))
          coulhlp1 = 0

          ic2 = 0
          DO itype = 1, atoms%ntype
             DO ineq = 1, atoms%neq(itype)
                DO l = 0, hybrid%lcutm1(itype)
                   DO M = -l, l
                      DO n = 1, hybrid%nindxm1(l, itype) - 1
                         ic2 = ic2 + 1
                      END DO
                   END DO
                END DO
             END DO
          END DO

          ic = 0
          ic1 = 0
          DO itype = 1, atoms%ntype
             DO ineq = 1, atoms%neq(itype)
                DO l = 0, hybrid%lcutm1(itype)
                   DO M = -l, l
                      DO n = 1, hybrid%nindxm1(l, itype) - 1
                         coulhlp1(ic + n, ic + 1:ic + hybrid%nindxm1(l, itype) - 1) &
                            = coulhlp(ic1 + n, ic1 + 1:ic1 + hybrid%nindxm1(l, itype) - 1)
                      END DO

                      coulhlp1(ic2 + 1, ic + 1:ic + hybrid%nindxm1(l, itype) - 1) &
                         = coulhlp(ic1 + hybrid%nindxm1(l, itype), ic1 + 1:ic1 + hybrid%nindxm1(l, itype) - 1)

                      ic = ic + hybrid%nindxm1(l, itype) - 1
                      ic1 = ic1 + hybrid%nindxm1(l, itype)
                      ic2 = ic2 + 1

                   END DO
                END DO
             END DO
          END DO

          IF (ikpt == 1) THEN
             ! add contributions from s channel and G=0 component

             ic = 0; ic1 = 0
             DO itype = 1, atoms%ntype
                DO ineq = 1, atoms%neq(itype)
                   WRITE (*, *) ic + 1, ic + hybrid%nindxm1(0, itype) - 1, ic1 + 1, &
                      ic1 + hybrid%nindxm1(0, itype) - 1
                   coulhlp1(ic + 1:ic + hybrid%nindxm1(0, itype) - 1, nbasp + 1) &
                      = coulhlp(ic1 + 1:ic1 + hybrid%nindxm1(0, itype) - 1, nbasp + 1)
                   coulhlp1(nbasp + 1, ic + 1:ic + hybrid%nindxm1(0, itype) - 1) &
                      = coulhlp(nbasp + 1, ic1 + 1:ic1 + hybrid%nindxm1(0, itype) - 1)
                   ic = ic + SUM((/((2*l + 1)*(hybrid%nindxm1(l, itype) - 1), &
                                    l=0, hybrid%lcutm1(itype))/))
                   ic1 = ic1 + SUM((/((2*l + 1)*hybrid%nindxm1(l, itype), &
                                      l=0, hybrid%lcutm1(itype))/))
                END DO
             END DO

             ic2 = 0
             DO itype = 1, atoms%ntype
                DO ineq = 1, atoms%neq(itype)
                   DO l = 0, hybrid%lcutm1(itype)
                      DO M = -l, l
                         DO n = 1, hybrid%nindxm1(l, itype) - 1
                            ic2 = ic2 + 1
                         END DO
                      END DO
                   END DO
                END DO
             END DO

             ic1 = 0
             DO itype = 1, atoms%ntype
                DO ineq = 1, atoms%neq(itype)
                   DO l = 0, hybrid%lcutm1(itype)
                      DO M = -l, l
                         ic2 = ic2 + 1

                         ic1 = ic1 + hybrid%nindxm1(l, itype)

                         IF (l /= 0) CYCLE
                         WRITE (900, *) ic2, ic1, itype, ineq

                         ic3 = 0
                         ic4 = 0
                         DO itype1 = 1, atoms%ntype
                            ishift = SUM((/((2*l2 + 1)*hybrid%nindxm1(l2, itype1), &
                                            l2=0, hybrid%lcutm1(itype1))/))
                            ishift1 = SUM((/((2*l2 + 1)*(hybrid%nindxm1(l2, itype1) - 1), &
                                             l2=0, hybrid%lcutm1(itype1))/))
                            DO ineq1 = 1, atoms%neq(itype1)
                               ic5 = ic3 + (ineq1 - 1)*ishift + 1
                               ic6 = ic5 + hybrid%nindxm1(0, itype1) - 2

                               ic7 = ic4 + (ineq1 - 1)*ishift1 + 1
                               ic8 = ic7 + hybrid%nindxm1(0, itype1) - 2
                               WRITE (901, *) ic2, ic7, ic8, ic1, ic5, ic6, itype, itype1

                               coulhlp1(ic2, ic7:ic8) = coulhlp(ic1, ic5:ic6)
                               IF (sym%invs) THEN
                                  coulhlp1(ic7:ic8, ic2) = coulhlp1(ic2, ic7:ic8)
                               ELSE
                                  coulhlp1(ic7:ic8, ic2) = CONJG(coulhlp1(ic2, ic7:ic8))
                               ENDIF
                            END DO
                            ic3 = ic3 + ishift*atoms%neq(itype1)
                            ic4 = ic4 + ishift1*atoms%neq(itype1)
                         END DO

                      END DO
                   END DO
                END DO
             END DO

          END IF

          ic2 = 0
          DO itype = 1, atoms%ntype
             DO ineq = 1, atoms%neq(itype)
                DO l = 0, hybrid%lcutm1(itype)
                   DO M = -l, l
                      DO n = 1, hybrid%nindxm1(l, itype) - 1
                         ic2 = ic2 + 1
                      END DO
                   END DO
                END DO
             END DO
          END DO

          ic1 = 0
          DO itype = 1, atoms%ntype
             DO ineq = 1, atoms%neq(itype)
                DO l = 0, hybrid%lcutm1(itype)
                   DO M = -l, l
                      ic2 = ic2 + 1

                      ic1 = ic1 + hybrid%nindxm1(l, itype)

                      ic3 = 0
                      ic4 = ic2
                      DO itype1 = 1, atoms%ntype
                         DO ineq1 = 1, atoms%neq(itype1)
                            DO l1 = 0, hybrid%lcutm1(itype1)
                               DO m1 = -l1, l1
                                  ic3 = ic3 + hybrid%nindxm1(l1, itype1)
                                  IF (ic3 < ic1) CYCLE
                                  WRITE (300, '(4i6,2f15.10)') ic2, ic4, ic1, ic3, coulhlp(ic1, ic3)
                                  coulhlp1(ic2, ic4) = coulhlp(ic1, ic3)
                                  IF (sym%invs) THEN
                                     coulhlp1(ic4, ic2) = coulhlp1(ic2, ic4)
                                  ELSE
                                     coulhlp1(ic4, ic2) = CONJG(coulhlp1(ic2, ic4))
                                  ENDIF
                                  ic4 = ic4 + 1
                               END DO
                            END DO
                         END DO
                      END DO

                      DO igpt = 1, hybrid%ngptm(ikpt)
                         coulhlp1(ic2, ic4) = coulhlp(ic1, nbasp + igpt)
                         IF (sym%invs) THEN
                            coulhlp1(ic4, ic2) = coulhlp1(ic2, ic4)
                         ELSE
                            coulhlp1(ic4, ic2) = CONJG(coulhlp1(ic2, ic4))
                         ENDIF
                         ic4 = ic4 + 1
                      END DO

                   END DO
                END DO
             END DO
          END DO

          coulhlp1(nbasp + 1:, nbasp + 1:) = coulhlp(nbasp + 1:, nbasp + 1:)

          ic = 0
          DO itype = 1, atoms%ntype
             DO ineq = 1, atoms%neq(itype)
                DO l = 0, hybrid%lcutm1(itype)

                   DO M = -l, l
                      WRITE (800 + ikpt, *) l, M
                      DO n = 1, hybrid%nindxm1(l, itype) - 1
                         WRITE (800 + ikpt, '(16f8.4)') coulhlp1(ic + n, ic + 1:ic + hybrid%nindxm1(l, itype) - 1)
                      END DO
                      ic = ic + hybrid%nindxm1(l, itype) - 1
                   ENDDO
                END DO
             END DO
          END DO

          ic = 0
          DO i = 1, nbasm1(ikpt)
             DO j = 1, i
                IF (ABS(coulhlp1(i, j)) > 1E-8) THEN
                   WRITE (850 + ikpt, '(2i6)') i, j
                   WRITE (850 + ikpt, '(2i6)') j, i
                END IF
             END DO
          END DO
          DEALLOCATE (coulhlp, coulhlp1)
       END DO
       STOP
    END IF ! lplot
