      MODULE m_rwnoco
      USE m_juDFT
!---------------------------------------------------------------------
!     read or write the nocoinp-file
!---------------------------------------------------------------------
      CONTAINS
        SUBROUTINE rw_noco_read(atoms,noco,input)
          USE m_types_atoms
          USE m_types_noco
          USE m_types_input
          USE m_constants
          IMPLICIT NONE
          ! ..
          ! ..  Arguments ..
          TYPE(t_atoms),INTENT(IN)   :: atoms
          TYPE(t_noco),INTENT(INOUT) :: noco
          TYPE(t_input),INTENT(INOUT):: input


          ! ..
          INTEGER :: itype, j, fileend
          REAL    :: qsc(3),rdum
          CHARACTER(len=8) :: inpchar
          logical          :: ldum



          DO itype = 1,atoms%ntype
             READ (24,'(22x,l1)') ldum

             inpchar(1:2)= 'XX'
             READ (24,fmt='(19x,a2)') inpchar(1:2)
             BACKSPACE (24)
             IF (inpchar(1:2)/='pi') THEN
                READ (24,fmt='(7x,f14.10,11x,f14.10)') &
                     &      noco%alph_inp(itype),rdum
             ELSE
                READ (24,fmt='(7x,f12.8,13x,f14.10)') &
                     &      noco%alph_inp(itype),rdum
                noco%alph_inp(itype)= noco%alph_inp(itype)*pi_const
             ENDIF
             inpchar(1:2)= 'XX'
             READ (24,fmt='(19x,a2)') inpchar(1:2)
             BACKSPACE (24)
             IF (inpchar(1:2)/='pi') THEN
                READ (24,fmt='(7x,f14.10,11x,f14.10)') &
                     &      noco%beta_inp(itype),rdum
             ELSE
                READ (24,fmt='(7x,f12.8,13x,f14.10)') &
                     &      noco%beta_inp(itype),rdum
                noco%beta_inp(itype)= noco%beta_inp(itype)*pi_const
             ENDIF

             READ (24,*)
                WRITE (oUnit,8010)        itype,.false.
                WRITE (oUnit,8020)        noco%alph_inp(itype)
                WRITE (oUnit,8025)         noco%beta_inp(itype)
                WRITE (oUnit,*)
          ENDDO
8010      FORMAT ('atom-type',i4,',l_relax=',l1,',l_magn=',l1,&
               &',M=',f8.5,',magtype=',i4)
8020      FORMAT ('alpha =',f14.10,',b_cons_x =',f14.10)
8025      FORMAT ('beta  =',f14.10,',b_cons_y =',f14.10)
8026      FORMAT ('Total number of magnetic atoms',i4,',magnetic types',i4)

          READ (24,*)
          READ (24,8036) noco%l_ss,noco%l_mperp,ldum
!!$          IF ( (inpchar(1:8)=='sso_opt=') .OR. (noco%l_ss .AND. noco%l_soc) ) THEN
!!$             BACKSPACE (24)
!!$             READ (24,fmt='(45x,2l1)') input%sso_opt(1),input%sso_opt(2)
!!$          ELSE
!!$             input%sso_opt(1)= .FALSE.
!!$             input%sso_opt(2)= .FALSE.
!!$          ENDIF
          READ (24,8046) noco%mix_b

          WRITE (oUnit,fmt='(5(A,l1),l1)') &
               & 'l_ss=',noco%l_ss,',l_mperp=',noco%l_mperp,',l_constr=',any(noco%l_constrained)
          WRITE (oUnit,8040) noco%mix_b
8030      FORMAT ('l_ss=',l1,',l_mperp=',l1,',l_constr=',l1,',l_disp=',l1)
8035      FORMAT (5x,l1,9x,l1,10x,l1,8x,l1)
8036      FORMAT (5x,l1,9x,l1,10x,l1)
8040      FORMAT ('mix_b=',f6.3,',thetaJ=',f14.10,',nsh=',i4)
8045      FORMAT (6x,f6.3,8x,f14.10,5x,i4)
8046      FORMAT (6x,f6.3)

          IF (noco%l_ss) THEN
             READ (24,fmt='(5x,3(f14.10,1x))') &
                  &    noco%qss_inp(1),noco%qss_inp(2),noco%qss_inp(3)
             inpchar(1:3)= 'XXX'
             READ (24,fmt='(a4)',iostat=fileend) inpchar(1:4)
             IF (fileend==0) THEN
                IF (inpchar(1:4)=='qsc=') THEN
                   BACKSPACE (24)
                   READ (24,fmt='(5x,3(f14.10,1x))') &
                        &      qsc(1),qsc(2),qsc(3)
                   DO j= 1,3
                      IF ( ABS(qsc(j)) < 1.e-6 ) THEN
                         WRITE (oUnit,fmt='(A,i1,A,1x,f14.10)')&
                              &          'Error reading nocoinp: qsc(',j,') =',qsc(j)
                         CALL juDFT_error("Error reading nocoinp",calledby&
                              &              ="rw_noco")
                      ENDIF
                      noco%qss_inp(j)= noco%qss_inp(j)/qsc(j)
                   ENDDO
                ENDIF
             ENDIF
             WRITE (oUnit,*) 'This is a Spin-Spiral (SS) calculation. The'
             WRITE (oUnit,*) 'q-vector of the Spin-Spiral is:'
             WRITE (oUnit,8060) noco%qss_inp
8060         FORMAT ('qss=(',f14.10,',',f14.10,',',f14.10,')')
          ENDIF
        END SUBROUTINE rw_noco_read

      SUBROUTINE rw_noco_write(atoms,noco,input)
      USE m_types_atoms
      USE m_types_noco
      USE m_types_input
      USE m_constants
      IMPLICIT NONE
! ..
! ..  Arguments ..
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_noco),INTENT(INOUT) :: noco
      TYPE(t_input),INTENT(IN)   :: input

! ..
      INTEGER :: itype
      CHARACTER(len=8) :: inpchar

      DO itype = 1,atoms%ntype
         WRITE (24,8010) itype,.false.
         WRITE (24,8020) noco%alph_inp(itype),0.0
         WRITE (24,8025)  noco%beta_inp(itype),0.0
         WRITE (24,*)
      ENDDO

      WRITE (24,*) '-- logical parameters --'
      WRITE (24,fmt='(5(A,l1),l1)') &
     & 'l_ss=',noco%l_ss,',l_mperp=',noco%l_mperp,',l_constr=',any(noco%l_constrained)
      WRITE (24,8040) noco%mix_b

      IF (noco%l_ss) THEN
       WRITE (24,8060) noco%qss_inp
       WRITE (oUnit,8060) noco%qss_inp
      ENDIF
8010      FORMAT ('atom-type',i4,',l_relax=',l1,',l_magn=',l1,&
               &',M=',f8.5,',magtype=',i4)
8020      FORMAT ('alpha =',f14.10,',b_cons_x =',f14.10)
8025      FORMAT ('beta  =',f14.10,',b_cons_y =',f14.10)
8040      FORMAT ('mix_b=',f6.3,',thetaJ=',f14.10,',nsh=',i4)
8060         FORMAT ('qss=(',f14.10,',',f14.10,',',f14.10,')')
      END SUBROUTINE rw_noco_write
      END MODULE m_rwnoco
