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
          REAL    :: qsc(3)
          CHARACTER(len=8) :: inpchar 



          DO itype = 1,atoms%ntype
             READ (24,'(22x,l1)') noco%l_relax(itype)
             
             inpchar(1:2)= 'XX'
             READ (24,fmt='(19x,a2)') inpchar(1:2)
             BACKSPACE (24)
             IF (inpchar(1:2)/='pi') THEN
                READ (24,fmt='(7x,f14.10,11x,f14.10)') &
                     &      noco%alph(itype),noco%b_con(1,itype)
             ELSE
                READ (24,fmt='(7x,f12.8,13x,f14.10)') &
                     &      noco%alph(itype),noco%b_con(1,itype)
                noco%alph(itype)= noco%alph(itype)*pi_const
             ENDIF
             inpchar(1:2)= 'XX'
             READ (24,fmt='(19x,a2)') inpchar(1:2)
             BACKSPACE (24)
             IF (inpchar(1:2)/='pi') THEN
                READ (24,fmt='(7x,f14.10,11x,f14.10)') &
                     &      noco%beta(itype),noco%b_con(2,itype)
             ELSE
                READ (24,fmt='(7x,f12.8,13x,f14.10)') &
                     &      noco%beta(itype),noco%b_con(2,itype)
                noco%beta(itype)= noco%beta(itype)*pi_const
             ENDIF

             READ (24,*)
                WRITE (6,8010)        itype,noco%l_relax(itype)
                WRITE (6,8020)        noco%alph(itype),noco%b_con(1,itype)
                WRITE (6,8025)         noco%beta(itype),noco%b_con(2,itype)
                WRITE (6,*)
          ENDDO
8010      FORMAT ('atom-type',i4,',l_relax=',l1,',l_magn=',l1,&
               &',M=',f8.5,',magtype=',i4)
8020      FORMAT ('alpha =',f14.10,',b_cons_x =',f14.10)
8025      FORMAT ('beta  =',f14.10,',b_cons_y =',f14.10)
8026      FORMAT ('Total number of magnetic atoms',i4,',magnetic types',i4)

          READ (24,*)
          READ (24,8036) noco%l_ss,noco%l_mperp,noco%l_constr
!!$          IF ( (inpchar(1:8)=='sso_opt=') .OR. (noco%l_ss .AND. noco%l_soc) ) THEN
!!$             BACKSPACE (24)
!!$             READ (24,fmt='(45x,2l1)') input%sso_opt(1),input%sso_opt(2)
!!$          ELSE
!!$             input%sso_opt(1)= .FALSE. 
!!$             input%sso_opt(2)= .FALSE.
!!$          ENDIF
          READ (24,8046) noco%mix_b

          WRITE (6,fmt='(5(A,l1),l1)') &
               & 'l_ss=',noco%l_ss,',l_mperp=',noco%l_mperp,',l_constr=',noco%l_constr
          WRITE (6,8040) noco%mix_b
8030      FORMAT ('l_ss=',l1,',l_mperp=',l1,',l_constr=',l1,',l_disp=',l1)
8035      FORMAT (5x,l1,9x,l1,10x,l1,8x,l1)
8036      FORMAT (5x,l1,9x,l1,10x,l1)
8040      FORMAT ('mix_b=',f6.3,',thetaJ=',f14.10,',nsh=',i4)
8045      FORMAT (6x,f6.3,8x,f14.10,5x,i4)
8046      FORMAT (6x,f6.3)

          IF (noco%l_ss) THEN
             READ (24,fmt='(5x,3(f14.10,1x))') &
                  &    noco%qss(1),noco%qss(2),noco%qss(3)
             inpchar(1:3)= 'XXX'
             READ (24,fmt='(a4)',iostat=fileend) inpchar(1:4)
             IF (fileend==0) THEN
                IF (inpchar(1:4)=='qsc=') THEN
                   BACKSPACE (24)
                   READ (24,fmt='(5x,3(f14.10,1x))') &
                        &      qsc(1),qsc(2),qsc(3)
                   DO j= 1,3
                      IF ( ABS(qsc(j)) < 1.e-6 ) THEN
                         WRITE (6,fmt='(A,i1,A,1x,f14.10)')&
                              &          'Error reading nocoinp: qsc(',j,') =',qsc(j)
                         CALL juDFT_error("Error reading nocoinp",calledby&
                              &              ="rw_noco")
                      ENDIF
                      noco%qss(j)= noco%qss(j)/qsc(j)
                   ENDDO
                ENDIF
             ENDIF
             WRITE (6,*) 'This is a Spin-Spiral (SS) calculation. The'
             WRITE (6,*) 'q-vector of the Spin-Spiral is:'
             WRITE (6,8060) noco%qss(1),noco%qss(2),noco%qss(3)
8060         FORMAT ('qss=(',f14.10,',',f14.10,',',f14.10,')')
          ENDIF
        END SUBROUTINE rw_noco_read

      SUBROUTINE rw_noco_write(atoms,noco,input)
      USE m_types_atoms
      USE m_types_noco
      USE m_types_input
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
         IF ( noco%l_relax(itype) ) noco%b_con(:,itype) = 0.0
         WRITE (24,8010) itype,noco%l_relax(itype)
         WRITE (24,8020) noco%alph(itype),noco%b_con(1,itype)
         WRITE (24,8025)  noco%beta(itype),noco%b_con(2,itype)
         WRITE (24,*)
      ENDDO
 
      WRITE (24,*) '-- logical parameters --'
      WRITE (24,fmt='(5(A,l1),l1)') &
     & 'l_ss=',noco%l_ss,',l_mperp=',noco%l_mperp,',l_constr=',noco%l_constr
      WRITE (24,8040) noco%mix_b

      IF (noco%l_ss) THEN
       WRITE (24,8060) noco%qss(1),noco%qss(2),noco%qss(3)
       WRITE (6,8060) noco%qss(1),noco%qss(2),noco%qss(3)
      ENDIF
8010      FORMAT ('atom-type',i4,',l_relax=',l1,',l_magn=',l1,&
               &',M=',f8.5,',magtype=',i4)
8020      FORMAT ('alpha =',f14.10,',b_cons_x =',f14.10)
8025      FORMAT ('beta  =',f14.10,',b_cons_y =',f14.10)
8040      FORMAT ('mix_b=',f6.3,',thetaJ=',f14.10,',nsh=',i4)
8060         FORMAT ('qss=(',f14.10,',',f14.10,',',f14.10,')')
      END SUBROUTINE rw_noco_write
      END MODULE m_rwnoco
