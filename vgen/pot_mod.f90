MODULE m_potmod
  USE m_juDFT
CONTAINS
  SUBROUTINE pot_mod(&
       &                   atoms,sphhar,vacuum,stars,&
       &                   input,&
       &                   vr,vxy,vz,vpw,vpw_w)


    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)  :: vacuum
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms
    REAL,    INTENT (INOUT) :: vr(atoms%jmtd,0:sphhar%nlhd,atoms%ntypd,input%jspins)
    REAL,    INTENT (INOUT) :: vz(vacuum%nmzd,2,input%jspins)
    COMPLEX, INTENT (INOUT) :: vpw(stars%n3d,input%jspins),vpw_w(stars%n3d,input%jspins)
    COMPLEX, INTENT (INOUT) :: vxy(vacuum%nmzxyd,stars%n2d-1,2,input%jspins)

    INTEGER i,j,n,ivac ,typmag,bxcflag,bxc_r,bxc_c,nat

    ! --- modify mag.-pot.  >> ---------------------------------------------

    ! what do you want ?
    bxcflag= 0
    typmag= atoms%ntype-1
    !  0 : potential is not changed
    !  1 : B_xc is kept in MTs and set to zero elsewhere
    !  2 : B_xc is kept in first typmag MTs and set to zero elsewhere
    !  3 : B_xc is read from file
    !  4 : B_xc is written to file

    IF (bxcflag/=0) THEN

       IF ( (bxcflag<0) .OR. (bxcflag>5) ) THEN
          CALL juDFT_error("bxcflag out of bounds",calledby="pot_mod")
       ENDIF
       IF ( (bxcflag==2) .AND. ((typmag<0).OR.(typmag>atoms%ntype)) ) THEN
          CALL juDFT_error("typmag out of bounds",calledby="pot_mod")
       ENDIF
       IF (input%jspins/=2) THEN
          CALL juDFT_error("no B-field as input%jspins==1",calledby="pot_mod")
       ENDIF

       IF (bxcflag/=4) THEN

          IF (bxcflag==3) THEN
             OPEN(201,file='bxc',form='unformatted',action='read')
             DO j= 1,atoms%ntype+2
                READ(201)
             ENDDO
          ENDIF
          IF (bxcflag/=1) THEN
             nat= 1
             DO n= 1,atoms%ntype
                IF ( (bxcflag==3) .OR. (n>typmag) ) THEN
                   DO j= 0,sphhar%nlh(atoms%ntypsy(nat))
                      DO i= 1,atoms%jri(n)
                         vr(i,j,n,1)= ( vr(i,j,n,1)+vr(i,j,n,2) )/2.
                         vr(i,j,n,2)= vr(i,j,n,1)
                         IF (bxcflag==3) THEN
                            READ(201) bxc_r
                            vr(i,j,n,1)= vr(i,j,n,1) + bxc_r
                            vr(i,j,n,2)= vr(i,j,n,2) - bxc_r
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
                nat= nat + atoms%neq(n)
             ENDDO
          ENDIF
          DO j= 1,stars%ng3
             vpw(j,1)= ( vpw(j,1)+vpw(j,2) )/2.
             vpw(j,2)= vpw(j,1)
             vpw_w(j,1)= ( vpw_w(j,1)+vpw_w(j,2) )/2.
             vpw_w(j,2)= vpw_w(j,1)
             IF (bxcflag==3) THEN
                READ(201) bxc_c
                vpw(j,1)= vpw(j,1) + bxc_c
                vpw(j,2)= vpw(j,2) - bxc_c
                READ(201) bxc_c
                vpw_w(j,1)= vpw_w(j,1) + bxc_c
                vpw_w(j,2)= vpw_w(j,2) - bxc_c
             ENDIF
          ENDDO
          IF (input%film) THEN
             DO ivac= 1,vacuum%nvac
                DO i= 1,vacuum%nmz
                   vz(i,ivac,1)= ( vz(i,ivac,1)+vz(i,ivac,2) )/2.
                   vz(i,ivac,2)= vz(i,ivac,1)
                   IF (bxcflag==3) THEN
                      READ(201) bxc_r
                      vz(i,ivac,1)= vz(i,ivac,1) + bxc_r
                      vz(i,ivac,2)= vz(i,ivac,2) - bxc_r
                   ENDIF
                ENDDO
                DO i= 1,vacuum%nmzxy
                   DO j= 1,stars%ng2-1
                      vxy(i,j,ivac,1)=&
                           &             ( vxy(i,j,ivac,1)+vxy(i,j,ivac,2) )/2.
                      vxy(i,j,ivac,2)= vxy(i,j,ivac,1)
                      IF (bxcflag==3) THEN
                         READ(201) bxc_c
                         vxy(i,j,ivac,1)= vxy(i,j,ivac,1) + bxc_c
                         vxy(i,j,ivac,2)= vxy(i,j,ivac,2) - bxc_c
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          WRITE(6,fmt='(1x)')
          SELECT CASE (bxcflag)
          CASE(1)
             WRITE(6,fmt='(A)') 'B_xc outside MTs is set to zero !!'
          CASE(2)
             WRITE(6,fmt='(A,i3,1x,A)')&
                  &         'B_xc outside the first',typmag,'MTs is set to zero !!'
          CASE(3)
             CLOSE(201)
             WRITE(6,fmt='(A)') 'B_xc is read from file "bxc" !!'
          END SELECT
          WRITE(6,fmt='(1x)')

       ELSE ! (bxc==4), write B_xc it file

          OPEN(201,file='bxc',form='unformatted',status='replace')
          WRITE(201) stars%ng3, atoms%ntype, input%film
          nat= 1
          DO n= 1,atoms%ntype
             WRITE(201) sphhar%nlh(atoms%ntypsy(nat)), atoms%jri(n), atoms%neq(n)
             nat= nat + atoms%neq(n)
          ENDDO
          IF (input%film) THEN
             WRITE(201) vacuum%nvac, vacuum%nmz, vacuum%nmzxy, stars%ng2
          ELSE
             WRITE(201) 0,1,1,1
          ENDIF
          nat= 1
          DO n= 1,atoms%ntype
             DO j= 0,sphhar%nlh(atoms%ntypsy(nat))
                DO i= 1,atoms%jri(n)
                   WRITE(201) ( vr(i,j,n,1)-vr(i,j,n,2) )/2.
                ENDDO
             ENDDO
             nat= nat + atoms%neq(n)
          ENDDO
          DO j= 1,stars%ng3
             WRITE(201) ( vpw(j,1)-vpw(j,2) )/2.
             WRITE(201) ( vpw_w(j,1)-vpw_w(j,2) )/2.
          ENDDO
          IF (input%film) THEN
             DO ivac= 1,vacuum%nvac
                DO i= 1,vacuum%nmz
                   WRITE(201) ( vz(i,ivac,1)-vz(i,ivac,2) )/2.
                ENDDO
                DO i= 1,vacuum%nmzxy
                   DO j= 1,stars%ng2-1
                      WRITE(201) ( vxy(i,j,ivac,1)-vxy(i,j,ivac,2) )/2.
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
          CLOSE(201)
          CALL juDFT_end("B_xc is written to 'bxc'")

       ENDIF

    ENDIF ! bxcflag/=0


  END SUBROUTINE pot_mod
END MODULE m_potmod
