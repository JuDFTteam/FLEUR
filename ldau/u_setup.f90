!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_usetup
  USE m_juDFT
  !-------------------------------------------------------------------+
  ! Sets up the quantities needed for the LDA+U subroutines:          |
  !     radial integrals: us,dus,uds,duds                             |
  !     overlap of dot u: ddn                                         |
  !     potential matrix: vs_mmp                                      |
  !     total energy contribution: e_ldau                             |
  !                                                  G.B. Oct. 2000   |
  !                                                                   |
  !     Extension to multiple U per atom type  G.M. 2017              |
  !-------------------------------------------------------------------+
CONTAINS
  SUBROUTINE u_setup(sym,atoms,lmaxb,sphhar, input,el,vr,mpi, vs_mmp,results,number)
    USE m_umtx
    USE m_uj2f
    USE m_nmat_rot
    USE m_vmmp
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_results),INTENT(INOUT)   :: results
    TYPE(t_mpi),INTENT(IN)          :: mpi
    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_sphhar),INTENT(IN)       :: sphhar
    TYPE(t_atoms),INTENT(INOUT)     :: atoms !n_u might be modified if no density matrix is found

    ! ... Arguments ...
    INTEGER, INTENT (IN) :: lmaxb  

    REAL,    INTENT (IN) :: el(0:,:,:) !(0:atoms%lmaxd,ntype,jspd)
    REAL,    INTENT (IN) :: vr(:,0:,:,:) !(atoms%jmtd,0:sphhar%nlhd,ntype,jspd)
    COMPLEX,INTENT (OUT)::vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins)
    INTEGER,INTENT(IN),OPTIONAL::number
    ! ... Local Variables ...
    INTEGER itype,ispin,j,k,l,jspin,urec,i_u
    INTEGER noded,nodeu,ios,lty(atoms%n_u)
    REAL wronk
    LOGICAL n_mmp_exist,n_exist
    CHARACTER*8 l_type*2,l_form*9
    CHARACTER*12 ::filename
    REAL f(atoms%jmtd,2),g(atoms%jmtd,2),theta(atoms%n_u),phi(atoms%n_u),zero(atoms%n_u)
    REAL f0(atoms%n_u,input%jspins),f2(atoms%n_u,input%jspins),f4(atoms%n_u,input%jspins),f6(atoms%n_u,input%jspins)
    REAL, ALLOCATABLE :: u(:,:,:,:,:,:)
    COMPLEX, ALLOCATABLE :: ns_mmp(:,:,:,:)
    !
    ! look, whether density matrix exists already:
    !
    filename="n_mmp_mat"
    IF (PRESENT(number)) THEN
       WRITE(filename,"('n_mmp_mat.',i0)") number
    ENDIF
    INQUIRE (file=filename,exist=n_mmp_exist)
    IF (n_mmp_exist.AND.atoms%n_u>0) THEN
       !
       ! calculate slater integrals from u and j
       !
       CALL uj2f(input%jspins,atoms,&
                 f0,f2,f4,f6)
       !
       ! set up e-e- interaction matrix
       !
       ALLOCATE ( u(-lmaxb:lmaxb,-lmaxb:lmaxb,&
            -lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins) )
       DO ispin = 1, 1 ! input%jspins
          f0(:,1) = (f0(:,1) + f0(:,input%jspins) ) / 2
          f2(:,1) = (f2(:,1) + f2(:,input%jspins) ) / 2
          f4(:,1) = (f4(:,1) + f4(:,input%jspins) ) / 2
          f6(:,1) = (f6(:,1) + f6(:,input%jspins) ) / 2
          CALL umtx(atoms,lmaxb,f0(1,ispin),&
                    f2(1,ispin),f4(1,ispin),f6(1,ispin),&
                    u(-lmaxb,-lmaxb,-lmaxb,-lmaxb,1,ispin) )
       END DO
       !
       ! read density matrix
       !
       ALLOCATE (ns_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins))
       OPEN (69,file=filename,status='unknown',form='formatted')
       READ (69,9000) ns_mmp
9000   FORMAT(7f20.13)
       CLOSE (69)
       !
       ! check for possible rotation of n_mmp
       !
       INQUIRE (file='n_mmp_rot',exist=n_exist)
       IF (n_exist) THEN
          OPEN (68,file='n_mmp_rot',status='old',form='formatted')
          DO i_u = 1, atoms%n_u
             itype = atoms%lda_u(i_u)%atomType
             l = atoms%lda_u(i_u)%l
             READ(68,*,iostat=ios) theta(i_u),phi(i_u)
             IF (ios == 0) THEN
                lty(i_u) = l
             ELSE
                IF (i_u == 1)  CALL juDFT_error("ERROR reading n_mmp_rot",calledby ="u_setup")
                theta(i_u) = theta(i_u-1) ; phi(i_u) = phi(i_u-1)
                lty(i_u) = lty(i_u-1)
             END IF
          END DO
          CLOSE (68)
          zero = 0.0
          CALL nmat_rot(phi,theta,zero,3,atoms%n_u,input%jspins,lty,ns_mmp)
       ENDIF
       !
       ! calculate potential matrix and total energy correction
       !
       CALL v_mmp(sym,atoms,input%jspins,lmaxb,ns_mmp,u,f0,f2,&
                  vs_mmp,results)
       IF (mpi%irank.EQ.0) THEN
          DO jspin = 1,input%jspins
             WRITE (6,'(a7,i3)') 'spin #',jspin
             DO i_u = 1, atoms%n_u
                itype = atoms%lda_u(i_u)%atomType
                l = atoms%lda_u(i_u)%l
                WRITE (l_type,'(i2)') 2*(2*l+1)
                l_form = '('//l_type//'f12.7)'
                WRITE (6,'(a20,i3)') 'n-matrix for atom # ',itype
                WRITE (6,l_form) ((ns_mmp(k,j,i_u,jspin),k=-l,l),j=-l,l)
                WRITE (6,'(a20,i3)') 'V-matrix for atom # ',itype
                IF (atoms%lda_u(i_u)%l_amf) THEN
                   WRITE (6,*) 'using the around-mean-field limit '
                ELSE
                   WRITE (6,*) 'using the atomic limit of LDA+U '
                ENDIF
                WRITE (6,l_form) ((vs_mmp(k,j,i_u,jspin),k=-l,l),j=-l,l)
             END DO
          END DO
          WRITE (6,*) results%e_ldau
       ENDIF
       DEALLOCATE ( u,ns_mmp )
    ELSE
       IF (mpi%irank.EQ.0) THEN
          WRITE (*,*) 'no density matrix found ... skipping LDA+U'
          WRITE(*,*) "File:",filename
       ENDIF
       vs_mmp(:,:,:,:) = CMPLX(0.0,0.0)
       results%e_ldau = 0.0
       atoms%n_u = 0
    ENDIF

    RETURN
  END SUBROUTINE u_setup
END MODULE m_usetup
