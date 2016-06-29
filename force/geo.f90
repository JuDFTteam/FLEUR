!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_geo
  USE m_juDFT
CONTAINS
  SUBROUTINE geo(&
       &               atoms,sym,cell,oneD,input_in,tote,&
       &               forcetot)

    !    *********************************************************************
    !    * calculates the NEW atomic positions after the results%force calculation   *
    !    * SUBROUTINE is based on a BFGS method implemented by jij%M. Weinert    *
    !    *                                [cf. PRB 52 (9) p. 6313 (1995)]    * 
    !    *                                                                   *
    !    * as a first step we READ in the file 'inp' WITH some additional    *
    !    * information (SUBROUTINE rw_inp)                                   *
    !    * THEN recover the old geometry optimisation information from file  *
    !    * 'forces.dat' (SUBROUTINE bfsg0)                                   *
    !    * this input together WITH the NEW forces (forcetot) are now used   *
    !    * to calculate the NEW atomic positions (SUBROUTINE bfsg)           *
    !    * finally the NEW 'inp' file is written (SUBROUTINE rw_inp)         *
    !    *                                                           Gustav  *
    !
    ! input: 
    !        ntype .... total number of atom types
    !        thetad ... approx. debye temperature
    !        zat(ntype) mass number of the atom (or atomic number)
    !        xa ....... mixing factor 
    !        epsdisp .. limit for displacement to be converged
    !        epsforce . the same for force
    !        istepnow . steps to be done in this run
    !
    !    *********************************************************************
    USE m_rwinp
    USE m_bfgs
    USE m_bfgs0
    USE m_types
    IMPLICIT NONE
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_input),INTENT(IN):: input_in
    ! ..
    ! ..  Scalar Arguments ..
    REAL,    INTENT (IN) :: tote
    ! ..
    ! ..  Array Arguments ..
    REAL,    INTENT (INOUT) :: forcetot(3,atoms%ntypd)
    ! ..
    ! ..  Local Scalars ..
    INTEGER i,j,na ,istep0,istep,itype,jop,ieq
    LOGICAL lconv       
    TYPE(t_atoms)  :: atoms_new
    ! ..
    ! ..  Local Arrays ..
    REAL xold(3*atoms%ntypd),y(3*atoms%ntypd),h(3*atoms%ntypd,3*atoms%ntypd),zat(atoms%ntypd)
    REAL tau0(3,atoms%ntypd),tau0_i(3,atoms%ntypd) 

    TYPE(t_input):: input

    input=input_in

    CALL bfgs0(&
         &           atoms%ntype,&
         &           istep0,xold,y,h)

    DO itype=1,atoms%ntype
       IF (atoms%l_geo(itype)) THEN
          WRITE (6,'(6f10.5)') (tau0(j,itype),j=1,3),&
               &                         (forcetot(i,itype),i=1,3)
          DO i = 1,3
             forcetot(i,itype)=forcetot(i,itype)*REAL(atoms%relax(i,itype))
          ENDDO
          WRITE (6,'(6f10.5,a,3i2)') (tau0(j,itype),j=1,3),&
               &             (forcetot(i,itype),i=1,3),' atoms%relax: ',&
               &             (atoms%relax(i,itype),i=1,3)
       ELSE
          DO i = 1,3
             forcetot(i,itype)=0.0
          ENDDO
       ENDIF
    ENDDO

    istep = 1
    CALL bfgs(&
         &          atoms%ntype,istep,istep0,forcetot,&
         &          zat,input%xa,input%thetad,input%epsdisp,input%epsforce,tote,&
         &          xold,y,h,tau0,&
         &          lconv)
    IF (lconv) THEN
       WRITE (6,'(a)') "Des woars!"
       CALL juDFT_end(" GEO Des woars ", 1) ! The 1 is temporarily. Should be mpi%irank.
    ELSE

       atoms_new=atoms
       na = 0
       DO itype=1,atoms%ntype
          tau0_i(:,itype)=MATMUL(cell%bmat,tau0(:,itype))
          DO ieq = 1,atoms%neq(itype)
             na = na + 1
             jop = sym%invtab(atoms%ngopr(na))
             IF (oneD%odi%d1) jop = oneD%ods%ngopr(na)
             DO i = 1,3
                atoms_new%taual(i,na) = 0.0
                DO j = 1,3
                   IF (.NOT.oneD%odi%d1) THEN
                      atoms_new%taual(i,na) = atoms_new%taual(i,na) +&
                           &                   sym%mrot(i,j,jop) * tau0_i(j,itype)
                   ELSE
                      atoms_new%taual(i,na) = atoms_new%taual(i,na) +&
                           &                   oneD%ods%mrot(i,j,jop) * tau0_i(j,itype)
                   END IF
                ENDDO
                IF (oneD%odi%d1) THEN
                   atoms_new%taual(i,na) = atoms_new%taual(i,na) +&
                        &              oneD%ods%tau(i,jop)/cell%amat(3,3)
                ELSE
                   atoms_new%taual(i,na) = atoms_new%taual(i,na) + sym%tau(i,jop)
                END IF
             ENDDO
          ENDDO
       ENDDO

       CALL judft_error("Writing on new input file not implemented in geo")
       !        input%l_f = .false.
       !        CALL rw_inp(&
       !     &            'W',atoms_new,obsolete,vacuum,input,stars,sliceplot,banddos,&
       !     &            cell,sym,xcpot,noco,jij,oneD,hybrid,kpts,&
       !     &            noel,namex,relcor,a1,a2,a3,scale,dtild,name)

    ENDIF

    RETURN
  END SUBROUTINE geo
END MODULE m_geo
