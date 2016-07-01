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
    USE m_rinpXML
    USE m_winpXML
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

    ! temporary variables for XML IO
    TYPE(t_input)                 :: input_temp
    TYPE(t_dimension)             :: dimension_temp
    TYPE(t_atoms)                 :: atoms_temp
    TYPE(t_sphhar)                :: sphhar_temp
    TYPE(t_cell)                  :: cell_temp
    TYPE(t_stars)                 :: stars_temp
    TYPE(t_sym)                   :: sym_temp
    TYPE(t_noco)                  :: noco_temp
    TYPE(t_vacuum)                :: vacuum_temp
    TYPE(t_sliceplot)             :: sliceplot_temp
    TYPE(t_banddos)               :: banddos_temp
    TYPE(t_obsolete)              :: obsolete_temp
    TYPE(t_enpara)                :: enpara_temp
    TYPE(t_xcpot)                 :: xcpot_temp
    TYPE(t_results)               :: results_temp
    TYPE(t_jij)                   :: jij_temp
    TYPE(t_kpts)                  :: kpts_temp
    TYPE(t_hybrid)                :: hybrid_temp
    TYPE(t_oneD)                  :: oneD_temp
    LOGICAL                       :: l_opti_temp
    INTEGER                       :: numSpecies
    INTEGER                       :: div(3)
    INTEGER, ALLOCATABLE          :: xmlElectronStates(:,:)
    INTEGER, ALLOCATABLE          :: atomTypeSpecies(:)
    INTEGER, ALLOCATABLE          :: speciesRepAtomType(:)
    REAL, ALLOCATABLE             :: xmlCoreOccs(:,:,:)
    LOGICAL, ALLOCATABLE          :: xmlPrintCoreStates(:,:)
    CHARACTER(len=3), ALLOCATABLE :: noel_temp(:)
    CHARACTER(len=4)              :: namex_temp
    CHARACTER(len=12)             :: relcor_temp
    CHARACTER(LEN=20)             :: filename
    REAL                          :: a1_temp(3),a2_temp(3),a3_temp(3)
    REAL                          :: scale_temp, dtild_temp

    input=input_in
    atoms_new=atoms

    na = 1
    DO i = 1,atoms_new%ntype
       IF (input%film) atoms_new%taual(3,na) = atoms_new%taual(3,na)/cell%amat(3,3)
       tau0_i(:,i) = atoms_new%taual(:,na)
       tau0(:,i)=MATMUL(cell%amat,tau0_i(:,i))
       na = na + atoms_new%neq(i)
    END DO

    CALL bfgs0(atoms%ntype, istep0,xold,y,h)

    DO itype=1,atoms%ntype
       IF (atoms%l_geo(itype)) THEN
          WRITE (6,'(6f10.5)') (tau0(j,itype),j=1,3), (forcetot(i,itype),i=1,3)
          DO i = 1,3
             forcetot(i,itype)=forcetot(i,itype)*REAL(atoms%relax(i,itype))
          ENDDO
          WRITE (6,'(6f10.5,a,3i2)') (tau0(j,itype),j=1,3),&
               &     (forcetot(i,itype),i=1,3),' atoms%relax: ', (atoms%relax(i,itype),i=1,3)
       ELSE
          DO i = 1,3
             forcetot(i,itype)=0.0
          ENDDO
       ENDIF
    ENDDO

    istep = 1
    CALL bfgs(atoms%ntype,istep,istep0,forcetot,&
         &          zat,input%xa,input%thetad,input%epsdisp,input%epsforce,tote,&
         &          xold,y,h,tau0, lconv)
    IF (lconv) THEN
       WRITE (6,'(a)') "Des woars!"
       CALL juDFT_end(" GEO Des woars ", 1) ! The 1 is temporarily. Should be mpi%irank.
    ELSE
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
                      atoms_new%taual(i,na) = atoms_new%taual(i,na) + sym%mrot(i,j,jop) * tau0_i(j,itype)
                   ELSE
                      atoms_new%taual(i,na) = atoms_new%taual(i,na) + oneD%ods%mrot(i,j,jop) * tau0_i(j,itype)
                   END IF
                ENDDO
                IF (oneD%odi%d1) THEN
                   atoms_new%taual(i,na) = atoms_new%taual(i,na) + oneD%ods%tau(i,jop)/cell%amat(3,3)
                ELSE
                   atoms_new%taual(i,na) = atoms_new%taual(i,na) + sym%tau(i,jop)
                END IF
             ENDDO
          ENDDO
       ENDDO

!       CALL judft_error("Writing on new input file not implemented in geo")
       input%l_f = .false.
!       CALL rw_inp('W',atoms_new,obsolete,vacuum,input,stars,sliceplot,banddos,&
!                   cell,sym,xcpot,noco,jij,oneD,hybrid,kpts,&
!                   noel,namex,relcor,a1,a2,a3,scale,dtild,name)
       IF(input%l_inpXML) THEN
          ALLOCATE(noel_temp(1),atomTypeSpecies(1),speciesRepAtomType(1))
          ALLOCATE(xmlElectronStates(1,1),xmlPrintCoreStates(1,1))
          ALLOCATE(xmlCoreOccs(1,1,1))
          CALL r_inpXML(&
                        atoms_temp,obsolete_temp,vacuum_temp,input_temp,stars_temp,sliceplot_temp,&
                        banddos_temp,dimension_temp,cell_temp,sym_temp,xcpot_temp,noco_temp,Jij_temp,&
                        oneD_temp,hybrid_temp,kpts_temp,enpara_temp,sphhar_temp,l_opti_temp,noel_temp,&
                        namex_temp,relcor_temp,a1_temp,a2_temp,a3_temp,scale_temp,dtild_temp,xmlElectronStates,&
                        xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType)
          numSpecies = SIZE(speciesRepAtomType)
          filename = 'inp_new.xml'
          input_temp%l_f = input%l_f
          div(:) = MIN(kpts_temp%nmop(:),1)
          CALL w_inpXML(&
                        atoms_new,obsolete_temp,vacuum_temp,input_temp,stars_temp,sliceplot_temp,&
                        banddos_temp,cell_temp,sym_temp,xcpot_temp,noco_temp,jij_temp,oneD_temp,hybrid_temp,&
                        kpts_temp,kpts_temp%nmop,kpts_temp%l_gamma,noel_temp,namex_temp,relcor_temp,a1_temp,a2_temp,a3_temp,&
                        scale_temp,dtild_temp,input_temp%comment,xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
                        atomTypeSpecies,speciesRepAtomType,.FALSE.,filename,numSpecies,enpara_temp)
          DEALLOCATE(noel_temp,atomTypeSpecies,speciesRepAtomType)
          DEALLOCATE(xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs)
       END IF
    ENDIF

    RETURN
  END SUBROUTINE geo
END MODULE m_geo
