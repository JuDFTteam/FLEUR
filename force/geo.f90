!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_geo
  USE m_juDFT
CONTAINS
  SUBROUTINE geo(atoms,sym,cell,oneD,vacuum,input_in,tote,&
                 forcetot)

    !    *********************************************************************
    !    * calculates the NEW atomic positions after the results%force calculation   *
    !    * SUBROUTINE is based on a BFGS method implemented by M. Weinert    *
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
    USE m_types
    USE m_rwinp
    USE m_bfgs
    USE m_bfgs0
    USE m_constants
    USE m_rinpXML
    USE m_winpXML
    USE m_init_wannier_defaults
    USE m_xsf_io
    IMPLICIT NONE
    TYPE(t_oneD),INTENT(IN)   :: oneD
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_vacuum),INTENT(IN) :: vacuum
    TYPE(t_input),INTENT(IN)  :: input_in
    ! ..
    ! ..  Scalar Arguments ..
    REAL,    INTENT (IN) :: tote
    ! ..
    ! ..  Array Arguments ..
    REAL,    INTENT (INOUT) :: forcetot(3,atoms%ntype)
    ! ..
    ! ..  Local Scalars ..
    INTEGER i,j,na ,istep0,istep,itype,jop,ieq
    LOGICAL lconv
    TYPE(t_atoms)  :: atoms_new
    ! ..
    ! ..  Local Arrays ..
    REAL xold(3*atoms%ntype),y(3*atoms%ntype),h(3*atoms%ntype,3*atoms%ntype),zat(atoms%ntype)
    REAL tau0(3,atoms%ntype),tau0_i(3,atoms%ntype) 

    TYPE(t_input):: input

    ! temporary variables for XML IO
    TYPE(t_input)                 :: input_temp
    TYPE(t_dimension)             :: dimension_temp
    TYPE(t_atoms)                 :: atoms_temp
    TYPE(t_cell)                  :: cell_temp
    TYPE(t_stars)                 :: stars_temp
    TYPE(t_sym)                   :: sym_temp
    TYPE(t_noco)                  :: noco_temp
    TYPE(t_vacuum)                :: vacuum_temp
    TYPE(t_sliceplot)             :: sliceplot_temp
    TYPE(t_banddos)               :: banddos_temp
    TYPE(t_enpara)                :: enpara_temp
    CLASS(t_xcpot),ALLOCATABLE    :: xcpot_temp
    TYPE(t_results)               :: results_temp
    TYPE(t_kpts)                  :: kpts_temp
    TYPE(t_hybrid)                :: hybrid_temp
    TYPE(t_oneD)                  :: oneD_temp
    TYPE(t_coreSpecInput)         :: coreSpecInput_temp
    TYPE(t_wann)                  :: wann_temp
    LOGICAL                       :: l_kpts_temp
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
    REAL                          :: forceAllAtoms(3,atoms%nat)
    CLASS(t_forcetheo),ALLOCATABLE:: forcetheo
    TYPE(t_field)                 :: field
    input=input_in
    atoms_new=atoms

    istep0 = 0
    xold = 0.0
    y = 0.0
    h= 0.0

    na = 1
    DO i = 1,atoms_new%ntype
       zat(i)=real(atoms%nz(i))
       tau0(:,i)=atoms%pos(:,na)
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
    CALL bfgs(atoms%ntype,istep,istep0,forcetot,zat,input%xa,input%thetad,input%epsdisp,&
              input%epsforce,tote,xold,y,h,tau0,lconv)

    !write out struct_force.xsf file
    forceAllAtoms = 0.0
    na = 0
    DO itype=1,atoms%ntype
       forcetot(:,itype)=MATMUL(cell%bmat,forcetot(:,itype))/tpi_const ! to inner coordinates
       DO ieq = 1,atoms%neq(itype)
          na = na + 1
          jop = sym%invtab(sym%ngopr(na))
          IF (oneD%odi%d1) jop = oneD%ods%ngopr(na)
          DO i = 1,3
             DO j = 1,3
                IF (.NOT.oneD%odi%d1) THEN
                   forceAllAtoms(i,na) = forceAllAtoms(i,na) + sym%mrot(i,j,jop) * forcetot(j,itype)
                ELSE
                   forceAllAtoms(i,na) = forceAllAtoms(i,na) + oneD%ods%mrot(i,j,jop) * forcetot(j,itype)
                END IF
             END DO
          END DO
          forceAllAtoms(:,na) = MATMUL(cell%amat,forceAllAtoms(:,na)) ! to external coordinates
       END DO
    END DO
    OPEN (55,file="struct_force.xsf",status='replace')
    CALL xsf_WRITE_atoms(55,atoms,input%film,.false.,cell%amat,forceAllAtoms)
    CLOSE (55)

    IF (lconv) THEN
       WRITE (6,'(a)') "Des woars!"
       CALL juDFT_end(" GEO Des woars ", 0) ! The 0 is temporarily. Should be mpi%irank.
    ELSE
       na = 0
       DO itype=1,atoms%ntype
          tau0_i(:,itype)=MATMUL(cell%bmat,tau0(:,itype))/tpi_const
          DO ieq = 1,atoms%neq(itype)
             na = na + 1
             jop = sym%invtab(sym%ngopr(na))
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

       input%l_f = .FALSE.

       IF(.NOT.input%l_inpXML) THEN
          CALL judft_error("Relaxation no longer supported with old inp-file") 
       ELSE
          kpts_temp%ntet = 1
          kpts_temp%numSpecialPoints = 1
          ALLOCATE(kpts_temp%specialPoints(3,kpts_temp%numSpecialPoints))
          ALLOCATE(noel_temp(1),atomTypeSpecies(1),speciesRepAtomType(1))
          ALLOCATE(xmlElectronStates(1,1),xmlPrintCoreStates(1,1))
          ALLOCATE(xmlCoreOccs(1,1,1))
          CALL initWannierDefaults(wann_temp)
          CALL r_inpXML(atoms_temp,vacuum_temp,input_temp,stars_temp,sliceplot_temp,&
                        banddos_temp,dimension_temp,forcetheo,field,cell_temp,sym_temp,xcpot_temp,noco_temp,&
                        oneD_temp,hybrid_temp,kpts_temp,enpara_temp,coreSpecInput_temp,wann_temp,noel_temp,&
                        namex_temp,relcor_temp,a1_temp,a2_temp,a3_temp,dtild_temp,xmlElectronStates,&
                        xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType,l_kpts_temp)
          numSpecies = SIZE(speciesRepAtomType)
          filename = 'inp_new.xml'
          input_temp = input
          !input_temp%l_f = input%l_f
          input_temp%gw_neigd = dimension_temp%neigd
          div(:) = MIN(kpts_temp%nkpt3(:),1)
          stars_temp%gmax = stars_temp%gmaxInit
          CALL w_inpXML(atoms_new,vacuum_temp,input_temp,stars_temp,sliceplot_temp,forcetheo,&
                        banddos_temp,cell_temp,sym_temp,xcpot_temp,noco_temp,oneD_temp,hybrid_temp,&
                        kpts_temp,kpts_temp%nkpt3,kpts_temp%l_gamma,noel_temp,namex_temp,relcor_temp,a1_temp,a2_temp,a3_temp,&
                        dtild_temp,input_temp%comment,xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
                        atomTypeSpecies,speciesRepAtomType,.FALSE.,filename,.TRUE.,numSpecies,enpara_temp)
          DEALLOCATE(atomTypeSpecies,speciesRepAtomType)
          DEALLOCATE(xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs)
       END IF

    ENDIF
    RETURN
  END SUBROUTINE geo
END MODULE m_geo
