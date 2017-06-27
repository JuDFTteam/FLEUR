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
    USE m_constants
    USE m_rinpXML
    USE m_winpXML
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
    CALL bfgs(atoms%ntype,istep,istep0,forcetot,&
         &          zat,input%xa,input%thetad,input%epsdisp,input%epsforce,tote,&
         &          xold,y,h,tau0, lconv)

    IF (lconv) THEN
       WRITE (6,'(a)') "Des woars!"
       CALL juDFT_end(" GEO Des woars ", 0) ! The 0 is temporarily. Should be mpi%irank.
    ELSE
       na = 0
       DO itype=1,atoms%ntype
          tau0_i(:,itype)=MATMUL(cell%bmat,tau0(:,itype))/tpi_const
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

       input%l_f = .FALSE.

       IF(.NOT.input%l_inpXML) THEN
          ALLOCATE(atoms_temp%nz(atoms%ntype))
          ALLOCATE(atoms_temp%zatom(atoms%ntype))
          ALLOCATE(atoms_temp%jri(atoms%ntype))
          ALLOCATE(atoms_temp%dx(atoms%ntype))
          ALLOCATE(atoms_temp%lmax(atoms%ntype))
          ALLOCATE(atoms_temp%nlo(atoms%ntype))
          ALLOCATE(atoms_temp%ncst(atoms%ntype))
          ALLOCATE(atoms_temp%lnonsph(atoms%ntype))
          ALLOCATE(atoms_temp%nflip(atoms%ntype))
          ALLOCATE(atoms_temp%l_geo(atoms%ntype))
          ALLOCATE(atoms_temp%lda_u(atoms%ntype))
          ALLOCATE(atoms_temp%bmu(atoms%ntype))
          ALLOCATE(atoms_temp%relax(3,atoms%ntype))
          ALLOCATE(atoms_temp%neq(atoms%ntype))
          ALLOCATE(atoms_temp%taual(3,atoms%nat))
          ALLOCATE(atoms_temp%pos(3,atoms%nat))
          ALLOCATE(atoms_temp%rmt(atoms%ntype))

          ALLOCATE(atoms_temp%ncv(atoms%ntype))
          ALLOCATE(atoms_temp%ngopr(atoms%nat))
          ALLOCATE(atoms_temp%lapw_l(atoms%ntype))
          ALLOCATE(atoms_temp%invsat(atoms%nat))

          ALLOCATE(noco_temp%soc_opt(atoms%ntype+2),noco_temp%l_relax(atoms%ntype),noco_temp%b_con(2,atoms%ntype))
          ALLOCATE(noco_temp%alph(atoms%ntype),noco_temp%beta(atoms%ntype))

          ALLOCATE (Jij_temp%alph1(atoms%ntype),Jij_temp%l_magn(atoms%ntype),Jij_temp%M(atoms%ntype))
          ALLOCATE (Jij_temp%magtype(atoms%ntype),Jij_temp%nmagtype(atoms%ntype))

          ALLOCATE(atoms_temp%llo(atoms%nlod,atoms%ntype))
          ALLOCATE(atoms_temp%ulo_der(atoms%nlod,atoms%ntype))
          ALLOCATE(atoms_temp%l_dulo(atoms%nlod,atoms%ntype))

          ALLOCATE(vacuum_temp%izlay(vacuum%layerd,2))
          atoms_temp%ntype = atoms%ntype
          ALLOCATE(noel_temp(atoms%ntype))

          ALLOCATE (hybrid_temp%nindx(0:atoms%lmaxd,atoms%ntype))
          ALLOCATE (hybrid_temp%select1(4,atoms%ntype),hybrid_temp%lcutm1(atoms%ntype))
          ALLOCATE (hybrid_temp%select2(4,atoms%ntype),hybrid_temp%lcutm2(atoms%ntype),hybrid_temp%lcutwf(atoms%ntype))

          CALL rw_inp('r',atoms_temp,obsolete_temp,vacuum_temp,input_temp,stars_temp,sliceplot_temp,&
                      banddos_temp,cell_temp,sym_temp,xcpot_temp,noco_temp,Jij_temp,oneD_temp,hybrid_temp,&
                      kpts_temp,noel_temp,namex_temp,relcor_temp,a1_temp,a2_temp,a3_temp,scale_temp,dtild_temp,&
                      input_temp%comment)
          input_temp%l_f = input%l_f
          input_temp%tkb = input%tkb
          input_temp%delgau = input%tkb
          cell_temp = cell
          sym_temp = sym
          vacuum_temp = vacuum
          CALL rw_inp('W',atoms_new,obsolete_temp,vacuum_temp,input_temp,stars_temp,sliceplot_temp,&
               banddos_temp,cell_temp,sym_temp,xcpot_temp,noco_temp,Jij_temp,oneD_temp,hybrid_temp,&
               kpts_temp,noel_temp,namex_temp,relcor_temp,a1_temp,a2_temp,a3_temp,scale_temp,a3_temp(3),&
               input_temp%comment)
    
       ELSE
          kpts_temp%numSpecialPoints = 1
          ALLOCATE(kpts_temp%specialPoints(3,kpts_temp%numSpecialPoints))
          ALLOCATE(noel_temp(1),atomTypeSpecies(1),speciesRepAtomType(1))
          ALLOCATE(xmlElectronStates(1,1),xmlPrintCoreStates(1,1))
          ALLOCATE(xmlCoreOccs(1,1,1))
          CALL r_inpXML(atoms_temp,obsolete_temp,vacuum_temp,input_temp,stars_temp,sliceplot_temp,&
                        banddos_temp,dimension_temp,cell_temp,sym_temp,xcpot_temp,noco_temp,Jij_temp,&
                        oneD_temp,hybrid_temp,kpts_temp,enpara_temp,sphhar_temp,l_opti_temp,noel_temp,&
                        namex_temp,relcor_temp,a1_temp,a2_temp,a3_temp,scale_temp,dtild_temp,xmlElectronStates,&
                        xmlPrintCoreStates,xmlCoreOccs,atomTypeSpecies,speciesRepAtomType)
          numSpecies = SIZE(speciesRepAtomType)
          filename = 'inp_new.xml'
          input_temp%l_f = input%l_f
          input_temp%gw_neigd = dimension_temp%neigd
          div(:) = MIN(kpts_temp%nkpt3(:),1)
          stars_temp%gmax = stars_temp%gmaxInit
          CALL w_inpXML(atoms_new,obsolete_temp,vacuum_temp,input_temp,stars_temp,sliceplot_temp,&
                        banddos_temp,cell_temp,sym_temp,xcpot_temp,noco_temp,jij_temp,oneD_temp,hybrid_temp,&
                        kpts_temp,kpts_temp%nkpt3,kpts_temp%l_gamma,noel_temp,namex_temp,relcor_temp,a1_temp,a2_temp,a3_temp,&
                        scale_temp,dtild_temp,input_temp%comment,xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs,&
                        atomTypeSpecies,speciesRepAtomType,.FALSE.,filename,.TRUE.,numSpecies,enpara_temp)
          DEALLOCATE(atomTypeSpecies,speciesRepAtomType)
          DEALLOCATE(xmlElectronStates,xmlPrintCoreStates,xmlCoreOccs)
       END IF

    ENDIF
    RETURN
  END SUBROUTINE geo
END MODULE m_geo
