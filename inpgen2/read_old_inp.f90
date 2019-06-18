MODULE m_read_old_inp
  IMPLICIT NONE
CONTAINS
  SUBROUTINE read_old_inp(input,atoms,cell,stars,sym,noco,vacuum,forcetheo,&
       sliceplot,banddos,enpara,xcpot,kpts,hybrid,&
       oneD)
    USE m_types_input
    USE m_types_atoms
    USE m_types_cell
    USE m_types_sym
    USE m_types_noco
    USE m_types_vacuum
    USE m_types_banddos
    USE m_types_hybrid
    USE m_types_xcpot_inbuild_nofunction
    USE m_types_forcetheo
    USE m_types_kpts
    USE m_types_enpara
    USE m_types_oneD
    USE m_types_sliceplot
    USE m_types_stars
    USE m_fleur_init_old
    USE m_atompar
    USE m_constants


    TYPE(t_input),INTENT(OUT):: input
    TYPE(t_atoms),INTENT(OUT):: atoms
    TYPE(t_cell) ,INTENT(OUT):: cell
    TYPE(t_sym)  ,INTENT(OUT):: sym
    TYPE(t_noco) ,INTENT(OUT):: noco
    TYPE(t_vacuum) ,INTENT(OUT):: vacuum
    TYPE(t_banddos),INTENT(OUT):: banddos
    TYPE(t_hybrid) ,INTENT(OUT):: hybrid
    TYPE(t_xcpot_inbuild_nf),INTENT(OUT)::xcpot
    TYPE(t_oned) ,INTENT(OUT):: oned
    TYPE(t_sliceplot),INTENT(OUT):: sliceplot
    TYPE(t_stars),INTENT(OUT):: stars
    TYPE(t_enpara) ,INTENT(OUT):: enpara

    TYPE(t_forcetheo):: forcetheo
    TYPE(t_kpts) ,INTENT(OUT):: kpts
   
    !local only
    TYPE(t_dimension)    ::dimension
    TYPE(t_sphhar)       ::sphhar
    TYPE(t_atompar)      :: ap
    INTEGER              :: n,grid(3)

    CALL fleur_init_old(&
         input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,&
         sliceplot,banddos,enpara,xcpot,kpts,hybrid,&
         oneD,grid) !kpt grid not used...

    CALL sym%init(cell,input%film)
    ALLOCATE(enpara%qn_el(0:3,atoms%ntype,2),enpara%qn_ello(size(enpara%ello0),atoms%ntype,2))
    enpara%qn_el=0.0
    enpara%qn_ello=0.0
 
    ALLOCATE(atoms%label(atoms%nat)); atoms%label=""
    ALLOCATE(atoms%speciesname(atoms%ntype))
    DO n=1,atoms%ntype
       WRITE(atoms%speciesname(n),"(a,a,i0)") ADJUSTL(TRIM(namat_const(atoms%nz(n)))),"-",n
    END DO

    DO n=1,atoms%ntype
       ap=find_atompar(atoms%nz(n),atoms%rmt(n))
       IF (atoms%nlo(n).NE.LEN_TRIM(ap%lo)/2) CALL judft_warn("Number of LOs changed in new version.")
       call atoms%econf(n)%init(ap%econfig)
    END DO

  
  END SUBROUTINE read_old_inp
 END MODULE m_read_old_inp
