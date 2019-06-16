  TYPE(t_input)    :: input
      TYPE(t_atoms)    :: atoms
      TYPE(t_cell)     :: cell
      TYPE(t_sym)      :: sym
      TYPE(t_noco)     :: noco
      TYPE(t_vacuum)   :: vacuum
      TYPE(t_banddos)  :: banddos
      TYPE(t_hybrid)   :: hybrid
      TYPE(t_xcpot_inbuild_nf)::xcpot
      TYPE(t_oned)     :: oned
      TYPE(t_sliceplot):: sliceplot
      TYPE(t_stars)    :: stars


        TYPE(t_enpara)   :: enpara


        TYPE(t_forcetheo):: forcetheo
      TYPE(t_kpts)     :: kpts
    


      call fleur_init_old(&
       input,DIMENSION,atoms,sphhar,cell,stars,sym,noco,vacuum,forcetheo,&
       sliceplot,banddos,enpara,xcpot,kpts,hybrid,&
       oneD,coreSpecInput)
