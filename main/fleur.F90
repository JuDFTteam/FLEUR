!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur
      IMPLICIT NONE
      CONTAINS
        SUBROUTINE fleur_execute(mpi_comm)

          !     ***************************************************************
          !
          !     based on flapw7 (c.l.fu, m.weinert, e.wimmer):
          !     full potential linearized augmented plane wave method for thin
          !     films and superlattices (version 7 ---- general symmetry)
          !     symmetry part       ---  e.wimmer
          !     potential generator ---  c.l.fu,r.podloucky
          !     matrix elements     ---  m.weinert
          !     charge density      ---  c.l.fu
          !                                c.l.fu        1987
          !     2nd variation diagon.  --- r.-q. wu      1992
          !     forces a la Yu et al   --- r.podloucky   1995
          !     full relativistic core --- a.shick       1996
          !     broyden mixing         --- r.pentcheva   1996
          !     gga (pw91, pbe)        --- t.asada       1997
          !     local orbitals         --- p.kurz        1997
          !     automatic symmetry     --- w.hofer       1997
          !     core tails & start     --- r.abt         1998
          !     spin orbit coupling    --- a.shick,x.nie 1998
          !     non-colinear magnet.   --- p.kurz        1999
          !     one-dimensional        --- y.mokrousov   2002
          !     exchange parameters    --- m.lezaic      2004
          !
          !                       g.bihlmayer, s.bluegel 1999
          !     ***************************************************************
          !----------------------------------------
          ! this routine is the main PROGRAM

          USE m_types
          USE m_fleur_init
          USE m_pldngen
          USE m_optional
          USE m_vgen
          USE m_rhodirgen
          USE m_writexcstuff
          USE m_vmatgen
          USE m_icorrkeys
          USE m_eigen
          USE m_eigenso
          USE m_fermie
          USE m_force0
          USE m_cdngen
          USE m_totale
          USE m_potdis
          USE m_mix
          USE m_xmlOutput
          USE m_juDFT_time
          !          USE m_jcoff
          !          USE m_jcoff2
          !          USE m_ssomat
#ifdef CPP_WANN
          USE m_wann_optional
          USE m_wannier
#endif
          USE m_mixedbasis
          USE m_io_hybrid
          USE m_coulomb
          USE m_gen_map
          USE m_dwigner
          !          USE m_generate_pntgpt
          !          USE m_rotate_eig
          USE m_ylm
#ifdef CPP_MPI
          USE m_mpi_bc_all,  ONLY : mpi_bc_all
#endif
          USE m_eig66_io,   ONLY : open_eig, close_eig
          IMPLICIT NONE

          INTEGER,INTENT(IN) :: mpi_comm

          !     Types, these variables contain a lot of data!
          TYPE(t_input)    :: input
          TYPE(t_dimension):: dimension
          TYPE(t_atoms)    :: atoms
          TYPE(t_sphhar)   :: sphhar
          TYPE(t_cell)     :: cell
          TYPE(t_stars)    :: stars
          TYPE(t_sym)      :: sym
          TYPE(t_noco)     :: noco
          TYPE(t_vacuum)   :: vacuum
          TYPE(t_sliceplot):: sliceplot
          TYPE(t_banddos)  :: banddos
          TYPE(t_obsolete) :: obsolete
          TYPE(t_enpara)   :: enpara
          TYPE(t_xcpot)    :: xcpot
          TYPE(t_results)  :: results
          TYPE(t_jij)      :: jij
          TYPE(t_kpts)     :: kpts
          TYPE(t_hybrid)   :: hybrid
          TYPE(t_oneD)     :: oneD
          TYPE(t_mpi)      :: mpi

          TYPE(t_potden)   :: v,vx

          !     .. Local Scalars ..
          INTEGER:: eig_id
          INTEGER:: it,ithf,itype,l
          LOGICAL:: stop80,reap,l_endit,l_opti,l_cont
          !--- J<
          INTEGER             :: phn
          REAL, PARAMETER     :: tol = 1.e-8
          INTEGER             :: qcount ,imt,i_J,j_J
          !--- J>
          !     HF/hybrid-functionals/EXX
          LOGICAL               ::  l_restart
#ifdef CPP_MPI
          include 'mpif.h'
          integer:: ierr(2)
#endif
          mpi%mpi_comm=mpi_comm
         
         CALL timestart("Initialization")
         CALL fleur_init(mpi,input,dimension,atoms,sphhar,cell,stars,sym,noco,vacuum,&
                 sliceplot,banddos,obsolete,enpara,xcpot,results,jij,kpts,hybrid,&
                 oneD,l_opti)
         CALL timestop("Initialization")

          hybrid%l_hybrid   = (&
               xcpot%icorr == icorr_pbe0 .OR.&
               xcpot%icorr == icorr_hse  .OR.&
               xcpot%icorr == icorr_vhse .OR.&
               xcpot%icorr == icorr_hf   .OR.&
               xcpot%icorr == icorr_exx)

          IF (l_opti) THEN
             IF (sliceplot%iplot .AND. (mpi%irank==0) ) THEN
                IF (noco%l_noco) THEN
                   CALL pldngen(sym,stars,atoms,sphhar,vacuum,&
                        cell,input,noco,oneD,sliceplot)
                ENDIF
             ENDIF
             CALL OPTIONAL(mpi,atoms,sphhar,vacuum,dimension,&
                  stars,input,sym,cell,sliceplot,obsolete,xcpot,noco,oneD)
          ENDIF
          !
          IF (sliceplot%iplot)      CALL juDFT_end("density plot o.k.",mpi%irank)
          IF (input%strho)          CALL juDFT_end("starting density generated",mpi%irank)
          IF (input%swsp)           CALL juDFT_end("spin polarised density generated",mpi%irank)
          IF (input%lflip)          CALL juDFT_end("magnetic moments flipped",mpi%irank)
          IF (input%l_bmt)          CALL juDFT_end('"cdnbmt" written',mpi%irank)


#ifdef CPP_WANN
          input%l_wann = .FALSE.
          INQUIRE (file='wann_inp',exist=input%l_wann)
          IF (input%l_wann .AND. (mpi%irank == 0))THEN
             CALL wann_optional(input,atoms, sym,cell,oneD,noco)
          ENDIF
#endif

          l_restart = .TRUE.

          it     = 0
          ithf   = 0
          l_cont = ( it < input%itmax )
          results%last_distance = -1.0
          IF (mpi%irank.EQ.0) CALL openXMLElementNoAttributes('scfLoop')
          DO 80 WHILE ( l_cont )
             it = it + 1
             !+t3e
             IF (input%alpha.LT.10.0) THEN
                !
                IF (it.GT.1) THEN
                   obsolete%pot8 = .FALSE.
                   input%alpha = input%alpha - NINT(input%alpha)
                END IF
                !
                CALL resetIterationDependentTimers()
                CALL timestart("Iteration")
                IF (mpi%irank.EQ.0) THEN
                   !-t3e
                   WRITE (6,FMT=8100) it
                   WRITE (16,FMT=8100) it
8100               FORMAT (/,10x,'   it=    ',i5)
                   !
                   IF (.NOT.obsolete%pot8) THEN
                      !
                      !      ----> potential generator
                      !
                      !---> pk non-collinear
                      !--->        reload the density matrix from file rhomat_in
                      !--->        calculate spin-up and -down density for USE in the
                      !--->        potential generator and store the direction of
                      !--->        magnetization on file dirofmag
                      IF (noco%l_noco) THEN
                         CALL timestart("gen. spin-up and -down density")
                         CALL rhodirgen(dimension,sym,stars,atoms,sphhar,vacuum,22,cell,input,oneD)
                         CALL timestop("gen. spin-up and -down density")
                      ENDIF
                      !---> pk non-collinear

                      reap=.NOT.obsolete%disp
                      input%total = .TRUE.
                   ENDIF!(obsolete%pot8)
                ENDIF !mpi%irank.eq.0
#ifdef CPP_MPI
                CALL MPI_BCAST(input%total,1,MPI_LOGICAL,0,mpi%mpi_comm,ierr)
#endif

                !--- J<
                IF(jij%l_jenerg) GOTO 234

                jij%alph1(:)=noco%alph(:)
                stop80= .FALSE.
                IF ( (noco%l_soc .AND. noco%l_ss) ) THEN
                   IF ( (jij%l_J).OR.(jij%nqpt/=1).OR.(jij%nmagn/=1).OR.(jij%phnd/=1) ) THEN
                      CALL juDFT_error("fleur: J-loop with ss+soc" ,calledby ="fleur")
                   ENDIF
                ENDIF
                DO qcount=1,jij%nqpt
                   IF (jij%l_J) THEN
                      noco%qss(:)=jij%qj(:,qcount)
                      jij%qn = ( noco%qss(1)**2 + noco%qss(2)**2 + noco%qss(3)**2 )
                   ENDIF
                   IF ( jij%l_J.AND.(mpi%irank.EQ.0) ) THEN
                      WRITE(6,*) 'qss=(',noco%qss(1),',',noco%qss(2),',',noco%qss(3),')'
                      CALL timestart("Q-point for J_ij(total)")

                   ENDIF
                   !HF
                   hybrid%l_subvxc = ( hybrid%l_hybrid.and.xcpot%icorr /= icorr_exx )

                   IF ( noco%l_soc ) THEN
                      dimension%neigd2 = dimension%neigd*2
                   ELSE
                      dimension%neigd2 = dimension%neigd
                   END IF
                   IF( .NOT. ALLOCATED(results%w_iks) )&
                        ALLOCATE ( results%w_iks(dimension%neigd2,kpts%nkpt,dimension%jspd) )

#ifdef CPP_NEVER
                   IF(  hybrid%l_hybrid .AND. it == 1 ) THEN
                      CALL juDFT_WARN ("Hybrid functionals not working in this version")
                      CALL timestart("generation of mixedbasis and coulombmatrix")

                      IF ( mpi%irank == 0 ) WRITE(*,'(/A)',advance='no') ' calculation of mixedbasis...'
                      print *,"symcheck:",sym%invs,noco%l_noco 
                      eig_id=open_eig(&
                      mpi%mpi_comm,dimension%nbasfcn,dimension%neigd,kpts%nkpt,dimension%jspd,atoms%lmaxd,atoms%nlod,atoms%ntype,atoms%nlotot&
                      ,noco%l_noco,.FALSE.,sym%invs.AND..NOT.noco%l_noco,noco%l_soc,.FALSE.,&
         mpi%n_size,layers=vacuum%layers,nstars=vacuum%nstars,ncored=DIMENSION%nstd,&
         nsld=atoms%nat,nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,l_mcd=banddos%l_mcd,&
         l_orb=banddos%l_orb)
                      if (kpts%nkptf==0) call judft_error("kpoint-set of full BZ not available",hint="to generate kpts in the full BZ you should specify a k-mesh in inp.xml")
                      !construct the mixed-basis
                      CALL mixedbasis(atoms,kpts, dimension,input,cell,sym,xcpot,hybrid, eig_id,mpi,v,l_restart)
                      IF ( mpi%irank == 0 ) WRITE(*,'(A)')'...done'
                      hybrid%maxlmindx = MAXVAL((/ ( SUM( (/ (hybrid%nindx(l,itype)*(2*l+1), l=0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )
                    
                      call open_hybrid_io(hybrid,dimension,atoms,sym%invs)
                      
                      IF ( mpi%irank == 0 ) WRITE(*,'(A)',advance='no') ' calculation of coulomb matrix ...'
                      CALL coulombmatrix(mpi,atoms,kpts,cell,sym,hybrid,xcpot,l_restart)
                      IF ( mpi%irank == 0 ) WRITE(*,'(A)')'...done'

#ifdef CPP_MPI
                      CALL MPI_Bcast( hybrid%maxbasm1,1,MPI_INTEGER4,0, mpi%mpi_comm,ierr(1) )
                      CALL MPI_Bcast( hybrid%radshmin,1,MPI_REAL8,   0, mpi%mpi_comm,ierr(1) )
#endif
                      CALL timestop("generation of mixedbasis and coulombmatrix")


                      IF ( noco%l_soc ) THEN
                         input%zelec = input%zelec * 2
                      END IF

                      IF ( mpi%irank == 0 ) WRITE(*,'(A)',advance='no') ' start fermie....'
                      CALL fermie(eig_id, mpi,kpts,obsolete,&
                           input,noco,enpara%epara_min,jij,cell,results)

                      IF ( noco%l_soc ) THEN
                         input%zelec = input%zelec / 2
                      END IF

                      IF ( mpi%irank == 0 ) WRITE(*,'(A)') '...done'

                   ELSEIF ( it == 1 ) THEN ! allocate some dummy arrays
                   IF (it==1) THEN ! temporary until HF is not excluded by #if any more
                      IF ( noco%l_soc ) THEN
                         dimension%neigd2 = dimension%neigd*2
                      ELSE
                         dimension%neigd2 = dimension%neigd
                      END IF
                      kpts%nkptf = 0; hybrid%maxindx = 0; hybrid%gptmd = 0; hybrid%maxgptm = 0
                      hybrid%maxgptm1 = 0; hybrid%maxgptm2 = 0; hybrid%maxindxm1 = 0; hybrid%maxindxm2 = 0
                      hybrid%maxlcutm1 = 0; hybrid%maxlcutm2 = 0; hybrid%maxindxp1 = 0; hybrid%maxindxp2 = 0
                      ALLOCATE(hybrid%gptm(0,0),hybrid%ngptm(0),hybrid%pgptm(0,0),hybrid%ngptm1(0),&
                           hybrid%pgptm1(0,0),hybrid%ngptm2(0),hybrid%pgptm2(0,0),hybrid%basm1(0,0,0,0),&
                           hybrid%basm2(0,0,0,0),hybrid%nindxm1(0,0),hybrid%nindxm2(0,0))
                   END IF ! first iteration hybrids
#endif
                   !HF
                   !#endif
                   IF (.NOT.obsolete%pot8) THEN
                      CALL timestart("generation of potential")
                      CALL vgen(reap,input,xcpot,dimension, atoms,sphhar,stars,vacuum,&
                           sym,obsolete,cell, oneD,sliceplot,mpi ,results,noco,v,vx)
                      CALL timestop("generation of potential")

                      IF (mpi%irank.EQ.0) THEN
                         !---> pk non-collinear
                         !--->          generate the four component matrix potential from spin up
                         !--->          and down potentials and direction of the magnetic field
                         IF (noco%l_noco) THEN
                            CALL timestart("generation of potential-matrix")
                            CALL vmatgen(stars, atoms,sphhar,vacuum,sym,input,oneD,8,22,26)
                            CALL timestop("generation of potential-matrix")
                         ENDIF
                         !---> end pk non-collinear
                         !---> do some output for the tddft calculations:
                         IF (input%gw /= 0) THEN
                            CALL write_xcstuff(sphhar,atoms,dimension,sym, stars,vacuum,input)
                         ENDIF
                         !
                      ENDIF ! mpi%irank.eq.0

                      !
                      !+t3e
                   ENDIF ! .not.obsolete%pot8
                   IF(  hybrid%l_hybrid .AND. it == 1 ) THEN
                      CALL juDFT_WARN ("Hybrid functionals not working in this version")
                      CALL timestart("generation of mixedbasis and coulombmatrix")

                      IF ( mpi%irank == 0 ) WRITE(*,'(/A)',advance='no') ' calculation of mixedbasis...'
                      eig_id=open_eig(&
                      mpi%mpi_comm,dimension%nbasfcn,dimension%neigd,kpts%nkpt,dimension%jspd,atoms%lmaxd,atoms%nlod,atoms%ntype,atoms%nlotot&
                      ,noco%l_noco,.FALSE.,sym%invs.AND..NOT.noco%l_noco,noco%l_soc,.FALSE.,&
         mpi%n_size,layers=vacuum%layers,nstars=vacuum%nstars,ncored=DIMENSION%nstd,&
         nsld=atoms%nat,nat=atoms%nat,l_dos=banddos%dos.OR.input%cdinf,l_mcd=banddos%l_mcd,&
         l_orb=banddos%l_orb)
                      if (kpts%nkptf==0) call judft_error("kpoint-set of full BZ not available",hint="to generate kpts in the full BZ you should specify a k-mesh in inp.xml")
                      !construct the mixed-basis
                      CALL mixedbasis(atoms,kpts, dimension,input,cell,sym,xcpot,hybrid, eig_id,mpi,v,l_restart)
                      IF ( mpi%irank == 0 ) WRITE(*,'(A)')'...done'
                      hybrid%maxlmindx = MAXVAL((/ ( SUM( (/ (hybrid%nindx(l,itype)*(2*l+1), l=0,atoms%lmax(itype)) /) ),itype=1,atoms%ntype) /) )
                      
                      call open_hybrid_io(hybrid,dimension,atoms,sym%invs)
                      IF ( mpi%irank == 0 ) WRITE(*,'(A)',advance='no') ' calculation of coulomb matrix ...'
                      CALL coulombmatrix(mpi,atoms,kpts,cell,sym,hybrid,xcpot,l_restart)
                      IF ( mpi%irank == 0 ) WRITE(*,'(A)')'...done'

#ifdef CPP_MPI
                      CALL MPI_Bcast( hybrid%maxbasm1,1,MPI_INTEGER4,0, mpi%mpi_comm,ierr(1) )
                      CALL MPI_Bcast( hybrid%radshmin,1,MPI_REAL8,   0, mpi%mpi_comm,ierr(1) )
#endif
                      CALL timestop("generation of mixedbasis and coulombmatrix")
                   ENDIF
                   
#ifdef CPP_MPI
                   CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

                   !
                   !          ----> eigenvalues and eigenfunctions
                   !
                   !--- J<
                   IF(jij%l_disp)THEN
                      jij%mtypes=1
                      jij%nmagn=1
                      jij%phnd=1
                   ENDIF
                   i_J=1
                   DO imt=1,jij%mtypes
                      DO j_J=i_J,jij%nmagn
                         DO phn=1,jij%phnd


                            input%eigvar(1)= .TRUE.
                            input%eigvar(2)= .TRUE.
                            input%eigvar(3)= .TRUE.

                            input%eigvar(2)= input%eigvar(2) .AND. ( noco%l_soc .AND. (.NOT.noco%l_noco) )
                            ! eigvar(1/2)= 1st/2nd var. ; eigvar(3)= calc density,etc

                            IF ( noco%l_soc .AND. (.NOT.noco%l_noco) ) THEN
                               input%evonly(1)= .FALSE.
                               input%evonly(2)= input%eonly
                            ELSE
                               input%evonly(1)= input%eonly
                               input%evonly(2)= .FALSE.
                            ENDIF

                            IF ( input%eigvar(1).OR.input%eigvar(2) ) THEN
                               IF (jij%l_J) THEN
                                  input%tkb=0.
#ifdef CPP_NEVER
                                  CALL jcoff(i_J,j_J,phn,mpi,atoms,atoms, noco,jij)
#endif
                               ENDIF
                               IF (input%eigvar(1)) THEN
                                  CALL timestart("generation of hamiltonian and diagonalization (total)")

                                  ! WRITE(6,fmt='(A)') 'Starting 1st variation ...'
                                  CALL timestart("eigen")
                                  CALL eigen(mpi,stars,sphhar,atoms,obsolete,xcpot,&
                                       sym,kpts,dimension,vacuum,input,cell,enpara,banddos,noco,jij,oneD,hybrid,&
                                       it,eig_id, results,v,vx)
                                  CALL timestop("eigen")
                                  !
                                  !                   add all contributions to total energy
                                  !
                                  IF( hybrid%l_subvxc) THEN
                                     DEALLOCATE( results%w_iks )
#ifdef CPP_MPI
                                     ! send all result of local total energies to the r
                                     IF (mpi%irank==0) THEN
                                        CALL MPI_Reduce(MPI_IN_PLACE,results%te_hfex%valence,&
                                             1,MPI_REAL8,MPI_SUM,0,mpi%mpi_comm,ierr(1))
                                        CALL MPI_Reduce(MPI_IN_PLACE,results%te_hfex%core,&
                                             1,MPI_REAL8,MPI_SUM,0,mpi%mpi_comm,ierr(1))
                                     ELSE
                                        CALL MPI_Reduce(results%te_hfex%valence,MPI_IN_PLACE,&
                                             1,MPI_REAL8,MPI_SUM,0, mpi%mpi_comm,ierr(1))
                                        CALL MPI_Reduce(results%te_hfex%core,MPI_IN_PLACE,&
                                             1,MPI_REAL8,MPI_SUM,0, mpi%mpi_comm,ierr(1))
                                     ENDIF
                                     !                                  END IF
#endif
                                  END IF ! xcpot%icorr = any hybrid



                            ENDIF
                            IF (input%eigvar(2)) THEN
                               ! RS: open unit for SOC vectors for GW
                               IF(noco%l_soc.AND.input%gw.EQ.2) THEN
                                  WRITE(6,'(A)') 'RS: open SOCVEC unit 4649'
                                  OPEN(4649,file='SOCVEC',form='unformatted')
                               ENDIF
                               ! WRITE(6,fmt='(A)') 'Starting 2nd variation ...'
                               CALL eigenso(eig_id,mpi,dimension,stars,vacuum,atoms,sphhar,&
                                    obsolete,sym,cell,noco,input,kpts, oneD)
                               IF(noco%l_soc.AND.input%gw.EQ.2) THEN
                                  CLOSE(4649)
                                  INQUIRE(1014,opened=l_endit)
                                  IF(l_endit) CLOSE(1014)
                                  INQUIRE(667,opened=l_endit)
                                  IF(l_endit) CLOSE(667)
                                  CALL juDFT_end("GW+SOC finished",mpi%irank)
                               ENDIF
                            ENDIF
                            CALL timestop("generation of hamiltonian and diagonalization (total)")

#ifdef CPP_MPI
                            CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif

                         ENDIF ! ( input%eigvar(1) .OR. input%eigvar(2) )

                         IF ( input%eigvar(3) .AND. noco%l_soc .AND. noco%l_ss ) THEN
#ifdef CPP_NEVER
                            CALL ssomat(eig_id, mpi,dimension,stars,vacuum,atoms,sphhar,&
                                 sym,cell,noco,input,obsolete,kpts,oneD,MPI_DOUBLE_PRECISION )
#endif
                            stop80= .TRUE.
                         ENDIF
                         !-t3e
                         !
                         !              ----> fermi level and occupancies
                         !

                         IF ( input%eigvar(3) .AND. ( .NOT.(noco%l_soc .AND. noco%l_ss) ) ) THEN
                            IF (jij%l_J) THEN
                               CALL timestart("determination of fermi energy")
                               ALLOCATE ( results%w_iks(dimension%neigd,kpts%nkpt,dimension%jspd) )
                               CALL fermie(eig_id, mpi,kpts,obsolete,input,&
                                    noco,enpara%epara_min,jij,cell,results)
                               DEALLOCATE ( results%w_iks )
                               CALL timestop("determination of fermi energy")
                            ENDIF
                            IF ( noco%l_soc .AND. (.NOT. noco%l_noco) ) dimension%neigd = 2*dimension%neigd
                            IF( .NOT. ALLOCATED(results%w_iks) )&
                                 ALLOCATE ( results%w_iks(dimension%neigd,kpts%nkpt,dimension%jspd) )
                            IF ( (mpi%irank.EQ.0).AND.(.NOT.jij%l_J) ) THEN
                               CALL timestart("determination of fermi energy")

                               IF ( noco%l_soc .AND. (.NOT. noco%l_noco) ) THEN
                                  input%zelec = input%zelec*2
                                  CALL fermie(eig_id,mpi,kpts,obsolete,&
                                       input,noco,enpara%epara_min,jij,cell,results)
                                  results%seigscv = results%seigscv/2
                                  results%ts = results%ts/2
                                  input%zelec = input%zelec/2
                               ELSE
                                  CALL fermie(eig_id,mpi,kpts,obsolete,&
                                       input,noco,enpara%epara_min,jij,cell,results)
                               ENDIF
                               CALL timestop("determination of fermi energy")

                            ENDIF

                            IF (input%eonly) THEN
                               CALL close_eig(eig_id)

                               IF (.NOT. jij%l_J) THEN
                                  DEALLOCATE( results%w_iks )
#ifdef CPP_MPI
                                  CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif
                                  IF (mpi%irank==0) THEN
                                     WRITE (*,fmt='(A)') 'eigenvalues written, program stops'
                                  ENDIF
                                  stop80= .TRUE.
                               ENDIF
                            ENDIF ! input%eonly

                         ENDIF ! ( input%eigvar(3) .and. .not.(noco%l_soc .and. noco%l_ss) )
                         !--- J<
                         IF(jij%l_J) THEN
                            IF (.NOT. input%eonly) THEN
                               DEALLOCATE ( results%w_iks )
                            ENDIF
                            IF (((i_J.EQ.j_J)).OR.(sym%invs.AND.(jij%qn.GT.tol))) GOTO 33
                         ENDIF
                      ENDDO !phn
33                    CONTINUE
                   ENDDO !j_J
                   i_J=i_J+jij%nmagtype(imt)
                ENDDO !imt
                IF ((mpi%irank.EQ.0).AND.(jij%l_J)) THEN
                   CALL timestop("Q-point for J_ij(total)")
                ENDIF
             ENDDO !qcount
             IF (stop80) THEN
                IF ((mpi%irank.EQ.0).AND.(isCurrentXMLElement("iteration"))) THEN
                   CALL closeXMLElement('iteration')
                END IF
                EXIT ! it
             ENDIF

234          CONTINUE

             IF (mpi%irank.EQ.0) THEN
                IF(jij%l_J) THEN
                   IF(.NOT.jij%l_disp)THEN
                      REWIND(113)
                      REWIND(114)
#ifdef CPP_NEVER
                      CALL jcoff2(atoms,sym,cell,jij,input)
#endif
                   ENDIF
                   CLOSE(113)
                   CLOSE(114)
                ENDIF
             ENDIF

             IF (.NOT.jij%l_J) THEN
                !--- J>
#ifdef CPP_MPI
                CALL MPI_BCAST(results%ef,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
                CALL MPI_BCAST(results%w_iks,SIZE(results%w_iks),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
#endif
                !-t3e
                !
                !        ----> initialise force_old
                !
                CALL force_0(results)
                !
                !        ----> charge density
                !
                !+Wannier functions
#ifdef CPP_WANN
                input%l_wann = .FALSE.
                INQUIRE (file='wann_inp',exist=input%l_wann)
                IF (input%l_wann) THEN
                   CALL wannier(mpi,atoms,noco, dimension,sym,obsolete,cell,kpts,&
                        stars,oneD,vacuum,sphhar,input, sliceplotresults)
                ENDIF
#endif
                !-Wannier
                CALL timestart("generation of new charge density (total)")

                CALL cdngen(eig_id,mpi,input,banddos,sliceplot,vacuum,&
                     dimension,kpts,atoms,sphhar,stars,sym,obsolete,&
                     enpara,cell,noco,jij,results,oneD)
                !
                ! the w_iks are needed for hybrid functionals so do not
                ! deallocate them in that case

                IF ( hybrid%l_hybrid ) THEN
                   DEALLOCATE ( results%w_iks )
                END IF

                IF ( noco%l_soc .AND. (.NOT. noco%l_noco) ) dimension%neigd=dimension%neigd/2
                !+t3e
#ifdef CPP_MPI
                CALL MPI_BCAST(enpara%evac0,SIZE(enpara%evac0),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
                CALL MPI_BCAST(enpara%el0,SIZE(enpara%el0),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
                CALL MPI_BCAST(enpara%ello0,SIZE(enpara%ello0),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)

                IF (noco%l_noco) THEN
                   DO n= 1,atoms%ntype
                      IF (noco%l_relax(n)) THEN
                         CALL MPI_BCAST(noco%alph(n),1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
                         CALL MPI_BCAST(noco%beta(n),1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
                      ENDIF
                   ENDDO
                   IF (noco%l_constr) THEN
                      CALL MPI_BCAST(noco%b_con,SIZE(noco%b_con),MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
                   ENDIF
                ENDIF
#endif
                CALL timestop("generation of new charge density (total)")
                IF (mpi%irank.EQ.0) THEN
                   !-t3e

                   IF (banddos%ndir.GT.0) THEN
                      CALL juDFT_end("NDIR",mpi%irank)
                   END IF
                   !          ----> output potential and potential difference
                   IF (obsolete%disp) THEN
                      reap = .FALSE.
                      input%total = .FALSE.
                      CALL timestart("generation of potential (total)")
                      CALL vgen(reap,input,xcpot,dimension, atoms,sphhar,stars,vacuum,sym,&
                           obsolete,cell,oneD,sliceplot,mpi, results,noco,v,vx)
                      CALL timestop("generation of potential (total)")

                      CALL potdis(stars,vacuum,atoms,sphhar, input,cell,sym)
                   END IF
                   !
                   !i         ----> total energy
                   !


                   CALL timestart('determination of total energy')
                   CALL totale(atoms,sphhar,stars,vacuum,dimension,&
                        sym,input,noco,cell,oneD,xcpot,hybrid,it,results)

                   CALL timestop('determination of total energy')


                   ! in case of parallel processing, the total energy calculation is done
                   ! only for irank.eq.0, since no parallelization is required here. once
                   ! a force calculation is applied, however, the irank.eq.0 process is
                   ! led into a MPI_FINALIZE after convergence, while the other processes
                   ! are not, resulting in fleur not terminating despite having finished
                   ! the calculation. the next 7 lines correct that issue.
                   ! (other files subject to this correction: geo.F, force_w.F)
                   ! Schlipf/Klueppelberg Jun 2012
#ifdef CPP_MPI
                ELSEIF (input%l_f) THEN ! forces, but mpi%irank.ne.0
                   !This does not work, you can not call MPI_BCAST within a
                   !else part of irank==0 as PE=0 will not call this!
                   !CALL MPI_BCAST(lconv,1,MPI_LOGICAL,0,mpi_comm,ierr)
                   !IF (lconv) THEN
                   !  CALL MPI_FINALIZE(ierr)
                   !END IF
#endif

                ENDIF ! mpi%irank.EQ.0
                !Close file if not a hybrid calculation
                IF ( hybrid%l_hybrid ) CALL close_eig(eig_id)

             ENDIF !(if not jij%l_J)
          ELSE
             input%alpha = input%alpha - 10.
          END IF !(if input%alpha <10.)
          IF (.NOT.jij%l_J) THEN

             IF (mpi%irank.EQ.0) THEN
                !-t3e
                !
                !          ----> mix input and output densities
                !
                CALL timestart("mixing")
                CALL mix(stars,atoms,sphhar,vacuum,input,sym,cell,it,noco,oneD,hybrid,results)
                !
                CALL timestop("mixing")
                WRITE (6,FMT=8130) it
                WRITE (16,FMT=8130) it
8130            FORMAT (/,5x,'******* it=',i3,'  is completed********',/,/)
                write(*,*) "Iteration:",it," Distance:",results%last_distance
                CALL timestop("Iteration")
                !+t3e
             ENDIF ! mpi%irank.EQ.0
             ! hybrid functionals - ddist is needed on all processes
#        ifdef CPP_MPI
             CALL MPI_BCAST(hybrid%ddist,dimension%jspd,MPI_REAL8,0,mpi%mpi_comm,ierr)
#        endif

             !--- J<
          ELSE
          ENDIF !(if not jij%l_J)
          !--- J>

#ifdef CPP_MPI
          CALL MPI_BCAST(results%last_distance,1,MPI_DOUBLE_PRECISION,0,mpi%mpi_comm,ierr)
          CALL MPI_BARRIER(mpi%mpi_comm,ierr)
#endif
          !-t3e
          ! Delete the broyden files after the fifth iteration
          ! in the case of a HF or hybrid functional calculation
          IF (it.EQ. 5 .AND. hybrid%l_subvxc .AND. input%imix .LE. 10) THEN
             CALL system('rm -f broyd*')
          END IF
          !+fo
          call priv_geo_end(mpi)

          !-fo
          IF ( hybrid%l_calhf ) ithf = ithf + 1
          IF ( hybrid%l_subvxc ) THEN
             l_cont = ( ithf < input%itmax )
             results%te_hfex%core    = 0
             results%te_hfex%valence = 0
          ELSE
             l_cont = ( it < input%itmax )
          END IF
          CALL writeTimesXML()
          CALL check_time_for_next_iteration(it,l_cont)
          l_cont=l_cont.AND.((input%mindistance<=results%last_distance).OR.input%l_f)
          IF ((mpi%irank.EQ.0).AND.(isCurrentXMLElement("iteration"))) THEN
             CALL closeXMLElement('iteration')
          END IF
80        CONTINUE
          IF (mpi%irank.EQ.0) CALL closeXMLElement('scfLoop')
          CALL juDFT_end("all done",mpi%irank)

        contains
          subroutine priv_geo_end(mpi)
            TYPE(t_mpi),INTENT(IN)::mpi
            LOGICAL :: l_exist
            !Check if a new input was generated
            INQUIRE (file='inp_new',exist=l_exist)
            IF (l_exist) THEN
               CALL juDFT_end(" GEO new inp created ! ",mpi%irank)
            END IF
            !check for inp.xml
            INQUIRE (file='inp_new.xml',exist=l_exist)
            IF (.NOT.l_exist) return
            IF (mpi%irank==0) then
               CALL system('mv inp.xml inp_old.xml')
               CALL system('mv inp_new.xml inp.xml')
               INQUIRE (file='qfix',exist=l_exist)
               IF (l_exist) THEN
                  OPEN(2,file='qfix')
                  WRITE(2,*)"F"
                  CLOSE(2)
                  print *,"qfix set to F"
               ENDIF
               INQUIRE(file='broyd',exist=l_exist)
               IF (l_exist) THEN
                  CALL system('rm broyd')
                  print *,"broyd file deleted"
               ENDIF
               INQUIRE(file='broyd.7',exist=l_exist)
               IF (l_exist) THEN
                  CALL system('rm broyd.7')
                  print *,"broyd.7 file deleted"
               ENDIF
            ENDIF
            CALL juDFT_end(" GEO new inp.xml created ! ",mpi%irank)
          end subroutine priv_geo_end

        END SUBROUTINE fleur_execute
      END MODULE
