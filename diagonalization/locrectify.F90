!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!*************************************************************
!     Find linear combinations of local orbitals that are
!     eigenfunctions of the z-reflection operation.
!     Used to exploit z-reflection symmetry when solving
!     the eigenvalue problem (film calculations).
!     Frank Freimuth, January 2007
!*************************************************************
#include"cpp_double.h"
MODULE m_locrectify
  USE m_juDFT
  CONTAINS
    SUBROUTINE locrectify(jsp,input,lapw,bkpt,atoms, kveclo, sym,cell,&
         kindlocrec,evenlocs,oddlocs,evenlocindex, oddlocindex,locrec_r,locrec_c)

      USE m_ylm
      USE m_constants,ONLY:tpi_const
      USE m_loccoeff
      USE m_types
      IMPLICIT NONE
      INTEGER,INTENT(IN)         :: jsp
      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_sym),INTENT(IN)     :: sym
      TYPE(t_cell),INTENT(IN)    :: cell
      TYPE(t_atoms),INTENT(IN)   :: atoms
      TYPE(t_lapw),INTENT(IN)    :: lapw
      REAL,INTENT(in)::bkpt(3)
      INTEGER,INTENT(in)::kveclo(atoms%nlotot)

      REAL,OPTIONAL,INTENT(out)::locrec_r(atoms%nlotot,atoms%nlotot)
      COMPLEX,OPTIONAL,INTENT(out)::locrec_c(atoms%nlotot,atoms%nlotot)

      LOGICAL,INTENT(out)::kindlocrec(atoms%nlotot)
      INTEGER,INTENT(out)::evenlocs
      INTEGER,INTENT(out)::oddlocs
      INTEGER,INTENT(out)::evenlocindex(atoms%nlotot)
      INTEGER,INTENT(out)::oddlocindex(atoms%nlotot)

      INTEGER numorec
      INTEGER invatom,zrefatom
      LOGICAL l_planetwo,l_spacetwo,l_spacefour,l_spacetwoinv
      LOGICAL l_crosstwo
      INTEGER pos,loc,nn,angmom,bas,locm,loopi,j,i
      INTEGER lm,locmm,lmm
      INTEGER  ind,vec,natom,nap,inap,locnum,ll1
      INTEGER ipiv(28),info
      COMPLEX loccoffs(28,28),reccoffs(28)
      COMPLEX ylm(16),phase
      REAL fk(3),tmk,fkp(3)
      !      complex coeffi(1:28,1:28,0:13)     transferred to mod_loccoeff.f
      !      logical coffkind(1:28,0:13)        import with USE m_loccoeff
      INTEGER jump(atoms%nat)
      INTEGER mqnum,nalo,naup,basindex(atoms%nat,atoms%nlod)
      INTEGER zrefnap,zrefinap,coffindex
      LOGICAL l_real,l_invsatswap

      !      ci=cmplx(0.0,1.0) transferred to mod_loccoeff.f

      l_real=PRESENT(locrec_r)

      jump(1:atoms%nat)=0
      IF (l_real) THEN
         locrec_r=0.0
      ELSE
         locrec_c=0.0
      ENDIF
      numorec=0

      !***********************************************
      !   basindex holds first index of subset of locs
      !***********************************************
      bas=1
      nalo=1
      DO ind=1,atoms%ntype
         naup=nalo+atoms%neq(ind)-1
         DO natom=nalo,naup
            IF(atoms%invsat(natom).EQ.2)CYCLE
            DO loc=1,atoms%nlo(ind)
               basindex(natom,loc)=bas
               locnum=(atoms%invsat(natom)+1)*(2*atoms%llo(loc,ind)+1)
               bas=bas+locnum
            ENDDO !loc
         ENDDO !natom
         nalo=naup+1
      ENDDO !ind

      bas=0
      nalo=1
      DO ind=1,atoms%ntype
         naup=nalo+atoms%neq(ind)-1
         DO natom=nalo,naup
            IF(atoms%invsat(natom).EQ.2)CYCLE
            IF(jump(natom).NE.0)THEN
               bas=bas+jump(natom)
               CYCLE
            ENDIF

            !*******************************************************
            !       find out by what kind of symmetry the atoms
            !       are interrelated
            !*******************************************************
            !two atoms related by inversion symmetry are in xy-plane
            l_planetwo=.FALSE.    
            !two atoms are related by z-reflection, but not by inversion
            l_spacetwo=.FALSE.
            !two atoms are related by z-reflection and by inversion
            l_spacetwoinv=.FALSE.
            !four atoms are entangled, by both inversion and z-reflection
            l_spacefour=.FALSE.
            !two atoms are related by inversion, but not by z-reflection
            l_crosstwo=.FALSE.
            !order of locs corresponds to setup of coeffi-array => l_invsatswap=.false.
            !order of locs is reversed with respect to coeffi-array => l_invsatswap=.true.
            l_invsatswap=.FALSE.
            IF((ABS(ABS(atoms%taual(3,natom))-0.5).LT.1e-6).AND..NOT.input%film)THEN !atom on face
               IF(atoms%invsat(natom).EQ.1)THEN ! two atoms
                  invatom=sym%invsatnr(natom)
                  l_planetwo=.TRUE.
               ENDIF
            ELSEIF(ABS(atoms%taual(3,natom)).GT.1e-6) THEN !atom not in xy-plane
               DO zrefatom=nalo,naup !find z-reflection image
                  IF(ALL(ABS(atoms%taual(1:2,natom)-atoms%taual(1:2,zrefatom))-&
                       NINT(ABS(atoms%taual(1:2,natom)-atoms%taual(1:2,zrefatom))).LT.1e-6)) THEN
                     IF(ABS(atoms%taual(3,natom)+atoms%taual(3,zrefatom)).LT.1e-6) GOTO 543
                  ENDIF
               ENDDO !zrefatom
               if (l_real) THEN
                  !if there is no z-reflected counterpart there must be an inversion 
                  !counterpart instead
                  IF(atoms%invsat(natom).EQ.1)THEN
                     zrefatom=0 !this switches on l_crosstwo later on
                  ELSE
                     CALL juDFT_error("nozref-invsat",calledby ="locrectify")
                  ENDIF
               else
                  !in the complex version there must be a z-reflected counterpart
                  CALL juDFT_error("zrefatom not found",calledby ="locrectify")
               endif
543            CONTINUE

               IF(atoms%invsat(natom).EQ.0) THEN !no inversion-image
                  l_spacetwo=.TRUE.
                  if (l_real)   CALL juDFT_error("spacetwo & INVERSION",calledby ="locrectify")
               ELSE                        !inversion-image
                  invatom=sym%invsatnr(natom)
                  IF(zrefatom.EQ.0)THEN
                     l_crosstwo=.TRUE.
                  ELSEIF(invatom.EQ.zrefatom)THEN    !inversion = z-reflection
                     l_spacetwoinv=.TRUE.
                  ELSE                           !inversion /= z-reflection
                     if (l_real) THEN
                        l_spacefour=.TRUE. !four entangled locs
                     else
                        CALL juDFT_error("spacefour & no INVERSION",calledby ="locrectify")
                     endif
                     IF(atoms%invsat(zrefatom).EQ.1)THEN

                     ELSEIF(atoms%invsat(zrefatom).EQ.2)THEN
                        l_invsatswap=.TRUE.
                        DO info=1,atoms%nat
                           IF(sym%invsatnr(info).EQ.zrefatom)GOTO 723
                        ENDDO
                        CALL juDFT_error("no atoms for zref",calledby ="locrectify")
723                     CONTINUE
                        zrefatom=info
                        IF(atoms%invsat(zrefatom)/=1) CALL juDFT_error("atoms%invsat" ,calledby ="locrectify")
                        IF(.NOT.ALL(ABS(atoms%taual(1:2,invatom)- atoms%taual(1:2,zrefatom))<1e-6)) &
                             CALL juDFT_error ("invsat-zref2",calledby ="locrectify")
                        IF(.NOT.ABS(atoms%taual(3,invatom)+atoms%taual(3,zrefatom)) <1e-6)  &
                             CALL juDFT_error("invsat-zref3", calledby="locrectify")

                     ELSEIF(atoms%invsat(zrefatom).EQ.0)THEN
                        CALL juDFT_error("spacefour: atoms-zref",calledby ="locrectify")
                     ENDIF
                  ENDIF
               ENDIF
            ELSE !atom lies in the xy-plane
               IF(atoms%invsat(natom).EQ.1)THEN ! two atoms
                  invatom=sym%invsatnr(natom)
                  l_planetwo=.TRUE.
               ENDIF
            ENDIF !symmetry partners

            DO loc=1,atoms%nlo(ind)
               angmom=atoms%llo(loc,ind)
               IF(angmom.EQ.0)THEN
                  IF(l_planetwo)THEN
                     kindlocrec(numorec+1)=.TRUE.
                     kindlocrec(numorec+2)=.TRUE.
                     IF (l_real) THEN
                        locrec_r(bas+1,numorec+1)=1.0
                        locrec_r(bas+2,numorec+2)=1.0
                     ELSE
                        locrec_c(bas+1,numorec+1)=1.0
                        locrec_c(bas+2,numorec+2)=1.0
                     ENDIF
                     bas=bas+2
                     numorec=numorec+2
                     CYCLE
                  ELSEIF(l_spacetwoinv.OR.l_crosstwo)THEN
                     coffindex=0
                     locnum=2
                  ELSEIF(l_spacetwo) THEN
                     locnum=1
                     coffindex=0
                  ELSEIF(l_spacefour)THEN
                     locnum=2
                     coffindex=6
                  ELSE
                     kindlocrec(numorec+1)=.TRUE.
                     IF (l_real) THEN
                        locrec_r(bas+1,numorec+1)=1.0
                     ELSE
                        locrec_c(bas+1,numorec+1)=1.0
                     END IF
                     bas=bas+1
                     numorec=numorec+1
                     CYCLE
                  ENDIF !l_planetwo,l_crosstwo,....
               ELSE
                  IF(l_planetwo)THEN
                     locnum=2*(2*angmom+1)
                     IF(angmom.EQ.1)THEN
                        coffindex=3
                     ELSEIF(angmom.EQ.2)THEN
                        coffindex=9
                     ELSEIF(angmom.EQ.3)THEN
                        coffindex=12
                     ELSE
                        CALL juDFT_error("angmom1",calledby ="locrectify")
                     ENDIF
                  ELSEIF(l_spacetwoinv.OR.l_crosstwo)THEN
                     locnum=2*(2*angmom+1)
                     IF(angmom.EQ.1)THEN
                        coffindex=4
                     ELSEIF(angmom.EQ.2)THEN
                        coffindex=8
                     ELSEIF(angmom.EQ.3)THEN
                        coffindex=11
                     ELSE
                        CALL juDFT_error("angmom5",calledby ="locrectify")
                     ENDIF
                  ELSEIF(l_spacetwo) THEN
                     locnum=2*angmom+1
                     IF(angmom.EQ.1)THEN
                        coffindex=4
                     ELSEIF(angmom.EQ.2)THEN
                        coffindex=8
                     ELSEIF(angmom.EQ.3)THEN
                        coffindex=11
                     ELSE
                        CALL juDFT_error("angmom2",calledby="locrectify")
                     ENDIF
                  ELSEIF(l_spacefour) THEN
                     locnum=2*(2*angmom+1)
                     IF(angmom.EQ.1)THEN
                        coffindex=7
                     ELSEIF(angmom.EQ.2)THEN
                        coffindex=10
                     ELSEIF(angmom.EQ.3)THEN
                        coffindex=13
                     ELSE
                        CALL juDFT_error("angmom3",calledby ="locrectify")
                     ENDIF
                  ELSE
                     locnum=2*angmom+1
                     IF(angmom.EQ.1)THEN
                        coffindex=1
                     ELSEIF(angmom.EQ.2)THEN
                        coffindex=2
                     ELSEIF(angmom.EQ.3)THEN
                        coffindex=5
                     ELSE
                        CALL juDFT_error("angmom",calledby ="locrectify")
                     ENDIF
                  ENDIF! l_planetwo
               ENDIF !angmom

               mqnum=2*angmom+1
               ll1=angmom*angmom-1
               loccoffs(:,:)=CMPLX(0.0,0.0)
               !**************************************************
               !        Write the expansion of the locs
               !        in terms of spherical harmonics into the
               !        loccoffs-array.
               !        First index of loccoffs: m-quantum-number.
               !        Second index of loccoffs: index of loc.
               !**************************************************
               DO locm=1,locnum
                  pos=bas+locm
                  vec=kveclo(pos)
                  !The vector of the planewave onto which the loc is matched.
                  fk(:)=bkpt(:)+(/lapw%k1(vec,jsp),lapw%k2(vec,jsp),lapw%k3(vec,jsp)/)

                  tmk=tpi_const*DOT_PRODUCT(lapw%k1(:,jsp),atoms%taual(:,natom))
                  phase = CMPLX(COS(tmk),SIN(tmk))
                  fkp=MATMUL(fk,cell%bmat)
                  CALL ylm4(3,fkp,ylm)
                  DO locmm=1,mqnum
                     lmm=ll1+locmm
                     loccoffs(locmm,locm)=ci**angmom*phase*CONJG(ylm(lmm+1))
                  ENDDO
               ENDDO !locm
               !*************************************************************
               !        If locs are entangled by symmetry, the number of locs
               !        that have to be treated at the same time is larger.
               !        => l_spacetwo: two locs, which are NOT in the
               !        xy-plane and which are not related by inversion, are
               !        entangled by z-reflexion.
               !        => l_spacefour: two pairs the partners of which are
               !        connected by inversion are interrelated by z-reflection.
               !*************************************************************
               IF(l_spacetwo.OR.l_spacefour)THEN
                  DO locm=1,locnum
                     pos=basindex(zrefatom,loc)-1+locm
                     vec=kveclo(pos)
                     fk(:)=bkpt(:)+(/lapw%k1(vec,jsp),lapw%k2(vec,jsp),lapw%k3(vec,jsp)/)

                     tmk=tpi_const*DOT_PRODUCT(lapw%k1(:,jsp),atoms%taual(:,zrefatom))
                     phase = CMPLX(COS(tmk),SIN(tmk))
                     fkp=MATMUL(fk,cell%bmat)
                     CALL ylm4(3,fkp,ylm)
                     DO locmm=1,mqnum
                        lmm=ll1+locmm
                        loccoffs(locnum+locmm,locnum+locm)= ci**angmom*phase*CONJG(ylm(lmm+1))
                     ENDDO
                  ENDDO !locm
                  jump(zrefatom)=jump(zrefatom)+locnum
                  IF(l_spacetwo)locnum=locnum*2
               ENDIF !l_spacetwo
               !*************************************************************
               !        If locs are entangled by symmetry, loccoffs has a larger
               !        number of components for a given loc.
               !        l_planetwo => two locs in xy-plane are entangled
               !        by inversion symmetry.
               !        l_spacetwoinv => two z-reflexion-symmetric atoms
               !        are entangled also by inversion.
               !        l_crosstwo => two inversion-symmetric atoms
               !        are not directly connected by z-reflexion.
               !        l_spacefour =>  two pairs the partners of which are
               !        connected by inversion are interrelated by z-reflection.
               !*************************************************************
               IF(l_planetwo.OR.l_crosstwo.OR.l_spacetwoinv .OR.l_spacefour) THEN
                  DO locm=1,locnum
                     DO locmm=1,mqnum
                        loccoffs(mqnum+locmm,locm)=(-1)**(locmm-1)* CONJG(loccoffs(mqnum+1-locmm,locm))
                     ENDDO! locmm
                  ENDDO ! locm
               ENDIF ! l_planetwo

               !******************************************************************
               !        If four atoms are entangled, something remains to be added
               !        here:
               !******************************************************************
               IF(l_spacefour)THEN
                  DO locm=1,locnum
                     DO locmm=1,mqnum
                        loccoffs(locnum+mqnum+locmm,locnum+locm)=(-1)**(locmm-1)*&
                             CONJG(loccoffs(locnum+mqnum+1-locmm,locnum+locm))
                     ENDDO! locmm
                  ENDDO ! locm
                  locnum=locnum*2
               ENDIF
               !******************************************************************
               !        Find linear combinations of locs that are eigenfunctions
               !        of the z-reflection operation by solving linear system
               !        of equations. Write the transformation into locrec-matrix.
               !******************************************************************
               CALL CPP_LAPACK_cgetrf(locnum,locnum,loccoffs,28,ipiv,info)
               IF(info /= 0)  CALL juDFT_error("cgetrf",calledby ="locrectify")
               DO loopi=1,locnum
                  numorec=numorec+1
                  kindlocrec(numorec)=coffkind(loopi,coffindex)
                  IF(l_invsatswap)THEN
                     reccoffs(1:locnum/2)=coeffi(1:locnum/2,loopi,coffindex)
                     reccoffs(locnum/2+1:locnum/4*3)= coeffi(locnum/4*3+1:locnum,loopi,coffindex)
                     reccoffs(locnum/4*3+1:locnum)= coeffi(locnum/2+1:locnum/4*3,loopi,coffindex)
                  ELSE
                     reccoffs(1:locnum)=coeffi(1:locnum,loopi,coffindex)
                  ENDIF
                  CALL CPP_LAPACK_cgetrs('N',locnum,1,loccoffs,28,ipiv, reccoffs,28,info)
                  IF(info /= 0)  CALL juDFT_error("cgetrs", calledby="locrectify")
                  IF (l_real) THEN
                     IF(ANY(ABS(AIMAG(reccoffs(1:locnum))).GT.1e-6)) THEN
                        !In the real version of the code, the hamiltonian and overlap
                        !matrices are real. Hence the transformation matrices that make
                        !the locs become eigenfunctions of the z-reflection operation
                        !have to be real also. Check this here. 
                        PRINT*,"Sorry!"
                        PRINT*,"The transformation is not purely real."
                        PRINT*,"=> Something is wrong here."
                        CALL juDFT_error("unfortunately complex",calledby ="locrectify")
                     ENDIF
                     IF(l_spacefour)THEN
                        locrec_r(bas+1:bas+locnum/2,numorec)=REAL(reccoffs(:locnum/2))
                        DO locmm=1,locnum/2
                           locrec_r(basindex(zrefatom,loc)-1+locmm,numorec)= REAL(reccoffs(locmm+locnum/2))
                        ENDDO

                     ELSE
                        locrec_r(bas+1:bas+locnum,numorec)=REAL(reccoffs(:locnum))

                     ENDIF
                  ELSE
                     IF(l_spacetwo)THEN
                        locrec_c(bas+1:bas+locnum/2,numorec)=reccoffs(1:locnum/2)
                        locrec_c(basindex(zrefatom,loc):basindex(zrefatom,loc)+locnum/2-1,numorec)=&
                             reccoffs(1+locnum/2:locnum)
                     ELSE
                        locrec_c(bas+1:bas+locnum,numorec)=reccoffs(1:locnum)
                     ENDIF
                  ENDIF
               ENDDO !loopi
               IF(l_spacetwo) locnum=locnum/2
               IF(l_spacefour) locnum=locnum/2
               bas=bas+locnum
            ENDDO !loc
         ENDDO !natom
         nalo=naup+1
      ENDDO !ind


      !*********************************************************
      !     Determine number of even and odd locs.
      !     Prepare sort-arrays.
      !*********************************************************
      evenlocs=0
      oddlocs=0
      DO loopi=1,atoms%nlotot
         IF(kindlocrec(loopi)) THEN
            evenlocs=evenlocs+1
            evenlocindex(evenlocs)=loopi
         ELSE
            oddlocs=oddlocs+1
            oddlocindex(oddlocs)=loopi
         ENDIF
      ENDDO
    END SUBROUTINE locrectify
END MODULE m_locrectify

