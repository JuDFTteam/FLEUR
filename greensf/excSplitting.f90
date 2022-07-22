MODULE m_excSplitting

   use m_types
   use m_constants
   use m_trapz
   use m_xmlOutput
   use m_intgr

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE excSplitting(gfinp,atoms,input,scalarGF,greensfImagPart,ef,vTot)

      TYPE(t_gfinp),             INTENT(IN)  :: gfinp
      TYPE(t_atoms),             INTENT(IN)  :: atoms
      TYPE(t_input),             INTENT(IN)  :: input
      TYPE(t_scalarGF),          INTENT(IN)  :: scalarGF(:)
      TYPE(t_greensfImagPart),   INTENT(IN)  :: greensfImagPart
      REAL,                      INTENT(IN)  :: ef
      type(t_potden),            intent(in)  :: vTot

      INTEGER :: i_gf,i_elem,i_elemLO,ispin,m,l,atomType,nLO,iLO,iLOp
      LOGICAL :: l_sphavg
      REAL    :: excSplit,del,delta
      REAL, ALLOCATABLE :: eMesh(:)
      COMPLEX, ALLOCATABLE :: imag(:,:,:)
      REAL, ALLOCATABLE :: intCOM(:,:), intNorm(:,:)
      REAL, ALLOCATABLE :: bxc(:)
      CHARACTER(LEN=20) :: attributes(4)
      CALL openXMLElementNoAttributes('onSiteExchangeSplitting')

      !Calculate the onsite exchange splitting from Bxc
      allocate(bxc(atoms%jmtd), source=0.0)
      do atomType = 1, atoms%ntype
         !Get the Bxc part of the potential
         !L=0 of potential has an additional rescaling of r/sqrt(4pi)
         !The sqrt(4pi) is the part of the integral over the angular part
         bxc(:atoms%jri(atomType)) = (vTot%mt(:atoms%jri(atomType),0,atomType,1) - vTot%mt(:atoms%jri(atomType),0,atomType,2))/2.0 &
                                    * sfp_const*atoms%rmsh(:atoms%jri(atomType),atomType)
         
         CALL intgr3(bxc,atoms%rmsh(:,atomType),atoms%dx(atomType),atoms%jri(atomType),delta)

         attributes = ''
         WRITE(attributes(1),'(i0)') atomType
         WRITE(attributes(2),'(f12.7)') delta * hartree_to_ev_const
         WRITE(attributes(3),'(a2)') 'eV'
         CALL writeXMLElementForm('bxcIntegral',['atomType','Delta   ','units   '],&
                                  attributes(:3),reshape([8,1,5,5,6,1,12,2],[3,2]))

      enddo

      IF(.NOT.gfinp%checkOnsite()) RETURN

      CALL gfinp%eMesh(ef,del=del,eMesh=eMesh)

      ALLOCATE(imag(SIZE(eMesh),-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const))
      ALLOCATE(intNorm(SIZE(eMesh),input%jspins),source=0.0)
      ALLOCATE(intCOM(SIZE(eMesh),input%jspins),source=0.0)


      WRITE(oUnit,9000)
9000  FORMAT(/,'Onsite Exchange Splitting (from imag. part of GF)')
      WRITE(oUnit,'(A)') '---------------------------------------------------'

      DO i_gf = 1, gfinp%n

         l  = gfinp%elem(i_gf)%l
         atomType = gfinp%elem(i_gf)%atomType
         l_sphavg = gfinp%elem(i_gf)%l_sphavg
         nLO = gfinp%elem(i_gf)%countLOs(atoms)
         !Only onsite exchange splitting
         IF(gfinp%elem(i_gf)%isOffDiag()) CYCLE
         IF(gfinp%elem(i_gf)%l_kresolved_int) CYCLE
         IF(.NOT.gfinp%isUnique(i_gf, distinct_kresolved_int=.TRUE.)) CYCLE

         i_elem = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg)
         i_elemLO = gfinp%uniqueElements(atoms,max_index=i_gf,l_sphavg=l_sphavg,lo=.TRUE.)

         !-------------------------------------------------
         ! Evaluate center of mass of the bands in question
         ! and take the difference between spin up/down
         ! We use the same cutoffs as for the kk-integration
         !-------------------------------------------------
         ! E_COM = 1/(int dE N_l(E)) * int dE E * N_l(E)
         !-------------------------------------------------
         excSplit = 0.0
         DO ispin = 1, input%jspins
            intCOM = 0.0
            intNorm = 0.0
            IF(l_sphavg) THEN
               imag = greensfImagPart%applyCutoff(i_elem,i_gf,ispin,l_sphavg)
            ELSE
               imag = greensfImagPart%applyCutoff(i_elem,i_gf,ispin,l_sphavg,imat=1)
               imag = imag + greensfImagPart%applyCutoff(i_elem,i_gf,ispin,l_sphavg,imat=2) * scalarGF(i_gf)%ddn(ispin,ispin)
               IF(nLO>0) THEN
                  DO iLO = 1, nLO
                     imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,ispin,l_sphavg,imat=1,iLO=iLO) * scalarGF(i_gf)%uulon(ilo,ispin,ispin)
                     imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,ispin,l_sphavg,imat=2,iLO=iLO) * scalarGF(i_gf)%uloun(ilo,ispin,ispin)
                     imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,ispin,l_sphavg,imat=3,iLO=iLO) * scalarGF(i_gf)%dulon(ilo,ispin,ispin)
                     imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,ispin,l_sphavg,imat=4,iLO=iLO) * scalarGF(i_gf)%ulodn(ilo,ispin,ispin)
                     DO iLOp = 1, nLO
                        imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,ispin,l_sphavg,iLO=iLO,iLOp=iLOp) * scalarGF(i_gf)%uloulopn(ilo,ilop,ispin,ispin)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
            DO m = -l, l
               intCOM(:,ispin) = intCOM(:,ispin) + eMesh*REAL(imag(:,m,m))
               intNorm(:,ispin) = intNorm(:,ispin) + REAL(imag(:,m,m))
            ENDDO
            excSplit = excSplit + (-1)**(ispin) * 1.0/(trapz(intNorm(:,ispin),del,SIZE(eMesh))) &
                                                     * trapz(intCOM(:,ispin),del,SIZE(eMesh))
         ENDDO
         WRITE(oUnit,'(A,I4,A,I4,A,f10.4,A)') '  atom: ', atomType, '   l: ', l,&
                                            '    DeltaExc: ',excSplit * hartree_to_ev_const, ' eV'

         attributes = ''
         WRITE(attributes(1),'(i0)') atomType
         WRITE(attributes(2),'(i0)') l
         WRITE(attributes(3),'(f12.7)') excSplit * hartree_to_ev_const
         WRITE(attributes(4),'(a2)') 'eV'
         CALL writeXMLElementForm('excSplit',['atomType','l       ','Delta   ','units   '],&
                                  attributes,reshape([8,1,5,5,6,1,12,2],[4,2]))

      ENDDO
      WRITE(oUnit,'(/)')
      CALL closeXMLElement('onSiteExchangeSplitting')

   END SUBROUTINE excSplitting

END MODULE m_excSplitting
