MODULE m_excSplitting

   USE m_types
   USE m_constants
   USE m_trapz
   USE m_xmlOutput

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE excSplitting(gfinp,atoms,input,usdus,greensfImagPart,ef)

      TYPE(t_gfinp),             INTENT(IN)  :: gfinp
      TYPE(t_atoms),             INTENT(IN)  :: atoms
      TYPE(t_input),             INTENT(IN)  :: input
      TYPE(t_usdus),             INTENT(IN)  :: usdus
      TYPE(t_greensfImagPart),   INTENT(IN)  :: greensfImagPart
      REAL,                      INTENT(IN)  :: ef

      INTEGER :: i_gf,i_elem,i_elemLO,indUnique,ispin,m,l,lp,atomType,atomTypep,nLO,iLO,iLOp
      LOGICAL :: l_sphavg
      REAL    :: excSplit,del,atomDiff(3)
      REAL, ALLOCATABLE :: eMesh(:), imag(:)
      REAL, ALLOCATABLE :: intCOM(:,:), intNorm(:,:)
      CHARACTER(LEN=20) :: attributes(4)

      IF(.NOT.gfinp%checkOnsite()) RETURN

      CALL gfinp%eMesh(ef,del=del,eMesh=eMesh)

      ALLOCATE(intNorm(SIZE(eMesh),input%jspins),source=0.0)
      ALLOCATE(intCOM(SIZE(eMesh),input%jspins),source=0.0)


      WRITE(oUnit,9000)
9000  FORMAT(/,'Onsite Exchange Splitting (from imag. part of GF)')
      WRITE(oUnit,'(A)') '---------------------------------------------------'

      CALL openXMLElementNoAttributes('onSiteExchangeSplitting')

      DO i_gf = 1, gfinp%n

         l  = gfinp%elem(i_gf)%l
         lp = gfinp%elem(i_gf)%lp
         atomType = gfinp%elem(i_gf)%atomType
         atomTypep = gfinp%elem(i_gf)%atomTypep
         l_sphavg = gfinp%elem(i_gf)%l_sphavg
         atomDiff = gfinp%elem(i_gf)%atomDiff
         nLO = gfinp%elem(i_gf)%countLOs(atoms)
         !Only onsite exchange splitting
         IF(l /= lp) CYCLE
         IF(atomType /= atomTypep) CYCLE
         IF(ANY(ABS(atomDiff).GT.1e-12)) CYCLE

         i_elem = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,indUnique=indUnique)
         i_elemLO = gfinp%uniqueElements(atoms,ind=i_gf,l_sphavg=l_sphavg,lo=.TRUE.)

         IF(i_gf /= indUnique) CYCLE
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
            DO m = -l, l
               IF(l_sphavg) THEN
                  imag = greensfImagPart%applyCutoff(i_elem,i_gf,m,m,ispin,l_sphavg)
               ELSE
                  imag = greensfImagPart%applyCutoff(i_elem,i_gf,m,m,ispin,l_sphavg,imat=1)
                  imag = imag + greensfImagPart%applyCutoff(i_elem,i_gf,m,m,ispin,l_sphavg,imat=2) * usdus%ddn(l,atomType,ispin)
                  IF(nLO>0) THEN
                     DO iLO = 1, nLO
                        imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,m,m,ispin,l_sphavg,imat=1,iLO=iLO) * usdus%uulon(ilo,atomType,ispin)
                        imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,m,m,ispin,l_sphavg,imat=2,iLO=iLO) * usdus%uulon(ilo,atomType,ispin)
                        imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,m,m,ispin,l_sphavg,imat=3,iLO=iLO) * usdus%dulon(ilo,atomType,ispin)
                        imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,m,m,ispin,l_sphavg,imat=4,iLO=iLO) * usdus%dulon(ilo,atomType,ispin)
                        DO iLOp = 1, nLO
                           imag = imag + greensfImagPart%applyCutoff(i_elemLO,i_gf,m,m,ispin,l_sphavg,iLO=iLO,iLOp=iLOp) * usdus%uloulopn(ilo,ilop,atomType,ispin)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
               intCOM(:,ispin) = intCOM(:,ispin) + eMesh*imag
               intNorm(:,ispin) = intNorm(:,ispin) + imag
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
         CALL writeXMLElementForm('excSplit',['atomType','l       ','Delta   ','unit    '],&
                                  attributes,reshape([8,1,5,4,6,1,12,2],[4,2]))

      ENDDO
      WRITE(oUnit,'(/)')
      CALL closeXMLElement('onSiteExchangeSplitting')

   END SUBROUTINE excSplitting

END MODULE m_excSplitting