MODULE m_tetrahedronInit

   CONTAINS

   SUBROUTINE tetrahedronInit(ikpt,kpts,input,neig,eig,gCoeffs,ef,resWeights,dosWeights,boundInd)

      !This Subroutine sets up the weights for the brillouin zone integration using
      !the tetrahedron method for the greens function calculations
      !
      !If we are in a non-magnetic or collinear case we can decompose the greens function
      !into real and imaginary part and obtain the real part from the imaginary part
      !via Kramers-Kronig transformation. Then we only need the weights for a dos calculations
      !
      !Otherwise we need the weights for the resolvent function

      USE m_types
      USE m_juDFT
      USE m_dosWeights
      USE m_resWeights

      IMPLICIT NONE

      INTEGER,                INTENT(IN)  :: ikpt !Current k-point
      TYPE(t_kpts),           INTENT(IN)  :: kpts
      TYPE(t_input),          INTENT(IN)  :: input
      INTEGER,                INTENT(IN)  :: neig
      REAL,                   INTENT(IN)  :: eig(:,:)
      TYPE(t_greensfCoeffs),  INTENT(IN)  :: gCoeffs
      REAL,                   INTENT(IN)  :: ef

      REAL,                   INTENT(INOUT) :: resWeights(:,:)
      REAL,                   INTENT(INOUT) :: dosWeights(:,:)
      INTEGER,                INTENT(INOUT) :: boundInd(:,:) !Indices in between which the weights are non zero
                                                                !to save computation time

      resWeights = 0.0
      dosWeights = 0.0
      boundInd = 0

      IF(input%film) THEN
         IF(input%l_resolvent) THEN
            CALL juDFT_error("Not yet ready")
            !CALL resWeightsCalc(ikpt,kpts,neig,eig,ez,nz,resWeights,boundInd)
         ENDIF
         CALL dosWeightsCalcTria(ikpt,kpts,neig,eig,gCoeffs,ef,dosWeights,boundInd)
      ELSE
         IF(input%l_resolvent) THEN
            CALL resWeightsCalc(ikpt,kpts,neig,eig,gCoeffs,resWeights,boundInd)
         ENDIF
         CALL dosWeightsCalc(ikpt,kpts,neig,eig,gCoeffs,ef,dosWeights,boundInd)
      ENDIF

      IF(input%l_resolvent) THEN
         boundInd(:,1) = 1
         boundInd(:,2) = gCoeffs%ne
      ENDIF

   END SUBROUTINE tetrahedronInit

END MODULE m_tetrahedronInit