#include "kplib/src/kPointLattice.h"
#include "kplib/src/kPointLatticeGenerator.h"
static KPointLattice kp;

extern "C"{
int kplib_generate(const double primVectorsArray[3][3],
  const int latticePointOperatorsArray[][3][3],
  const int numOperators,
  const bool isConventionalHexagonal,
  const bool Gamma,
  const double minDistance,
  const int minSize
){
  KPointLatticeGenerator generator=KPointLatticeGenerator(primVectorsArray,primVectorsArray,latticePointOperatorsArray,numOperators,isConventionalHexagonal);
  if (Gamma){
    generator.includeGamma(TRUE);
  }else{
    generator.includeGamma(FALSE);
  }
  kp=generator.getKPointLattice(minDistance,minSize);
  return kp.numTotalKPoints();
}


void kplib_setkpoints(double kpoints[][3],double weights[]){
  Tensor<double> coords = kp.getKPointCoordinates();
  std::vector<int> w=kp.getKPointWeights();
  for (int i = 0; i < (int) coords.size(); ++i) {
    kpoints[i][0]=coords[i][0];
    kpoints[i][1]=coords[i][1];
    kpoints[i][2]=coords[i][2];
    weights[i]=w[i];
  }
}
}
