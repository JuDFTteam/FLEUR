/* -*- Mode: C++; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
// This file is a part of ChASE.
// Copyright (c) 2015-2018, Simulation Laboratory Quantum Materials,
//   Forschungszentrum Juelich GmbH, Germany
// and
// Copyright (c) 2016-2018, Aachen Institute for Advanced Study in Computational
//   Engineering Science, RWTH Aachen University, Germany All rights reserved.
// License is 3-clause BSD:
// https://github.com/SimLabQuantumMaterials/ChASE/

#include <complex.h>
#include <mpi.h>
#include <random>

#include "ChASE-MPI/chase_mpi.hpp"
#include "ChASE-MPI/chase_mpi_properties.hpp"
#include "ChASE-MPI/impl/chase_mpihemm_blas.hpp"
#include "ChASE-MPI/impl/chase_mpihemm_blas_seq_inplace.hpp"

using namespace chase;
using namespace chase::mpi;

class ChASE_State {
 public:
  template <typename T>
  static ChaseMpiProperties<T>* constructProperties(std::size_t N, std::size_t nev,
                                             std::size_t nex, MPI_Comm comm);
  template <typename T>
  static ChaseMpiProperties<T>* getProperties();

  static ChaseMpiProperties<double>* double_prec;
  static ChaseMpiProperties<std::complex<double>>* complex_double_prec;
};

ChaseMpiProperties<double>* ChASE_State::double_prec = nullptr;
ChaseMpiProperties<std::complex<double>>* ChASE_State::complex_double_prec =
    nullptr;

template <>
 ChaseMpiProperties<double>* ChASE_State::constructProperties(std::size_t N,
                                                             std::size_t nev,
                                                             std::size_t nex,
                                                             MPI_Comm comm) {
  double_prec = new ChaseMpiProperties<double>(N, nev, nex, comm);
  return double_prec;
}

template <>
ChaseMpiProperties<std::complex<double>>* ChASE_State::constructProperties(
    std::size_t N, std::size_t nev, std::size_t nex, MPI_Comm comm) {
  complex_double_prec =
      new ChaseMpiProperties<std::complex<double>>(N, nev, nex, comm);
  return complex_double_prec;
}

template <>
ChaseMpiProperties<double>* ChASE_State::getProperties() {
  return double_prec;
}
template <>
ChaseMpiProperties<std::complex<double>>* ChASE_State::getProperties() {
  return complex_double_prec;
}

template <typename T>
void call_chase(T* H, int* N, T* V, Base<T>* ritzv, int* nev, int* nex,
                int* deg, double* tol, char* mode, char* opt) {
  typedef ChaseMpi<ChaseMpiHemmBlasSeqInplace, T> SEQ_CHASE;
  std::cerr << "entering chase" << std::endl;

  std::mt19937 gen(2342.0);
  std::normal_distribution<> d;

  SEQ_CHASE single(*N, *nev, *nex, V, ritzv, H);

  ChaseConfig<T>& config = single.GetConfig();
  config.SetTol(*tol);
  config.SetDeg(*deg);
  config.SetOpt(*opt == 'S');
  config.SetApprox(*mode == 'A');

  if (!config.UseApprox())
    for (std::size_t k = 0; k < *N * (*nev + *nex); ++k)
      V[k] = getRandomT<T>([&]() { return d(gen); });

  chase::Solve(&single);
}

template <typename T>
void chase_setup(MPI_Fint* fcomm, int* N, int* nev, int* nex, std::size_t* xoff,
                 std::size_t* yoff, std::size_t* xlen, std::size_t* ylen) {
  MPI_Comm comm = MPI_Comm_f2c(*fcomm);
  auto props = ChASE_State::constructProperties<T>(*N, *nev, *nex, comm);
  props->get_off(xoff, yoff, xlen, ylen);
}

template <typename T>
void chase_solve(T* H, T* V, Base<T>* ritzv, int* deg, double* tol, char* mode,
                 char* opt) {
  typedef ChaseMpi<ChaseMpiHemmBlas, T> CHASE;
  std::cerr << "entering chase" << std::endl;

  std::mt19937 gen(2342.0);
  std::normal_distribution<> d;
  ChaseMpiProperties<T>* props = ChASE_State::getProperties<T>();

  CHASE single(props, V, ritzv);

  T* H_ = single.GetMatrixPtr();
  std::size_t xlen, ylen, xoff, yoff;
  props->get_off(&xoff, &yoff, &xlen, &ylen);
  ChaseConfig<T>& config = single.GetConfig();
  auto N = config.GetN();
  auto nev = config.GetNev();
  auto nex = config.GetNex();

  if (!config.UseApprox())
    for (std::size_t k = 0; k < N * (nev + nex); ++k)
      V[k] = getRandomT<T>([&]() { return d(gen); });

  for (std::size_t k = 0; k < xlen * ylen; ++k) H_[k] = H[k];

  config.SetTol(*tol);
  config.SetDeg(*deg);
  config.SetOpt(*opt == 'S');
  config.SetApprox(*mode == 'A');

  chase::Solve(&single);
}

extern "C" {

void zchase_(std::complex<double>* H, int* N, std::complex<double>* V,
             double* ritzv, int* nev, int* nex, int* deg, double* tol,
             char* mode, char* opt) {
  call_chase<std::complex<double>>(H, N, V, ritzv, nev, nex, deg, tol, mode,
                                   opt);
}

void dchase_(double* H, int* N, double* V, double* ritzv, int* nev, int* nex,
             int* deg, double* tol, char* mode, char* opt) {
  call_chase<double>(H, N, V, ritzv, nev, nex, deg, tol, mode, opt);
}

void zchase_init(MPI_Fint* fcomm, int* N, int* nev, int* nex, int* xoff,
                 int* yoff, int* xlen, int* ylen) {
  std::size_t xoff_, yoff_, xlen_, ylen_;

  chase_setup<std::complex<double>>(fcomm, N, nev, nex, &xoff_, &yoff_, &xlen_,
                                    &ylen_);

  *xoff = static_cast<int>(xoff_);
  *yoff = static_cast<int>(yoff_);
  *xlen = static_cast<int>(xlen_);
  *ylen = static_cast<int>(xlen_);
}

void dchase_init(MPI_Fint* fcomm, int* N, int* nev, int* nex, int* xoff,
                 int* yoff, int* xlen, int* ylen) {
  std::size_t xoff_, yoff_, xlen_, ylen_;

  chase_setup<double>(fcomm, N, nev, nex, &xoff_, &yoff_, &xlen_, &ylen_);

  *xoff = static_cast<int>(xoff_);
  *yoff = static_cast<int>(yoff_);
  *xlen = static_cast<int>(xlen_);
  *ylen = static_cast<int>(xlen_);
}

void zchase_solve(std::complex<double>* H, std::complex<double>* V,
                  double* ritzv, int* deg, double* tol, char* mode, char* opt) {
  chase_solve<std::complex<double>>(H, V, ritzv, deg, tol, mode, opt);
}

void dchase_solve(double* H, double* V, double* ritzv, int* deg, double* tol,
                  char* mode, char* opt) {
  chase_solve<double>(H, V, ritzv, deg, tol, mode, opt);
}

}  // extern C
