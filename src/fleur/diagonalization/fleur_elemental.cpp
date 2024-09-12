/*
   Copyright (c) 2014, Daniel Wortmann
   All rights reserved.

   This file provides an interface from FLEUR to Elemental 
*/
#include "elemental.hpp"
using namespace std;
using namespace elem;

// Typedef our real or complex types to 'C' for convenience
#ifdef CPP_INVERSION
typedef double C;
#else
typedef Complex<double> C;
#endif

class global_data {
 public:
  Grid *g;
  mpi::Comm mpi_comm;
  int matrix_dimension;
  DistMatrix<C> *H_mat,*S_mat;
  DistMatrix<C,STAR,VC> eigenvectors;
  DistMatrix<double,VC,STAR> eigenvalues;

  
};

//Global variables
  

static global_data *gd;


DistMatrix<C>* fleur_matrix(int n,C* buffer)
{
  // Create the distributed matrix
  DistMatrix<C,STAR,VC> *mat;
  // The Matrix should be a n x n matrix in 1-D cyclic distribution
  mat= new DistMatrix<C,STAR,VC>(n,n,*(gd->g));
  C* localbuffer=mat->Buffer(); // this is the local buffer of the matrix
  const int buffersize1=mat->LocalWidth();
  const int buffersize2=mat->LocalHeight();
  int localindex=0;
  int fleur_index=0; 
  // Now copy all the data into local buffer to initialize the matrix
  // Loop over columns of local data
  for (int i=mpi::CommRank(gd->mpi_comm);i<n;i+=mpi::CommSize(gd->mpi_comm))
    {
      for (int j=0;j<=i;j++)
	{
	  localbuffer[localindex++]=buffer[fleur_index++];
	}
      // further Off-diagonal elements are set to zero !Probably not needed
      for (int j=0;j<n-i-1;j++)
	{
	  localbuffer[localindex++]=0.0;
	}
    }
  DistMatrix<C> *mat2=new DistMatrix<C>(*mat);
  delete mat;
  //*mat2=*mat;
  return mat2;
}
extern "C"{ 
void fl_el_initialize(int n, C* hbuf, C* sbuf, int mpi_used_comm)
// Set the two matrices
{
  // Initialize the Library
  int argc=0; char** argv;
  Initialize( argc, argv );
  
  //Store the matrix dimension& the mpi_communicator
  gd = new global_data;
  gd->mpi_comm=MPI_Comm_f2c(mpi_used_comm);
  gd->matrix_dimension=n;
  
  
  // First we need a mpi-grid
  gd->g= new Grid(gd->mpi_comm);
  // Store the Matrices
  gd->H_mat=fleur_matrix(n,hbuf);
  gd->S_mat=fleur_matrix(n,sbuf);
  
}

void fl_el_diagonalize(int no_of_eigenpairs)
// Diagonalize the Matrix and return the number of local eigenvalues
{
  
  /* this is for the development version
  // The subset determines the no of eigenvalues found
  HermitianEigSubset<double> subset;
  subset.indexSubset=true;
  subset.lowerIndex=0;
  subset.upperIndex=no_of_global_eigenpairs;
  
  // Space for eigenvalues
  DistMatrix<double> eigenval(g);
  DistMatrix<C> eigenvec(g);
  //default sorting
  const SortType sort = static_cast<SortType>(0);

  //call diagonalization
  HermitianGenDefEig( AXBX, LOWER, H_mat, S_mat, eigenval, eigenvec, sort, subset );
  */
  DistMatrix<double, VR, STAR> eigenval(*(gd->g));
  DistMatrix<C> evec(*(gd->g));
  HermitianGenDefiniteEigType eigtype=AXBX;
  UpperOrLower uplo=UPPER;
  if (mpi::CommRank(gd->mpi_comm)==0) {
    cout<<"H/S-matrix of size "<<gd->matrix_dimension<<endl;
  }
  Display(*(gd->H_mat));
  Display(*(gd->S_mat));
  HermitianGenDefiniteEig(eigtype, uplo, *(gd->H_mat), *(gd->S_mat), eigenval, evec,0,no_of_eigenpairs);
  //redistribute matrices
  //eigenvalues are of type DistMatrix(C,STAR,VC);
  gd->eigenvalues=eigenval;
  gd->eigenvectors=evec;

  no_of_eigenpairs=gd->eigenvectors.LocalWidth();
    }
  

void fl_el_eigenvalues(int neig, double* eig){
  //Return the eigenvalues

  double* buf=gd->eigenvalues.Buffer();
  if (neig > gd->eigenvalues.LocalWidth()*gd->eigenvalues.LocalHeight())
    {
      cerr<<"Error in dimensions in fleur_elemental\n";
    }
  for (int i=0; i<neig;i++){
    eig[i]=buf[i];
  }
}

void fl_el_eigenvectors(int neig, double* eig, C* eigvec){
  //Return all the local eigenvectors&eigenvalues
  double* eigbuf=gd->eigenvalues.Buffer();
  Display(gd->eigenvalues);
  Display(gd->eigenvectors);
  C* eigbuff=gd->eigenvectors.Buffer();
  int local_index=0;
  for (int i=0; i<neig;i++)
    {
      //Copy eigenvalue
      int pe=mpi::CommRank(gd->mpi_comm);
      int in=i*mpi::CommSize(gd->mpi_comm)+pe;
      cout<< "PE:"<<pe<<":"<<i<<"->"<<in<<endl;
      eig[i]=eigbuf[i];
      //Copy eigenvector
      for (int j=0;j<gd->matrix_dimension; j++){
	eigvec[local_index]=eigbuff[local_index];
        local_index++;
      }
    }
}
}	
  

