/* This is a sample program which will explore the use of the Eigen
 * package for linear algebra.
 * Specifically of interest are the array and matrix classes
 * and the SVD routines provided.
 */

/* Update: so far added is the use of the Eigen package for
 * non-square matrices and an example of how to use template
 * functions such that arguments and returns can be generic
 * matrix objects.
 */


#include "Eigen/Dense"
#include <iostream>
using namespace Eigen;

#define complex std::complex<double>

/* Function to multiply two matrices in canonical form
 * for a given site with Nsig states and ni x nj submatrices
 */

//void canonmult(<matrixbase> &matrix){
//}

template <typename DerivedA,typename DerivedB>
DerivedA squaredist(const MatrixBase<DerivedA>& p1,const MatrixBase<DerivedB>& p2)
{
  return p1*p2.transpose();
}

int main(void) {

	int Nsigma = 3;		//Quantum dimension of a site
	int Ni = 2;			//Dimension for matrix indices i
	int Nj = 2;			//and j
	int Nsigmai = Ni * Nsigma; //Canonical i dimension
	typedef complex typeword;  //type for all matrices/arrays

	Matrix<typeword, Dynamic, Dynamic> M(Nsigmai, Nj);
	M = Matrix<typeword, Dynamic, Dynamic>::Random(Nsigmai, Nj);
	std::cout << M << std::endl;
	std::cout << std::endl;

	std::cout << squaredist(M,2*M) << std::endl;

	std::cout << std::endl;

	Matrix<typeword, Dynamic, Dynamic> N = M.transpose();
	std::cout << squaredist(N,N) << std::endl;



}
