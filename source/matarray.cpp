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

//HEY MY NAMES MICHAEL !!!

#include "Eigen/Dense"
#include "Eigen/SVD"
#include <iostream>
/* selected definitions from std*/
#define complex std::complex<double>
#define cout std::cout
#define endl std::endl

using namespace Eigen;

/* Function to multiply two matrices in canonical form
 * for a given site with Nsig states and ni x nj submatrices
 */

template<typename DerivedA, typename DerivedB>
DerivedA transmult(const MatrixBase<DerivedA>& p1,
		const MatrixBase<DerivedB>& p2) {
	return p1 * p2.transpose();
}

/* Takes a square matrix SV [from SVD of M = A S V] of side n
 * and post-multiplies it by a canonical matrix Min
 * [i.e. a column of Nsigma n x m matrices] to give a
 * new canonical matrix Mout of the same dimensions as Min
 * <see Eq.(136) of Schollwoeck>
 */
template<typename DerivedA, typename DerivedB>
DerivedB canonmult(const MatrixBase<DerivedA>& SV,
		const MatrixBase<DerivedB>& Min, int Nsigma) {

	DerivedB Mout(Min.rows(), Min.cols());
	//Mout = DerivedB::Zero(Min.rows(), Min.cols());

	/*cout << endl;
	 cout << "Min is" <<endl;
	 cout << Min << endl;
	 cout << endl;
	 cout << "Mout is" <<endl;
	 cout << Mout << endl;
	 cout << endl;
*/

	for (int k = 0; k < Nsigma; k++) {
		Mout.block(SV.cols() * k, 0, SV.cols(), Min.cols()) = SV
				* Min.block(SV.cols() * k, 0, SV.cols(), Min.cols());
	}

	return Mout;
}



int main(void) {

	int Nsigma = 3;		//Quantum dimension of a site
	int Ni = 3;			//Dimension for matrix indices i
	int Nj = 3;			//and j
	int Nsigmai = Ni * Nsigma; //Canonical i dimension
	typedef complex typeword;  //type for all matrices/arrays


	Matrix<typeword, Dynamic, Dynamic> M(Nsigmai, Nj), N(Nsigmai, Nj);
	for (int k = 0; k < Nsigma; k++) {
		M.block(Ni * k, 0, Ni, Nj) = (k + 1)
				* Matrix<typeword, Dynamic, Dynamic>::Ones(Ni, Nj);
	}
	Matrix<typeword, Dynamic, Dynamic> P(Ni, Ni);
	P = 2 * Matrix<typeword, Dynamic, Dynamic>::Identity(Ni, Ni);

	cout<< "This is an example of an M matrix (stack of Nsigma ="<<Nsigma<< "matrices)"<<endl;
	cout<< M << endl;
	cout<< endl;

	cout<< "This is a sample SV to be multiplied in canonical form"<<endl;
	cout<< P << endl;
	cout<< endl;

	cout<< "We call a canonical matrix multiplication without copying:" << endl;
	N = canonmult(P, M, Nsigma);
	cout<< N << endl<< endl;

	cout<< "We can easily to SVD too:" << endl;
	JacobiSVD<Matrix<typeword, Dynamic, Dynamic> > svd(M, ComputeThinU | ComputeThinV);
	//cout << "The singular values are" << endl << svd.singularValues() << endl << endl;
	cout << "This U matrix will be pointed to by eg. A^{sigma=1}" << endl << svd.matrixU() << endl<< endl;
	cout << "This V will go on to make the next M^{sigma=2}" << endl << svd.matrixV() << endl;





}
