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

#include <iostream>
#include "matrixdefs.h"
#include "mps.h"
#include "qstate.h"
#include "hushdebug.h"

/* selected definitions from std*/
#define complex std::complex<double>
#define cout std::cout
#define endl std::endl

using namespace Eigen;



/* Function to multiply two matrices in canonical form
 * for a given site with Nsig states and ni x nj submatrices
 */
/*template<typename DerivedA, typename DerivedB>
DerivedA transmult(const MatrixBase<DerivedA>& p1,
		const MatrixBase<DerivedB>& p2) {
	return p1 * p2.transpose();
}
*/

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

	for (int k = 0; k < Nsigma; k++) {
		Mout.block(SV.cols() * k, 0, SV.cols(), Min.cols()) = SV
				* Min.block(SV.cols() * k, 0, SV.cols(), Min.cols());
	}

	return Mout;
}

/* A mockup of the *TYPE* of thing we need to do taking canonical
 * matrices and then SVDing and re-pointing
 */
/*
template<typename Derived>
void svdmult(const MatrixBase<Derived>& M1, const MatrixBase<Derived>& M2) {

	JacobiSVD<Derived> svd(M1, ComputeThinU | ComputeThinV);

	const_cast<MatrixBase<Derived>&> (M1) = svd.matrixU();
	const_cast<MatrixBase<Derived>&> (M2) = svd.matrixV();

}

typedef complex datatype;*/


/*
typedef Matrix<datatype, Dynamic, Dynamic> CanonMat;
class Svd: public JacobiSVD<CanonMat>{
public:
	Svd(CanonMat M): JacobiSVD<CanonMat>(M, ComputeThinU | ComputeThinV){}
	~Svd(){}
};
*/

int main( int argc, char* argv[] ) {

	/*Split the code into parts so we each have our own effective mains
	 * where we can run test code relevant to what we are developing
	 * Feel free to add whatever you want
	 * Just feed a the appropriate character to run to run your code
	 *
	 */

	if (argc > 1) {

		if (strcmp(argv[1],"h") == 0) {
			HushDebug::hushmain();

			return 0;
		}

		if (strcmp(argv[1],"g") == 0) {
			//Sam debug

			return 0;
		}

		if (strcmp(argv[1],"l") == 0) {
			//Levi debug

			return 0;
		}


	}


	//cout<< "Hey here is number of the parameters you've entered" <<endl;
	//cout<< argc <<endl;
	//cout<< "and what they said" <<endl;
	//cout<< argv[0] <<endl;
	//cout<< argv[1] <<endl;



	int Nsigma = 3;		//Quantum dimension of a site
	int Ni = 3;			//Dimension for matrix indices i
	int Nj = 3;			//and j
	int Nsigmai = Ni * Nsigma; //Canonical i dimension
	typedef complex typeword;  //type for all matrices/arrays

	//Define and allocate array of pointers to type Matrix
	//NOTE const is the Matrix not the pointer!!!
	const Matrix<typeword, Dynamic, Dynamic> **objarray;
	objarray = new const Matrix<typeword, Dynamic, Dynamic>*[2];

	//Allocate some matrices to play with
	Matrix<typeword, Dynamic, Dynamic> M(Nsigmai, Nj), N(Nsigmai, Nj);
	for (int k = 0; k < Nsigma; k++) {
		M.block(Ni * k, 0, Ni, Nj) = (k + 1)
				* Matrix<typeword, Dynamic, Dynamic>::Ones(Nj,Nj);
	}
	Matrix<typeword, Dynamic, Dynamic> P(Ni, Ni);
	P = 2 * Matrix<typeword, Dynamic, Dynamic>::Identity(Nj,Nj);


	//Point array elements to these matrix types now
	objarray[0] = &M;
	objarray[1] = &P;

	//Now we'll deal with only the pointers as far as possible

	cout<< "This is an example of an M matrix (stack of Nsigma ="<<Nsigma << " matrices)"<<endl;
	cout<< "It is the deference of the first element of an array of pointers"<<endl;
	cout<< "*objarray[0] ="<<endl;
	cout<< *objarray[0] << endl;
	cout<< endl;

	cout<< "This is a sample SV to be multiplied in canonical form"<<endl;
	cout<< "*objarray[1] ="<<endl;
	cout<< *objarray[1] << endl;
	cout<< endl;

	cout<< "We call a canonical matrix multiplication without copying:" << endl;
	N = canonmult(*objarray[1], *objarray[0], Nsigma);
	cout<< N << endl;


	cout<< "We can easily to SVD too:" << endl;

	//We will need to do this via function calls and the templates are
	//the way to do this in Eigen:
	/*cout<< "This can be done by a function call - *take note* of the syntax"<<endl;
	cout<< "Now calling svdmult:"<<endl<<endl;
	svdmult(*objarray[0],*objarray[1]);
	cout << "This U matrix IS pointed to by objarray[0]" << endl;
	cout << "*objarray[0] =" << endl;
	cout<< *objarray[0]<< endl<< endl;
	cout << "This V matrix IS pointed to by objarray[1]" << endl;
	cout << "*objarray[1] =" << endl;
	cout<< *objarray[1]<< endl<<endl;*/
	JacobiSVD<CanonMat> KK(M);
	//cout<<KK.singularValues()<<endl;
	//cout<<KK.matrixU()<<endl;
	//cout<<KK.matrixV()<<endl;

cout << "making an mps object"<<endl;
Mps pp(3,3,'i');

cout<<"mps before"<<endl;
cout<<pp.return_matrix_at_site(0)<<endl<<endl;
cout<<pp.return_matrix_at_site(1)<<endl<<endl;
cout<<pp.return_matrix_at_site(2)<<endl<<endl;

pp.change_mps_storage_to_right();
cout<<"left to right storage"<<endl;


cout<<"mps right storage"<<endl;
cout<<pp.return_matrix_at_site(0)<<endl<<endl;
cout<<pp.return_matrix_at_site(1)<<endl<<endl;
cout<<pp.return_matrix_at_site(2)<<endl<<endl;

//pp.sweep_from_right_at(1);
pp.right_compress(1);

cout<<"mps after right canonicalising"<<endl;
cout<<pp.return_matrix_at_site(0)<<endl<<endl;
cout<<pp.return_matrix_at_site(1)<<endl<<endl;
cout<<pp.return_matrix_at_site(2)<<endl<<endl;



exit(0);


//cout<<pp.stored_matrix_dimensions[0]<<endl;
//cout<<pp.stored_matrix_dimensions[1]<<endl;
//cout<<pp.stored_matrix_dimensions[2]<<endl;
//cout<<pp.stored_matrix_dimensions[3]<<endl;
//cout<<pp.stored_matrix_dimensions[4]<<endl;
//cout<<pp.stored_matrix_dimensions[5]<<endl<<endl;
//const CanonMat* MM=pp.mps_pointers[0];
//cout<<*MM<<endl<<endl;
//pp.sweep_from_left_at(0, pp.mps_pointers);
//cout<<pp.mps_pointers[0]<<endl;

//cout<<*(pp.mps_pointers[0])<<endl;
//cout<<*MM<<endl;


//Trying to get the resizing of matricies to work/understand it.
//----------------------------------------------------------

CanonMat original_matrix;

original_matrix.setRandom(3,3);
cout<< "---------------------" << endl;
cout<< "Random 3 by 3 matrix:" << endl;
cout<<original_matrix<<endl;

CanonMat copied_matrix;

copied_matrix = original_matrix;

copied_matrix = copied_matrix*copied_matrix;

cout<< "Original matrix after copy:" << endl;
cout<<original_matrix<<endl;
cout<< "Copied matrix after being multipled by 100:" << endl;
cout<<copied_matrix<<endl;

original_matrix.resize(9,1);

cout<< "Original matrix after resize to vector:" << endl;
cout<<original_matrix<<endl;


CanonVec copied_vector;

copied_vector = original_matrix;

cout<< "Attempt to copy the matrix onto a vector:" << endl;
cout<<copied_vector<<endl;

//Testing out QState

int n_states = 4;
QState simple_state(n_states,2);

cout<< "Print out q state that is stored" << endl;
cout<< simple_state.return_qstate() <<endl;

Mps simple_mps(simple_state);

cout<< "Print out mps state that was made from qstate" << endl;
for (int i=0; i<n_states; i++) {
	cout<< simple_mps.return_matrix_at_site(i) <<endl;
}

CanonMat matrix_1, matrix_2, store_matrix;

matrix_1.setRandom(3,3);
matrix_2.setRandom(2,2);

store_matrix = matrix_1;

cout<< "---------------------" << endl;
cout<< "Stored matrix before:" << endl;
cout<<store_matrix<<endl;

store_matrix = matrix_2;

cout<< "Stored matrix after:" << endl;
cout<<store_matrix<<endl;







}
