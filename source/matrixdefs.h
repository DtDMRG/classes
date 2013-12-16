/*
 * matrixdefs.h
 *
 *  Created on: 9 Dec 2013
 *      Author: genway
 */

#ifndef MATRIXDEFS_H_
#define MATRIXDEFS_H_

#include "Eigen/Dense"
#include "Eigen/SVD"

//Need to check that we are using c++11 if not, we have to get shared pointers from an older reource
#if __cplusplus<201103L
	#include <memory>
	#define shared_ptr std::shared_ptr
#else
	#include <tr1/memory>
	#define shared_ptr std::tr1::shared_ptr
#endif

#include <vector>

#define vector std::vector

using namespace Eigen;

#define complex std::complex<double>

//Type all matricies are set to by default
typedef complex datatype;

//Definition of the canonical matrix.
typedef Matrix<datatype, Dynamic, Dynamic> CanonMat;

//A diagonal version os the canonical matrix
typedef DiagonalMatrix<datatype, Dynamic, Dynamic> DiagonalCanonMat;

//Vector used to store all values
typedef Matrix<datatype, Dynamic, 1> CanonVec;

typedef shared_ptr<CanonMat> CanonMat_ptr;
typedef shared_ptr<CanonVec> CanonVec_ptr;

typedef vector<CanonMat_ptr>::iterator CanonMat_itr;

typedef JacobiSVD<CanonMat> Svd;


#endif /* MATRIXDEFS_H_ */
