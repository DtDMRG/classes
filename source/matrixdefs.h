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

#include <vector>

#define vector std::vector

using namespace Eigen;

#define complex std::complex<double>

//Type all matrices are set to by default
typedef complex datatype;

//Definition of the canonical matrix.
typedef Matrix<datatype, Dynamic, Dynamic> CanonMat;

//A diagonal version of the canonical matrix
typedef DiagonalMatrix<datatype, Dynamic, Dynamic> DiagonalCanonMat;

//Vector used to store all values
typedef Matrix<datatype, Dynamic, 1> CanonVec;

typedef JacobiSVD<CanonMat> Svd;


#endif /* MATRIXDEFS_H_ */
