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
using namespace Eigen;

#define complex std::complex<double>
typedef complex datatype;

typedef Matrix<datatype, Dynamic, Dynamic> CanonMat;


class Svd: public JacobiSVD<CanonMat>{
public:
	Svd(CanonMat M): JacobiSVD<CanonMat>(M, ComputeThinU | ComputeThinV){}
	~Svd(){}
};







#endif /* MATRIXDEFS_H_ */
