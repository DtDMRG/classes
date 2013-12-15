/*
 * matrixdefs.h
 *
 *  Created on: 9 Dec 2013
 *      Author: genway
 */

#ifndef MATRIXDEFS_H_
#define MATRIXDEFS_H_

//Define the index type if needed
//#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int

#include "Eigen/Dense"
#include "Eigen/SVD"
#include <memory>
#include <vector>

#define shared_ptr std::shared_ptr
#define vector std::vector

using namespace Eigen;

#define complex std::complex<double>

typedef complex datatype;
typedef Matrix<datatype, Dynamic, Dynamic> CanonMat;
typedef shared_ptr<CanonMat> CanonMat_ptr;
typedef vector<CanonMat_ptr>::iterator CanonMat_itr;

class Svd: public JacobiSVD<CanonMat>{
public:
	Svd(CanonMat M): JacobiSVD<CanonMat>(M, ComputeThinU | ComputeThinV){}
	~Svd(){}
};







#endif /* MATRIXDEFS_H_ */
