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
#include <tr1/memory>
#include <vector>
#define shared_ptr std::tr1::shared_ptr
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
