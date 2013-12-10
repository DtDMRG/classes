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


#define complex std::complex<double>
typedef complex datatype;
typedef Matrix<datatype, Dynamic, Dynamic> CanonMat;

class Svd {
public:
	CanonMat U,V, Sigma;
	Svd(CanonMat* M){
		Svd=JacobiSVD<CanonMat>(M, ComputeThinU | ComputeThinV);
		U=svd.matrixU();
		V=svd.matrixV();
		int nz=svd.nonzeroSingularValues();
blabla;
		CanonMat S(nz,nz);
		//for (int i=0;i<nz;i++){S(i,1)=}
	}
	virtual ~Svd();
};


#endif /* MATRIXDEFS_H_ */
