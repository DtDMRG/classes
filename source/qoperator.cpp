/*
 * qoperator.cpp
 *
 *  Created on: Dec 20, 2013
 *      Author: michaelhush
 */

#include "qoperator.h"

QOperator::QOperator(unsigned td) {

	basic_contructor_code(td);

}

QOperator::QOperator(unsigned ns, unsigned db) {

	unsigned total_dim = pow(db,ns);

	basic_contructor_code(total_dim);

}

QOperator::~QOperator(){}

void QOperator::basic_contructor_code(unsigned td) {

	//Store matrix dimensions, which will be of size td td
	matrix_dim = td;

	//Create Operator initially filled with Zeros
	operator_matrix = CanonMat::Identity(matrix_dim,matrix_dim);

}

int QOperator::return_dim(){
	return matrix_dim;
}

CanonMat QOperator::return_qoperator_matrix() {
	return operator_matrix;
}


