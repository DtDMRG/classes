/*
 *
 */

#ifndef QOPERATOR_H_
#define QOPERATOR_H_


#include "matrixdefs.h"



/*

*/


class QOperator {

	int matrix_dim;

	CanonMat operator_matrix;

	/* Private Functions */

	void basic_contructor_code(unsigned);

public:

	/* Constructors */

	/*
	* Given a int prouduces a idenity matrix state with that specific size
	* */
	QOperator(unsigned);

	/*
	* Given a int for the number of sites
	* and an int for the Hilbert space dimension for each site
	* produces an operator which is the identity
	* */
	QOperator(unsigned, unsigned);

	/* Destructor
	 *  Destructor NOT IMPLEMENTED!!!
	 */
	~QOperator();

	int return_dim();

	CanonMat return_qoperator_matrix();

	/*
	 * Check data is valid and consistent
	 * */
	void validate_QOperator();

};

#endif /* QOperator_H_ */
