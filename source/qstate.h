/*
 *
 */

#ifndef QSTATE_H_
#define QSTATE_H_


#include "matrixdefs.h"



/*

*/


class QState {

	int n_sites, hilbert_dim, full_state_dim;

	CanonVec_ptr ptr_state;

public:

	/* Constructor
	 * Given a int for the number of sites
	 * and an int for the Hilbert space dimension for each site
	 * produce a state in the form 000...
	 * */
	QState(int, int);

	/* Destructor
	 *  Destructor NOT IMPLEMENTED!!!
	 */
	~QState();


	int return_n_sites();

	int return_hilbert_dim();

	int return_full_state_dim();

	CanonVec return_qstate();

	/*
	 * Check data is valid and consistent
	 * */
	void validate_QState();

};

#endif /* QSTATE_H_ */
