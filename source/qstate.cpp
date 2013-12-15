/*

 */

#include "qstate.h"

QState::QState(int ns, int db) {

	//Set the value of n_sites and hilbert_dim
	n_sites = ns;
	hilbert_dim = db;

	//Derived dull state dimension.
	full_state_dim = hilbert_dim^n_sites;

	//Our aim is to define simple a 1 by hilbert_dim^n_sites vector.
	//By default we simple set it to be the 000... state

	//Create state initially filled with Zeros
	//CanonMat init_state = CanonMat::Zero( full_state_dim, 1);


	CanonMat init_state = CanonMat::Zero( full_state_dim, 1);

	//Set the 00... state to 1.0
	init_state(0,0) = 1.0;

	//Store the state in the shared pointer
	*ptr_state = init_state;

}

QState::~QState(){}

int QState::return_n_sites(){
	return n_sites;
}


int QState::return_hilbert_dim(){
	return hilbert_dim;
}


int QState::return_full_state_dim(){
	return full_state_dim;
}

CanonMat QState::return_qstate() {
	return *ptr_state;
}

