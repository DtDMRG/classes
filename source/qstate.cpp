/*

 */

#include "qstate.h"

QState::QState(unsigned ns, unsigned db) {

	basic_contructor_code(ns,db);

	//Set the 00... state to 1.0
	state_vector(0) = 1.0;

}

QState::QState(unsigned ns, unsigned db, unsigned state) {

	basic_contructor_code(ns,db);

	//Set the store... state to 1.0
	state_vector(state) = 1.0;

}

QState::QState(unsigned ns, unsigned db, char state_type) {

	//Special states

	//Currently only one option 'r' which gives random state

	n_sites = ns;
	hilbert_dim = db;

	//Derived dull state dimension.
	full_state_dim = pow(hilbert_dim,n_sites);

	//Our aim is to define simple a 1 by hilbert_dim^n_sites vector.
	//By default we simple set it to be the 000... state

	//Create state initially filled with Zeros
	state_vector = CanonVec::Random( full_state_dim);

}

QState::~QState(){}

void QState::basic_contructor_code(unsigned ns,unsigned db) {

	//Set the value of n_sites and hilbert_dim
	n_sites = ns;
	hilbert_dim = db;

	//Derived dull state dimension.
	full_state_dim = pow(hilbert_dim,n_sites);

	//Our aim is to define simple a 1 by hilbert_dim^n_sites vector.
	//By default we simple set it to be the 000... state

	//Create state initially filled with Zeros
	state_vector = CanonVec::Zero( full_state_dim);


}



int QState::return_n_sites(){
	return n_sites;
}


int QState::return_hilbert_dim(){
	return hilbert_dim;
}


int QState::return_full_state_dim(){
	return full_state_dim;
}

CanonVec QState::return_qstate() {
	return state_vector;
}

