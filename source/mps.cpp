/*
 * mps.cpp
 *
 */

#include "mps.h"

#include <iostream>
#include <stdio.h>      //fputs
#include <stdlib.h>     //abort
#define cout std::cout
#define endl std::endl


/*CONSTRUCTOR
 * Given number of sites and dimension of the quantum basis initializes the array dimensions
 * to all 1, and the array mmps to (dimbasis,1) matrices corresponding to the "lowest" basis
 * element.
 * */

Mps::Mps(unsigned ns, unsigned db) {

	//Set the value of numsites and dimbasis
	n_sites = ns;
	hilbert_dim = db;

	//Initialize memory for pointers and array of dimensions
	constructor_common_memory_init();

	//Define the column matrix of the down state
	CanonMat init= CanonMat( db ,1 );

	stored_matrix_dimensions[0] = 1;

	for (unsigned i=0; i<n_sites;i++) {

		stored_matrix_dimensions[i+1] = 1;

		mps_matrices[i] = CanonMat::Identity( hilbert_dim,1);

	}

}

//To make a random MPS for testing with variable matrix dimensions
Mps::Mps(unsigned ns, unsigned db, char ifpresent) {

	//Set the value of numsites and dimbasis
	n_sites = ns;
	hilbert_dim = db;

	//Initialize memory for pointers and array of dimensions
	constructor_common_memory_init();

	//Define the column matrix of the down state
	stored_matrix_dimensions[0] = 1;

	for (unsigned i=1; i<n_sites/2+1;i++) {

		stored_matrix_dimensions[i] = i+1;

	}
	for (unsigned i=n_sites; i>n_sites/2;i--) {

		stored_matrix_dimensions[i] = n_sites-i+1;

	}
	for (unsigned i=0; i<n_sites; i++){

		mps_matrices[i] = CanonMat::Random( hilbert_dim*stored_matrix_dimensions[i],stored_matrix_dimensions[i+1]);

	}

}

Mps::Mps(QState input_qstate) {

	//Extract and store common variables
	n_sites = input_qstate.return_n_sites();
	hilbert_dim = input_qstate.return_hilbert_dim();

	//Initialize memory for pointers and array of dimensions
	constructor_common_memory_init();


	//Code that actually changes an quantum state (QState) into an MPS (Mps).
	//Implementation of Equation (35)
	//--------------------------------

	//Variables used for the loop

	//Matrix which we use to store what part of the state hasn't yet formed part of the Mps
	CanonMat current_state_residue;
	//To start we copy in the whole state stored passed at the input
	current_state_residue = input_qstate.return_qstate();

	//Store the rank of the first dimension which by definition should be 1
	stored_matrix_dimensions[0] = 1;

	for (unsigned i=0; i<n_sites;i++) {

		unsigned lead_mat_index =stored_matrix_dimensions[i]*hilbert_dim;

		current_state_residue.resize(lead_mat_index, current_state_residue.size()/lead_mat_index);

		Svd temp_svd(current_state_residue, ComputeThinU | ComputeThinV);

		unsigned current_rank = temp_svd.nonzeroSingularValues();

		stored_matrix_dimensions[i+1] = current_rank;

		CanonMat temp_u = temp_svd.matrixU().leftCols(current_rank);

		mps_matrices[i] = temp_u;

		current_state_residue = temp_svd.singularValues().head(current_rank).asDiagonal() * temp_svd.matrixV().adjoint().topRows(current_rank);
	}



}



void Mps::constructor_common_memory_init() {

	    stored_matrix_dimensions= new unsigned[n_sites+1];

	    mps_matrices.resize(n_sites);



}



/*DESTRUCTOR*/
Mps::~Mps(){}



void Mps::sweep_from_left_at(unsigned position) {

	Svd temp_svd(mps_matrices[position], ComputeThinU | ComputeThinV);

	unsigned current_rank = temp_svd.nonzeroSingularValues();

	stored_matrix_dimensions[position+1] = current_rank;

	CanonMat temp_u = temp_svd.matrixU().leftCols(current_rank);

	mps_matrices[position] = temp_u;

	CanonMat temp_new_matrix = leftcanonmult(
			temp_svd.singularValues().head(current_rank).asDiagonal() * temp_svd.matrixV().adjoint().topRows(current_rank),
			temp_new_matrix, current_rank,stored_matrix_dimensions[position+2]);

	mps_matrices[position+1] = temp_new_matrix;
}


//To be made private
//Multiplies a matrix "sv_matrix" by canonical matrix "next_matrix" as hilbert_dim blocks
template<typename DerivedA, typename DerivedB>
DerivedB Mps::leftcanonmult(const MatrixBase<DerivedA>& sv_matrix,
		const MatrixBase<DerivedB>& next_matrix, unsigned rank, unsigned next_matrix_dimension) {

	DerivedB new_next_matrix(next_matrix_dimension*hilbert_dim, rank);

	for (unsigned k = 0; k < hilbert_dim; k++) {
		new_next_matrix.block(rank * k, 0,  next_matrix_dimension, rank) = sv_matrix
				* next_matrix.block(rank * k, 0,  next_matrix_dimension, rank);
	}

	return new_next_matrix;
}



//Must template
CanonMat Mps::return_matrix_at_site(unsigned site) {

	return mps_matrices[site];
}





//Compares the elements of stored_matrix_dimensions with the dimensions of the matrices contained in mps_matrices
int Mps::validate_MPS(){
	for(unsigned int i=0; i<mps_matrices.size(); i++){
		unsigned left =stored_matrix_dimensions[i]*hilbert_dim;


		if( left !=mps_matrices[i].rows() || stored_matrix_dimensions[i+1]!=mps_matrices[i].cols() ){
			fputs ("Dimensions of matrices do not match.\n",stderr);
			abort();
		}
//		else	cout<<mps_matrices[i].rows()/hilbert_dim<<" "<<mps_matrices[i].cols()<<endl;
	}

	return 0;
}





