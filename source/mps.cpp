/*
 * mps.cpp
 *
 *  Created on: Dec 5, 2013
 *      Author: emanuelelevi
 */

#include "mps.h"

#include <iostream>
#define cout std::cout
#define endl std::endl


/*CONSTRUCTOR
 * Given number of sites and dimension of the quantum basis initializes the array dimensions
 * to all 1, and the array mmps to (dimbasis,1) matrices corresponding to the "lowest" basis
 * element.
 * */

Mps::Mps(int ns, int db) {

	//Set the value of numsites and dimbasis
	n_sites = ns;
	hilbert_dim = db;

	//Initialize memory for pointers and array of dimensions
	constructor_common_memory_init();


	//Define the column matrix of the down state
	CanonMat init= CanonMat( db ,1 );

    for (int i=0;i<db-1;i++){
    	init(i,0)= complex(double(i),0.);
    }
    init(db-1,0)=1.0;

    //Deference vector (1star) and deference shared_ptr (second star) to set
    for ( CanonMat_itr it = mps_pointers.begin(); it != mps_pointers.end(); it++ ){
    	**it=init;
    }

    //Initialize the dimensions vector to one for each site
	for (int i=0;i<n_sites;i++){
		stored_matrix_dimensions[i]=1;
	//	mps_pointers[i]=&init;
	}
	stored_matrix_dimensions[n_sites]=1;

	 //This is to test it is working here --- TO BE REMOVED
	    for ( CanonMat_itr it = mps_pointers.begin(); it != mps_pointers.end(); it++ ){
	        	//**it=CanonMat::Ones(db,1);
	        	cout<< "allocated CanonMat:"<<endl;
	        	cout<< **it <<endl <<endl;
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
	CanonMat_ptr current_state_residue;
	//To start we copy in the whole state stored passed at the input
	*current_state_residue = input_qstate.return_qstate();
	//Possible
	//(*current_state_residue).conservativeResize( input_qstate.return_hilbert_dim(), input_qstate.return_full_state_dim()/input_qstate.return_hilbert_dim());

	int last_rank = 1;

	stored_matrix_dimensions[0] = last_rank;

	//Iterator used to access all the pointers I need to get to
	CanonMat_itr current_mps_itr;
	current_mps_itr = mps_pointers.begin();

	for (int i=0; i<n_sites;i++) {
		Svd temp_svd(*current_state_residue);
		stored_matrix_dimensions[i+1] = temp_svd.rank();
		//(*stored_matrix_dimensions[1])stored_matrix_dimensions[i+1]
		//**current_mps_itr = temp_svd.matrixU().conservativeResize(hilbert_dim, hilbert_dim);

		current_mps_itr++;
	}

}

/*
MPS = Range[Ns];
 rankList = Range[Ns + 1];
 currentState = inputState;
 rlast = 1;
 rankList[[1]] = rlast;
 Do[
   currentMat =
    Partition[currentState, Length[currentState]/(rlast d)];
   {U, S, V, r} = CustomSVD[currentMat];
   rankList[[i + 1]] = r;
   MPS[[i]] = Transpose[ArrayReshape[U, {rlast, d, r}]];
   currentState = Flatten[ S.V\[ConjugateTranspose]];
   rlast = r;
   , {i, Ns} ];
 {MPS, rankList}]
*/

void Mps::constructor_common_memory_init() {

		//Resize the arrays to the right size

	    stored_matrix_dimensions= new int[n_sites+1];
	    //mps_pointers = new const CanonMat*[ns];
	    //mps_pointers.resize(ns);
	    //myVector.push_back(Base_p(new Derived));

	    //Push back to set a whole load of new Canon_Mat_ptrs of the right size
	    for (int i=0; i<n_sites; ++i){
	    	mps_pointers.push_back(CanonMat_ptr(new CanonMat(hilbert_dim, 1)) );
	    }



}



/*DESTRUCTOR*/
Mps::~Mps(){}



void Mps::sweep_from_left_at(int position, const CanonMat**& pointers){
	//cout<<pointers[position]<<endl;
	Svd localsweep(*pointers[position]);
	pointers[position]= & localsweep.matrixU();
	//cout<<pointers[position]<<endl;
	//cout<<*pointers[position]<<endl;
}





