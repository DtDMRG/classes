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

	//Define the column matrix of the down state
	CanonMat init= CanonMat( db ,1 );

    for (int i=0;i<db-1;i++){
    	init(i,0)= complex(double(i),0.);
    }
    init(db-1,0)=1.0;


    //Resize the arrays to the right size
    stored_matrix_dimensions= new int[ns+1];
    //mps_pointers = new const CanonMat*[ns];
    //mps_pointers.resize(ns);
    //myVector.push_back(Base_p(new Derived));

    //Push back to set a whole load of new Canon_Mat_ptrs of the right size
    for (int i=0; i<ns; ++i){
    	mps_pointers.push_back(CanonMat_ptr(new CanonMat(db, 1)) );
    }

    //Deference vector (1star) and deference shared_ptr (second star) to set
    for ( CanonMat_itr it = mps_pointers.begin(); it != mps_pointers.end(); it++ ){
    	**it=init;
    }

    //This is to test it is working here --- TO BE REMOVED
    for ( CanonMat_itr it = mps_pointers.begin(); it != mps_pointers.end(); it++ ){
        	//**it=CanonMat::Ones(db,1);
        	cout<< "allocated CanonMat:"<<endl;
        	cout<< **it <<endl <<endl;
    }



    //Initialize the dimensions vector to one for each site
	for (int i=0;i<ns;i++){
		stored_matrix_dimensions[i]=1;
	//	mps_pointers[i]=&init;
	}
	stored_matrix_dimensions[ns]=1;

}

Mps::Mps(QState input_qstate) {

	//Extract and store common variables
	n_sites = input_qstate.return_n_sites();
	hilbert_dim = input_qstate.return_hilbert_dim();

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





