/*
 * mps.cpp
 *
 *  Created on: Dec 5, 2013
 *      Author: emanuelelevi
 */

#include "mps.h"
#define complex std::complex<double>



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
	MatrixXc M(db,1);
    for (int i=0;i<db-1;i++){M(i,1)=complex(0.,0.);}
    M(db-1,1)=1.0;


    //Resize the arrays to the right size
    stored_matrix_dimensions= new int[ns+1];
    stored_mps = new MatrixXc[ns];



    //Initialize the dimensions vector to one for each site
	for (int i=0;i<ns;i++){
		stored_matrix_dimensions[i]=1;
		stored_mps[i]=M;
	}
	stored_matrix_dimensions[ns]=1;


}






/*DESTRUCTOR*/
Mps::~Mps() {}

