/*
 * mps.h
 *
 *  Created on: Dec 5, 2013
 *      Author: emanuelelevi
 */

#ifndef MPS_H_
#define MPS_H_


#include "Eigen/Dense"
#include "Eigen/SVD"
#include "matrixdefs.h"
/* selected definitions from std*/
#define complex std::complex<double>


using namespace Eigen;



/*
This is the type definition of a matrix of complex entries.
Must be defined on a separate file.
*/
typedef Matrix<complex, Dynamic, Dynamic> CanonMat;


class Mps {


public:
	//Private members
	int n_sites, hilbert_dim;
	CanonMat init;



	//private function
	void sweep_from_left_at(int, const CanonMat** &); //Perform an SVD and then change the matrix at site and also pass the residue to site plus 1
	void sweep_from_right_at(int); //*
	void trunc_sweep_from_left_at(int,int); //Perform an SVD and then change the matrix at site and also pass the residue to site plus 1
	void trunc_sweep_from_right_at(int,int); //*



	int* stored_matrix_dimensions;
	//const CanonMat** mps_pointers;
	vector<CanonMat_ptr> mps_pointers;


	//Constructor, Destructor
	//Mps(int, int);
	Mps(int, int);
	~Mps();


	//Some functions for debugging
	void return_matrix(int,int); //returns matrix at site and dim
	void return_array_list(); //Simply return the array list
	void validate_MPS(); //Check all the data stored in memory is valid and consistent. For example check that the stored_MPS is consistent with the stored_matrix_dimensions




	//Functions for manipulating MPS
	void make_left_canonincal(); //Make the matrix left canonical (calls private function sweep left many times)
	void make_right_canonical(); //Make the matrix
	void compress(int);

};

#endif /* MPS_H_ */
