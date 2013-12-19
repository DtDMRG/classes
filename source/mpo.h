/*
 * mps.h
 *
 *  Created on: Dec 5, 2013
 *      Author: emanuelelevi
 */

#ifndef MPO_H_
#define MPO_H_


/* selected definitions from std
#define complex std::complex<double>

/*
This is the type definition of a matrix of complex entries.
Must be defined on a separate file.

typedef Matrix<complex, Dynamic, Dynamic> CanonMat;

/*
 * Define the matrix product operator class
 * Essentially an array of matricies for each site of the system
 * At each site this matrix is a composed of the four indexes
 * two related to the dimension of the hilbert space
 * two related to size of the matrix which can change

class Mpo {

	//Private variables
	int n_sites, hilbert_dim;

	int *matrix_dimensions;
	CanonMat *stored_mps;

public:

	//Constructors
	/*
	 * Create an MPO given first the number of sites then the Hilbert space dimension

	Mpo(int, int);


	//Destructor

	~Mpo();


	//Some functions for debugging
	void return_matrix(int,int,int); //returns matrix at site and 2 dims
	void return_array_list(); //Simply return the array list
	void validate_MPO(); //Check all the data stored in memory is valid and consistent. For example check that the stored_MPS is consistent with the stored_matrix_dimensions

};
*/

#endif /* MPS_H_ */
