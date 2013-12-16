#ifndef MPS_H_
#define MPS_H_


#include "matrixdefs.h"

#include "qstate.h"


class Mps {

	/*
	 * Private memory where things are stored
	 */

	int n_sites, hilbert_dim;

	int* stored_matrix_dimensions;
	//const CanonMat** mps_pointers;

	vector<CanonMat> mps_matricies;

	/*
	* Private functions
	* After running a private function the validate MPS state is _not_ required to return true
	* SO ONLY USE ONE WHEN YOU NOW WHAT YOU'RE DOING
	*/

	/*
	 * Initializes memory for the matrix dimensions and all the pointers
	 * Requires the number of sites and Hilbert dimension have already been defined (bad programming my bad)
	 */
	void constructor_common_memory_init();

public:

	/*
	 * Constructors
	 */

	/*
	 * Default creator, states MPS in all high state
	 * given the 1 the number of sites
	 * and 2 the hilbert dimension of each site
	 */
	Mps(int, int);

	/*
	 * Create an MPS from a state
	 * only requires a QState everything else is simply extracted from that object
	 */
	Mps(QState);

	/*
	 * Destructor NOT IMPLEMENTED!!!
	 */
	~Mps();

	/*
	* Public functions
	* After running a public function the validate MPS state function MUST ALWAYS RETURN TRUE
	* Code taking that into account
	*/

	//Simple functions for returning internal variables

	void sweep_from_left_at(int); //Perform an SVD and then change the matrix at site and also pass the residue to site plus 1
	void sweep_from_right_at(int); //*
	void trunc_sweep_from_left_at(int,int); //Perform an SVD and then change the matrix at site and also pass the residue to site plus 1
	void trunc_sweep_from_right_at(int,int); //*

	//Functions for manipulating MPS
		void make_left_canonincal(); //Make the matrix left canonical (calls private function sweep left many times)
		void make_right_canonical(); //Make the matrix
		void compress(int);




	//Some functions for debugging
	void return_matrix(int,int); //returns matrix at site and dim
	void return_array_list(); //Simply return the array list
	void validate_MPS(); //Check all the data stored in memory is valid and consistent. For example check that the stored_MPS is consistent with the stored_matrix_dimensions

	//Temporary functions for debugging
	CanonMat return_matrix_at_site(int);

};

#endif /* MPS_H_ */
