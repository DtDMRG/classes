#ifndef MPS_H_
#define MPS_H_


#include "matrixdefs.h"

#include "qstate.h"


class Mps {

	/*
	 * Private memory where things are stored
	 */

	unsigned n_sites, hilbert_dim;

	unsigned* stored_matrix_dimensions;
	//const CanonMat** mps_pointers;

	vector<CanonMat> mps_matrices;

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
	Mps(unsigned, unsigned);


	Mps(unsigned, unsigned, char);
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

	void sweep_from_left_at(unsigned); //Perform an SVD and then change the matrix at site and also pass the residue to site plus 1
	void sweep_from_right_at(unsigned); //*
	void trunc_sweep_from_left_at(unsigned,unsigned); //Perform an SVD and then change the matrix at site and also pass the residue to site plus 1
	void trunc_sweep_from_right_at(unsigned,unsigned); //*

	//Functions for manipulating MPS
	void make_left_canonincal(); //Make the matrix left canonical (calls private function sweep left many times)
	void make_right_canonical(); //Make the matrix
	void compress(unsigned);




	//Some functions for debugging
	void return_matrix(unsigned,unsigned); //returns matrix at site and dim
	void return_array_list(); //Simply return the array list
	int validate_MPS(); //Check all the data stored in memory is valid and consistent. For example check that the stored_MPS is consistent with the stored_matrix_dimensions

	//Temporary functions for debugging
	CanonMat return_matrix_at_site(unsigned);

	//Multiplying canonical matrices for sweep from left
	template<typename DerivedA, typename DerivedB>
	DerivedB leftcanonmult(const MatrixBase<DerivedA>& ,
			const MatrixBase<DerivedB>&, unsigned, unsigned);


};

#endif /* MPS_H_ */
