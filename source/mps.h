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
	void make_left_canonical(void); //Make the matrix left canonical (calls private function sweep left many times)
	void make_right_canonical(void); //Make the matrix
	void left_compress(unsigned);
	void right_compress(unsigned);

	//Change storage of an mps from left-storage to right
	void change_mps_storage_to_right(void);
	//Change storage of an mps from right-storage to left
	void change_mps_storage_to_left(void);

	//Some functions for debugging
	void return_matrix(unsigned,unsigned); //returns matrix at site and dim
	void return_array_list(); //Simply return the array list
	int validate_MPS(void); //Check all the data stored in memory is valid and consistent. For example check that the stored_MPS is consistent with the stored_matrix_dimensions


	//Temporary functions for debugging
	CanonMat return_matrix_at_site(unsigned);

	/* FUNCTIONS required for more than one method
	 * within the class.
	 */

	//Multiplying canonical matrices for sweep from left
	template<typename DerivedA, typename DerivedB>
	DerivedB left_canon_mult(const MatrixBase<DerivedA>& ,
			const MatrixBase<DerivedB>&);

	//Multiplying canonical matrices for sweep from right
	template<typename DerivedA, typename DerivedB>
	DerivedB right_canon_mult(const MatrixBase<DerivedA>& ,
			const MatrixBase<DerivedB>&);

	//Change from left to right storage for canonical matrices
	template<typename Derived>
	Derived left_storage_to_right(const MatrixBase<Derived>&);

	//Change from left to right storage for canonical matrices
	template<typename Derived>
	Derived right_storage_to_left(const MatrixBase<Derived>&);


};

#endif /* MPS_H_ */
