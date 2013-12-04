/* This is a sample program which will explore the use of the Eigen
 * package for linear algebra.
 * Specifically of interest are the array and matrix classes
 * and the SVD routines provided.
 */

#include "Eigen/Dense"
#include <iostream>
using namespace Eigen;

int main()
{
  MatrixXd m(2,5);
  m.resize(4,3);
  std::cout << "The matrix m is of size "
            << m.rows() << "x" << m.cols() << std::endl;
  std::cout << "It has " << m.size() << " coefficients";
  std::cout  << std::endl;
  VectorXd v(2);
  v.resize(5);
  std::cout << "The vector v is of size " << v.size();
  std::cout << std::endl;
  std::cout << "As a matrix, v is of size "
            << v.rows() << "x" << v.cols() << std::endl;
}
