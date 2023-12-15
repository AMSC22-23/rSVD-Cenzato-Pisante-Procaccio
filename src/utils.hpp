#ifndef HH__UTILS__HH
#define HH__UTILS__HH

#ifdef EIGEN
#include <Eigen/Dense>
using Matrix=Eigen::MatrixXd;
using Vector=Eigen::VectorXd;
#else
#include "fullMatrix.hpp"
using Matrix=FullMatrix<double,ORDERING::ROWMAJOR>;
using Vector=FullMatrix<double,ORDERING::ROWMAJOR>;
#endif

#endif