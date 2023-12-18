#ifndef HH__UTILS__HH
#define HH__UTILS__HH

#include <iostream>
#include <cmath>
#include <omp.h>
#include <tuple>
#include <chrono>
#include <random>
#include <fstream>
#include <sstream> 
#include <iomanip>
//@note use at least namespaces to avoid polluting the global namespace with 
//      the names defined in the headers
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