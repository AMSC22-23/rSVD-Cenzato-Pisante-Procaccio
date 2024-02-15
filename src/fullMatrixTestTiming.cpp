#include "fullMatrix.hpp"
#include <chrono>
#include <cmath>

using Real=double;

#ifdef BYROWS
using Matrix=FullMatrix<Real,ORDERING::ROWMAJOR>;
#else
using Matrix=FullMatrix<Real,ORDERING::COLMAJOR>;
#endif

#ifdef EIGEN
#include <Eigen/Dense>
using EigenMatrix=Eigen::MatrixXd;
void constructHilbertMatrix(EigenMatrix& h){
	for(auto i=0;i<h.rows();++i)
		for(auto j=0;j<h.cols();++j)
			h(i,j)=1/static_cast<double>(i+j+1);
}
Real checkMatrixDiff(const Matrix& a,const EigenMatrix& e){
	Real diff=0.;
	for(size_t i=0;i<a.rows();++i)
		for(size_t j=0;j<a.cols();++j)
			diff+=(a(i,j)-e(i,j))*(a(i,j)-e(i,j));
	return std::sqrt(diff);
}
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/*
Method to construct a MxN Hilbert Matrix (can also not be squared)
https://en.wikipedia.org/wiki/Hilbert_matrix
*/
void constructHilbertMatrix(Matrix& h){
	for(size_t i=0;i<h.rows();++i)
		for(size_t j=0;j<h.cols();++j)
			h(i,j)=1/static_cast<double>(i+j+1);
}

int main(int argc, char**argv){
	//@note its ok, but I suggest to use tools like getpot to parse command line arguments
	size_t m1=(argc>=2)?std::stoul(argv[1]):1000;
	size_t n1=(argc>=3)?std::stoul(argv[2]):1000;
	size_t m2=n1;
	size_t n2=(argc>=4)?std::stoul(argv[3]):1000;

	std::cout<<"-----------------------------------------------------------"<<std::endl;
	std::cout<<"|                                                         |"<<std::endl;
	std::cout<<"| Welcome to the Official Test for Matrix Multiplication! |"<<std::endl;
	std::cout<<"|                                                         |"<<std::endl;
	std::cout<<"-----------------------------------------------------------"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"  Based on your inputs and our sophisticated AI model"<<std::endl;
	std::cout<<"  we decided that these are the dimensions for the test matrices: "<<std::endl;
	std::cout<<std::endl;
	std::cout<<" Matrix 1:\t"<<"Matrix 2:\t"<<std::endl;
	std::cout<<" M: "<<m1<<"\t"<<"M: "<<m2<<"\t"<<std::endl;
	std::cout<<" N: "<<n1<<"\t"<<"N: "<<n2<<"\t"<<std::endl;
	std::cout<<std::endl;
	std::cout<<" Some other information for nerds: "<<std::endl;
	std::cout<<std::endl;
	std::cout<<"   Ordering: ";
	#ifdef BYROWS
	std::cout<<"ROWMAJOR"<<std::endl;
	#else
	std::cout<<"COLMAJOR"<<std::endl;
	#endif
	std::cout<<std::endl;
	std::cout<<"     OpenMP: ";
	#ifdef _OPENMP
	std::cout<<"ON "<<std::endl;
	std::cout<<"           # of threads active : "<<omp_get_num_threads()<<std::endl;
	std::cout<<"           # of maximum threads: "<<omp_get_max_threads()<<std::endl;
	#else
	std::cout<<"OFF"<<std::endl;
	#endif
	std::cout<<std::endl;
	std::cout<<"      Eigen: ";
	#ifdef EIGEN
	std::cout<<"ON "<<std::endl;
	#else
	std::cout<<"OFF"<<std::endl;
	#endif
	std::cout<<std::endl;
	std::cout<<"-----------------------------------------------------------"<<std::endl;
	std::cout<<std::endl;

	std::cout<<"  Creating the matrices... ";
	Matrix mat1(m1,n1);
	Matrix mat2(m2,n2);
	constructHilbertMatrix(mat1);
	constructHilbertMatrix(mat2);
	std::cout<<"Completed"<<std::endl;

	#ifdef EIGEN
	std::cout<<"  Creating Eigen matrices for test reference... ";
	EigenMatrix mat1e(m1,n1);
	EigenMatrix mat2e(m2,n2);
	constructHilbertMatrix(mat1e);
	constructHilbertMatrix(mat2e);
	std::cout<<"Completed"<<std::endl;
	std::cout<<"  Difference between the first  matrices: "<<checkMatrixDiff(mat1,mat1e)<<std::endl;
	std::cout<<"  Difference between the second matrices: "<<checkMatrixDiff(mat2,mat2e)<<std::endl;
	#endif

	std::cout<<std::endl;
	std::cout<<"  Multiplicating the two matrices... ";
	const auto t0=std::chrono::high_resolution_clock::now();
	//NB: they must not be auto as they will be returned as expressions
	Matrix ret=mat1*mat2;
	const auto t1=std::chrono::high_resolution_clock::now();
	std::cout<<"Completed"<<std::endl;
	const auto dt=std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();

	//std::cout<<ret<<std::endl;

	#ifdef _OPENMP
	std::cout<<"  Multiplicating the two matrices with just one thread... ";
	omp_set_num_threads(1);
	const auto t0p=std::chrono::high_resolution_clock::now();
	Matrix retp=mat1*mat2;
	const auto t1p=std::chrono::high_resolution_clock::now();
	std::cout<<"Completed"<<std::endl;
	const auto dtp=std::chrono::duration_cast<std::chrono::microseconds>(t1p-t0p).count();
	#endif

	#ifdef EIGEN
	std::cout<<"  Multiplicating the two Eigen matrices... ";
	const auto t0e=std::chrono::high_resolution_clock::now();
	EigenMatrix rete=mat1e*mat2e;
	const auto t1e=std::chrono::high_resolution_clock::now();
	std::cout<<"Completed"<<std::endl;
	std::cout<<std::endl;
	const auto dte=std::chrono::duration_cast<std::chrono::microseconds>(t1e-t0e).count();
	#endif
	
	std::cout<<std::endl;
	std::cout<<"-----------------------------------------------------------"<<std::endl;
	std::cout<<std::endl;

	std::cout<<"  Elapsed time         :\t"<< dt<<"\t[mus]"<<std::endl;
	
	#ifdef _OPENMP
	std::cout<<std::endl;
	std::cout<<"  Elapsed time 1 thread:\t"<<dtp<<"\t[mus]"<<std::endl;
	std::cout<<std::endl;
	std::cout<<"      Ratio Mine/OpenMP:\t"<<((double)dtp/dt)<<std::endl;
	#endif

	#ifdef EIGEN
	std::cout<<std::endl;
	std::cout<<"  Elapsed time Eigen   :\t"<<dte<<"\t[mus]"<<std::endl;
	
	std::cout<<std::endl;
	std::cout<<"       Ratio Mine/Eigen:\t"<<((double)dt/dte)<<std::endl;
	std::cout<<std::endl;
	std::cout<<"  Difference between results:\t"<<checkMatrixDiff(ret,rete)<<std::endl;
	#endif

	std::cout<<std::endl;
	std::cout<<"-----------------------------------------------------------"<<std::endl;

	return 0;
}
