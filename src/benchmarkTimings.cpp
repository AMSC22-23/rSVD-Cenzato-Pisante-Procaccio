#include "fullMatrix.hpp"
#include <chrono>
#include <algorithm>
#include <tuple>
#include <iomanip>

/*
	 This benchmark performs some operations for different
	 sizes of the matrix.

	 The output is in the csv format and it indicates the size,the time with my class (one thread),
	 my class with the maximum parallelization with openmp and eigen
 */

#define _CRT_NONSTDC_NO_WARNINGS

using Real=double;

#ifdef BYROWS
using Matrix=FullMatrix<Real,ORDERING::ROWMAJOR>;
#else
using Matrix=FullMatrix<Real,ORDERING::ROWMAJOR>;
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef EIGEN
#include <Eigen/Dense>
using EigenMatrix=Eigen::MatrixXd;
#endif

template<class M>
inline void constructHilbertMatrix(M &h){
	for(size_t i=0;i<h.rows();++i)
		for(size_t j=0;j<h.cols();++j)
			h(i,j)=1/static_cast<double>(i+j+1);
}

template<class M>
void operation(M &a, M &b, M &c){

#if defined LAZYPRODLONG

	c=a*2.*2.*2.*2.*2.*2.*2.*2.*2.;

#elif defined LAZYPRODSHORT

	c=a*2.;

#elif defined LAZYSUMLONG

	c=a+b+a+b+a+b+a+b+a+b;

#elif defined LAZYSUMSHORT

	c=a+b;

#elif defined MIXED 

	c=a*(a+b*2.)*a.transpose();

#else

	c=a*b;

#endif

}

std::tuple<double,double,double,double> 
getInterpolationCoeff(const size_t& _s, const size_t& _e, const size_t& _f){

	//Cast
	const auto s=static_cast<double>(_s);
	const auto e=static_cast<double>(_e);
	const auto f=static_cast<double>(_f);

	//Some useful numbers
	const auto s3=s*s*s;
	const auto s2=s*s;
	const auto e3=e*e*e;
	const auto e2=e*e;

	//Entries of the original matrix
	const auto a11=s3;
	const auto a12=s2;
	const auto a13=s;
	const auto a14=1.;
	const auto a21=e3/8;
	const auto a22=e2/4;
	const auto a23=e/2;
	const auto a24=1.;
	const auto a31=6*e/*27/64*e3*/;
	const auto a32=2./*9/16*e2*/;
	const auto a33=0./*3/4*e*/;
	const auto a34=0./*1.*/;
	const auto a41=e3;
	const auto a42=e2;
	const auto a43=e;
	const auto a44=1.;

	//Construct rhs
	const auto b1=e;
	const auto b2=e/8;
	const auto b3=0./*e/64*/;
	const auto b4=f;
	
	//Determinant
	const auto det=+a11*(
									  a22*a33*a44+a23*a34*a42+a24*a32*a43
									 -a24*a33*a42-a23*a32*a44-a22*a34*a43 )
								 -a21*(
								    a12*a33*a44+a13*a34*a42+a14*a32*a43
									 -a14*a33*a42-a13*a32*a44-a12*a34*a43)
								 +a31*(
								 		a12*a23*a44+a13*a24*a42+a14*a22*a43
									 -a14*a23*a42-a13*a22*a44-a12*a24*a43)
								 -a41*(
								 		a12*a23*a34+a13*a24*a32+a14*a22*a33
									 -a14*a23*a32-a13*a22*a34-a12*a24*a33) ;

	//Entries of the inverted matrix
	const auto A11=+a22*a33*a44+a23*a34*a42+a24*a32*a43
								 -a24*a33*a42-a23*a32*a44-a22*a34*a43 ;
	const auto A12=-a12*a33*a44-a13*a34*a42-a14*a32*a43
								 +a14*a33*a42+a13*a32*a44+a12*a34*a43 ;
	const auto A13=+a12*a23*a44+a13*a24*a42+a14*a22*a43
								 -a14*a23*a42-a13*a22*a44-a12*a24*a43 ;
	const auto A14=-a12*a23*a34-a13*a24*a32-a14*a22*a33
								 +a14*a23*a32+a13*a22*a34+a12*a24*a33 ;

	const auto A21=-a21*a33*a44-a23*a34*a41-a24*a31*a43
								 +a24*a33*a41+a23*a31*a44+a21*a34*a43 ;
	const auto A22=+a11*a33*a44+a13*a34*a41+a14*a31*a43
								 -a14*a33*a41-a13*a31*a44-a11*a34*a43 ;
	const auto A23=-a11*a23*a44-a13*a24*a41-a14*a21*a43
								 +a14*a23*a41+a13*a21*a44+a11*a24*a43 ;
	const auto A24=+a11*a23*a34+a13*a24*a31+a14*a21*a33
								 -a14*a23*a31-a13*a21*a34-a11*a24*a33 ;

	const auto A31=+a21*a32*a44+a22*a34*a41+a24*a31*a42
								 -a24*a32*a41-a22*a31*a44-a21*a34*a42 ;
	const auto A32=-a11*a32*a44-a12*a34*a41-a14*a31*a42
								 +a14*a32*a41+a12*a31*a44+a11*a34*a42 ;
	const auto A33=+a11*a22*a44+a12*a24*a41+a14*a21*a42
								 -a14*a22*a41-a12*a21*a44-a11*a24*a42 ;
	const auto A34=-a11*a22*a34-a12*a24*a31-a14*a21*a32
								 +a14*a22*a31+a12*a21*a34+a11*a24*a32 ;

	const auto A41=-a21*a32*a43-a22*a33*a41-a23*a31*a42
								 +a23*a32*a41+a22*a31*a43+a21*a33*a42 ;
	const auto A42=+a11*a32*a43+a12*a33*a41+a13*a31*a42
								 -a13*a32*a41-a12*a31*a43-a11*a33*a42 ;
	const auto A43=-a11*a22*a43-a12*a23*a41-a13*a21*a42
								 +a13*a22*a41+a12*a21*a43+a11*a23*a42 ;
	const auto A44=+a11*a22*a33+a12*a23*a31+a13*a21*a32
								 -a13*a22*a31-a12*a21*a33-a11*a23*a32 ;

	//Now construct the results of the matrix vector operations

	const auto ret1=A11/det*b1+A12/det*b2+A13/det*b3+A14/det*b4;
	const auto ret2=A21/det*b1+A22/det*b2+A23/det*b3+A24/det*b4;
	const auto ret3=A31/det*b1+A32/det*b2+A33/det*b3+A34/det*b4;
	const auto ret4=A41/det*b1+A42/det*b2+A43/det*b3+A44/det*b4;

	return std::make_tuple(ret1,ret2,ret3,ret4);

}

int main(int argc, char** argv){

	std::cout<<std::setprecision(10);

	size_t start=(argc>=2)?std::stoul(argv[1]):2;
	size_t end  =(argc>=3)?std::stoul(argv[2]):1024;
	size_t step =(argc>=4)?std::stoul(argv[3]):100;

	Matrix a,b,c;

#ifdef EIGEN
	EigenMatrix e,f,g;
#endif

	std::cout<<"N,Times,Mine";

#ifdef _OPENMP
	std::cout<<",OpenMP";
#endif

#ifdef EIGEN
	std::cout<<",Eigen";
#endif

	std::cout<<std::endl;

	std::vector<Real> timings;

	double tau=0.;

	//Get the interpolant values
	//auto [c1,c2,c3,c4]=getInterpolationCoeff(start,end,2U);

#if defined LAZYPRODLONG

	tau=10.;	

#elif defined LAZYPRODSHORT

	tau=5.;	

#elif defined LAZYSUMLONG

	tau=10.;	

#elif defined LAZYSUMSHORT

	tau=5.;	

#elif defined MIXED 

	tau=100.;	

#else

	tau=100.;	

#endif

	for(size_t i=start;i<=end;i+=step){

			/*
				In order to reduce the misuration error with small matrices
				we perform the test multiple times and get the minimum
			*/

			timings.clear();

			a.resize(i,i);
			b.resize(i,i);

			constructHilbertMatrix(a);
			constructHilbertMatrix(b);

			//NB: the linear interpolation is good for something that grows linearly
			//size_t times=(end*end+2*i-i*end-2*start)/(end-start);

			const auto i_=static_cast<double>(i);
			
			//size_t times=static_cast<size_t>(c1*i_*i_*i_+c2*i_*i_+c3*i_+c4);

			//We can try an exponential 
			//The exponential is best
			size_t times=static_cast<size_t>(2.+(end-2.)*std::exp((tau/end)*(start-i_)));

			timings.reserve(times);

			std::cout<<i<<","<<times;

#ifdef _OPENMP

			omp_set_num_threads(1);
			
#endif

			for(size_t j=0;j<times;++j){

				const auto t0=std::chrono::high_resolution_clock::now();
				operation(a,b,c);
				const auto t1=std::chrono::high_resolution_clock::now();

				const auto dt=std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
				
				timings.emplace_back(dt);

			}

			//From all the tests i did take the minimum

			const auto iter=std::min_element(timings.cbegin(),timings.cend());
			std::cout<<","<<*iter;

#ifdef _OPENMP
			
			//This are my optimal threads, i can achieve more than that (up to 20) but the performance decreases
			omp_set_num_threads(12);

			timings.clear();
			timings.shrink_to_fit();

			timings.reserve(times);

			for(size_t j=0;j<times;++j){

				const auto t0=std::chrono::high_resolution_clock::now();
				operation(a,b,c);
				const auto t1=std::chrono::high_resolution_clock::now();

				const auto dt=std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
				
				timings.emplace_back(dt);

			}

			//From all the tests i did take the minimum

			const auto itero=std::min_element(timings.cbegin(),timings.cend());
			std::cout<<","<<*itero;


#endif

#ifdef EIGEN

			e.resize(i,i);
			f.resize(i,i);

			constructHilbertMatrix(e);
			constructHilbertMatrix(f);

			timings.clear();
			timings.shrink_to_fit();

			timings.reserve(times);

			for(size_t j=0;j<times;++j){

				const auto t0=std::chrono::high_resolution_clock::now();
				operation(e,f,g);	
				const auto t1=std::chrono::high_resolution_clock::now();

				const auto dt=std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();
				
				timings.emplace_back(dt);

			}

			//From all the tests i did take the minimum

			const auto itere=std::min_element(timings.cbegin(),timings.cend());
			std::cout<<","<<*itere;

#endif

			std::cout<<std::endl;
			
	}

	return 0;
}
