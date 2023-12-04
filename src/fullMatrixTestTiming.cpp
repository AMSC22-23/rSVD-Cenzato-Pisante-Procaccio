#include "fullMatrix.hpp"
#include <chrono>

int main(int argc, char**argv){

	int n=(argc==2)?std::atoi(argv[1]):1000;

	#ifdef BYROWS
	FullMatrix<double,ORDERING::ROWMAJOR> m1(n,n,1.),m2(n,n,2.);
	#else
	FullMatrix<double,ORDERING::COLMAJOR> m1(n,n,1.),m2(n,n,2.);
	#endif

	const auto t0=std::chrono::high_resolution_clock::now();
	auto ret=m1*m2;
	const auto t1=std::chrono::high_resolution_clock::now();

	const auto dt=std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();

	std::cout<<"Elapsed time for matrix-matrix multiplication: "<<dt<<" [ms]"<<std::endl;

	return 0;
}
