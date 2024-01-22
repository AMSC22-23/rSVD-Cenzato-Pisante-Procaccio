#include "fullMatrix.hpp"
#include <chrono>

int main(int argc, char* argv[]){

	int n=argc>1?std::atoi(argv[1]):100;
	int m=argc>2?std::atoi(argv[2]):100;

	std::cout<<"Testing expression templates with matrices of dimension "<<n<<"x"<<m<<std::endl;
	std::cout<<std::endl;
	
	//It does not matter the order
	using Real=double;
	using Matrix=FullMatrix<Real,ORDERING::ROWMAJOR>;

	Matrix a;

	a.resize(n,m,1.);

	std::cout<<"Doing some operations with scalar multiplication."<<std::endl<<std::endl;
	std::cout<<"In particular doing 2*2*2*2*2*2*2*2*2*A where A is only ones."<<std::endl<<std::endl;

	std::cout<<"Starting...";

	const auto t0=std::chrono::high_resolution_clock::now();
	//NB: these have to be floating points
	Matrix res=2.*2.*2.*2.*2.*2.*2.*2.*2.*a;
	const auto t1=std::chrono::high_resolution_clock::now();

	const auto dt=std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count();

	std::cout<<" done"<<std::endl<<std::endl;

	std::cout<<"Timing: "<<dt<< "[mus]"<<std::endl;
	std::cout<<"Norm:   "<<res.norm()<<std::endl<<std::endl;

	return 0;
}
