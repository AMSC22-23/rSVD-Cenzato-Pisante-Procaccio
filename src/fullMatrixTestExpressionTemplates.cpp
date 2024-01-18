#include "fullMatrix.hpp"

int main(int argc, char* argv[]){

	int n=argc>1?std::atoi(argv[1]):10;
	int m=argc>2?std::atoi(argv[2]):5;

	std::cout<<"Testing expression templates with matrices of dimension "<<n<<"x"<<m<<std::endl;
	std::cout<<std::endl;
	
	//It does not matter the order
	using Real=double;
	using Matrix=FullMatrix<Real,ORDERING::ROWMAJOR>;

	Matrix a;

	a.resize(n,m,1.);

	std::cout<<"Matrix: "<<std::endl;
	std::cout<<std::endl;
	std::cout<<a<<std::endl;
	std::cout<<"Matrix transposed: "<<std::endl;
	std::cout<<std::endl;
	std::cout<<a.transpose()<<std::endl;

	Matrix b=a*a.transpose();

	std::cout<<"Multiplying the matrix by it transposed: "<<std::endl;
	std::cout<<std::endl;
	std::cout<<b<<std::endl;
	
	std::cout<<"Multiplying the matrix by it transposed and view its norm: "<<std::endl;
	std::cout<<std::endl;
	std::cout<<(a*a.transpose()).norm()<<std::endl<<std::endl;

	a=a*a.transpose();

	std::cout<<"Multiplying the matrix by it transposed but now it gives an error: "<<std::endl;
	std::cout<<std::endl;
	std::cout<<a<<std::endl;

	return 0;
}
