#include "fullMatrix.hpp"

#include <unordered_map>
#include <functional>

/*
	Everything in this file is a test for the correct implementation of the class FullMatrix.
	It can (and it should) be compiled indipendently from the main file.
*/
int main(){

	//Create the dimensions of the test matrix
	constexpr size_t n=5;
	constexpr size_t m=3;
	constexpr size_t new_n=4;
	constexpr size_t new_m=6;
	bool correct=true;

	//Using an alias for faster change if needed 
	using Real=double;

	
	//Create the test matrix of dimensions n x m
	FullMatrix<Real> testMatrix=FullMatrix<Real>::Zero(n,m);
	
	std::cout<<"Executing tests on FullMatrix:"<<std::endl;
	testMatrix.print(std::cout);
	std::cout<<std::endl;

	std::cout<<"Test: Correctness of creation: ";
	for(size_t i=0;i<n;++i)
		for(size_t j=0;j<m;++j)
			if(testMatrix[i][j]!=0.)
				correct=false;
	std::cout<<correct<<std::endl;
	correct=true;

	//Change the value of one element and check if it is still correct
	testMatrix[n-1][0]=1.;

	std::cout<<std::endl;
	std::cout<<"Now matrix changed to:"<<std::endl;
	testMatrix.print(std::cout);
	std::cout<<std::endl;

	std::cout<<"Test: Changed just one value: "<<std::endl;
	for(size_t i=0;i<n;++i)
		for(size_t j=0;j<m;++j){
			if(i!=n-1||j!=0){
				if(testMatrix(i,j)!=0.)
					correct=false;
			}
			else{
				if(testMatrix(i,j)!=1.)
					correct=false;
			}
		}
	std::cout<<correct<<std::endl;
	correct=true;

	//Resize the matrix to be of new dimensions 
	testMatrix.resize(new_n,new_m);
	//Now the matrix should be again full zeros

	std::cout<<std::endl;
	std::cout<<"Now matrix changed to:"<<std::endl;
	testMatrix.print(std::cout);
	std::cout<<std::endl;

	std::cout<<"Test: Correctness of resize: ";
	for(size_t i=0;i<new_n;++i)
		for(size_t j=0;j<new_m;++j)
			if(testMatrix[i][j]!=0.)
				correct=false;
	std::cout<<correct<<std::endl;
	correct=true;

	return 0;
}
