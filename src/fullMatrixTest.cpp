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

	std::cout<<"Test: Changed just one value: ";
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

	std::cout<<std::endl;
	std::cout<<"Now matrix changed to:"<<std::endl;
	testMatrix.print(std::cout);
	std::cout<<std::endl;

	std::cout<<"Test: Correctness of resize: ";
	if(testMatrix.rows()!=new_n||testMatrix.cols()!=new_m)
		correct=false;
	std::cout<<correct<<std::endl;
	correct=true;

	//Refilling the matrix with numbers
	for(size_t i=0;i<new_n;++i)
		for(size_t j=0;j<new_m;++j)
			testMatrix[i][j]=i+j;

	std::cout<<std::endl;
	std::cout<<"Now matrix changed to:"<<std::endl;
	testMatrix.print(std::cout);
	std::cout<<std::endl;

	//Create a vector for which the matrix multiplication is not valid
	std::vector<Real> testVector(new_m-1,1.);

	std::cout<<"Using vector:"<<std::endl;
	for(size_t i=0;i<testVector.size();++i)
		std::cout<<testVector[i]<<std::endl;
	std::cout<<std::endl;

	std::cout<<"Test: Correctness of erroneous input in matrix-vector multiplication: ";
	std::vector<Real> testResult=testMatrix*testVector;
	if(testResult.size()!=0)
		correct=false;
	std::cout<<correct<<std::endl;
	correct=true;

	//Now correct the vector
	testVector.emplace_back(1.);

	std::cout<<std::endl;
	std::cout<<"Now using vector;"<<std::endl;
	for(size_t i=0;i<testVector.size();++i)
		std::cout<<testVector[i]<<std::endl;
	std::cout<<std::endl;

	std::cout<<"Test: Correctness of result in matrix-vector multiplication: ";
	testResult=testMatrix*testVector;
	Real correctVal=new_m*(new_m-1)/2;
	for(size_t i=0;i<testResult.size();++i)
		if(testResult[i]!=correctVal+i*new_m)
			correct=false;
	std::cout<<correct<<std::endl;

	//Now test the matrix-matrix multiplication
	FullMatrix<Real> testMultMatrix(new_m-1, 1, 1.);

	std::cout<<std::endl;
	std::cout<<"New matrix:"<<std::endl;
	testMultMatrix.print(std::cout);
	std::cout<<std::endl;

	std::cout<<"Test: Correctness of erroneous input in matrix-matrix multiplication: ";
	FullMatrix<Real> testMultResult=testMatrix*testMultMatrix;
	if(testMultResult.rows()==1&&testMultResult.cols()==1&&testMultResult[0][0]!=0.){
			correct=false;
	}
	std::cout<<correct<<std::endl;

	//Now test the correct matrix-matrix multiplication
	testMultMatrix.resize(testMultMatrix.rows()+1,1);
	testMultMatrix[testMultMatrix.rows()-1][0]=1.;

	std::cout<<std::endl;
	std::cout<<"Matrix updated to:"<<std::endl;
	testMultMatrix.print(std::cout);
	std::cout<<std::endl;

	std::cout<<"Test: Correctness of result in matrix-matrix multiplication: ";
	testMultResult=testMatrix*testMultMatrix;
	if(testMultResult.rows()!=testResult.size()&&testMultResult.cols()!=1)
		correct=false;
	else{
		for(size_t i=0;i<testMultResult.rows();++i){
			if(testMultResult[i][0]!=testResult[i])
				correct=false;
		}
	}
	std::cout<<correct<<std::endl;
	correct=true;

	//Noew check the correct transposed
	FullMatrix<Real> testTransposed=testMatrix.transpose();

	std::cout<<std::endl;
	std::cout<<"Matrix transposed:"<<std::endl;
	testTransposed.print(std::cout);
	std::cout<<std::endl;

	std::cout<<"Test: Correctness of the transposed: ";
	for(size_t i=0;i<testTransposed.rows();++i)
		for(size_t j=0;j<testTransposed.cols();++j)
			if(testTransposed[i][j]!=testMatrix[j][i])
				correct=false;
	std::cout<<correct<<std::endl;
	correct=true;

	return 0;
}









