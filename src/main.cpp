#include <iostream>
#include "fullMatrix.hpp"

int main(){

	std::cout<<"Randomized SVD!"<<std::endl;
	
	FullMatrix<double> obj(3,3,1);

	std::cout<<obj[0][0]<<std::endl;
	obj[0][0]=2.;
	std::cout<<obj[0][0]<<std::endl;

	for(auto i:obj[0])
		std::cout<<i<<" ";

	std::cout<<std::endl;

	return 0;
}
