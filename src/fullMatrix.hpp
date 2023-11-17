//This is to prevent the .hpp file to be reused 
//multiple times when compiling (it is an header guard)
#pragma once

#include <vector>

/*
	T is the type of the element contained in the matrix.
	Can be int, float, double...
*/
template <typename T>
/*
	Class that contains all the informations and methods for 
	the management of a generic full matrix
*/
class FullMatrix{

/*
	Alias so that the type of the element is more clear
*/
using Real=T;

public:

	/*
		Constructor of a n x m matrix with a value initVal
	*/
	FullMatrix(const size_t n=1,const size_t m=1,const Real initVal=0){
		entries.reserve(n);
		for(size_t i=0;i<m;i++){
			entries.emplace_back(m,initVal);
		}
	}

	/*
		Override of the operator [] so that the element at i,j 
		can be accessed as: 
		auto elem=A[i][j]
	*/
	std::vector<Real>& operator[](size_t index){
		return entries[index];
	}

private:

	/*
		Elements of the matrix
	*/
	std::vector<std::vector<Real>> entries;
};
