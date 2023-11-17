/*
	 This is to prevent the .hpp file to be reused 
	 multiple times when compiling (it is an header guard)
 */
#pragma once

#include <vector>
#include <iostream>

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

	public:

		/*
			 Alias so that the type of the element is more clear
		 */
		using Real=T;

		/*
			 Constructor of a n x m matrix with a value initVal
		 */
		FullMatrix(const size_t n=1,const size_t m=1,const Real initVal=0){
			m_entries.reserve(n);
			for(size_t i=0;i<m;i++){
				m_entries.emplace_back(m,initVal);
			}
		}

		/*
			 Override of the operator [] so that the element at i,j 
			 can be accessed as: 
			 auto elem=A[i][j];
		 */
		std::vector<Real>& operator[](const size_t index){
			return m_entries[index];
		}

		/*
			 Override of the operator [] for read only operations
		 */
		const std::vector<Real>& operator[](const size_t index) const{
			return m_entries[index];
		}

		/*
			 Override of the operator () so that the element in position i,j
			 can be accessed as:
			 auto elem=A(i,j);
		 */
		Real& operator()(const size_t i, const size_t j){
			return m_entries[i][j];
		}

		/*
			 Override of the operator () for read only operations
		 */
		const Real& operator()(const size_t i, const size_t j) const{
			return m_entries[i][j];
		}

		/*
			 Method for printing all elements in the matrix in a generic stream
			 (same as the professor did in one of the labs)
		 */
		void print(std::ostream &os=std::cout) const{
			for(auto row:m_entries){
				for(auto obj:row){
					os<<obj<<" ";
				}
				os<<std::endl;
			}
		}

	private:

		/*
			 Elements of the matrix
		 */
		std::vector<std::vector<Real>> m_entries;
};
