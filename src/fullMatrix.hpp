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
			for(size_t i=0;i<n;++i){
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
			Overload of the multiplication operator for matrix-vector
		*/
		friend std::vector<Real> operator*(const FullMatrix& A, const std::vector<Real>& x){
			
			std::vector<Real> toReturn;

			//Check if matrix-vector multiplication is feasible
			if(A.cols()!=x.size())
				return toReturn;

			toReturn.reserve(A.rows());
			Real sum=0.;

			for(size_t i=0;i<A.rows();++i){
				sum=0.;
				for(size_t j=0;j<x.size();++j){
					sum+=A.m_entries[i][j]*x[j];
				}
				toReturn.emplace_back(sum);
			}

			return toReturn;
		}

		/*
			Overload of the multiplication operator for matrix-matrix
		*/
		friend FullMatrix operator*(const FullMatrix& A, const FullMatrix& B){
			
			FullMatrix toReturn;
			
			//Check if matrix-matrix multiplication is feasible
			if(A.cols()!=B.rows())
				return toReturn;

			toReturn.resize(A.rows(),B.cols());

			Real sum=0.;

			for(size_t i=0;i<A.rows();++i){
				for(size_t j=0;j<B.cols();++j){
					sum=0.;

					for(size_t k=0;k<B.rows();++k){
						sum+=A[i][k]*B[k][j];
					}

					toReturn[i][j]=sum;
				}
			}

			return toReturn;
		}

		/*
			Method to access the number of rows of the matrix
		*/
		const size_t rows() const{
			return m_entries.size();
		}

		/*
			Method to access the number of cols of the matrix
		*/
		const size_t cols() const{
			return m_entries[0].size();
		}

		/*
			Static method that constructs a new matrix filled with zeros
		*/
		static FullMatrix& Zero(const size_t n, const size_t m){
			static FullMatrix toReturn(n,m);
			return toReturn;
		}

		/*
			Resize of the matrix
			Note that all the elements added will be defaulted, and all elements in excess are truncated 
		*/
		void resize(const size_t n, const size_t m){
			m_entries.resize(n);
			for(size_t i=0;i<n;++i)
				m_entries[i].resize(m);
		}

		/*
			Method to return the transposed of a matrix
		*/
		FullMatrix transpose() const{
			FullMatrix toReturn(cols(),rows());

			for(size_t i=0;i<rows();++i)
				for(size_t j=0;j<cols();++j)
					toReturn[j][i]=m_entries[i][j];

			return toReturn;
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
