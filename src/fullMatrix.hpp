/*
	 This is to prevent the .hpp file to be reused 
	 multiple times when compiling (it is an header guard)
 */
#ifndef HH__FULL_MATRIX__HH
#define HH__FULL_MATRIX__HH

#include <vector>
#include <iostream>

/*
	 Enum for the storage of the ordering information of the matrix
 */
enum class ORDERING{
	ROWMAJOR=0,
	COLMAJOR=1
};

/*
	 T is the type of the element contained in the matrix.
	 Can be int, float, double...
	 ORDER is the type of ordering that you want for matrix storage
 */
template <typename T, ORDERING ORDER=ORDERING::COLMAJOR>
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
		FullMatrix(const size_t n=1,const size_t m=1,const Real initVal=0) : m_rows(n), m_cols(m) {
			m_entries.reserve(n*m);
			for(size_t i=0;i<n*m;++i){
				m_entries.emplace_back(initVal);
			}
		}

		/*
			 Override of the operator [] so that the element at i,j 
			 can be accessed as: 
			 auto elem=A[i][j];

			 std::vector<Real>& operator[](const size_t index){
			 return m_entries[index];
			 }
		 */

		/*
			 Override of the operator [] for read only operations

			 const std::vector<Real>& operator[](const size_t index) const{
			 return m_entries[index];
			 }
		 */

		/*
			 Override of the operator () so that the element in position i,j
			 can be accessed as:
			 auto elem=A(i,j);
		 */
		Real& operator()(const size_t i, const size_t j){
			if constexpr(ORDER==ORDERING::ROWMAJOR)
				return m_entries[i*m_cols+j];
			else
				return m_entries[i+m_rows*j];
		}
		/*
			 Override of the operator () for read only operations
		 */
		const Real& operator()(const size_t i, const size_t j) const{
			if constexpr(ORDER==ORDERING::ROWMAJOR)
				return m_entries[i*m_cols+j];
			else
				return m_entries[i+m_cols*j];
		}
		/*
			Another operator for accessing an element of the matrix
		*/
		Real& coeffRef(const size_t i, const size_t j){
			return this->operator()(i,j);
		}
		/*
			Another operator for accessing an element in read only operations
		*/
		const Real& coeffref(const size_t i, const size_t j) const{
			return this->operator()(i,j);
		}
		
		/*
			 Override of the multiplication for matrix-scalar
			 Note that it returns a new matrix
		 */
		friend FullMatrix operator*(const FullMatrix& A, const Real k){
			FullMatrix toReturn(A.rows(),A.cols());

			for(size_t i=0;i<A.m_entries.size();++i)
				toReturn.m_entries[i]=A.m_entries[i]*k;

			return toReturn;
		}
		/*
			 Override of the associative multiplication scalar-matrix
		 */
		friend FullMatrix operator*(const Real k, const FullMatrix& A){
			return A*k;
		}
		/*
			 Override of the operator for the summation (?) of two matrices
		 */
		friend FullMatrix operator+(const FullMatrix& A,const FullMatrix& B){
			FullMatrix toReturn;
			if(A.rows()!=B.rows() || A.cols()!=B.cols())
				return toReturn;

			toReturn.resize(A.rows(),A.cols());

			for(size_t i=0;i<A.m_entries.size();++i)
				toReturn.m_entries[i]=A.m_entries[i]+B.m_entries[i];

			return toReturn;
		}
		/*
			 Override of the operator for the subtraction of two matrices
		 */
		friend FullMatrix operator-(const FullMatrix& A,const FullMatrix& B){
			FullMatrix toReturn;
			if(A.rows()!=B.rows() || A.cols()!=B.cols())
				return toReturn;

			toReturn.resize(A.rows(),A.cols());

			for(size_t i=0;i<A.m_entries.size();++i)
				toReturn.m_entries[i]=A.m_entries[i]-B.m_entries[i];

			return toReturn;
		}
		/*
			 Overload of the multiplication operator for matrix-vector
		 */
		friend std::vector<Real> operator*(const FullMatrix& A, const std::vector<Real>& x){

			std::vector<Real> toReturn;

			//Check if matrix-vector multiplication is feasible
			if(A.cols()!=x.size())
				return toReturn;

			toReturn.resize(A.rows(),0.);

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				for(size_t i=0;i<A.rows();++i){
					for(size_t j=0;j<x.size();++j){
						toReturn[i]+=A.m_entries[i*A.cols()+j]*x[j];
					}
				}
			}
			else{
				for(size_t j=0;j<A.cols();++j){
					for(size_t i=0;i<A.rows();++i){
						toReturn[i]+=A.m_entries[i+A.rows()*j]*x[j];
					}
				}
			}

			return toReturn;
		}
		/*
			 Override of the multiplication operator for matrix-matrix
			 It is cache friendly: reference: https://siboehm.com/articles/22/Fast-MMM-on-CPU
		 */
		friend FullMatrix operator*(const FullMatrix& A, const FullMatrix& B){

			FullMatrix toReturn;

			//Check if matrix-matrix multiplication is feasible
			if(A.cols()!=B.rows())
				return toReturn;

			toReturn.resize(A.rows(),B.cols());

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				for(size_t i=0;i<A.rows();++i){
					for(size_t k=0;k<B.rows();++k){
						for(size_t j=0;j<B.cols();++j){
							toReturn.m_entries[i*B.cols()+j]+=A.m_entries[i*A.cols()+k]*B.m_entries[k*B.cols()+j];
						}
					}
				}
			}
			else{
				for(size_t j=0;j<B.cols();++j){
					for(size_t k=0;k<A.cols();++k){
						for(size_t i=0;i<A.rows();++i){
							toReturn.m_entries[i+B.rows()*j]+=A.m_entries[i+A.rows()*k]*B.m_entries[k+B.rows()*j];
						}
					}
				}
			}
			/*
				 for(size_t i=0;i<A.rows();++i){
				 for(size_t j=0;j<B.cols();++j){
				 for(size_t k=0;k<B.rows();++k){
				 toReturn[i][j]+=A[i][k]*B[k][j];
				 }
				 }
				 }
			 */

			return toReturn;
		}
		/*
			 Method to access the number of rows of the matrix
		 */
		const size_t& rows() const{
			return m_rows;
		}

		/*
			 Method to access the number of cols of the matrix
		 */
		const size_t& cols() const{
			return m_cols;
		}

		/*
			Add a vector to a specified row.
			Hyphothesis: row<m_rows && toInsert.size()<=m_cols
		*/
		void row(const size_t row, const std::vector<Real>& toInsert){
			if(row>=m_rows || toInsert.size()>m_cols)
				return;

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				for(size_t j=0;j<toInsert.size();++j){
					m_entries[row*m_cols+j]=toInsert[j];
				}
			}
			else{
				for(size_t j=0;j<toInsert.size();++j){
					m_entries[row+m_rows*j]=toInsert[j];
				}
			}
		}
		/*
			Add a vector to a specified column.
			Hyphothesis: col<m_cols && toInsert.size()<=m_rows
		*/
		void col(const size_t col, const std::vector<Real>& toInsert){
			if(col>=m_cols || toInsert.size()>m_rows)
				return;

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				for(size_t i=0;i<toInsert.size();++i){
					m_entries[i*m_cols+col]=toInsert[i];
				}
			}
			else{
				for(size_t i=0;i<toInsert.size();++i){
					m_entries[i+m_rows*col]=toInsert[i];
				}
			}
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
			m_entries.resize(n*m,0.);
			m_rows=n;
			m_cols=m;
		}
		/*
			 Method to multiply the current matrix by a constant
		 */
		void scale(const Real k){
			for(size_t i=0;i<m_entries.size();++i)
				m_entries[i]*=k;
		}
		/*
			 Method to return the transposed of a matrix
		 */
		FullMatrix transpose() const{
			FullMatrix toReturn(m_cols,m_rows);

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				for(size_t i=0;i<rows();++i)
					for(size_t j=0;j<cols();++j)
						toReturn.m_entries[j*m_rows+i]=m_entries[i*m_cols+j];
			}
			else{
				for(size_t j=0;j<cols();++j)
					for(size_t i=0;i<rows();++i)
						toReturn.m_entries[i*m_cols+j]=m_entries[j*m_rows+i];
			}

			return toReturn;
		}
		/*
			 Method for printing all elements in the matrix in a generic stream
			 (same as the professor did in one of the labs)
		 */
		void print(std::ostream &os=std::cout) const{
			
			for(size_t i=0;i<m_rows;++i){
				for(size_t j=0;j<m_cols;++j){
					if constexpr(ORDER==ORDERING::ROWMAJOR)
						os<<m_entries[i*m_cols+j]<<" ";
					else
						os<<m_entries[i+m_rows*j]<<" ";
				}
				os<<std::endl;
			}
		}

	private:

		/*
			 Elements of the matrix
		 */
		std::vector<Real> m_entries;
		/*
			 Rows of the matrix
		 */
		size_t m_rows;
		/*
			 Columns of the matrix
		 */
		size_t m_cols;
};

#endif
