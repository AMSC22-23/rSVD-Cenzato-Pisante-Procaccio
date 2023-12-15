/*
	 This is to prevent the .hpp file to be reused 
	 multiple times when compiling (it is an header guard)
 */
#ifndef HH__FULL_MATRIX__HH
#define HH__FULL_MATRIX__HH

#include <vector>
#include <iostream>
#include <cmath>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"

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
		FullMatrix(const size_t n,const size_t m,const Real initVal=0) : m_rows(n), m_cols(m) {
			m_entries.reserve(n*m);
			for(size_t i=0;i<n*m;++i){
				m_entries.emplace_back(initVal);
			}
		}
		FullMatrix(const size_t n): FullMatrix(n,1)
		{}

		FullMatrix() =default;

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
		const Real& coeffRef(const size_t i, const size_t j) const{
			return this->operator()(i,j);
		}
		/*
			 Another operator for accessing an element directly
		 */
		Real& coeffRef(const size_t i) {
			return m_entries[i];
		}
		/*
			 Another operator for accessing an element in read only operations
		 */
		const Real& coeffRef(const size_t i) const{
			return m_entries[i];
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
					auto ofs=A.cols()*i;
					Real sum=0.;
#pragma omp parallel for shared(ofs,x,A) reduction(+:sum)
					for(size_t j=0;j<x.size();++j){
						sum+=A.m_entries[ofs+j]*x[j];
					}
					toReturn[i]=sum;
				}
			}
			else{
				for(size_t j=0;j<A.cols();++j){
					auto ofs=A.rows()*j;
#pragma omp parallel for shared(ofs,x,A,toReturn)
					for(size_t i=0;i<A.rows();++i){
						toReturn[i]+=A.m_entries[i+ofs]*x[i];
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
#pragma omp parallel for shared(A,B,toReturn)
				for(size_t i=0;i<A.rows();++i){
					for(size_t k=0;k<B.rows();++k){
						auto ofs1=i*B.cols();
						auto ofs2=i*A.cols();
						auto ofs3=k*B.cols();
						//#pragma omp parallel for shared(ofs1,ofs2,ofs3,k,A,B,toReturn)
						for(size_t j=0;j<B.cols();++j){
							toReturn.m_entries[ofs1+j]+=A.m_entries[ofs2+k]*B.m_entries[ofs3+j];
						}
					}
				}
			}
			else{
#pragma omp parallel for shared(A,B,toReturn)
				for(size_t j=0;j<B.cols();++j){
					for(size_t k=0;k<A.cols();++k){
						auto ofs1=A.rows()*k;
						auto ofs2=B.rows()*j;
						//#pragma omp parallel for shared(ofs1,ofs2,A,B,toReturn) 
						for(size_t i=0;i<A.rows();++i){
							toReturn.m_entries[i+ofs2]+=A.m_entries[i+ofs1]*B.m_entries[k+ofs2];
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
		void row(const size_t _row, const std::vector<Real>& toInsert){
			if(_row>=m_rows || toInsert.size()>m_cols)
				return;

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				for(size_t j=0;j<toInsert.size();++j){
					m_entries[_row*m_cols+j]=toInsert[j];
				}
			}
			else{
				for(size_t j=0;j<toInsert.size();++j){
					m_entries[_row+m_rows*j]=toInsert[j];
				}
			}
		}
		/*
			 Overload of the row method for adding a one dimensional matrix
			 We need to check that it is a 1D vector
		 */
		void row(const size_t _row, const FullMatrix& toInsert){
			if(toInsert.m_rows!=1&&toInsert.m_cols!=1)
				return;
			row(_row,toInsert.m_entries);
		}
		/*
			 Add a vector to a specified column.
Hyphothesis: col<m_cols && toInsert.size()<=m_rows
		 */
		void col(const size_t _col, const std::vector<Real>& toInsert){
			if(_col>=m_cols || toInsert.size()>m_rows)
				return;

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				for(size_t i=0;i<toInsert.size();++i){
					m_entries[i*m_cols+_col]=toInsert[i];
				}
			}
			else{
				for(size_t i=0;i<toInsert.size();++i){
					m_entries[i+m_rows*_col]=toInsert[i];
				}
			}
		}
		/*
			 Overload of the col method for adding a one dimensional matrix
			 We need to check that it is a 1D vector
		 */
		void col(const size_t _col, const FullMatrix& toInsert){
			if(toInsert.m_rows!=1&&toInsert.m_cols!=1)
				return;
			col(_col,toInsert.m_entries);
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
			 Clear the matrix
			 It will completely empty the matrix (indeed m_rows==0 and m_cols==0) (USE WITH CARE)
		 */
		void clear(){
			m_entries.clear();
			m_entries.shrink_to_fit();
			m_rows=0;
			m_cols=0;
		}
		/*
			 Eliminate the last n columns of the matrix
			 Return true if operation was successful (we want n<m_cols)
		 */
		void trimCols(const size_t n){
			if(n>=m_cols)
				return;

			const auto toShift=m_cols-n;

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				for(auto i=m_entries.begin()+toShift;i<=m_entries.end();i+=toShift){
					m_entries.erase(i,i+n);
				}
			}
			else{
				m_entries.resize(m_rows*toShift);
			}
			m_cols-=n;

		}
		/*
			 Eliminate the last n rows of the matrix
			 Return true if operation was successful (we want n<m_rows)
		 */
		void trimRows(const size_t n){
			if(n>=m_rows)
				return;

			const auto toShift=m_rows-n;

			if constexpr(ORDER==ORDERING::ROWMAJOR){
				m_entries.resize(toShift*m_cols);
			}
			else{
				for(auto i=m_entries.begin()+toShift;i<=m_entries.end();i+=toShift){
					m_entries.erase(i,i+n);
				}
			}
			m_rows-=n;

		}
		/*
			 Method to multiply the current matrix by a constant
		 */
		void scale(const Real k){
			for(size_t i=0;i<m_entries.size();++i)
				m_entries[i]*=k;
		}

		/*
			 Method for retrieving the Frobenius norm squared
		 */
		const Real normSquared() const {
			Real norm=0.;
			for(auto el: m_entries)
				norm+=el*el;

			return norm;
		}
		/*
			 Method for retrieving the Frobenius norm
		 */
		const Real norm() const {
			return std::sqrt(normSquared());
		}
		/*
			 Method to return the transposed of a matrix
TODO: you can return a wrapper class AdjointMatrix with a reference to the class A but
that overloads the operator () in a COLMAJOR way (if matrix is saved as ROWMAJOR)
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
		/*
			 Another method for printing a matrix
		 */
		friend std::ostream& operator<<(std::ostream& os, FullMatrix const & mat){

			for(size_t i=0;i<mat.m_rows;++i){
				for(size_t j=0;j<mat.m_cols;++j){
					if constexpr(ORDER==ORDERING::ROWMAJOR)
						os<<mat.m_entries[i*mat.m_cols+j]<<" ";
					else
						os<<mat.m_entries[i+mat.m_rows*j]<<" ";
				}
				os<<std::endl;
			}

			return os;
		}
		/*
			 I can work on the assumption that the matrix is squared
		 */
		void setIdentity(){
			resize(m_rows,m_cols);
			for(size_t i=0;i<m_rows;++i)
				this->operator()(i,i)=1.;
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

#pragma GCC diagnostic pop

#endif
