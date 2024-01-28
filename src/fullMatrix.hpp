/*
	 This is to prevent the .hpp file to be reused 
	 multiple times when compiling (it is an header guard)
 */
#ifndef HH__FULLMATRIX__HH
#define HH__FULLMATRIX__HH

#include "expressionWrapper.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

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
	 Real is the type of the element contained in the matrix.
	 Can be float, double... (better not int)
	 ORDER is the type of ordering that you want for matrix storage
 */
template <typename Real, ORDERING ORDER=ORDERING::COLMAJOR>
/*
	 Class that contains all the informations and methods for 
	 the management of a generic full matrix

	 It derives from the expression class in order to activate expression templating
 */
class FullMatrix : public Expr<FullMatrix<Real,ORDER>> {

	public:

		/*
			 Constructor of a n x m matrix with a value initVal
		 */
		FullMatrix(const size_t n,const size_t m,const Real initVal=0): m_rows(n), m_cols(m) {
			m_entries.reserve(n*m);
			for(size_t i=0;i<n*m;++i){
				m_entries.emplace_back(initVal);
			}
		}
		/*
			 Constructor with just one dimension, useful for vector operations
		 */
		FullMatrix(const size_t n): FullMatrix(n,1) {}

		/*
			 Default constructor
		 */
		FullMatrix() =default;
		/*
			 Copy constructor
		 */
		FullMatrix(const FullMatrix &) = default;
		/*
			 Move constructor 
		 */
		FullMatrix(FullMatrix &&) = default;
		/*
			 Copy assignment
		 */
		FullMatrix &operator=(const FullMatrix &) = default;
		/*
			 Move assignment
		 */
		FullMatrix &operator=(FullMatrix &&) = default;

		/*
			 Constructor from an expression
			 NB:this is called also when creating a matrix from a transposed
		 */
		template <class T> FullMatrix(const Expr<T> &e): m_rows(e.rows()), m_cols(e.cols()) {
			//Casting
			const T &el(e);
			m_entries.resize(m_rows*m_cols);
			if constexpr(ORDER==ORDERING::ROWMAJOR){
#pragma omp parallel for shared(m_entries,m_rows,m_cols,el) 
				for(size_t i=0;i<m_rows;++i){
					const auto ofs=i*m_cols;
					for(size_t j=0;j<m_cols;++j)
						m_entries[ofs+j]=el(i,j);
				}
			}
			else{
#pragma omp parallel for shared(m_entries,m_rows,m_cols,el) 
				for(size_t j=0;j<m_cols;++j){
					const auto ofs=j*m_rows;
					for(size_t i=0;i<m_rows;++i)
						m_entries[ofs+i]=el(i,j);
				}
			}
		}
		/*
			 Assignment from an expression
		 */
		template <class T> FullMatrix& operator=(const Expr<T> &e){
			//Casting
			const T &el(e);
			this->resize(e.rows(),e.cols());
			if constexpr(ORDER==ORDERING::ROWMAJOR){
#pragma omp parallel for shared(m_entries,m_rows,m_cols,el) 
				for(size_t i=0;i<m_rows;++i){
					const auto ofs=i*m_cols;
					for(size_t j=0;j<m_cols;++j)
						m_entries[ofs+j]=el(i,j);
				}
			}
			else{
#pragma omp parallel for shared(m_entries,m_rows,m_cols,el) 
				for(size_t j=0;j<m_cols;++j){
					const auto ofs=j*m_rows;
					for(size_t i=0;i<m_rows;++i)
						m_entries[ofs+i]=el(i,j);
				}
			}
			return *this;
		}

		/*
			 When i found a matrix i collapse it by using the constructor of the FullMatrix class
		 */
		template <class T>
			const Real normExpr(const Expr<T> &e) const {
				FullMatrix toReturn(e);
				return toReturn.norm();
			}

		/*
			 Override of the operator []
		 */
		Real operator[](const size_t index){
			return m_entries[index];
		}

		/*
			 Override of the operator [] for read only operations
		 */
		const Real operator[](const size_t index) const{
			return m_entries[index];
		}

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
				return m_entries[i+m_rows*j];
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
			return this->operator[](i);
		}
		/*
			 Another operator for accessing an element in read only operations
		 */
		const Real& coeffRef(const size_t i) const{
			return this->operator[](i);
		}

		/*
			 This class stores a reference to a FullMatrix 
			 which is either stored in rowmajor or colmajor, then it is accessed in the opposite way.

			 It also need to enable expression templating
		 */

		class FullMatrixAdjoint: public Expr<FullMatrixAdjoint>{

			public:

				/*
					 Constructor of the class
				 */
				FullMatrixAdjoint(const FullMatrix<Real,ORDER>& matrix) : m_matrix(matrix) {}

				/*
					 The only thing i need is an overload to the access operator in the opposite way
				 */
				Real operator[](const size_t index){
					if constexpr(ORDER==ORDERING::ROWMAJOR)
						return m_matrix(index%cols(),index/cols());
					else
						return m_matrix(index/rows(),index%rows());
				}
				const Real operator[](const size_t index) const{
					if constexpr(ORDER==ORDERING::ROWMAJOR)
						return m_matrix(index%cols(),index/cols());
					else
						return m_matrix(index/rows(),index%rows());
					//return m_matrix[index];
				}
				//This is the only important
				const Real& operator()(const size_t i, const size_t j) const{
					if constexpr(ORDER==ORDERING::ROWMAJOR)
						return m_matrix.m_entries[i+rows()*j];
					else
						return m_matrix.m_entries[i*cols()+j];
				}
				const Real& coeffRef(const size_t i, const size_t j) const{
					return this->operator()(i,j);
				}
				const Real& coeffRef(const size_t i) const{
					return this->operator[](i);
				}
				/*
					 Printing methods
				 */
				void print(std::ostream &os=std::cout) const{
					os<<rows()<<" "<<cols()<<std::endl;
					for(size_t i=0;i<rows();++i){
						for(size_t j=0;j<cols();++j){
							os<<this->operator()(i,j)<<" ";
						}
						os<<std::endl;
					}
				}
				friend std::ostream& operator<<(std::ostream& os, FullMatrixAdjoint const & mat){
					os<<mat.rows()<<" "<<mat.cols()<<std::endl;
					for(size_t i=0;i<mat.rows();++i){
						for(size_t j=0;j<mat.cols();++j){
							os<<mat(i,j)<<" ";
						}
						os<<std::endl;
					}

					return os;
				}
				/*
					 Some other useful methods
				 */
				const size_t& rows() const{
					return m_matrix.m_cols;
				}
				const size_t& cols() const{
					return m_matrix.m_rows;
				}

			private:

				/*
					 Reference to the matrix
				 */
				const FullMatrix<Real,ORDER>& m_matrix;

		};

		/*
			 This class is declared as friend so that it can also access the matrix private
			 attributes directly.
		 */
		friend class FullMatrixAdjoint;

#ifndef LAZY
		/*
			 Override of the multiplication for matrix-scalar
			 Note that it returns a new matrix
		 */
		friend FullMatrix operator*(const FullMatrix& A, const Real k){
			FullMatrix toReturn(A.rows(),A.cols());

#pragma omp parallel for shared(A)
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

#pragma omp parallel for shared(A)
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

#pragma omp parallel for shared(A)
			for(size_t i=0;i<A.m_entries.size();++i)
				toReturn.m_entries[i]=A.m_entries[i]-B.m_entries[i];

			return toReturn;
		}
#endif
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

			//https://stackoverflow.com/questions/2324658/how-to-determine-the-version-of-the-c-standard-used-by-the-compiler
			//#if __cplusplus == 202002L

			//			std::cout<<"HELLO BABY!!!"<<std::endl;

			//#else

			if constexpr(ORDER==ORDERING::ROWMAJOR){
#pragma omp parallel for shared(B,toReturn)
				for(size_t i=0;i<A.rows();++i){
					for(size_t k=0;k<B.rows();++k){
						const auto ofs1=i*B.cols();
						const auto ofs2=i*A.cols();
						const auto ofs3=k*B.cols();
						//#pragma omp parallel for shared(ofs1,ofs2,ofs3,k,B,toReturn)
						for(size_t j=0;j<B.cols();++j){
							toReturn.m_entries[ofs1+j]+=A.m_entries[ofs2+k]*B.m_entries[ofs3+j];
						}
					}
				}
			}
			else{
#pragma omp parallel for shared(B,toReturn)
				for(size_t j=0;j<B.cols();++j){
					for(size_t k=0;k<A.cols();++k){
						const auto ofs1=A.rows()*k;
						const auto ofs2=B.rows()*j;
						//#pragma omp parallel for shared(ofs1,ofs2,B,toReturn) 
						for(size_t i=0;i<A.rows();++i){
							toReturn.m_entries[i+ofs2]+=A.m_entries[i+ofs1]*B.m_entries[k+ofs2];
						}
					}
				}
			}

			//#endif

			return toReturn;
		}
		/*
			 Here i am missing also the case when i am doing Expression*Matrix and vice-versa

			 NB: i have to do the same for the transposed

Careful: if i do not specify the overload for the multiplication between a matrix and a 
transposed it enters here, hence creating a temporary for nothing
		 */
		template <is_expr T>
			friend FullMatrix operator*(const T &e, const FullMatrix& B){
				//I need to collapse the expression
				FullMatrix A=e;
				return A*B;
			}
		template <is_expr T>
			friend FullMatrix operator*(const FullMatrix& A, const T &e){
				//I need to collapse the expression
				FullMatrix B=e;
				return A*B;
			}
		template <is_expr T>
			friend FullMatrix operator*(const T &e, const FullMatrixAdjoint& B){
				//I need to collapse the expression
				FullMatrix A=e;
				return A*B;
			}
		template <is_expr T>
			friend FullMatrix operator*(const FullMatrixAdjoint& A, const T &e){
				//I need to collapse the expression
				FullMatrix B=e;
				return A*B;
			}
		/*
			 Also for matrix and adjoint 
TODO: think more about this
		 */
		friend FullMatrix operator*(const FullMatrixAdjoint& A, const FullMatrix& B){

			FullMatrix toReturn;

			if(A.cols()!=B.rows())
				return toReturn;

			toReturn.resize(A.rows(),B.cols());

			if constexpr(ORDER==ORDERING::ROWMAJOR){
#pragma omp parallel for shared(B,toReturn)
				for(size_t i=0;i<A.rows();++i){
					for(size_t k=0;k<B.rows();++k){
						const auto ofs1=i*B.cols();
						const auto ofs2=i*A.cols();
						const auto ofs3=k*B.cols();
						for(size_t j=0;j<B.cols();++j){
							toReturn.m_entries[ofs1+j]+=A[ofs2+k]*B.m_entries[ofs3+j];
						}
					}
				}
			}
			else{
#pragma omp parallel for shared(B,toReturn)
				for(size_t j=0;j<B.cols();++j){
					for(size_t k=0;k<A.cols();++k){
						const auto ofs1=A.rows()*k;
						const auto ofs2=B.rows()*j;
						for(size_t i=0;i<A.rows();++i){
							toReturn.m_entries[i+ofs2]+=A[i+ofs1]*B.m_entries[k+ofs2];
						}
					}
				}
			}

			//#endif

			return toReturn;
		}
		friend FullMatrix operator*(const FullMatrix& A, const FullMatrixAdjoint& B){

			FullMatrix toReturn;

			if(A.cols()!=B.rows())
				return toReturn;

			toReturn.resize(A.rows(),B.cols());

			if constexpr(ORDER==ORDERING::ROWMAJOR){
#pragma omp parallel for shared(B,toReturn)
				for(size_t i=0;i<A.rows();++i){
					for(size_t k=0;k<B.rows();++k){
						const auto ofs1=i*B.cols();
						const auto ofs2=i*A.cols();
						const auto ofs3=k*B.cols();
						for(size_t j=0;j<B.cols();++j){
							toReturn.m_entries[ofs1+j]+=A.m_entries[ofs2+k]*B[ofs3+j];
						}
					}
				}
			}
			else{
#pragma omp parallel for shared(B,toReturn)
				for(size_t j=0;j<B.cols();++j){
					for(size_t k=0;k<A.cols();++k){
						const auto ofs1=A.rows()*k;
						const auto ofs2=B.rows()*j;
						for(size_t i=0;i<A.rows();++i){
							toReturn.m_entries[i+ofs2]+=A.m_entries[i+ofs1]*B[k+ofs2];
						}
					}
				}
			}

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
			this->resize(n,m,0.);
		}
		/*
			 Resize of the matrix with a parameter also
		 */
		void resize(const size_t n, const size_t m, const Real val){
			m_entries.resize(n*m,val);
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
			Method to return a block of the matrix starting from position i,j and of dimension nxm
		*/
		FullMatrix block(const size_t i, const size_t j, const size_t n, const size_t m) const{
			
			FullMatrix toReturn;
		
			//Do some checks on the bound
			if(i+n>this->m_rows || j+m>this->m_cols)
				return toReturn;

			toReturn.resize(n,m,0.);

			for(size_t h=0;h<n;++h)
				for(size_t k=0;k<m;++k)
					toReturn(h,k)=this->operator()(i+h,j+k);
			
			return toReturn;
		}
		/*
			 Method to multiply the current matrix by a constant
		 */
		void scale(const Real k){
			for(size_t i=0;i<m_entries.size();++i)
				m_entries[i]*=k;
		}

		/*
			 Method for retrieving the Frobenius norm
		 */
		const Real norm() const {
			Real norm=0.;
			for(auto el: m_entries)
				norm+=el*el;

			return std::sqrt(norm);
		}
		/*
			 Method to return the transposed of a matrix

			 Now it returns a FullMatrixAdjoint object
		 */
		FullMatrixAdjoint transpose() const{

			FullMatrixAdjoint toReturn(*this);

			return toReturn;
			/*
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
			*/
		}
		/*
			 Method for printing all elements in the matrix in a generic stream
			 (same as the professor did in one of the labs)
		 */
		void print(std::ostream &os=std::cout) const{
			os<<m_rows<<" "<<m_cols<<std::endl;
			for(size_t i=0;i<m_rows;++i){
				for(size_t j=0;j<m_cols;++j){
					os<<this->operator()(i,j)<<" ";
				}
				os<<std::endl;
			}
		}
		/*
			 Another method for printing a matrix
		 */
		friend std::ostream& operator<<(std::ostream& os, FullMatrix const & mat){
			os<<mat.rows()<<" "<<mat.cols()<<std::endl;
			for(size_t i=0;i<mat.m_rows;++i){
				for(size_t j=0;j<mat.m_cols;++j){
          os<<std::setw(8)<<std::fixed<<std::setprecision(4)<<mat(i,j)<<" ";
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

		/*
			Method to return the data of the matrix

			Useful for MPI
		*/
		Real* data(){
			return this->m_entries.data();
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


//EXPRESSION TEMPLATES

#ifdef LAZY

/*
	 Using some useful aliases to alleviate the notation for user functions
 */
template<class L, class R> using AddExpr=CWiseBinaryOperator<L,R,Add>;
template<class L, class R> using SubExpr=CWiseBinaryOperator<L,R,Sub>;
template<class L, class R> using MultExpr=CWiseBinaryOperator<L,R,Mult>;

/*
	 These are the operators that the user will use
	 They are defined inline so that they can be substituted at compile time

	 Moreover they are specialized only for my class 
 */

template<is_expr L, is_expr R>
inline AddExpr<L,R>
operator+(const L &l, const R &r){
	return AddExpr<L,R>(l,r);
}

template<is_expr L, is_expr R>
inline SubExpr<L,R>
operator-(const L &l, const R &r){
	return SubExpr<L,R>(l,r);
}

template<is_expr L, std::floating_point R>
inline MultExpr<L,R> 
operator*(const L &l, const R &r){
	return MultExpr<L,R>(l,r);
}

template<std::floating_point L, is_expr R>
inline MultExpr<L,R> 
operator*(const L &l, const R &r){
	return MultExpr<L,R>(l,r);
}

#endif

#pragma GCC diagnostic pop

#endif
