#ifndef HH__EXPRESSION_WRAPPER__HH
#define HH__EXPRESSION_WRAPPER__HH

#include "fullMatrix.hpp"
#include <functional>
#include <iostream>
#include <type_traits>
#include <utility>

using Real=double;

/*
	 Generic class that represents some expression 

	 Any class that incapsulates an expression should inherit from it.

	 It contains the (static) cast to the derived class (both const and not)
 */
template <class E> struct Expr{
	
	public:
	
		/*
			Note: the class should be never instantiated alone, so this contains a pointer to the derived class
		 */
		operator const E&() const{
			return static_cast<const E&>(*this);
		}	
		/*
			 Non const cast operator
		 */
		operator E&(){
			return static_cast<E&>(*this);
		}
		//Also as methods for more clarity
		constexpr const E& asDerived() const{
			return static_cast<const E&>(*this);
		}
		constexpr E& asDerived(){
			return static_cast<E&>(*this);
		}

		/*
			Operator to access an element directly (and constant version)
		*/
		Real operator()(const std::size_t i, const std::size_t j) const{
			return asDerived().operator()(i,j);
		}
		Real& operator()(const std::size_t i, const std::size_t j){
			return asDerived().operator()(i,j);
		}

		/*
			Methods to interrogate on the dimensions of the matrix
		*/
		std::size_t rows() const{
			return asDerived().rows();
		}	
		std::size_t cols() const{
			return asDerived().cols();
		}	

		/*
			Method to return the norm of the expression.
			It keeps track of the current expression so that it can be passed to a FullMatrix
		*/
		const Real norm() const{
			return asDerived().normExpr(*this);
		}

};

#include <concepts>
#include <type_traits>

/*
	Class that encapsulates component wise operations
*/
template<class L, class R, class OP>
class CWiseBinaryOperator: public Expr<CWiseBinaryOperator<L,R,OP>>{

	public:
		/*
			Create a constructor since attributes are const
		*/
		CWiseBinaryOperator(const L &l,const R &r): m_l(l), m_r(r) {};

		/*
			Apply the specified operator to the left and right classes
		*/
		Real operator()(const size_t i, const size_t j) const{
			return OP()(m_l(i,j), m_r(i,j));
		}

		std::size_t rows() const{
			return m_l.rows();
		}
		std::size_t cols() const{
			return m_r.cols();
		}

		/*	
			This class has two options, it either contains an Expression or a Matrix
			So i just need to call one of the two classes L or R recursively.

			Note that it does not matter which since i will arrive to a Matrix anyway
		*/
		template <class T>
		const Real normExpr(const Expr<T> &e) const {
			return m_l.normExpr(e);
		}

	private:
		/*
		Keep a reference to the current expressions
		*/
		const L &m_l;
		const R &m_r;

};

/*
	Specialization for scalar operations (right)

	NB: with scalar operations we can access the element also using (,) operator.
*/
template<class L, class OP>
class CWiseBinaryOperator<L,Real,OP>: public Expr<CWiseBinaryOperator<L,Real,OP>>{

	public:
		/*
			Create a constructor since attributes are const
		*/
		CWiseBinaryOperator(const L &l,const Real &r): m_l(l), m_r(r) {};

		/*
			Apply the specified operator to the left and right classes
		*/
		Real operator()(const size_t i, const size_t j) const{
			return OP()(m_l(i,j), m_r);
		}
		
		std::size_t rows() const{
			return m_l.rows();
		}
		std::size_t cols() const{
			return m_l.cols();
		}

		/*	
			Here i call the non-double operator
		*/
		template <class T>
		const Real normExpr(const Expr<T> &e) const {
			return m_l.normExpr(e);
		}

	private:
		/*
		Keep a reference to the current operators
		*/
		const L &m_l;
		const Real &m_r;

};

/*
	Another specialization (left)
*/
template< class R, class OP>
class CWiseBinaryOperator<Real,R,OP>: public Expr<CWiseBinaryOperator<Real,R,OP>>{

	public:
		/*
			Create a constructor since attributes are const
		*/
		CWiseBinaryOperator(const Real &l,const R &r): m_l(l), m_r(r) {};

		/*
			Apply the specified operator to the left and right classes
		*/
		Real operator()(const size_t i, const size_t j) const{
			return OP()(m_l, m_r(i,j));
		}
		
		std::size_t rows() const{
			return m_r.rows();
		}
		std::size_t cols() const{
			return m_r.cols();
		}
	
		/*	
			Here i call the non-double operator
		*/
		template <class T>
		const Real normExpr(const Expr<T> &e) const {
			return m_r.normExpr(e);
		}

	private:
		/*
		Keep a reference to the current operators
		*/
		const Real &m_l;
		const R &m_r;

};

/*
	Use std library multiplication for matrix-scalar
*/

using Mult=std::multiplies<Real>; 
using Add=std::plus<Real>;
using Sub=std::minus<Real>;

/*
	Create a concept for the expression
*/
template<class T>
concept is_expr=std::is_base_of<Expr<T>, T>::value;

#endif
