#ifndef INCLUDED_ObjexxFCL_DimensionExpressionCon_hh
#define INCLUDED_ObjexxFCL_DimensionExpressionCon_hh


// DimensionExpressionCon: Constant-Valued DimensionExpression
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/DimensionExpression.hh>


namespace ObjexxFCL {


/// @brief DimensionExpressionCon: Constant-Valued DimensionExpression
class DimensionExpressionCon :
	public DimensionExpression
{


private: // Types


	typedef  DimensionExpression  Super;


public: // Creation


	/// @brief Copy Constructor
	inline
	DimensionExpressionCon( DimensionExpressionCon const & exp ) :
		Super(),
		value_( exp.value_ ),
		integer_( exp.integer_ )
	{}


	/// @brief int Constructor
	inline
	explicit
	DimensionExpressionCon( int const value_a ) :
		value_( static_cast< double >( value_a ) ),
		integer_( true )
	{}


	/// @brief double Constructor
	inline
	explicit
	DimensionExpressionCon( double const value_a ) :
		value_( value_a ),
		integer_( false )
	{}


	/// @brief Clone
	inline
	DimensionExpressionCon *
	clone() const override
	{
		return new DimensionExpressionCon( *this );
	}


	/// @brief Clone with Dimension Substitution
	inline
	DimensionExpressionCon *
	clone( Dimension const & ) const override
	{
		return new DimensionExpressionCon( *this );
	}


	/// @brief Destructor
	inline
	~DimensionExpressionCon() override
	{}


public: // Inspector


	/// @brief Initialized?
	inline
	bool
	initialized() const override
	{
		return true;
	}


	/// @brief Integer?
	inline
	bool
	integer() const override
	{
		return integer_;
	}


	/// @brief Constant?
	inline
	bool
	constant() const override
	{
		return true;
	}


	/// @brief Reference?
	inline
	bool
	reference() const override
	{
		return false;
	}


	/// @brief Reducible?
	inline
	bool
	reducible() const override
	{
		return false;
	}


	/// @brief Value
	inline
	double
	operator ()() const override
	{
		return value_;
	}


	/// @brief Value
	inline
	double
	value() const override
	{
		return value_;
	}


	/// @brief Insert an Observer
	inline
	void
	insert_observer( Observer & ) const override
	{}


	/// @brief Remove an Observer
	inline
	void
	remove_observer( Observer & ) const override
	{}


public: // Modifier


	/// @brief Update for Destruction of a Subject
	inline
	void
	destructed( Subject const & ) override
	{}


private: // Data


	/// @brief Value
	double value_;

	/// @brief Integer-valued?
	bool integer_;


}; // DimensionExpressionCon


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_DimensionExpressionCon_HH
