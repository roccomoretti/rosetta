// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/KofNConstraint.hh.cxxtest.hh
/// @brief  test suite for KofNConstraint
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/scoring/constraints/KofNConstraint.hh> // DO NOT AUTO-REMOVE
// basic headers

// unit headers

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class KofNConstraintTests : public CxxTest::TestSuite
{
public:


	void test_serialize_KofNConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		ConstraintCOPs csts;
		csts.push_back( utility::pointer::make_shared< AtomPairConstraint >( AtomID( 1, 2 ), AtomID( 2, 3 ), some_func ));
		csts.push_back( utility::pointer::make_shared< AtomPairConstraint >( AtomID( 3, 4 ), AtomID( 4, 5 ), some_func ));
		csts.push_back( utility::pointer::make_shared< AtomPairConstraint >( AtomID( 4, 5 ), AtomID( 5, 6 ), some_func ));
		csts.push_back( utility::pointer::make_shared< AtomPairConstraint >( AtomID( 6, 7 ), AtomID( 7, 8 ), some_func ));
		ConstraintOP instance( new KofNConstraint( csts, 2 ) ); // serialize this through a pointer to the base class

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		ConstraintOP instance2; // deserialize also through a pointer to the base class
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		// make sure the deserialized base class pointer points to a KofNConstraint
		TS_ASSERT( utility::pointer::dynamic_pointer_cast< KofNConstraint > ( instance2 ));
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

