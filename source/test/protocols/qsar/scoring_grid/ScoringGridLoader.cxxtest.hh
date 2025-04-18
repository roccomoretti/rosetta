// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/ScoringGridLoader.cxxtest.hh
/// @brief  test suite for protocols::moves::ScoringGridLoader
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package headers
#include <protocols/qsar/scoring_grid/ScoringGridLoader.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.hh>

// Test utility headers
#include <test/util/schema_utilities.hh>

using namespace protocols::qsar::scoring_grid;
using namespace utility::tag;

class ScoringGridLoaderTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_scoring_grid_loader_defines_valid_xsd() {
		XMLSchemaDefinition xsd;
		ScoringGridLoader::provide_xml_schema( xsd );
		ensure_rosetta_scripts_like_XSD_validates_w_ct( xsd, ScoringGridLoader::scoring_grid_loader_ct_namer( ScoringGridLoader::loader_name() ) );
	}

};

