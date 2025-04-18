// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/methods/PairEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::PairEnergy.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers


// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>


//Auto Headers


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class OccludedHbondSolEnergyTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_occ_sol_fitted_start_score_start_func_match_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( occ_sol_fitted, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		//std::cout.precision( 16 );
		//std::cout << sfxn.score( pose ) << std::endl;
		adv.validate_start_func_matches_start_score( 29.77499608096777, false, 1e-6 );
	}

	void test_occ_sol_fitted_deriv_check_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( occ_sol_fitted, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( false, 1e-6 );
	}

	void test_occ_sol_fitted_start_score_start_func_match_w_partial_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( occ_sol_fitted, 0.5 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		//std::cout.precision( 16 );
		//std::cout << sfxn.score( pose ) << std::endl;
		adv.validate_start_func_matches_start_score( 29.77499608096777, false, 1e-6 );

	}

	void test_occ_sol_fitted_deriv_check_w_partial_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( occ_sol_fitted, 0.5 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.simple_deriv_check( false, 1e-6 );
	}

};


