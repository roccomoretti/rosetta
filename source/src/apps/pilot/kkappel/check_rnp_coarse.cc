// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/kkappel/calculate_distances.cc
/// @brief see what's going on with protein centroid with RNA FA
/// @author Kalli Kappel kappel@stanford.edu


//Unit Headers
//Package Headers
//Project Headers
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/util/SwitchResidueTypeSet.hh>
//Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
//Numeric Headers
//C++ Headers
#include <iostream>

#include <devel/init.hh>



// New options for this application
using namespace basic::options::OptionKeys;

//OPT_KEY( String, mutfile )

static basic::Tracer TR( "apps.pilot.kkappel.check_RNP_coarse" );

void check_structures() {

	// Make a pose from a pdb file
	core::pose::Pose pose;
	//core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	//pose = *(protocols::stepwise::setup::initialize_pose_and_other_poses_from_command_line( rsd_set ));
	utility::vector1< std::string > input_structs = basic::options::option[ basic::options::OptionKeys::in::file::s ]();
	core::import_pose::pose_from_file( pose, input_structs[ 1 ] );
	// If I want to try to import it as a centroid pose
	//core::import_pose::centroid_pose_from_pdb( pose, input_structs[ 1 ] );
	pose.dump_pdb("input_pose.pdb");

	// Ok now what if I want to change one residue to be coarse grained

	// Get the centroid residue type set
	core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
	//core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID_ROT );

	core::scoring::ScoreFunctionOP sfxn;
	sfxn = core::scoring::get_score_function();
	sfxn->show( pose );

	// Make a new residue
	core::conformation::Residue rsd = pose.residue( 70 ); // totally random for testing
	core::conformation::ResidueOP new_rsd( nullptr );
	core::chemical::ResidueType const & new_rsd_type( rsd_set->name_map(rsd.name()) );
	new_rsd = core::conformation::ResidueFactory::create_residue( new_rsd_type, rsd, pose.conformation() );

	// Now try to replace it
	pose.replace_residue(70, *new_rsd, false );

	// Is the pose ok with this?
	pose.dump_pdb("after_replace_residue.pdb");

	// Try scoring
	core::scoring::ScoreFunctionOP sfxn_coarse;
	sfxn_coarse = core::scoring::ScoreFunctionFactory::create_score_function( "cen_std.wts" );
	// what about scoring with rna_lores?
	core::scoring::ScoreFunctionOP sfxn_rna;
	sfxn_rna = core::scoring::ScoreFunctionFactory::create_score_function( "rna/denovo/rna_lores.wts" );

	std::cout << "trying to score with rna_lores" << std::endl;

	sfxn->show( pose );
	sfxn_rna->show( pose );
	// If you try to score with the coarse score function, it seg faults!
	//sfxn_coarse->show( pose );
	// what if I totally switch residue type sets first?
	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
	std::cout << "Just switched to centroid mode" << std::endl;
	sfxn_coarse->show( pose );
	pose.dump_pdb("centroid_pose.pdb");
	sfxn_rna->show( pose );
	// RNA looks just the same after switching to the CENTROID residue type set
	// Is the switch_to_residue_type_set actually putting the "CEN" atom of the
	// centroid residues at the correct place??
	// Yes, b/c this is defined for each residue type by averaging over observed
	// side-chain conformations in known protein structures

}

int main( int argc, char ** argv ) {

	try {
		using namespace basic::options;
		//NEW_OPT( mutfile, "File listing mutations", "test_mutfile" );

		devel::init( argc, argv );
		check_structures();
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
