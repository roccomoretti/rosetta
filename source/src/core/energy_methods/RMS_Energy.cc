// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RMS_Energy.cc
/// @brief  RMS Energy function. Used to optimize the RMSD between two structures.
/// @author James Thompson


// Unit headers
#include <core/energy_methods/RMS_Energy.hh>
#include <core/energy_methods/RMS_EnergyCreator.hh>

// Package headers

#include <core/scoring/rms_util.hh>
#include <basic/options/option.hh>
#include <utility/exit.hh>

// option key includes
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/io/pdb/build_pose_as_is.hh>
#include <core/scoring/EnergyMap.hh>


namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the RMS_Energy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
RMS_EnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< RMS_Energy >();
}

core::scoring::ScoreTypes
RMS_EnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( rms );
	return sts;
}

/// c-tor
RMS_Energy::RMS_Energy() :
	parent( utility::pointer::make_shared< RMS_EnergyCreator >() )
{
	// configure native pose. Not perfect, ideally we should have some way of
	// guaranteeing that the native pose and pose provided for scoring have
	// the same ResidueTypeSet.
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		//core::import_pose::pose_from_file( native_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		core::io::pdb::build_pose_from_pdb_as_is(native_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ]() );
	} else {
		utility_exit_with_message( "Error: must provide native pose when scoring with RMS_Energy!\n" );
	}

	// configure target RMSD.
	rms_target_ = basic::options::option[ basic::options::OptionKeys::score::rms_target ]();
}


/// clone
core::scoring::methods::EnergyMethodOP
RMS_Energy::clone() const
{
	return utility::pointer::make_shared< RMS_Energy >();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// @brief Calculate the RMS difference between native_pose_ (provided by
/// the option -in::file::native and the given Pose. The actual energy calculation
/// is the difference between the RMSD and the target RMSD. Target RMSD is specified
/// the option -score::rms_target.

void
RMS_Energy::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals
) const {

	// PROF_START( basic::RMS );

	///////////////////////////////////////
	//
	// RMS SCORE
	core::Real rmsd;
	if ( basic::options::option[ basic::options::OptionKeys::evaluation::rms_type ].value() == "RNP" ) {
		rmsd = core::scoring::CA_or_equiv_rmsd( native_pose_, pose );
	} else {
		rmsd = core::scoring::CA_rmsd( native_pose_, pose );
	}
	totals[ core::scoring::rms ]  = std::abs( rms_target_ - rmsd );

	// PROF_STOP( basic::RMS );
}
core::Size
RMS_Energy::version() const
{
	return 1; // Initial versioning
}

} // namespace energy_methods
} // namespace core
