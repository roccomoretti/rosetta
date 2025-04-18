// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    mp_mutate_relax.cc
/// @brief   Mutate a residue, then do range relax for a membrane protein
/// @author  JKLeman (julia.koehler1982@gmail.com)

// Unit Headers
#include <protocols/relax/membrane/MPMutateRelaxMover.hh>
#include <protocols/relax/membrane/MPMutateRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/kinematics/FoldTree.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/relax/membrane/MPRangeRelaxMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/chemical/AA.hh>

// Package Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.hh>
#include <protocols/membrane/util.hh>
#include <utility/io/util.hh>
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <core/pose/util.hh>

// C++ Headers

#include <protocols/relax/membrane/util.hh> // AUTO IWYU For add_mutant_to_vectors, check_mutants_ok
#include <core/pack/task/PackerTask.hh> // AUTO IWYU For PackerTask
#include <utility/stream_util.hh> // AUTO IWYU For operator<<

static basic::Tracer TR( "protocols.membrane.MPMutateRelaxMover" );

namespace protocols {
namespace membrane {

using core::Size;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Mutations will be read in by commandline
MPMutateRelaxMover::MPMutateRelaxMover() : protocols::moves::Mover()
{
	register_options();
	init_from_cmd();
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
MPMutateRelaxMover::MPMutateRelaxMover( MPMutateRelaxMover const & ) = default;

/// @brief Assignment Operator
MPMutateRelaxMover & MPMutateRelaxMover::operator = ( MPMutateRelaxMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new MPMutateRelaxMover( *this ) ); // ??????????????????????
}

/// @brief Destructor
MPMutateRelaxMover::~MPMutateRelaxMover() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPMutateRelaxMover::clone() const {
	return ( utility::pointer::make_shared< MPMutateRelaxMover >( *this ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPMutateRelaxMover::fresh_instance() const {
	return utility::pointer::make_shared< MPMutateRelaxMover >();
}

/// @brief Parse Rosetta Scripts Options for this Mover
void
MPMutateRelaxMover::parse_my_tag(
	utility::tag::TagCOP /*tag*/,
	basic::datacache::DataMap &
) {

	// TODO: implement this

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
MPMutateRelaxMoverCreator::create_mover() const {
	return utility::pointer::make_shared< MPMutateRelaxMover >();
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MPMutateRelaxMoverCreator::keyname() const {
	return MPMutateRelaxMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MPMutateRelaxMoverCreator::mover_name() {
	return "MPMutateRelaxMover";
}

/// @brief Get the name of this Mover (MPMutateRelaxMover)
std::string
MPMutateRelaxMover::get_name() const {
	return "MPMutateRelaxMover";
}

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Mutate residue and then range relax the membrane protein
void MPMutateRelaxMover::apply( core::pose::Pose & pose ) {

	using namespace basic::options;
	using namespace utility;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::membrane;
	using namespace protocols::relax::membrane;
	using namespace protocols::simple_moves;
	using namespace core::scoring;
	using namespace core::pack::task;

	// finalize setup
	finalize_setup( pose );

	TR << "Running MPMutateRelax protocol with SCOREFUNCTION " << utility::CSI_Green() << sfxn_->get_name() << utility::CSI_Reset() << std::endl;

	// final foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// WHAT IS HAPPENING HERE:
	// I am dumping the poses outside of JD2
	// Why would I do that?
	// because I want the mutant to be part of the output filename

	// make a deepcopy of the pose to start from
	Pose original_pose = Pose( pose );
	Pose working_pose;

	// construct counter
	core::Size counter(0);

	// for debugging, iterating over constructs and mutants
	for ( core::Size c = 1; c <= wt_res_.size(); ++c ) {
		TR.Debug << "construct " << c << std::endl;

		for ( core::Size m = 1; m <= wt_res_[c].size(); ++m ) {
			TR.Debug << "mutant " << m << " wt_res " << wt_res_[c][m] << " resid " << resn_[c][m] << " mut_res " << mut_res_[c][m] << std::endl;
		}
	}

	// go through each construct (i.e. outer vector)
	for ( core::Size c = 1; c <= wt_res_.size(); ++c ) {

		TR << "going through construct " << c << std::endl;

		// counter
		counter = 1;
		//std::string mutations;

		// iterate over nstruct
		while ( counter <= nstruct_ ) {

			TR << "working on nstruct " << counter << std::endl;

			// get original starting pose
			working_pose = Pose( original_pose );

			// make mutations
			std::string mutations = make_mutations( working_pose, c );

			// if repacking
			if ( repack_mutation_only_ == true || repack_radius_ > 0 ) {

				// create task factory and repack
				TR << "Repacking only..." << std::endl;
				PackerTaskOP repack = TaskFactory::create_packer_task( working_pose );
				TR.Debug << "repacking vector size: " << repack_residues_[c].size() << std::endl;
				TR.Debug << "pose length: " << working_pose.size() << std::endl;
				repack->restrict_to_residues( repack_residues_[c] );
				repack->restrict_to_repacking();
				core::pack::pack_rotamers( working_pose, *sfxn_, repack );

				// if relaxing
			} else if ( relax_ == true ) {

				// do range relax
				TR << "Running MPRangeRelax..." << std::endl;
				MPRangeRelaxMoverOP relax( new MPRangeRelaxMover() );
				//    relax->optimize_membrane( false );
				relax->apply( working_pose );
			}

			// create output filename
			std::string output;

			// if model exists, increment counter
			// this means that the app should be able to start from an existing file number
			core::Size a = counter;
			while ( a <= nstruct_ ) {
				output = output_filename( mutations, a );
				if ( utility::file::file_exists( output ) ) {
					++a;
				} else {
					break;
				}
			}

			// dump pose
			output = output_filename( mutations, a );
			working_pose.dump_scored_pdb( output, *sfxn_ );

			// increment counter
			++counter;

		}// nstruct
	}// iterate over construct

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void MPMutateRelaxMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::in::file::s );
	option.add_relevant( OptionKeys::mp::mutate_relax::mutation );
	option.add_relevant( OptionKeys::mp::mutate_relax::mutant_file );
	option.add_relevant( OptionKeys::mp::mutate_relax::nmodels );
	option.add_relevant( OptionKeys::mp::mutate_relax::repack_mutation_only );
	option.add_relevant( OptionKeys::mp::mutate_relax::repack_radius );
	option.add_relevant( OptionKeys::mp::mutate_relax::relax );

}

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize from commandline
void MPMutateRelaxMover::init_from_cmd() {

	using namespace basic::options;

	// making sure models won't be overwritten
	if ( option[ OptionKeys::out::nstruct ]() > 1 ) {
		utility_exit_with_message( "Set '-out:nstruct 1'! Otherwise models will be overwritten!" );
	}

	// input checking
	if ( ! option[ OptionKeys::mp::mutate_relax::mutant_file ].user() &&
			! option[ OptionKeys::mp::mutate_relax::mutation ].user() ) {

		utility_exit_with_message( "Too few inputs: You must EITHER specify the mutant file with -mp:mutate_relax:mutant_file OR specify a single mutation with -mp:mutate_relax:mutation! Quitting..." );
	}

	// more input checking
	if ( option[ OptionKeys::mp::mutate_relax::mutant_file ].user() &&
			option[ OptionKeys::mp::mutate_relax::mutation ].user() ) {

		utility_exit_with_message( "Too many inputs: You must EITHER specify the mutant file with -mp:mutate_relax:mutant_file OR specify a single mutation with -mp:mutate_relax:mutation! Quitting..." );
	}

	// get mutant file
	if ( option[ OptionKeys::mp::mutate_relax::mutant_file ].user() ) {
		mutant_file_ = option[ OptionKeys::mp::mutate_relax::mutant_file ]();
	}

	// Number of iterations to run
	nstruct_ = 1;
	if ( option[ OptionKeys::mp::mutate_relax::nmodels ].user() ) {
		nstruct_ = option[ OptionKeys::mp::mutate_relax::nmodels ]();
	}

	// get protein name for dumping PDBs
	if ( option[ OptionKeys::in::file::s ].user() ) {
		protein_ = option[ OptionKeys::in::file::s ](1);
	} else {
		utility_exit_with_message("No PDB given, please use -in:file:s to provide PDB.");
	}

	// repack options, overwrite iterations to 1
	repack_mutation_only_ = false;
	if ( option[ OptionKeys::mp::mutate_relax::repack_mutation_only ].user() ) {
		repack_mutation_only_ = option[ OptionKeys::mp::mutate_relax::repack_mutation_only ]();
		nstruct_ = 1;
		TR << "Repacking only: setting number of iterations to 1." << std::endl;
	}

	// repacking within radius, overwrite iterations to 1
	repack_radius_ = 0;
	if ( option[ OptionKeys::mp::mutate_relax::repack_radius ].user() ) {
		repack_radius_ = option[ OptionKeys::mp::mutate_relax::repack_radius ]();
		nstruct_ = 1;
		TR << "Repacking only: setting number of iterations to 1." << std::endl;
	}

	// only running relax
	relax_ = false;
	if ( option[ OptionKeys::mp::mutate_relax::relax ].user() ) {
		relax_ = option[ OptionKeys::mp::mutate_relax::relax ]();
		TR << "Setting relax option to " << relax_ << std::endl;
	}

	// check number of iterations
	if ( relax_ == true && nstruct_ == 0 ) {
		nstruct_ = 100;
		TR << "Relaxing structures: setting number of iterations to 100." << std::endl;
	} else if ( relax_ == true ) {
		TR << "Relaxing structures: number of iterations is " << nstruct_ << " as per user-input." << std::endl;
		TR << "I assume you know what you are doing, a good number of iterations is 100." << std::endl;
	} else if ( relax_ == false && repack_mutation_only_ == false && repack_radius_ == 0 ) {
		nstruct_ = 1;
		TR << "Neither repacking nor relaxing structures: setting number of iterations to 1." << std::endl;
	}

	// checking inputs
	if ( repack_mutation_only_ == true && repack_radius_ > 0 ) {
		utility_exit_with_message("Set EITHER repack_mutation_only OR repack_radius, not both!");
	}
	if ( ( repack_mutation_only_ == true || repack_radius_ > 0 ) && relax_ == true ) {
		utility_exit_with_message("Set EITHER repacking option OR relax option, not both!");
	}

}// init from cmdline

////////////////////////////////////////////////////////////////////////////////

/// @brief Get repack residues
void MPMutateRelaxMover::get_repack_residues( Pose & pose ) {

	using namespace core::pose;

	// go through number of constructs
	for ( core::Size i = 1; i <= resn_.size(); ++i ) {

		// initialize boolean vector with false
		utility::vector1< bool > repack_res( pose.size(), false );

		// if repacking
		if ( repack_mutation_only_ == true ) {

			// go through number of mutations within construct
			for ( core::Size j = 1; j <= resn_[i].size(); ++j ) {

				// set residue number of the mutation to true
				core::Size resn = resn_[i][j];
				repack_res[ resn ] = true;
			}
		} else if ( repack_radius_ > 0 ) {

			// go through number of mutations within construct
			for ( core::Size j = 1; j <= resn_[i].size(); ++j ) {

				core::Size resn = resn_[i][j];
				core::Vector xyz_mut = pose.residue( resn ).xyz( "CA" );

				// go through residues in pose
				for ( core::Size k = 1; k <= nres_protein( pose ); ++k ) {

					// check whether the residue is within repack radius
					core::Vector xyz_k = pose.residue( k ).xyz( "CA" );
					core::Real dist = ( xyz_k - xyz_mut ).length();

					// if yes, set repack flag of this residue to true
					if ( dist <= repack_radius_ ) {
						repack_res[ k ] = true;
					}
				}
			}
		}

		// add vector describing construct to total constructs
		repack_residues_.push_back( repack_res );
	}

} // get repack residues

////////////////////////////////////////////////////////////////////////////////

/// @brief Finalize setup
void MPMutateRelaxMover::finalize_setup( Pose & pose ){

	using namespace utility;
	using namespace utility::io;
	using namespace basic::options;
	using namespace protocols::relax::membrane;

	// one line per input file: several mutants within a single construct
	utility::vector1< std::string > mutations;

	// mutants given in command line
	// input format A163F in pose numbering, can give multiple mutations
	// split multiple mutations by whitespace
	if ( option[ OptionKeys::mp::mutate_relax::mutation ].user() ) {
		utility::vector1< std::string > muts = option[ OptionKeys::mp::mutate_relax::mutation ]();
		mutations.push_back( join( muts, " " ) );
	}
	// mutants given in mutant file
	if ( option[ OptionKeys::mp::mutate_relax::mutant_file ].user() ) {
		//  mutations = read_mutant_file_constructs( mutant_file_ );
		mutations = get_lines_from_file_data( mutant_file_ );
	}

	TR << "Looking at mutations " << mutations << std::endl;

	// loop over constructs
	for ( core::Size i = 1; i <= mutations.size(); ++i ) {

		// add mutants in that line to vectors
		add_mutant_to_vectors( pose, mutations[i], wt_res_, resn_, mut_res_ );

	}

	// error checking
	if ( ! check_mutants_ok( pose, wt_res_, resn_ ) ) {
		utility_exit_with_message( "Residue identity in input file doesn't match the pose!" );
	}

	// call AddMembraneMover
	if ( ! pose.conformation().is_membrane() ) {
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );
	}

	// get repack residues
	get_repack_residues( pose );

	// create scorefunction
	sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );

}// finalize setup

////////////////////////////////////////////////////////////////////////////////

/// @brief Make mutations, returns string of output file
std::string MPMutateRelaxMover::make_mutations( Pose & pose, core::Size num_construct ) {

	using namespace utility;
	using namespace core::chemical;
	using namespace protocols::simple_moves;
	using namespace protocols::relax::membrane;

	// mutations
	std::string mutations = "";
	core::Size c = num_construct;

	// go through each mutation for this construct
	for ( core::Size m = 1; m <= wt_res_[ c ].size(); ++m ) {

		// mutate residue
		TR << "Mutating residue " << wt_res_[c][m] << resn_[c][m] << " to " << mut_res_[c][m] << std::endl;
		std::string aa3 = name_from_aa( aa_from_oneletter_code( mut_res_[c][m] ) );
		MutateResidueOP mutate( new MutateResidue( resn_[c][m], aa3 ) );
		mutate->apply( pose );

		// get mutation as output tag
		mutations += "_" + to_string( wt_res_[c][m] ) + to_string( resn_[c][m] ) + to_string( mut_res_[c][m] );

	}// iterate over mutations in construct

	return mutations;
} // make mutations

////////////////////////////////////////////////////////////////////////////////

/// @brief Create output filename
std::string MPMutateRelaxMover::output_filename( std::string mutation_tag, core::Size counter ) {

	using namespace basic::options;
	using namespace utility;

	// create output filename to /my/path/myprotein.pdb to myprotein
	const std::string tmp( file_basename( protein_ ) );
	std::string output;

	// get pathname
	if ( option[ OptionKeys::out::path::pdb ].user() ) {
		output = option[ OptionKeys::out::path::pdb ]().path() + trim( tmp, ".pdb");
	} else {
		output = trim( tmp, ".pdb");
	}

	// add mutation to output filename
	output += mutation_tag;

	// add counter to output filename
	std::string cnt;
	if ( counter < 10 ) {
		cnt = "000" + to_string( counter );
	} else if ( counter > 9 && counter < 100 ) {
		cnt = "00" + to_string( counter );
	} else if ( counter > 99 && counter < 1000 ) {
		cnt = "0" + to_string( counter );
	} else if ( counter > 999 && counter < 10000 ) {
		cnt = to_string( counter );
	} else {
		utility_exit_with_message( "Please choose an -nstruct < 9999." );
	}
	output += "_" + cnt + ".pdb";

	return output;

}
} // membrane
} // protocols

