// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief

//MPI headers typically have to go first
#ifdef USEMPI
#include <mpi.h>
#endif

// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>


#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>
#include <core/scoring/rms_util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>

#include <devel/init.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/FastRelax.fwd.hh>
#include <protocols/match/Hit.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <utility>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/excn/Exceptions.hh>

#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashSampler.hh>
#include <protocols/loophash/LocalInserter.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/loops/Loops.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

// C++ headers
//#include <cstdlib>

#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>


static basic::Tracer TR( "main" );

using namespace protocols::moves;
using namespace core::scoring;
using namespace core;
using namespace core::pose;
using namespace conformation;
using namespace kinematics;
using namespace protocols::match;
using namespace protocols::frag_picker;
using namespace protocols::loophash;
using namespace numeric::geometry::hashing;

class LoopHashRelax_Sampler;
using LoopHashRelax_SamplerOP = utility::pointer::shared_ptr<LoopHashRelax_Sampler>;
using LoopHashRelax_SamplerCOP = utility::pointer::shared_ptr<const LoopHashRelax_Sampler>;

class LoopHashRelax_Sampler: public protocols::moves::Mover {
public:

	LoopHashRelax_Sampler(
		LoopHashLibraryOP library
	):
		library_(std::move(library))

	{
	}

	void apply( core::pose::Pose& pose ) override;

	protocols::moves::MoverOP clone() const override {
		return utility::pointer::make_shared< LoopHashRelax_Sampler >( *this );
	}


	std::string get_name() const override {
		return "LoopHashRelax_Sampler";
	}

	protocols::moves::MoverOP fresh_instance() const override {
		return utility::pointer::make_shared< LoopHashRelax_Sampler >( library_ );
	}

private:
	LoopHashLibraryOP library_;

};

void
LoopHashRelax_Sampler::apply( core::pose::Pose& pose )
{
	if ( !library_ ) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string prefix = option[ out::prefix ]();
	core::Size skim_size = option[ lh::skim_size ]();   // default 100
	core::Real mpi_metropolis_temp_ = option[ lh::mpi_metropolis_temp ]();


	LocalInserter_SimpleMinOP simple_inserter( new LocalInserter_SimpleMin() );
	LoopHashSampler  lsampler( library_, simple_inserter );
	lsampler.set_min_bbrms( option[ lh::min_bbrms ]()   );
	lsampler.set_max_bbrms( option[ lh::max_bbrms ]() );
	lsampler.set_min_rms( option[ lh::min_rms ]() );
	lsampler.set_max_rms( option[ lh::max_rms ]() );
	lsampler.set_max_struct( skim_size );
	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_file( native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
	} else {
		utility_exit_with_message("This app requires specifying the -in:file:native flag.");
	}

	// Set up contraints
	ScoreFunctionOP fascorefxn = core::scoring::get_score_function();

	// convert pose to centroid pose:
	if ( !pose.is_fullatom() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
	}

	// pre relax!
	//protocols::relax::FastRelax *pre_relax = new protocols::relax::FastRelax( fascorefxn,  option[ OptionKeys::relax::sequence_file ]() );
	//pre_relax->apply( pose );


	core::Real last_accepted_score = (*fascorefxn)(pose);

	// See if a loopfile was defined - if so restrict sampling to those loops

	protocols::loops::Loops loops(true);
	utility::vector1< core::Size > selection;
	loops.get_residues( selection );

	TR << "Userdefined Loopregions: " << loops.size() << std::endl;
	TR << loops << std::endl;
	TR << "Residues: ";
	for ( core::Size i=1; i <= selection.size(); ++i ) TR <<  selection[i] << " ";
	TR << std::endl;

	//read_coord_cst(); //include this function later !

	for ( int round = 1; round <= option[ OptionKeys::lh::rounds ]; round ++ ) {
		core::Size total_starttime = time(nullptr);

		//static int casecount = 0;
		core::pose::Pose opose = pose;
		std::vector< core::io::silent::SilentStructOP > lib_structs;

		TR.Info << "Loophash apply function ! " << std::endl;

		//protocols::relax::FastRelax *qrelax = new protocols::relax::FastRelax( fascorefxn, 1 );
		protocols::relax::FastRelaxOP relax( new protocols::relax::FastRelax( fascorefxn,  option[ OptionKeys::relax::sequence_file ]() ) );

		// convert pose to centroid pose:
		core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
		core::pose::set_ss_from_phipsi( pose );


		// Generate alternate structures
		core::Size starttime2 = time(nullptr);
		core::Size sampler_chunk_size = 1;
		core::Size start_res = std::max( core::Size(2), core::Size(rand()%(pose.size() - sampler_chunk_size - 2 )) );
		core::Size stop_res  = std::min( core::Size(pose.size()), core::Size(start_res + sampler_chunk_size - 1 )  );

		// If a loopfile was given choose your insertion site from there
		TR.Info << "Selection size: " << selection.size() << std::endl;
		if ( selection.size() > 0 ) {
			utility::vector1< core::Size > temp_selection = selection;
			//std::random__shuffle( temp_selection.begin(), temp_selection.end());
			numeric::random::random_permutation(temp_selection.begin(), temp_selection.end(), numeric::random::rg());

			start_res = std::max( core::Size(2), core::Size( temp_selection[1] ) );
			stop_res  = std::min( core::Size(pose.size()), core::Size(start_res + sampler_chunk_size - 1)  );
			TR.Info << "SubselectionSample: " << start_res << " - " << stop_res << std::endl;
		}

		lsampler.set_start_res( start_res );
		lsampler.set_stop_res(  stop_res );
		lsampler.build_structures( pose, lib_structs );
		core::Size endtime2 = time(nullptr);
		core::Size loophash_time = endtime2 - starttime2;
		TR.Info << "FOUND (" << start_res << " to " << stop_res << "): "
			<< lib_structs.size() << " states in time: "
			<< endtime2 - starttime2 << " s " << std::endl;

		// try again if we have failed to find structures
		if ( lib_structs.size() == 0 ) continue;

		// choose up to "skim_size" of them
		numeric::random::random_permutation(lib_structs.begin(), lib_structs.end(), numeric::random::rg());

		std::vector< core::io::silent::SilentStructOP > select_lib_structs;
		for ( core::Size k=0; k< std::min(skim_size, lib_structs.size() ) ; k++ ) {
			select_lib_structs.push_back( lib_structs[k] );
		}

		core::pose::Pose ref_pose;
		if ( option[ lh::refstruct].user() ) {
			core::import_pose::pose_from_file( ref_pose, option[ lh::refstruct ]() , core::import_pose::PDB_file);
		}

		core::Real bestcenscore = MAXIMAL_FLOAT;
		core::Size bestcenindex = 0;
		for ( core::Size h = 0; h < select_lib_structs.size(); h++ ) {
			core::pose::Pose rpose;
			select_lib_structs[h]->fill_pose( rpose );

			//rpose.dump_pdb("struct_" + string_of(h) + ".pdb" );

			core::Real refrms = 0;
			if ( option[ lh::refstruct].user() ) {
				refrms = scoring::CA_rmsd( ref_pose, rpose );
			}
			core::Real rms_factor = 10.0;
			core::Real decoy_score = select_lib_structs[h]->get_energy("lh_censcore") + refrms * rms_factor;

			select_lib_structs[h]->add_energy( "refrms",     refrms,      1.0 );
			select_lib_structs[h]->add_energy( "comb_score", decoy_score, 1.0 );
			TR.Info << "refrms: " << refrms << "  Energy: " << decoy_score << std::endl;
			if ( decoy_score < bestcenscore ) {
				bestcenscore = decoy_score;
				bestcenindex = h;
			}
		}
		TR.Info << "Best:" << "  Energy: " << bestcenscore << std::endl;


		if ( (  option[ OptionKeys::lh::write_centroid_structs ]() ) ||
				(  option[ OptionKeys::lh::centroid_only ]() ) ) {

			core::io::silent::SilentFileOptions opts; // initialized from the command line
			core::io::silent::SilentFileData sfd(opts);
			std::string silent_file_ = option[ OptionKeys::out::file::silent ]();
			if ( option[ OptionKeys::lh::centroid_only ]() ) {
				silent_file_ += ".centroid.out" ;
			}

			for ( core::Size h = 0; h < select_lib_structs.size(); h++ ) {
				core::pose::Pose rpose;
				select_lib_structs[h]->fill_pose( rpose );
				core::Real rms = scoring::CA_rmsd( native_pose, rpose );
				select_lib_structs[h]->add_energy( "round", round, 1.0 );
				select_lib_structs[h]->add_energy( "rms", rms, 1.0 );
				select_lib_structs[h]->set_decoy_tag( "S_" + ObjexxFCL::string_of( round ) + "_" + ObjexxFCL::string_of(  h )  );
				sfd.write_silent_struct( *(select_lib_structs[h]) , silent_file_ );
			}

		}

		/// In centroid cleanup mode this is IT
		if ( option[ OptionKeys::lh::centroid_only ]() ) {
			select_lib_structs[bestcenindex]->fill_pose( pose );
			TR.Info << "Centroid mode. Skipping relax" << std::endl;
			continue;
		}


		/// For fullatom goodness, continue

		core::Size starttime = time(nullptr);
		relax->batch_apply( select_lib_structs );
		core::Size endtime = time(nullptr);
		core::Size batchrelax_time = endtime - starttime;
		TR.Info << "Batchrelax time: " << endtime - starttime << " for " << select_lib_structs.size() << " structures " << std::endl;

		core::Real bestscore = MAXIMAL_FLOAT;
		core::Size bestindex = 0;
		core::pose::Pose relax_winner;
		for ( core::Size h = 0; h < select_lib_structs.size(); h++ ) {
			TR.Info << "DOING: " << h << " / " << select_lib_structs.size() << std::endl;
			core::pose::Pose rpose;

			select_lib_structs[h]->fill_pose( rpose );

			//core::Real score = scoring::CA_rmsd( native_pose, rpose );
			core::Real score = (*fascorefxn)(rpose);
			TR.Info << "score: " << h << "  " << score << std::endl;

			select_lib_structs[h]->add_energy("lh_score_new", score );
			select_lib_structs[h]->add_energy("lh_score_old", last_accepted_score );
			select_lib_structs[h]->add_energy("lh_score_diff", score - last_accepted_score );

			if ( score < bestscore ) {
				bestscore = score;
				bestindex = h;
				relax_winner = rpose;
			}
		}
		//casecount++;


		TR.Info << "Metropolis decision: " << std::endl;
		// apply metropolis criterion
		core::Real new_energy = bestscore;
		core::Real old_energy = last_accepted_score;

		bool metropolis_replace = false;
		core::Real energy_diff_T = 0;
		if ( mpi_metropolis_temp_ > 0.0 ) energy_diff_T = old_energy - new_energy;

		if ( ( energy_diff_T >= 0.0 ) ) metropolis_replace = true; // energy of new is simply lower
		else if ( (energy_diff_T/mpi_metropolis_temp_) > (-10.0) ) {
			core::Real random_float = numeric::random::rg().uniform();
			if ( random_float < exp( energy_diff_T ) )  metropolis_replace = true;
		}

		TR.Info << "Metropolis: " << new_energy <<  ( (new_energy<old_energy)?" < ":" > ") << old_energy << " :" << (metropolis_replace?"ACC":"REJ") << std::endl;

		if ( metropolis_replace ) {
			pose = relax_winner;
			// Ok, pose has new content
			core::Real bestrms = scoring::CA_rmsd( native_pose, pose );
			TR.Info << "BESTSCORE: " << bestscore << "BESTRMS" << bestrms << std::endl;
			//pose.dump_pdb( "lhb_" + prefix + "_" + utility::to_string( round ) + ".pdb" );
			core::io::silent::SilentFileOptions opts; // initialized from the command line
			core::io::silent::SilentFileData sfd( opts );
			std::string silent_file_ = option[ OptionKeys::out::file::silent ]();
			for ( core::Size h = 0; h < select_lib_structs.size(); h++ ) {

				if ( h == bestindex || option[ OptionKeys::lh::write_all_fa_structs ]()  ) {
					core::pose::Pose rpose;
					select_lib_structs[h]->fill_pose( rpose );
					core::Real rms = scoring::CA_rmsd( native_pose, rpose );
					select_lib_structs[h]->add_energy( "round", round, 1.0 );
					select_lib_structs[h]->add_energy( "rms", rms, 1.0 );
					select_lib_structs[h]->set_decoy_tag( "S_" + ObjexxFCL::string_of( round ) + "_" + ObjexxFCL::string_of(  h )  );
					select_lib_structs[h]->sort_silent_scores();
					select_lib_structs[h]->print_score_header( std::cout );

					sfd.write_silent_struct( *(select_lib_structs[h]) , silent_file_ );
				}
			}
			last_accepted_score = new_energy;
		}  // end of metropolis acceptance if-block

		core::Size total_endtime = time(nullptr);

		TR.Info << "------------------------------------------------------------------------------------" << std::endl;
		TR.Info << " Energy: " << round << "  " << old_energy << "  " << new_energy << std::endl;
		TR.Info << " Timing: " << total_endtime - total_starttime << " L: " << loophash_time << " B: " << batchrelax_time << std::endl;
		TR.Info << "------------------------------------------------------------------------------------" << std::endl;


	}

}

void run_sandbox( LoopHashLibraryOP /*loop_hash_library*/ ){

	core::pose::Pose tgtpose, srcpose;
	core::import_pose::pose_from_file( tgtpose, "input/S_00001_0000001_0_0001.pdb" , core::import_pose::PDB_file);
	core::import_pose::pose_from_file( srcpose, "input/S_00001_0000001_0_1_0001.pdb" , core::import_pose::PDB_file);

	// test silent store class
	/*
	protocols::wum::SilentStructStore mystore;

	mystore.add( tgtpose );
	mystore.add( srcpose );


	std::string sdata;
	mystore.serialize( sdata );
	TR.Info << "----" << std::endl;
	mystore.print( TR.Info );


	protocols::wum::SilentStructStore mystore2;

	mystore2.read_from_cmd_line();
	TR.Info << "-------" << std::endl;
	mystore2.print( TR.Info );
	TR.Info << "-AA-AA-" << std::endl;

	mystore.add( mystore2 );

	TR.Info << "-------" << std::endl;

	mystore.print( TR.Info );

	std::string serial_form;
	mystore.serialize( serial_form );
	mystore.serialize_to_file( "testoutput.out" );


	protocols::wum::SilentStructStore recovered_store;
	recovered_store.read_from_string( serial_form );

	recovered_store.print( TR.Info );

	*/
	// Loopgraft
	// Test 1
	//37 63
	//loop_hash_library->graft_loop( srcpose, tgtpose, protocols::loops::Loop( 37, 63 ) );
	//tgtpose.dump_pdb( "output_graft.pdb" );


	// Test 2

	// results
	//21396 Fullatom protein
	//52499     "    binary
	//13879 Centroid protein
	//23918     "    binary
	// 2424 BackboneSegment float
	// 1212 BackboneSegment word

	{
		core::io::silent::SilentFileOptions opts; // initialized from the command line
		core::io::silent::ProteinSilentStruct pss(opts);
		pss.fill_struct(tgtpose, "lala" );

		std::ostringstream ss;
		pss.print_scores( ss );
		pss.print_conformation( ss );
		// TR.Info << ss.str() << std::endl;
		TR.Info << pss.mem_footprint() << "  " << ss.str().length() << std::endl;
	}

	//  {
	//    core::io::silent::BinarySilentStruct pss;
	//    pss.fill_struct(tgtpose, "lala" );
	//
	//    std::ostringstream ss;
	//    pss.print_scores( ss );
	//    pss.print_conformation( ss );
	//   // TR.Info << ss.str() << std::endl;
	//    TR.Info << pss.mem_footprint() << "  " << ss.str().length() << std::endl;
	//  }

	core::util::switch_to_residue_type_set( tgtpose, core::chemical::CENTROID_t );


	{
		core::io::silent::SilentFileOptions opts; // initialized from the command line
		core::io::silent::ProteinSilentStruct pss(opts);
		pss.fill_struct(tgtpose, "lala" );

		std::ostringstream ss;
		pss.print_scores( ss );
		pss.print_conformation( ss );
		// TR.Info << ss.str() << std::endl;
		TR.Info <<pss.mem_footprint() << "  " <<  ss.str().length() << std::endl;
	}

	{
		core::io::silent::SilentFileOptions opts; // initialized from the command line
		core::io::silent::BinarySilentStruct pss(opts);
		pss.fill_struct(tgtpose, "lala" );

		std::ostringstream ss;
		pss.print_scores( ss );
		pss.print_conformation( ss );
		// TR.Info << ss.str() << std::endl;
		//TR.Info <<pss.mem_footprint() << "  " <<  ss.str().length() << std::endl;
	}
}


int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols;
		using namespace protocols::jd2;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core;


		// initialize core
		devel::init(argc, argv);

#ifdef USEMPI
	int mpi_rank_, mpi_npes_;
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank_ ) );
	MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &mpi_npes_ ) );

	// Ok the point of the following is to make sure each worker starts at a different
	// point in the pseudo random sequence. Offsetting each by 1 is enough, since the same random
	// numbers will fall on different actions and the individual trajectories will diverge.
	// this still allows us to make overall reproducible runs when settign a constant overall start seed.
	core::Size random_sum=0;
	for ( int i = 0; i < mpi_rank_; i ++ ){
	  random_sum+=numeric::random::rg().random_range(1,65536)+rand()%65536;
	}
	TR << "Random Sum: " << random_sum    << std::endl;
	TR << "Random:     " << numeric::random::rg().uniform()  << std::endl;
	TR << "Random2:     " << rand() << std::endl;
#endif


		TR.Info << "SIZEOF short: " << sizeof( short ) << std::endl;
		TR.Info << "SIZEOF short*: " << sizeof( short * ) << std::endl;
		TR.Info << "SIZEOF core::Size: " << sizeof( core::Size ) << std::endl;
		TR.Info << "SIZEOF core::Size: " << sizeof( unsigned int ) << std::endl;

		utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
		LoopHashLibraryOP loop_hash_library( new LoopHashLibrary( loop_sizes ) );


		// sandbox mode ?
		if ( option[lh::sandbox]() ) { ;
			run_sandbox( loop_hash_library );
			return 0;
		}

		// Run simple sampling run test or create the db ?
		if ( option[lh::create_db]() ) { ;
			loop_hash_library->create_db();
			loop_hash_library->save_db();
			return 0;
		}

		LoopHashRelax_SamplerOP lh_sampler( new LoopHashRelax_Sampler( loop_hash_library ) );

		// Normal mode with external loophash library
		loop_hash_library->load_db();
		//protocols::jd2::JobDistributor::get_instance()->go( loop_hash_library );
		protocols::jd2::JobDistributor::get_instance()->go( lh_sampler );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
