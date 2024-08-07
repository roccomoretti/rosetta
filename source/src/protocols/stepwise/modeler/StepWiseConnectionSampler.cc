// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/StepWiseConnectionSampler.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/StepWiseConnectionSampler.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/packer/StepWiseMasterPacker.hh>
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_Closer.hh>
#include <protocols/stepwise/modeler/protein/loop_close/util.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosureChecker.hh>
#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/rna/rigid_body/util.hh>
#include <protocols/stepwise/modeler/rna/sugar/util.hh>
#include <protocols/stepwise/modeler/rna/sugar/VirtualSugarJustInTimeInstantiator.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/legacy/screener/RNA_AtrRepScreener.hh>
#include <protocols/stepwise/screener/BaseCentroidScreener.hh>
#include <protocols/stepwise/screener/BulgeApplier.hh>
#include <protocols/stepwise/screener/PartitionContactScreener.hh>
#include <protocols/stepwise/screener/RNA_ChainClosableGeometryScreener.hh>
#include <protocols/stepwise/screener/RNA_ChainClosableGeometryResidueBasedScreener.hh>
#include <protocols/stepwise/screener/RNA_ChainClosureScreener.hh>
#include <protocols/stepwise/screener/FastForwardToNextRigidBody.hh>
#include <protocols/stepwise/screener/FastForwardToNextResidueAlternative.hh>
#include <protocols/stepwise/screener/IntegrationTestBreaker.hh>
#include <protocols/stepwise/screener/AlignRMSD_Screener.hh>
#include <protocols/stepwise/screener/PoseSelectionScreener.hh>
#include <protocols/stepwise/screener/ProteinCCD_ClosureScreener.hh>
#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/StubApplier.hh>
#include <protocols/stepwise/screener/StubDistanceScreener.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <protocols/stepwise/screener/TagDefinition.hh>
#include <protocols/stepwise/screener/VDW_BinScreener.hh>
#include <protocols/stepwise/StepWiseSampleAndScreen.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/sampler/StepWiseSampler.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerComb.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerOneTorsion.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerOneDOF.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerRingConformer.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSamplerComb.hh>
#include <protocols/stepwise/sampler/protein/util.hh>
#include <protocols/stepwise/sampler/rna/util.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>
#include <protocols/toolbox/rigid_body/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/types.hh>

#include <utility/tools/make_vector1.hh>

#include <basic/Tracer.hh>

//Req's on WIN32

#include <protocols/stepwise/modeler/align/StepWiseClusterer.hh> // AUTO IWYU For StepWise...
#include <numeric/conversions.hh> // AUTO IWYU For radians

static basic::Tracer TR( "protocols.stepwise.modeler.StepWiseConnectionSampler" );

using namespace core;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::stepwise::modeler::rna::rigid_body;
using namespace protocols::stepwise::sampler::copy_dofs;
using namespace protocols::stepwise::sampler;

using core::pose::PoseOP;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Enumerates over:
//    main-chain torsions that build specified residues (for protein)
//     or...
//    suite backbone, sugar pucker, and base chi angle; and packs 2' hydroxyls (for RNA)
//  or...
// Enumerate rigid body degrees of freedom for a chunk of the pose.
//
// Decides which one to do based on what residue is upstream of the supplied moving_residue.
//
//
//        Reference residue        Moving residue                 Distal residue
//                                (virtual*)   (virtual*)        (virtual*)
//     5'   --Sugar -- ...       -- Phosphate -- Sugar -- ... -- Phosphate -- Sugar -- Phosphate -- 3'
//              |                       |                           |
//          Reference                Moving                      Distal
//             Base                   Base!                       Base
//              |_______________________________________________|
//                         jump or suite connection
//
//  * may be virtual -- the exceptions are if connected through chain closure to reference and/or distal residue.
//
// Note that sugar of the floating base will be instantiated if it needs to close chain to
//  anchor or distal residues. See 'to do' note below on how this could be improved.
//  In other cases, code below carries out a geometric 'sanity check' that involves temporarily
//  instantiating the sugar (within the screening_pose).
//
// The Sample-and-screen uses several poses and checkers to reduce the number of variant changes &
// energy recomputations -- faster code but somewhat complicated. Note that some information on, say, multiple
//  options for a sugar can be passed into  the sampler in 'residue_alternative_set'; some of the first screens
//  are based on base centroid distances and are independent of sugar conformation, and so some computation
//  that would be redundant across anchor sugar conformations can be avoided.
//
// As in other stepwise code,
//   StepWiseWorkingParameters holds information about this particular modeling job.
//   StepWiseModelerOptions holds information that is const across all of stepwise monte carlo.
//
// And, there are some 'useful info' variables derived below that are inferred from pose and above information -- use of KIC, etc. --
//   they are set at the beginning and should not change again (perhaps should store them together as a COP).
//
// Did not carry over:
//
//  protein KIC (which should be refactored anyway to be inside sampler, as Fang nicely carried out for RNA.)
//  build_pose_from_scratch (for RNA). BaseBinMap (for RNA rigid body).
//
// To do:
//
//  It should be possible to get rid of ResidueAlternativeSets -- this is a deeply refactored carryover from
//   parin's original SWA code. Could instead (1) sample base & sugar for residues that are at closure cutpoints but
//   are otherwise free (there's already a helper function for identifying these residues), and (2) have a screener
//   that checks for closure at virtual sugars. This might also deprecate "fast-forward" stuff, which is complicated.
//
//   -- Rhiju, 2014.
//    [after massive refactoring of code from Parin Sripakdeevong & Rhiju Das, 2009-2013]
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace modeler {


//Constructor
StepWiseConnectionSampler::StepWiseConnectionSampler( working_parameters::StepWiseWorkingParametersCOP & working_parameters ):
	working_parameters_( working_parameters ), // may deprecate
	moving_res_list_( working_parameters->working_moving_res_list() ),
	moving_res_( working_parameters->working_moving_res() ),
	moving_partition_res_( working_parameters->working_moving_partition_res() ),
	rigid_body_modeler_( false ),  // "useful_info" -- updated below
	reference_res_( 0 ),            // "useful_info" -- updated below
	kic_modeler_( false ),         // "useful_info" -- updated below
	protein_connection_( false ),   // "useful_info" -- updated below
	max_distance_squared_( 0.0 ),   // updated below
	virt_sugar_atr_rep_screen_( false )
{
	// not necessarily native -- just used for alignment & rmsd calcs.
	set_native_pose( working_parameters->working_native_pose() );
}

//Destructor
StepWiseConnectionSampler::~StepWiseConnectionSampler() = default;

/////////////////////
std::string
StepWiseConnectionSampler::get_name() const {
	return "StepWiseConnectionSampler";
}

///////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::apply( core::pose::Pose & pose ){

	clock_t const time_start( clock() );
	initialize_useful_info( pose );
	if ( !initialize_pose( pose ) ) return;
	if ( !initialize_sampler( pose ) ) return;
	initialize_screeners( pose );

	StepWiseSampleAndScreen sample_and_screen( sampler_, screeners_ );
	sample_and_screen.set_verbose( !options_->choose_random() || options_->integration_test_mode() || options_->verbose_sampler() );
	// Since get_max_ntries doesn't have access to the Pose, I modify it further here
	// based on whether the residue in question is a ligand.
	if ( moving_ligand_ ) {
		sample_and_screen.set_max_ntries( get_max_ntries() * 10 );
		sample_and_screen.set_num_random_samples( options_->num_random_samples() * 4 );
	} else {
		sample_and_screen.set_max_ntries( get_max_ntries() );
		sample_and_screen.set_num_random_samples( options_->num_random_samples() );
	}
	sample_and_screen.run();

	pose_list_ = clusterer_->pose_list(); // held in the last screener.

	TR.Debug << "Sampling time: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::initialize_useful_info( pose::Pose const & pose ){
	core::Size reference_res = toolbox::rigid_body::figure_out_reference_res_for_jump( pose, moving_res_ );
	// Complicate the rigid body modeler logic. This works iff the moving_res is NOT bonded to its connections
	// (i.e. JUMP_DOCK) but we need to manage to ignore cases where this is in fact the reference res for a ligand
	// or vrt but the edge is 'written backwards'
	rigid_body_modeler_ = reference_res > 0 && pose.aa( reference_res ) != core::chemical::aa_vrt;

	figure_out_reference_res( pose );
	protein_connection_ = is_protein( pose, moving_res_list_ ); // argh.
	// sets up rna_cutpoints_closed_, five_prime_chain_break_res_, three_prime_chain_break_res_ and chainbreak_gaps_
	rna::figure_out_moving_rna_chain_breaks( pose, moving_partition_res_,
		rna_cutpoints_closed_,
		rna_five_prime_chain_breaks_, rna_three_prime_chain_breaks_, rna_chain_break_gap_sizes_ );
	kic_modeler_ = ( !rigid_body_modeler_ && options_->kic_modeler_if_relevant() && ( rna_cutpoints_closed_.size() > 0 ) );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::figure_out_reference_res( pose::Pose const & pose ){
	// reference res is only defined when there is a single moving_res (for the downstream res)
	if ( rigid_body_modeler_ ) {
		figure_out_reference_res_with_rigid_body_rotamer( pose );
	} else if ( moving_res_list_.size() == 1 ) { // reference_res only makes sense with single residue.
		reference_res_ = figure_out_reference_res_for_suite( pose, moving_res_ );
		runtime_assert( moving_partition_res_ == figure_out_moving_partition_res_for_suite( pose, moving_res_, reference_res_ ) );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::figure_out_reference_res_with_rigid_body_rotamer( pose::Pose const & pose ){
	rigid_body_rotamer_ = utility::pointer::make_shared< sampler::rigid_body::RigidBodyStepWiseSampler >( pose, moving_res_ );
	reference_res_ = rigid_body_rotamer_->reference_res();
	runtime_assert( moving_partition_res_ == rigid_body_rotamer_->moving_partition_res() );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::initialize_screeners( pose::Pose & pose )
{
	initialize_checkers( pose );
	screeners_.clear();
	if ( rigid_body_modeler_ && pose.residue_type( moving_res_ ).is_polymer() ) initialize_residue_level_screeners( pose );
	initialize_pose_level_screeners( pose );
}

////////////////////////////////////////////////////////////////////////////////
// geometry checks that require base conformation (but not sugar) at moving_res
////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::initialize_residue_level_screeners( pose::Pose & pose ) {

	using namespace screener;
	using namespace core::conformation;
	using namespace protocols::stepwise::screener;
	using utility::tools::make_vector1;

	runtime_assert( rigid_body_rotamer_ != nullptr );

	screeners_.push_back( utility::pointer::make_shared< StubApplier >( moving_res_base_stub_ ) ); // will pull stub out of the sampler

	screeners_.push_back( utility::pointer::make_shared< StubDistanceScreener >( moving_res_base_stub_, rigid_body_rotamer_->reference_stub(), max_distance_squared_ ) );

	if ( base_centroid_checker_ && !protein_connection_ && pose.residue_type( moving_res_ ).is_NA() ) {
		screeners_.push_back( utility::pointer::make_shared< BaseCentroidScreener >( base_centroid_checker_,
			moving_res_base_stub_ ) );
	}

	tag_definition_ = utility::pointer::make_shared< TagDefinition >( pose, screeners_[ screeners_.size() ] );
	screeners_.push_back( tag_definition_ );

	for ( core::Size n = 1; n <= rna_five_prime_chain_breaks_.size(); n++ ) {
		screeners_.push_back( utility::pointer::make_shared< RNA_ChainClosableGeometryResidueBasedScreener >( rna_chain_closable_geometry_checkers_[ n ] ) );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// clash checks
	ResidueCOP screening_moving_rsd_at_origin = rigid_body_rotamer_->get_residue_at_origin( moving_res_ ).clone();
	if ( VDW_bin_checker_ != nullptr && !protein_connection_ && pose.residue_type( moving_res_ ).is_polymer() ) {
		screeners_.push_back( utility::pointer::make_shared< VDW_BinScreener >( VDW_bin_checker_, *virt_sugar_screening_pose_, moving_res_,
			screening_moving_rsd_at_origin, moving_res_base_stub_ ) );
	}
	// User-input VDW: Does not work for chain_closure move and is_internal_ move yet, since the checker does not know that
	// moving residue atoms can bond to previous or next residues.
	if ( user_input_VDW_bin_checker_ && user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() ) {
		screeners_.push_back( utility::pointer::make_shared< VDW_BinScreener >( user_input_VDW_bin_checker_, *virt_sugar_screening_pose_, moving_res_,
			screening_moving_rsd_at_origin, moving_res_base_stub_ ) );
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// 'Standard screens' that look at entire pose.
//////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::initialize_pose_level_screeners( pose::Pose & pose ) {

	using namespace screener;
	using namespace core::conformation;
	using namespace protocols::stepwise::screener;
	using utility::tools::make_vector1;

	//////////////////////////////////////////////////////////////////////////////////////////////
	// finally copy in the base conformation into the pose for RMSD and VDW screens
	// at first, don't copy in dofs for sugar/backbone (i.e., residue alternative), because those copy_dofs take extra computation;
	//  just apply rigid_body transformation.
	screeners_.push_back( utility::pointer::make_shared< SampleApplier >( *screening_pose_, false /*apply_residue_alternative_sampler*/ ) );

	AlignRMSD_ScreenerOP align_rmsd_screener;
	if ( ( get_native_pose() != nullptr ) && (moving_res_ != 0) &&
			( options_->rmsd_screen() > 0.0 || options_->integration_test_mode() ) ) {
		bool do_screen = ( ( options_->rmsd_screen() > 0.0 ) && !options_->integration_test_mode() ); // gets toggled to true in integration tests.
		// not necessarily native -- just used for alignment & rmsd calcs.

		// AMW: consider, here, expanding rmsd_screen for ligand sampling.
		align_rmsd_screener = utility::pointer::make_shared< AlignRMSD_Screener >( *get_native_pose(), *screening_pose_,
			working_parameters_->working_moving_partition_res() /*used to be working_parameters_->working_moving_res_list()*/,
			options_->rmsd_screen(), do_screen );
		screeners_.push_back( align_rmsd_screener );
	}

	// For KIC closure, immediate check that a closed loop solution was actually found.
	if ( kic_modeler_ ) {
		screeners_.push_back( utility::pointer::make_shared< RNA_ChainClosureScreener >( rna_chain_closure_checkers_[ 1 ], *screening_pose_, true /*just do closure check*/ ) );
	}

	// Following may still work, but has not been tested.
	// if ( options_->combine_long_loop_mode()  && ( rna_cutpoints_closed_.size() == 0 ) ) {
	//  screeners_.push_back( new ResidueContactScreener( *screening_pose_, last_append_res_,  last_prepend_res_, atom_atom_overlap_dist_cutoff_ ) );
	// }

	if ( !rigid_body_modeler_ && base_centroid_checker_ && !protein_connection_ && pose.residue_type( moving_res_ ).is_NA() ) {
		bool const force_centroid_interaction = ( rigid_body_modeler_ || options_->force_centroid_interaction()
			|| ( rna_cutpoints_closed_.size() == 0 ) );
		screeners_.push_back( utility::pointer::make_shared< BaseCentroidScreener >( base_centroid_checker_, screening_pose_, force_centroid_interaction ) );
	}

	if ( VDW_bin_checker_ ) {
		screeners_.push_back( utility::pointer::make_shared< VDW_BinScreener >( VDW_bin_checker_, *screening_pose_, moving_res_ ) );
	}
	if ( user_input_VDW_bin_checker_ && user_input_VDW_bin_checker_->user_inputted_VDW_screen_pose() ) {
		screeners_.push_back( utility::pointer::make_shared< VDW_BinScreener >( user_input_VDW_bin_checker_, *screening_pose_, moving_res_ ) );
	}

	// This one can actually filters out models compared to the other atr-rep screener below -- since it
	// involves a virtualized sugar, can fail atr screen.
	bool const use_loose_rep_cutoff = ( kic_modeler_ || moving_partition_res_.size() > 1 /* is_internal */ );
	if ( rigid_body_modeler_ && virt_sugar_atr_rep_screen_ &&
			moving_res_list_.size() > 0 && options_->atr_rep_screen() ) {
		screeners_.push_back( utility::pointer::make_shared< SampleApplier >( *virt_sugar_screening_pose_, false /*apply_residue_alternative_sampler*/ ) );
		screeners_.push_back( utility::pointer::make_shared< PartitionContactScreener >( *virt_sugar_screening_pose_, working_parameters_, use_loose_rep_cutoff, scorefxn_->energy_method_options() ) );
		//  screeners_.push_back( new RNA_AtrRepScreener( virt_sugar_atr_rep_checker_, *virt_sugar_screening_pose_ ) );
	}

	screeners_.push_back( utility::pointer::make_shared< SampleApplier >( *screening_pose_, true /*apply_residue_alternative_sampler*/ ) );
	for ( core::Size n = 1; n <= rna_cutpoints_closed_.size(); n++ ) screeners_.push_back( utility::pointer::make_shared< RNA_ChainClosableGeometryScreener >( rna_chain_closable_geometry_checkers_[ n ], screening_pose_ ) );

	for ( core::Size n = 1; n <= rna_five_prime_chain_breaks_.size(); n++ ) {
		bool strict = rigid_body_modeler_ && rna_cutpoints_closed_.has_value( rna_five_prime_chain_breaks_[n] ) && (rna_cutpoints_closed_.size() < 3);
		// if strict == false, no need to create a redundant screener, and its confusing. Remove?
		screeners_.push_back( utility::pointer::make_shared< RNA_ChainClosableGeometryScreener >( rna_chain_closable_geometry_checkers_[ n ], screening_pose_, strict  /*strict*/ ) );
	}

	// AMW: why no PCS for carbohydrate? don't care.
	// AMW: also disable for DNA. Not QUITE sure why.
	// AMW: ditto TNA. There is something up here. Maybe not necessary
	// AMW: also protein -- contacts are not common in (say) cyclic peptoids
	// but just hard to get results in some integration tests, or something?
	PartitionContactScreenerOP atr_rep_screener;
	if ( options_->atr_rep_screen() && moving_res_list_.size() > 0 &&
			( options_->atr_rep_screen_for_docking() || !rigid_body_modeler_ )
			&& pose.residue_type( moving_res_ ).is_RNA()
			&& !pose.residue_type( moving_res_ ).has_variant_type( core::chemical::DEOXY_O2PRIME ) ) {
		atr_rep_screener = utility::pointer::make_shared< PartitionContactScreener >( *screening_pose_, working_parameters_, use_loose_rep_cutoff, scorefxn_->energy_method_options() );
		screeners_.push_back( atr_rep_screener );
	}

	for ( core::Size n = 1; n <= rna_cutpoints_closed_.size(); n++ ) {
		if ( kic_modeler_ && n == 1 ) continue; // if KIC, first one is screened above, actually.
		if ( screening_pose_->residue_type( rna_cutpoints_closed_[n] ).is_NA() && screening_pose_->residue_type( rna_cutpoints_closed_[n] + 1 ).is_NA() ) {
			screeners_.push_back( utility::pointer::make_shared< RNA_ChainClosureScreener >( rna_chain_closure_checkers_[ n ] ) );
		}
	}

	for ( core::Size n = 1; n <= protein_ccd_closers_.size(); n++ )  screeners_.push_back( utility::pointer::make_shared< ProteinCCD_ClosureScreener >( protein_ccd_closers_[n], *protein_ccd_poses_[n] ) );

	// phosphate screener, o2prime screener, should be here.
	if ( !options_->lores() ) {
		master_packer_->add_packer_screeners( screeners_, pose, screening_pose_ );
	}

	//  Want this legacy BulgeApplier to run in very specific SWA RNA cases.
	// the atr_rep_checker has been replaced with a more general and efficient PartitionContactScreener, but
	//  still in use by BulgeApplier.
	if ( !rigid_body_modeler_ && !working_parameters_->rebuild_bulge_mode() &&
			options_->allow_bulge_at_chainbreak() && moving_partition_res_.size() == 1 &&
			( rna_cutpoints_closed_.size() > 0 )  && !protein_connection_ && !kic_modeler_ ) {
		runtime_assert( rna_atr_rep_checker_ != nullptr );
		legacy::screener::RNA_AtrRepScreenerOP rna_atr_rep_screener( new legacy::screener::RNA_AtrRepScreener( rna_atr_rep_checker_, *screening_pose_ ) );
		screeners_.push_back( StepWiseScreenerOP( rna_atr_rep_screener ) );  // will actually carry out the legacy RNA atr/rep check.
		screeners_.push_back( utility::pointer::make_shared< BulgeApplier >( rna_atr_rep_checker_ /* it turns out that this is not even used*/ , base_centroid_checker_, moving_res_ ) ); // apply bulge at the last minute.
	}

	/////////////////////////////////////////////////////
	screeners_.push_back( utility::pointer::make_shared< SampleApplier >( pose ) );

	if ( !tag_definition_ ) { // may have been defined above in residue level modeler.
		tag_definition_ = utility::pointer::make_shared< TagDefinition >( pose, screeners_[1], options_->sampler_include_torsion_value_in_tag(),
			moving_res_, reference_res_, "" /* extra_tag_ */ );
		screeners_.push_back( tag_definition_ );
	}

	screeners_.push_back( utility::pointer::make_shared< PoseSelectionScreener >( pose, scorefxn_, clusterer_ ) );

	if ( rigid_body_modeler_ ) {
		// As long as we're not instantiating any sugars in the final pose, break once we find any valid sugar rotamer [really?]
		if ( rna_chain_closure_checkers_.size() == 0 && !options_->sampler_try_sugar_instantiation() ) {
			screeners_.push_back( utility::pointer::make_shared< FastForwardToNextRigidBody >() ); // generalize to non-rigid-body rotamer case.
		} else {
			screeners_.push_back( utility::pointer::make_shared< FastForwardToNextResidueAlternative >( moving_res_ ) ); // generalize -- should fast forward past all virtual residue alternatives.
		}
	}

	if ( options_->integration_test_mode() ) {
		screeners_.insert( screeners_.begin() /*right at beginning!*/,
			utility::pointer::make_shared< IntegrationTestBreaker >( atr_rep_screener, screeners_[ screeners_.size() ], align_rmsd_screener ) );
	}

}

/////////////////////////////////////////////////////////////////////////////////////////
// initialization of variants for actual pose.
// also currently responsible for rna_chainbreaks_ & kic_modeler_ -- perhaps move this out.
bool
StepWiseConnectionSampler::initialize_pose( pose::Pose & pose  ){

	bool ready_for_more = presample_virtual_sugars( pose ); //defines alternatives for residues at chainbreaks with a virtual sugar.
	if ( !ready_for_more ) return false;

	// connection_sampler cannot handle waters - but might be better to handle this by virtualizing in screening pose,
	//  and then allowing StepWisePacker to do Mg(2+) hydration -- currently that's happening later in StepWiseMinimizer.
	// magnesium::remove_mg_bound_waters( pose, magnesium::get_mg_res( pose ), false /*leave other waters*/ );
	for ( core::Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue_type( n ).aa() == core::chemical::aa_h2o ) add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RESIDUE_VARIANT, n );
	}

	// note that constraint addition must happen after all variant changes (or the atom numbering will be off).
	// AMW: why is this not being handled by chainbreak scoreterm alone, if applicable??
	for ( core::Size n = 1; n <= rna_cutpoints_closed_.size(); n++ ) rna::add_harmonic_chain_break_constraint( pose, rna_cutpoints_closed_[ n ] );

	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::initialize_checkers( pose::Pose const & pose  ){

	using namespace rna::checker;

	// RNA Base Centroid stuff -- could soon generalize to find protein partners too. That would be cool.
	if ( working_parameters_->working_moving_partition_res().size() > 0 ) {
		base_centroid_checker_ = RNA_BaseCentroidCheckerOP(new RNA_BaseCentroidChecker ( pose, working_parameters_,
			( working_parameters_->floating_base() /*rigid body*/ && options_->tether_jump() ) ) );
		base_centroid_checker_->set_floating_base( working_parameters_->floating_base() &&
			working_parameters_->working_moving_partition_res().size() == 1  );
		base_centroid_checker_->set_allow_base_pair_only_screen( options_->allow_base_pair_only_centroid_screen() );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// set up screening pose -- do not change pose itself.
	// get rid of stuff that will be CCD-closed or packed (2'-OH, terminal phosphates) at last stages.
	// The idea is that if the screening pose fails basic stereochemistry checks, then we don't
	// have to carry out expensive CCD closure or packing.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////
	screening_pose_ = pose.clone();
	if ( options_->o2prime_legacy_mode() || options_->virtualize_packable_moieties_in_screening_pose() ) {
		core::pose::rna::add_virtual_O2Prime_hydrogen( *screening_pose_ );
	}
	if ( options_->virtualize_packable_moieties_in_screening_pose() ) phosphate::remove_terminal_phosphates( *screening_pose_ );
	for ( core::Size n = 1; n <= rna_cutpoints_closed_.size(); n++ ) {
		if ( screening_pose_->residue_type( rna_cutpoints_closed_[n] + 1 ).is_NA() ) {
			add_variant_type_to_pose_residue( *screening_pose_, core::chemical::VIRTUAL_PHOSPHATE,
				rna_cutpoints_closed_[n] + 1 ); // PS May 31, 2010 -- updated to all cutpoints by rhiju, feb. 2014
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////
	// VDW bin checker can take a while to set up... becomes rate-limiting in random, i.e.
	//  stepwise monte carlo runs.
	// Could in principle instantiate for proteins too.
	if ( !options_->choose_random() && !protein_connection_ && moving_res_ > 0 && !pose.residue_type( moving_res_ ).is_carbohydrate() ) {
		TR << TR.Magenta << "Creating VDW Bin Checker " << TR.Reset << std::endl;
		VDW_bin_checker_ = utility::pointer::make_shared< checker::RNA_VDW_BinChecker >();
		VDW_bin_checker_->setup_using_working_pose( *screening_pose_, working_parameters_ );
	}
	if ( !user_input_VDW_bin_checker_ /* could be externally defined for speed */ && options_->VDW_rep_screen_info().size() > 0 ) {
		TR << TR.Magenta << "Creating USER VDW Bin Checker " << TR.Reset << std::endl;
		user_input_VDW_bin_checker_ = utility::pointer::make_shared< RNA_VDW_BinChecker >( pose );
		options_->setup_options_for_VDW_bin_checker( user_input_VDW_bin_checker_ );
		user_input_VDW_bin_checker_->setup_using_user_input_VDW_pose( options_->VDW_rep_screen_info(), pose, working_parameters_ );
	}

	virt_sugar_screening_pose_ = screening_pose_->clone(); //Hard copy. Used for trying out sugar at moving residue.
	// virtual sugars even at residues that have instantiated sugars -- we can quickly screen this pose,
	// and it provides the appropriate baseline atr/rep for checking contacts and clashes.
	for ( core::Size n = 1; n <= residue_alternative_sets_.size(); n++ ) {
		if ( residue_alternative_sets_[ n ].size() > 1 ) {
			pose::add_variant_type_to_pose_residue( *virt_sugar_screening_pose_,
				core::chemical::VIRTUAL_RIBOSE, residue_alternative_sets_[ n ].representative_seqpos() );
		}
	}

	// following is to check atr/rep even on sugars that will remain virtualized.
	//  may now be deprecated due to development of PartitionContactScreener, which handles both protein & RNA.
	bool const use_loose_rep_cutoff = ( kic_modeler_ || moving_partition_res_.size() > 1 /* is_internal */ );
	if ( !rigid_body_modeler_ && options_->allow_bulge_at_chainbreak() && moving_partition_res_.size() == 1 && !protein_connection_ ) {  // need this for legacy bulge application code.
		rna_atr_rep_checker_ = utility::pointer::make_shared< RNA_AtrRepChecker >( *screening_pose_, working_parameters_, use_loose_rep_cutoff, scorefxn_->energy_method_options().clone() );
		// not in use anymore -- delete after MAR 2015 if SWA looks OK:
		//   rna_virt_sugar_atr_rep_checker_ = utility::pointer::make_shared< RNA_AtrRepChecker >( *virt_sugar_screening_pose_, working_parameters_, use_loose_rep_cutoff );
	}

	// we will be checking clashes of even virtual sugars compared to no-sugar baseline.
	for ( core::Size n = 1; n <= residue_alternative_sets_.size(); n++ ) {
		if ( residue_alternative_sets_[ n ].size() > 1 ) {
			pose::remove_variant_type_from_pose_residue( *screening_pose_,
				core::chemical::VIRTUAL_RIBOSE, residue_alternative_sets_[ n ].representative_seqpos() );
		}
	}

	for ( core::Size n = 1; n <= rna_five_prime_chain_breaks_.size(); n++ ) {
		rna_chain_closable_geometry_checkers_.push_back( utility::pointer::make_shared< RNA_ChainClosableGeometryChecker >(
			rna_five_prime_chain_breaks_[n], rna_three_prime_chain_breaks_[n], rna_chain_break_gap_sizes_[n] ) );
	}

	screening_pose_->remove_constraints(); // chain closure actually assumes no constraints in pose -- it sets up its own constraints.
	for ( core::Size n = 1; n <= rna_cutpoints_closed_.size(); n++ ) {
		// AMW TODO: I am not sure if this condition is safe to the rare case that we're looking at a cyclic cutpoint, in which case
		// instead of looking at rna_cutpoints_closed_[n] + 1 we either need to go into cyclize_res or we need to find the
		// similarly indexed rna_three_prime_chain_breaks_ for every rna_five_prime_chain_breaks_ or something.
		if ( screening_pose_->residue_type( rna_cutpoints_closed_[n] ).is_NA() && screening_pose_->residue_type( rna_cutpoints_closed_[n] + 1 ).is_NA() ) {
			rna_chain_closure_checkers_.push_back( utility::pointer::make_shared< RNA_ChainClosureChecker >( *screening_pose_, rna_cutpoints_closed_[n] ) );
			rna_chain_closure_checkers_[rna_chain_closure_checkers_.size()]->set_reinitialize_CCD_torsions( options_->reinitialize_CCD_torsions() );
		}
	}

	// AMW TODO: carbohydrate loop closers?

	// protein loop closers
	utility::vector1< core::Size > all_moving_res = get_moving_res_including_domain_boundaries( pose, moving_res_list_ );
	protein_cutpoints_closed_ = protein::just_protein( figure_out_moving_cutpoints_closed_from_moving_res( pose, all_moving_res ), pose );
	if ( !options_->kic_modeler_if_relevant() ) {
		// CCD closure -- heuristic closer but will accept fewer than 6 torsions (and does not require
		//  the torsions to be triaxial like kinematic closer).
		// Involves up to 5 residues: takeoff - bridge_res1 - bridge_res2 -bridge_res3 - landing
		for ( core::Size n = 1; n <= protein_cutpoints_closed_.size(); n++ ) {
			protein::loop_close::StepWiseProteinCCD_CloserOP protein_ccd_closer( new protein::loop_close::StepWiseProteinCCD_Closer( working_parameters_ ) );
			protein_ccd_closer->set_ccd_close_res( protein_cutpoints_closed_[ n ] );
			protein_ccd_closer->set_working_moving_res_list( moving_res_list_ );
			protein_ccd_closers_.push_back( protein_ccd_closer );

			PoseOP protein_ccd_pose =  pose.clone();
			protein_ccd_closer->init( *protein_ccd_pose );
			protein_ccd_poses_.push_back( protein_ccd_pose ); // does this really need to be an independent pose?
		}
	}

	master_packer_->reset( pose );

	clusterer_ = utility::pointer::make_shared< align::StepWiseClusterer >( options_ );
	clusterer_->set_max_decoys( get_num_pose_kept() );
	clusterer_->set_calc_rms_res( moving_res_list_ );
	if ( !options_->choose_random() ) {
		clusterer_->set_do_checks( false ); // for speed!
		clusterer_->set_assume_atom_ids_invariant( true ); // for speed!
	}
}

/////////////////////////////////////////////////////////////////////////////////
Size
StepWiseConnectionSampler::get_max_ntries() {
	core::Size max_ntries( 0 );
	if ( rigid_body_modeler_ ) {
		max_ntries = std::max( 100000, 1000 * int( options_->num_random_samples() ) );
		if ( options_->rmsd_screen() && !options_->integration_test_mode() ) max_ntries *= 10;
	} else {
		//  if ( options_->lores() ) return 100;
		max_ntries = std::max( 10000, 100 * int( options_->num_random_samples() ) );
		if ( rna_chain_closure_checkers_.size() > 0 ) {
			max_ntries *= options_->max_tries_multiplier_for_ccd() /* default 10 */;
		}
		if ( kic_modeler_ ) max_ntries = 5 * options_->num_random_samples(); // some chains just aren't closable.
		if ( protein_connection_ ) max_ntries = 500;
	}
	return max_ntries;
}

/////////////////////////////////////////////////////////////////////////////////
Size
StepWiseConnectionSampler::get_num_pose_kept() {
	core::Size num_pose_kept( 108 );
	if ( options_->sampler_num_pose_kept() > 0 ) num_pose_kept = options_->sampler_num_pose_kept();
	if ( !rigid_body_modeler_ && working_parameters_ && working_parameters_->sample_both_sugar_base_rotamer() ) num_pose_kept *= 4; //12;
	if ( rigid_body_modeler_ && base_centroid_checker_->allow_base_pair_only_screen() ) num_pose_kept *= 4;
	return num_pose_kept;
}

/////////////////////////////////////////////////////////////////////////////////////
Size
StepWiseConnectionSampler::which_residue_alternative_set_is_moving_residue() const {
	core::Size which_set_is_moving_res( 0 );
	for ( core::Size n = 1; n <= residue_alternative_sets_.size(); n++ ) {
		if ( residue_alternative_sets_[n].representative_seqpos() == moving_res_ ) which_set_is_moving_res = n;
	}
	return which_set_is_moving_res;
}

/////////////////////////////////////////////////////////////////////////////////////
// sets up a 'bare-bones' sugar modeling object with moving_residue information.
void
StepWiseConnectionSampler::initialize_moving_residue_pose_list( pose::Pose const & pose ){
	if ( which_residue_alternative_set_is_moving_residue() > 0 ) return; // already initialized.
	utility::vector1< pose::PoseOP > pose_list;
	if ( pose.residue_type( moving_res_ ).is_RNA() && // can't be a ligand, or a protein at the moment for that matter...
			rigid_body_modeler_ && moving_partition_res_.size() == 1 ) { // single floating base [classic]
		pose_list = setup_pose_with_moving_residue_alternative_list( pose, moving_res_, options_->extra_chi(), options_->use_phenix_geo() );
	} else {
		pose_list = utility::tools::make_vector1( pose.clone() ); // no alternatives.
	}
	sampler::copy_dofs::ResidueAlternativeSet residue_alternative_set( pose_list, moving_res_ );
	residue_alternative_sets_.push_back( residue_alternative_set );
}

/////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::initialize_euler_angle_grid_parameters(){
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// definition of euler angle grid search parameters
	// Following are set by #define in StepWiseRNA_FloatingBase_Samper_util_hh
	sampler::rigid_body::RigidBodyStepWiseSamplerValueRange & value_range = rigid_body_rotamer_->value_range();
	value_range.set_euler_angle_bin_size( STANDARD_EULER_ANGLE_BIN_SIZE );
	value_range.set_euler_z_bin_size( STANDARD_EULER_Z_BIN_SIZE );
	value_range.set_centroid_bin_size( STANDARD_CENTROID_BIN_SIZE );
	if ( options_->rmsd_screen() && moving_ligand_ ) {
		TR << "Centroid bin size is " << STANDARD_CENTROID_BIN_SIZE * ( options_->rmsd_screen() /  1.8 ) / 8.0  << std::endl;
		value_range.set_centroid_bin_size( STANDARD_CENTROID_BIN_SIZE * ( options_->rmsd_screen() /  1.8 ) / 8.0 );
	}
	if ( options_->integration_test_mode() ) { // use coarser search for speed
		value_range.set_euler_angle_bin_size( STANDARD_EULER_ANGLE_BIN_SIZE * 4 );
		value_range.set_centroid_bin_size( STANDARD_CENTROID_BIN_SIZE * 4);
	}
}

//////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::initialize_xyz_grid_parameters(){
	Distance max_distance = options_->sampler_max_centroid_distance(); // if unspecified (0.0), will be replaced
	// Reducing by 1.8 (rounding up sqrt(3))

	if ( options_->rmsd_screen() && moving_ligand_  /*&& max_distance == 0.0*/ ) {
		TR << "Max distance is now " << options_->rmsd_screen() / 1.8 << std::endl;
		max_distance = options_->rmsd_screen() /  1.8;
	}
	if ( options_->tether_jump() && max_distance == 0.0 ) max_distance = 8.0;
	int centroid_bin_min, centroid_bin_max;
	initialize_xyz_parameters( max_distance, max_distance_squared_,
		centroid_bin_min, centroid_bin_max,
		get_moving_rsd_list(), truly_floating_base() );
	sampler::rigid_body::RigidBodyStepWiseSamplerValueRange & value_range = rigid_body_rotamer_->value_range();
	value_range.set_max_distance( max_distance );
	value_range.set_centroid_bin_min( centroid_bin_min );
	value_range.set_centroid_bin_max( centroid_bin_max );
}

//////////////////////////////////////////////////////////////////////
// used in setting max_distance -- can use a tighter tether if moving res is covalently connected to the reference (anchor) residue.
Size  // core::Size instead of bool for historical reasons.
StepWiseConnectionSampler::truly_floating_base() {
	runtime_assert( reference_res_ > 0 ); // don't be in this function unless there is a jump from reference to moving!
	if ( rna_cutpoints_closed_.has_value( moving_res_    ) && reference_res_ == moving_res_+1 ) return 0;
	if ( rna_cutpoints_closed_.has_value( reference_res_ ) && reference_res_ == moving_res_-1 ) return 0;
	return 1;
}

//////////////////////////////////////////////////////////////////////
// kind of here for historical reasons -- set max distance based on
// list of all the options for moving_residue. This will probably
// be unnecessary soon when we shift to a mode where max_distance
// is a fixed number.
utility::vector1< core::conformation::ResidueOP >
StepWiseConnectionSampler::get_moving_rsd_list() const {
	runtime_assert( rigid_body_modeler_ );
	core::Size const n = which_residue_alternative_set_is_moving_residue();
	runtime_assert( n > 0 );
	utility::vector1< core::conformation::ResidueOP > moving_rsd_list;
	sampler::copy_dofs::ResidueAlternativeSet const & moving_res_modeling = residue_alternative_sets_[ n ];
	for ( core::Size n = 1; n <= moving_res_modeling.pose_list().size(); n++ ) {
		moving_rsd_list.push_back( moving_res_modeling.pose( n )->residue( moving_res_ ).clone() );
	}
	return moving_rsd_list;
}

//////////////////////////////////////////////////////////////////////
bool
StepWiseConnectionSampler::initialize_sampler( pose::Pose const & pose ){
	if ( moving_res_ != 0 && !pose.residue_type( moving_res_ ).is_polymer() ) {
		moving_ligand_ = true;
	}
	if ( rigid_body_modeler_ ) {
		initialize_euler_angle_grid_parameters();
		initialize_xyz_grid_parameters();
		initialize_full_rigid_body_sampler();
		if ( moving_ligand_ ) {
			TR << "About to initialize ligand bond sampler" << std::endl;
			sampler_ = initialize_ligand_bond_sampler( pose );
		}
	} else {
		// AMW: We now default to a 'bonded ligand' sampler.
		// This means that we DO NOT call initialize_rna_bond_sampler if moving_res_ == 0 anymore.
		if ( protein_connection_ ) {
			sampler_ = initialize_protein_bond_sampler( pose );
		} else if ( moving_res_ != 0 && pose.residue_type( moving_res_ ).is_carbohydrate() ) {
			TR << "Sampling carbohydrate" << std::endl;
			sampler_ = initialize_carbohydrate_bond_sampler( pose );
		} else if ( moving_res_ != 0 && pose.residue_type( moving_res_ ).is_TNA() ) {
			TR << "Sampling TNA" << std::endl;
			sampler_ = initialize_tna_bond_sampler( pose );
		} else if ( moving_res_ == 0 || pose.residue_type( moving_res_ ).is_RNA() ) {
			sampler_ = initialize_rna_bond_sampler( pose );
		} else { // if ( moving_res_ != 0 ) {
			runtime_assert( pose.residue_type( moving_res_ ).is_polymer() );
			sampler_ = initialize_generic_polymer_bond_sampler( pose );
		}
	}

	if ( sampler_ == nullptr ) {
		pose_list_.push_back( pose.clone() );
		( *scorefxn_ )( *( pose_list_[ 1 ]) );
		return false;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
StepWiseSamplerOP
StepWiseConnectionSampler::initialize_protein_bond_sampler( pose::Pose const & pose ){
	using namespace protocols::stepwise::sampler;
	using namespace protocols::stepwise::sampler::protein;

	utility::vector1< Real > allowed_values{ -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180 };
	utility::vector1< Real > dof_vals{
		numeric::conversions::radians( -160.0 ),
		numeric::conversions::radians( -140.0 ),
		numeric::conversions::radians( -120.0 ),
		numeric::conversions::radians( -100.0 ),
		numeric::conversions::radians(  -80.0 ),
		numeric::conversions::radians(  -60.0 ),
		numeric::conversions::radians(  -40.0 ),
		numeric::conversions::radians(  -20.0 ),
		numeric::conversions::radians(    0.0 ),
		numeric::conversions::radians(   20.0 ),
		numeric::conversions::radians(   40.0 ),
		numeric::conversions::radians(   60.0 ),
		numeric::conversions::radians(   80.0 ),
		numeric::conversions::radians(  100.0 ),
		numeric::conversions::radians(  120.0 ),
		numeric::conversions::radians(  140.0 ),
		numeric::conversions::radians(  160.0 ),
		numeric::conversions::radians(  180.0 ) };
	auto sampler = utility::pointer::make_shared< StepWiseSamplerSizedComb >();
	if ( moving_res_ == 0 ) { return nullptr; }
	auto e = pose.fold_tree().get_residue_edge( moving_res_ );

	if ( e.label() == kinematics::Edge::CHEMICAL ) {
		TR << "Setting up RNA-bonded protein sampler for " << pose.residue_type( moving_res_ ).name() << std::endl;
		using namespace core::id;
		// residue is *built by* the branch
		sampler->add_external_loop_rotamer(
			utility::pointer::make_shared< StepWiseSamplerOneDOF >(
			core::id::DOF_ID( AtomID( pose.residue_type( moving_res_ ).atom_index( "C" ), moving_res_ ), id::PHI ), dof_vals ) );
		sampler->add_external_loop_rotamer(
			utility::pointer::make_shared< StepWiseSamplerOneDOF >(
			core::id::DOF_ID( AtomID( pose.residue_type( moving_res_ ).atom_index( "CA" ), moving_res_ ), id::PHI ), dof_vals ) );
		sampler->add_external_loop_rotamer(
			utility::pointer::make_shared< StepWiseSamplerOneDOF >(
			core::id::DOF_ID( AtomID( pose.residue_type( moving_res_ ).atom_index( "N" ), moving_res_ ), id::PHI ), dof_vals ) );
		// these are actually the last two BB!
	} else {

		auto protein_sampler = get_basic_protein_sampler( pose, moving_res_list_,
			working_parameters_, options_, input_streams_ );
		// protein_sampler->init();

		if ( protein_cutpoints_closed_.size() > 0 && options_->kic_modeler_if_relevant() ) {
			runtime_assert( protein_cutpoints_closed_.size() == 1 );
			protein::loop_close::kic_close_loops_in_samples( protein_sampler, pose, working_parameters_, options_ );
		}

		// need to provide BB-omitting flag to this IMO.
		sampler->add_external_loop_rotamer( protein_sampler );
	}

	sampler->set_random( options_->choose_random() );
	sampler->init();

	return sampler;
}

//////////////////////////////////////////////////////////////////////
sampler::StepWiseSamplerOP
StepWiseConnectionSampler::initialize_rna_bond_sampler( pose::Pose const & pose ){
	using namespace sampler;
	if ( moving_res_list_.size() == 0 ) return nullptr;
	sampler::StepWiseSamplerOP sampler_ = sampler::rna::setup_sampler( pose, options_,
		working_parameters_, false /*build_pose_from_scratch_*/,
		kic_modeler_, (rna_cutpoints_closed_.size() > 0) );
	ResidueAlternativeStepWiseSamplerCombOP rsd_alternatives_rotamer = get_rsd_alternatives_rotamer();
	if ( rsd_alternatives_rotamer == nullptr ) return sampler_;

	StepWiseSamplerCombOP sampler( new StepWiseSamplerComb );
	sampler->add_external_loop_rotamer( rsd_alternatives_rotamer );
	sampler->add_external_loop_rotamer( sampler_ );
	sampler->set_random( options_->choose_random() );
	sampler->init();
	return sampler;
}


//////////////////////////////////////////////////////////////////////
sampler::StepWiseSamplerOP
StepWiseConnectionSampler::initialize_tna_bond_sampler( pose::Pose const & pose ){
	using namespace sampler;
	if ( moving_res_list_.size() == 0 ) return nullptr;
	sampler::StepWiseSamplerOP sampler_ = sampler::rna::setup_sampler( pose, options_,
		working_parameters_, false /*build_pose_from_scratch_*/,
		kic_modeler_, (rna_cutpoints_closed_.size() > 0) );
	ResidueAlternativeStepWiseSamplerCombOP rsd_alternatives_rotamer = get_rsd_alternatives_rotamer();
	if ( rsd_alternatives_rotamer == nullptr ) return sampler_;

	StepWiseSamplerCombOP sampler( new StepWiseSamplerComb );
	sampler->add_external_loop_rotamer( rsd_alternatives_rotamer );
	sampler->add_external_loop_rotamer( sampler_ );
	sampler->set_random( options_->choose_random() );
	sampler->init();
	return sampler;
}


//////////////////////////////////////////////////////////////////////
sampler::StepWiseSamplerOP
StepWiseConnectionSampler::initialize_generic_polymer_bond_sampler( pose::Pose const & pose ) {
	using namespace sampler;
	if ( moving_res_list_.size() == 0 ) return nullptr;

	// AMW: Note that this won't work if the residue is LINKed to another residue, but it does not
	// have the polymer property...
	utility::vector1< Real > allowed_values{ -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180 };
	utility::vector1< Real > dof_vals{
		numeric::conversions::radians( -160.0 ),
		numeric::conversions::radians( -140.0 ),
		numeric::conversions::radians( -120.0 ),
		numeric::conversions::radians( -100.0 ),
		numeric::conversions::radians(  -80.0 ),
		numeric::conversions::radians(  -60.0 ),
		numeric::conversions::radians(  -40.0 ),
		numeric::conversions::radians(  -20.0 ),
		numeric::conversions::radians(    0.0 ),
		numeric::conversions::radians(   20.0 ),
		numeric::conversions::radians(   40.0 ),
		numeric::conversions::radians(   60.0 ),
		numeric::conversions::radians(   80.0 ),
		numeric::conversions::radians(  100.0 ),
		numeric::conversions::radians(  120.0 ),
		numeric::conversions::radians(  140.0 ),
		numeric::conversions::radians(  160.0 ),
		numeric::conversions::radians(  180.0 ) };
	utility::vector1< Real > chi_values{ -150, -120, -90, -60, -30, 0, 30, 60, 90, 120, 150, 180 };
	utility::vector1< Real > symm_chi_values{ -150, -120, -90, -60, -30, 0 };
	utility::vector1< Real > rotchi{ -60.0, 60.0, 180.0 };
	StepWiseSamplerCombOP sampler( new StepWiseSamplerComb );
	TR << "Setting up bonded polymer sampler for " << pose.residue_type( moving_res_ ).name() << std::endl;
	// AMW: be very careful here because it's not clear all of these will be very useful.

	auto e = pose.fold_tree().get_residue_edge( moving_res_ );

	bool omit_first_two_bb = false;
	if ( e.label() == kinematics::Edge::CHEMICAL ) {
		// we must add BRANCH torsions.
		//Size other_seqpos = ( e.start() == moving_res_ ) ? e.stop() : e.start();
		//TR.Warning << "pose ft is " << pose.fold_tree() << std::endl;
		//utility_exit();

		/*
		id::DOF_ID
		AtomTree::torsion_angle_dof_id(
		AtomID const & atom1_in_id,
		AtomID const & atom2_in_id,
		AtomID const & atom3_in_id,
		AtomID const & atom4_in_id,
		Real & offset,
		bool const quiet
		) const
		{*/
		using namespace core::id;
		// residue is *built by* the branch
		sampler->add_external_loop_rotamer(
			utility::pointer::make_shared< StepWiseSamplerOneDOF >(
			core::id::DOF_ID( AtomID( pose.residue_type( moving_res_ ).atom_index( "C" ), moving_res_ ), id::PHI ), dof_vals ) );
		sampler->add_external_loop_rotamer(
			utility::pointer::make_shared< StepWiseSamplerOneDOF >(
			core::id::DOF_ID( AtomID( pose.residue_type( moving_res_ ).atom_index( "CA" ), moving_res_ ), id::PHI ), dof_vals ) );
		// these are actually the last two BB!
		omit_first_two_bb = true;
	}

	for ( core::Size ii = pose.residue_type( moving_res_ ).mainchain_atoms().size() - 1; ii >= 1; --ii ) {
		if ( omit_first_two_bb && ii >= pose.residue_type( moving_res_ ).mainchain_atoms().size() - 1 ) continue;
		// skip if lower connect not there. Might be bad if it exists but is unfulfilled?
		// ideally skip if it's there but unfulfilled
		//if ( !pose.residue_type( moving_res_ ).lower_connect_id() && ii == 1 ) continue;
		// This was a good idea... but... what if you naturally don't have an upper
		//if ( !pose.residue_type( moving_res_ ).upper_connect_id() && ii >= pose.residue_type( moving_res_ ).mainchain_atoms().size() - 1 ) continue;

		// AMW TODO: use functions in conformation or something to only add sensible ones.
		sampler->add_external_loop_rotamer( utility::pointer::make_shared< StepWiseSamplerOneTorsion >( core::id::TorsionID( moving_res_, id::BB, ii ), allowed_values ) );
	}
	for ( core::Size ii = 1; ii <= pose.residue_type( moving_res_ ).nchi() - pose.residue_type( moving_res_ ).n_proton_chi(); ++ii ) {
		// nrchi
		if ( ii == 2 && (
				pose.residue_type( moving_res_ ).aa() == chemical::aa_tyr ||
				pose.residue_type( moving_res_ ).aa() == chemical::aa_phe ||
				pose.residue_type( moving_res_ ).aa() == chemical::aa_asp ) ) {
			sampler->add_external_loop_rotamer( utility::pointer::make_shared< StepWiseSamplerOneTorsion >( core::id::TorsionID( moving_res_, id::CHI, ii ), symm_chi_values ) );
		} else if ( ii == 3 && pose.residue_type( moving_res_ ).aa() == chemical::aa_glu ) {
			sampler->add_external_loop_rotamer( utility::pointer::make_shared< StepWiseSamplerOneTorsion >( core::id::TorsionID( moving_res_, id::CHI, ii ), symm_chi_values ) );
		} else if ( ii == pose.residue_type( moving_res_ ).nchi() ) {
			// final asymmetric nrchi: trp, gln, asn
			sampler->add_external_loop_rotamer( utility::pointer::make_shared< StepWiseSamplerOneTorsion >( core::id::TorsionID( moving_res_, id::CHI, ii ), allowed_values ) );
		} else {
			sampler->add_external_loop_rotamer( utility::pointer::make_shared< StepWiseSamplerOneTorsion >( core::id::TorsionID( moving_res_, id::CHI, ii ), rotchi ) );
		}
	}

	sampler->set_random( options_->choose_random() );
	sampler->init();
	return sampler;
}

//////////////////////////////////////////////////////////////////////
sampler::StepWiseSamplerOP
StepWiseConnectionSampler::initialize_ligand_bond_sampler( pose::Pose const & pose ){
	using namespace sampler;
	if ( moving_res_list_.size() == 0 ) return nullptr;

	utility::vector1< Real > allowed_values{ -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180 };
	StepWiseSamplerCombOP sampler( new StepWiseSamplerComb );
	TR << "Setting up ligand sampler for " << pose.residue_type( moving_res_ ).name() << std::endl;
	TR << "\t" << pose.residue_type( moving_res_ ) << std::endl;
	for ( core::Size ii = 1; ii <= pose.residue_type( moving_res_ ).nchi(); ++ii ) {
		StepWiseSamplerOneTorsionOP foo( new StepWiseSamplerOneTorsion( core::id::TorsionID( moving_res_, id::CHI, ii ), allowed_values ) );
		sampler->add_external_loop_rotamer( foo );
	}
	sampler->add_external_loop_rotamer( rigid_body_rotamer_ );
	sampler->set_random( options_->choose_random() );
	sampler->init();
	return sampler;
}

//////////////////////////////////////////////////////////////////////
sampler::StepWiseSamplerOP
StepWiseConnectionSampler::initialize_carbohydrate_bond_sampler( pose::Pose const & pose ){
	using namespace sampler;
	if ( moving_res_list_.size() == 0 ) return nullptr;

	StepWiseSamplerCombOP sampler( new StepWiseSamplerComb );
	// Really should be sampling ring conformations or something. We'd need a special class for that.
	for ( core::Size ii = 1; ii <= pose.residue_type( moving_res_ ).n_ring_conformer_sets(); ++ii ) {
		StepWiseSamplerRingConformerOP rc( new StepWiseSamplerRingConformer( ii, moving_res_, pose.residue_type( moving_res_ ).ring_conformer_set( ii ) ) );
		sampler->add_external_loop_rotamer( rc );
	}

	//for ( core::Size ii = 1; ii <= pose.residue_type( moving_res_ ).nchi(); ++ii ) {
	utility::vector1< Real > allowed_values{ -160, -140, -120, -100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180 };
	// Chi 1 is also the connection to lower. As a pseudo-MC torsion, it should be sampled more finely
	StepWiseSamplerOneTorsionOP foo( new StepWiseSamplerOneTorsion( core::id::TorsionID( moving_res_, id::CHI, 1 ), allowed_values ) );
	sampler->add_external_loop_rotamer( foo );

	utility::vector1< Real > sp3_rots{ -60, 60, 180 };
	for ( core::Size ii = 2; ii <= pose.residue_type( moving_res_ ).nchi(); ++ii ) {
		// There will (usually?) be one chi that is the hydroxyl that is "really a mainchain torsion" in the connected form.
		// that isn't necessarily cause not to sample it, but be aware: if chi atom 3 is final mainchain atom, be wary.
		StepWiseSamplerOneTorsionOP foo( new StepWiseSamplerOneTorsion( core::id::TorsionID( moving_res_, id::CHI, ii ), sp3_rots ) );
		sampler->add_external_loop_rotamer( foo );
	}
	sampler->init();
	return sampler;
}


//////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::initialize_full_rigid_body_sampler(){

	using namespace sampler::rigid_body;
	using namespace sampler::copy_dofs;

	ResidueAlternativeStepWiseSamplerCombOP rsd_alternatives_rotamer = get_rsd_alternatives_rotamer();
	sampler_ = utility::pointer::make_shared< RigidBodyStepWiseSamplerWithResidueAlternatives >( rsd_alternatives_rotamer, rigid_body_rotamer_ );
	sampler_->set_random( options_->choose_random() );
	sampler_->init();
}

//////////////////////////////////////////////////////////////////////
sampler::copy_dofs::ResidueAlternativeStepWiseSamplerCombOP
StepWiseConnectionSampler::get_rsd_alternatives_rotamer(){

	if ( residue_alternative_sets_.size() == 0 ) return nullptr;
	ResidueAlternativeStepWiseSamplerCombOP rsd_alternatives_rotamer( new ResidueAlternativeStepWiseSamplerComb() );
	// note that following will include moving_res_ for sampler, as well as any other chunks that might move...
	for ( ResidueAlternativeSet const & residue_alternative_set : residue_alternative_sets_ ) {
		ResidueAlternativeStepWiseSamplerOP rsd_alt_rotamer;
		if ( rigid_body_rotamer_ != nullptr ) {
			rsd_alt_rotamer = utility::pointer::make_shared< ResidueAlternativeStepWiseSampler >( residue_alternative_set,
				*rigid_body_rotamer_->pose_at_origin() /*take representative residues from this pose after applying copy_dofs*/);
		} else {
			rsd_alt_rotamer = utility::pointer::make_shared< ResidueAlternativeStepWiseSampler >( residue_alternative_set );
		}
		rsd_alternatives_rotamer->add_residue_alternative_rotamer( rsd_alt_rotamer );
	}
	return rsd_alternatives_rotamer;
}

void
StepWiseConnectionSampler::add_residue_alternative_set( sampler::copy_dofs::ResidueAlternativeSet const & residue_alternative_set ){
	residue_alternative_sets_.push_back( residue_alternative_set );
}

////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< PoseOP > &
StepWiseConnectionSampler::get_pose_list(){ return pose_list_; }

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker ){ user_input_VDW_bin_checker_ = user_input_VDW_bin_checker; }

//////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::set_pose_list( utility::vector1< pose::PoseOP > & pose_list ){
	pose_list_ = pose_list;
}

//////////////////////////////////////////////////////////////////
void
StepWiseConnectionSampler::set_options( options::StepWiseModelerOptionsCOP options ){
	options_ = options;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
StepWiseConnectionSampler::presample_virtual_sugars( pose::Pose & pose ){

	using namespace modeler::rna::sugar;
	if ( moving_res_ == 0 ) return true;
	VirtualSugarJustInTimeInstantiatorOP virtual_sugar_just_in_time_instantiator =
		instantiate_any_virtual_sugars( pose, working_parameters_, scorefxn_, options_ );
	if ( !virtual_sugar_just_in_time_instantiator->success() ) return false;
	virtual_sugar_just_in_time_instantiator->instantiate_sugars_at_cutpoint_closed( pose );

	// in backwards order to match some old runs.
	for ( core::Size n = virtual_sugar_just_in_time_instantiator->num_sets(); n >= 1; n-- ) {
		add_residue_alternative_set( virtual_sugar_just_in_time_instantiator->residue_alternative_set( n ) );
	}

	// sets up residue alternatives for moving_res in rigid_body runs.
	if ( rigid_body_modeler_ ) initialize_moving_residue_pose_list( pose );
	return true;
}


} //modeler
} //stepwise
} //protocols
