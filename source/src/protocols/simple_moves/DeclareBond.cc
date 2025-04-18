// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song
/// @author Modified by Vikram K. Mulligan (vmullig@uw.edu)

#include <protocols/simple_moves/DeclareBond.hh>
#include <protocols/simple_moves/DeclareBondCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>

#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>


#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.simple_moves.DeclareBond" );

namespace protocols {
namespace simple_moves {

namespace {
///@brief run the residue selector on the pose, assert that only one residue is selected, return the resid for that selected residue
core::Size
get_selected_residue(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSelector const & selector
) {
	utility::vector1< bool > const sele = selector.apply( pose );
	utility::vector1< core::Size > const positions =
		core::select::residue_selector::selection_positions( sele );

	if ( positions.size() != 1 ) {
		utility_exit_with_message( "DeclareBond expected a residue selector to select exactly one residues. Instead, it selected " + std::to_string( positions.size() ) );
	}

	return positions[ 1 ];
}
}

DeclareBond::DeclareBond():
	res1_(0),
	atom1_(""),
	res2_(0),
	atom2_(""),
	add_termini_(false),
	run_kic_(false),
	kic_res1_(0),
	kic_res2_(0),
	rebuild_fold_tree_(false)
{}
DeclareBond::~DeclareBond()= default;

void
DeclareBond::set( core::Size const res1,
	std::string const & atom1,
	core::Size const res2,
	std::string const & atom2,
	bool const add_termini,
	bool const run_kic,
	core::Size const kic_res1,
	core::Size const kic_res2,
	bool const rebuild_fold_tree
)
{
	res1_ = res1;
	atom1_ = atom1;
	res2_ = res2;
	atom2_ = atom2;
	add_termini_ = add_termini;
	run_kic_ = run_kic;
	kic_res1_ = kic_res1;
	kic_res2_ = kic_res2;
	rebuild_fold_tree_ = rebuild_fold_tree;
}

void DeclareBond::apply( core::pose::Pose & pose )
{
	using namespace core::chemical;

	if ( selector1_ != nullptr ) {
		res1_ = get_selected_residue( pose, *selector1_ );
	}
	runtime_assert( res1_ != 0 );

	if ( selector2_ != nullptr ) {
		res2_ = get_selected_residue( pose, *selector2_ );
	}
	runtime_assert( res2_ != 0 );

	//printf("Stripping termini.\n"); fflush(stdout); //DELETE ME
	if ( atom1_=="N" && ( pose.residue_type(res1_).is_alpha_aa() || pose.residue_type(res1_).is_beta_aa() || pose.residue_type(res1_).is_peptoid() ) && pose.residue(res1_).has_variant_type(LOWER_TERMINUS_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue(pose, LOWER_TERMINUS_VARIANT, res1_);
	}
	if ( atom2_=="N" && ( pose.residue_type(res2_).is_alpha_aa() || pose.residue_type(res2_).is_beta_aa() || pose.residue_type(res2_).is_peptoid() ) && pose.residue(res2_).has_variant_type(LOWER_TERMINUS_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue(pose, LOWER_TERMINUS_VARIANT, res2_);
	}
	if ( atom1_=="C" && ( pose.residue_type(res1_).is_alpha_aa() || pose.residue_type(res1_).is_beta_aa() || pose.residue_type(res1_).is_peptoid() ) && pose.residue(res1_).has_variant_type(UPPER_TERMINUS_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue(pose, UPPER_TERMINUS_VARIANT, res1_);
	}
	if ( atom2_=="C" && ( pose.residue_type(res2_).is_alpha_aa() || pose.residue_type(res2_).is_beta_aa() || pose.residue_type(res1_).is_peptoid() ) && pose.residue(res2_).has_variant_type(UPPER_TERMINUS_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue(pose, UPPER_TERMINUS_VARIANT, res2_);
	}

	//printf("Declaring bond.\n"); fflush(stdout); //DELETE ME
	pose.conformation().declare_chemical_bond(res1_, atom1_, res2_, atom2_);

	//Rebuild the polymer bond dependent atoms:
	//printf("Rebuilding bond-dependent atoms.\n"); fflush(stdout); //DELETE ME
	pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(res1_);
	pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(res2_);
	pose.conformation().rebuild_residue_connection_dependent_atoms( res1_, pose.residue_type(res1_).residue_connection_id_for_atom( pose.residue_type(res1_).atom_index(atom1_) ) );
	pose.conformation().rebuild_residue_connection_dependent_atoms( res2_, pose.residue_type(res2_).residue_connection_id_for_atom( pose.residue_type(res2_).atom_index(atom2_) ) );

	if ( rebuild_fold_tree_ ) {
		core::pose::Pose const pose_copy(pose); //Make a reference copy of pose (const to prevent accidentally altering it).

		pose.clear();

		for ( core::Size ires=1; ires<=pose_copy.size(); ++ires ) {
			if ( ires==1 ) {
				pose.append_residue_by_jump(pose_copy.residue(ires),1);
			} else {
				core::Size anchor_rsd(0);
				core::Size anchor_conid(0);
				core::Size icon=1;
				for ( ; icon<=pose_copy.residue_type(ires).n_possible_residue_connections(); ++icon ) {
					if ( pose_copy.residue(ires).connected_residue_at_resconn(icon) != 0 ) {
						if ( pose_copy.residue(ires).connected_residue_at_resconn(icon) < ires ) {
							anchor_rsd = pose_copy.residue(ires).connected_residue_at_resconn(icon);
							anchor_conid = pose_copy.residue(ires).connect_map(icon).connid();
							break;
						}
					}
				}
				if ( anchor_rsd != 0 ) {
					pose.append_residue_by_bond(pose_copy.residue(ires),false, icon, anchor_rsd, anchor_conid);
				} else {
					if ( pose_copy.chain(ires-1)!=pose_copy.chain(ires) ) { //If this is a new chain, connect this by a jump to residue 1, starting a new chain.
						pose.append_residue_by_jump(pose_copy.residue(ires), 1, "", "", true);
					} else {
						pose.append_residue_by_jump(pose_copy.residue(ires), ires-1);
					}
				}
			}
		}

		// add back all the connections
		for ( core::Size ires=1; ires<=pose_copy.size(); ++ires ) {
			for ( core::Size icon=1; icon<=pose_copy.residue_type(ires).n_possible_residue_connections(); ++icon ) {
				if ( pose_copy.residue(ires).connected_residue_at_resconn(icon) != 0 ) {
					core::Size anchor_rsd = pose_copy.residue(ires).connected_residue_at_resconn(icon);
					core::Size anchor_conid = pose_copy.residue(ires).connect_map(icon).connid();

					if ( pose.residue(ires).connected_residue_at_resconn(icon) == 0 ) {
						//if (pose.residue_type(ires).name3() != "CYS") {
						pose.conformation().set_noncanonical_connection(ires, icon, anchor_rsd, anchor_conid);
						TR << "Adding connection Res " << ires << " to Residue " << anchor_rsd << std::endl;
						//}
					}
				}
			}
		}
	}

	if ( add_termini_ ) {
		for ( core::Size ir=1; ir<=pose.size() ; ++ir ) {
			if ( pose.residue_type(ir).lower_connect_id() != 0 ) {
				if ( pose.residue(ir).connected_residue_at_lower() == 0 ) {
					if ( pose.residue(ir).has_variant_type(CUTPOINT_UPPER) ) {
						core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, ir );
					}
					core::pose::add_variant_type_to_pose_residue(pose, LOWER_TERMINUS_VARIANT, ir);
				}
			}
			if ( pose.residue_type(ir).upper_connect_id() != 0 ) {
				if ( pose.residue(ir).connected_residue_at_upper() == 0 ) {
					if ( pose.residue(ir).has_variant_type(CUTPOINT_LOWER) ) {
						core::pose::remove_variant_type_from_pose_residue(pose, CUTPOINT_LOWER, ir);
					}
					core::pose::add_variant_type_to_pose_residue(pose, UPPER_TERMINUS_VARIANT, ir);
				}
			}
		}
	}

	//Kinematic closure to build the rest of the peptide:
	if ( run_kic_ ) {
		protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kinmover( new protocols::loops::loop_closure::kinematic_closure::KinematicMover );
		core::scoring::ScoreFunctionOP sfxn;
		sfxn = core::scoring::get_score_function();
		if ( kic_res1_ == 0 ) {
			kic_res1_ = 1;
		}
		if ( kic_res2_ == 0 ) {
			kic_res2_ = pose.size();
		}

		kinmover->set_temperature( 1.0 );
		kinmover->set_vary_bondangles( false );
		kinmover->set_sample_nonpivot_torsions( true );
		kinmover->set_rama_check( true );
		kinmover->set_idealize_loop_first( true );
		kinmover->set_sfxn(sfxn);
		kinmover->set_pivots(kic_res1_, (int)( ( (core::Real)(kic_res1_+kic_res2_) )/2.0 ), kic_res2_);

		pose.update_residue_neighbors();
		kinmover->apply(pose);
	}
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
DeclareBond::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data
)
{
	atom1_ = tag->getOption< std::string >( "atom1" );
	atom2_ = tag->getOption< std::string >( "atom2" );
	add_termini_ = tag->getOption< bool >( "add_termini", true );
	rebuild_fold_tree_ = tag->getOption< bool >( "rebuild_fold_tree", false );
	run_kic_ = tag->getOption< bool >( "run_KIC", false);

	kic_res1_ = tag->getOption< core::Size >( "KIC_res1", core::Size(0) );
	kic_res2_ = tag->getOption< core::Size >( "KIC_res2", core::Size(0) );

	if ( tag->hasOption( "res1_selector" ) ) {
		selector1_ = protocols::rosetta_scripts::parse_residue_selector( tag, data, "res1_selector" );
		res1_ = 0;
		runtime_assert_msg( ! tag->hasOption( "res1" ), "Please only use one of 'res1_selector' or 'res1', not both" );
	} else if ( tag->hasOption( "res1" ) ) {
		selector1_ = nullptr;
		res1_ = tag->getOption< core::Size >( "res1" );
	} else {
		utility_exit_with_message( "Please provide either res1 or res1_selector to DeclareBond" );
	}

	if ( tag->hasOption( "res2_selector" ) ) {
		selector2_ = protocols::rosetta_scripts::parse_residue_selector( tag, data, "res2_selector" );
		res2_ = 0;
		runtime_assert_msg( ! tag->hasOption( "res2" ), "Please only use one of 'res2_selector' or 'res2', not both" );
	} else if ( tag->hasOption( "res2" ) ) {
		selector2_ = nullptr;
		res2_ = tag->getOption< core::Size >( "res2" );
	} else {
		utility_exit_with_message( "Please provide either res2 or res2_selector to DeclareBond" );
	}

}

moves::MoverOP DeclareBond::clone() const { return utility::pointer::make_shared< DeclareBond >( *this ); }
moves::MoverOP DeclareBond::fresh_instance() const { return utility::pointer::make_shared< DeclareBond >(); }





std::string DeclareBond::get_name() const {
	return mover_name();
}

std::string DeclareBond::mover_name() {
	return "DeclareBond";
}

void DeclareBond::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "res1", xsct_non_negative_integer, "Residue containing first atom" )
		+ XMLSchemaAttribute( "res2", xsct_non_negative_integer, "Residue containing second atom" )
		+ XMLSchemaAttribute::required_attribute( "atom1", xs_string, "Name of first atom" )
		+ XMLSchemaAttribute::required_attribute( "atom2", xs_string, "Name of second atom" )
		+ XMLSchemaAttribute::attribute_w_default( "add_termini", xsct_rosetta_bool, "Add termini to pose?", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "rebuild_fold_tree", xsct_rosetta_bool, "Rebuild the fold tree after declaring this bond?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "run_kic", xsct_rosetta_bool, "Run KIC to close any chainbreak caused by the declared chemical bond?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "KIC_res1", xsct_non_negative_integer, "First residue to use in KIC", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "KIC_res2", xsct_non_negative_integer, "Second residue to use in KIC", "0" );

	core::select::residue_selector::attributes_for_parse_residue_selector(
		attlist, "res1_selector",
		"Alternative to using the res1 option. This residue selector must select exactly one residue!" );

	core::select::residue_selector::attributes_for_parse_residue_selector(
		attlist, "res2_selector",
		"Alternative to using the res2 option. This residue selector must select exactly one residue!" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Declares a chemical bond between two atoms", attlist );

}

std::string DeclareBondCreator::keyname() const {
	return DeclareBond::mover_name();
}

protocols::moves::MoverOP
DeclareBondCreator::create_mover() const {
	return utility::pointer::make_shared< DeclareBond >();
}

void DeclareBondCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DeclareBond::provide_xml_schema( xsd );
}


} // moves
} // protocols
