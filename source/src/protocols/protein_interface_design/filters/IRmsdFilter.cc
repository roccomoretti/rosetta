// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/filters/RmsdFilter.cc
/// @brief rmsd filtering
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/IRmsdFilter.hh>
#include <protocols/protein_interface_design/filters/IRmsdFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/rosetta_scripts/util.hh>


#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

IRmsdFilter::IRmsdFilter() :
	protocols::filters::Filter( "IRmsd" ),
	threshold_( 5.0 ),
	reference_pose_( /* NULL */ )
{}

IRmsdFilter::IRmsdFilter(protocols::docking::DockJumps const movable_jumps,
	core::Real const threshold,
	core::pose::PoseOP reference_pose)
: protocols::filters::Filter( "IRmsd" ),
	threshold_(threshold),
	reference_pose_(std::move(reference_pose)),
	movable_jumps_(movable_jumps)
{}

IRmsdFilter::~IRmsdFilter() = default;

protocols::filters::FilterOP
IRmsdFilter::clone() const {
	return utility::pointer::make_shared< IRmsdFilter >( *this );
}

static basic::Tracer TR( "protocols.protein_interface_design.filters.IRmsdFilter" );
core::Real
IRmsdFilter::compute( core::pose::Pose const & pose ) const
{
	//ScoreFunctionOP temp = new ScoreFuntion(scorefxn_)
	return protocols::docking::calc_Irmsd(pose, *reference_pose_, scorefxn_, movable_jumps_);
}

bool
IRmsdFilter::apply( core::pose::Pose const & pose ) const {
	TR << "Beginning compute" << std::endl;
	core::Real const rmsd( compute( pose ));
	TR << "I_rmsd: " << rmsd ;
	if ( rmsd <= threshold_ ) {
		TR<<" passing."<<std::endl;
		return( true );
	} else TR<<" failing." << std::endl;
	return( false );
}

void
IRmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	out<<"RMSD: " << rmsd <<'\n';
}

core::Real
IRmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	return( (core::Real) rmsd );
}

void
IRmsdFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map )
{
	reference_pose_ = protocols::rosetta_scripts::legacy_saved_pose_or_input( tag, data_map, class_name() );

	threshold_ = tag->getOption<core::Real>( "threshold", 5.0 );

	//TODO: support multiple jumps
	auto jump_num = tag->getOption<core::Size>( "jump", 1);

	movable_jumps_.push_back(jump_num);

	scorefxn_ = rosetta_scripts::parse_score_function( tag, data_map )->clone();

	TR<<"Built IRmsdFilter with threshold " << threshold_ << ", scorefxn " <<
		rosetta_scripts::get_score_function_name(tag) <<", jump "<< jump_num << std::endl;

}



std::string IRmsdFilter::name() const {
	return class_name();
}

std::string IRmsdFilter::class_name() {
	return "IRmsd";
}

void IRmsdFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_saved_reference_pose( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "RMSD threshold above which we would fail the filter", "5" )
		+ XMLSchemaAttribute::attribute_w_default( "jump", xsct_non_negative_integer, "Jump for calculating docking RMSD, numbered sequentially from 1", "1" );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Calculates an interface rmsd. Rmsd is calculated over all backbone atoms for those residues found in the interface of the reference structure. Interface residues are those residues which are within 8 Angstroms of any residue on the other side of the interface.", attlist );
}

std::string IRmsdFilterCreator::keyname() const {
	return IRmsdFilter::class_name();
}

protocols::filters::FilterOP
IRmsdFilterCreator::create_filter() const {
	return utility::pointer::make_shared< IRmsdFilter >();
}

void IRmsdFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	IRmsdFilter::provide_xml_schema( xsd );
}



} // filters
} // protein_interface_design
} // devel


