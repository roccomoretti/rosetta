// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/protein_interface_design/filters/DisulfideFilter.hh
/// @brief Filters for interfaces which could form a disulfide bond between
/// docking partners.
/// @author Sarel Fleishman (sarelf@uw.edu)

#include <protocols/simple_filters/AtomicContactFilter.hh>
#include <protocols/simple_filters/AtomicContactFilterCreator.hh>


// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

//parsing
#include <utility/tag/Tag.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/selection.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueSpanSelector.hh>
#include <core/select/util.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.filters.AtomicContactFilter" );

/// @brief default ctor
AtomicContactFilter::AtomicContactFilter() :
	parent( "AtomicContact" ),
	protocols::moves::ResId( 0u )
{}

/// @brief Constructor with a single target residue
AtomicContactFilter::AtomicContactFilter( core::Size const res1, core::Size const res2, core::Real const distance, bool const sidechain, bool const backbone, bool const protons ) :
	parent( "AtomicContact" ),
	protocols::moves::ResId( res2 ),
	side1_( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( res1 ) ),
	distance_( distance ),
	sidechain_( sidechain ),
	backbone_( backbone ),
	protons_( protons )
{}

/// @return Whether a disulfide bond is possible between any of the targets
bool AtomicContactFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	report( TR.Debug, pose );
	if ( dist <= distance_ ) return true;
	return false;
}

core::Real
AtomicContactFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core::conformation;

	if ( !get_resid(pose) ) {
		TR.Error << "residue2 has not been defined"<<std::endl;
		runtime_assert( get_resid(pose) );
	}
	core::Real nearest_distance( 10000 );
	Residue const res2( pose.residue( get_resid(pose) ) );
	debug_assert( side1_ );
	for ( core::Size residue1 : core::select::get_residues_from_subset( side1_->apply( pose ) ) ) {
		Residue const res1( pose.residue( residue1 ) );

		auto atom1_begin( res1.atom_begin() ), atom1_end( res1.atom_end() ), atom2_begin( res2.atom_begin() ), atom2_end( res2.atom_end() );
		if ( sidechain_ && !backbone_ ) {
			atom1_begin = res1.sidechainAtoms_begin();
			atom2_begin = res2.sidechainAtoms_begin();
		}
		if ( !sidechain_ && backbone_ ) {
			atom1_end = res1.sidechainAtoms_begin();
			atom2_end = res2.sidechainAtoms_begin();
		}
		if ( !protons_ ) {
			atom1_end = res1.heavyAtoms_end();
			atom2_end = res2.heavyAtoms_end();
		}
		for ( auto atom1=atom1_begin; atom1!=atom1_end; ++atom1 ) {
			for ( auto atom2=atom2_begin; atom2!=atom2_end; ++atom2 ) {
				core::Real const dist( atom1->xyz().distance( atom2->xyz() ) );
				if ( dist <= nearest_distance ) nearest_distance = dist;
			}//foreach atom2
		}//foreach atom1
	}
	return( nearest_distance );
}

core::Real
AtomicContactFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	return( dist );
}

void AtomicContactFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	core::Real const dist( compute( pose ) );
	out<<"Minimal distance between residues is "<<dist<<std::endl; // Unfortunately, no show() on ResidueSelectors
}

void AtomicContactFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	distance_ = tag->getOption< core::Real >( "distance", 4.0 );
	if ( tag->hasOption("range1") ) {
		std::istringstream range_str( tag->getOption< std::string >( "range1" ) );
		core::Size num1, num2;
		range_str >> num1 >> num2;
		if ( !range_str.good() ) TR << "cannot read parameter range1" << std::endl;
		side1_ = utility::pointer::make_shared< core::select::residue_selector::ResidueSpanSelector >( num1, num2 );
	} else {
		std::string const res1( tag->getOption< std::string >( "residue1" ) );
		side1_ = utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( res1 );
	}
	if ( tag->hasOption( "residue2" ) ) {
		std::string const res2( tag->getOption< std::string >( "residue2" ) );
		set_resid( core::pose::parse_resnum( res2 ) );
		modifiable( false );
	} else {
		modifiable( true );
		TR<<"AtomicContact: residue2 was not defined. A mover/filter will have to set it during the protocol\n";
	}
	sidechain_ = tag->getOption< bool >( "sidechain", true );
	backbone_  = tag->getOption< bool >( "backbone",  false );
	protons_   = tag->getOption< bool >( "protons",   false );

	TR<<"AtomicContact filter between residues with distance cutoff of "<<distance_<<std::endl; // Unfortunately, no show() on ResidueSelectors
}



std::string AtomicContactFilter::name() const {
	return class_name();
}

std::string AtomicContactFilter::class_name() {
	return "AtomicContact";
}

void AtomicContactFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "distance", xsct_real, "Distance below which a contact counts", "4.0" )
		+ XMLSchemaAttribute( "range1", xs_string, "Space-separated pair of integers indicating the first residue range over which to evaluate the filter" )
		+ XMLSchemaAttribute( "residue1", xsct_refpose_enabled_residue_number, "Residue number that may be making atomic contacts with residue2; only read if range1 is not provided" )
		+ XMLSchemaAttribute( "residue2", xsct_refpose_enabled_residue_number, "Residue number that may be making atomic contacts with the range of residues range1" )
		+ XMLSchemaAttribute( "sidechain", xsct_rosetta_bool, "Count sidechain contacts" )
		+ XMLSchemaAttribute( "backbone", xsct_rosetta_bool, "Count contacts to the backbone" )
		+ XMLSchemaAttribute( "protons", xsct_rosetta_bool, "Count contacts made by protons" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string AtomicContactFilterCreator::keyname() const {
	return AtomicContactFilter::class_name();
}

protocols::filters::FilterOP
AtomicContactFilterCreator::create_filter() const {
	return utility::pointer::make_shared< AtomicContactFilter >();
}

void AtomicContactFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AtomicContactFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
