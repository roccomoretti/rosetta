// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/ContactMolecularSurfaceFilter.cc
/// @brief  Filter structures by weighted contact molecular surface area
/// @author Longxing Cao (longxing@uw.edu)

// Unit Headers
#include <protocols/simple_filters/ContactMolecularSurfaceFilter.hh>
#include <protocols/simple_filters/ContactMolecularSurfaceFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueVector.hh>
#include <core/select/residue_selector/util.hh>
#include <protocols/simple_filters/ShapeComplementarityFilter.hh>

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


//// C++ headers
static basic::Tracer tr( "protocols.simple_filters.ContactMolecularSurfaceFilter" );

namespace protocols {
namespace simple_filters {

using core::Size;

// @brief default constructor
ContactMolecularSurfaceFilter::ContactMolecularSurfaceFilter():
	Filter( "ContactMolecularSurface" ),
	filtered_area_( 250 ),
	distance_weight_( 1 ),
	near_squared_size_( 0 ),
	apolar_target_( false ),
	quick_( false ),
	verbose_( false ),
	use_rosetta_radii_( false ),
	selector1_(),
	selector2_()
{}


// @brief constructor with arguments
ContactMolecularSurfaceFilter::ContactMolecularSurfaceFilter( Real const & filtered_area, Real const & distance_weight,
	bool const quick, bool const verbose, bool const use_rosetta_radii ):
	Filter( "ContactMolecularSurface" ),
	filtered_area_( filtered_area ),
	distance_weight_( distance_weight ),
	near_squared_size_( 0 ),
	apolar_target_( false ),
	quick_( quick ),
	verbose_( verbose ),
	use_rosetta_radii_( use_rosetta_radii ),
	selector1_(),
	selector2_()
{}


// @brief constructor with arguments and selectors
ContactMolecularSurfaceFilter::ContactMolecularSurfaceFilter( Real const & filtered_area, Real const & distance_weight,
	bool const quick, bool const verbose,
	core::select::residue_selector::ResidueSelectorCOP const selector1,
	core::select::residue_selector::ResidueSelectorCOP const selector2,
	bool const use_rosetta_radii ):
	Filter( "ContactMolecularSurface" ),
	filtered_area_( filtered_area ),
	distance_weight_( distance_weight ),
	near_squared_size_( 0 ),
	apolar_target_( false ),
	quick_( quick ),
	verbose_( verbose ),
	use_rosetta_radii_( use_rosetta_radii ),
	selector1_( selector1 ),
	selector2_( selector2 )
{}


// accessors
bool ContactMolecularSurfaceFilter::quick() const { return quick_; }
bool ContactMolecularSurfaceFilter::verbose() const { return verbose_; }
bool ContactMolecularSurfaceFilter::use_rosetta_radii() const { return use_rosetta_radii_; }
bool ContactMolecularSurfaceFilter::apolar_target() const { return apolar_target_; }
core::Real ContactMolecularSurfaceFilter::filtered_area() const { return filtered_area_; }
core::Real ContactMolecularSurfaceFilter::distance_weight() const { return distance_weight_; }
core::Real ContactMolecularSurfaceFilter::near_squared_size() const { return near_squared_size_; }
core::select::residue_selector::ResidueSelectorCOP ContactMolecularSurfaceFilter::selector1() const { return selector1_; }
core::select::residue_selector::ResidueSelectorCOP ContactMolecularSurfaceFilter::selector2() const { return selector2_; }

// mutators
void ContactMolecularSurfaceFilter::quick( bool const quick ) { quick_ = quick; }
void ContactMolecularSurfaceFilter::verbose( bool const verbose ) { verbose_ = verbose; }
void ContactMolecularSurfaceFilter::use_rosetta_radii( bool const use_rosetta_radii ) { use_rosetta_radii_ = use_rosetta_radii; }
void ContactMolecularSurfaceFilter::apolar_target( bool const apolar_target ) { apolar_target_ = apolar_target; }
void ContactMolecularSurfaceFilter::filtered_area( core::Real const filtered_area ) { filtered_area_ = filtered_area; }
void ContactMolecularSurfaceFilter::distance_weight( core::Real const distance_weight ) { distance_weight_ = distance_weight; }
void ContactMolecularSurfaceFilter::near_squared_size( core::Real const near_squared_size ) { near_squared_size_ = near_squared_size; }
void ContactMolecularSurfaceFilter::selector1( core::select::residue_selector::ResidueSelectorCOP const selector ) { selector1_ = selector; }
void ContactMolecularSurfaceFilter::selector2( core::select::residue_selector::ResidueSelectorCOP const selector ) { selector2_ = selector; }

filters::FilterOP
ContactMolecularSurfaceFilter::clone() const {
	return utility::pointer::make_shared<ContactMolecularSurfaceFilter>( *this );
}

filters::FilterOP
ContactMolecularSurfaceFilter::fresh_instance() const {
	return utility::pointer::make_shared<ContactMolecularSurfaceFilter>();
}


core::Real
ContactMolecularSurfaceFilter::compute( Pose const & pose ) const
{
	ContactMolecularSurfaceCalculator scc;

	if ( !scc.Init() ) {
		throw CREATE_EXCEPTION(EXCN_InitFailed, "");
	}

	scc.settings.weight = distance_weight_;
	scc.settings.near_squared_size = near_squared_size_;
	if ( quick_ ) {
		scc.settings.density = 5.0;
	}
	if ( use_rosetta_radii_ ) {
		scc.settings.use_rosetta_radii = true;
	}
	scc.Reset(); // this may not be needed anymore, but I'm leaving it here for safety

	core::Real result( 0.0 );

	if ( selector1_ && selector2_ ) {
		// selector-based
		setup_from_selectors( pose, scc );
		result = scc.CalcContactArea();
		if (  -1 == result ) {
			throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
		}
	} else {
		// jump-based
		throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
	}
	if ( verbose_ ) tr << "Weighted contact molecular surface is " << result;

	return result;
}

/// @brief
core::Real
ContactMolecularSurfaceFilter::report_sm( Pose const & pose ) const
{
	core::Real r;
	try {
		r = compute( pose );
	} catch( EXCN_InitFailed const & ) {
		tr.Error << "Issue initializing shape complementarity calculator - returning -2 instead." << std::endl;
		return -2;
	} catch( EXCN_CalcFailed const & ) {
		tr.Error << "Issue running shape complementarity calculator - returning -1 instead." << std::endl;
		return -1;
	}
	return r;
	//return r.sc;
}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose has high enough shape
// complementarity.
bool ContactMolecularSurfaceFilter::apply( Pose const & pose ) const
{
	core::Real r;
	try {
		r = compute( pose );
	} catch( EXCN_InitFailed const & ) {
		tr.Error << "Issue initializing shape complementarity calculator - failing filter." << std::endl;
		return false;
	} catch( EXCN_CalcFailed const & ) {
		tr.Error << "Issue running shape complementarity calculator - failing filter." << std::endl;
		return false;
	}

	if ( r < filtered_area_ ) {
		tr << "Filter failed current < threshold interface area: " << r << " < " << filtered_area_ << std::endl;
		return false;
	}

	tr << "Successfully filtered: " << r << std::endl;
	return true;
} // apply_filter

/// @brief parse xml
void
ContactMolecularSurfaceFilter::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data
) {
	filtered_area_ = tag->getOption<Real>( "min_interface", 1.0 );
	distance_weight_ = tag->getOption<Real>( "distance_weight", 1.0 );
	near_squared_size_ = tag->getOption<Real>( "near_squared_size", 0.0 );
	verbose_ = tag->getOption<bool>( "verbose", false );
	quick_ = tag->getOption<bool>( "quick", false );
	use_rosetta_radii_ = tag->getOption<bool>( "use_rosetta_radii", false );
	apolar_target_ = tag->getOption<bool>( "apolar_target", false );

	std::string const selector1name = tag->getOption< std::string >( "target_selector", "" );
	if ( !selector1name.empty() ) selector1_ = core::select::residue_selector::get_residue_selector( selector1name, data );

	std::string const selector2name = tag->getOption< std::string >( "binder_selector", "" );
	if ( !selector2name.empty() ) selector2_ = core::select::residue_selector::get_residue_selector( selector2name, data );


	if ( quick_ ) {
		tr.Info << "Calculating contact moleulcar surface in quick mode with less accuracy." << std::endl;
	}

	if ( !selector1_ || !selector2_ ) {
		tr.Info << "Target and binder are not properly set by the residue selectors, Quit!!!" << std::endl;
		utility_exit_with_message( "Fail to parse residue selectors for ContactMolecular surface." );
	}
}

/// @brief Uses residue selectors to set up the ContactMolecularSurfaceCalculator
/// @param[in]  pose Pose to be analyzed
/// @param[out] scc Initialized, empty ContactMolecularSurfaceCalculator, to which pose residues are added
void
ContactMolecularSurfaceFilter::setup_from_selectors( Pose const & pose, ContactMolecularSurfaceCalculator & scc ) const
{
	using core::select::residue_selector::ResidueVector;

	ResidueVector const residues1( selector1_->apply( pose ) );
	ResidueVector const residues2( selector2_->apply( pose ) );

	// Dump information about residues
	if ( tr.Info.visible() ) {
		tr.Info << "Using residues for molecule surface (rosetta numbering):" << std::endl;
		tr.Info << "  Target Surface: ";
		for ( auto r = residues1.begin(); r != residues1.end(); ++r ) {
			tr.Info << (r == residues1.begin() ? "" : ", ") << *r;
		}
		tr.Info << std::endl;
		tr.Info << "  Binder Surface: ";
		for ( auto r = residues2.begin(); r != residues2.end(); ++r ) {
			tr.Info << (r == residues2.begin() ? "" : ", ") << *r;
		}
		tr.Info << std::endl;
	}

	for ( core::Size r : residues1 ) {
		scc.AddResidue( 0, pose.residue( r ), apolar_target_ );
	}

	for ( core::Size r : residues2 ) {
		scc.AddResidue( 1, pose.residue( r ) );
	}
}


std::string ContactMolecularSurfaceFilter::name() const {
	return class_name();
}

std::string ContactMolecularSurfaceFilter::class_name() {
	return "ContactMolecularSurface";
}

void ContactMolecularSurfaceFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "min_interface" , xsct_real , "The filter fails is the calculated interface area is less than the given value." , "0" )
		+ XMLSchemaAttribute::attribute_w_default( "distance_weight" , xsct_real , "The weight factor of the cloest distance betweent the distance that is multiplied by the area by each surface dot." , "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "near_squared_size" , xsct_real , "The weight factor of the cloest distance betweent the distance that is multiplied by the area by each surface dot." , "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "verbose" , xsct_rosetta_bool , "If true, print extra calculation details to the tracer." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "quick" , xsct_rosetta_bool , "If true, do a quicker, less accurate calculation by reducing the density." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "use_rosetta_radii" , xsct_rosetta_bool , "Use rosetta radii instead of the historical shape complementarity radii." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "apolar_target" , xsct_rosetta_bool , "Only look at non-polar atoms on the target side." , "false" )
		+ XMLSchemaAttribute( "target_selector" , xs_string , "Explicitly set which residues are on each side of the interface using residue_selectors." )
		+ XMLSchemaAttribute( "binder_selector" , xs_string , "Explicitly set which residues are on each side of the interface using residue_selectors." ) ;

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Calculates the contact molecular surface area on the target defined by the target_selector, and the area is weighted by the cloest distance between the target and the binder. The binder is defined by the binder residue selector", attlist );
}

std::string ContactMolecularSurfaceFilterCreator::keyname() const {
	return ContactMolecularSurfaceFilter::class_name();
}

protocols::filters::FilterOP
ContactMolecularSurfaceFilterCreator::create_filter() const {
	return utility::pointer::make_shared<ContactMolecularSurfaceFilter>();
}

void ContactMolecularSurfaceFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ContactMolecularSurfaceFilter::provide_xml_schema( xsd );
}


} // filters
} // protocols
