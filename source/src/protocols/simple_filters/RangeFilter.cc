// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/RangeFilter.cc
/// @brief
/// @details
/// @author Javier Castellanos (javiercv@uw.edu)

// Unit Headers
#include <protocols/simple_filters/RangeFilter.hh>
#include <protocols/simple_filters/RangeFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility>
#include <utility/tag/Tag.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

//// C++ headers
static basic::Tracer tr( "protocols.filters.RangeFilter" );

namespace protocols {
namespace simple_filters {

// @brief default constructor
RangeFilter::RangeFilter() : protocols::filters::Filter()
{}

RangeFilter::RangeFilter(Real lower_bound, Real upper_bound, FilterOP  filter ) : protocols::filters::Filter(),
	filter_(std::move(filter)),
	lower_bound_( lower_bound ),
	upper_bound_( upper_bound )
{}


/// @brief
void
RangeFilter::report( std::ostream & out, Pose const & pose ) const
{
	Real value = filter_->apply( pose );
	out << value << " in range " << lower_bound_ << " - " << upper_bound_ << std::endl;
}


// @brief returns true if the given pose passes the filter, false otherwise.
bool RangeFilter::apply( Pose const & pose ) const
{
	Real value = filter_->report_sm( pose );
	if ( value > lower_bound_ && value < upper_bound_ ) {
		tr << "Successfully filtered: " << value << " in range " << lower_bound_ << " - " << upper_bound_ << std::endl;
		return true;
	} else {
		tr << "Filter failed: value = " << value << " range = "<< lower_bound_ << " - " << upper_bound_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
RangeFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data
)
{
	std::string const filter_name( tag->getOption< std::string >( "filter") );
	filter_ = protocols::rosetta_scripts::parse_filter( filter_name, data );

	lower_bound_ = tag->getOption<Real>( "lower_bound");
	upper_bound_ = tag->getOption<Real>( "upper_bound");
	debug_assert(lower_bound_ < upper_bound_);
}



std::string RangeFilter::name() const {
	return class_name();
}

std::string RangeFilter::class_name() {
	return "Range";
}

void RangeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"filter", xs_string,
		"Filter to be examined")
		+ XMLSchemaAttribute::required_attribute(
		"lower_bound", xsct_real,
		"minimal filter score allowed")
		+ XMLSchemaAttribute::required_attribute(
		"upper_bound", xsct_real,
		"maximal filter score allowed");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"This filter returns true if the return value of the given filter is between lower_bound and upper_bound.",
		attlist );
}

std::string RangeFilterCreator::keyname() const {
	return RangeFilter::class_name();
}

protocols::filters::FilterOP
RangeFilterCreator::create_filter() const {
	return utility::pointer::make_shared< RangeFilter >();
}

void RangeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RangeFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
