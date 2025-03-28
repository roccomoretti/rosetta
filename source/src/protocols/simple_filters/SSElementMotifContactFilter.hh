// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/SSElementMotifContactFilter.hh
/// @brief  Reports either the average ss-contact or the worst ss-contacting element
/// @author TJ Brunette (tjbrunette@gmail.com)


#ifndef INCLUDED_protocols_simple_filters_SSElementMotifContactFilter_hh
#define INCLUDED_protocols_simple_filters_SSElementMotifContactFilter_hh

// Unit Headers
#include <protocols/simple_filters/SSElementMotifContactFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/scoring/motif/motif_hash_stuff.fwd.hh>
#include <protocols/loops/Loops.hh>
// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <set>

//// C++ headers

namespace protocols {
namespace simple_filters {

class SSElementMotifContactFilter : public protocols::filters::Filter{
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef std::string String;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;


public:// constructor/destructor


	// @brief default constructor
	SSElementMotifContactFilter();

	// @brief copy constructor
	SSElementMotifContactFilter( SSElementMotifContactFilter const & rval );

	~SSElementMotifContactFilter() override;


public:// virtual constructor


	// @brief make clone
	filters::FilterOP clone() const override { return utility::pointer::make_shared< SSElementMotifContactFilter >(*this);}
	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return utility::pointer::make_shared< SSElementMotifContactFilter >();}


public:// mutator


	// @brief
	void filtered_value( Real const & value );


public:// accessor


	// @brief get name of this filter


public:// virtual main operation


	Real report_sm(const Pose & pose ) const override;
	void report( std::ostream & out,const Pose & pose ) const override;
	protocols::loops::Loops get_ss_elements(const Pose & pose) const;
	core::Size which_ssElement(core::Size res,protocols::loops::Loops ssElements) const;
	core::Size get_SSelements_in_contact(core::Size element,protocols::loops::Loops ssElements, const Pose & pose) const;
	core::Size get_ssElements_in_contact_w_threshold(std::multiset<core::Size> ssElements_in_contact) const;
	Real compute( const Pose & pose ) const;
	bool apply(const Pose & pose ) const override;


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap &
	) override;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:

	Real filtered_value_;
	core::Size ignore_n_terminal_SS_;
	core::Size ignore_c_terminal_SS_;
	bool only_n_term_;
	bool only_c_term_;
	Real threshold_;
	core::Size contacts_between_ssElement_threshold_;
	bool report_avg_;
	bool use_atomic_contact_filter_; //normally use motifs
	core::scoring::motif::MotifHashManager *mman_;
};

} // filters
} // protocols

#endif
