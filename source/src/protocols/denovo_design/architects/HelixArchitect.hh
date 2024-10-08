// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/HelixArchitect.hh
/// @brief Architect for helices
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_denovo_design_architects_HelixArchitect_hh
#define INCLUDED_protocols_denovo_design_architects_HelixArchitect_hh

// Unit headers
#include <protocols/denovo_design/architects/HelixArchitect.fwd.hh>
#include <protocols/denovo_design/architects/DeNovoArchitect.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>

#include <protocols/denovo_design/components/Segment.fwd.hh> // AUTO IWYU For SegmentCOPs

// Utility headers

namespace protocols {
namespace denovo_design {
namespace architects {

///@brief Architect for helices
class HelixArchitect : public DeNovoArchitect {
public:
	HelixArchitect( std::string const & id_value );

	~HelixArchitect() override;

	static std::string
	class_name() { return "HelixArchitect"; }

	std::string
	type() const override;

	DeNovoArchitectOP
	clone() const override;

	StructureDataOP
	design( core::pose::Pose const & pose, core::Real & random ) const override;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

public:
	void
	set_lengths( std::string const & lengths_str );

	void
	set_lengths( Lengths const & lengths );

	/// @brief Returns the list of allowed segments; One of these will be chosen during design
	components::SegmentCOPs const &
	motifs() const;

	components::SegmentCOPs::const_iterator
	motifs_begin() const;

	components::SegmentCOPs::const_iterator
	motifs_end() const;

private:
	components::SegmentCOPs motifs_;
	Lengths lengths_;

};

} //protocols
} //denovo_design
} //architects

#endif //INCLUDED_protocols_denovo_design_architects_HelixArchitect_hh

