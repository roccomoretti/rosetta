// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/SetResidueAliasMover.hh
/// @brief Sets a residue alias in the StructureData
/// @author Tom Linsky (tlinsky@gmail.com)

#ifndef INCLUDED_protocols_denovo_design_movers_SetResidueAliasMover_hh
#define INCLUDED_protocols_denovo_design_movers_SetResidueAliasMover_hh

// Unit headers
#include <protocols/denovo_design/movers/SetResidueAliasMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

///@brief Sets a residue alias in the StructureData
class SetResidueAliasMover : public protocols::moves::Mover {

public:

	SetResidueAliasMover();

	// copy constructor (not needed unless you need deep copies)
	//SetResidueAliasMover( SetResidueAliasMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	~SetResidueAliasMover() override;


public:
	// mover virtual API
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;


	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//SetResidueAliasMover & operator=( SetResidueAliasMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Sets the name of the alias
	void
	set_alias_name( std::string const & alias ) { alias_name_ = alias; }

	/// @brief Sets the name of the segment that contains the residue to which the alias refers
	void
	set_segment_name( std::string const & segment ) { segment_name_ = segment; }

	/// @brief Sets the residue within the segment to which the alias refers
	void
	set_resid( core::Size const resid ) { resid_ = resid; }

private:
	std::string alias_name_;
	std::string segment_name_;
	core::Size resid_;

};

std::ostream &
operator<<( std::ostream & os, SetResidueAliasMover const & mover );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_SetResidueAliasMover_hh
