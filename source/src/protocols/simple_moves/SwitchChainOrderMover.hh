// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SwitchChainOrderMover.hh
/// @brief switch the chain order

#ifndef INCLUDED_protocols_simple_moves_SwitchChainOrderMover_hh
#define INCLUDED_protocols_simple_moves_SwitchChainOrderMover_hh

#include <protocols/simple_moves/SwitchChainOrderMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMapObj.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

class SwitchChainOrderMover : public moves::Mover {
public:
	SwitchChainOrderMover();
	// Undefinded, commenting out to fix PyRosetta build  SwitchChainOrderMover( std::string const & );

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void chain_order( std::string const & co );
	utility::vector1< core::Size > chain_order( core::pose::Pose const & ) const;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP s );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	utility::vector1< std::string > chain_ids_;
	bool invert_chains_ = false;

	utility::pointer::shared_ptr< basic::datacache::DataMapObj< utility::vector1< core::Size > > > residue_numbers_; /// dflt NULL; a vector of residue numbers placed on the basic::datacache::DataMap which need to be changed due to the chain order switch
	core::scoring::ScoreFunctionOP scorefxn_; // dflt NULL; needed to score the pose for the energy graph to behave well
};


} // simple_moves
} // protocols

#endif
