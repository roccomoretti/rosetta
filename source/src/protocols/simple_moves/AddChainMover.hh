// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AddChainMover.hh
/// @brief switch the chain order

#ifndef INCLUDED_protocols_simple_moves_AddChainMover_hh
#define INCLUDED_protocols_simple_moves_AddChainMover_hh

#include <protocols/simple_moves/AddChainMover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// C++ Headers
namespace protocols {
namespace simple_moves {

class AddChainMover : public moves::Mover {
public:
	AddChainMover();

	void apply( core::pose::Pose & pose ) override;

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	void fname( std::string const & f ){ fname_ = f; }
	std::string fname() const{ return fname_; }

	void new_chain( bool const n ){ new_chain_ = n; }
	bool new_chain() const{ return new_chain_; }

	void swap_chain_number( core::Size const wc ){ swap_chain_number_ = wc; }
	core::Size swap_chain_number() const { return swap_chain_number_; }

	void scorefxn( core::scoring::ScoreFunctionOP s );
	core::scoring::ScoreFunctionOP scorefxn() const;

	bool random_access() const{ return random_access_; }
	void random_access( bool const b ){ random_access_ = b; }

	void update_PDBInfo( bool const n ){ update_PDBInfo_ = n; }
	bool update_PDBInfo() const{ return update_PDBInfo_; }

	void add_new_chain( core::pose::Pose & pose ) const; // Adds new chain to pose
	void swap_chain( core::pose::Pose & pose ) const; // Adds new chain to pose

	// There's no point in a setter for this. It has to happen during parse_my_tag
	// void spm_reference_name( std::string const & name ) { spm_reference_name_ = name; }
	std::string spm_reference_name( ) const { return spm_reference_name_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	void
	load_pose( core::pose::Pose & pose, std::string & name ) const;

private:
	std::string fname_; //dflt ""; pdb names to load (can accept a comma-separated list)
	bool new_chain_; //dflt true; add as a new chain?
	bool random_access_; //dflt false; if true randomly choose one file name from a list and work with that throughout the run.
	bool update_PDBInfo_; //dflt true; update chain ids.
	core::Size swap_chain_number_; //dflt 2; swap chain with specified chain number
	core::scoring::ScoreFunctionOP scorefxn_; //dflt score12; used to score the new pose
	std::string spm_reference_name_; // This is only used for display and adding to the pose comments
	core::pose::PoseCOP spm_reference_pose_;

};


} // simple_moves
} // protocols

#endif
