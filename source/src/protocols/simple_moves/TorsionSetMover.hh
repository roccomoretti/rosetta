// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/TorsionSetMover.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_simple_moves_TorsionSetMover_HH
#define INCLUDED_protocols_simple_moves_TorsionSetMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/TorsionSetMover.fwd.hh>
#include <core/id/TorsionID.hh>

#include <utility/vector1.hh> // AUTO IWYU For vector1

//////////////////////////////////////////////////////////////////////
// Bare bones mover -- just applied predefined set of torsions
//  to pose.
//////////////////////////////////////////////////////////////////////

namespace protocols {
namespace simple_moves {

class TorsionSetMover: public moves::Mover {

public:

	/// @brief Default constructor -- explicitly deleted.
	TorsionSetMover() = delete;

	/// @brief Constructor.
	TorsionSetMover(
		utility::vector1< core::id::TorsionID >  torsion_ids,
		utility::vector1< core::Real > torsion_values );

	/// @brief Destructor.
	~TorsionSetMover() override;

	/// @brief Clone function: create a copy of this object and return an owning pointer to the copy.
	protocols::moves::MoverOP clone() const override;

public:

	/// @brief Apply function: apply this mover to a pose.
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override { return "TorsionSetMover"; }

private:

	utility::vector1< core::id::TorsionID >  torsion_ids_;
	utility::vector1< core::Real > torsion_values_;

};

} //simple_moves
} //protocols

#endif
