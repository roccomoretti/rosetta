// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/BasePairStep.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_BasePairStep_HH
#define INCLUDED_core_pose_rna_BasePairStep_HH

#include <utility/VirtualBase.hh>
#include <core/pose/rna/BasePairStep.fwd.hh>
#include <core/types.hh>
#include <iosfwd>

namespace core {
namespace pose {
namespace rna {


///////////////////////////////////
//  5'-- ( i )     --   ( i+1) -- 3'
//         |               |
//         |               |
//  3'-- (j+q+1)         ( j ) -- 5'
//           \           /
//             n - n - n
//       allowing bulges (n's) on the second strand
///////////////////////////////////
class BasePairStep: public utility::VirtualBase {

public:

	//constructor
	// This is more information than is required -- i_next should be i+1, j_next should be j+q+1.
	//  but I wanted this to be explicit to prevent confusion for which numbers correspond to what.
	BasePairStep( Size const i, Size const i_next,
		Size const j, Size const j_next );

	//destructor
	~BasePairStep() override;

	Size const & i() const { return base_pair_step_.first.first; };
	Size const & i_next() const { return base_pair_step_.first.second; };
	Size const & j() const { return base_pair_step_.second.first; };
	Size const & j_next() const { return base_pair_step_.second.second; };

	friend
	std::ostream &
	operator <<( std::ostream & os, BasePairStep const & bps );

	friend
	bool operator == (BasePairStep const & lhs, BasePairStep const & rhs )
	{
		//There must be a more elegant way to do this...
		if ( lhs.base_pair_step_.first.first   == rhs.base_pair_step_.first.first &&
				lhs.base_pair_step_.first.second  == rhs.base_pair_step_.first.second &&
				lhs.base_pair_step_.second.first  == rhs.base_pair_step_.second.first &&
				lhs.base_pair_step_.second.second == rhs.base_pair_step_.second.second ) {
			return true;
		}
		return false;
	}

private:

	typedef std::pair< Size, Size > DinucleotideStrand;
	std::pair< DinucleotideStrand, DinucleotideStrand > base_pair_step_;

};

} //base_pairs
} //denovo
} //rna

#endif
