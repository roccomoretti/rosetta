// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file NWAligner.hh
/// @brief class definition for a class that aligns two Sequence objects using a
/// the Needleman-Wunsch alignment algorithm.
/// @author James Thompson

#ifndef INCLUDED_core_sequence_NWAligner_hh
#define INCLUDED_core_sequence_NWAligner_hh

#include <core/types.hh>

#include <core/sequence/ScoringScheme.fwd.hh>

#include <core/sequence/Aligner.hh>
#include <core/sequence/DP_Matrix.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>

namespace core {
namespace sequence {

class NWAligner : public Aligner {

public:

	/// @brief ctor
	NWAligner() {}

	/// @brief dtor
	~NWAligner() override = default;


	SequenceAlignment align(
		SequenceOP seq_y,
		SequenceOP seq_x,
		ScoringSchemeOP ss
	) override;

	static void init_matrix(
		Size y_len,
		Size x_len,
		ScoringSchemeOP ss,
		DP_Matrix & scores
	);

}; // class NWAligner

} // sequence
} // core

#endif
