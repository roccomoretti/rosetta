// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/DownstreamAlgorithm.fwd.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_downstream_DownstreamAlgorithm_fwd_hh
#define INCLUDED_protocols_match_downstream_DownstreamAlgorithm_fwd_hh

// Package headers

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {
namespace downstream {

class DownstreamAlgorithm;
typedef utility::pointer::shared_ptr< DownstreamAlgorithm > DownstreamAlgorithmOP;
typedef utility::pointer::shared_ptr< DownstreamAlgorithm const > DownstreamAlgorithmCOP;

}
}
}

#endif
