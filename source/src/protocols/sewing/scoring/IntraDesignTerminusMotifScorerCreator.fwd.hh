// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   RequirementCreator.hh
/// @brief  Forward Header for base class for RequirementCreators for the Requirement load-time factory registration scheme
/// @author Tim Jacobs


#ifndef INCLUDED_protocols_sewing_scoring_IntraDesignTerminusMotifScorerCreator_fwd_hh
#define INCLUDED_protocols_sewing_scoring_IntraDesignTerminusMotifScorerCreator_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

/// @brief Abstract base class for a Requirement factory; the Creator class is responsible for
/// creating a particular Requirement class.
class IntraDesignTerminusMotifScorerCreator;
typedef utility::pointer::shared_ptr< IntraDesignTerminusMotifScorerCreator > IntraDesignTerminusMotifScorerCreatorOP;
typedef utility::pointer::shared_ptr< IntraDesignTerminusMotifScorerCreator const > IntraDesignTerminusMotifScorerCreatorCOP;

} //namespace
} //namespace
} //namespace

#endif
