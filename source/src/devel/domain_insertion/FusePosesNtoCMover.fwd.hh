// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/domain_insertion/FusePosesNtoCMover.fwd.hh
/// @brief  fwd hh file for FusePosesNtoCMover
/// @author Florian Richter, flosopher@gmail.com, february 2013

#ifndef INCLUDED_devel_domain_insertion_FusePosesNtoCMover_fwd_hh
#define INCLUDED_devel_domain_insertion_FusePosesNtoCMover_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace devel {
namespace domain_insertion {


class FusePosesNtoCMover;

typedef utility::pointer::shared_ptr< FusePosesNtoCMover > FusePosesNtoCMoverOP;
typedef utility::pointer::shared_ptr< FusePosesNtoCMover const > FusePosesNtoCMoverCOP;

class SetupCoiledCoilFoldTreeMover;
typedef utility::pointer::shared_ptr< SetupCoiledCoilFoldTreeMover > SetupCoiledCoilFoldTreeMoverOP;
typedef utility::pointer::shared_ptr< SetupCoiledCoilFoldTreeMover const > SetupCoiledCoilFoldTreeMoverCOP;

}
}

#endif
