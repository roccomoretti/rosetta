// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/ShortestPathInFoldTree.fwd.hh
/// @brief  kinematics::ShortestPathInFoldTree forward declarations header
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_KinematicTaskControl_fwd_hh
#define INCLUDED_protocols_abinitio_KinematicTaskControl_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace abinitio {

// Forward
class KinematicTaskControl;

// Types
typedef  utility::pointer::shared_ptr< KinematicTaskControl >  KinematicTaskControlOP;
typedef  utility::pointer::shared_ptr< KinematicTaskControl const >  KinematicTaskControlCOP;

// Forward
class LoopJumpFoldCst;

// Types
typedef  utility::pointer::shared_ptr< LoopJumpFoldCst >  LoopJumpFoldCstOP;
typedef  utility::pointer::shared_ptr< LoopJumpFoldCst const >  LoopJumpFoldCstCOP;

} // namespace kinematics
} // namespace core

#endif
