// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/SCMinMultifunc.fwd.hh
/// @brief  Forward declaration of the interface class to the minimizer during sidechain minimization.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_scmin_SCMinMultifunc_fwd_hh
#define INCLUDED_core_pack_scmin_SCMinMultifunc_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace scmin {

class SCMinMultifunc;

typedef utility::pointer::shared_ptr< SCMinMultifunc > SCMinMultifuncOP;
typedef utility::pointer::shared_ptr< SCMinMultifunc > SCMinMultifuncCOP;

} // namespace scmin
} // namespace pack
} // namespace core

#endif
