// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/CenHBPotential.cc
/// @brief  Smooth, differentiable version of centroid hbond term
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_CenHBPotential_fwd_hh
#define INCLUDED_core_scoring_CenHBPotential_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class CenHBPotential;

typedef utility::pointer::shared_ptr< CenHBPotential > CenHBPotentialOP;
typedef utility::pointer::shared_ptr< CenHBPotential const > CenHBPotentialCOP;

} // ns scoring
} // ns core

#endif
