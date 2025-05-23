// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/electron_density_atomwise/ElecDensAtomwiseEnergy.fwd.hh
/// @brief  elec_dens_atomwise energy method forward declaration
/// @author Fang-Chieh Chou


#ifndef INCLUDED_core_scoring_electron_density_atomwise_ElecDensAtomwiseEnergy_FWD_HH
#define INCLUDED_core_scoring_electron_density_atomwise_ElecDensAtomwiseEnergy_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace energy_methods {

class ElecDensAtomwiseEnergy;

typedef utility::pointer::shared_ptr< ElecDensAtomwiseEnergy > ElecDensAtomwiseEnergyOP;
typedef utility::pointer::shared_ptr< ElecDensAtomwiseEnergy const > ElecDensAtomwiseEnergyCOP;

} // scoring
} // core


#endif
