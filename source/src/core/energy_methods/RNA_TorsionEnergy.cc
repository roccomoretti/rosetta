// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RNA_TorsionEnergy.cc
/// @brief  RNA_Torsion energy method class implementation
/// @author Rhiju Das

// Unit Headers
#include <core/energy_methods/RNA_TorsionEnergy.hh>
#include <core/energy_methods/RNA_TorsionEnergyCreator.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/conformation/Residue.hh>

// Utility headers

#include <core/scoring/EnergyMap.hh>
#include <utility>
#include <utility/vector1.hh>

// C++
using namespace core::chemical;

namespace core {
namespace energy_methods {


/// @details This must return a fresh instance of the RNA_TorsionEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
RNA_TorsionEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< RNA_TorsionEnergy >( options.rna_options() );
}

core::scoring::ScoreTypes
RNA_TorsionEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( rna_torsion );
	sts.push_back( rna_torsion_sc );
	return sts;
}


/// ctor
RNA_TorsionEnergy::RNA_TorsionEnergy( core::scoring::rna::RNA_EnergyMethodOptions const & options,
	core::scoring::rna::RNA_TorsionPotentialOP rna_torsion_potential /* = 0 */ ) :
	parent( utility::pointer::make_shared< RNA_TorsionEnergyCreator >() ),
	options_( options ),
	rna_torsion_potential_(std::move( rna_torsion_potential ))
{
	if ( rna_torsion_potential_ == nullptr ) rna_torsion_potential_ = utility::pointer::make_shared< core::scoring::rna::RNA_TorsionPotential >( options );
}

/// clone
core::scoring::methods::EnergyMethodOP
RNA_TorsionEnergy::clone() const
{
	return utility::pointer::make_shared< RNA_TorsionEnergy >( options_, rna_torsion_potential_ );
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {
	if ( rsd1.has_variant_type( REPLONLY ) ) return;
	if ( rsd2.has_variant_type( REPLONLY ) ) return;

	emap[ core::scoring::rna_torsion ] += rna_torsion_potential_->residue_pair_energy( rsd1, rsd2, pose );
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const {
	if ( rsd.has_variant_type( REPLONLY ) ) return;

	emap[ core::scoring::rna_torsion ]    += rna_torsion_potential_->eval_intrares_energy( rsd, pose );
	emap[ core::scoring::rna_torsion_sc ] += rna_torsion_potential_->intrares_side_chain_score(); // evaluated at same time as above.
}


///////////////////////////////////////////////////////////////////////////////
void
RNA_TorsionEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	rna_torsion_potential_->eval_atom_derivative( id, pose, weights, F1, F2 );
}


void RNA_TorsionEnergy::indicate_required_context_graphs( utility::vector1< bool > & /*context_graphs_required */ ) const {}

/// @brief RNA_PairwiseLowResolutionEnergy distance cutoff
Distance
RNA_TorsionEnergy::atomic_interaction_cutoff() const
{
	return 0.0; /// Uh, I don't know.
}

core::Size
RNA_TorsionEnergy::version() const
{
	return 2; // A new torsion potential (integration from Das lab branch -- Aug 2011)
}


} //scoring
} //core

