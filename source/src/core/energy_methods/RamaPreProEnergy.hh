// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/RamaPreProEnergy.hh
/// @brief  A variation on the Ramachandran scorefunction that has separate probability tables for residues that precede prolines.
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- modified this to work with D-amino acids, BACKBONE_AA amino acids, and cyclic geometry.

#ifndef INCLUDED_core_energy_methods_RamaPreProEnergy_hh
#define INCLUDED_core_energy_methods_RamaPreProEnergy_hh

// Unit headers
#include <core/energy_methods/RamaPreProEnergy.fwd.hh>
#include <core/scoring/RamaPrePro.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// C++ headers

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class RamaPreProEnergy : public core::scoring::methods::ContextIndependentLRTwoBodyEnergy {
public:
	typedef core::scoring::methods::ContextIndependentLRTwoBodyEnergy  parent;

public:
	RamaPreProEnergy( );

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size rsd1,
		Size rsd2
	) const override;


	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override { }

	utility::vector1< id::PartialAtomID >
	atoms_with_dof_derivatives( conformation::Residue const & res, pose::Pose const & pose ) const override;

	Real
	eval_intraresidue_dof_derivative(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const & /*min_data*/,
		id::DOF_ID const & /*dof_id*/,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & /*sfxn*/,
		core::scoring::EnergyMap const & weights
	) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return false; }

	bool
	defines_intrares_dof_derivatives( pose::Pose const & ) const override { return true; }

	virtual
	Distance
	atomic_interaction_cutoff() const { return 0.0; }

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override { }

	//fpd  use the new minimizer interface
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	core::scoring::methods::LongRangeEnergyType
	long_range_type() const override;

	bool is_allowed_type( core::chemical::ResidueType const & rt ) const;

private:
	core::scoring::RamaPrePro const & potential_;

	core::Size version() const override;

};


} // namespace energy_methods
} // namespace core


#endif // INCLUDED_core_energy_methods_RamaPreProEnergy_HH
