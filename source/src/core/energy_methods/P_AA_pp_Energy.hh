// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/P_AA_pp_Energy.hh
/// @brief  Probability of observing an amino acid, given its phi/psi energy method declaration
/// @author Andrew Leaver-Fay


#ifndef INCLUDED_core_energy_methods_P_AA_pp_Energy_hh
#define INCLUDED_core_energy_methods_P_AA_pp_Energy_hh

// Unit headers
#include <core/energy_methods/P_AA_pp_Energy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/P_AA.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/MinimizationData.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



class P_AA_pp_Energy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	/// ctor
	P_AA_pp_Energy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////


	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::EnergyMap & emap
	) const override;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const override { return false; }

	/// @brief The P_AA_pp energy function describes derivatives wrt phi and psi.
	bool
	defines_dof_derivatives( pose::Pose const & p ) const override;

	/// @brief Return a vector of the atoms that determine the phi/psi DOFs
	utility::vector1< id::PartialAtomID >
	atoms_with_dof_derivatives( conformation::Residue const & res, pose::Pose const & pose ) const override;

	/// @brief Evaluate the P_AA_pp DOF derivative for a particular residue.
	Real
	eval_residue_dof_derivative(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights
	) const override;


	/// @brief APL Deprecated 6.29.2010
	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights
	) const;

	/// @brief P_AA_pp_Energy is context independent; indicates that no
	/// context graphs are required
	void indicate_required_context_graphs( utility::vector1< bool > & ) const override;


	// data
private:
	core::scoring::P_AA const & p_aa_;
	core::Size version() const override;

};

} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
