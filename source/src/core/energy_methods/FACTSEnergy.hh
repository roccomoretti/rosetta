// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// @file   core/scoring/facts/FACTSEnergy.hh
// @brief
// @author Hahnbeom Park

#ifndef INCLUDED_devel_facts_FACTSEnergy_HH
#define INCLUDED_devel_facts_FACTSEnergy_HH

// Unit Headers
#include <core/scoring/facts/FACTSPotential.fwd.hh>
#include <core/energy_methods/FACTSEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/EnergyMap.fwd.hh>
//#include <core/pack/task/PackerTask.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
//#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace core {
namespace energy_methods {


class FACTSEnergy : public core::scoring::methods::ContextDependentTwoBodyEnergy  {

public:
	typedef core::scoring::methods::ContextDependentTwoBodyEnergy parent;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::facts::FACTSPotential FACTSPotential;
	typedef core::scoring::methods::EnergyMethod EnergyMethod;
	typedef core::scoring::methods::EnergyMethodOP EnergyMethodOP;
	typedef core::scoring::methods::EnergyMethodOptions EnergyMethodOptions;

public:
	FACTSEnergy( core::scoring::methods::EnergyMethodOptions const & options );

	FACTSEnergy( FACTSEnergy const & src );

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override {
		return utility::pointer::make_shared< FACTSEnergy >(*this);
	}

	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const & residues_repacking,
		utility::vector1< bool > const &
	) const override;

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & ) const override;

	void
	prepare_rotamers_for_packing(
		pose::Pose const & pose,
		conformation::RotamerSetBase & set
	) const override;

	void
	update_residue_for_packing(
		pose::Pose &,
		Size resid ) const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;


	void
	evaluate_rotamer_intrares_energies(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		utility::vector1< core::PackerEnergy > & energies
	) const override;


	void
	evaluate_rotamer_intrares_energy_maps(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		utility::vector1< core::scoring::EnergyMap > & emaps
	) const override;


	/// @brief Batch computation of rotamer pair energies.  Need not be overriden in
	/// derived class -- by default, iterates over all pairs of rotamers,
	/// and calls derived class's residue_pair_energy method.  Since short range rotamer pairs
	/// may not need calculation, the default method looks at blocks of residue type pairs
	/// and only calls the residue_pair_energy method if the rotamer pairs are within range
	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const override;


	/// @brief Batch computation of rotamer/background energies.  Need not be overriden
	/// in derived class -- by default, iterates over all rotamers in the set, and calls
	/// derived class's residue_pair_energy method for each one against the background rotamer
	/// Since short range rotamer pairs may not need calculation, the default method
	/// looks at blocks of residue type pairs and only calls the residue_pair_energy method
	/// if the rotamer pairs are within range
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const override;


	/// @brief Batch computation of rotamer/background energies.  Need not be overriden
	/// in derived class -- by default, iterates over all rotamers in the set, and calls
	/// derived class's residue_pair_energy method for each one against the background rotamer
	/// Since short range rotamer pairs may not need calculation, the default method
	/// looks at blocks of residue type pairs and only calls the residue_pair_energy method
	/// if the rotamer pairs are within range

	void
	evaluate_rotamer_background_energy_maps(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::EnergyMap > & emaps
	) const override;


	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {}

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return true; }

	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override;

	/// this is our own special function
	Real
	packing_interaction_cutoff() const { return 5.5; }

private:
	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	// const-ref to scoring database
	FACTSPotential const & potential_;

	bool const exclude_DNA_DNA_;
	Real max_dis_;

	core::Size version() const override;

};


}
}

#endif
