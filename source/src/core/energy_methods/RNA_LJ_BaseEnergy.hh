// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/RNA_LJ_BaseEnergy.hh
/// @brief  van der waals over aromatics
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_LJ_BaseEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_LJ_BaseEnergy_hh

// Unit Headers
#include <core/energy_methods/RNA_LJ_BaseEnergy.fwd.hh>

// Package headers
#include <core/conformation/Atom.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>


namespace core {
namespace energy_methods {


class RNA_LJ_BaseEnergy : public core::scoring::methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentTwoBodyEnergy  parent;

public:

	RNA_LJ_BaseEnergy( core::scoring::etable::Etable const & etable_in );

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	RNA_LJ_BaseEnergy( RNA_LJ_BaseEnergy const & src );

	void
	setup_for_derivatives(
		pose::Pose & pose,
		core::scoring::ScoreFunction const & scfxn
	) const override;

	/// called during gradient-based minimization inside dfunc
	/**
	F1 and F2 are not zeroed -- contributions from this atom are
	just summed in
	**/
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;


	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return false; }

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override {}

	Distance
	atomic_interaction_cutoff() const override;


	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

	Real
	eval_atom_energy(
		id::AtomID const & id,
		pose::Pose const & pose
	) const;

private:

	void
	eval_lj(
		conformation::Atom const & atom1,
		conformation::Atom const & atom2,
		Real const & d2,
		Real & fa_atr_score,
		Real & fa_rep_score,
		Real & deriv_atr,
		Real & deriv_rep
	) const;

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &// totals
	) const override;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:
	core::scoring::etable::Etable const & etable_; // shouldn't this be a pointer? Reference count information is (dangerously) lost when
	//a reference is taken, instead of a smart pointer.  There's the potential for a dangling reference with this.


	/// these guys are taken from the etable
	ObjexxFCL::FArray3D < Real > const & ljatr_;
	ObjexxFCL::FArray3D < Real > const & ljrep_;

	ObjexxFCL::FArray3D < Real > const & dljatr_;
	ObjexxFCL::FArray3D < Real > const & dljrep_;

	Real const safe_max_dis2_;
	Real const get_bins_per_A2_;

	bool const verbose_;
	core::Size version() const override;

};

} //scoring
} //core

#endif // INCLUDED_core_energy_methods_RNA_LJ_BaseEnergy_HH
