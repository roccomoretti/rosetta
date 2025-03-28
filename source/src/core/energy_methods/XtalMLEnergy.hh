// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/cryst/XtalMLEnergy.hh
/// @brief  ML target
/// @author Frank DiMaio


#ifndef INCLUDED_core_scoring_cryst_XtalMLEnergy_hh
#define INCLUDED_core_scoring_cryst_XtalMLEnergy_hh

// Package headers

#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <numeric/xyzVector.hh>

#include <string>

namespace core {
namespace energy_methods {

class XtalMLEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;

public:


	XtalMLEnergy();


	XtalMLEnergy( XtalMLEnergy const & src ) :
		parent(src),
		dml_dx(src.dml_dx),
		ml(src.ml)
	{}

	core::scoring::methods::EnergyMethodOP clone() const override;

	void finalize_total_energy( pose::Pose & pose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & totals ) const override;

	void setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & scorefxn ) const override;

	void setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & sf) const override;

	void
	setup_for_minimizing( pose::Pose & pose, core::scoring::ScoreFunction const & sf, kinematics::MinimizerMapBase const & ) const override;

	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const &, // domain_map,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const override;

	void indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

private:
	//////
	// PRIVATE DATA
	//////
	// precomputed derivatives
	mutable utility::vector1< utility::vector1< numeric::xyzVector< core::Real > > > dml_dx;

	// ml target function
	mutable core::Real ml;

	core::Size version() const override;
};

}
}

#endif
