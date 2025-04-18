// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/rna/data/RNA_ChemicalMappingEnergy.hh
/// @brief  Statistically derived chemical mapping score
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_ChemicalMappingEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_ChemicalMappingEnergy_hh

// Unit Headers
#include <core/energy_methods/RNA_ChemicalMappingEnergy.fwd.hh>
#include <core/scoring/rna/data/RNA_DMS_Potential.fwd.hh>
#include <core/scoring/rna/data/RNA_DMS_LowResolutionPotential.fwd.hh>

// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

namespace core {
namespace energy_methods {


class RNA_ChemicalMappingEnergy : public core::scoring::methods::WholeStructureEnergy {
public:

	typedef core::scoring::methods::WholeStructureEnergy parent;


	RNA_ChemicalMappingEnergy();

	~RNA_ChemicalMappingEnergy() override;


	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	Real
	calculate_energy( pose::Pose & pose,
		bool const use_low_res = false,
		bool const rna_base_pair_computed = false ) const;

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &// totals
	) const override;

	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {}

	Size version() const override { return 0; }

private:

	core::scoring::rna::data::RNA_DMS_Potential & DMS_potential_;
	core::scoring::rna::data::RNA_DMS_LowResolutionPotential & DMS_low_resolution_potential_;

};


} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
