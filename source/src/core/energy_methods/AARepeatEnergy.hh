// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/aa_repeat_energy/AARepeatEnergy.hh
/// @brief Headers for an EnergyMethod that penalizes stretches of a repeating amino acid (e.g. poly-Q sequences).
/// @details This energy method is inherently not pairwise decomposible.  However, it is intended for very rapid calculation,
/// and has been designed to plug into Alex Ford's modifications to the packer that permit it to work with non-pairwise scoring
/// terms.
/// @author Vikram K. Mulligan (vmullig@uw.edu).



#ifndef INCLUDED_core_scoring_aa_repeat_energy_AARepeatEnergy_hh
#define INCLUDED_core_scoring_aa_repeat_energy_AARepeatEnergy_hh

// Unit headers
#include <core/energy_methods/AARepeatEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/annealing/ResidueArrayAnnealableEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {

/// @brief AARepeatEnergy, an energy function to penalize stretches of the same residue,
/// derived from base class for EnergyMethods, which are meaningful only on entire structures.
/// These EnergyMethods do all of their work in the "finalize_total_energy" section of score
/// function evaluation.
class AARepeatEnergy : public core::scoring::methods::WholeStructureEnergy, public core::scoring::annealing::ResidueArrayAnnealableEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent1;
	typedef core::scoring::annealing::ResidueArrayAnnealableEnergy parent2;
	typedef core::scoring::ScoreFunction ScoreFunction;
	typedef core::scoring::EnergyMap EnergyMap;

public:

	/// @brief Default constructor.
	///
	AARepeatEnergy();

	/// @brief Copy constructor.
	///
	AARepeatEnergy( AARepeatEnergy const &src );

	/// @brief Default destructor.
	///
	~AARepeatEnergy() override;

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	core::scoring::methods::EnergyMethodOP clone() const override;

	/// @brief AARepeatEnergy is context-independent and thus indicates that no context graphs need to be maintained by
	/// class Energies.
	void indicate_required_context_graphs( utility::vector1< bool > &context_graphs_required ) const override;

	/// @brief AARepeatEnergy is version 1.0 right now.
	///
	core::Size version() const override;

	/// @brief Actually calculate the total energy
	/// @details Called by the scoring machinery.
	void finalize_total_energy( core::pose::Pose & pose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & totals ) const override;

	/// @brief Calculate the total energy given a vector of const owning pointers to residues.
	/// @details Called by finalize_total_energy().
	core::Real calculate_energy(
		utility::vector1< core::conformation::ResidueCOP > const & resvect,
		utility::vector1< core::Size > const & rotamer_ids,
		core::Size const substitution_position = 0
	) const override;

private:

	/******************
	Private functions:
	******************/

	/// @brief Return a penalty for N residues in a row, from a lookup table.
	/// @details The last entry in the lookup table is the penalty to return for more residues
	/// than there are entries in the lookup table.  Returns 0 if nres is zero.
	inline core::Real penalties(core::Size const nres) const;

	/// @brief Read a penalty data file and load the penalties into the energy method.  Called once by the constructor.
	/// @details Comment lines are ignored in the penalty data file.  The file should have one line that's a whitespace-separated
	/// row of numbers.  The numbers represent the penalty for having a stretch of 1, 2, 3, 4, etc. of the same residue as a repeating
	/// stretch.  The function looks in the working directory and in database/scoring/score_functions/aa_repeat_energy/ for the penalty
	/// file.
	void load_penalty_table_from_file( std::string const &filename );

	/// @brief Parse a series of floats from a string and store them in the penalties vector.  Return true if successful and false otherwise.
	/// @details Only modifies the penalties vector if successful.
	bool parse_line( std::string const &line, utility::vector1< core::Real > &penalties ) const;

	/******************
	Private variables:
	******************/

	utility::vector1 < core::Real > penalties_;

};

} // energy_methods
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
