// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs/PseudocontactShiftEnergy.hh
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs_PseudocontactShiftEnergy_hh
#define INCLUDED_protocols_scoring_methods_pcs_PseudocontactShiftEnergy_hh

// Package headers
#include <protocols/scoring/methods/pcs/PseudocontactShiftData.fwd.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftTensor.fwd.hh>
#include <protocols/scoring/methods/pcs/PseudocontactShiftEnergy.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Utility headers
#include <utility/SingletonBase.hh>

// Numeric headers
#include <numeric/xyzVector.fwd.hh>

#include <utility/vector1.hh>

// Objexx headers

// C++ headers

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs {

class PCS_Energy : public core::scoring::methods::WholeStructureEnergy {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent;

public:

	PCS_Energy();

	~PCS_Energy() override;

	PCS_Energy &
	operator=(PCS_Energy const & other);

	PCS_Energy(PCS_Energy const & other);

	core::scoring::methods::EnergyMethodOP
	clone() const override;

	void
	indicate_required_context_graphs( utility::vector1< bool > & ) const override;

	void
	finalize_total_energy(
		core::pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;

	core::Real
	calculate_pcs_score( core::pose::Pose & pose, bool print_to_tracer ) const;

	core::Real
	calculate_scores_and_tensors_from_pose_and_PCS_data(
		utility::vector1<core::Real> & vec_best_score,
		utility::vector1<PCS_tensor> & vec_best_tensor,
		numeric::xyzVector< core::Real > & best_coo,
		core::pose::Pose const & pdb,
		PCS_data & pcs_d
	) const;

	core::Real
	minimize_tensors_from_PCS_data(
		utility::vector1<PCS_tensor> & vec_best_tensor,
		numeric::xyzVector< core::Real > & best_coo,
		PCS_data const & pcs_d
	) const;

	/*core::Real
	calculate_single_score_and_tensor_from_PCS_data_per_lanthanides(
	PCS_tensor & PCS_t,
	PCS_data_per_lanthanides & pcs_d_p_l
	) const; */


	PCS_data &
	PCS_data_from_pose(core::pose::Pose & pose) const;

	void
	dump_PCS_info(
		utility::vector1<PCS_tensor> const & vec_tensor,
		numeric::xyzVector< core::Real > const & best_coo,
		PCS_data const &pcs_d
	) const;

	void
	show_additional_info(std::ostream & out, core::pose::Pose & pose, bool verbose=false) const override;

	core::Size version() const override;
};


class PCS_Energy_parameters_manager : public utility::SingletonBase< PCS_Energy_parameters_manager > {
public:

	friend utility::SingletonBase< PCS_Energy_parameters_manager >;

	PCS_Energy_parameters_manager(PCS_Energy_parameters_manager const &) = delete;

	PCS_Energy_parameters_manager &
	operator=( PCS_Energy_parameters_manager const &) = delete;

private:

	PCS_Energy_parameters_manager();

	core::Real grid_edge_;
	core::Real grid_step_;
	core::Real grid_small_cutoff_;
	core::Real grid_large_cutoff_;
	core::Real grid_cone_angle_cutoff_;
	std::string grid_atom_name_1_;
	std::string grid_atom_name_2_;
	core::Size grid_residue_num_1_;
	core::Size grid_residue_num_2_;
	core::Real grid_k_vector_;
	bool minimize_best_tensor_;
	core::Real pcs_weight_;
	utility::vector1<std::string> vec_filename_;
	utility::vector1<core::Real> vec_individual_weight_;

	utility::vector1< bool > vec_exclude_residues_;
	bool vec_exclude_residues_exists_;
	bool vec_exclude_residues_changed_;

public:

	//rvernon -> partial PCS score machinery in development
	void
	set_vector_exclude_residues(utility::vector1< core::Size > const vec_exclude);

	void
	remove_vector_exclude_residues();

	bool
	has_exclude_residues_vector();

	bool
	has_exclude_residues_vector_changed();

	void
	exclude_residues_vector_is_current();

	utility::vector1< bool >
	get_vector_exclude_residues();
	//rvernon


	void
	set_vector_name_and_weight(utility::vector1<std::string> const vec_filename,
		utility::vector1<core::Real> const vec_individual_weight);

	void
	set_grid_param(core::Real const grid_edge,
		core::Real const grid_step,
		core::Real const grid_small_cutoff,
		core::Real const grid_large_cutoff,
		core::Real const grid_cone_angle_cutoff,
		std::string const grid_atom_name_1,
		std::string const grid_atom_name_2,
		core::SSize const grid_residue_num_1,
		core::SSize const grid_residue_num_2,
		core::Real const grid_k_vector,
		bool const minimize_best_tensor,
		core::Real const pcs_weight
	);
	void
	print_grid_param() const;

	core::Real
	get_grid_edge() const;

	core::Real
	get_grid_step() const;

	core::Real
	get_grid_small_cutoff() const;

	core::Real
	get_grid_large_cutoff() const;

	core::Real
	get_grid_cone_angle_cutoff() const;

	std::string
	get_grid_atom_name_1() const;

	std::string
	get_grid_atom_name_2() const;

	core::Size
	get_grid_residue_num_1() const;

	core::Size
	get_grid_residue_num_2() const;

	core::Real
	get_grid_k_vector() const;

	bool
	get_minimize_best_tensor() const;

	core::Real
	get_pcs_weight() const;

	utility::vector1<std::string> const &
	get_vector_filename() const;

	utility::vector1<core::Real> const &
	get_vector_weight() const;

};

} //pcs
} //methods
} //scoring
} //protocols

#endif // INCLUDED_protocols_scoring_methods_PseudocontactShiftEnergy_HH
