// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/FA_ElecEnergyCD.hh
/// @brief  Electrostatic energy with a distance-dependant dielectric
/// @author Hahnbeom Park

#ifndef INCLUDED_core_energy_methods_FA_ElecEnergyCD_hh
#define INCLUDED_core_energy_methods_FA_ElecEnergyCD_hh

// Unit Headers
#include <core/energy_methods/FA_GrpElecEnergy.fwd.hh>
#include <core/scoring/elec/GroupElec.hh>

// Package headers

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers


#include <basic/datacache/CacheableData.hh> // AUTO IWYU For CacheableData

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace energy_methods {

class FAElecContextData : public basic::datacache::CacheableData {

public:
	FAElecContextData();
	~FAElecContextData() override;

	void initialize( Size const nres );

	basic::datacache::CacheableDataOP clone() const override {
		return utility::pointer::make_shared< FAElecContextData >( *this );
	}

	Real &n( core::Size i ){ return n_[i]; };
	Real get_n( core::Size i ) const { return n_[i]; };
	Vector &dw_dr( core::Size i ){ return dw_dr_[i]; };
	Vector get_dw_dr( core::Size i ) const { return dw_dr_[i]; };
	utility::vector1< Size > &boundary_neighs( core::Size i ) { return boundary_neighs_[i]; };
	utility::vector1< Size > get_boundary_neighs( core::Size i ) const { return boundary_neighs_[i]; };

private:
	utility::vector1< Real > n_;
	utility::vector1< Vector > dw_dr_;
	utility::vector1< utility::vector1< Size > > boundary_neighs_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


class FA_GrpElecEnergy : public core::scoring::methods::ContextDependentTwoBodyEnergy  {
public:
	typedef core::scoring::methods::ContextDependentTwoBodyEnergy  parent;
public:


	FA_GrpElecEnergy( core::scoring::methods::EnergyMethodOptions const & options );


	FA_GrpElecEnergy( FA_GrpElecEnergy const & src );

	/// @brief Initilize constants.
	void
	initialize();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/// stashes nblist if use_nblist is true
	void
	setup_for_minimizing(
		pose::Pose & pose,
		core::scoring::ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map
	) const override;

	void
	setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const &scfxn ) const override;

	void
	setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const &scfxn ) const override;


	void
	setup_for_packing(
		pose::Pose & pose,
		utility::vector1< bool > const &,
		utility::vector1< bool > const &
	) const override;

	// Creates a rotamer trie for the input set of rotamers and stores the trie
	// in the rotamer set.
	void
	prepare_rotamers_for_packing(
		pose::Pose const &,
		conformation::RotamerSetBase & ) const override;

	// Updates the cached rotamer trie for a residue if it has changed during the course of
	// a repacking
	void
	update_residue_for_packing( pose::Pose &, Size ) const override;

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

	/// @brief Returns true if we're using neighborlist-autoupdate
	bool
	minimize_in_whole_structure_context( pose::Pose const & pose ) const override;

	bool
	defines_score_for_residue_pair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		bool res_moving_wrt_eachother
	) const override;

	bool
	use_extended_residue_pair_energy_interface() const override;

	void
	residue_pair_energy_ext(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResPairMinimizationData const & ,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::scoring::EnergyMap & emap
	) const override;

	void
	setup_for_minimizing_for_residue_pair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		kinematics::MinimizerMapBase const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData &
	) const override;

	/// @brief Evaluate the atom derivative f1/f2 vectors for all atoms on rsd1
	/// in response to the atoms on rsd2, and all the atoms on rsd2 as they
	/// in response to the atoms on rsd1.  This method is used with the
	/// MinimizationGraph and when nblist_autoupdate is not in use.
	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResSingleMinimizationData const &,
		core::scoring::ResPairMinimizationData const &,
		pose::Pose const & pose, // provides context
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs,
		utility::vector1< core::scoring::DerivVectorPair > & r2_atom_derivs
	) const override;

	/// @brief Evaluate the derivative vectors for a particular atom in a given
	/// (asymmetric) pose when nblist_autoupdate is being used.  nblist_autoupdate
	/// cannot be used with symmetric poses, in rtmin, or in minpack.
	void
	eval_atom_derivative(
		id::AtomID const &,
		pose::Pose const &,
		kinematics::DomainMap const &,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const &,
		Vector &,
		Vector &
	) const override;

	void
	finalize_total_energy(
		pose::Pose & ,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap &
	) const override;


	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & emap
	) const override;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		core::scoring::ResSingleMinimizationData const &,
		pose::Pose const & pose,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::scoring::DerivVectorPair > & atom_derivs
	) const override;

	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const override;


	//@brief overrides default rotamer/background energy calculation and uses
	// the trie-vs-trie algorithm instead
	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set,
		conformation::Residue const & residue,
		pose::Pose const & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap const & weights,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const override;


	bool
	defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const override { return false; }

	/// @brief Interface function for class core::scoring::NeighborList.
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_intrares_countpair(
		conformation::Residue const &,
		pose::Pose const &,
		core::scoring::ScoreFunction const &
	) const;

	/// @brief Interface function for class core::scoring::NeighborList.
	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size const,
		Size const,
		pose::Pose const &,
		core::scoring::ScoreFunction const &
	) const;

	core::scoring::etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;


	Distance
	atomic_interaction_cutoff() const override;

	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const override;

	/// Private methods
private:

	void
	set_nres_mono(
		core::pose::Pose const & pose
	) const;

	void
	precalc_context( pose::Pose & pose,
		FAElecContextDataOP data
	) const;


	Real
	eval_n( Real const cendist,
		Real &dn_drij,
		bool const eval_deriv ) const;

	void
	eval_context_derivatives(
		conformation::Residue const & rsd1,
		FAElecContextDataCOP data,
		core::scoring::EnergyMap const &,
		utility::vector1< core::scoring::DerivVectorPair > & r1_atom_derivs
	) const;

	core::Real
	burial_weight( core::Real const nb ) const;

	core::Real
	burial_deriv( core::Real const nb ) const;

	bool
	monomer_test(
		Size irsd,
		Size jrsd
	) const;

protected:

	inline
	core::scoring::etable::coulomb::Coulomb const &
	coulomb() const {return coulomb_; }

	inline
	core::scoring::elec::GroupElec const &
	groupelec() const {return groupelec_; }

private:

	core::scoring::etable::coulomb::Coulomb coulomb_;
	core::scoring::elec::GroupElec groupelec_;

	bool exclude_protein_protein_;
	bool exclude_RNA_RNA_;
	bool exclude_monomer_;
	bool exclude_DNA_DNA_;
	//Real intrares_scale_;
	bool context_dependent_;
	Real context_minstrength_;
	Real ncb_maxburial_, ncb_minburial_;


	mutable Size nres_monomer_;

	// temporary for derivative stuffs; turn off for now
	//mutable utility::vector1< Real > Eres_;

	core::Size version() const override;

};

} // namespace energy_methods
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_energy_methods_FA_GrpElecEnergy )
#endif // SERIALIZATION


#endif
