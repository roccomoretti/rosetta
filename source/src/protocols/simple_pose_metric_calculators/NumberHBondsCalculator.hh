// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/NumberHBondsCalculator.hh
/// @brief calculates h-bonds in pose
/// @author Florian Richter

#ifndef INCLUDED_protocols_toolbox_simple_pose_metric_calculators_NumberHBondsCalculator_hh
#define INCLUDED_protocols_toolbox_simple_pose_metric_calculators_NumberHBondsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>
//#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>


// utility headers
#include <utility/vector1.hh>

// C++ headers
#include <set>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace simple_pose_metric_calculators {

/// @brief use the "score12_params" set if the -restore_pre_talaris_2013_behavior flag is on the command line
/// and otherwise use the new and improved sp2_elec_params parameter set
std::string
choose_hbond_parameter_set();

class NumberHBondsCalculator : public core::pose::metrics::EnergyDependentCalculator {
public:

	NumberHBondsCalculator();
	NumberHBondsCalculator( bool const generous, std::set< core::Size > special_region=std::set< core::Size >() );
	~NumberHBondsCalculator() override;

	//core::pose::metrics::PoseMetricCalculatorOP clone() const;

	core::pose::metrics::PoseMetricCalculatorOP clone() const override {
		return utility::pointer::make_shared< NumberHBondsCalculator >( *this );
	};

	static core::Real
	sum_Hbond_terms( core::scoring::EnergyMap const & emap );

	// getters and setters
	void set_generous( bool generous ){ use_generous_hbonds_ = generous; };
	void set_special_region( std::set< core::Size > selection ){ special_region_ = selection; };
	core::Size get_total_hbonds() { return all_Hbonds_; };
	std::set< core::Size > get_special_region(){ return special_region_; };
	//core::id::AtomID_Map< core::Size > get_atom_honds(){ return atom_Hbonds_; };
	utility::vector1< core::Size > get_residue_hbonds(){ return residue_Hbonds_; };
	//core::scoring::hbonds::HBondSetCOP get_hbond_set() { return hbond_set_; };
	void set_max_hb_energy( core::Real energy ) { max_hb_energy_ = energy; }
	void set_cumul_hbonds( bool cumul_hbonds ) { cumul_hbonds_ = cumul_hbonds; }

protected:

	void lookup( std::string const & key, basic::MetricValueBase * valptr ) const override;
	std::string print( std::string const & key ) const override;
	void recompute( core::pose::Pose const & this_pose ) override;

	void
	determine_res_to_recompute(
		core::pose::Pose const & pose,
		utility::vector1< bool > & res_to_recompute
	);

	void
	compute_Hbonds_for_residue(
		core::pose::Pose const & pose,
		core::Size i,
		utility::vector1< bool > const & res_to_recompute,
		core::scoring::hbonds::HBondSet & hb_set
	);

private:
	bool use_generous_hbonds_;
	std::set< core::Size > special_region_;
	core::Size all_Hbonds_;
	core::Size special_region_Hbonds_;
	core::id::AtomID_Map< core::Size > atom_Hbonds_;
	utility::vector1< core::Size > residue_Hbonds_;
	utility::vector1< core::Real > ref_residue_total_energies_; //holds the calculated energies to prevent unnecessary recalculation
	core::Real max_hb_energy_;
	bool cumul_hbonds_;
	//core::scoring::hbonds::HBondSetOP hbond_set_;
	//core::scoring::hbonds::HBondDatabaseCOP hb_database_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace simple_pose_metric_calculators
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_simple_pose_metric_calculators_NumberHBondsCalculator )
#endif // SERIALIZATION


#endif
