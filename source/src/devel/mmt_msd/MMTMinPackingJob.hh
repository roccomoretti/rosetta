// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/mmt_msd/MMTMinPackingJob.hh
/// @brief  declaration for class MMTMinPackingJob
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_devel_mmt_msd_MMTMinPackingJob_HH
#define INCLUDED_devel_mmt_msd_MMTMinPackingJob_HH

// Unit headers
#include <devel/mmt_msd/MMTMinPackingJob.fwd.hh>
#include <devel/mmt_msd/MMTPackingJob.hh>

// Core headers
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/pack/scmin/AtomTreeCollection.fwd.hh>
#include <core/pack/scmin/SCMinMinimizerMap.fwd.hh>
#include <core/pack/scmin/SidechainStateAssignment.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>

#include <core/scoring/annealing/RotamerSets.fwd.hh> // AUTO IWYU For RotamerSetsOP


namespace devel {
namespace mmt_msd {

class MMTMinPackingJob : public MMTPackingJob
{
public:
	typedef core::pack::scmin::SidechainStateAssignment SidechainStateAssignment;
	typedef core::pack::scmin::SidechainStateAssignmentOP SidechainStateAssignmentOP;

public:
	MMTMinPackingJob();
	~MMTMinPackingJob() override;

	void cartesian( bool setting );
	void nonideal( bool setting );

	void setup() override;
	void optimize() override;
	void update_pose( core::pose::Pose & pose ) override;

	bool best_assignment_exists() const;

	SidechainStateAssignment const &
	get_best_assignment() const;

	core::Real final_energy() const override;

protected:
	void clean_up() override;

private:

	bool cartesian_;
	bool nonideal_;

	core::pack::rotamer_set::RotamerSetsOP  rotsets_;
	core::pack::scmin::SCMinMinimizerMapOP  scminmap_;
	core::scoring::MinimizationGraphOP      mingraph_;
	core::pack::scmin::AtomTreeCollectionOP atc_;
	core::optimization::MinimizerOptionsOP  min_options_;
	core::pack::scmin::SidechainStateAssignmentOP best_assignment_;

};

}
}

#endif
