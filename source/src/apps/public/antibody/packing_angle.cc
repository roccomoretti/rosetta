// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/nmarze/tilt_angle.cc
/// @brief pilot app for extracting VL/VH packing angle from pdb/list of pdbs
/// @author Nick Marze (nickmarze@gmail.com)


#include <protocols/jd2/internal_util.hh>
#include <utility/excn/Exceptions.hh>
//#include <utility/exit.hh>

#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <protocols/antibody/metrics.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

static basic::Tracer TR( "apps.pilot.nmarze.packing_angle" );

using utility::vector1;

class PackingAngle : public protocols::moves::Mover {
public:
	// default constructor
	PackingAngle()= default;
	// destructor
	~PackingAngle() override= default;

	void apply( core::pose::Pose & pose_in ) override
	{
		TR << "Applying Packing Angle Calculator" << std::endl;

		using namespace core;
		using namespace protocols;
		using namespace core::pose;
		using namespace protocols::moves;
		using namespace protocols::antibody;

		protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

		antibody::AntibodyInfoCOP ab_info_ = antibody::AntibodyInfoCOP ( new AntibodyInfo( pose_in ) );
		PoseCOP new_pose = utility::pointer::make_shared< Pose >( pose_in );

		vector1< Real > orientation_coords_ = vl_vh_orientation_coords( *new_pose , *ab_info_ );

		job->add_string_real_pair( "VL_VH_distance", orientation_coords_[1] );
		job->add_string_real_pair( "VL_VH_opening_angle", orientation_coords_[2] );
		job->add_string_real_pair( "VL_VH_opposite_opening_angle", orientation_coords_[3] );
		job->add_string_real_pair( "VL_VH_packing_angle", orientation_coords_[4] );
		job->add_string_real_pair( "H1_length", ab_info_->get_CDR_length( h1 ) );
		job->add_string_real_pair( "H2_length", ab_info_->get_CDR_length( h2 ) );
		job->add_string_real_pair( "H3_length", ab_info_->get_CDR_length( h3 ) );
		job->add_string_real_pair( "L1_length", ab_info_->get_CDR_length( l1 ) );
		job->add_string_real_pair( "L2_length", ab_info_->get_CDR_length( l2 ) );
		job->add_string_real_pair( "L3_length", ab_info_->get_CDR_length( l3 ) );

		TR << "Finished applying Packing Angle Calculator" << std::endl;

		return;
	} // PackingAngle::apply()

	std::string get_name() const override { return "PackingAngle"; }


	protocols::moves::MoverOP
	fresh_instance() const override {
		return utility::pointer::make_shared< PackingAngle >();
	}


	bool
	reinitialize_for_each_job() const override { return false; }


	bool
	reinitialize_for_new_input() const override { return false; }

private:

};

using PackingAngleOP = utility::pointer::shared_ptr<PackingAngle>;


/////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		protocols::jd2::register_options();

		// initialize core
		devel::init(argc, argv);

		PackingAngleOP packing_angle = utility::pointer::make_shared< PackingAngle >();
		protocols::jd2::JobDistributor::get_instance()->go( packing_angle );

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
