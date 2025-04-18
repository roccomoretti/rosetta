// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/TopologyBrokerMover.cc
/// @brief
/// @author Lei Shi (shilei@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/TopologyBrokerMover.hh>
#include <protocols/protein_interface_design/movers/TopologyBrokerMoverCreator.hh>

// Package headers

#include <core/pose/Pose.hh>
#include <core/pose/init_id_map.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>


//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/calc_taskop_movers/DesignRepackMover.hh>

//#include <protocols/topology_broker/TopologyBroker.hh>
//#include <protocols/topology_broker/util.hh>
#include <protocols/abinitio/BrokerMain.hh>
#include <protocols/abinitio/AbrelaxMover.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>

#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <protocols/protein_interface_design/movers/SaveAndRetrieveSidechains.hh>
#include <core/pose/PDBInfo.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.TopologyBrokerMover" );




TopologyBrokerMover::TopologyBrokerMover() : calc_taskop_movers::DesignRepackMover( TopologyBrokerMover::mover_name() ) {}
TopologyBrokerMover::~TopologyBrokerMover() = default;

void
TopologyBrokerMover::apply( pose::Pose & pose )
{
	using namespace core::scoring;
	using namespace core::pack::task::operation;

	if ( pose.conformation().num_chains()!=2 ) {
		utility_exit_with_message( "In TopologyBrokerMover, pose must contain contain exactly two chains: A is target, B is design to be folded" );
	}

	if ( align_ ) {
		if ( start_ < pose.conformation().chain_begin( 2 ) || end_ < pose.conformation().chain_begin( 2 )
				|| start_ > pose.conformation().chain_end( 2 ) || end_ > pose.conformation().chain_end( 2 )
				|| start_ >= end_ ) {
			TR.Error << "rigid_start " << start_ << " and rigid_end " << end_ << " is not in range " << pose.conformation().chain_begin( 2 ) <<"-" << pose.conformation().chain_end( 2 ) << std::endl;
			utility_exit_with_message( "To realign in TopologyBrokerMover, please provide start and end values in range." );
		}
	}

	//separate the two chains
	core::pose::Pose pose1 = pose;
	core::pose::Pose pose2 = pose;
	core::Size chain1primet(pose.conformation().chain_begin( 1 ));
	core::Size chain1end(pose.conformation().chain_end( 1 ) );
	core::Size chain2primet(pose.conformation().chain_begin( 2 ));
	core::Size chain2end(pose.conformation().chain_end( 2 ) );
	pose1.conformation().delete_residue_range_slow( chain2primet, chain2end);
	pose2.conformation().delete_residue_range_slow( chain1primet, chain1end);
	core::pose::Pose ref_pose = pose2;

	protocols::protein_interface_design::movers::SaveAndRetrieveSidechainsOP srsc( new protocols::protein_interface_design::movers::SaveAndRetrieveSidechains(ref_pose, true, false, 0) );

	protocols::abinitio::AbrelaxApplication::register_options();
	protocols::abinitio::register_options_broker();
	protocols::abinitio::AbrelaxMoverOP abrelax( new protocols::abinitio::AbrelaxMover() );
	//abrelax->set_b_return_unrelaxed_fullatom(true);
	abrelax->apply( pose2 );
	//switch to fa output, relax is set in commond line
	//if ( !pose.is_fullatom() ) {
	//  core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
	//}
	srsc->apply(pose2);

	// superimpose
	core::id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map_AtomID ( atom_map, pose2, id::AtomID::BOGUS_ATOM_ID() );
	for ( core::Size i=start_; i<=end_; ++i ) {
		core::id::AtomID const id1( pose2.residue(i-chain1end).atom_index("CA"), i-chain1end );
		core::id::AtomID const id2( ref_pose.residue(i-chain1end).atom_index("CA"), i-chain1end);
		atom_map[ id1 ] = id2;
	}
	core::scoring::superimpose_pose( pose2, ref_pose, atom_map );
	Real rms=core::scoring::CA_rmsd(pose2,ref_pose,start_-chain1end,end_-chain1end);
	if ( rms > 0.001 ) {
		utility_exit_with_message( "Something is wrong, the rigid part is not rigid, check your rigid_start, rigid_end and setup.tpb for broker !" );
	}
	//runtime_assert( rms < 0.001 );

	//coordinate constrained relax

	// append target
	pose1.append_pose_by_jump(pose2,chain1end);
	pose=pose1;
	//copy pdbinfo
	core::pose::PDBInfoOP pdbinfo( new core::pose::PDBInfo( pose) );
	pose.pdb_info(pdbinfo);

}


void
TopologyBrokerMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & )
{
	TR<<"Setup TopologyBrokerMover Mover " << std::endl;

	align_ = tag->getOption< bool >( "realign", false );
	if ( align_ ) {
		if ( tag->hasOption("rigid_start") && tag->hasOption("rigid_end") ) {
			start_ = tag->getOption< core::Size >( "rigid_start", 1 );
			end_ = tag->getOption< core::Size >( "rigid_end", 1 );
		} else {
			utility_exit_with_message( "If you want to realign, please provide motif rigid_start and rigid_end in pose numbering!" );
		}
	} else  {
		utility_exit_with_message( "Expect to use realign option for futher design purpose" );
	}
}

protocols::moves::MoverOP
TopologyBrokerMover::clone() const {
	return( utility::pointer::make_shared< TopologyBrokerMover >( *this ));
}

std::string TopologyBrokerMover::get_name() const {
	return mover_name();
}

std::string TopologyBrokerMover::mover_name() {
	return "TopologyBrokerMover";
}

void TopologyBrokerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "realign", xsct_rosetta_bool, "Realign must be true, but its default is false. We are all children in the eyes of the Broker.", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rigid_start", xsct_positive_integer, "Rigid segment startend; must be on the second chain.", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "rigid_end", xsct_positive_integer, "Rigid segment end; must be on the second chain", "1" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string TopologyBrokerMoverCreator::keyname() const {
	return TopologyBrokerMover::mover_name();
}

protocols::moves::MoverOP
TopologyBrokerMoverCreator::create_mover() const {
	return utility::pointer::make_shared< TopologyBrokerMover >();
}

void TopologyBrokerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TopologyBrokerMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
