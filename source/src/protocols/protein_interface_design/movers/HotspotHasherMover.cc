// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/protein_interface_design/movers/HotspotHasherMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/HotspotHasherMover.hh>
#include <protocols/protein_interface_design/movers/HotspotHasherMoverCreator.hh>

// Package headers

// Project headers
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility>
#include <utility/file/file_sys_util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <utility/vector0.hh>
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

static basic::Tracer TR( "protocols.protein_interface_design.movers.HotspotHasherMover" );




HotspotHasherMover::HotspotHasherMover() : protocols::moves::Mover( HotspotHasherMover::mover_name() ) { }
HotspotHasherMover::HotspotHasherMover(
	std::vector<std::string> const & resnames,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size const n_stubs,
	std::string const & target_resnum,
	protocols::filters::FilterOP hotspot_filter,
	core::Real const target_distance,
	std::string const & hashin_fname,
	std::string const & hashout_fname
) :
	protocols::moves::Mover( HotspotHasherMover::mover_name() ),
	scorefxn_(std::move(scorefxn)),
	resnames_(resnames),
	n_stubs_(n_stubs),
	target_resnum_(target_resnum),
	target_distance_(target_distance),
	hashin_fname_(hashin_fname),
	hashout_fname_(hashout_fname)
{
	hotspot_filter_ = hotspot_filter;
}

HotspotHasherMover::~HotspotHasherMover() = default;

protocols::moves::MoverOP
HotspotHasherMover::clone() const {
	return( utility::pointer::make_shared< HotspotHasherMover >( *this ) );
}

void HotspotHasherMover::apply( core::pose::Pose & pose ) {

	// finding a request for ALL adds the main 18 aa's to
	if ( std::find( resnames_.begin(), resnames_.end(), "ALL" ) != resnames_.end() ) {
		resnames_.emplace_back("ALA" );
		resnames_.emplace_back("ARG" );
		resnames_.emplace_back("ASN" );
		resnames_.emplace_back("ASP" );
		resnames_.emplace_back("GLU" );
		resnames_.emplace_back("GLN" );
		resnames_.emplace_back("HIS" );
		resnames_.emplace_back("ILE" );
		resnames_.emplace_back("LEU" );
		resnames_.emplace_back("LYS" );
		resnames_.emplace_back("MET" );
		resnames_.emplace_back("PHE" );
		resnames_.emplace_back("PRO" );
		resnames_.emplace_back("SER" );
		resnames_.emplace_back("THR" );
		resnames_.emplace_back("TRP" );
		resnames_.emplace_back("TYR" );
		resnames_.emplace_back("VAL" );
	}

	protocols::hotspot_hashing::HotspotStubSet stubset;

	// read existing hashes
	if ( utility::file::file_exists( hashin_fname_ ) ) {
		stubset.read_data( hashin_fname_ );
		TR << "Found hash file " << hashin_fname_ << std::endl;
	}
	if ( utility::file::file_exists( hashout_fname_ ) ) { // useful for interrupted runs
		stubset.read_data( hashout_fname_ );
		TR << "Found hash file " << hashout_fname_ << std::endl;
	}

	core::Size target_resnum( 0 );
	if ( ! target_resnum_.empty() ) {
		target_resnum = core::pose::parse_resnum( target_resnum_, pose );
	}

	// for each residue requested
	for ( std::vector< std::string >::const_iterator it=resnames_.begin() ; it!=resnames_.end(); ++it ) {
		std::string resname = *it;

		TR << "Hash contains " << stubset.size(resname) << " " << resname << " stubs." << std::endl;

		// check to see if we've already finished our hash
		core::Size stubs_left = n_stubs_;
		stubs_left -= stubset.size( resname );
		if ( stubs_left <= 0 ) {
			// perform a scorecut
			if ( basic::options::option[ basic::options::OptionKeys::out::scorecut ].user() ) {
				Real score_cutoff = basic::options::option[ basic::options::OptionKeys::out::scorecut ]();
				std::stringstream i;
				i.str("");
				i << score_cutoff;
				stubset.clear();
				stubset.read_data( hashout_fname_ );
				protocols::hotspot_hashing::HotspotStubSetOP cut_stubs = stubset.subset( score_cutoff );
				std::string newfname = i.str() + "cut_" + hashout_fname_;
				cut_stubs->write_all( newfname );
			}
			return;
		}

		// do hashing in 10-stub cycles to minimize file i/o
		core::Size n_per(10);

		core::Size n_cycles = n_stubs_ / n_per;
		// make sure we do at least one cycle
		if ( n_cycles <= 0 ) n_cycles = 1;
		// PERFORM HASHING
		for ( core::Size i = 1; i <= n_cycles; ++i ) {
			stubset.clear();
			stubset.score_threshold( score_threshold_ );
			TR << "Finding " << n_per*i << "/" << n_stubs_ << " " << resname << " stubs" ;
			if ( ! target_resnum_.empty() ) {
				TR << " " << target_distance_ << "A from " << target_resnum_ << std::endl;
				stubset.fill( pose, scorefxn_, target_resnum, target_distance_, resname, n_per );
			} else {
				TR << "." << std::endl;
				stubset.fill( pose, scorefxn_, resname, n_per );
			}
			stubset.write_all( hashout_fname_ );
		}
	} // for each residue

	// perform a scorecut
	if ( basic::options::option[ basic::options::OptionKeys::out::scorecut ].user() ) {
		Real score_cutoff = basic::options::option[ basic::options::OptionKeys::out::scorecut ]();
		std::stringstream i;
		i.str("");
		i << score_cutoff;
		stubset.clear();
		stubset.read_data( hashout_fname_ );
		protocols::hotspot_hashing::HotspotStubSetOP cut_stubs = stubset.subset( score_cutoff );
		std::string newfname = i.str() + "cut_" + hashout_fname_;
		cut_stubs->write_all( newfname );
	}
} // HotspotHasherMover::apply



void
HotspotHasherMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data )
{

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();

	n_stubs_ = tag->getOption<core::Size>( "nstubs", 1000 );

	// target residue
	target_resnum_ = core::pose::get_resnum_string( tag, "target_residue_", "" );

	target_distance_ = tag->getOption<core::Real>( "target_distance", 15.0 );

	score_threshold_ = tag->getOption<core::Real>( "threshold", -1.0 );

	// hash in/out
	hashin_fname_ = tag->getOption<std::string>( "in", "");
	hashout_fname_ = tag->getOption<std::string>( "out", "hash.stubs");

	// filter
	std::string const hotspot_filter_name( tag->getOption<std::string>( "hotspot_filter", "true_filter" ) );
	protocols::filters::FilterOP hotspot_filter = protocols::rosetta_scripts::parse_filter_or_null( hotspot_filter_name, data );
	if ( hotspot_filter ) {
		hotspot_filter_ = hotspot_filter->clone();
	} else {
		if ( hotspot_filter_name != "true_filter" ) {
			TR<<"***WARNING WARNING! Filter defined for HotspotHasher not found in filter_list!!!!***"<<std::endl;
			utility_exit_with_message("Filter " + hotspot_filter_name + " not found.");
		} else {
			hotspot_filter_ = utility::pointer::make_shared< protocols::filters::TrueFilter >();
		}
	}

	// residues
	utility::vector0< TagCOP > const & hasher_tags( tag->getTags() );
	for ( auto hash_tag_ptr : hasher_tags ) {
		std::string tag_name = hash_tag_ptr->getName();
		if ( tag_name == "residue" ) {
			std::string resname( hash_tag_ptr->getOption< std::string >( "type", "" ) );
			resnames_.push_back( resname );
		}
	}
	runtime_assert( resnames_.size() > 0 );

	TR<<"hashing mover finding residues: ";
	for ( std::vector< std::string >::const_iterator it=resnames_.begin() ; it!=resnames_.end(); ++it ) TR<<*it<<" ";
	if ( ! target_resnum_.empty() ) TR << target_distance_ << "A away from residue " << target_resnum_ << std::endl;
	TR<<std::endl;
} // HotspotHasherMover::parse_my_tag

std::string HotspotHasherMover::get_name() const {
	return mover_name();
}

std::string HotspotHasherMover::mover_name() {
	return "HotspotHasher";
}

std::string subtag_for_hotspothasher( std::string const & tag ) {
	return "stfhsh_" + tag;
}

void HotspotHasherMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "nstubs", xsct_non_negative_integer, "Number of hotspot stubs to generate", "1000 ")
		+ XMLSchemaAttribute( "target_residue_pdb_num", xsct_refpose_enabled_residue_number, "Residue number from which to generate hot spot stubs; provide this or target_residue_res_num" )
		+ XMLSchemaAttribute( "target_residue_res_num", xsct_refpose_enabled_residue_number, "Residue number from which to generate hot spot stubs; provide this or target_residue_pdb_num" )
		+ XMLSchemaAttribute::attribute_w_default( "target_distance", xsct_real, "Interface distance definition", "15.0" )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Threshold below which stubs are favorable enough in energy to dump", "-1.0" )
		+ XMLSchemaAttribute( "in", xs_string, "Input for a checkpoint of the hotspot hashing effort" )
		+ XMLSchemaAttribute::attribute_w_default( "out", xs_string, "Output where hashed stubs should be stored", "hash.stubs" )
		+ XMLSchemaAttribute::attribute_w_default( "hotspot_filter", xs_string, "Filter to apply to all hotspots before storage", "true_filter" );

	AttributeList subtag_attributes;
	/*** Parse Sub-Tags ***/
	subtag_attributes + XMLSchemaAttribute::required_attribute( "type", xs_string, "Residue name for consideration" );

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "residue", subtag_attributes, "Tags describing individual residues hotspots of which are to be hashed"/*, 1 minoccurs*/ )
		.complex_type_naming_func( & subtag_for_hotspothasher );

	//protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_mover )
		.element_name( mover_name() )
		.description( "XRW TODO" )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( ssl, 1, xsminmax_unbounded )
		.write_complex_type_to_schema( xsd );
}

std::string HotspotHasherMoverCreator::keyname() const {
	return HotspotHasherMover::mover_name();
}

protocols::moves::MoverOP
HotspotHasherMoverCreator::create_mover() const {
	return utility::pointer::make_shared< HotspotHasherMover >();
}

void HotspotHasherMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HotspotHasherMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
