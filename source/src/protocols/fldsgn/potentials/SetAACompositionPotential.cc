// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/potentials/SetAACompositionEnergy.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// unit headers
#include <protocols/fldsgn/potentials/SetAACompositionPotential.hh>
#include <protocols/fldsgn/potentials/SetAACompositionPotentialCreator.hh>
#include <protocols/fldsgn/potentials/AACompositionEnergy.hh>

// C++ headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

// utility header
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

// boost
#include <boost/lexical_cast.hpp>

// parser
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {


static basic::Tracer TR( "protocols.fldsgn.potentials.SetAACompositionPotential" );


/// @brief


/// @brief


/// @brief


/// @brief default constructor
SetAACompositionPotential::SetAACompositionPotential() :
	Super( "SetAACompositionPotential" ),
	loaded_( false )
{}


/// @Brief copy constructor
SetAACompositionPotential::SetAACompositionPotential( SetAACompositionPotential const & /*rval*/ ) = default;


/// @brief default destructor
SetAACompositionPotential::~SetAACompositionPotential() = default;


/// @brief clone this object
SetAACompositionPotential::MoverOP
SetAACompositionPotential::clone() const
{
	return utility::pointer::make_shared< SetAACompositionPotential >( *this );
}


/// @brief create this type of object
SetAACompositionPotential::MoverOP
SetAACompositionPotential::fresh_instance() const
{
	return utility::pointer::make_shared< SetAACompositionPotential >();
}

/// @brief set parameters
bool
SetAACompositionPotential::set_parameters( String const & file )
{
	utility::io::izstream data( file );
	if ( !data ) {
		TR.Error << "can not open file " << file << std::endl;
		return false;
	}

	String line;
	while ( getline( data, line ) ) {
		utility::vector1< String > tokens ( utility::split( line ) );
		if ( tokens[1][0] == '#' ) continue;    // skip reading line that is commented out
		runtime_assert( tokens.size() == 3 );

		core::chemical::AA const aa = core::chemical::aa_from_name( tokens[1] );
		auto const lower_threshold = boost::lexical_cast<Real>( tokens[2] );
		auto const upper_threshold = boost::lexical_cast<Real>( tokens[3] );

		runtime_assert( lower_threshold <= upper_threshold );

		std::pair< Real, Real > thresholds( lower_threshold, upper_threshold );

		comp_constraint_aas_.insert( std::map< core::chemical::AA, std::pair< Real, Real >  >::value_type( aa, thresholds ) );

	}


	Real checkL( 0.0 ), checkU( 0.0 );
	auto it = comp_constraint_aas_.begin();
	while ( it != comp_constraint_aas_.end() ) {
		AA aa( it->first );
		checkL += comp_constraint_aas_[ aa ].first;
		checkU += comp_constraint_aas_[ aa ].second;
		++it;
	}
	runtime_assert( checkL <= 1.0 && checkU <= 1.0 );


	return true;
} // set_parameters


/// @brief apply defined moves to given Pose
void
SetAACompositionPotential::apply( Pose & )
{
	using core::scoring::ScoreType;

	if ( loaded_ ) return;

	protocols::fldsgn::potentials::AACompositionEnergy cce( comp_constraint_aas_ );

	std::map< ScoreType, Real > new_weights;
	new_weights.insert( std::map< ScoreType, Real >::value_type( core::scoring::aa_cmp, weight_ ) );

	sfx_->set_weight( core::scoring::aa_cmp, 0.0 );
	sfx_->add_extra_method( new_weights, cce );

	//sfx_->show( TR, pose );
}


/// @brief
/// @brief parse xml
void
SetAACompositionPotential::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data
)
{

	std::string const file( tag->getOption<String>( "file", "" ) );
	if ( file == "" ) {
		TR.Fatal << "No input of file ! " << std::endl;
		utility_exit_with_message( "No input of file ! " );
	}
	set_parameters( file );

	weight_ = ( tag->getOption<Real>( "weight", 1.0 ) );

	// set scorefxn
	String const sfxn ( tag->getOption<String>( "scorefxn", "" ) );
	sfx_ = data.get_ptr<ScoreFunction>( "scorefxns", sfxn );
	if ( sfxn == "" ) {
		TR.Fatal << "No input of sfxn ! " << std::endl;
		utility_exit_with_message( "No input of sfxn ! " );
	}

}

std::string SetAACompositionPotential::get_name() const {
	return mover_name();
}

std::string SetAACompositionPotential::mover_name() {
	return "SetAACompositionPotential";
}

void SetAACompositionPotential::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "file", xs_string, "File specifying desired amino acid composition" )
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "Weight for amino acid composition potential", "1.0" )
		+ XMLSchemaAttribute::required_attribute( "scorefxn", xs_string, "Score function to add this potential to" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Sets the target amino acid composition of a pose and adds the relevant term to the indicated score function", attlist );
}

std::string SetAACompositionPotentialCreator::keyname() const {
	return SetAACompositionPotential::mover_name();
}

protocols::moves::MoverOP
SetAACompositionPotentialCreator::create_mover() const {
	return utility::pointer::make_shared< SetAACompositionPotential >();
}

void SetAACompositionPotentialCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetAACompositionPotential::provide_xml_schema( xsd );
}



} // namespace potentials
} // namespace fldsgn
} // namespace protocols
