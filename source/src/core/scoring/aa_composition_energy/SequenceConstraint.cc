// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/aa_composition_energy/SequenceConstraint.cc
/// @brief A base class for constraining sequences, analogous to a geometric constraint.
/// @details
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <core/scoring/aa_composition_energy/SequenceConstraint.hh>

#include <basic/Tracer.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION



namespace core {
namespace scoring {
namespace aa_composition_energy {

static basic::Tracer TR( "core.scoring.aa_composition_energy.SequenceConstraint" );

/// @brief Constructor
///
SequenceConstraint::SequenceConstraint( core::scoring::ScoreType const & t ):
	core::scoring::constraints::Constraint( t ),
	dummy_atomid_(0,0)
	//TODO -- initialize variables here.
{}

/// @brief Copy constructor
///
SequenceConstraint::SequenceConstraint( SequenceConstraint const &src ):
	core::scoring::constraints::Constraint( src.score_type() ),
	dummy_atomid_(0,0)
	//TODO -- copy variables here.
{
}

/// @brief Destructor
///
SequenceConstraint::~SequenceConstraint() = default;

} // aa_composition_energy
} // scoring
} // core

#ifdef    SERIALIZATION
/// @details Default constructor that's needed by Cereal in order to deserialize a constraint.
/// This should ONLY be used by derived classes in their default constructors and then only
/// for the sake of deserializing a constraint.
core::scoring::aa_composition_energy::SequenceConstraint::SequenceConstraint() :
	core::scoring::constraints::Constraint(),
	dummy_atomid_(0,0)
	//TODO -- initialize variables here.
{}

template< class Archive >
void
core::scoring::aa_composition_energy::SequenceConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::constraints::Constraint >( this ) );
	arc( CEREAL_NVP( dummy_atomid_ ) );
}

template< class Archive >
void
core::scoring::aa_composition_energy::SequenceConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::constraints::Constraint >( this ) );
	arc( dummy_atomid_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::aa_composition_energy::SequenceConstraint );
CEREAL_REGISTER_TYPE( core::scoring::aa_composition_energy::SequenceConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_aa_composition_energy_SequenceConstraint )
#endif // SERIALIZATION
