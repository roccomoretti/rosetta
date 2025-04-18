// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/id/JumpID.cc
/// @brief  Jump identifier
/// @author Phil Bradley

// Unit headers
#include <core/id/JumpID.hh>

#include <iostream>

namespace core {
namespace id {

std::ostream &
operator <<(
	std::ostream & os,
	JumpID const & a
)
{
	os << "JumpID " << a.rsd1_ << ' ' << a.rsd2_;
	return os;
}

} // namespace id
} // namespace core

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

/// @brief Automatically generated serialization method
template< class Archive >
void
core::id::JumpID::save( Archive & arc ) const {
	arc( CEREAL_NVP( rsd1_ ) ); // Size
	arc( CEREAL_NVP( rsd2_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::id::JumpID::load( Archive & arc ) {
	arc( rsd1_ ); // Size
	arc( rsd2_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::id::JumpID );

#endif // SERIALIZATION
