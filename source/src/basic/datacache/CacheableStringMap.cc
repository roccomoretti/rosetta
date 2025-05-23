// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableStringMap.cc
/// @brief
/// @author Phil Bradley

// unit headers
#include <basic/datacache/CacheableStringMap.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

CacheableStringMap::CacheableStringMap() : CacheableData() {}

CacheableStringMap::~CacheableStringMap() = default;

CacheableDataOP
CacheableStringMap::clone() const {
	return utility::pointer::make_shared< CacheableStringMap >(*this);
}

std::map< std::string, std::string > &
CacheableStringMap::map() {
	return map_;
}

std::map< std::string, std::string > const &
CacheableStringMap::map() const {
	return map_;
}

CacheableStringMapOP
CacheableStringMap::shared_from_this(){
	return utility::pointer::static_pointer_cast<CacheableStringMap>( CacheableData::shared_from_this() );
}


} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::CacheableStringMap::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( map_ ) ); // std::map<std::string, std::string>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::CacheableStringMap::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( map_ ); // std::map<std::string, std::string>
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::CacheableStringMap );
CEREAL_REGISTER_TYPE( basic::datacache::CacheableStringMap )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_CacheableStringMap )
#endif // SERIALIZATION
