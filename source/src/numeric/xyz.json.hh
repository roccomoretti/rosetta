// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/xyz.json.hh
/// @author Sam DeLuca

#ifndef INCLUDED_numeric_xyz_json_hh
#define INCLUDED_numeric_xyz_json_hh

#include <numeric/xyzVector.hh>

#include <utility/tools/make_vector.hh>

#include <json.hpp>


namespace numeric {


template<typename T>
void to_json(nlohmann::json& j, const xyzVector<T>& v) {
	j.push_back(v.x());
	j.push_back(v.y());
	j.push_back(v.z());
}

template<typename T>
void from_json(const nlohmann::json& j, xyzVector<T>& v) {
	if ( !j.is_array() || j.size() != 3 ) {
		throw std::invalid_argument("Input json object of incorrect size/type.");
	}

	v.x() = j.at(0).get<T>();
	v.y() = j.at(1).get<T>();
	v.z() = j.at(2).get<T>();
}


/// @brief Convert vector to a json array
/// @note Format is a list in the form [x,y,z]
template<typename T>
inline
nlohmann::json serialize(xyzVector<T> const & coords)
{
	nlohmann::json j = nlohmann::json::array();
	j.push_back(coords.x());
	j.push_back(coords.y());
	j.push_back(coords.z());
	return j;
}

template<typename T>
inline
xyzVector<T> deserialize(nlohmann::json const & data)
{
	xyzVector<T> coords;
	coords.x(data.at(0).get<T>());
	coords.y(data.at(1).get<T>());
	coords.z(data.at(2).get<T>());
	return coords;
}


}


#endif /* XYZ_JSON_HH_ */
