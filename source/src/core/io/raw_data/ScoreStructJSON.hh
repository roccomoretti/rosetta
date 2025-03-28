// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/raw_data/ScoreStructJSON.hh
///
/// @brief Write out only the scores in JSON format
/// @author Luki Goldschmidt

#ifndef INCLUDED_core_io_raw_data_ScoreStructJSON_hh
#define INCLUDED_core_io_raw_data_ScoreStructJSON_hh

// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/raw_data/RawStruct.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <string>

namespace core {
namespace io {
namespace raw_data {

class ScoreStructJSON : public RawStruct {

public:
	// Undefined, commenting out to fix PyRosetta build  ScoreStructJSON();
	ScoreStructJSON( std::string const & tag );

	// Not implemented for this format
	void fill_pose( core::pose::Pose & ) override {}
	void fill_pose( core::pose::Pose &, core::chemical::ResidueTypeSet const & ) override {}
	void print_conformation( std::ostream & ) const override {}
	Real get_debug_rmsd() override { return 0; }

	void print_header(
		std::ostream& /* out */,
		std::map < std::string, core::Real > const & /* score_map */,
		std::map < std::string, std::string > const & /* string_map */,
		bool /* print_sequence */ ) const override
	{
		// No header in JSON
	}

	void print_scores(
		std::ostream& out,
		std::map < std::string, core::Real > const & score_map,
		std::map < std::string, std::string > const & string_map = ( std::map < std::string, std::string > () ) ) const override;

}; // class ScoreStructJSON

} // namespace silent
} // namespace io
} // namespace core

#endif
