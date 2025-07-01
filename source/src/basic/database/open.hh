// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/database/open.hh
/// @brief  Functions for opening database files
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_basic_database_open_hh
#define INCLUDED_basic_database_open_hh


// Utility headers
#include <utility/io/izstream.fwd.hh>

// C++ headers
#ifdef WIN32
#include <string>
#else
#include <iosfwd>
#endif

namespace basic {
namespace database {

/// @brief Open a file from the database on a provided stream
/// File is a location relative to the paths on -in:path:database
/// Local directories will not be consulted.
/// Throws a utility::excn::IOError exception if the file can't be loaded.
void
open(
	utility::io::izstream & db_stream,
	std::string const & db_file
);

/// @brief Open a file from the database and returns its contents
/// File is a location relative to the paths on -in:path:database
/// Local directories will not be consulted.
/// Throws a utility::excn::IOError exception if the file can't be loaded.
std::string
open(
	std::string const & db_file
);

/// @brief Open a file from the database and returns its contents
/// Uses the utility::io::GeneralFileManager to make sure that the specified file is only ever read from disk once.
/// File is a location relative to the paths on -in:path:database
/// Local directories will not be consulted.
/// Throws a utility::excn::IOError exception if the file can't be loaded.
std::string
cached_open(
	std::string const & db_file
);

/// @brief Full-path database file name
/// The paths of -in:path:database will be consulted to find the first existing option
/// @details NOTE: This function is provided primarily for debug output messages & interfacing with third party code which need the filename itself
/// Due to caching/lazy downloads/alternate database file handling,
/// don't use it for file existence checks -- just try to open (with warn=false) and
/// look at the `izstream::good()` results.
std::string
full_name(
	std::string const & db_file,
	bool warn = true
);

/// @brief Open a file from the local directory or the database on a provided stream
/// File is a location relative to the working dirtectory or the paths on -in:path:database
/// Local directory files are preferred.
/// Throws a utility::excn::IOError exception if the file can't be loaded.
void
open_with_local(
	utility::io::izstream & db_stream,
	std::string const & db_file
);


/// @brief Find a path to a file.
///
/// Try various combinations to locate the specific file being requested by the user.
/// (inspired by core::scoring::ScoreFunction::find_weights_file())
///
/// mute_error_if_failure returns "" if the file does not exist
/// instead of throwing a loud exception
///
/// dir like chemical/carbohydrates/linkage_conformers/
std::string
find_database_path(
	std::string const & dir,
	std::string const & filename,
	bool mute_error_if_failure = false
);

/// @brief Find a path to a file.
///
/// Try various combinations to locate the specific file being requested by the user.
/// (inspired by core::scoring::ScoreFunction::find_weights_file())
///
/// dir like chemical/carbohydrates/linkage_conformers/
std::string
find_database_path( std::string const & dir, std::string const & filename, std::string const & ext );


/// @brief Does cache file (absolute path) exist?
/// if dir_only is true, will return true if the cache file could be created.
bool
find_cache_file(
	std::string const & cache_file,
	bool dir_only
);

/// @brief Get the (absolute) path to a given cached file.
/// If source_file is given, it's the full path to the source database file that's being cached.
/// If for_writing is true, will only check that the given file would be creatable.
/// Will return an empty string if it can't find a cache file.
std::string
full_cache_name(
	std::string const & short_name,
	std::string const & source_file,
	bool for_writing
);


/// @brief Utility function for when the settings to automatically download database files is set.
/// Handle the downloading of the file in cases where it's missing.
/// Returns true on success and false on failure
bool
handle_database_download( std::string const & db_file, std::string const & db_file_full );


} // namespace database
} // namespace basic


#endif // INCLUDED_basic_io_database_open_HH
