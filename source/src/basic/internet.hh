// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/internet.hh
/// @brief  Implements utility functions which deal with the internet
/// Note that the general approach of Rosetta is to NOT touch the internet,
/// unless the user explicitly opts in to internet-connection functionality.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_basic_internet_hh
#define INCLUDED_basic_internet_hh

#include <string>

namespace basic {

/// @brief Download a file from the internet at location URL, saving it to the file dest
/// Returns true on success and false on failure
bool
download_file( std::string const & URL, std::string const & dest );

/// @brief Takes a function pointer to a function which downloads a file from a URL to a destination
/// (Provided mainly so that PyRosetta can set its own downloader.)
void
set_file_downloader( bool func(std::string const &, std::string const &) );

/////////////////////////////////////////////////////////////

/// @brief Use a system call to wget to download the file
bool wget_downloader( std::string const & URL, std::string const & dest );


}

#endif // INCLUDED_basic_internet_hh
