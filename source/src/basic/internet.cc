// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/internet.cc
/// @brief  Implements utility functions which deal with the internet
/// Note that the general approach of Rosetta is to NOT touch the internet,
/// unless the user explicitly opts in to internet-connection functionality.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <basic/internet.hh>
#include <basic/execute.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>
#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>

namespace basic {

static basic::Tracer TR("basic.internet");

bool (*DOWNLOADER_FUNCTION)( std::string const &, std::string const & ) = nullptr;

bool
download_file( std::string const & URL, std::string const & dest ) {
	if ( DOWNLOADER_FUNCTION ) {
		return DOWNLOADER_FUNCTION( URL, dest );
	} else {
#ifndef WIN32
		return wget_downloader( URL, dest );
#else
		utility_exit_with_message("Default internet downloader not supported on Windows.").
#endif
	}
}

void
set_file_downloader( bool func(std::string const &, std::string const &) ) {
	DOWNLOADER_FUNCTION = func;
}

bool
wget_downloader( std::string const & URL, std::string const & dest ) {
	// wget doesn't automatically create directories, so we need to do that manually if they don't already exist
	std::string dirname = utility::file::PathName( dest ).parent();
	if ( !dirname.empty() && ! utility::file::is_directory( dirname ) ) {
		TR << "Creating directory for downloaded file at " << dirname << std::endl;
		utility::file::create_directory_recursive( dirname );
	}

	std::string command = "wget";
	std::string message = "Downloading " + URL + " with wget";
	// From original EsmPerplexityTensorflowProtocol location: TODO: remove certificate flag after updating ssl stuff on our gitlab
	std::vector<std::string> args = {URL, "--progress=bar:force", "--no-check-certificate", "-O", dest};
	basic::ExecutionResult result = basic::execute(message, command, args, false, true );
	return result.result == 0;
}

}
