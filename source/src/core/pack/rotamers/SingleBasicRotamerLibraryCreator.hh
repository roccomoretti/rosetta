// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamers/SingleBasicRotamerLibraryCreator.hh
/// @brief  Class for instantiating a Basic rotlib
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_pack_rotamers_SingleBasicRotamerLibraryCreator_HH
#define INCLUDED_core_pack_rotamers_SingleBasicRotamerLibraryCreator_HH

// Package headers
#include <core/pack/rotamers/SingleResidueRotamerLibraryCreator.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.fwd.hh>

// Program header
#include <core/chemical/ResidueType.fwd.hh>

// Utility headers

// C++ headers
#include <string>

namespace core {
namespace pack {
namespace rotamers {

class SingleBasicRotamerLibraryCreator : public SingleResidueRotamerLibraryCreator {
public:
	core::pack::rotamers::SingleResidueRotamerLibraryCOP
	create( core::chemical::ResidueType const & ) const override;

	std::string keyname() const override;
};


typedef utility::pointer::shared_ptr< SingleBasicRotamerLibraryCreator > SingleBasicRotamerLibraryCreatorOP;
typedef utility::pointer::shared_ptr< SingleBasicRotamerLibraryCreator const > SingleBasicRotamerLibraryCreatorCOP;

} //namespace rotamers
} //namespace pack
} //namespace core


#endif
