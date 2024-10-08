// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/palette/PackerPaletteCreator.hh
/// @brief  Declaration of the base class for PackerPalette factory registration and creation
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_pack_palette_PackerPaletteCreator_hh
#define INCLUDED_core_pack_palette_PackerPaletteCreator_hh

// Unit headers
#include <core/pack/palette/PackerPaletteCreator.fwd.hh>

// Package headers
#include <core/pack/palette/PackerPalette.fwd.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

#include <string>


namespace core {
namespace pack {
namespace palette {

/// @brief The PackerPaletteCreator class's responsibilities are to create
/// on demand a new instance of a PackerPalette class.
/// @details The PackerPaletteCreator must register itself with the PackerPaletteFactory
/// at load time (before main() begins) so that the PackerPaletteFactory is ready
/// to start creating PackerPalettes by the time any protocol
/// requests one.
class PackerPaletteCreator : public utility::VirtualBase
{
public:
	~PackerPaletteCreator() override {}

	/// @brief Instantiate a new PackerPalette
	virtual PackerPaletteOP create_packer_palette() const = 0;
	virtual std::string keyname() const = 0;

	/// @brief Describe the allowed XML options for a particular PackerPalette subclass.
	/// @details Pure virtual.  Must be implemented by subclasses.
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;
};

}
}
}

#endif
