// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/FragmentStoreProvider.fwd.hh
/// @brief Abstract fragment store backing interface.
/// @author Brian D. Weitzner (bweitzner@lyell.com)
//

#ifndef INCLUDED_protocols_indexed_structure_store_FragmentStoreProvider_FWD_HH
#define INCLUDED_protocols_indexed_structure_store_FragmentStoreProvider_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols
{
namespace indexed_structure_store
{

class FragmentStoreProvider;

typedef utility::pointer::shared_ptr<FragmentStoreProvider>         FragmentStoreProviderOP;
typedef utility::pointer::shared_ptr<FragmentStoreProvider const>   FragmentStoreProviderCOP;

}
}

#endif
