// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/RNA_Fragments/RNA_Fragments.fwd.hh
/// @brief  RNA_Fragments forward declarations header
/// @author Rhiju Das

#ifndef INCLUDED_core_fragment_rna_RNA_Fragments_FWD_HH
#define INCLUDED_core_fragment_rna_RNA_Fragments_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace fragment {
namespace rna {

enum SYN_ANTI_RESTRICTION {
	SYN,
	ANTI,
	ANY
};

class RNA_Fragments;

typedef utility::pointer::shared_ptr< RNA_Fragments > RNA_FragmentsOP;
typedef utility::pointer::shared_ptr< RNA_Fragments const > RNA_FragmentsCOP;

} //fragments
} //rna
} //protocols

#endif
