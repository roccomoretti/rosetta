// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/RNA_SecStruct.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_pose_rna_RNA_SecStruct_FWD_HH
#define INCLUDED_core_pose_rna_RNA_SecStruct_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pose {
namespace rna {

class RNA_SecStruct;
typedef utility::pointer::shared_ptr< RNA_SecStruct > RNA_SecStructOP;
typedef utility::pointer::shared_ptr< RNA_SecStruct const > RNA_SecStructCOP;

} //secstruct
} //farna
} //protocols

#endif
