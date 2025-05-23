// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/options/keys/OptionKeys.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_basic_options_keys_OptionKeys_hh
#define INCLUDED_basic_options_keys_OptionKeys_hh


// Utility headers

#include <platform/types.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/ResidueChainVectorOption.fwd.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringVectorOption.fwd.hh>

// This file is included by all the auto-generated options files and
// so REQUIRES all these includes.
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ResidueChainVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.hh>

namespace basic {
namespace options {


typedef  utility::options::OptionKey  OptionKey;
typedef  utility::options::BooleanOptionKey  BooleanOptionKey;
typedef  utility::options::IntegerOptionKey  IntegerOptionKey;
typedef  utility::options::RealOptionKey  RealOptionKey;
typedef  utility::options::StringOptionKey  StringOptionKey;
typedef  utility::options::FileOptionKey  FileOptionKey;
typedef  utility::options::PathOptionKey  PathOptionKey;
typedef  utility::options::BooleanVectorOptionKey  BooleanVectorOptionKey;
typedef  utility::options::IntegerVectorOptionKey  IntegerVectorOptionKey;
typedef  utility::options::RealVectorOptionKey  RealVectorOptionKey;
typedef  utility::options::ResidueChainVectorOptionKey  ResidueChainVectorOptionKey;
typedef  utility::options::StringVectorOptionKey  StringVectorOptionKey;
typedef  utility::options::FileVectorOptionKey  FileVectorOptionKey;
typedef  utility::options::PathVectorOptionKey  PathVectorOptionKey;


namespace OptionKeys {


/// @brief Lookup functors
typedef  OptionKey  KeyType;
#include <utility/keys/KeyLookup.functors.hh> // DO NOT AUTO-REMOVE (not an include, code injection)


} // namespace OptionKeys
} // namespace options
} // namespace basic


#endif // INCLUDED_basic_options_keys_OptionKeys_HH
