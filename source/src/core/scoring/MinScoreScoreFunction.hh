// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/symmetry/DockingScoreFunction.hh
/// @brief  Symmetric Score function class
/// @author Ingemar Andre


#ifndef INCLUDED_core_scoring_MinScoreScoreFunction_hh
#define INCLUDED_core_scoring_MinScoreScoreFunction_hh

// Unit headers
#include <core/scoring/MinScoreScoreFunction.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// Project headers

#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

class MinScoreScoreFunction : public ScoreFunction
{
public:
	typedef ScoreFunction parent;

public:

	/// ctor
	MinScoreScoreFunction();

	/// ctor
	MinScoreScoreFunction( utility::options::OptionCollection const & options );

private:

	MinScoreScoreFunction &
	operator=( MinScoreScoreFunction const & );

	MinScoreScoreFunction( MinScoreScoreFunction const & );

public:

	MinScoreScoreFunction( ScoreFunction const & src, core::Real const );
	MinScoreScoreFunction( ScoreFunction const & src, core::Real const, utility::options::OptionCollection const & options );

	MinScoreScoreFunction( core::Real const );
	MinScoreScoreFunction( core::Real const, utility::options::OptionCollection const & options );

	/// @brief INTERNAL USE ONLY
	void
	assign( ScoreFunction const & src) override;

	/// @brief INTERNAL USE ONLY
	virtual void
	assign( MinScoreScoreFunction const & src);

	ScoreFunctionOP clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// score
	/////////////////////////////////////////////////////////////////////////////

	Real
	operator ()( pose::Pose & pose ) const override;

	static
	void
	list_options_read( utility::options::OptionKeyList & options_read );

private:
	core::Real min_score_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_MinScoreScoreFunction )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_MinScoreScoreFunction_HH
