// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/spline/CompoundInterpolator.cc
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#include <numeric/interpolation/spline/CompoundInterpolator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <utility/tools/make_vector.hh>
#include <utility/pointer/owning_ptr.hh>

#include <algorithm>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace numeric {
namespace interpolation {
namespace spline {

struct compare_interp_range {
	bool operator()( interp_range const & a, interp_range const & b ) {
		return a.lb < b.ub;
	}
};

CompoundInterpolator::CompoundInterpolator( CompoundInterpolator const & other ):
	Interpolator( other )
{
	// Make deep copy of the sub-interpolater ranges.
	for ( interp_range const & other_range: other.interpolators_ ) {
		interp_range new_range( other_range );
		new_range.interp = other_range.interp->clone();
		interpolators_.push_back( new_range );
	}
}

InterpolatorOP CompoundInterpolator::clone() const {
	return utility::pointer::make_shared< CompoundInterpolator >( *this );
}

void
CompoundInterpolator::add_range(
	InterpolatorOP interp,
	Real lb,
	Real ub
) {
	interp_range ir;
	ir.lb = lb;
	ir.ub = ub;
	ir.interp = interp;
	interpolators_.push_back( ir );
	std::sort( interpolators_.begin(), interpolators_.end(), compare_interp_range() );
	for ( size_t i = 1; i < interpolators_.size(); ++i ) {
		assert( interpolators_[i].ub <= interpolators_[i+1].lb );
	}
}


void
CompoundInterpolator::interpolate(
	Real x,
	Real & y,
	Real & dy
) const {

	if ( has_lb_function() && x < get_lb_function_cutoff() ) {
		return compute_lb_function_solution(x,y);
	}
	if ( has_ub_function() && x > get_ub_function_cutoff() ) {
		return compute_ub_function_solution(x,y);
	}

	for ( size_t i = 1; i <= interpolators_.size(); ++i ) {
		if ( interpolators_[i].lb <= x && x <= interpolators_[i].ub ) {
			return interpolators_[i].interp->interpolate(x,y,dy);
		}
	}
	assert(false);
}

/// @brief serialize the Interpolator to a json object
nlohmann::json CompoundInterpolator::serialize() const
{
	std::vector< nlohmann::json > interpolator_data;
	for ( auto const & it : interpolators_ ) {
		nlohmann::json inner {
			{ "ub", it.ub },
			{ "lb", it.lb },
			{ "interp", it.interp->serialize() }
		};
		interpolator_data.emplace_back( std::move(inner) );
	}

	nlohmann::json j {
		{ "interp_list", interpolator_data },
		{ "base_data", Interpolator::serialize() }
	};
	return j;

}

/// @brief deserialize a json object to a Interpolator
void CompoundInterpolator::deserialize(nlohmann::json const & data)
{
	interpolators_.clear();
	for ( auto & interpolator_record : data["interp_list"] ) {
		InterpolatorOP current_interpolator( new SimpleInterpolator() );
		current_interpolator->deserialize(interpolator_record["interp"]);
		add_range(current_interpolator,data["lb"].get<Real>(),data["ub"].get<Real>());
	}

	Interpolator::deserialize(data["base_data"]);
}

bool CompoundInterpolator::operator == ( Interpolator const & other ) const
{
	if ( ! Interpolator::operator==( other ) ) return false;

	auto const & other_downcast( static_cast< CompoundInterpolator const & > ( other ) );
	if ( interpolators_.size() != other_downcast.interpolators_.size() ) return false;
	for ( platform::Size ii = 1; ii <= interpolators_.size(); ++ii ) {
		if ( interpolators_[ii].lb != other_downcast.interpolators_[ii].lb ) return false;
		if ( interpolators_[ii].ub != other_downcast.interpolators_[ii].ub ) return false;
		if ( ! (*interpolators_[ii].interp == *other_downcast.interpolators_[ii].interp) ) return false;
	}
	return true;
}

bool CompoundInterpolator::same_type_as_me( Interpolator const & other ) const
{
	return dynamic_cast< CompoundInterpolator const * > (&other);
}

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::interpolation::spline::interp_range::save( Archive & arc ) const {
	arc( CEREAL_NVP( lb ) ); // Real
	arc( CEREAL_NVP( ub ) ); // Real
	arc( CEREAL_NVP( interp ) ); // InterpolatorOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::interpolation::spline::interp_range::load( Archive & arc ) {
	arc( lb ); // Real
	arc( ub ); // Real
	arc( interp ); // InterpolatorOP
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::interpolation::spline::interp_range );

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::interpolation::spline::CompoundInterpolator::save( Archive & arc ) const {
	arc( cereal::base_class< class numeric::interpolation::spline::Interpolator >( this ) );
	arc( CEREAL_NVP( interpolators_ ) ); // utility::vector1<interp_range>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::interpolation::spline::CompoundInterpolator::load( Archive & arc ) {
	arc( cereal::base_class< class numeric::interpolation::spline::Interpolator >( this ) );
	arc( interpolators_ ); // utility::vector1<interp_range>
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::interpolation::spline::CompoundInterpolator );
CEREAL_REGISTER_TYPE( numeric::interpolation::spline::CompoundInterpolator )

CEREAL_REGISTER_DYNAMIC_INIT( numeric_interpolation_spline_CompoundInterpolator )
#endif // SERIALIZATION

