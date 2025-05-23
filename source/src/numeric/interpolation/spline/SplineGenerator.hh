// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/SplineGenerator.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#ifndef INCLUDED_numeric_interpolation_spline_SplineGenerator_hh
#define INCLUDED_numeric_interpolation_spline_SplineGenerator_hh

#include <numeric/interpolation/spline/SplineGenerator.fwd.hh>

#include <numeric/interpolation/spline/Interpolator.hh>

#include <numeric/types.hh>
#include <utility/vector1.hh>
#include <map>

namespace numeric {
namespace interpolation {
namespace spline {

struct Point {
	Point( Real xin, Real yin            ) : x(xin), y(yin), dy(-12345.0), has_dy(false) {}
	Point( Real xin, Real yin, Real dyin ) : x(xin), y(yin), dy(  dyin  ), has_dy(true ) {}
	Real x;
	Real y;
	Real dy;
	bool has_dy;
};

struct LinearFunction{
	LinearFunction() : cutoff(0.0), slope(0.0),intercept(0.0) {}
	LinearFunction(Real cutoff_in, Real slope_in, Real intercept_in) : cutoff(cutoff_in), slope(slope_in),intercept(intercept_in) {}

	Real cutoff;
	Real slope;
	Real intercept;

};

class SplineGenerator {
public:

	SplineGenerator(
		Real lbx, Real lby, Real lbdy,
		Real ubx, Real uby, Real ubdy
	);

	SplineGenerator();

	~SplineGenerator();
	void add_known_value( Real x, Real y );

	void add_known_value( Real x, Real y, Real dy );

	void add_boundary_function(std::string const & tag, Real const & cutoff, Real const & slope, Real const & intercept);

	InterpolatorOP get_interpolator();

	//getters for the rest of the private data
	Real get_lbx() const
	{
		return lbx_;
	}

	Real get_lby() const
	{
		return lby_;
	}

	Real get_lbdy() const
	{
		return lbdy_;
	}

	Real get_ubx() const
	{
		return ubx_;
	}

	Real get_uby() const
	{
		return uby_;
	}

	Real get_ubdy() const
	{
		return ubdy_;
	}

	numeric::Size get_num_points() const
	{
		return points_.size();
	}

	utility::vector1<Point> const & get_points() const
	{
		return points_;
	}

	std::map<std::string,LinearFunction> const & get_boundary_functions() const
	{
		return boundary_functions_;
	}

private:

	Real lbx_,lby_,lbdy_,ubx_,uby_,ubdy_;

	utility::vector1<Point> points_;

	std::map<std::string,LinearFunction> boundary_functions_;

	InterpolatorOP interpolator_;

};

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#endif
