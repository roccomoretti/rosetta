// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RG_Energy_RNA.cc
/// @brief  Radius of gyration energy function definition.
/// @author Rhiju Das


// Unit headers
#include <core/energy_methods/RG_Energy_RNA.hh>
#include <core/energy_methods/RG_Energy_RNACreator.hh>

// Package headers
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>


// Utility headers

#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


// C++

using namespace core::chemical;
using namespace core::chemical::rna;

namespace core {
namespace energy_methods {


/// @details This must return a fresh instance of the RG_Energy_RNA class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
RG_Energy_RNACreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< RG_Energy_RNA >();
}

core::scoring::ScoreTypes
RG_Energy_RNACreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( rna_rg );
	return sts;
}


/// c-tor
RG_Energy_RNA::RG_Energy_RNA() :
	parent( utility::pointer::make_shared< RG_Energy_RNACreator >() )
{}

/// clone
core::scoring::methods::EnergyMethodOP
RG_Energy_RNA::clone() const
{
	return utility::pointer::make_shared< RG_Energy_RNA >();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
RG_Energy_RNA::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	core::scoring::rna::RNA_ScoringInfo  & rna_scoring_info( core::scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	//Doesn't recalculate stuff if already updated:
	rna_centroid_info.update( pose );

	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
	Size const nres( pose.size() );

	// calculate center of mass -- mutable.

	center_of_mass_ = 0.0;
	for ( Size i = 1; i <= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue( i ) );
		if ( !rsd.is_RNA() ) continue;
		Vector const v( base_centroids[i] );
		center_of_mass_ += v;
	}
	center_of_mass_ /= nres;

	///////////////////////////////////////
	//
	// RG SCORE
	// calculate RG based on distance from center of mass
	Real rg_squared = 0;
	for ( Size i = 1; i <= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue( i ) );
		if ( !rsd.is_RNA() ) continue;
		Vector const v( base_centroids[i] );
		rg_squared += ( v - center_of_mass_ ).length_squared();
	}

	// This definition of rg differs with the conventional definition which
	// divides by nres and not by nres-1.  For the sake of matching r++, it's
	// being left at nres-1 for now, but is a candidate for change in the near
	// future.
	rg_squared /= ( nres - 1 );
	rg_ = sqrt( rg_squared ); //Save in this class
}


/////////////////////////////////////////////////////////////////////////////
void
RG_Energy_RNA::setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & scorefxn ) const
{
	//score the pose.... that should generate the center_of_mass position, and the Rg.
	setup_for_scoring( pose, scorefxn );
}


///////////////////////////////////////////////////////////////////////////////
void
RG_Energy_RNA::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals
) const {
	using namespace conformation;
	totals[ core::scoring::rna_rg ] = rg_;

	core::scoring::rna::RNA_ScoringInfo  & rna_scoring_info( core::scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	rna_centroid_info.set_calculated( false );
}


///////////////////////////////////////////////////////////////////////////////
// Following makes the approximation that center of mass of the RNA
// will stay fixed during minimimize!
void
RG_Energy_RNA::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	using namespace conformation;

	Size const i ( atom_id.rsd() );
	Size const atom_num_i ( atom_id.atomno() );

	conformation::Residue const & rsd( pose.residue( i ) );
	if ( !rsd.has_variant_type( REPLONLY ) ) return;
	if ( !rsd.is_RNA() ) return;

	core::scoring::rna::RNA_ScoringInfo  const & rna_scoring_info( core::scoring::rna::rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_CentroidInfo const & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );

	Size const nres( pose.size() );

	//Apply force at base sidechain atom?
	//
	// Alternatively spread over all atoms in base. That's perhaps more kosher, but
	// equivalent if base is rigid.
	//
	// Or alternatively, if I set up a base centroid virtual atom with a fixed
	// geometry relative to the base first sidechain atom -- apply it there.

	if ( atom_num_i != first_base_atom_index( rsd.type() ) ) return;

	Vector const v( base_centroids[i] );
	Vector f2 = ( v - center_of_mass_ )/ ( ( nres - 1 ) * rg_ );
	Vector f1 = cross( f2, v );

	F1 += weights[ core::scoring::rna_rg ] * f1;
	F2 += weights[ core::scoring::rna_rg ] * f2;
} // eval atom derivative

core::Size
RG_Energy_RNA::version() const
{
	return 1; // Initial versioning
}


} //scoring
} //core
