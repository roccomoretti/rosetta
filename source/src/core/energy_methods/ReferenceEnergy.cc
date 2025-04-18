// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/ReferenceEnergy.hh
/// @brief  Reference energy method implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/energy_methods/ReferenceEnergy.hh>
#include <core/energy_methods/ReferenceEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the ReferenceEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
ReferenceEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	if ( options.has_method_weights( core::scoring::ref ) ) {
		return utility::pointer::make_shared< ReferenceEnergy >( options.method_weights( core::scoring::ref ), options.ordered_wat_penalty() );
	} else {
		return utility::pointer::make_shared< ReferenceEnergy >( options.ordered_wat_penalty() );
	}
}

core::scoring::ScoreTypes
ReferenceEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( ref );
	return sts;
}

ReferenceEnergy::ReferenceEnergy( Real ordered_wat_penalty ) :
	parent( utility::pointer::make_shared< ReferenceEnergyCreator >() ),
	ordered_wat_penalty_ (ordered_wat_penalty)
{}

ReferenceEnergy::ReferenceEnergy( utility::vector1< Real > const & aa_weights_in, Real ordered_wat_penalty ):
	parent( utility::pointer::make_shared< ReferenceEnergyCreator >() ),
	aa_weights_( aa_weights_in ),
	ordered_wat_penalty_ (ordered_wat_penalty)
{}

ReferenceEnergy::~ReferenceEnergy() = default;

core::scoring::methods::EnergyMethodOP
ReferenceEnergy::clone() const
{
	return utility::pointer::make_shared< ReferenceEnergy >( aa_weights_, ordered_wat_penalty_ );
}


/// This is a terrible terrible terrible hack that will do for now.
void
ReferenceEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	core::scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) || rsd.is_virtual_residue() ) {
		return;
	}

	using namespace chemical;

	// water
	if ( rsd.aa() == core::chemical::aa_h2o ) {
		if ( rsd.name() == "HOH" ) {
			emap[ core::scoring::ref ] += ordered_wat_penalty_;
		}
	}

	if ( !aa_weights_.empty() ) {
		///
		AA const & aa( rsd.aa() );
		AA const aa2 ( is_d_aminoacid(aa) ? get_l_equivalent(aa) : aa);
		if ( Size(aa2) > aa_weights_.size() ) return;

		if ( rsd.type().base_analogue() == na_rad ) {
			emap[ core::scoring::ref ] += aa_weights_[ na_rad ];
		} else if ( rsd.type().base_analogue() == na_rcy ) {
			emap[ core::scoring::ref ] += aa_weights_[ na_rcy ];
		} else if ( rsd.type().base_analogue() == na_rgu ) {
			emap[ core::scoring::ref ] += aa_weights_[ na_rgu ];
		} else if ( rsd.type().base_analogue() == na_ura ) {
			emap[ core::scoring::ref ] += aa_weights_[ na_ura ];
		} else if ( aa2 == na_lra ) { // Catch L-RNA
			emap[ core::scoring::ref ] += aa_weights_[ na_rad ];
		} else if ( aa2 == na_lrc ) {
			emap[ core::scoring::ref ] += aa_weights_[ na_rcy ];
		} else if ( aa2 == na_lrg ) {
			emap[ core::scoring::ref ] += aa_weights_[ na_rgu ];
		} else if ( aa2 == na_lur ) {
			emap[ core::scoring::ref ] += aa_weights_[ na_ura ];
		} else {
			emap[ core::scoring::ref ] += aa_weights_[ aa2 ];
		}
		//   if ( rsd.is_DNA() ) {
		//    std::cout << "using dna refE " << aa_weights_[aa] << std::endl;
		//   }
		return;
	}

	// reference weights for RNA...
	if ( rsd.is_RNA() && ( aa_weights_.size() > num_canonical_aas ) ) {
		emap[ core::scoring::ref ] += aa_weights_[ rsd.aa() ];
		return;
	}

	/// else -- use the default reference weights from r++
	if ( rsd.type().aa() > num_canonical_aas ) return;

	switch ( rsd.type().aa() ) {
	case aa_ala : emap[ core::scoring::ref ] +=  0.16; break;
	case aa_cys : emap[ core::scoring::ref ] +=  1.70; break;
	case aa_asp :
		emap[ core::scoring::ref ] += (rsd.type().has_variant_type( chemical::PROTONATED )) ? -0.262 : -0.67;
		break;
	case aa_glu :
		emap[ core::scoring::ref ] += -0.81; // (rsd.type().has_variant_type( chemical::PROTONATED )) ? -0.81 : -0.81;
		break;
	case aa_phe : emap[ core::scoring::ref ] +=  0.63; break;
	case aa_gly : emap[ core::scoring::ref ] +=  -0.17; break;
	case aa_his :
		emap[ core::scoring::ref ] += (rsd.type().has_variant_type( chemical::PROTONATED )) ? 0.288 : 0.56;
		break;
	case aa_ile : emap[ core::scoring::ref ] +=  0.24; break;
	case aa_lys :
		emap[ core::scoring::ref ] += -0.65; //(rsd.type().has_variant_type( chemical::DEPROTONATED )) ? -0.65 : -0.65;
		break;
	case aa_leu : emap[ core::scoring::ref ] +=  -0.10; break;
	case aa_met : emap[ core::scoring::ref ] +=  -0.34; break;
	case aa_asn : emap[ core::scoring::ref ] +=  -0.89; break;
	case aa_pro : emap[ core::scoring::ref ] +=  0.02; break;
	case aa_gln : emap[ core::scoring::ref ] +=  -0.97; break;
	case aa_arg : emap[ core::scoring::ref ] +=  -0.98; break;
	case aa_ser : emap[ core::scoring::ref ] +=  -0.37; break;
	case aa_thr : emap[ core::scoring::ref ] +=  -0.27; break;
	case aa_val : emap[ core::scoring::ref ] +=  0.29; break;
	case aa_trp : emap[ core::scoring::ref ] +=  0.91; break;
	case aa_tyr :
		emap[ core::scoring::ref ] += (rsd.type().has_variant_type( chemical::DEPROTONATED )) ? 0.238 : 0.51;
		break;
	default :
		utility_exit_with_message("Data consistency error in ReferenceEnergy");
		break;
	}

}

///////////////////////////////////////////////////////////////////////////////

/// @brief Returns true if passed a core::chemical::AA corresponding to a
/// D-amino acid, and false otherwise.
bool
ReferenceEnergy::is_d_aminoacid(
	core::chemical::AA const res_aa
) const {
	using namespace core::chemical;
	if ( res_aa >= aa_dal && res_aa <= aa_dty ) return true;
	return false;
}

///////////////////////////////////////////////////////////////////////////////

/// @brief When passed a d-amino acid, returns the l-equivalent.  Returns
/// aa_unk otherwise.
core::chemical::AA
ReferenceEnergy::get_l_equivalent(
	core::chemical::AA const d_aa
) const {
	using namespace core::chemical;
	if ( d_aa==aa_dal ) return aa_ala;
	else if ( d_aa==aa_dcs ) return aa_cys;
	else if ( d_aa==aa_das ) return aa_asp;
	else if ( d_aa==aa_dgu ) return aa_glu;
	else if ( d_aa==aa_dph ) return aa_phe;
	else if ( d_aa==aa_dhi ) return aa_his;
	else if ( d_aa==aa_dil ) return aa_ile;
	else if ( d_aa==aa_dly ) return aa_lys;
	else if ( d_aa==aa_dle ) return aa_leu;
	else if ( d_aa==aa_dme ) return aa_met;
	else if ( d_aa==aa_dan ) return aa_asn;
	else if ( d_aa==aa_dpr ) return aa_pro;
	else if ( d_aa==aa_dgn ) return aa_gln;
	else if ( d_aa==aa_dar ) return aa_arg;
	else if ( d_aa==aa_dse ) return aa_ser;
	else if ( d_aa==aa_dth ) return aa_thr;
	else if ( d_aa==aa_dva ) return aa_val;
	else if ( d_aa==aa_dtr ) return aa_trp;
	else if ( d_aa==aa_dty ) return aa_tyr;

	return aa_unk;
}


Real
ReferenceEnergy::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const &
) const
{
	return 0.0;
}

/// @brief ReferenceEnergy is context independent; indicates that no
/// context graphs are required
void
ReferenceEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
ReferenceEnergy::version() const
{
	return 1; // Initial versioning
}


} // scoring
} // core

