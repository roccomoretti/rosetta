// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/rmoretti/resscore.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <devel/init.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>

#include <core/id/AtomID.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <utility/vector0.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>

#include <string>
#include <map>

static basic::Tracer TR("apps.resscore_featurize");


static const std::vector< std::string > FEATURE_NAMES = {
	"HYDRO",
	"ELEM_C",
	"ELEM_N",
	"ELEM_O",
	"ELEM_P",
	"ELEM_S",
	"ELEM_F",
	"ELEM_Cl",
	"ELEM_Br",
	"ELEM_I",
	"ELEM_X",
	"GEOM_tetra",
	"GEOM_trig",
	"GEOM_lin",
	"NUMH_H0",
	"NUMH_H1",
	"NUMH_H2",
	"NUMH_H3",
	"NUMH_H4+",
	"NBOND_DX",
	"NBOND_D1",
	"NBOND_D2",
	"NBOND_D3",
	"NBOND_D4"
};

class ResidueFeaturizer {
public:

	class FeatureSet {
	public:
		FeatureSet() = default; // Needed for std::map default construction

		FeatureSet(core::chemical::ResidueType const & restype) {
			for ( core::Size atm(1); atm <= restype.natoms(); ++atm ) {
				if ( restype.is_virtual(atm) ) { continue; }

				utility::vector0<bool> af; // atom features

				af.push_back( restype.atom_is_hydrogen(atm) );
				af.append( get_element(restype, atm) );
				af.append( get_geom(restype, atm) );
				af.append( get_conn(restype, atm) );
				af.append( get_bonded(restype, atm) );

				features_[atm] = af;
			}
		}

		std::map< core::Size, utility::vector0<bool> > const &
		get_feature_vector() const {
			return features_;
		}

		utility::vector0<bool>
		get_element(core::chemical::ResidueType const & restype, core::Size atm) {
			using namespace core::chemical::element;

			auto elem = restype.element(atm);
			if ( restype.atom_is_hydrogen(atm) ) { // The element for a hydrogen is the attached element type
				elem = restype.element( restype.atom_base(atm) );
			}

			utility::vector0<bool> ef;
			ef.push_back( elem == C );
			ef.push_back( elem == N );
			ef.push_back( elem == O );
			ef.push_back( elem == P );
			ef.push_back( elem == S );
			ef.push_back( elem == F );
			ef.push_back( elem == Cl );
			ef.push_back( elem == Br );
			ef.push_back( elem == I );
			ef.push_back( std::none_of( ef.begin(), ef.end(), [](bool b){ return b; } ) ); //

			return ef;
		}

		utility::vector0<bool>
		get_geom(core::chemical::ResidueType const & restype, core::Size atm_in) {

			core::Size atm = atm_in;
			if ( restype.atom_is_hydrogen(atm_in) ) {
				atm = restype.atom_base(atm_in);
			}

			auto types = restype.bonded_neighbor_types(atm);

			core::Size n_double = 0, n_aro = 0, n_triple = 0;
			for ( auto type: restype.bonded_neighbor_types(atm) ) {
				switch ( type ) {
				case core::chemical::BondName::TripleBond:
					++n_triple;
					break;
				case core::chemical::BondName::DoubleBond:
					++n_double;
					break;
				case core::chemical::BondName::AromaticBond:
					++n_aro;
					break;
				default:
					break;
				}
			}

			if ( n_triple > 0 ) {
				return {false, false, true}; // lin
			}
			if ( n_aro > 0 ) {
				return {false, true, false}; // tri
			}
			if ( n_double == 0 ) {
				// TODO: Deal with COO, CON, etc. (and the connections)
				return {true, false, false}; // tet -- all single bonds
			}

			auto elem = restype.element(atm); // base for hydrogens.
			if ( elem == core::chemical::element::P || elem == core::chemical::element::S ) {
				if ( n_double >= 2 ) {
					return {true, false, false}; // Phoshpate, sulfate
				} else {
					return {false, true, false}; // tri
				}
			} else {
				if (n_double >= 2) {
					return {false, false, true}; // lin: C=C=C
				} else {
					return {false, true, false}; // tri
				}
			}
			return {true, false, false}; // tet -- should never get here
		}

		utility::vector0<bool>
		get_conn(core::chemical::ResidueType const & restype, core::Size atm_in) {
			core::Size atm = atm_in;
			if ( restype.atom_is_hydrogen(atm_in) ) {
				atm = restype.atom_base(atm_in);
			}
			core::Size nhydro = restype.number_bonded_hydrogens(atm);

			utility::vector0<bool> cf = {false, false, false, false, false};

			if ( nhydro >= cf.size() ) {
				cf[ cf.size() - 1 ] = true;
			} else {
				cf[ nhydro ] = true;
			}
			return cf;
		}

		utility::vector0<bool>
		get_bonded(core::chemical::ResidueType const & restype, core::Size atm_in) {
			core::Size atm = atm_in;
			if ( restype.atom_is_hydrogen(atm_in) ) {
				atm = restype.atom_base(atm_in);
			}
			auto bonded_to = restype.bonded_neighbor(atm);

			core::Size nbonded = 0;
			for ( auto b: bonded_to ) {
				if ( ! restype.is_virtual(b) ) {
					++nbonded;
				}
			}

			utility::vector0<bool> bf = {false, false, false, false, false};

			if ( nbonded >= bf.size() ) {
				bf[ 0 ] = true;
			} else {
				bf[ nbonded ] = true;
			}

			return bf;
		}

	private:
		std::map< core::Size, utility::vector0<bool> > features_;
	};

	ResidueFeaturizer() = default;

	FeatureSet const & featurize(core::chemical::ResidueType const & restype) {
		std::string rtn = restype.name();
		if ( features_.count(rtn) == 0 ) {
			features_[rtn] = FeatureSet(restype);
		}
		return features_[rtn];
	}

private:

	std::map< std::string, FeatureSet > features_;

};


class PoseFeaturizer {
	core::Real const DISTANCE_CUTOFF = 10.0;

public:
	PoseFeaturizer(core::pose::Pose const & pose, ResidueFeaturizer & featurizer) {
		for ( core::Size jj(1); jj <= pose.size(); ++jj ) {
			core::conformation::Residue const & jj_res = pose.residue(jj);

			collect_features( jj_res, featurizer );

			// ii in inner loop such that the featurization above happens before.
			for ( core::Size ii(1); ii < jj; ++ii ) {
				core::conformation::Residue const & ii_res = pose.residue(ii);

				// Quick out if residues are too far apart
				if ( (ii_res.nbr_atom_xyz().distance( jj_res.nbr_atom_xyz() ) - ii_res.nbr_radius() - jj_res.nbr_radius()) > DISTANCE_CUTOFF ) {
					break;
				}

				core::Real respair_idx = grab_energies(pose, ii, jj);

				for ( core::Size ai(1); ai <= ii_res.natoms(); ++ai ) {
					core::Real ai_idx = get_atom_index(ii, ai);
					if ( ai_idx < 0 ) { continue; }

					for ( core::Size aj(1); aj <= jj_res.natoms(); ++aj ) {
						core::Real aj_idx = get_atom_index(jj, aj);
						if ( aj_idx < 0 ) { continue; }

						// DIST info
						core::Real dist = ii_res.xyz(ai).distance( jj_res.xyz(aj) );
						if ( dist > DISTANCE_CUTOFF ) { continue; } // Too far

						// Bit better coherence -- always put hydrogens second
						if ( ii_res.atom_is_hydrogen(ai) && ! jj_res.atom_is_hydrogen(aj) ) {
							bonds_.push_back( {respair_idx, aj_idx, ai_idx, dist} );
						} else {
							// TODO: Should we permute the order here?
							bonds_.push_back( {respair_idx, ai_idx, aj_idx, dist} );
						}

						add_angles( respair_idx, ii_res, ai, jj_res, aj, true );
						add_angles( respair_idx, jj_res, aj, ii_res, ai, false );
					}
				}
			}
		}
	}

	void
	collect_features(core::conformation::Residue const & res, ResidueFeaturizer & featurizer ) {
		ResidueFeaturizer::FeatureSet const & feat = featurizer.featurize( res.type() );
		for ( auto const & pair: feat.get_feature_vector() ) {
			id_to_atomno_[ core::id::AtomID(pair.first, res.seqpos()) ] = features_.size(); // pre-push gives zero-indexed
			features_.push_back( pair.second );
		}
	}

	core::Size
	grab_energies(core::pose::Pose const & pose, core::Size ii, core::Size jj) {
		core::scoring::Energies const & pose_energies( pose.energies() );
		core::scoring::EnergyGraph const & egraph( pose_energies.energy_graph() );

		core::Size idx = energies_.size();

		respair_to_idx_[ std::make_pair(ii, jj) ] = idx;
		respair_to_idx_[ std::make_pair(jj, ii) ] = idx;

		auto edge = egraph.find_energy_edge(ii, jj);

		core::Real score = 0;
		if ( edge ) {
			score = edge->fill_energy_map().dot( pose_energies.weights() );
		}

		energies_.push_back( score );
		return idx;
	}

	core::Real
	get_atom_index(core::Size resnum, core::Size atomnum) {
		core::id::AtomID id(atomnum, resnum);
		if ( id_to_atomno_.count(id) == 0 ) {
			return -1; // Not featurized, for some reason
		} else {
			return core::Real( id_to_atomno_[id] );
		}
	}

	void
	add_angles( core::Real respair_idx, core::conformation::Residue const & res, core::Size aa, core::conformation::Residue const & res2, core::Size bb, bool sym = false) {
		for ( core::Size base: res.bonded_neighbor(aa) ) {
			if ( res.atom_is_hydrogen(base) || res.is_virtual(base) ) { continue; }
			core::Real base_idx = get_atom_index(res.seqpos(), base);
			if ( base_idx < 0 ) { continue; }

			core::Real angle = numeric::angle_degrees( res.xyz(base), res.xyz(aa), res2.xyz(bb) );

			core::Real aa_idx = get_atom_index(res.seqpos(), aa); // Already validated
			core::Real bb_idx = get_atom_index(res2.seqpos(), bb);
			angles_.push_back( {respair_idx, base_idx, aa_idx, bb_idx, angle} );

			// Asymmetric torsions
			for ( core::Size base2: res.bonded_neighbor(base) ) {
				if ( base2 == aa ) { continue; }
				if ( res.atom_is_hydrogen(base2) || res.is_virtual(base2) ) { continue; }
				core::Real base2_idx = get_atom_index(res.seqpos(), base2);
				if ( base2_idx < 0 ) { continue; }

				core::Real torsion = numeric::dihedral_degrees( res.xyz(base2), res.xyz(base), res.xyz(aa), res2.xyz(bb) );
				asym_torsions_.push_back( {respair_idx, base2_idx, base_idx, aa_idx, bb_idx, torsion} );
			}

			if ( sym ) {
				// Symmetric torsions -- one on each residue
				for ( core::Size baseb: res2.bonded_neighbor(bb) ) {
					if ( res2.atom_is_hydrogen(baseb) || res2.is_virtual(baseb) ) { continue; }
					core::Real baseb_idx = get_atom_index(res2.seqpos(), baseb);
					if ( baseb_idx < 0 ) { continue; }

					core::Real torsion = numeric::dihedral_degrees( res.xyz(base), res.xyz(aa), res2.xyz(bb), res2.xyz(baseb) );
					if ( res.atom_is_hydrogen(aa) && ! res2.atom_is_hydrogen(bb) ) {
						sym_torsions_.push_back( {respair_idx, baseb_idx, bb_idx, aa_idx, base_idx, torsion} );
					} else {
						sym_torsions_.push_back( {respair_idx, base_idx, aa_idx, bb_idx, baseb_idx, torsion} );
					}
				}

			}
		}
	}

	template< class T >
	void
	dump_table( utility::io::ozstream & out, utility::vector0< utility::vector0< T > > const & vector ) {
		for ( auto const & inner: vector ) {
			for ( auto const & value: inner ) {
				out << value << "\t";
			}
			out << "\n";
		}
	}

	template< class T >
	void
	dump_vector( utility::io::ozstream & out, utility::vector0< T > const & vector ) {
		for ( auto const & value: vector ) {
			out << value << "\n";
		}
	}

	void
	dump( std::string const & filename ) {
		utility::io::ozstream out(filename);

		out << "##FEATURES";
		for (std::string const & name: FEATURE_NAMES ) {
			out << "\t" << name;
		}
		out << "\n";
		dump_table( out, features_ );

		out << "##RESPAIR_ENERGIES\n";
		dump_vector( out, energies_ );

		out << "##BONDS\n";
		dump_table( out, bonds_ );

		out << "##ANGLES\n";
		dump_table( out, angles_ );

		out << "##ASYMTORSIONS\n";
		dump_table( out, asym_torsions_ );

		out << "##SYMTORSIONS\n";
		dump_table( out, sym_torsions_ );
	}

private:
	std::map< core::id::AtomID, core::Size > id_to_atomno_;
	std::map< std::pair< core::Size, core::Size >, core::Size > respair_to_idx_;

	utility::vector0< utility::vector0< bool > > features_;
	utility::vector0< utility::vector0< core::Real > > bonds_; // respair, atm1, atm2, dist
	utility::vector0< utility::vector0< core::Real > > angles_; // respair, atm1, atm2, atm3, angle
	utility::vector0< utility::vector0< core::Real > > asym_torsions_; // respair, atm1, atm2, atm3, atm4, angle
	utility::vector0< utility::vector0< core::Real > > sym_torsions_; // respair, atm1, atm2, atm3, atm4, angle

	utility::vector0< core::Real > energies_;
};

int
main( int argc, char* argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try {

		devel::init( argc, argv );

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		core::scoring::methods::EnergyMethodOptionsOP emopts( new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		scorefxn->set_energy_method_options( *emopts );
		// Post-loading commandline scorefunction alterations.
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scorefxn );
		// TODO: Do we need other scorefunction modifications here?

		using namespace core::import_pose::pose_stream;
		MetaPoseInputStream input = streams_from_cmd_line();
		core::chemical::ResidueTypeSetCOP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
			basic::options::option[ in::file::residue_type_set ]()
		);

		ResidueFeaturizer featurizer;

		core::pose::Pose pose;

		while ( input.has_another_pose() ) {
			clock_t starttime = clock();

			input.fill_pose( pose, *rsd_set );

			// Load commandline pose modifications
			core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose );
			// TODO: Others that should go here? Disulfides, rotamer bonuses, etc.?

			(*scorefxn)(pose);

			PoseFeaturizer pose_featurizer(pose, featurizer);

			std::string filename = std::string(basic::options::option[ out::path::all ].value()) + utility::file::FileName( core::pose::tag_from_pose(pose) ).base() + ".dat.gz";
			TR << "Saving data to `" << filename << "`" << std::endl;

			pose_featurizer.dump(filename);

			clock_t stoptime = clock();
			TR << "Processing file " << core::pose::tag_from_pose(pose) << " took " << ((double) stoptime - starttime) / CLOCKS_PER_SEC << " seconds." << std::endl;
		} // while ( input.has_another_pose() )
		return 0;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

}
