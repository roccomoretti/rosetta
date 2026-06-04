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

#include <core/chemical/AtomTypeSet.hh>

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

#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/util.hh>

#include <utility/vector0.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>

#include <json.hpp>

#include <fstream>
#include <string>
#include <map>

static basic::Tracer TR("apps.resscore_featurize");

OPT_KEY( String, ligand_chain )
OPT_KEY( String, ligand_name3 )
OPT_KEY( File, ligand_file ) // A file listing the ligands for each PDB
OPT_KEY( Real, dist_max )
OPT_KEY( Boolean, use_hydro )

class ResidueFeaturizer {
public:

	ResidueFeaturizer()
	{}

	utility::vector1< std::string >
	get_feature_names() const {
		return {
			"hydro",
			"element",
			"geom",
			"nhydro",
			"nbonded",
			"rtype",
			"pcharge"
		};
	}

	utility::vector1< std::string > const &
	get_feature_vector( core::chemical::ResidueType const & restype, core::Size atm ) {
		std::string const & name = restype.name();

		if ( atom_features_.count( name ) == 0 || atom_features_[name].count(atm) == 0 ) {
			utility::vector1< std::string > features;

			features.push_back( get_hydro(restype, atm) );
			features.push_back( get_element(restype, atm) );
			features.push_back( get_geom(restype, atm) );
			features.push_back( std::to_string( get_nhydro(restype, atm) ) );
			features.push_back( std::to_string( get_nbonded(restype, atm) ) );
			features.push_back( get_rtype(restype, atm) );
			features.push_back( get_pcharge(restype, atm) );

			atom_features_[name][atm] = std::move(features);
		}
		return atom_features_[name][atm];
	}

	static
	std::string
	get_hydro(core::chemical::ResidueType const & restype, core::Size atm) {
		if ( restype.atom_is_hydrogen(atm) ) {
			return "T";
		} else {
			return "F";
		}
	}

	static
	std::string
	get_element(core::chemical::ResidueType const & restype, core::Size atm) {
		using namespace core::chemical::element;

		auto elem = restype.element(atm);
		if ( restype.atom_is_hydrogen(atm) ) { // The element for a hydrogen is the attached element type
			elem = restype.element( restype.atom_base(atm) );
		}

		return core::chemical::element::name_from_elements( elem );
	}

	static
	std::string
	get_geom(core::chemical::ResidueType const & restype, core::Size atm_in) {
		core::Size atm = atm_in;
		if ( restype.atom_is_hydrogen(atm_in) ) {
			atm = restype.atom_base(atm_in);
		}

		auto types = restype.bonded_neighbor_types(atm);
		if ( types.size() == 0 ) {
			return "UNK";
		}

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
			return "LIN";
		}
		if ( n_aro > 0 ) {
			return "TRI";
		}

		if ( n_double == 0 ) {
			// Check for amide, carboxylate, aromatic amine, etc.
			if ( has_lone_pair_and_attached_to_delocalizable_pi(restype, atm) ) {
					return "TRI";
			}
			return "TET"; // all single bondes
		} else if ( n_double == 1 ) {
			return "TRI";
		} else { // n_double >= 2
			auto elem = restype.element(atm); // base for hydrogens.
			if ( elem == core::chemical::element::P || elem == core::chemical::element::S ) {
				return "TET"; // Phoshpate, sulfate
			} else {
				return "LIN"; // C=C=C
			}
		}
		return "TET"; // should never get here
	}

	static
	bool
	has_lone_pair_and_attached_to_delocalizable_pi(core::chemical::ResidueType const & restype, core::Size atm) {
		auto elem = restype.element(atm);
		if ( elem != core::chemical::element::O && elem != core::chemical::element::N && elem != core::chemical::element::S ) {
			return false; // Does not have lone pair
		}
		for ( core::Size nbr: restype.bonded_neighbor(atm) ) {
			auto nbr_elem = restype.element(nbr);
			if ( nbr_elem != core::chemical::element::C && elem != core::chemical::element::N ) {
				continue;
			}
			for ( auto nbr_bond_type: restype.bonded_neighbor_types(nbr) ) {
				if ( nbr_bond_type == core::chemical::BondName::DoubleBond || nbr_bond_type == core::chemical::BondName::AromaticBond ) {
					return true; // At least one delocalizable pi system
				}
			}
		}
		return false; // No delocalizable pi found
	}

	static
	int
	get_nhydro(core::chemical::ResidueType const & restype, core::Size atm_in) {
		static constexpr int MAX_HYDRO = 4;

		core::Size atm = atm_in;
		if ( restype.atom_is_hydrogen(atm_in) ) {
			atm = restype.atom_base(atm_in);
		}
		core::Size nhydro = restype.number_bonded_hydrogens(atm);

		if ( nhydro >= MAX_HYDRO ) {
			return MAX_HYDRO;
		} else {
			return nhydro;
		}
	}

	static
	int
	get_nbonded(core::chemical::ResidueType const & restype, core::Size atm_in) {
		static constexpr int MAX_BONDED = 4;
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

		if ( nbonded >= MAX_BONDED ) {
			return 0;
		} else {
			return nbonded;
		}
	}

	static
	std::string
	get_rtype(core::chemical::ResidueType const & restype, core::Size atm_in) {
		core::Size atm = atm_in;
		if ( restype.atom_is_hydrogen(atm_in) ) {
			// While hydrogens have their own Rosetta types, it's entirely dependent on what they're bonded to, which is a richer feature data
			atm = restype.atom_base(atm_in);
		}

		return restype.atom_type(atm).atom_type_name();
	}

	static
	std::string // String to handle rounding issues.
	get_pcharge(core::chemical::ResidueType const & restype, core::Size atm_in) {
		// For hydrogens, we're returning the partial charge of the hydrogen itself

		core::Real charge = restype.atom_charge(atm_in);
		charge = std::round( charge*10 ) / 10.0; // Be sure to round, rather than truncate
		std::stringstream out;
		out << std::fixed << std::setprecision(1) << charge;
		return out.str();
	}

private:

	// Indexed by restype name, atomnum & position
	std::map< std::string,
		std::map< core::Size,
			utility::vector1< std::string >
		>
	> atom_features_;

};


class PoseFeaturizer {

public:
	PoseFeaturizer(ResidueFeaturizer & rf, core::Real dist_max=10.0, bool use_hydro=false):
		residue_featurizer_(rf),
		dist_max_(dist_max),
		use_hydro_(use_hydro)
	{}

	void
	featurize( core::pose::Pose const & pose, utility::vector1< core::Size > const & focus, utility::vector1< core::Size > const & other ) {

		for ( core::Size ii: focus ) {
			core::conformation::Residue const & ii_res = pose.residue(ii);
			core::chemical::ResidueType const & ii_type = ii_res.type();

			// Compute residue-residue interactions
			for( core::Size jj: other ) {
				core::conformation::Residue const & jj_res = pose.residue(jj);
				core::chemical::ResidueType const & jj_type = jj_res.type();
				if ( jj_type.is_virtual_residue() ) { continue; }

				if ( ii_res.nbr_atom_xyz().distance( jj_res.nbr_atom_xyz() ) > (dist_max_ + ii_res.nbr_radius() + jj_res.nbr_radius()) ) {
					continue;
				}

				for ( core::Size ai(1); ai <= ii_res.natoms(); ++ai ) {
					if ( ii_res.is_virtual(ai) ) { continue; }
					if ( ! use_hydro_ && ai > ii_res.nheavyatoms() ) { break; }
					utility::vector1< std::string > const & ii_atom_feat = residue_featurizer_.get_feature_vector(ii_type,ai);

					for ( core::Size aj(1); aj <= jj_res.natoms(); ++aj ) {
						if ( jj_res.is_virtual(aj) ) { continue; }
						if ( ! use_hydro_ && aj > jj_res.nheavyatoms() ) { break; }

						core::Real dist = ii_res.xyz(ai).distance( jj_res.xyz(aj) );
						if ( dist > dist_max_ ) { continue; } // Too far

						utility::vector1< std::string > jj_atom_feat = residue_featurizer_.get_feature_vector(jj_type,aj);

						jj_atom_feat.append( ii_atom_feat );

						dist_features_.emplace_back( std::move(jj_atom_feat) );
						dist_values_.push_back( dist );
					}
				}

			}

			// Compute internal (e.g. score) values
			if ( pose.energies().energies_updated() ) {
				auto res_energies = pose.energies().residue_total_energies(ii) * pose.energies().weights();
				for ( int ii = 1; ii <= core::scoring::n_score_types; ++ii ) {
					core::scoring::ScoreType const type = core::scoring::ScoreType(ii);
					if ( res_energies[ type ] == 0.0 ) continue; // Avoid adding entry for zeros
					score_values_[ core::scoring::name_from_score_type(type) ] += res_energies[ type ];
				}
				score_values_["total"] += pose.energies().residue_total_energy(ii);
			}
		}
	}

public:

	void
	dump( std::string const & filename ) const {
		utility::io::ozstream out(filename);
		json output;
		json settings;
		settings["dist_max"] = dist_max_;
		settings["use_hydro"] = use_hydro_;
		output["settings"] = settings;

		utility::vector1<std::string> feature_names = residue_featurizer_.get_feature_names();
		feature_names.append( residue_featurizer_.get_feature_names() );
		output["dist_feature_names"] = feature_names;

		output["dist_features"] = dist_features_;
		output["dist_values"] = dist_values_;
		if ( ! score_values_.empty() ) {
			output["scores"] = score_values_;
		}
		out << output.dump( /*2*/ );
	}

private:
	ResidueFeaturizer & residue_featurizer_;

	core::Real dist_max_=10.0;
	bool use_hydro_ = false;

	utility::vector1< utility::vector1< std::string > > dist_features_;
	utility::vector1< float > dist_values_; // float because we don't need that much precision

	// By score type name
	std::map< std::string, core::Real > score_values_;

};

json
load_config(std::string const & filename) {
	std::ifstream stream(filename);
	return json::parse(stream);
}

int
main( int argc, char* argv [] ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( ligand_chain, "Which chain letter to use for analysis", "X" );
	NEW_OPT( ligand_name3, "Which three letter code to use for analysis", "" );
	NEW_OPT( ligand_file, "A file listing the ligand position for each PDB", "" );
	NEW_OPT( dist_max, "The maximum distance to consider the interaction", 10.0 );
	NEW_OPT( use_hydro, "Include interactions with hydrogen atoms?", false );

	try {

		devel::init( argc, argv );

		core::select::residue_selector::ResidueSelectorOP ligand_selector;
		if ( basic::options::option[ basic::options::OptionKeys::ligand_name3 ].user() ) {
			ligand_selector = utility::pointer::make_shared< core::select::residue_selector::ResidueNameSelector >( basic::options::option[ basic::options::OptionKeys::ligand_name3 ], true );
		} else {
			ligand_selector = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( basic::options::option[ basic::options::OptionKeys::ligand_chain ] );
		}

		std::map< std::string, std::string > ligand_mapping;
		if ( basic::options::option[ basic::options::OptionKeys::ligand_file ].user() ) {
			utility::io::izstream ligfile( basic::options::option[ basic::options::OptionKeys::ligand_file ] );
			std::string line;
			while ( ligfile.getline(line) ) {
				utility::vector1< std::string > splitline = utility::string_split( line );
				if ( splitline.empty() ) { continue; }
				if ( splitline.size() >= 2 ) {
					ligand_mapping[ splitline[1] ] = ligand_mapping[ splitline[2] ];
				} else {
					ligand_mapping[ splitline[1] ] = "";
				}
			}
		}

		using namespace core::import_pose::pose_stream;
		MetaPoseInputStream input = streams_from_cmd_line();
		core::chemical::ResidueTypeSetCOP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
			basic::options::option[ in::file::residue_type_set ]()
		);

		ResidueFeaturizer featurizer;
		core::Real dist_max = basic::options::option[ basic::options::OptionKeys::dist_max ];
		core::Real use_hydro = basic::options::option[ basic::options::OptionKeys::use_hydro ];

		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );

		core::pose::Pose pose;

		while ( input.has_another_pose() ) {
			clock_t starttime = clock();
			TR << "Processing " << core::pose::tag_from_pose(pose) << std::endl;

			input.fill_pose( pose, *rsd_set );
			std::string pose_tag = utility::file::FileName( core::pose::tag_from_pose(pose) ).base();
			sfxn->score(pose);

			core::select::residue_selector::ResidueSelectorOP my_ligand_selector = ligand_selector;
			if ( ligand_mapping.count( pose_tag ) ) {
				if ( ligand_mapping[pose_tag].empty() ) {
					TR.Warning << "Ligand correspondence for " << pose_tag << " is missing -- skipping" << std::endl;
					continue;
				}
				TR << "Using residue " << ligand_mapping[pose_tag] << " for " << pose_tag << std::endl;
				my_ligand_selector = utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >(ligand_mapping[pose_tag]);
			}

			PoseFeaturizer pose_featurizer(featurizer, dist_max, use_hydro); // New one for each input

			utility::vector1< core::Size > ligand_residues = core::select::residue_selector::selection_positions( my_ligand_selector->apply(pose) );
			ligand_residues.resize(1); // Just the first residue
			utility::vector1< core::Size > other_residues = core::select::residue_selector::unselection_positions( my_ligand_selector->apply(pose) );

			pose_featurizer.featurize( pose, ligand_residues, other_residues );

			std::string filename = std::string(basic::options::option[ out::path::all ].value()) + pose_tag + ".json.gz";
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
