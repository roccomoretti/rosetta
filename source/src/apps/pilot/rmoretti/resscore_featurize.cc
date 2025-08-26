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

#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
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
OPT_KEY( File, config )

enum class FeatureType {
	HYDRO,
	ELEM,
	GEOM,
	NUMH,
	NBOND
};

FeatureType
FeatureType_from_string(std::string const & type_str) {
	if ( type_str == "HYDRO" ) {
		return FeatureType::HYDRO;
	} else if ( type_str == "ELEM" ) {
		return FeatureType::ELEM;
	} else if ( type_str == "GEOM" ) {
		return FeatureType::GEOM;
	} else if ( type_str == "NUMH" ) {
		return FeatureType::NUMH;
	} else if ( type_str == "NBOND" ) {
		return FeatureType::NBOND;
	} else {
		utility_exit_with_message("Cannot interpret `"+type_str+"` as a feature type name");
	}
}

std::string
FeatureType_to_string(FeatureType ft) {
	switch( ft ) {
	case FeatureType::HYDRO:
		return "HYDRO";
	case FeatureType::ELEM:
		return "ELEM";
	case FeatureType::GEOM:
		return "GEOM";
	case FeatureType::NUMH:
		return "NUMH";
	case FeatureType::NBOND:
		return "NBOND";
	}
	utility_exit_with_message("Can't understand feature specification");
}

utility::vector1< int >
parse_feature_value(std::string const & value_str, std::string const & colname) {
	if ( value_str.empty() ) {
		if (colname == "ELEM") {
			return {1,2,3,4,5,6,7,8,9,0};
		} else if ( colname == "GEOM" ) {
			return {0,1,2,3};
		} else if ( colname == "NUMH" ) {
			return {0,1,2,3,4};
		} else if ( colname == "NBOND" ) {
			return {1,2,3,4,0};
		} else {
			utility_exit_with_message("Cannot column `"+colname+"`");
		}
	}
	if ( colname == "ELEM" ) {
		if ( value_str == "X" ) { return {0}; }
		else if ( value_str == "C" ) { return {1}; }
		else if ( value_str == "N" ) { return {2}; }
		else if ( value_str == "O" ) { return {3}; }
		else if ( value_str == "S" ) { return {4}; }
		else if ( value_str == "P" ) { return {5}; }
		else if ( value_str == "F" ) { return {6}; }
		else if ( value_str == "Cl" ) { return {7}; }
		else if ( value_str == "Br" ) { return {8}; }
		else if ( value_str == "I" ) { return {9}; }
	} else if ( colname == "GEOM" ) {
		if ( value_str == "TET" ) { return {0}; }
		else if ( value_str == "TRI" ) { return {1}; }
		else if ( value_str == "LIN" ) { return {2}; }
		else if ( value_str == "UNK" ) { return {3}; }
	}
	// Direct numeric
	try {
		return {std::stoi(value_str)};
	} catch ( std::invalid_argument const & e ) {
		utility_exit_with_message("Cannot parse value `"+value_str+"` for column `"+colname+"`");
	}
}

std::string
feature_value_to_string(int value, FeatureType ft) {
	switch( ft ) {
	case FeatureType::ELEM:
		switch( value ) {
			case 0: return "X";
			case 1: return "C";
			case 2: return "N";
			case 3: return "O";
			case 4: return "S";
			case 5: return "P";
			case 6: return "F";
			case 7: return "Cl";
			case 8: return "Br";
			case 9: return "I";
		}
		break;
	case FeatureType::GEOM:
		switch( value ) {
			case 0: return "TET";
			case 1: return "TRI";
			case 2: return "LIN";
			case 3: return "UNK";
		}
		break;
	default:
		// The rest are numeric
		break;
	}
	return std::to_string(value);
}

utility::vector1< core::Size >
parse_pos( std::string const & pos_designation, std::string const & type = "DIST" ) {
	if ( ! pos_designation.empty() ) {
		return utility::vector1< core::Size >{ core::Size( std::stol( pos_designation ) ) };
	} else if ( type == "DIST" ) {
		return {1, 2};
	} else {
		utility_exit_with_message("Unable to parse position designation `" + pos_designation + "` of type " + type );
	}
}

class AtomFeaturizer {

public:
	AtomFeaturizer() = default;

	static
	bool
	skip( core::chemical::ResidueType const & restype, core::Size atm ) {
		// For now, just skip virtual atoms
		return restype.is_virtual(atm);
	}

	static
	int
	get_value( core::chemical::ResidueType const & restype, core::Size atm, FeatureType ft ) {
		switch ( ft ) {
		case FeatureType::HYDRO:
			return int(restype.atom_is_hydrogen(atm));
		case FeatureType::ELEM:
			return get_element(restype, atm);
		case FeatureType::GEOM:
			return get_geom(restype,atm);
		case FeatureType::NUMH:
			return get_nhydro(restype,atm);
		case FeatureType::NBOND:
			return get_bonded(restype,atm);
		}
		utility_exit_with_message("Can't understand feature type");
	}

	static
	int
	get_element(core::chemical::ResidueType const & restype, core::Size atm) {
		using namespace core::chemical::element;

		auto elem = restype.element(atm);
		if ( restype.atom_is_hydrogen(atm) ) { // The element for a hydrogen is the attached element type
			elem = restype.element( restype.atom_base(atm) );
		}

		switch (elem) {
		case C:
			return 1;
		case N:
			return 2;
		case O:
			return 3;
		case S:
			return 4;
		case P:
			return 5;
		case F:
			return 6;
		case Cl:
			return 7;
		case Br:
			return 8;
		case I:
			return 9;
		default:
			return 0;
		}
	}

	static
	int
	get_geom(core::chemical::ResidueType const & restype, core::Size atm_in) {
		static constexpr int TET = 0;
		static constexpr int TRI = 1;
		static constexpr int LIN = 2;
		static constexpr int UNK = 3;

		core::Size atm = atm_in;
		if ( restype.atom_is_hydrogen(atm_in) ) {
			atm = restype.atom_base(atm_in);
		}

		auto types = restype.bonded_neighbor_types(atm);
		if ( types.size() == 0 ) {
			return UNK;
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
			return LIN;
		}
		if ( n_aro > 0 ) {
			return TRI;
		}

		if ( n_double == 0 ) {
			// Check for amide, carboxylate, aromatic amine, etc.
			if ( has_lone_pair_and_attached_to_delocalizable_pi(restype, atm) ) {
					return TRI;
			}
			return TET; // all single bondes
		} else if ( n_double == 1 ) {
			return TRI;
		} else { // n_double >= 2
			auto elem = restype.element(atm); // base for hydrogens.
			if ( elem == core::chemical::element::P || elem == core::chemical::element::S ) {
				return TET; // Phoshpate, sulfate
			} else {
				return LIN; // C=C=C
			}
		}
		return TET; // should never get here
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
	get_bonded(core::chemical::ResidueType const & restype, core::Size atm_in) {
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

private:
};

struct Condition {
	core::Size pos = 0;
	FeatureType feature_type = FeatureType::HYDRO;
	int value = -1;
	bool inv = false;

	std::string to_string() const {
		return FeatureType_to_string(feature_type) + ":" + std::to_string(pos) + (inv?"!=":"=") + feature_value_to_string(value,feature_type);
	}

	bool
	is_same(Condition const & other, bool check_pos=true) const {
		if ( check_pos && pos != other.pos ) { return false; }
		return feature_type == other.feature_type && value == other.value && inv == other.inv;
	}

};

class FeatureSpec {

public:

	FeatureSpec( int group=999 ):
		group_(group)
	{}

	void
	add_condition( Condition const & cond ) {
		for ( auto const & curr_cond: conditions_by_pos_[cond.pos] ) {
			if ( cond.is_same( curr_cond ) ) {
				return; // Don't double-add a condition.
			}
		}
		conditions_by_pos_[cond.pos].push_back( cond );
	}

// Format of feature specifications
// features: [
// 	{
// 		restrict: ["COLUMN"] ## All values for all positions
// 	}
// 	{
// 		restrict: ["COLUMN:POS"] ## All values for the given pos
// 	}
// 	{
// 		restrict: ["COLUMN=VAL"] ## Specific value for all columns
// 	}
//	{ restrict: ["COLUMN:POS=VAL","COLUMN:POS!=VAL",...] }
//	...
// ]

	static
	utility::vector1< FeatureSpec >
	parse_json( json const & config ) {
		utility::vector1< FeatureSpec > feature_specs;

		for ( auto const & entry: config["features"] ) {
			// Other non-restrict entries reserved for future use

			int group = entry.value("group",999);

			utility::vector1< FeatureSpec > condition_sets;
			condition_sets.emplace_back( group );

			for ( auto const & cond: entry["restrict"] ) {
				std::string column_str = "";
				std::string pos_str = "";
				std::string value_str = "";
				bool inv = false;

				utility::vector1<std::string> split_eq = utility::string_split(cond, '=');
				if ( split_eq.size() == 0 ) { continue; }
				if ( split_eq.size() > 2 ) { utility_exit_with_message("Too many '=' in `restrict`"); }
				if ( split_eq.size() == 2) {
					value_str = split_eq[2];
				}
				std::string colpos = split_eq[1];
				if ( colpos[ colpos.size()-1 ] == '!' ) {
					inv = true;
					colpos = colpos.substr(0, colpos.size()-1);
				}
				utility::vector1<std::string> split_colon = utility::string_split(colpos, ':');
				if ( split_colon.size() > 2 ) { utility_exit_with_message("Too many ':' in `restrict`"); }
				if ( split_colon.size() == 2 ) {
					pos_str = split_colon[2];
				}
				column_str = split_colon[1];

				// Now convert the designations to Conditions
				FeatureType ft = FeatureType_from_string(column_str);
				for (core::Size pos: parse_pos(pos_str) ) {
					utility::vector1< Condition > conditions_for_position;
					for ( int val: parse_feature_value(value_str, column_str) ) {
						Condition new_cond;
						new_cond.pos = pos;
						new_cond.feature_type = ft;
						new_cond.value = val;
						new_cond.inv = inv;

						conditions_for_position.push_back( new_cond );
					}
					utility::vector1< FeatureSpec > new_condition_sets;
					for ( auto const & old_conds: condition_sets ) {
						for ( auto const & new_cond: conditions_for_position ) {
							FeatureSpec new_conds( old_conds );
							new_conds.add_condition( new_cond );
							new_condition_sets.emplace_back( std::move(new_conds) );
						}
					}
					if ( new_condition_sets.size() >= condition_sets.size() ) {
						condition_sets = std::move(new_condition_sets);
					}
				}
			}

			feature_specs.append( condition_sets );
		}

		return feature_specs;
	}

	json
	to_json() const {
		auto retval = json::object();
		retval["group"] = group_;

		std::vector< std::string > restrict;
		for ( auto const & entry: conditions_by_pos_ ) {
			for ( auto const & cond: entry.second ) {
				restrict.push_back( cond.to_string() );
			}
		}
		retval["restrict"] = restrict;

		return retval;
	}

	bool
	matches( core::chemical::ResidueType const & restype, core::Size atm, int pos ) const {
		if ( AtomFeaturizer::skip(restype, atm) ) { return false; }
		if ( conditions_by_pos_.count(pos) == 0 ) { return true; } // No conditions mean pass
		for ( Condition const & cond: conditions_by_pos_.at(pos) ) {
			int atom_val = AtomFeaturizer::get_value(restype, atm, cond.feature_type);
			if ( cond.inv ) {
				if ( atom_val == cond.value ) { return false; }
			} else {
				if ( atom_val != cond.value ) { return false; }
			}
		}
		return true;
	}

	bool
	compare_positions( FeatureSpec const & other, core::Size my_pos, core::Size other_pos ) const {
		if ( other.conditions_by_pos_.count(other_pos) == 0 && conditions_by_pos_.count(my_pos) == 0 ) {
			return true; // Not present matches
		}
		if ( ( other.conditions_by_pos_.count(other_pos) == 0 && conditions_by_pos_.count(my_pos) != 0 ) ||
			(conditions_by_pos_.count(my_pos) == 0 && other.conditions_by_pos_.count(other_pos) != 0 )
		) {
			return false;
		}
		utility::vector1< Condition > const & my_cond = conditions_by_pos_.at(my_pos);
		utility::vector1< Condition > const & o_cond = other.conditions_by_pos_.at(other_pos);
		// We make a bit of an assumption that we don't have any duplicate conditions in the list
		// so if we're the same size and have a match for each, we're matching.
		if ( my_cond.size() != o_cond.size() ) {
			return false;
		}
		for ( auto const & cond: my_cond ) {
			bool found = false;
			for ( auto const & other: o_cond ) {
				if ( cond.is_same(other, my_pos==other_pos) ) {
					found = true;
					break;
				}
			}
			if ( !found ) {
				return false;
			}
		}
		return true;
	}

	bool
	is_same( FeatureSpec const & other ) const {
		bool exact_same = true;
		for ( auto const & entry: conditions_by_pos_ ) {
			core::Size pos = entry.first;
			if ( ! compare_positions(other, pos, pos) ) {
				exact_same = false;
				break;
			}
		}
		if ( exact_same ) { return true; }

		// Now we handle distance inversions
		if ( compare_positions(other, 1, 2) && compare_positions(other, 2, 1) ) {
			return true;
		} else {
			return false;
		}
	}


private:

	int group_ = 999;

	std::map< core::Size, utility::vector1< Condition > > conditions_by_pos_;

};

class FeatureDescriber {
public:

	FeatureDescriber( json const & config ):
		original_config_( config )
	{
		for ( auto const & new_fs: FeatureSpec::parse_json( config ) ) {
			add_feature_spec(new_fs);
		}
		parse_binning(config);
	}

	void parse_binning( json const & config ) {
		if ( config.count("binnning") ) {
			auto const & binning = config["binning"];
			dist_min_ = binning.value("dist_min",dist_min_);
			dist_max_ = binning.value("dist_max",dist_max_);
			dist_width_ = binning.value("dist_width",dist_width_);
		}
	}

	void dump_config( std::string const & filename ) const {
		json output = original_config_; // Copy over everything we're not going to reset.

		output["binning"] = json::object();
		output["binning"]["dist_min"] = dist_min_;
		output["binning"]["dist_max"] = dist_max_;
		output["binning"]["dist_width"] = dist_width_;
		output["binning"]["dist_nbins"] = ndistbin();

		json features_out = json::array();
		for ( auto const & fs: features_ ) {
			features_out.push_back( fs.to_json() );
		}
		output["features"] = features_out;

		utility::io::ozstream f( filename );
		f << std::setw(4) << output << std::endl;
	}

	core::Real dist_max() const { return dist_max_; }
	core::Size nfeatures() const { return features_.size(); }
	core::Size ndistbin() const {
		return (dist_max_ - dist_min_)/dist_width_ + 0.99; // Round up
	}

	// Calculate return 0 to ndistbin()-1
	core::Size distbin(core::Real dist) const {
		return core::Size( (dist - dist_min_)/dist_width_ ); // Round down.
	}

	utility::vector1< FeatureSpec > const &
	get_features() const {
		return features_;
	}

	void
	add_feature_spec(FeatureSpec const & fs) {
		for ( auto const & old_fs: features_ ) {
			if ( fs.is_same(old_fs) ) {
				//TR << "Feature \n\t" << fs.to_json() << " is the same as\n\t" << old_fs.to_json() << std::endl;
				return;
			}
		}
		features_.push_back(fs);
	}

private:
	json original_config_;

	core::Real dist_min_ = 0, dist_max_ = 10, dist_width_ = 0.2;

	utility::vector1< FeatureSpec > features_;

};


class ResidueFeaturizer {
public:

	ResidueFeaturizer(FeatureDescriber & feature_describer):
		feature_describer_(feature_describer)
	{}

	utility::vector1< bool > const &
	get_feature_vector( core::chemical::ResidueType const & restype, core::Size atm, int pos ) {
		std::string const & name = restype.name();
		if ( pos_features_.count( name ) == 0 || pos_features_[name].count(atm) == 0 || pos_features_[name][atm].count(pos) == 0 ) {
			utility::vector1< bool > & feat_vect = pos_features_[name][atm][pos];

			for ( FeatureSpec const & spec: feature_describer_.get_features() ) {
				feat_vect.push_back( spec.matches(restype, atm, pos) );
			}
		}
		return pos_features_[name][atm][pos];
	}

private:

	FeatureDescriber & feature_describer_;

	// Indexed by restype name, atomnum & position
	std::map< std::string,
		std::map< core::Size,
			std::map< int,
				utility::vector1< bool >
			>
		>
	> pos_features_;

};


class PoseFeaturizer {

public:
	PoseFeaturizer(FeatureDescriber & fd, ResidueFeaturizer & rf):
		feature_describer_(fd),
		residue_featurizer_(rf)
	{
		features_.resize( feature_describer_.nfeatures() );
		for ( core::Size ii(1); ii <= feature_describer_.nfeatures(); ++ii ) {
			features_[ii].resize( feature_describer_.ndistbin() );
		}
	}

	void
	featurize( core::pose::Pose const & pose, utility::vector1< core::Size > const & focus, utility::vector1< core::Size > const & other ) {

		core::Real dist_max = feature_describer_.dist_max();

		for ( core::Size ii: focus ) {
			core::conformation::Residue const & ii_res = pose.residue(ii);
			core::chemical::ResidueType const & ii_type = ii_res.type();

			for( core::Size jj: other ) {
				core::conformation::Residue const & jj_res = pose.residue(jj);
				core::chemical::ResidueType const & jj_type = jj_res.type();
				if ( ii_res.nbr_atom_xyz().distance( jj_res.nbr_atom_xyz() ) > (dist_max + ii_res.nbr_radius() + jj_res.nbr_radius()) ) {
					continue;
				}

				for ( core::Size ai(1); ai <= ii_res.natoms(); ++ai ) {
					utility::vector1< bool > const & ii_atom_feat0 = residue_featurizer_.get_feature_vector(ii_type,ai,0);//feature_describer_.matches( fs_ii, ai, 0 );
					utility::vector1< bool > const & ii_atom_feat1 = residue_featurizer_.get_feature_vector(ii_type,ai,1);//feature_describer_.matches( fs_ii, ai, 1 );

					for ( core::Size aj(1); aj <= jj_res.natoms(); ++aj ) {
						core::Real dist = ii_res.xyz(ai).distance( jj_res.xyz(aj) );
						if ( dist > dist_max ) { continue; } // Too far

						utility::vector1< bool > const & jj_atom_feat0 = residue_featurizer_.get_feature_vector(jj_type,aj,0);//feature_describer_.matches( fs_jj, aj, 0 );
						utility::vector1< bool > const & jj_atom_feat1 = residue_featurizer_.get_feature_vector(jj_type,aj,1);//feature_describer_.matches( fs_jj, aj, 1 );
						core::Size distbin = feature_describer_.distbin(dist);

						for ( core::Size ff(1); ff <= feature_describer_.nfeatures(); ++ff ) {
							if ( ii_atom_feat0[ff] && jj_atom_feat1[ff] ) {
								++features_[ff][distbin];
							}
							if ( jj_atom_feat0[ff] && ii_atom_feat1[ff] ) {
								++features_[ff][distbin];
							}
						}

					}
				}

			}
		}
	}

public:

	void
	dump( std::string const & filename ) {
		utility::io::ozstream out(filename);
		json outvalues{ features_ };
		out << outvalues.dump( /*2*/ );
	}

private:
	FeatureDescriber & feature_describer_;
	ResidueFeaturizer & residue_featurizer_;

	// Indexed by feature number, then by distance binning, storing counts of interactions.
	utility::vector1< utility::vector0< int > > features_;

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
	NEW_OPT( config, "What (JSON-formatted) configuration file to use", "config.json" );

	try {

		devel::init( argc, argv );

		json config = load_config( basic::options::option[ basic::options::OptionKeys::config ] );

		FeatureDescriber feature_describer( config );
		feature_describer.dump_config( "config_out.json" );

		core::select::residue_selector::ResidueSelectorOP ligand_selector;
		if ( basic::options::option[ basic::options::OptionKeys::ligand_name3 ].user() ) {
			ligand_selector = utility::pointer::make_shared< core::select::residue_selector::ResidueNameSelector >( basic::options::option[ basic::options::OptionKeys::ligand_name3 ], true );
		} else {
			ligand_selector = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( basic::options::option[ basic::options::OptionKeys::ligand_chain ] );
		}

		using namespace core::import_pose::pose_stream;
		MetaPoseInputStream input = streams_from_cmd_line();
		core::chemical::ResidueTypeSetCOP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
			basic::options::option[ in::file::residue_type_set ]()
		);

		ResidueFeaturizer featurizer(feature_describer);

		core::pose::Pose pose;

		while ( input.has_another_pose() ) {
			clock_t starttime = clock();

			input.fill_pose( pose, *rsd_set );

			TR << "Processing " << core::pose::tag_from_pose(pose) << std::endl;

			PoseFeaturizer pose_featurizer(feature_describer, featurizer); // New one for each input

			utility::vector1< core::Size > ligand_residues = core::select::residue_selector::selection_positions( ligand_selector->apply(pose) );
			ligand_residues.resize(1); // Just the first residue
			utility::vector1< core::Size > other_residues = core::select::residue_selector::unselection_positions( ligand_selector->apply(pose) );

			pose_featurizer.featurize( pose, ligand_residues, other_residues );

			std::string filename = std::string(basic::options::option[ out::path::all ].value()) + utility::file::FileName( core::pose::tag_from_pose(pose) ).base() + ".json.gz";
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
