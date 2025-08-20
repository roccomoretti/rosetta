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

static const std::vector< std::string > FEATURE_NAMES = {
	"HYDRO",
	"ELEM",
	"GEOM",
	"NUMH",
	"NBOND",
};

class ResidueFeaturizer {
public:
	static
	core::Size
	column_index(std::string const & colname) {
		if ( colname == "HYDRO" ) {
			return 0;
		} else if ( colname == "ELEM" ) {
			return 1;
		} else if ( colname == "GEOM" ) {
			return 2;
		} else if ( colname == "NUMH" ) {
			return 3;
		} else if ( colname == "NBOND" ) {
			return 4;
		} else {
			utility_exit_with_message("Cannot interpret `"+colname+"` as a column name");
		}
	}

	static
	std::string
	column_from_index(core::Size index) {
		switch( index ) {
		case 0:
			return "HYDRO";
		case 1:
			return "ELEM";
		case 2:
			return "GEOM";
		case 3:
			return "NUMH";
		case 4:
			return "NBOND";
		default:
			return std::to_string(index);
		}
	}

	static
	utility::vector1< int >
	parse_value(std::string const & value_str, std::string const & colname) {
		if ( value_str.empty() ) {
			if (colname == "ELEM") {
				return {1,2,3,4,5,6,7,8,9};
			} else if ( colname == "GEOM" ) {
				return {0,1,2};
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
		}
		// Direct numeric
		try {
			return {std::stoi(value_str)};
		} catch ( std::invalid_argument const & e ) {
			utility_exit_with_message("Cannot parse value `"+value_str+"` for column `"+colname+"`");
		}
	}

	static
	std::string
	value_to_string(int value, core::Size col_index) {
		switch( col_index ) {
		//case 0:
		//	return "HYDRO";
		case 1:
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
		case 2:
			switch( value ) {
				case 0: return "TET";
				case 1: return "TRI";
				case 2: return "LIN";
			}
			break;
		//case 3:
		//	return "NUMH";
		//case 4:
		//	return "NBOND";
		default:
			// fall off
			break;
		}
		return std::to_string(value);
	}

	class FeatureSet {

	public:
		FeatureSet() = default; // Needed for std::map default construction

		FeatureSet(core::chemical::ResidueType const & restype) {
			for ( core::Size atm(1); atm <= restype.natoms(); ++atm ) {
				if ( restype.is_virtual(atm) ) { continue; }

				std::vector<int> af; // atom features

				af.push_back( int(restype.atom_is_hydrogen(atm)) );
				af.push_back( get_element(restype, atm) );
				af.push_back( get_geom(restype, atm) );
				af.push_back( get_nhydro(restype, atm) );
				af.push_back( get_bonded(restype, atm) );

				features_[atm] = af;
			}
		}

		std::map< core::Size, utility::vector0<int> > const &
		get_feature_vector() const {
			return features_;
		}

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

		int
		get_geom(core::chemical::ResidueType const & restype, core::Size atm_in) {
			static constexpr int TET = 0;
			static constexpr int TRI = 1;
			static constexpr int LIN = 2;

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
				return LIN;
			}
			if ( n_aro > 0 ) {
				return TRI;
			}
			if ( n_double == 0 ) {
				// TODO: Deal with COO, CON, etc. (and the connections)
				return TET; // all single bondes
			}

			auto elem = restype.element(atm); // base for hydrogens.
			if ( elem == core::chemical::element::P || elem == core::chemical::element::S ) {
				if ( n_double >= 2 ) {
					return TET; // Phoshpate, sulfate
				} else {
					return TRI;
				}
			} else {
				if (n_double >= 2) {
					return LIN; // C=C=C
				} else {
					return TRI;
				}
			}
			return TET; // should never get here
		}

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
		// Atom number to feature list
		std::map< core::Size, utility::vector0<int> > features_;
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

class FeatureSpec {

public:
	struct Condition {
		core::Size which = 0;
		core::Size index = 0;
		int value = -1;
		bool inv = false;

		std::string to_string() const {
			return ResidueFeaturizer::column_from_index(index) + ":" + std::to_string(which) + (inv?"!=":"=") + ResidueFeaturizer::value_to_string(value,index);
		}
	};


	FeatureSpec( utility::vector1< Condition > const & conditions = utility::vector1< Condition >{} ):
		conditions_(conditions)
	{}

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

	static
	utility::vector1< FeatureSpec >
	parse_json( json const & config ) {
		utility::vector1< FeatureSpec > feature_specs;

		for ( auto const & entry: config["features"] ) {
			// Other non-restrict entries reserved for future use

			utility::vector1< FeatureSpec > condition_sets;
			condition_sets.push_back( FeatureSpec{} );

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
				core::Size column = ResidueFeaturizer::column_index(column_str);
				for (core::Size pos: parse_pos(pos_str) ) {
					utility::vector1< Condition > conditions_for_position;
					for ( int val: ResidueFeaturizer::parse_value(value_str, column_str) ) {
						Condition new_cond;
						new_cond.which = pos;
						new_cond.index = column;
						new_cond.value = val;
						new_cond.inv = inv;

						conditions_for_position.push_back( new_cond );
					}
					if ( condition_sets.empty() ) {
						condition_sets.push_back( FeatureSpec{} );
					}
					utility::vector1< FeatureSpec > new_condition_sets;
					for ( auto const & old_conds: condition_sets ) {
						for ( auto const & new_cond: conditions_for_position ) {
							FeatureSpec new_conds( old_conds );
							new_conds.conditions_.push_back( new_cond );
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

	utility::vector1< Condition >::const_iterator
	begin() const { return conditions_.begin(); }

	utility::vector1< Condition >::const_iterator
	end() const { return conditions_.end(); }

	std::vector< std::string >
	to_string_vector() const {
		std::vector< std::string > retval;
		for ( auto const & cond: conditions_ ) {
			retval.push_back( cond.to_string() );
		}
		return retval;
	}

	bool
	matches(ResidueFeaturizer::FeatureSet const & fs, core::Size atm, core::Size which) const {
		for ( Condition const & cond: conditions_ ) {
			if ( cond.which != which ) { continue; }
			auto const & feature_vector = fs.get_feature_vector();
			if ( ! feature_vector.count(atm) ) { return false; } // Non-featurized atoms, like virtual atoms

			if ( cond.inv ) {
				if ( feature_vector.at(atm)[cond.index] == cond.value ) { return false; }
			} else {
				if ( feature_vector.at(atm)[cond.index] != cond.value ) { return false; }
			}
		}
		return true;
	}

private:

	utility::vector1< Condition > conditions_;

};

/// Contain
class FeatureDescriber {
public:

	FeatureDescriber( json const & config ) {
		features_ = FeatureSpec::parse_json( config );
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
		json output;
		output["binning"] = json::object();
		output["binning"]["dist_min"] = dist_min_;
		output["binning"]["dist_max"] = dist_max_;
		output["binning"]["dist_width"] = dist_width_;

		json features_out = json::array();
		for ( auto const & fs: features_ ) {
			features_out.push_back( fs.to_string_vector() );
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

	utility::vector1< bool >
	matches(ResidueFeaturizer::FeatureSet const & fs, core::Size atm, core::Size which) const {
		utility::vector1< bool > retval;
		for ( FeatureSpec const & spec: features_ ) {
			retval.push_back( spec.matches(fs, atm, which) );
		}
		return retval;
	}

private:
	core::Real dist_min_ = 0, dist_max_ = 10, dist_width_ = 0.2;

	utility::vector1< FeatureSpec > features_;

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
			ResidueFeaturizer::FeatureSet const & fs_ii = residue_featurizer_.featurize( ii_res.type() );

			for( core::Size jj: other ) {
				core::conformation::Residue const & jj_res = pose.residue(jj);
				if ( ii_res.nbr_atom_xyz().distance( jj_res.nbr_atom_xyz() ) > (dist_max + ii_res.nbr_radius() + jj_res.nbr_radius()) ) {
					continue;
				}
				ResidueFeaturizer::FeatureSet const & fs_jj = residue_featurizer_.featurize( jj_res.type() );

				for ( core::Size ai(1); ai <= ii_res.natoms(); ++ai ) {
					utility::vector1< bool > ii_atom_feat0 = feature_describer_.matches( fs_ii, ai, 0 );
					utility::vector1< bool > ii_atom_feat1 = feature_describer_.matches( fs_ii, ai, 1 );

					for ( core::Size aj(1); aj <= jj_res.natoms(); ++aj ) {
						core::Real dist = ii_res.xyz(ai).distance( jj_res.xyz(aj) );
						if ( dist > dist_max ) { continue; } // Too far

						utility::vector1< bool > jj_atom_feat0 = feature_describer_.matches( fs_jj, aj, 0 );
						utility::vector1< bool > jj_atom_feat1 = feature_describer_.matches( fs_jj, aj, 1 );
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
		out << outvalues.dump( 2 );
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

		ResidueFeaturizer featurizer;

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
