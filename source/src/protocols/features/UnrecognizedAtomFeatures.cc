// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/UnrecognizedAtomFeatures.cc
/// @brief  report unrecognized atoms features to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/UnrecognizedAtomFeatures.hh>

//External

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.hh>

//Basic Headers
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/Tracer.hh>

#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>


// Utility Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>

// External Headers
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/UnrecognizedAtomFeaturesCreator.hh>

namespace protocols {
namespace features {

static basic::Tracer TR( "protocols.features.UnrecognizedAtomFeatures" );

using std::string;
using std::map;
using std::endl;
using std::pair;
using core::Size;
using core::Real;
using core::Distance;
using core::pose::UnrecognizedAtomRecord;
using core::pose::Pose;
using core::pose::PDBInfoCOP;
using core::conformation::Residue;
using basic::datacache::DataMap;
using utility::tag::TagCOP;
using utility::vector1;
using utility::tools::make_vector;
using utility::sql_database::sessionOP;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

UnrecognizedAtomFeatures::UnrecognizedAtomFeatures() :
	neighbor_distance_cutoff_(12.0)
{
	if ( !basic::options::option[basic::options::OptionKeys::in::remember_unrecognized_res]() || !basic::options::option[basic::options::OptionKeys::in::remember_unrecognized_water]() ) {
		TR.Warning << "Use -in:remember_unrecognized_res and -in:remember_unrecognized_water to locate unrecognized atoms." << endl;
	}
}

UnrecognizedAtomFeatures::UnrecognizedAtomFeatures(
	Real neighbor_distance_cutoff
) :
	neighbor_distance_cutoff_(neighbor_distance_cutoff)
{}

void
UnrecognizedAtomFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_unrecognized_residues_table_schema(db_session);
	write_unrecognized_atoms_table_schema(db_session);
	write_unrecognized_neighbors_table_schema(db_session);
}

void
UnrecognizedAtomFeatures::write_unrecognized_residues_table_schema(
	sessionOP db_session
) const {

	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", utility::pointer::make_shared< DbBigInt >());
	Column residue_number("residue_number", utility::pointer::make_shared< DbInteger >());
	Column name3("name3", utility::pointer::make_shared< DbText >());
	Column max_temperature("max_temperature", utility::pointer::make_shared< DbReal >());

	Columns residues_pkey_cols;
	residues_pkey_cols.push_back(struct_id);
	residues_pkey_cols.push_back(residue_number);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	vector1< string > reference_columns;
	reference_columns.push_back("struct_id");
	ForeignKey foreign_key(foreign_key_columns, "structures", reference_columns, true);

	Schema table("unrecognized_residues", PrimaryKey(residues_pkey_cols));
	table.add_column(name3);
	table.add_column(max_temperature);
	table.add_foreign_key(foreign_key);
	table.write(db_session);
}


void
UnrecognizedAtomFeatures::write_unrecognized_atoms_table_schema(
	sessionOP db_session
) const {

	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", utility::pointer::make_shared< DbBigInt >());
	Column residue_number("residue_number", utility::pointer::make_shared< DbInteger >());
	Column atom_name("atom_name", utility::pointer::make_shared< DbText >());
	Column coord_x("coord_x", utility::pointer::make_shared< DbReal >());
	Column coord_y("coord_y", utility::pointer::make_shared< DbReal >());
	Column coord_z("coord_z", utility::pointer::make_shared< DbReal >());
	Column temperature("temperature", utility::pointer::make_shared< DbReal >());

	Columns residues_pkey_cols;
	residues_pkey_cols.push_back(struct_id);
	residues_pkey_cols.push_back(residue_number);
	residues_pkey_cols.push_back(atom_name);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	vector1< string > reference_columns;
	reference_columns.push_back("struct_id");
	ForeignKey foreign_key(foreign_key_columns, "structures", reference_columns, true);

	Schema table("unrecognized_atoms", PrimaryKey(residues_pkey_cols));
	table.add_column(coord_x);
	table.add_column(coord_y);
	table.add_column(coord_z);
	table.add_column(temperature);
	table.add_foreign_key(foreign_key);
	table.write(db_session);
}

void
UnrecognizedAtomFeatures::write_unrecognized_neighbors_table_schema(
	sessionOP db_session
) const {

	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", utility::pointer::make_shared< DbBigInt >(), false);
	Column residue_number("residue_number", utility::pointer::make_shared< DbInteger >(), false);
	Column unrecognized_residue_number("unrecognized_residue_number", utility::pointer::make_shared< DbReal >(), false);
	Column closest_contact("closest_contact", utility::pointer::make_shared< DbReal >(), false);

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(residue_number);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(residue_number);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("resNum");
	ForeignKey foreign_key(foreign_key_columns, "residues", reference_columns, true);
	// No foreign key for neighbor_residue_number, because as there is
	// currently just unrecognized atoms, not unrecognized residues.

	Schema table("unrecognized_neighbors", primary_key);
	table.add_column(unrecognized_residue_number);
	table.add_column(closest_contact);
	table.add_foreign_key(foreign_key);
	table.write(db_session);
}


utility::vector1<std::string>
UnrecognizedAtomFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

void
UnrecognizedAtomFeatures::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/
) {
	neighbor_distance_cutoff_ = tag->getOption<Real>("neighbor_distance_cutoff", 12.0);
}

Size
UnrecognizedAtomFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	insert_unrecognized_residues_rows(pose, struct_id, db_session);
	insert_unrecognized_atoms_rows(pose, struct_id, db_session);
	insert_unrecognized_neighbors_rows(pose, relevant_residues, struct_id, db_session);
	return 0;
}


void
UnrecognizedAtomFeatures::insert_unrecognized_residues_rows(
	Pose const & pose,
	StructureID const struct_id,
	sessionOP db_session
){
	PDBInfoCOP pdb_info(pose.pdb_info());
	if ( !pdb_info ) return;


	InsertGenerator insert_generator("unrecognized_residues");
	insert_generator.add_column("struct_id");
	insert_generator.add_column("residue_number");
	insert_generator.add_column("name3");
	insert_generator.add_column("max_temperature");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id", struct_id) );

	map< core::Size, UnrecognizedAtomRecord const * > ur_found;

	for ( auto const & ua : pdb_info->get_unrecognized_atoms() ) {
		map< core::Size, UnrecognizedAtomRecord const * >::const_iterator
			i(ur_found.find(ua.res_num()));
		if ( i == ur_found.end() ) {
			ur_found[ua.res_num()] = &ua;
		} else {
			if ( i->second->temp() > ua.temp() ) {
				ur_found[ua.res_num()] = &ua;
			}
		}
	}

	for (
			map<core::Size, UnrecognizedAtomRecord const * >::const_iterator
			i = ur_found.begin(), ie = ur_found.end();
			i != ie; ++i ) {

		UnrecognizedAtomRecord const & ua(*(i->second));

		insert_generator.add_row(
			make_vector(
			struct_id_data,
			RowDataBaseOP(utility::pointer::make_shared< RowData<core::Size> >("residue_number", ua.res_num()) ),
			RowDataBaseOP(utility::pointer::make_shared< RowData<string> >("name3", ua.res_name()) ),
			RowDataBaseOP(utility::pointer::make_shared< RowData<Real> >("max_temperature", ua.temp()) )));
	}

	insert_generator.write_to_database(db_session);
}


void
UnrecognizedAtomFeatures::insert_unrecognized_atoms_rows(
	Pose const & pose,
	StructureID const struct_id,
	sessionOP db_session
){
	PDBInfoCOP pdb_info(pose.pdb_info());
	if ( !pdb_info ) return;


	InsertGenerator insert_generator("unrecognized_atoms");
	insert_generator.add_column("struct_id");
	insert_generator.add_column("residue_number");
	insert_generator.add_column("atom_name");
	insert_generator.add_column("coord_x");
	insert_generator.add_column("coord_y");
	insert_generator.add_column("coord_z");
	insert_generator.add_column("temperature");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id", struct_id) );

	for ( UnrecognizedAtomRecord const & ua : pdb_info->get_unrecognized_atoms() ) {

		insert_generator.add_row(
			make_vector(
			struct_id_data,
			RowDataBaseOP(utility::pointer::make_shared< RowData<core::Size> >("residue_number", ua.res_num()) ),
			RowDataBaseOP(utility::pointer::make_shared< RowData<string> >("atom_name", ua.atom_name()) ),
			RowDataBaseOP(utility::pointer::make_shared< RowData<Real> >("coord_x", ua.coords().x()) ),
			RowDataBaseOP(utility::pointer::make_shared< RowData<Real> >("coord_y", ua.coords().y()) ),
			RowDataBaseOP(utility::pointer::make_shared< RowData<Real> >("coord_z", ua.coords().z()) ),
			RowDataBaseOP(utility::pointer::make_shared< RowData<Real> >("temperature", ua.temp()) )));

	}

	insert_generator.write_to_database(db_session);
}

void
UnrecognizedAtomFeatures::insert_unrecognized_neighbors_rows(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	PDBInfoCOP pdb_info(pose.pdb_info());
	if ( !pdb_info ) return;


	InsertGenerator insert_generator("unrecognized_neighbors");
	insert_generator.add_column("struct_id");
	insert_generator.add_column("residue_number");
	insert_generator.add_column("unrecognized_residue_number");
	insert_generator.add_column("closest_contact");

	RowDataBaseOP struct_id_data( new RowData<StructureID>("struct_id", struct_id) );

	map< core::Size, pair< core::Size, Real > > closest_contact;


	for ( core::Size resNum=1; resNum <= pose.size(); ++resNum ) {
		if ( !check_relevant_residues(relevant_residues, resNum) ) continue;
		Residue const & res(pose.residue(resNum));

		core::Size closest_ua_resNum(0);
		Distance closest_ua_distance(10000000);
		for ( UnrecognizedAtomRecord const & ua : pdb_info->get_unrecognized_atoms() ) {
			Distance const ua_distance(res.actcoord().distance(ua.coords()));
			if (
					ua_distance < neighbor_distance_cutoff_ &&
					ua_distance < closest_ua_distance ) {
				closest_ua_resNum = ua.res_num();
				closest_ua_distance = ua_distance;
			}
		}
		if ( closest_ua_distance < neighbor_distance_cutoff_ ) {
			insert_generator.add_row(
				make_vector(
				struct_id_data,
				RowDataBaseOP( utility::pointer::make_shared< RowData<core::Size> >("residue_number", resNum) ),
				RowDataBaseOP( utility::pointer::make_shared< RowData<core::Size> >("unrecognized_residue_number", closest_ua_resNum) ),
				RowDataBaseOP( utility::pointer::make_shared< RowData<Distance> >("closest_contact", closest_ua_distance) )));
		}
	}

	insert_generator.write_to_database(db_session);
}


void
UnrecognizedAtomFeatures::delete_record(
	StructureID struct_id,
	sessionOP db_session) {

	delete_records_from_table("unrecognized_atoms", struct_id, db_session);
	delete_records_from_table("unrecognized_neighbors", struct_id, db_session);
}

std::string UnrecognizedAtomFeatures::type_name() const {
	return class_name();
}

std::string UnrecognizedAtomFeatures::class_name() {
	return "UnrecognizedAtomFeatures";
}

void UnrecognizedAtomFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("neighbor_distance_cutoff", xsct_real, "Neighbor distance cutoff", "12.0");

	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "report unrecognized atoms features to features Statistics Scientific Benchmark", attlist );
}

std::string UnrecognizedAtomFeaturesCreator::type_name() const {
	return UnrecognizedAtomFeatures::class_name();
}

protocols::features::FeaturesReporterOP
UnrecognizedAtomFeaturesCreator::create_features_reporter() const {
	return utility::pointer::make_shared< UnrecognizedAtomFeatures >();
}

void UnrecognizedAtomFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	UnrecognizedAtomFeatures::provide_xml_schema( xsd );
}



} // namesapce
} // namespace
