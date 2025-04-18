// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/constraints/FileConstraintGenerator.hh
///
/// @brief
/// @author Tom Linsky ( tlinsky@uw.edu ), Nov 2012

#ifndef INCLUDED_protocols_denovo_design_constraints_FileConstraintGenerator_hh
#define INCLUDED_protocols_denovo_design_constraints_FileConstraintGenerator_hh

// Unit Header
#include <protocols/denovo_design/constraints/FileConstraintGenerator.fwd.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>

// Package Header
#include <core/scoring/constraints/Constraint.fwd.hh>

// Proeject Header
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace constraints {

class FileConstraintGenerator : public protocols::constraint_generator::ConstraintGenerator {

public:
	FileConstraintGenerator();

	// Moves, so by-value
	FileConstraintGenerator( std::string const & filename );

	~FileConstraintGenerator() override;

	protocols::constraint_generator::ConstraintGeneratorOP
	clone() const override;

	core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & pose ) const override;

	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static std::string
	class_name();


protected:
	void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) override;

public:
	void
	set_cstfile( std::string const & filename );

protected:
	std::string
	clean_constraint_string( std::string const & cst_str ) const;

private:
	std::string filename_;

}; //class FileConstraintGenerator


} //namespace constraints
} //namespace denovo_design
} //namespace protocols


#endif // INCLUDED_protocols_denovo_design_constraints_FileConstraintGenerator_HH
