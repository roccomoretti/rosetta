// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IdentifyLigandMotifs.hh
/// @brief .hh fils for IdentifyLigandMotifs protocol. Protocol object reads in a pdb file(s) and outputs motifs for protein-ligand interactions in .pdb and .motifs format. App originally written by mdsmith, optimized and converted to protocol by Ari Ginsparg.

// libRosetta headers
//#include <core/options/option.hh>

#include <sstream>
#include <string>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/motifs/MotifLibrary.fwd.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <utility/vector1.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/motifs/Motif.fwd.hh>

class IdentifyLigandMotifs
{

public:

	// @brief default constructor
	IdentifyLigandMotifs();

	// @brief destructor
	~IdentifyLigandMotifs();

	// @brief process inputted pdb file(s), returns a MotifLibrary full of the motifs
	void
	process_file_list();

	// @brief function that looks at a residue (in our case a ligand) and returns vectors containing index lists of adjacent atoms within the residue
	void get_atom_trios(
		utility::vector1< utility::vector1< core::Size > > & motif_indices_list,
		utility::vector1< utility::vector1< core::Size > > & all_motif_indices,
		core::conformation::Residue const & ligres
	);

	// @brief returns the motif_library_ within the class, in case it is needed for additional usage beyond the scope of this protocol
	protocols::motifs::MotifLibrary
	get_motif_library();

	// @brief writes the motif_library_ to a file, using the motif_file_output_ as the file name prefix
	void
	write_motifs_to_disk();

private:
	//protected:

	// @brief get hbond score of ligand-residue interaction
	core::Real
	get_hbond_score(
		core::pose::Pose & pose,
		core::Size pos1,
		core::Size pos2,
		core::scoring::ScoreFunction & scorefxn
	);

	// @brief write a good motif interaction out to MotifLibrary
	void
	output_single_motif_to_MotifLibrary(
		core::pose::Pose & src_pose,
		protocols::motifs::MotifLibrary & motifs,
		std::string & pdb_name,
		core::Size prot_pos,
		core::Size lig_pos,
		utility::vector1< core::Size > & lig_atoms,
		core::Real pack_score,
		core::Real hb_score
	);

	// @brief write a good motif interaction out to pdb
	void
	output_single_motif_to_pdb(
		core::pose::Pose & src_pose,
		std::string & pdb_name,
		core::Size prot_pos,
		core::Size lig_pos,
		utility::vector1< core::Size > & lig_atoms
	);

	// @brief get packing score of ligand-residue interaction
	core::Real get_packing_score(
		core::pose::Pose & pose,
		core::Size pos1,
		core::Size pos2,
		core::scoring::ScoreFunction & scorefxn
	);

	// @brief main function, iterate over ligand (broken into adjacent 3 atom trios) and identify how they interact with nearby residues
	void
	process_for_motifs(
		core::pose::Pose & pose,
		std::string & pdb_name,
		protocols::motifs::MotifLibrary & motifs
	);

	// @brief function called in process_for_motifs() that collects any potential motifs between a ligand and amino acid residue, identified by index within passed pose; pdb_name is used for optional output of identified pdbs
	void
	ligand_to_residue_analysis(
		core::Size lig_pos,
		core::Size prot_pos,
		core::pose::Pose & pose,
		std::string & pdb_name,
		protocols::motifs::MotifLibrary & motifs,
		core::scoring::ScoreFunctionOP scorefxn,
		utility::vector1< utility::vector1< core::Size > > & motif_indices_list
	);

	// @brief string to control the directory that motif pdbs will be outputted to
	std::string motif_pdb_output_path_;
	// @brief sring to control the file that motifs will be written to
	std::string motif_file_output_;
	// @brief bool to control whether to output motifs to a MotifLibrary and .motifs file
	bool output_motifs_;
	// @brief bool to control whether to output motifs to a pdb
	bool output_motifs_as_pdb_;
	// @brief MotifLibrary to stay with the class
	protocols::motifs::MotifLibrary motif_library_;
};
