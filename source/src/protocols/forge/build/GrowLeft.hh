// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/GrowLeft.hh
/// @brief instruction to create an n-side extension
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_GrowLeft_hh
#define INCLUDED_protocols_forge_build_GrowLeft_hh

// unit headers
#include <protocols/forge/build/GrowLeft.fwd.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>

// project headers

// C++ headers



namespace protocols {
namespace forge {
namespace build {


/// @brief instruction to create an n-side extension
/// @details Use this for n-side insertions, but typically not n-terminal
///  extensions unless necessary.  It does not automatically cover the
///  additional residue on the right endpoint that needs to move during
///  n-terminal extensions due to invalid phi torsion.  For that case,
///  use the SegmentRebuild class replacing the n-terminal residue with
///  desired length+1.
class GrowLeft : public BuildInstruction {


private: // typedefs


	typedef BuildInstruction Super;


public: // typedefs


	typedef Super::Size Size;

	typedef Super::ResidueTypeSetCAP ResidueTypeSetCAP;
	typedef Super::LengthEvent LengthEvent;
	typedef Super::MoveMap MoveMap;
	typedef Super::Pose Pose;

	typedef Super::Positions Positions;
	typedef Super::String String;


public: // construct/destruct


	/// @brief default constructor
	GrowLeft();


	/// @brief constructor
	/// @param[in] pos grow an n-side extension prior to this position
	/// @param[in] ss the secondary structure desired, also defines length of extension
	/// @param[in] aa the annotated amino acid sequence desired, default is poly-alanine
	/// @param[in] rts the residue type set to use, default FA_STANDARD
	/// @remarks length of the *one-letter* aa must equal the length of ss
	GrowLeft(
		core::Size const pos,
		String const & ss,
		String const & aa = String(),
		ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
	);


	/// @brief copy constructor
	GrowLeft( GrowLeft const & rval );


	/// @brief default destructor
	~GrowLeft() override;


public: // assignment


	/// @brief copy assignment
	GrowLeft & operator =( GrowLeft const & rval );


public: // virtual constructors


	/// @brief clone this object
	BuildInstructionOP clone() const override;


public: // accessors


	/// @brief grow left from this anchor position
	/// @remarks this can change if listening to Conformation LengthEvents
	///  Use original_interval() to get the original anchor position.
	inline
	core::Size pos() const {
		return pos_;
	}


	/// @brief get secondary structure string
	inline
	String const & ss() const {
		return ss_;
	}


	/// @brief get annotated amino acid string
	inline
	String const & aa() const {
		return aa_;
	}


public: // virtual accessors


	/// @brief is the original interval storing valid information, or is empty
	///  or being used for something else?
	/// @return false, stores invalid interval
	inline
	bool original_interval_valid() const override {
		return false;
	}


	/// @brief a copy of the working range of residues specifying the modified region
	/// @remarks This can change if listening to Conformation LengthEvents
	inline
	Interval interval() const override {
		return Interval( pos_ - ss_.length(), pos_ - 1 );
	}


	/// @brief return a copy of the set of positions within the new region
	///  that were pre-existing in the original Pose prior to modify()
	/// @return An empty set -- no positions are pre-existing.
	Positions preexisting_positions() const override;


	/// @brief return a copy of the set of positions that are "new" and did
	///  not exist in the original Pose.
	/// @return A set of positions spanning the entire region -- all positions
	///  are new.
	Positions new_positions() const override;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has a defined conformation.  E.g. existing or copied residues.
	/// @return An empty set -- no positions are defined.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	Positions defined_positions() const override;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has an undefined conformation.  E.g. newly created residues.
	/// @return A set of positions spanning the entire region -- all positions
	///  are undefined.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	Positions undefined_positions() const override;


	/// @brief return a copy of the MoveMap that defines the moveable/fixed
	///  positions/dofs for this instruction
	/// @return a MoveMap with [interval.left, interval.right] bb & chi set to true
	///  at the MoveMapTorsionID level
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	MoveMap movemap() const override;


public: // virtual Conformation observer interface


	/// @brief update indexing on residue append
	void on_residue_append( LengthEvent const & event ) override;


	/// @brief update indexing on residue prepend
	void on_residue_prepend( LengthEvent const & event ) override;


	/// @brief update indexing on residue delete
	void on_residue_delete( LengthEvent const & event ) override;


public: // original positions


	/// @brief return the set of positions within the original interval that
	///  will be kept in this BuildInstruction
	/// @return An empty set -- no positions are kept.
	Positions original_kept_positions() const override;


	/// @brief return set of positions within the original interval that will
	///  be deleted in this BuildInstruction
	/// @return An empty set -- no positions are deleted.
	Positions original_deleted_positions() const override;


public: // instruction comparison


	/// @brief return set of any fixed positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	/// @return empty set if no fixed positions
	Positions original_fixed_positions() const override;


	/// @brief return set of any mutable positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	/// @return empty set if no mutable positions
	Positions original_mutable_positions() const override;


public: // virtual object descriptor


	/// @brief does this object create undefined backbone in the modified region?
	/// @return true
	inline
	bool creates_undefined_backbone() const override {
		return true;
	}


protected: // virtual Pose modification methods


	/// @brief are dependencies satisfied so that modify_impl() can complete
	///  successfully?
	/// @return always True, this BuildInstruction has no dependencies
	inline
	bool dependencies_satisfied() const override {
		return true;
	}


	/// @brief do the actual work of modifying the Pose
	void modify_impl( Pose & pose ) override;


protected: // virtual mutators


	/// @brief do the actual reset of intervals, positions, etc to initial state
	void reset_accounting_impl() override;

private: // data


	/// @brief make an n-terminal extension before this position
	/// @remarks this position can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	core::Size pos_;


	/// @brief secondary structure string, also defines length of extension
	String ss_;


	/// @brief annotated amino acid string, length of the one-letter version
	///  must be equal to length of ss
	String aa_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_GrowLeft_HH */
