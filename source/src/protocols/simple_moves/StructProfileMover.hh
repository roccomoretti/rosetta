// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/StructProfileMover.hh
/// @brief Quickly generates a structure profile

#ifndef INCLUDED_protocols_simple_moves_StructProfileMover_hh
#define INCLUDED_protocols_simple_moves_StructProfileMover_hh

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/StructProfileMover.fwd.hh>


#include <protocols/indexed_structure_store/SSHashedFragmentStore.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// C++ Headers
#include <string>
// Utility Headers
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

typedef  core::Real  Probability;

class StructProfileMover : public protocols::moves::Mover {
public:
	StructProfileMover();
	StructProfileMover(core::Real rmsThreshold,core::Real burialThreshold,core::Size consider_topN_frags,core::Real burialWt,bool only_loops,bool censorByBurial,core::Real allowed_deviation, core::Real allowed_deviation_loops,bool eliminate_background,bool psiblast_style_pssm,bool outputProfile,bool add_csts_to_pose,bool ignore_terminal_res,std::string fragment_store_path="",std::string fragment_store_format="",std::string fragment_store_compression="");
	core::Size ss_type_convert(char ss_type);
	void read_P_AA_SS_cen6();
	utility::vector1<std::string> get_closest_sequence_at_res(core::pose::Pose const & pose, core::Size res,utility::vector1<core::Real> cenList);
	utility::vector1<utility::vector1<std::string> > get_closest_sequences(core::pose::Pose const & pose,utility::vector1<core::Real> cenList, core::select::residue_selector::ResidueSubset const & subset);
	utility::vector1<utility::vector1<core::Size> >generate_counts(utility::vector1<utility::vector1<std::string> > top_frag_sequences,core::pose::Pose const & pose);
	utility::vector1<utility::vector1<core::Real> >generate_profile_score(utility::vector1<utility::vector1<core::Size> > res_per_pos,core::pose::Pose const & pose);
	utility::vector1<utility::vector1<core::Real> >generate_profile_score_wo_background(utility::vector1<utility::vector1<core::Size> > res_per_pos, utility::vector1<core::Real> cenList, core::pose::Pose const & pose);
	void save_MSAcst_file(utility::vector1<utility::vector1<core::Real> > profile_score,core::pose::Pose const & pose);
	void add_MSAcst_to_pose(utility::vector1<utility::vector1<core::Real> > profile_score,core::pose::Pose & pose);
	core::Real get_cen_deviation(std::vector<core::Real> cenListFrag,utility::vector1<core::Real> cenListModel);
	std::string censorFragByBurial(std::vector<core::Real> cenListFrag,utility::vector1<core::Real> cenListModel, std::string cenListFragSeq);
	void set_profile_save_name( std::string const & name ) { profile_save_filename_ = name; }
	utility::vector1< core::Real> calc_cenlist(core::pose::Pose const & pose);
	void apply( Pose & pose ) override;
	moves::MoverOP clone() const override { return utility::pointer::make_shared< StructProfileMover >( *this ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ) override;
	void set_residue_selector( core::select::residue_selector::ResidueSelector const & selector );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	core::Real rmsThreshold_;
	core::Real burialThreshold_;
	std::string aa_order_;
	core::Size consider_topN_frags_;
	core::Real burialWt_;
	bool psiblast_style_pssm_;
	bool outputProfile_;
	bool add_csts_to_pose_;
	core::Size cenType_;
	protocols::indexed_structure_store::SSHashedFragmentStoreOP SSHashedFragmentStoreOP_;
	typedef utility::vector1< utility::vector1< utility::vector1< Probability > > > Probability_AA_n_n;
	Probability_AA_n_n P_AA_SS_burial_;
	core::Real allowed_deviation_;
	core::Real allowed_deviation_loops_;
	bool only_loops_;
	bool censorByBurial_;
	bool eliminate_background_;
	bool ignore_terminal_res_;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;
	std::string profile_save_filename_;
	std::string fragment_store_path_;
	std::string fragment_store_format_;
	std::string fragment_store_compression_;

};


} // simple_moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_StructProfileMover_hh
