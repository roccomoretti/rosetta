// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @details
///
/// @author Oliver Lange, ported from dssp.cc in rosetta++ authored by Ben Blum (bblum)
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <core/scoring/dssp/StrandPairing.hh>

// Package Headers
#include <core/scoring/dssp/util.hh>
#include <core/scoring/dssp/Dssp.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <string>
#include <vector>
#include <iostream>

// Options
#include <basic/options/option.hh> // for quick-test from run:dry_run
#include <basic/options/keys/jumps.OptionKeys.gen.hh>

#include <core/scoring/dssp/PairingsList.hh>
#include <utility/vector1.hh>
#include <boost/functional/hash.hpp> // DO NOT AUTO-REMOVE

static basic::Tracer tr( "core.scoring.dssp" );

using core::Real;
using namespace core;
using namespace basic;
using namespace ObjexxFCL;

namespace core {
namespace scoring {
namespace dssp {

///////////////////////////////////////////////////////////////
///
/// @brief Constructor for set of StrandPairing objects
///
/// @details
/// The incoming hbonds matrix is indexed by (acceptor residue,
/// donor residue).  Residues with energy less than threshold are
/// considered paired, unless they are disallowed by the array
/// called allowed (this was included to prevent helical
/// residues from being considered paired, a problem which
/// occasionally arose).
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author bblum
///
///////////////////////////////////////////////////////////////

StrandPairingSet::StrandPairingSet( pose::Pose const& pose, Real threshold ) {
	ObjexxFCL::FArray2D_float hbond_bb_pair_score;
	core::scoring::dssp::fill_hbond_bb_pair_score_dssp( pose, hbond_bb_pair_score );
	compute( hbond_bb_pair_score, threshold, pose );
}

StrandPairingSet::StrandPairingSet( FArray2_float const &hbonds,
	float threshold, pose::Pose const& pose ) {
	compute( hbonds, threshold, pose);
}

StrandPairingSet::StrandPairingSet( core::scoring::dssp::PairingList const& in_pairings ) {
	for ( auto const & in_pairing : in_pairings ) {
		add_pairing( in_pairing );
	}
}

void StrandPairingSet::add_pairing( core::scoring::dssp::Pairing const& p ) {
	add_pairing( p.Pos1(), p.Pos2(), p.is_anti(), p.Pleating() );
}

void StrandPairingSet::compute( FArray2_float const &hbonds,
	float threshold, pose::Pose const& pose ) {

	Size const nres( pose.size() );
	if ( !nres ) {
		return;
	}
	debug_assert( hbonds.l1() );
	debug_assert( hbonds.l2() );
	debug_assert( hbonds.u1() <= static_cast< int >( pose.size() ) );
	debug_assert( hbonds.u2() <= static_cast< int >( pose.size() ) );

	for ( Size i = 2; i <= nres - 1; i++ ) {
		for ( Size j = i + 1; j <= nres - 1; j++ ) {
			if ( // antiparallel bridge
					( hbonds(i,j) < threshold
					&& hbonds(j,i) < threshold )
					|| ( hbonds(i-1,j+1) < threshold
					&& hbonds(j-1,i+1) < threshold ) ) {
				Size orientation, pleating;
				get_pleating(pose, i, j, orientation, pleating);
				//    tr.Debug << "found antiparallel pair at (" << i << "," << j << ") "<<std::endl;
				add_pairing(i, j, true , pleating );
			} else if ( // parallel bridge
					( hbonds(i-1,j) < threshold
					&& hbonds(j,i+1) < threshold )
					|| ( hbonds(j-1,i) < threshold
					&& hbonds(i,j+1) < threshold ) ) {
				Size orientation, pleating;
				get_pleating(pose, i, j, orientation, pleating);
				//    tr.Debug << "found parallel pair at (" << i << "," << j << ") "<<std::endl;
				//        tr.Trace << "ben: para phil: " << orientation << std::endl;
				add_pairing(i,j, false, pleating );
			}
		}
	}
}

std::istream & operator>>( std::istream &is, StrandPairingSet &set ) {
	std::string tag;
	Size nstrand;
	is >> tag >> nstrand;
	if ( tag != "STRAND_TOPOLOGY" ) {
		tr.Trace << "failed reading STRAND_TOPOLOGY --- found instead: " << tag << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}
	for ( Size ct = 1; ct <= nstrand && is.good(); ct++ ) {
		StrandPairing sp;
		is >> sp;
		if ( is.good() ) {
			set.pairings_.push_back( sp );
		}
	}
	return is;
}

std::ostream & operator<<(std::ostream & out, const StrandPairingSet &sp) {
	out << "STRAND_TOPOLGY " << sp.pairings_.size() << std::endl;
	for ( auto const & pairing : sp.pairings_ ) {
		out << pairing << std::endl;
	}
	return out;
}

///////////////////////////////////////////////////////////////
///
/// @brief Add a new pair of bonded residues to the set
///
/// @details
/// Look for a strand pairing to extend with the given pair of
/// residues; if none exists, create a new one.
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author bblum
///
///////////////////////////////////////////////////////////////
void StrandPairingSet::add_pairing( Size res1, Size res2, bool antiparallel, Size pleating) {
	bool addnew = true;
	// tr.Trace << "add pairing " << res1 << " " << res2 << (antiparallel ? " anti " : " para ") <<  std::endl;
	for ( auto & pairing : pairings_ ) {
		if ( pairing.extend(res1,res2,antiparallel,pleating) ) {
			addnew = false;
			if ( !pairing.range_check() ) {
				tr.Error << "just inconsistently added " << res1 << "-" << res2 << " to pairing " << pairing << std::endl;
				utility_exit_with_message("Inconsistent strand pairing set encountered.");
			}
			break;
		}
	}
	if ( ! addnew ) return;

	StrandPairing add(res1,res2, antiparallel, pleating);
	runtime_assert( add.range_check() );
	bool added = false;
	for ( auto it = pairings_.begin(), end = pairings_.end();
			it != end;
			++it ) {
		if ( add < *it ) {
			added = true;
			pairings_.insert( it, add );
			break;
		}
	}
	if ( !added ) {
		pairings_.push_back( add );
	}
}

bool StrandPairingSet::merge(const StrandPairingSet &other, bool domerge) {
	auto it = pairings_.begin(), end = pairings_.end();
	auto oit = other.pairings_.begin();
	while ( it != end ) {
		Size count = 0;
		while ( oit != other.pairings_.end() && it->merge(*oit, domerge) ) {
			++oit;
			++count;
		}
		if ( count == 0 ) {
			return false;
		}
		++it;
	}

	if ( oit != other.pairings_.end() ) {
		return false;
	}

	if ( domerge ) selfmerge();

	return true;
}

void StrandPairingSet::selfmerge() {
	StrandPairings goodpairings(pairings_);
	for (  auto it = pairings_.begin(), end = pairings_.end();
			it != end;
			++it ) {
		for ( auto oit = pairings_.begin();
				oit != it;
				++oit ) {
			while ( oit != it && it->merge(*oit, true) )
					oit = pairings_.erase( oit );
		}
	}
}

bool StrandPairingSet::check_pleat() const {
	for ( auto const & pairing : pairings_ ) {
		if ( !pairing.check_pleat() ) {
			return false;
		}
	}
	return true;
}

bool StrandPairing::check_pleat() const {
	for ( Size i = 1; i < (Size)pleating1_.size(); i++ ) {
		if ( pleating1_[i] == pleating1_[i-1] && pleating1_[i] != 0 ) {
			return false;
		}
	}
	return true;
}

core::Size StrandPairing::contact_order() const {
	if ( antiparallel() ) {
		return end2_-begin1_;
	} else {
		return begin2_-begin1_;
	}
	return 0;
}

///////////////////////////////////////////////////////////////
///
/// @brief If possible, extend this pairing by the given residues.
///
/// @details
/// If one of res1 or res2 is within 2 residues of the beginning
/// or end of one of the strands of the pairing, and the other
/// is within 5 residues, extend the pairing.  This is the dssp
/// definition of allowable beta bulges.  Return true if the
/// pairing was extended.
// Assumes we are running through res1 and res2 in an ordered
// way, so we extend at the beginning or end of the strand.
///
/// @param
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author bblum
///
///////////////////////////////////////////////////////////////
bool StrandPairing::extend( Size res1, Size res2, bool antiparallel, Size pleating ) {
	StrandPairing old_copy( *this );
	// Make sure res1 < res2
	Size temp = std::min(res1, res2);
	res2 = std::max(res1, res2);
	res1 = temp;
	if ( begin1_ == 0 ) { // uninitialized
		begin1_ = end1_ = res1;
		begin2_ = end2_ = res2;
		antipar_ = antiparallel;
		pairing1_.push_back( antipar_ ? end2_ : begin2_);
		pairing2_.push_back( antipar_ ? end1_ : begin1_ );
		pleating1_.push_back(pleating);
	}

	const Size SMALL_BULGE_LIMIT( basic::options::option[ basic::options::OptionKeys::jumps::max_strand_gap_allowed] );
	const Size BIG_BULGE_LIMIT( SMALL_BULGE_LIMIT + 3 );

	bool cando = false;
	if ( res1 >= begin1_ && res1 <= end1_ ) { // trivial--already in our strand
		tr.Trace <<"trivial case " << std::endl;
		cando = (pairing1_[res1 - begin1_] == res2);
	} else if ( res1 > end1_ && res1 <= end1_ + BIG_BULGE_LIMIT ) {
		tr.Trace << "case1 " << std::endl;
		if ( antiparallel ) {
			if ( res1 > end1_ + SMALL_BULGE_LIMIT ) {
				cando = (res2 < begin2_ && res2 + SMALL_BULGE_LIMIT >= begin2_ );
			} else {
				cando = (res2 < begin2_ && res2 + BIG_BULGE_LIMIT >= begin2_ );
			}
		} else {
			if ( res1 > end1_ + SMALL_BULGE_LIMIT ) { //bulge of size > 1
				cando = (res2 > end2_ && res2 <= end2_ + SMALL_BULGE_LIMIT);
			} else {
				cando = (res2 > end2_ && res2 <= end2_ + BIG_BULGE_LIMIT);
			}
		}
	} else if ( res1 < begin1_ && res1 + BIG_BULGE_LIMIT >= begin1_ ) {
		tr.Trace << "case2 "<< std::endl;
		if ( antiparallel ) {
			if (  res1 + SMALL_BULGE_LIMIT < begin1_ ) {
				cando = (res2 > end2_ && res2 <= end2_ + SMALL_BULGE_LIMIT);
			} else {
				cando = (res2 > end2_ && res2 <= end2_ + BIG_BULGE_LIMIT);
			}
		} else {
			if ( res1 + SMALL_BULGE_LIMIT < begin1_ ) { //bulge of size > 1
				cando = (res2 < begin2_ && res2 + SMALL_BULGE_LIMIT >= begin2_);
			} else {
				cando = (res2 < begin2_ && res2 + BIG_BULGE_LIMIT >= begin2_ );
			}
		}
	}
	runtime_assert( begin1_ <= end1_ );
	runtime_assert( begin2_ <= end2_ );

	// if extendable, insert this pairing in, adjust begins and ends
	if ( cando ) {
		tr.Trace << "extend " << *this << "   to residues " << res1 << " " << res2 << std::endl;
		if ( res1 < begin1_ ) {
			if ( res1 + 1 < begin1_ ) {
				pairing1_.insert(pairing1_.begin(),begin1_ - res1 - 1, 0);
				pleating1_.insert(pleating1_.begin(),begin1_ - res1 - 1, 0);
			}
			pairing1_.insert(pairing1_.begin(), res2);
			pleating1_.insert(pleating1_.begin(), pleating);
			begin1_ = res1;
		} else if ( res1 > end1_ ) {
			if ( res1 > end1_ + 1 ) {
				//add res1-end1_-1 0 to vector
				pairing1_.insert(pairing1_.end(), res1 - end1_ - 1, 0);
				pleating1_.insert(pleating1_.end(), res1 - end1_ - 1, 0);
			}
			pairing1_.push_back( res2 );
			pleating1_.push_back( pleating );
			end1_ = res1;
		}

		if ( res2 < begin2_ ) {
			if ( res2 + 1 < begin2_ ) {
				pairing2_.insert(pairing2_.begin(),begin2_ - res2 - 1, 0);
			}
			pairing2_.insert(pairing2_.begin(), res1);
			begin2_ = res2;
		} else if ( res2 > end2_ ) {
			if ( res2 > end2_ + 1 ) {
				pairing2_.insert(pairing2_.end(), res2 - end2_ - 1, 0);
			}
			pairing2_.insert(pairing2_.end(), res1);
			end2_ = res2;
		}
	} else {
		tr.Trace << " cannot extend "<< *this << "   to residues " << res1 << " " << res2 << std::endl;
	}

	show_internals(tr.Trace );
	if ( !valid_ends() ) {
		*this = old_copy;
		return false;
	}
	return cando;
}

void StrandPairing::extend_to(Size res) {
	Size res1, res2, pleat, diff = (antipar_ ? -1 : 1);
	if ( res < begin1_ ) {
		res1 = begin1_ - 1;
		res2 = pairing1_[0] - diff;
		pleat = 3 - pleating1_[0];
		while ( res1 >= res ) {
			extend(res1, res2, antipar_, pleat);
			res1 = res1 - 1;
			res2 = res2 - diff;
			pleat = 3 - pleat;
		}
	} else if ( res > end1_ ) {
		res1 = end1_ + 1;
		res2 = pairing1_[end1_ - begin1_] + diff;
		pleat = 3 - pleating1_[end1_ - begin1_];
		while ( res1 <= res ) {
			extend(res1, res2, antipar_, pleat);
			res1 += 1;
			res2 += diff;
			pleat = 3 - pleat;
		}
	}
}

void StrandPairing::show_internals( std::ostream& out ) const {
	out << "pairing1_: ";
	for ( unsigned long i : pairing1_ ) {
		out << i << " ";
	}
	out <<"\npairing2_: ";
	for ( unsigned long i : pairing2_ ) {
		out << i << " ";
	}
	out <<"\npairing2_: ";
	for ( unsigned long i : pleating1_ ) {
		out << i << " ";
	}
	out << std::endl;
}

std::size_t StrandPairing::hash_value() const {
	std::ostringstream str;
	str << (antipar_ ? 'A' : 'P') << '_';
	core::Size regA, regE;
	if ( antipar_ ) {
		regA=begin1_ + end2_;
		regE=end1_ + begin2_;
	} else {
		regA=begin2_-begin1_;
		regE=end2_-end1_;
	}
	str << regA << '_' << regE;
	return boost::hash_value(str.str());
}

bool StrandPairing::mergeable( const StrandPairing &other ) const {
	tr.Trace << "compare " << *this << " to " << other << std::endl;
	if ( antipar_ != other.antipar_ ) {
		tr.Trace << " not the same directionality " << std::endl;
		return false;
	}

	// Make sure both strands overlap (or at least almost overlap)
	const Size MARGIN( basic::options::option[ basic::options::OptionKeys::jumps::max_strand_gap_allowed] - 1 );
	if ( begin1_ > other.end1_ + MARGIN || begin2_ > other.end2_ + MARGIN
			|| other.begin1_ > end1_ + MARGIN || other.begin2_ > end2_ + MARGIN ) {
		tr.Trace << " no overlap between strands " << std::endl;
		tr.Trace << "begin1_ end1 begin2 end2 (repeat with other ) " <<
			begin1_ << " " << end1_ << " " << begin2_ << " " << end2_ << " " <<
			other.begin1_ << " " << other.end1_ << " " << other.begin2_ << " " << other.end2_ << " " << std::endl;
		return false;
	}

	// Make sure the merged strands won't overlap
	if ( end1_ >= other.begin2_ || other.end1_ >= begin2_ ) {
		tr.Trace << "merged strands will overlap " << std::endl;
		return false;
	}
	tr.Trace << "passed presecreen" << std::endl;
	// Make sure starting and ending registers match up
	// (redundant with later, rigorous test, but quick and easy
	// and gets rid of most pairs of unmergeable topologies)
	if ( antipar_ ) {
		if ( (begin1_ + end2_ != other.begin1_ + other.end2_) ||
				(end1_ + begin2_ != other.end1_ + other.begin2_) ) {
			tr.Trace << "register mismatch " << std::endl;
			return false;
		}
	} else {
		if ( (begin1_ + other.begin2_ != other.begin1_  + begin2_) ||
				(end1_ + other.end2_ != other.end1_ + end2_ ) ) {
			tr.Trace << "register mismatch " << std::endl;
			return false;
		}
	}

	runtime_assert( end1_ - begin1_ < pleating1_.size() );
	runtime_assert( end2_ - begin2_ < pairing2_.size() );
	StrandPairing myex(*this);
	myex.extend_to(other.begin1_);
	myex.extend_to(other.end1_);

	StrandPairing otherex(other);
	otherex.extend_to(begin1_);
	otherex.extend_to(end1_);

	if ( myex.end2_ != otherex.end2_ || myex.begin2_ != otherex.begin2_ ) {
		tr.Debug << "SURPRISE!\n" << myex << "\n" << otherex << "\n" << *this << "\n" << other << "\n" << std::endl;
		return false;
	}

	// Make sure pairings, bulges, and pleats match up.
	for ( Size i = 0; i <= myex.end1_ - myex.begin1_; i++ ) {
		runtime_assert( i < myex.pleating1_.size() && i < otherex.pleating1_.size() );
		if ( myex.pleating1_[i] != otherex.pleating1_[i] && myex.pleating1_[i] > 0 && otherex.pleating1_[i] > 0 ) {
			tr.Trace << " wrong pleating " << std::endl;
			return false;
		}
		runtime_assert( i < myex.pairing1_.size() && i < otherex.pairing1_.size() );
		if ( myex.pairing1_[i] != otherex.pairing1_[i] && myex.pairing1_[i] > 0 && otherex.pairing1_[i] > 0 ) {
			tr.Trace << "wrong pairing1_ " << std::endl;
			return false;
		}
	}
	for ( Size i = 0; i <= myex.end2_ - myex.begin2_; i++ ) {
		runtime_assert( i < myex.pairing2_.size() && i < otherex.pairing2_.size() );
		if ( myex.pairing2_[i] != otherex.pairing2_[i] && myex.pairing2_[i] > 0 && otherex.pairing2_[i] > 0 ) {
			tr.Trace << "wrong pairing2_ :" << myex.pairing2_[i] << " " << otherex.pairing2_[i] << std::endl;
			return false;
		}
	}
	tr.Trace << "EQUAL" << std::endl;
	return true;
}

bool StrandPairing::merge(const StrandPairing &other, bool domerge) {
	bool possible ( mergeable( other ) );
	if ( ! domerge || ! possible ) return possible;

	StrandPairing myex(*this);
	myex.extend_to(other.begin1_);
	myex.extend_to(other.end1_);

	StrandPairing otherex(other);
	otherex.extend_to(begin1_);
	otherex.extend_to(end1_);

	//bool changed = other.begin1_ < begin1_ || other.end1_ > end1_;
	// Now do actual merge
	// Add in holes in myex extension that are present in otherex
	for ( Size res = myex.begin1_; res <= myex.end1_; res++ ) {
		Size i = res - myex.begin1_;
		if ( (res < begin1_ || res > end1_) && otherex.pairing1_[i] == 0 ) {
			if ( myex.pairing1_[i] == 0 || otherex.pairing2_[myex.pairing1_[i] - otherex.begin2_] > 0 ) {
				std::cout << "SERIOUS PROBLEM.\n";
			}
			myex.pairing2_[myex.pairing1_[i] - myex.begin2_] = 0;
			myex.pairing1_[i] = myex.pleating1_[i] = 0;
		}
		// Fill in holes in myex overlap with values from otherex
		if ( res >= std::max(begin1_,other.begin1_) && res <= std::min(end1_,other.end1_) && myex.pairing1_[i] == 0 && otherex.pairing1_[i] > 0 ) {
			if ( myex.pairing2_[otherex.pairing1_[i]-myex.begin2_] > 0 ) {
				std::cout << "ANOTHER SERIOUS PROBLEM.\n";
			}
			myex.pairing1_[i] = otherex.pairing1_[i];
			myex.pleating1_[i] = otherex.pleating1_[i];
			myex.pairing2_[otherex.pairing1_[i] - myex.begin2_] = res;
			//changed = true;  // set but never used ~Labonte
		}
	}
	if ( !myex.valid_ends() ) return false;
	*this = myex;

	return possible;
}

Size StrandPairing::get_pleating( Size res ) const {
	if ( res >= begin1_ && res <= end1_ ) {
		return pleating1_[res - begin1_];
	} else if ( res >= begin2_ && res <= end2_ ) {
		if ( pairing2_[res - begin2_] != 0 ) {
			return pleating1_[pairing2_[res-begin2_]];
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}

// Returns the DSSP designation of the given residue (' ' if
// unpaired)
char StrandPairingSet::dssp_state( Size res ) const {
	char state = ' ';
	for ( auto const & pairing : pairings_ ) {
		if ( ! pairing.contains(res) ) continue;

		if ( pairing.is_ladder() ) {
			state = 'E';
		} else if ( state == ' ' ) {
			state = 'B';
		}
	}
	return state;
}

char StrandPairingSet::featurizer_state(Size res) const {
	char state = 'L';
	for ( auto const & pairing : pairings_ ) {
		if ( ! pairing.contains(res) ) continue;

		if ( pairing.is_bulge(res) ) {
			if ( state == 'e' ) {
				state = 'B';
			} else if ( state == 'b' ) {
				state = 'X';
			} else {
				state = 'b';
			}
		} else {
			if ( state == 'e' ) {
				state = 'E';
			} else if ( state == 'b' ) {
				state = 'B';
			} else {
				state = 'e';
			}
		}
	}
	return state;
}

bool StrandPairingSet::has_pairing( core::scoring::dssp::Pairing const& p ) const {
	return paired( p.Pos1(), p.Pos2(), p.is_anti() );
}

bool StrandPairingSet::has_pairing( StrandPairing const& p ) const {
	bool found ( false );
	for ( auto it = pairings_.begin(), end = pairings_.end();
			it != end && !found; ++it ) {
		found = p.mergeable( *it );
	}
	return found;
}

bool StrandPairing::has_pairing( core::scoring::dssp::Pairing const& p ) const {
	return paired( p.Pos1(), p.Pos2(), p.is_anti() );
}

bool StrandPairing::paired( Size res1, Size res2, bool antipar ) const {
	return ( antiparallel() == antipar && contains( res1 ) && get_pair( res1 ) == res2 );
}

bool StrandPairingSet::paired(Size res1, Size res2, bool antiparallel) const {
	for ( auto const & pairing : pairings_ ) {
		if ( pairing.paired( res1, res2, antiparallel ) ) {
			return true;
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////
// Handy function for getting out a list of easy-to-read beta pairings.
//////////////////////////////////////////////////////////////////////////
void StrandPairingSet::get_beta_pairs( core::scoring::dssp::PairingList & beta_pairs ) const {
	for ( auto const & pairing : pairings_ ) {
		pairing.get_beta_pairs( beta_pairs );
	}
}

StrandPairingSet::~StrandPairingSet() = default;

StrandPairing::StrandPairing(Size res1, Size res2, bool antiparallel, Size pleating) :
	begin1_( std::min( res1, res2 ) ), end1_(begin1_),
	begin2_( std::max( res1, res2 ) ), end2_(begin2_),
	antipar_(antiparallel)
{
	pairing1_.push_back( antipar_ ? end2_ : begin2_ );
	pairing2_.push_back( antipar_ ? begin1_ : end1_ );
	pleating1_.push_back(pleating);
}

StrandPairing::StrandPairing() :
	begin1_(0),
	end1_(0),
	begin2_(0),
	end2_(0),
	antipar_(true)
{ }

StrandPairing::~StrandPairing() = default;

Size StrandPairing::get_register() const {
	if ( antipar_ ) {
		return begin1_ + pairing1_[0];
	} else {
		return pairing1_[0] - begin1_;
	}
}

void StrandPairing::get_all_register_and_bulges( SizeList& regs, SizeList& bulges ) const {
	if ( antipar_ ) {
		Size reg = begin1_ + pairing1_[0];
		regs.push_back( reg );
		for ( Size i = 1; ( Size)i < pairing1_.size(); i++ ) {
			if ( pairing1_[i] != 0 && begin1_ + i + pairing1_[i] != reg ) {
				reg = begin1_ + i + pairing1_[i];
				regs.push_back( reg );
				bulges.push_back( begin1_ + i - 1 );
			}
		}
	} else {
		Size reg = pairing1_[0] - begin1_;
		regs.push_back( reg );
		for ( Size i = 1; ( Size)i < pairing1_.size(); i++ ) {
			if ( pairing1_[i] != 0 && pairing1_[i] - begin1_ - i != reg ) {
				reg = pairing1_[i] - begin1_ - i;
				regs.push_back( reg );
				bulges.push_back( begin1_ + i - 1 );
			}
		}
	}
}

bool StrandPairing::range_check() const {
	return (end1_-begin1_ + 1 == pairing1_.size() ) && ( end2_-begin2_ + 1 == pairing2_.size() );
}

std::istream & operator>>( std::istream &is, StrandPairing &sp ) {
	using namespace std;
	char direction;
	is >> direction;
	if ( direction == 'A' ) {
		sp.antipar_ = true;
	} else if ( direction == 'P' ) {
		sp.antipar_ = false;
	} else {
		tr.Trace << "failed reading A/P info --- found instead: " << direction << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}

	char minus;
	string to;
	if ( sp.antipar_ ) {
		is >> sp.begin1_ >> minus >> sp.end2_ >> to >> sp.end1_ >> minus >> sp.begin2_;
		tr.Trace << "read " << sp.begin1_ << "-" << sp.end2_ << " to " << sp.end1_ << "-" << sp.begin2_ << std::endl;
	} else {
		is >> sp.begin1_ >> minus >> sp.begin2_ >> to >> sp.end1_ >> minus >> sp.end2_;
		tr.Trace << "read " << sp.begin1_ << "-" << sp.begin2_ << " to " << sp.end1_ << "-" << sp.end2_ << std::endl;
	}
	if ( minus != '-' ) {
		tr.Trace << "failed reading start pairing --- found instead: " << sp.begin1_ << minus << sp.end2_ << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}

	string regstr;
	is >> regstr;
	if ( regstr != "reg:" ) {
		tr.Trace << "failed reading register tag --- found instead: " << regstr << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}

	runtime_assert( !is.fail() );

	utility::vector1< Size > regs;
	Size reg1;
	string comma  (",");
	while ( comma == "," && is >> reg1 >> comma ) {
		regs.push_back( reg1 );
	}

	if ( comma != "bulges:" && comma != "pleating:" ) {
		tr.Trace << "failed reading bulges tag --- found instead: " << comma << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}

	if ( tr.Trace.visible() ) {
		tr.Trace << " regs: ";
		for ( Size & reg : regs ) {
			tr.Trace << reg << " ";
		}
	}

	utility::vector1< Size > bulges;
	if ( comma == "bulges:" ) {
		comma = ",";
		Size bulge;
		while ( comma == "," && is >> bulge >> comma ) {
			bulges.push_back( bulge );
		}
	}

	if ( tr.Trace.visible() ) {
		tr.Trace << " bulges: ";
		for ( Size & bulge : bulges ) {
			tr.Trace << bulge << " ";
		}
	}

	if ( comma != "pleating:" ) {
		tr.Trace << "failed reading pleating tag --- found instead: " << comma << std::endl;
		is.setstate( std::ios_base::failbit );
		return is;
	}

	// work out how many pairings we have.. then read pleatings
	auto regit = regs.begin(), eregit = regs.end();
	auto bulgeit = bulges.begin(), ebulgeit = bulges.end();
	int dir = sp.antipar_ ? -1 : 1;
	Size bulge = 0;
	for ( Size pos1 = sp.begin1_,
			pos2 = sp.antipar_ ? sp.end2_ : sp.begin2_; pos1 <= sp.end1_; pos1++ ) {
		Size pleat;
		is >> pleat;
		runtime_assert( ( sp.antipar_ ? sp.end2_ - pos2 : pos2 - sp.begin2_ ) == sp.pairing2_.size() );
		Size reg = *regit;
		if ( bulgeit != ebulgeit ) bulge = *bulgeit;
		else bulge = 0;
		tr.Trace << " pos1 " << pos1 << " next bulge " << bulge << " pos2 " << pos2 << std::endl;

		if ( pleat == 0 ) {  //unsatisfied residue is in strand 1
			sp.pairing1_.push_back( 0 );
			sp.pleating1_.push_back( 0 );
			//work out if we have a bulge on strand1, ie. their is no corresponding pairing on strand2
			bool bulge2( false );
			{
				auto next_regit = regit;
				if ( next_regit != eregit ) ++next_regit;
				if ( next_regit != eregit ) {
					Size nex_reg = *next_regit;
					bulge2 = sp.antipar_ ? nex_reg > reg : nex_reg < reg;
				}
			}
			if ( !bulge2 && !bulge ) {
				sp.pairing2_.insert( sp.antipar_ ? sp.pairing2_.begin() : sp.pairing2_.end(), 0 );
				pos2+=dir;
			}
		} else { //
			sp.pairing1_.push_back( pos2 );
			sp.pleating1_.push_back( pleat );
			sp.pairing2_.insert( sp.antipar_ ? sp.pairing2_.begin() : sp.pairing2_.end(), pos1 );
			pos2+=dir;
		}

		if ( pos1 == bulge ) {
			if ( regit == eregit ) {
				tr.Error << "read bulge at " << bulge << " but no new register in list " << std::endl;
				is.setstate( std::ios_base::failbit );
				return is;
			}
			++regit;++bulgeit;

			Size new_pos2 = sp.antipar_ ? ( *regit - pos1 - 1 ) : ( *regit + pos1 + 1);
			tr.Trace << "jump to " << new_pos2 << " in strand2 due to new register " << *regit << " after bulge at " << bulge << std::endl;
			if ( sp.antipar_ && new_pos2 < pos2 ) {
				sp.pairing2_.insert( sp.pairing2_.begin(), pos2 - new_pos2, 0 );
			} else if ( !sp.antipar_ && new_pos2 > pos2 ) {
				sp.pairing2_.insert( sp.pairing2_.end(), new_pos2 - pos2, 0 );
			}
			pos2 = new_pos2;
		}
	} //for ... pos1 = begin1 .. end1

	if ( sp.begin1_ > sp.end1_ || sp.begin2_ > sp.end2_ ) {
		tr.Error  << "begin1 end1 begin2 end2 " << sp.begin1_ << " " << sp.end1_ << " " << sp.begin2_ << " " << sp.end2_ << std::endl;
		utility_exit_with_message( "error reading pairing from stream " );
	}

	tr.Trace << std::endl;
	runtime_assert( !is.fail() );
	sp.show_internals( tr.Trace );
	return is;
}

std::ostream & operator<<(std::ostream & out, const StrandPairing &sp) {
	runtime_assert( sp.begin1_ <= sp.end1_ );
	runtime_assert( sp.begin2_ <= sp.end2_ );
	out << (sp.antipar_ ? 'A' : 'P') << ' ' << sp.begin1_ << '-' << sp.pairing1_[0] << " to " << sp.end1_ << '-' << sp.pairing1_[sp.pairing1_.size()-1] << " reg: ";

	StrandPairing::SizeList regs,bulges;
	sp.get_all_register_and_bulges( regs, bulges );
	bool first( true );
	for ( StrandPairing::SizeList::const_iterator it = regs.begin(), eit = regs.end(); it != eit; ++it ) {
		if ( !first ) {
			out << ", ";
		}
		out << *it;
		first = false;
	}
	first = true;
	if ( bulges.size() ) {
		out << " bulges: ";
		for ( StrandPairing::SizeList::const_iterator it = bulges.begin(), eit = bulges.end(); it != eit; ++it ) {
			if ( !first ) {
				out << ", ";
			}
			out << *it;
			first = false;
		}
		out << " ";
	}
	out << " pleating: ";
	for ( unsigned long i : sp.pleating1_ ) {
		out << i << ' ';
	}
	return out;
}

Size StrandPairing::operator<( StrandPairing const& other) const {
	if ( antipar_ != other.antipar_ ) {
		return antipar_;
	}

	Size reg = antipar_ ? pairing1_[0] + begin1_ : pairing1_[0] - begin1_;
	Size otherreg = antipar_ ? other.pairing1_[0] + other.begin1_ : other.pairing1_[0] - other.begin1_;
	if ( reg == otherreg ) {
		if ( end1_ <= other.begin1_ ) {
			return true;
		} else if ( begin1_ >= other.end1_ ) {
			return false;
		} else {
			tr.Error << "DSSP error: strange strand pairing\n";
			tr.Error << begin1_ << ' ' << end1_ << ' ' << begin2_ << ' ' << end2_ <<' ' << std::endl;
			tr.Error << other.begin1_ << ' ' << other.end1_ << ' ' << other.begin2_ << ' ' << other.end2_ <<' ' << std::endl;
			return begin1_ < other.begin1_;
		}
	} else {
		return reg < otherreg;
	}
}

bool StrandPairing::has_common_pairing( const StrandPairing &other ) const {
	for ( Size i = begin1_; i <= end1_; i++ ) {
		if ( pairing1_[ i - begin1_ ] == 0 )  continue;

		Pairing pair;
		pair.Pos1(i);
		pair.Pos2(pairing1_[ i - begin1_ ]);
		pair.Orientation(antipar_ ? 1 : 2);
		pair.Pleating(pleating1_[ i - begin1_ ]);
		if ( other.has_pairing( pair ) ) return true;
	}
	return false;
}

Size StrandPairing::operator==(const StrandPairing &other) const {
	return (begin1_ == other.begin1_ && begin2_ == other.begin2_ &&
		end1_ == other.end1_ && end2_ == other.end2_ &&
		pairing1_ == other.pairing1_ && pairing2_ == other.pairing2_ &&
		pleating1_ == other.pleating1_);
}

// Return true if the given residue is part of a beta bulge
bool StrandPairing::is_bulge(Size res) const {
	if ( ! contains(res) ) return false;
	if ( get_pair(res) != 0 ) return false;
	return true;
}

// Return true if the given residue is part of this pairing
// (includes bulges)
bool StrandPairing::contains(Size res) const {
	return (res >= begin1_ && res <= end1_) || (res >= begin2_ && res <= end2_);
}

// Return the residue to which the given residue is paired
// (0 if the residue is unpaired, i.e. a bulge residue, or
// is not contained in the pairing at all)
Size StrandPairing::get_pair(Size res) const {
	if ( res >= begin1_ && res <= end1_ ) {
		return pairing1_[res-begin1_];
	} else if ( res >= begin2_ && res <= end2_ ) {
		return pairing2_[res-begin2_];
	} else {
		return 0;
	}
}

// Return true if this pairing is of length greater than 1.
bool StrandPairing::is_ladder() const {
	return end1_-begin1_ > 0 && end2_-begin2_ > 0;
}

bool StrandPairing::antiparallel() const {
	return antipar_;
}

void StrandPairing::get_beta_pairs( core::scoring::dssp::PairingList& beta_pairs ) const {
	for ( Size i = begin1_; i <= end1_; i++ ) {
		if ( pairing1_[ i - begin1_ ] == 0 ) continue;

		core::scoring::dssp::Pairing pair;
		pair.Pos1(i);
		pair.Pos2(pairing1_[ i - begin1_ ]);
		pair.Orientation(antipar_ ? 1 : 2);
		pair.Pleating(pleating1_[ i - begin1_ ]);
		beta_pairs.push_back( pair );
	}
}

bool StrandPairing::valid_ends() const {
	if ( antipar_ ) { ///bad pairing don't write
		Size end2( pairing1_[ 0 ] );
		Size begin2( pairing1_[ pairing1_.size()-1 ] );
		if ( begin2 > end2 ) return false;
	} else {
		Size begin2( pairing1_[ 0 ] );
		Size end2( pairing1_[ pairing1_.size()-1 ] );
		if ( begin2 > end2 ) return false;
	}
	return true;
}

}
}
}
