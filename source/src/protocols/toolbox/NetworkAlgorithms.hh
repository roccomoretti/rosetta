// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/toolbox/NetworkAlgorithms.hh
/// @brief Computes network properties of a pose
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)


#ifndef INCLUDED_protocols_toolbox_NetworkAlgorithms_hh
#define INCLUDED_protocols_toolbox_NetworkAlgorithms_hh

#include <protocols/toolbox/NetworkAlgorithms.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility>
#include <utility/VirtualBase.hh>

//// C++ headers
#include <list>
#include <string>

// Parser headers

namespace protocols {
namespace toolbox {

/// @brief residuenetwork base class
class ResidueNetwork : public utility::VirtualBase
{
public:
	/// @brief default constructor
	ResidueNetwork();

	/// @brief destructor
	~ResidueNetwork() override;

	/// @brief create a network from a pose
	void
	create_from_pose( core::pose::Pose const & pose );

	/// @brief generate a list of edges from the pose -- MUST be reimplemented for each type of network
	virtual void
	generate_edges( core::pose::Pose const & pose ) = 0;

	/// @brief empties edges
	void
	clear_edges();

	/// @brief run Dijkstra's shortest path algorithm on the given list of nodes
	/// after execution, the "distanceFromStart" variable of each node will contain the distance from residue resi
	void
	dijkstras( core::Size const resi ) const;

	/// @brief calculates the connectivity index of residue resi in the context of the network
	core::Real
	connectivity_index( core::Size const resi ) const;

	/// @brief calculates the average shortest path length of the network
	core::Real
	average_shortest_path_length() const;

	// accessors
public:
	// const accessor for nodes
	std::list< NodeOP > const & nodes() const { return nodes_; }

	// other data
private:
	std::list< NodeOP > nodes_;
};

// subclasses

/// @brief Creates networks based on residue proximity in space
class DistanceResidueNetwork : public ResidueNetwork
{
public:
	void
	generate_edges( core::pose::Pose const & pose ) override;
};

/// @brief Creates networks based on covalent connections between residues
class CovalentResidueNetwork : public ResidueNetwork
{
public:
	void
	generate_edges( core::pose::Pose const & pose ) override;
};

// Helper functions
NodeOP
ExtractSmallest( std::list< NodeOP > & nodes);

std::list< NodeOP >
AdjacentRemainingNodes( NodeOP node);

bool
Contains( std::list< NodeOP > const & nodes, NodeCOP node );

// Nodes (and edges) for network algorithms
struct Node : public utility::VirtualBase
{
	Node(
		std::string ID,
		core::Size const res
	) :
		resi( res ),
		id(std::move(ID)),
		distanceFromStart(9999),
		in_list( false )
	{
		neighbors.clear();
	}

	core::Size resi;
	std::string id;
	std::list< NodeOP > neighbors;
	int distanceFromStart;
	bool in_list;
};

} // toolbox
} // protocols

#endif
