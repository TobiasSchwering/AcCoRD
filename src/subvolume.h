/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 * 
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * subvolume.h - 	structure for storing subvolume properties. Simulation
 *					environment is partitioned into subvolumes
 *
 * Last revised for AcCoRD v0.5 (2016-04-15)
 *
 * Revision history:
 *
 * Revision v0.5 (2016-04-15)
 * - corrected memory allocation for subvolume helper arrays
 * - added surface subvolumes
 * - preliminary implementation of interaction between mesoscopic subvolumes and surfaces
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - modified use of unsigned long and unsigned long long to uint32_t and uint64_t
 * - added accommodation of spherical (microscopic) subvolumes and actors
 * - moved counting and IDing of neighbor subvolumes in different regions to a separate
 * function to simplify code reuse
 * - renamed region structure and array of structures to accommodate naming
 * of new boundary structure
 *
 * Revision v0.3.1
 * - header added
*/
#ifndef SUBVOLUME_H
#define SUBVOLUME_H

#include "region.h"
#include "randistrs.h" // For PRNGs
#include "base.h" // For subvolume adjacency
#include "global_param.h" // For region adjacency

//
// Constant definitions
//


//
// Subvolume data type declarations
//

/* The subvolume3D structure lists properties for a single individual subvolume,
* whether it is mesoscopic or microscopic. The entire simulation environment
* is defined by one array of these structures.
*/
struct subvolume3D {
	// ID of subvolume in global mesoscopic subvolume list
	uint32_t mesoID;
	
	// ID of the region containing this subvolume. This ID is important for
	// identifying transitions between regions (both mesoscopic and microscopic)
	unsigned short regionID;
	
	//double location[3]; // Coordinates of center of subvolume [ONLY NEEDED FOR INTERFACE]
	
	// The number of subvolumes that are neighbors to current one.
	// More important for mesoscopic subvolumes
	// To be neighbors, two subvolumes must (at least partially) share a common
	// face.
	unsigned short num_neigh;
	
	// Array of IDs of subvolumes that are neighbors to current subvolume.
	// Length is maximum value of num_neigh for all subvolumes in current region
	uint32_t * neighID;
	
	// Total number of reactions (including diffusion) that molecules in this
	// subvolume can participate in
	//unsigned short num_rxn;
	
	// Is subvolume along boundary of region?
	bool bBoundary;
	//uint32_t boundaryID; // (If applicable) Index of subvolume in boundary list
	
	// Diffusion transition rate from a (mesoscopic) boundary subvolume to all
	// of its neighbours
	// (whether or not each neighbour is in a different region)
	// Only allocated if bBoundary == true AND bMicro == false
	// Size is NUM_MOL_TYPES x num_neigh.
	// Each element gives the diffusion rate to the corresponding
	// neighbour in the neighID array.
	double ** diffRateNeigh;
		
	// Pointers to arrays
	uint64_t * num_mol;
	
	// FUTURE MEMBERS
};

//
// Function Declarations
//

// Initialize Subvolume Structure and Affiliated Arrays
void allocateSubvolArray(const uint32_t numSub,
	struct subvolume3D ** subvolArray);

void allocateSubvolHelper(const uint32_t numSub,
	uint32_t (** subCoorInd)[3],
	uint32_t ***** subID,
	uint32_t (*** subIDSize)[2],
	const short NUM_REGIONS,
	struct region regionArray[]);

// Free memory of subvolume structure
void deleteSubvolArray(const uint32_t numSub,
	struct subvolume3D subvolArray[],
	const unsigned short NUM_MOL_TYPES,
	const short NUM_REGIONS,
	struct region regionArray[]);
	
void deleteSubvolHelper(uint32_t subCoorInd[][3],
	uint32_t **** subID,
	uint32_t (** subIDSize)[2],
	const short NUM_REGIONS,
	struct region regionArray[]);

// Construct the array of structures with details of each subvolume
/* Each structure in the array subvol_spec defines a square/cube region of
* subvolumes with common properties */
void buildSubvolArray(const uint32_t numSub,
	uint32_t * numMesoSub,
	struct subvolume3D subvolArray[],
	const struct spec_region3D subvol_spec[],
	struct region regionArray[],
	const short NUM_REGIONS,
	const unsigned short NUM_MOL_TYPES,
	const unsigned short MAX_RXNS, 
	const double SUBVOL_BASE_SIZE,
	double DIFF_COEF[NUM_REGIONS][NUM_MOL_TYPES],
	uint32_t subCoorInd[numSub][3],
	uint32_t **** subID,
	uint32_t (** subIDSize)[2]);

// Determine whether two subvolumes in neighboring regions are neighbors themselves
// Assert that each subvolume is along its own region's boundary
// If subvolumes are neighbors, then the direction from curID towards neighID
// is written to adjDirection
bool checkSubvolNeigh(struct region regionArray[],
	const short NUM_REGIONS,
	short int curRegion,
	short int neighRegion,
	short int * sphRegion,
	short int * rectRegion,
	uint32_t curID,
	uint32_t curNeighID,
	uint32_t * sphSub,
	uint32_t * rectSub,
	uint32_t numSub,
	uint32_t subCoorInd[numSub][3],
	double boundAdjError,
	unsigned short * adjDirection,
	double curSubBound[6],
	double curNeighBound[6],
	unsigned short * numFaceSph);

// Calculate cartesian coordinates of rectangular subvolume
void findSubvolCoor(double subBound[6],
	const struct region regionSingle,
	uint32_t subCoorInd[3]);
	
#endif // SUBVOLUME_H
