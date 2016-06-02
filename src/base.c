/*
 * The AcCoRD Simulator
 * (Actor-based Communication via Reaction-Diffusion)
 *
 * Copyright 2016 Adam Noel. All rights reserved.
 *
 * For license details, read LICENSE.txt in the root AcCoRD directory
 * For user documentation, read README.txt in the root AcCoRD directory
 *
 * base.c - general utility functions that can apply to different simulation data
 * 			structures
 *
 * Last revised for AcCoRD v0.5 (2016-04-15)
 *
 * Revision history:
 *
 * Revision v0.5 (2016-04-15)
 * - filling in cases for 2D Rectangles
 * - added function to calculate boundary surface area. Renamed boundaryArea
 * function to boundaryVolume to avoid name confusion
 * - added function to return string of boundary name and integrated with error
 * messages
 *
 * Revision v0.4.1
 * - improved use and format of error messages
 *
 * Revision v0.4
 * - filled in cases for spherical boundaries to accommodate spherical regions and actors
 * - added clearance to bBoundaryIntersect to account for cases where we need to know
 * whether two shapes are within a specified distance of each other and not just
 * overlapping
 * - added functions for
 * 		- distance between 2 points
 *		- determining whether one shape is completely inside another shape
 * 		- squaring values
 *		- distance to nearest point on boundary
 * 		- determining intersection of line and boundary
 *
 * Revision v0.3.1.1
 * - added check in boundaryArea to identify invalid rectangle and box definitions. Will
 * 	 now return a zero value if boundary is invalid.
 *
 * Revision v0.3.1
 * - header added
 *
 * Created 2015-02-19
 */

#include "base.h" // for "Public" declarations
#include "float.h" //for DBL_MIN
//
// "Private" Declarations
//

//
// Definitions
//

// Is point inside of boundary?
bool bPointInBoundary(const double point[3], const int boundary1Type,
		const double boundary1[]) {
	switch (boundary1Type) {
	case RECTANGLE:
	case RECTANGULAR_BOX:
		return (point[0] >= boundary1[0] && point[0] <= boundary1[1]
				&& point[1] >= boundary1[2] && point[1] <= boundary1[3]
				&& point[2] >= boundary1[4] && point[2] <= boundary1[5]);
	case SPHERE:
		return (pointDistance(point, boundary1) <= boundary1[3]);
	case CYLINDER: //TODO changed, finish and test
		if (boundary1[4] == PLANE_XY) {
			return (point[2] >= boundary1[2]
					&& point[2] <= (boundary1[2] + boundary1[5])
					&& sqrt(
							squareDBL(point[0] - boundary1[0])
									+ squareDBL(point[1] - boundary1[1]))
							<= boundary1[3]);
		} else if (boundary1[4] == PLANE_XZ) {
			return (point[1] >= boundary1[1]
					&& point[1] <= (boundary1[1] + boundary1[5])
					&& sqrt(
							squareDBL(point[0] - boundary1[0])
									+ squareDBL(point[2] - boundary1[2]))
							<= boundary1[3]);
		} else if (boundary1[4] == PLANE_YZ) {
			return (point[0] >= boundary1[0]
					&& point[0] <= (boundary1[0] + boundary1[5])
					&& sqrt(
							squareDBL(point[1] - boundary1[1])
									+ squareDBL(point[2] - boundary1[2]))
							<= boundary1[3]);
		}
	default:
		fprintf(stderr, "ERROR: Cannot find point in shape type %s.\n",
				boundaryString(boundary1Type));
		return false;
	}
}

// Do two sets of boundaries overlap? TODO add full cylinder support
bool bBoundaryIntersect(const int boundary1Type, const double boundary1[],
		const int boundary2Type, const double boundary2[],
		const double clearance) {
	double d;
	bool lengthcheck = false;
	bool areacheck = false;
	switch (boundary1Type) {
	case RECTANGLE:
	case RECTANGULAR_BOX:
		switch (boundary2Type) {
		case RECTANGULAR_BOX:
			return (boundary1[2] < boundary2[3] && boundary1[3] > boundary2[2]
					&& boundary1[0] < boundary2[1]
					&& boundary1[1] > boundary2[0]
					&& boundary1[4] < boundary2[5]
					&& boundary1[5] > boundary2[4]);
		case SPHERE:
			d = 0;
			if (boundary2[0] < boundary1[0])
				d += squareDBL(boundary2[0] - boundary1[0]);
			else if (boundary2[0] > boundary1[1])
				d += squareDBL(boundary2[0] - boundary1[1]);
			if (boundary2[1] < boundary1[2])
				d += squareDBL(boundary2[1] - boundary1[2]);
			else if (boundary2[1] > boundary1[3])
				d += squareDBL(boundary2[1] - boundary1[3]);
			if (boundary2[2] < boundary1[4])
				d += squareDBL(boundary2[2] - boundary1[4]);
			else if (boundary2[2] > boundary1[5])
				d += squareDBL(boundary2[2] - boundary1[5]);

			return (d < squareDBL(boundary2[3] + clearance)
					&& !bBoundarySurround(RECTANGULAR_BOX, boundary1, SPHERE,
							boundary2, 0.)
					&& !bBoundarySurround(SPHERE, boundary2, RECTANGULAR_BOX,
							boundary1, 0.));
		case CYLINDER:
			; //dummy statement necessary to allow declarations after a label
			bool rectInCircle;
			bool circleInRect;
			int along = 0;
			int across1 = 0;
			int across2 = 0;

			if (boundary2[4] == PLANE_XY) {
				across1 = 0;
				across2 = 1;
				along = 2;
			} else if (boundary2[4] == PLANE_XZ) {
				across1 = 0;
				along = 1;
				across2 = 2;
			} else if (boundary2[4] == PLANE_YZ) {
				along = 0;
				across1 = 1;
				across2 = 2;
			} else {
				fprintf(stderr,
						"ERROR: Cannot determine the orientation of a %s.\n",
						boundaryString(boundary1Type));
				return false;
			}
			lengthcheck = boundary1[2 * along]
					<= boundary2[along] + boundary2[5] - clearance
					&& boundary1[2 * along + 1] >= boundary2[along] + clearance;
			areacheck =
					sqrt( //one of the edges of the rectangle inside the circle
							squareDBL(
									boundary1[2 * across1] - boundary2[across1])
									+ squareDBL(
											boundary1[2 * across2]
													- boundary2[across2]))
							<= boundary2[3] - clearance
							|| sqrt(
									squareDBL(
											boundary1[2 * across1]
													- boundary2[across1])
											+ squareDBL(
													boundary1[2 * across2]
															- boundary2[across2]))
									<= boundary2[3] - clearance
							|| sqrt(
									squareDBL(
											boundary1[2 * across1]
													- boundary2[across1])
											+ squareDBL(
													boundary1[2 * across2]
															- boundary2[across2]))
									<= boundary2[3] - clearance
							|| sqrt(
									squareDBL(
											boundary1[2 * across1]
													- boundary2[across1])
											+ squareDBL(
													boundary1[2 * across2]
															- boundary2[across2]))
									<= boundary2[3] - clearance
							|| (boundary2[across1] >= boundary1[across1 * 2] //or the center of the circle is inside the rectangle
							&& boundary2[across1] <= boundary1[across1 * 2 + 1]
									&& boundary2[across2]
											>= boundary1[across2 * 2]
									&& boundary2[across2]
											<= boundary1[across2 * 2 + 1]);
			return (lengthcheck && areacheck);
		default:
			fprintf(stderr,
					"ERROR: Cannot determine the intersection of a %s and a %s.\n",
					boundaryString(boundary2Type),
					boundaryString(boundary1Type));
			return false;
		}
	case SPHERE:
		switch (boundary2Type) {
		case SPHERE:
			d = pointDistance(boundary1, boundary2);
			return (d < boundary1[3] + boundary2[3] + clearance
					&& d > fabs(boundary1[3] - boundary2[3]));
		case RECTANGLE:
		case RECTANGULAR_BOX:
			d = 0;
			if (boundary1[0] < boundary2[0])
				d += squareDBL(boundary2[0] - boundary1[0]);
			else if (boundary1[0] > boundary2[1])
				d += squareDBL(boundary1[0] - boundary2[1]);
			if (boundary1[1] < boundary2[2])
				d += squareDBL(boundary1[1] - boundary2[2]);
			else if (boundary1[1] > boundary2[3])
				d += squareDBL(boundary1[1] - boundary2[3]);
			if (boundary1[2] < boundary2[4])
				d += squareDBL(boundary1[2] - boundary2[4]);
			else if (boundary1[2] > boundary2[5])
				d += squareDBL(boundary1[2] - boundary2[5]);

			return (d < squareDBL(boundary1[3] + clearance)
					&& !bBoundarySurround(RECTANGULAR_BOX, boundary2, SPHERE,
							boundary1, 0.)
					&& !bBoundarySurround(SPHERE, boundary1, RECTANGULAR_BOX,
							boundary2, 0.));
		default:
			fprintf(stderr,
					"ERROR: Cannot determine the intersection of a %s and a %s.\n",
					boundaryString(boundary2Type),
					boundaryString(boundary1Type));
			return false;
		}
	case CYLINDER:
		switch (boundary2Type) {
		case RECTANGULAR_BOX:
			; //dummy statement necessary to allow declarations after a label
			bool rectInCircle;
			bool circleInRect;
			int along = 0;
			int across1 = 0;
			int across2 = 0;

			if (boundary1[4] == PLANE_XY) {
				across1 = 0;
				across2 = 1;
				along = 2;
			} else if (boundary1[4] == PLANE_XZ) {
				across1 = 0;
				along = 1;
				across2 = 2;
			} else if (boundary1[4] == PLANE_YZ) {
				along = 0;
				across1 = 1;
				across2 = 2;
			} else {
				fprintf(stderr,
						"ERROR: Cannot determine the orientation of a %s.\n",
						boundaryString(boundary1Type));
				return false;
			}
			lengthcheck = boundary2[2 * along]
					<= boundary1[along] + boundary1[5] - clearance
					&& boundary2[2 * along + 1] >= boundary1[along] + clearance;
			areacheck =
					sqrt( //one of the edges of the rectangle inside the circle
							squareDBL(
									boundary2[2 * across1] - boundary1[across1])
									+ squareDBL(
											boundary2[2 * across2]
													- boundary1[across2]))
							<= boundary1[3] - clearance
							|| sqrt(
									squareDBL(
											boundary2[2 * across1]
													- boundary1[across1])
											+ squareDBL(
													boundary2[2 * across2]
															- boundary1[across2]))
									<= boundary1[3] - clearance
							|| sqrt(
									squareDBL(
											boundary2[2 * across1]
													- boundary1[across1])
											+ squareDBL(
													boundary2[2 * across2]
															- boundary1[across2]))
									<= boundary1[3] - clearance
							|| sqrt(
									squareDBL(
											boundary2[2 * across1]
													- boundary1[across1])
											+ squareDBL(
													boundary2[2 * across2]
															- boundary1[across2]))
									<= boundary1[3] - clearance
							|| (boundary1[across1] >= boundary2[across1 * 2] //or the center of the circle is inside the rectangle
							&& boundary1[across1] <= boundary2[across1 * 2 + 1]
									&& boundary1[across2]
											>= boundary2[across2 * 2]
									&& boundary1[across2]
											<= boundary2[across2 * 2 + 1]);
			return (lengthcheck && areacheck);
		default:
			fprintf(stderr,
					"ERROR: Cannot determine the intersection of a %s and a %s.\n",
					boundaryString(boundary2Type),
					boundaryString(boundary1Type));
			return false;
		}
	default:
		fprintf(stderr, "ERROR: Cannot find intersection with shape %s.\n",
				boundaryString(boundary1Type));
		return false;
	}
}

// Are two sets of boundaries adjacent? Intersections will not be detected.
// Both boundaries must be rectangular (either 2D or 3D)
bool bBoundaryAdjacent(const int boundary1Type, const double boundary1[],
		const int boundary2Type, const double boundary2[],
		const double distError, unsigned short * direction) {

	if ((boundary1Type == RECTANGULAR_BOX && boundary2Type == RECTANGULAR_BOX)
			|| (boundary1Type == RECTANGLE && boundary2Type == RECTANGULAR_BOX)
			|| (boundary1Type == RECTANGULAR_BOX && boundary2Type == RECTANGLE)) {
		if ( // Do boxes share face along xy-plane?
		(boundary1[1] > boundary2[0] + distError)
				&& (boundary2[1] > boundary1[0] + distError)
				&& (boundary1[3] > boundary2[2] + distError)
				&& (boundary2[3] > boundary1[2] + distError)) {
			if (fabs(boundary1[4] - boundary2[5]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower z
					{
				*direction = IN;
				return true;
			} else if (fabs(boundary2[4] - boundary1[5]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper z
					{
				*direction = OUT;
				return true;
			}
		} else if ( // Do boxes share face along zy-plane?
		(boundary1[3] > boundary2[2] + distError)
				&& (boundary2[3] > boundary1[2] + distError)
				&& (boundary1[5] > boundary2[4] + distError)
				&& (boundary2[5] > boundary1[4] + distError)) {
			if (fabs(boundary1[0] - boundary2[1]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower x
					{
				*direction = LEFT;
				return true;
			} else if (fabs(boundary2[0] - boundary1[1]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper x
					{
				*direction = RIGHT;
				return true;
			}
		} else if ( // Do boxes share face along zx-plane?
		(boundary1[1] > boundary2[0] + distError)
				&& (boundary2[1] > boundary1[0] + distError)
				&& (boundary1[5] > boundary2[4] + distError)
				&& (boundary2[5] > boundary1[4] + distError)) {
			if (fabs(boundary1[2] - boundary2[3]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower y
					{
				*direction = DOWN;
				return true;
			} else if (fabs(boundary2[2] - boundary1[3]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper y
					{
				*direction = UP;
				return true;
			}
		}
	} else if (boundary1Type == RECTANGLE && boundary2Type == RECTANGLE) // Boundaries are both rectangles. They must lie in same plane to have adjacency
	{
		if (boundary1[0] == boundary1[1]
				&& fabs(boundary1[0] - boundary2[0]) < distError
				&& fabs(boundary1[0] - boundary2[1]) < distError) // boundaries are both in YZ plane
						{
			if ((boundary1[3] > boundary2[2] + distError)
					&& (boundary2[3] > boundary1[2] + distError)) // There is overlap along Y
					{
				if (fabs(boundary1[4] - boundary2[5]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower z
						{
					*direction = IN;
					return true;
				} else if (fabs(boundary2[4] - boundary1[5]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper z
						{
					*direction = OUT;
					return true;
				}
			} else if ((boundary1[5] > boundary2[4] + distError)
					&& (boundary2[5] > boundary1[4] + distError)) // There is overlap along Z
					{
				if (fabs(boundary1[2] - boundary2[3]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower y
						{
					*direction = DOWN;
					return true;
				} else if (fabs(boundary2[2] - boundary1[3]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper y
						{
					*direction = UP;
					return true;
				}
			}
		} else if (boundary1[2] == boundary1[3]
				&& fabs(boundary1[2] - boundary2[2]) < distError
				&& fabs(boundary1[2] - boundary2[3]) < distError) // boundaries are both in XZ plane
						{
			if ((boundary1[1] > boundary2[0] + distError)
					&& (boundary2[1] > boundary1[0] + distError)) // There is overlap along X
					{
				if (fabs(boundary1[4] - boundary2[5]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower z
						{
					*direction = IN;
					return true;
				} else if (fabs(boundary2[4] - boundary1[5]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper z
						{
					*direction = OUT;
					return true;
				}
			} else if ((boundary1[5] > boundary2[4] + distError)
					&& (boundary2[5] > boundary1[4] + distError)) // There is overlap along Z
					{
				if (fabs(boundary1[0] - boundary2[1]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower x
						{
					*direction = LEFT;
					return true;
				} else if (fabs(boundary2[0] - boundary1[1]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper x
						{
					*direction = RIGHT;
					return true;
				}
			}
		} else if (boundary1[4] == boundary1[5]
				&& fabs(boundary1[4] - boundary2[4]) < distError
				&& fabs(boundary1[4] - boundary2[5]) < distError) // boundaries are both in XY plane
						{
			if ((boundary1[1] > boundary2[0] + distError)
					&& (boundary2[1] > boundary1[0] + distError)) // There is overlap along X
					{
				if (fabs(boundary1[2] - boundary2[3]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower y
						{
					*direction = DOWN;
					return true;
				} else if (fabs(boundary2[2] - boundary1[3]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper y
						{
					*direction = UP;
					return true;
				}
			} else if ((boundary1[3] > boundary2[2] + distError)
					&& (boundary2[3] > boundary1[2] + distError)) // There is overlap along Y
					{
				if (fabs(boundary1[0] - boundary2[1]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's lower x
						{
					*direction = LEFT;
					return true;
				} else if (fabs(boundary2[0] - boundary1[1]) < distError) // Boundary 2 is adjacent to boundary 1 along 1's upper x
						{
					*direction = RIGHT;
					return true;
				}
			}
		}
		//TODO changed, test and validate. The orientation sorting could be used and *direction calculated from along, but that's kinda hacky
	} else if (boundary1Type == CYLINDER && boundary2Type == CYLINDER) {
		if (boundary1[4] == boundary2[4]) {
			if (boundary1[4] == PLANE_XY) {
				if (sqrt(
						squareDBL(boundary1[0] - boundary2[0])
								+ squareDBL(boundary1[1] - boundary2[1]))
						< boundary1[3] + boundary2[3] + distError) { //radial overlap

					if (boundary1[2] > boundary2[2] + boundary2[5] - distError
							&& boundary1[2]
									< boundary2[2] + boundary2[5] + distError) {
						*direction = IN; // Boundary 2 is adjacent to boundary 1 along 1's lower z
						return true;
					} else if (boundary2[2]
							> boundary1[2] + boundary1[5] - distError
							&& boundary2[2]
									< boundary1[2] + boundary1[5] + distError) {
						*direction = OUT; // Boundary 2 is adjacent to boundary 1 along 1's upper z
						return true;
					}
				}
			} else if (boundary1[4] == PLANE_XZ) {
				if (sqrt(
						squareDBL(boundary1[0] - boundary2[0])
								+ squareDBL(boundary1[2] - boundary2[2]))
						< boundary1[3] + boundary2[3] + distError) { //radial overlap
					if (boundary1[1] > boundary2[1] + boundary2[5] - distError
							&& boundary1[1]
									< boundary2[1] + boundary2[5] + distError) {
						*direction = DOWN; // Boundary 2 is adjacent to boundary 1 along 1's lower y
						return true;
					} else if (boundary2[1]
							> boundary1[1] + boundary1[5] - distError
							&& boundary2[1]
									< boundary1[1] + boundary1[5] + distError) {
						*direction = UP; // Boundary 2 is adjacent to boundary 1 along 1's upper y
						return true;
					}
				}
			} else if (boundary1[4] == PLANE_YZ) {
				if (sqrt(
						squareDBL(boundary1[1] - boundary2[1])
								+ squareDBL(boundary1[2] - boundary2[2]))
						< boundary1[3] + boundary2[3] + distError) { //radial overlap

					if (boundary1[0] > boundary2[0] + boundary2[5] - distError
							&& boundary1[0]
									< boundary2[0] + boundary2[5] + distError) {
						*direction = LEFT; // Boundary 2 is adjacent to boundary 1 along 1's lower x
						return true;
					} else if (boundary2[0]
							> boundary1[0] + boundary1[5] - distError
							&& boundary2[0]
									< boundary2[0] + boundary1[5] + distError) {
						*direction = RIGHT; // Boundary 2 is adjacent to boundary 1 along 1's lower x
						return true;
					}
				} else {
					fprintf(stderr,
							"ERROR: Cannot determine the orientation of a %s.\n",
							boundaryString(boundary1Type));
					return false;
				}
			}
		} else {
			fprintf(stderr,
					"ERROR: Cannot determine whether 2 Cylinders are adjacent if they have different orientations.\n");
			return false;

		}
	} else {
		fprintf(stderr,
				"ERROR: Cannot determine whether a %s and a %s are adjacent.\n",
				boundaryString(boundary2Type), boundaryString(boundary1Type));
	}
	return false;
}

// Is first boundary entirely inside the second?
bool bBoundarySurround(const int boundary1Type, const double boundary1[],
		const int boundary2Type, const double boundary2[],
		const double clearance) {
	double p1[3];
	double p2[3];
	bool lengthcheck;
	bool areacheck;

	switch (boundary1Type) // Is boundary1 inside of boundary2?
	{
	case RECTANGLE:
	case RECTANGULAR_BOX:
		switch (boundary2Type) {
		case RECTANGLE:
		case RECTANGULAR_BOX:
			return (boundary1[0] >= boundary2[0] + clearance
					&& boundary1[1] <= boundary2[1] - clearance
					&& boundary1[2] >= boundary2[2] + clearance
					&& boundary1[3] <= boundary2[3] - clearance
					&& boundary1[4] >= boundary2[4] + clearance
					&& boundary1[5] <= boundary2[5] - clearance);
		case SPHERE:
			p1[0] = boundary1[0];
			p1[1] = boundary1[2];
			p1[2] = boundary1[4];
			if (boundary2[3] < pointDistance(p1, boundary2) + clearance)
				return false;
			p1[0] = boundary1[0];
			p1[1] = boundary1[2];
			p1[2] = boundary1[5];
			if (boundary2[3] < pointDistance(p1, boundary2) + clearance)
				return false;
			p1[0] = boundary1[0];
			p1[1] = boundary1[3];
			p1[2] = boundary1[4];
			if (boundary2[3] < pointDistance(p1, boundary2) + clearance)
				return false;
			p1[0] = boundary1[0];
			p1[1] = boundary1[3];
			p1[2] = boundary1[5];
			if (boundary2[3] < pointDistance(p1, boundary2) + clearance)
				return false;
			p1[0] = boundary1[1];
			p1[1] = boundary1[2];
			p1[2] = boundary1[4];
			if (boundary2[3] < pointDistance(p1, boundary2) + clearance)
				return false;
			p1[0] = boundary1[1];
			p1[1] = boundary1[2];
			p1[2] = boundary1[5];
			if (boundary2[3] < pointDistance(p1, boundary2) + clearance)
				return false;
			p1[0] = boundary1[1];
			p1[1] = boundary1[3];
			p1[2] = boundary1[4];
			if (boundary2[3] < pointDistance(p1, boundary2) + clearance)
				return false;
			p1[0] = boundary1[1];
			p1[1] = boundary1[3];
			p1[2] = boundary1[5];
			if (boundary2[3] < pointDistance(p1, boundary2) + clearance)
				return false;
			// All fail cases have been tried
			return true;
		case CYLINDER: //TODO changed, finish and validate
			if (boundary2[4] == PLANE_XY) {
				lengthcheck = boundary1[4] >= boundary2[2] + clearance
						&& boundary1[5]
								<= boundary2[2] + boundary2[5] - clearance;
				//check whether the distance of all 4 edges to the center of the circle is less than its radius
				areacheck = sqrt(
						squareDBL(boundary1[0] - boundary2[0])
								+ squareDBL(boundary1[2] - boundary2[1]))
						<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[1] - boundary2[0])
										+ squareDBL(
												boundary1[2] - boundary2[1]))
								<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[0] - boundary2[0])
										+ squareDBL(
												boundary1[3] - boundary2[1]))
								<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[1] - boundary2[0])
										+ squareDBL(
												boundary1[3] - boundary2[1]))
								<= boundary2[3] - clearance;
			} else if (boundary2[4] == PLANE_XZ) {
				lengthcheck = boundary1[2] >= boundary2[1] + clearance
						&& boundary1[3]
								<= boundary2[1] + boundary2[5] - clearance;
				areacheck = sqrt(
						squareDBL(boundary1[0] - boundary2[0])
								+ squareDBL(boundary1[4] - boundary2[2]))
						<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[1] - boundary2[0])
										+ squareDBL(
												boundary1[4] - boundary2[2]))
								<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[0] - boundary2[0])
										+ squareDBL(
												boundary1[5] - boundary2[2]))
								<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[1] - boundary2[0])
										+ squareDBL(
												boundary1[5] - boundary2[2]))
								<= boundary2[3] - clearance;
			} else if (boundary2[4] == PLANE_YZ) {
				lengthcheck = boundary1[0] >= boundary2[0] + clearance
						&& boundary1[1]
								<= boundary2[0] + boundary2[5] - clearance;
				areacheck = sqrt(
						squareDBL(boundary1[2] - boundary2[1])
								+ squareDBL(boundary1[4] - boundary2[2]))
						<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[3] - boundary2[1])
										+ squareDBL(
												boundary1[4] - boundary2[2]))
								<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[2] - boundary2[1])
										+ squareDBL(
												boundary1[5] - boundary2[2]))
								<= boundary2[3] - clearance
						&& sqrt(
								squareDBL(boundary1[3] - boundary2[1])
										+ squareDBL(
												boundary1[5] - boundary2[2]))
								<= boundary2[3] - clearance;
			}
			return areacheck && lengthcheck;
		default:
			fprintf(stderr,
					"ERROR: Cannot determine whether a %s is inside of a %s.\n",
					boundaryString(boundary2Type),
					boundaryString(boundary1Type));
			return false;
		}
	case SPHERE:
		switch (boundary2Type) {
		case RECTANGLE:
			return false; // A 3D object cannot be inside of a 2D object
		case RECTANGULAR_BOX:
			return (boundary1[3] <= (boundary1[0] - boundary2[0] - clearance)
					&& boundary1[3] <= (boundary2[1] - boundary1[0] - clearance)
					&& boundary1[3] <= (boundary1[1] - boundary2[2] - clearance)
					&& boundary1[3] <= (boundary2[3] - boundary1[1] - clearance)
					&& boundary1[3] <= (boundary1[2] - boundary2[4] - clearance)
					&& boundary1[3] <= (boundary2[5] - boundary1[2] - clearance));
		case SPHERE:
			return (boundary2[3]
					>= (boundary1[3] + pointDistance(boundary1, boundary2)
							+ clearance));
		default:
			fprintf(stderr,
					"ERROR: Cannot determine whether a %s is inside of a %s.\n",
					boundaryString(boundary2Type),
					boundaryString(boundary1Type));
			return false;
		}
	case CYLINDER:
		switch (boundary2Type) {
		case RECTANGULAR_BOX:
			if (boundary1[4] == PLANE_XY) {
				lengthcheck = boundary2[4] <= boundary1[2] - clearance
						&& boundary2[5]
								>= boundary1[2] + boundary1[5] + clearance;
				areacheck = boundary2[0]
						<= boundary1[0] - boundary1[3] - clearance
						&& boundary2[1]
								>= boundary1[0] + boundary1[3] + clearance
						&& boundary2[2]
								<= boundary1[1] - boundary1[3] - clearance
						&& boundary2[3]
								>= boundary1[1] + boundary1[3] + clearance;
			} else if (boundary1[4] == PLANE_XZ) {
				lengthcheck = boundary2[2] <= boundary1[1] - clearance
						&& boundary2[3]
								>= boundary1[1] + boundary1[5] + clearance;
				areacheck = boundary2[0]
						<= boundary1[0] - boundary1[3] - clearance
						&& boundary2[1]
								>= boundary1[0] + boundary1[3] + clearance
						&& boundary2[4]
								<= boundary1[2] - boundary1[3] - clearance
						&& boundary2[5]
								>= boundary1[2] + boundary1[3] + clearance;
			} else if (boundary1[4] == PLANE_YZ) {
				lengthcheck = boundary2[0] <= boundary1[0] - clearance
						&& boundary2[1]
								>= boundary1[0] + boundary1[5] + clearance;
				areacheck = boundary2[2]
						<= boundary1[1] - boundary1[3] - clearance
						&& boundary2[3]
								>= boundary1[1] + boundary1[3] + clearance
						&& boundary2[4]
								<= boundary1[2] - boundary1[3] - clearance
						&& boundary2[5]
								>= boundary1[2] + boundary1[3] + clearance;
			}
			return areacheck && lengthcheck;
		case CYLINDER: //TODO test and verify, perhaps add test with different orientations
			//TODO spread this definition to other calculations! far more typo resistant
			if (boundary1[4] == boundary2[4]) {
				int along = 0;
				int across1 = 0;
				int across2 = 0;
				if (boundary1[4] == PLANE_XY) {
					across1 = 0;
					across2 = 1;
					along = 2;
				} else if (boundary1[4] == PLANE_XZ) {
					across1 = 0;
					along = 1;
					across2 = 2;
				} else if (boundary1[4] == PLANE_YZ) {
					along = 0;
					across1 = 1;
					across2 = 2;
				} else {
					fprintf(stderr,
							"ERROR: Cannot determine the orientation of a %s.\n",
							boundaryString(boundary1Type));
					return false;
				}
				lengthcheck = boundary1[along] >= boundary2[along] + clearance
						&& boundary1[along] + boundary1[5]
								<= boundary2[along] + boundary2[5] - clearance;
				areacheck = sqrt(
						squareDBL(boundary1[across1] - boundary2[across1])
								+ squareDBL(
										boundary1[across2]
												- boundary2[across2]))
						<= boundary2[3] - boundary1[3] - clearance;
				return lengthcheck && areacheck;
			} else {

				fprintf(stderr,
						"ERROR: Cannot determine whether a %s is inside of a %s of a different orientation.\n",
						boundaryString(boundary2Type),
						boundaryString(boundary1Type));
				return false;
			}
		default:
			fprintf(stderr,
					"ERROR: Cannot determine whether a %s is inside of a %s.\n",
					boundaryString(boundary2Type),
					boundaryString(boundary1Type));
			return false;
		}
	default:
		fprintf(stderr,
				"ERROR: Cannot determine whether shape %s is inside another boundary.\n",
				boundaryString(boundary1Type));
		return false;
	}
}

// Does a point lie within box created by two other points?
bool bPointBetween(const double p1[3], const double p2[3],
		const double newPoint[3]) {
	int i;

	for (i = 0; i < 3; i++) {
		if (p1[i] > p2[i]) {
			if (newPoint[i] < p2[i] || newPoint[i] > p1[i])
				return false;
		} else {
			if (newPoint[i] > p2[i] || newPoint[i] < p1[i])
				return false;
		}
	}
	return true;
}

// Does a line segment intersect some boundary face? If so then which one and where?
// Returns the closest intersecting face from point p1 in positive direction along
// unit vector L
bool bLineHitBoundary(const double p1[3], const double L[3],
		const double length, const int boundary1Type, const double boundary1[],
		short * planeID, const short planeIDConst, const bool bInside,
		double * d, double intersectPoint[3]) {
	short curPlane;
	double minDist = INFINITY;
	double nearestIntersectPoint[3];
	bool bIntersect = false;

	switch (boundary1Type) {
	case RECTANGLE:
		if (bLineHitInfinitePlane(p1, L, length, RECTANGLE, boundary1,
				planeIDConst, false, d, intersectPoint)
				&& bPointOnFace(intersectPoint, RECTANGLE, boundary1,
						planeIDConst) && *d < minDist) {
			return true;
		}
		return false;
	case RECTANGULAR_BOX:
		for (curPlane = 0; curPlane < 6; curPlane++) {
			if (bLineHitInfinitePlane(p1, L, length, RECTANGULAR_BOX, boundary1,
					curPlane, false, d, intersectPoint)
					&& bPointOnFace(intersectPoint, RECTANGULAR_BOX, boundary1,
							curPlane) && *d < minDist) // Line does intersect this face at a valid distance and it is closest
							{
				bIntersect = true;
				*planeID = curPlane;
				minDist = *d;
				nearestIntersectPoint[0] = intersectPoint[0];
				nearestIntersectPoint[1] = intersectPoint[1];
				nearestIntersectPoint[2] = intersectPoint[2];
			}
		}
		if (bIntersect) {
			*d = minDist;
			intersectPoint[0] = nearestIntersectPoint[0];
			intersectPoint[1] = nearestIntersectPoint[1];
			intersectPoint[2] = nearestIntersectPoint[2];
			return true;
		}
		return false;

		//TODO: changed, test!
		//a check whether the intersection is on the boundary should not be necessary
		//as the nearest intersection is used
	case CYLINDER:
		; //dummy statement necessary to allow declarations after a label
		double centerToP1[3];
		double LDotCenterToP1;
		int along = 0;
		int across1 = 0;
		int across2 = 0;

		if (boundary1[4] == PLANE_XY) {
			across1 = 0;
			across2 = 1;
			along = 2;
		} else if (boundary1[4] == PLANE_XZ) {
			across1 = 0;
			along = 1;
			across2 = 2;
		} else if (boundary1[4] == PLANE_YZ) {
			along = 0;
			across1 = 1;
			across2 = 2;
		} else {
			fprintf(stderr,
					"ERROR: Cannot determine the orientation of a %s.\n",
					boundaryString(boundary1Type));
			return false;
		}

		//test against all planes
		for (curPlane = 0; curPlane < 6; curPlane++) {

			if (curPlane == along * 2) //lower circular face
				*d = (boundary1[along] - p1[along]) / L[along];
			else if (curPlane == along * 2 + 1) //upper circular face
				*d = (boundary1[along] + boundary1[5] - p1[along]) / L[along];
			else if (boundary1[5] > 0.) { //mantle face, length necessary
				//TODO update necessary, can still cause problems with tests against the outside
				centerToP1[across1] = p1[across1] - boundary1[across1];
				centerToP1[across2] = p1[across2] - boundary1[across2];
				LDotCenterToP1 = L[across1] * centerToP1[across1]
						+ L[across2] * centerToP1[across2];
				*d = sqrt(
						squareDBL(LDotCenterToP1) + squareDBL(boundary1[3])
								- squareDBL(centerToP1[across1])
								- squareDBL(centerToP1[across2]));
				if (bInside)
					*d = -LDotCenterToP1 + *d;
				else
					*d = -LDotCenterToP1 - *d;
			} else
				continue;

			intersectPoint[0] = p1[0] + L[0] * (*d);
			intersectPoint[1] = p1[1] + L[1] * (*d);
			intersectPoint[2] = p1[2] + L[2] * (*d);

			if (*d > 0. && *d <= length && *d < minDist) {
				bIntersect = true;
				*planeID = curPlane; //TODO does this still work?
				minDist = *d;
				nearestIntersectPoint[0] = intersectPoint[0];
				nearestIntersectPoint[1] = intersectPoint[1];
				nearestIntersectPoint[2] = intersectPoint[2];
			}
		}
		if (bIntersect) {
			//TODO: teststuff, improve or remove!
			if(boundary1[5] == 0. && sqrt(squareDBL(nearestIntersectPoint[across1] - boundary1[across1]) + squareDBL(nearestIntersectPoint[across2] - boundary1[across2])) > boundary1[3])
				return false;

			*d = minDist;
			intersectPoint[0] = nearestIntersectPoint[0];
			intersectPoint[1] = nearestIntersectPoint[1];
			intersectPoint[2] = nearestIntersectPoint[2];
			//TODO: added, test and validate
			//it can happen due to different calculation methods that the intersect point is not regarded as
			//inside the boundary by bPointInBoundary, which can lead to disastrous errors. so in this case
			//the intersectPoint is slightly pushed into the cylinders cross section
//			if (!bPointInBoundary(intersectPoint, boundary1Type, boundary1)) {
//				double pushFrac = 1e-6;
//				double vToCenter[3];
//
//				vToCenter[along] = 0.;
//				vToCenter[across1] = boundary1[across1]
//						- intersectPoint[across1];
//				vToCenter[across2] = boundary1[across2]
//						- intersectPoint[across2];
//
//				while (!bPointInBoundary(intersectPoint, boundary1Type,
//						boundary1) && pushFrac < 0.0001) {
//					pushFrac *= 2;
//					pushPoint(nearestIntersectPoint, intersectPoint, pushFrac,
//							vToCenter);
//				}
//				if (pushFrac >= 0.0001)
//					fprintf(stderr,
//							"ERROR: Failed to resolve inconsistency inside a %s.\n",
//							boundaryString(boundary1Type));
//			}
//			return bPointInBoundary(intersectPoint, boundary1Type, boundary1);
			return true;
		}
		return false;
	case SPHERE:
		return bLineHitInfinitePlane(p1, L, length, SPHERE, boundary1, curPlane,
				bInside, d, intersectPoint);
	default:
		fprintf(stderr,
				"ERROR: Cannot determine whether shape %s intersects another shape.\n",
				boundaryString(boundary1Type));
		return false;
	}
}

// Does a line segment hit an infinite plane? If so then where?
bool bLineHitInfinitePlane(const double p1[3], const double L[3],
		const double length, const int boundary1Type, const double boundary1[],
		const short planeID, const bool bInside, double * d,
		double intersectPoint[3]) {
	double centerToP1[3];
	double LDotCenterToP1;

	switch (boundary1Type) {
	case RECTANGLE:
		switch (planeID) {
		case PLANE_XY:
			*d = (boundary1[4] - p1[2]) / L[2];
			break;
		case PLANE_XZ:
			*d = (boundary1[2] - p1[1]) / L[1];
			break;
		case PLANE_YZ:
			*d = (boundary1[0] - p1[0]) / L[0];
			break;
		}
		intersectPoint[0] = (*d) * L[0] + p1[0];
		intersectPoint[1] = (*d) * L[1] + p1[1];
		intersectPoint[2] = (*d) * L[2] + p1[2];
		break;
	case RECTANGULAR_BOX:
		switch (planeID) {
		case 0:
			*d = (boundary1[0] - p1[0]) / L[0];
			break;
		case 1:
			*d = (boundary1[1] - p1[0]) / L[0];
			break;
		case 2:
			*d = (boundary1[2] - p1[1]) / L[1];
			break;
		case 3:
			*d = (boundary1[3] - p1[1]) / L[1];
			break;
		case 4:
			*d = (boundary1[4] - p1[2]) / L[2];
			break;
		case 5:
			*d = (boundary1[5] - p1[2]) / L[2];
			break;
		}
		intersectPoint[0] = (*d) * L[0] + p1[0];
		intersectPoint[1] = (*d) * L[1] + p1[1];
		intersectPoint[2] = (*d) * L[2] + p1[2];
		break;
	case SPHERE:

		centerToP1[0] = p1[0] - boundary1[0];
		centerToP1[1] = p1[1] - boundary1[1];
		centerToP1[2] = p1[2] - boundary1[2];

		LDotCenterToP1 = L[0] * centerToP1[0] + L[1] * centerToP1[1]
				+ L[2] * centerToP1[2];

		*d = sqrt(
				squareDBL(LDotCenterToP1) + squareDBL(boundary1[3])
						- squareDBL(centerToP1[0]) - squareDBL(centerToP1[1])
						- squareDBL(centerToP1[2]));

		if (bInside)
			*d = -LDotCenterToP1 + *d;
		else
			*d = -LDotCenterToP1 - *d;

		intersectPoint[0] = p1[0] + L[0] * (*d);
		intersectPoint[1] = p1[1] + L[1] * (*d);
		intersectPoint[2] = p1[2] + L[2] * (*d);
		break;
	default:
		fprintf(stderr,
				"ERROR: Cannot determine whether a line hits the plane of a %s.\n",
				boundaryString(boundary1Type));
		*d = 0;
		return false;
	}

	return *d > 0. && *d <= length;
}

// Is point that is in infinite plane also on boundary face?
// Assert that point is already on corresponding plane
bool bPointOnFace(const double p1[3], const int boundary1Type,
		const double boundary1[], const short planeID) {
	switch (boundary1Type) {
	case RECTANGLE:
		switch (planeID) {
		case PLANE_XY:
			return (p1[1] >= boundary1[2] && p1[1] <= boundary1[3]
					&& p1[0] >= boundary1[0] && p1[0] <= boundary1[1]);
		case PLANE_XZ:
			return (p1[0] >= boundary1[0] && p1[0] <= boundary1[1]
					&& p1[2] >= boundary1[4] && p1[2] <= boundary1[5]);
		case PLANE_YZ:
			return (p1[1] >= boundary1[2] && p1[1] <= boundary1[3]
					&& p1[2] >= boundary1[4] && p1[2] <= boundary1[5]);
		}
	case RECTANGULAR_BOX:
		switch (planeID) {
		case 0: // yz plane
		case 1:
			return (p1[1] >= boundary1[2] && p1[1] <= boundary1[3]
					&& p1[2] >= boundary1[4] && p1[2] <= boundary1[5]);
		case 2: // xz plane
		case 3:
			return (p1[0] >= boundary1[0] && p1[0] <= boundary1[1]
					&& p1[2] >= boundary1[4] && p1[2] <= boundary1[5]);
		case 4: // xy plane
		case 5:
			return (p1[1] >= boundary1[2] && p1[1] <= boundary1[3]
					&& p1[0] >= boundary1[0] && p1[0] <= boundary1[1]);
		}
	case SPHERE:
		// Trivially true
		return true;
	case CYLINDER: //TODO verify
		if (boundary1[4] == PLANE_XY) {
			if (planeID == 4 || planeID == 5)
				return (sqrt(
						squareDBL(p1[0] - boundary1[0])
								+ squareDBL(p1[1] - boundary1[1]))
						<= boundary1[5]
						&& (p1[2] == boundary1[2]
								|| p1[2] == boundary1[2] + boundary1[5]));
			else
				return (sqrt(
						squareDBL(p1[0] - boundary1[0])
								+ squareDBL(p1[1] - boundary1[1]))
						== boundary1[5]
						&& (boundary1[2] <= p1[2]
								&& p1[2] <= boundary1[2] + boundary1[5]));
		} else if (boundary1[4] == PLANE_XZ) {
			if (planeID == 2 || planeID == 3)
				return (sqrt(
						squareDBL(p1[0] - boundary1[0])
								+ squareDBL(p1[2] - boundary1[2]))
						<= boundary1[5]
						&& (p1[1] == boundary1[1]
								|| p1[1] == boundary1[1] + boundary1[5]));
			else
				return (sqrt(
						squareDBL(p1[0] - boundary1[0])
								+ squareDBL(p1[2] - boundary1[2]))
						== boundary1[5]
						&& (boundary1[1] <= p1[1]
								&& p1[2] <= boundary1[1] + boundary1[5]));
		} else if (boundary1[4] == PLANE_YZ) {
			if (planeID == 0 || planeID == 1)
				return (sqrt(
						squareDBL(p1[1] - boundary1[1])
								+ squareDBL(p1[2] - boundary1[2]))
						<= boundary1[5]
						&& (p1[0] == boundary1[0]
								|| p1[0] == boundary1[0] + boundary1[5]));
			else
				return (sqrt(
						squareDBL(p1[1] - boundary1[1])
								+ squareDBL(p1[2] - boundary1[2]))
						== boundary1[5]
						&& (boundary1[0] <= p1[0]
								&& p1[2] <= boundary1[0] + boundary1[5]));
		} else {
			fprintf(stderr,
					"ERROR: Cannot determine whether point is on the face of a %s.\n",
					boundaryString(boundary1Type));
			return false;
		}
	default:
		fprintf(stderr,
				"ERROR: Cannot determine whether point is on a face.\n");
		return false;
	}
}

// Do 2 boundaries share the same given surface?
// If so, faceShared specifies where they overlap
// This function is distinct from bBoundaryAdjacent because the shared
// face must be the same on both boundaries (e.g., lower x)
bool bSharedSurface(const int boundary1Type, const double boundary1[],
		const int boundary2Type, const double boundary2[], const short faceID,
		double faceShared[], const double error) {
	short i;
	short dim[2];
	switch (boundary1Type) {
	case RECTANGLE:
		switch (boundary2Type) {
		case RECTANGLE:
			// dim[0] will define the plane that the rectangles are on
			// dim[1] will define the shared varying coordinate

			// What plane are we on?
			if (boundary1[0] == boundary1[1])
				dim[0] = 0;
			else if (boundary1[2] == boundary1[3])
				dim[0] = 2;
			else if (boundary1[4] == boundary1[5])
				dim[0] = 4;

			// Is specified face normal to rectangle perimeter?
			// Are rectangles defined on the same plane?
			if (faceID == dim[0] || faceID == (dim[0] + 1)
					|| boundary2[dim[0]] != boundary2[dim[0] + 1])
				return false; // Shared face not possible

			switch (faceID) {
			case 0:
			case 1:
				if (dim[0] == 2)
					dim[1] = 4;
				if (dim[0] == 4)
					dim[1] = 2;
				break;
			case 2:
			case 3:
				if (dim[0] == 0)
					dim[1] = 4;
				if (dim[0] == 4)
					dim[1] = 0;
				break;
			case 4:
			case 5:
				if (dim[0] == 0)
					dim[1] = 2;
				if (dim[0] == 2)
					dim[1] = 0;
				break;
			default:
				fprintf(stderr,
						"ERROR: Specified face %d invalid for 2 Rectangles.\n",
						faceID);
				return false;
			}

			// Is the line actually shared?
			if (fabs(boundary1[faceID] - boundary2[faceID]) > error)
				return false; // Lines are different

			if (boundary1[dim[1]] >= boundary2[dim[1] + 1]
					|| boundary1[dim[1] + 1] <= boundary2[dim[1]])
				return false; // The segments do not overlap

			// We have overlap. Determine shared segment
			for (i = 0; i < 6; i++)
				faceShared[i] = boundary1[i];

			if (boundary1[dim[1]] < boundary2[dim[1]])
				faceShared[dim[1]] = boundary2[dim[1]];
			else
				faceShared[dim[1]] = boundary1[dim[1]];

			if (boundary1[dim[1] + 1] < boundary2[dim[1] + 1])
				faceShared[dim[1] + 1] = boundary1[dim[1] + 1];
			else
				faceShared[dim[1] + 1] = boundary2[dim[1] + 1];

			return true;
		default:
			fprintf(stderr,
					"ERROR: Boundary types %s and %s are not allowed to share a surface.\n",
					boundaryString(boundary1Type),
					boundaryString(boundary2Type));
			return false;
		}
	case RECTANGULAR_BOX:
		switch (boundary2Type) {
		case RECTANGULAR_BOX:
			// dim[0] and dim[1] will define the plane that the shared surface
			// would be on (if it exists)
			switch (faceID) {
			case 0:
			case 1:
				dim[0] = 2;
				dim[1] = 4;
				break;
			case 2:
			case 3:
				dim[0] = 0;
				dim[1] = 4;
				break;
			case 4:
			case 5:
				dim[0] = 0;
				dim[1] = 2;
				break;
			default:
				fprintf(stderr,
						"ERROR: Specified face %d invalid for 2 Rectanglular Boxes.\n",
						faceID);
				return false;
			}

			// Are the 2 faces on the same plane?
			if (fabs(boundary1[faceID] - boundary2[faceID]) > error)
				return false; // Planes are different

			// Do the 2 faces overlap?
			if (boundary1[dim[0]] >= boundary2[dim[0] + 1]
					|| boundary1[dim[0] + 1] <= boundary2[dim[0]]
					|| boundary1[dim[1]] >= boundary2[dim[1] + 1]
					|| boundary1[dim[1] + 1] <= boundary2[dim[1]])
				return false; // The faces do not overlap

			// We have overlap. Determine shared rectangle
			for (i = 0; i < 6; i++)
				faceShared[i] = boundary1[i];

			for (i = 0; i < 2; i++) {
				if (boundary1[dim[i]] < boundary2[dim[i]])
					faceShared[dim[i]] = boundary2[dim[i]];
				else
					faceShared[dim[i]] = boundary1[dim[i]];

				if (boundary1[dim[i] + 1] < boundary2[dim[i] + 1])
					faceShared[dim[i] + 1] = boundary1[dim[i] + 1];
				else
					faceShared[dim[i] + 1] = boundary2[dim[i] + 1];
			}

			return true;
		default:
			fprintf(stderr,
					"ERROR: Boundary types %s and %s are not allowed to share a surface.\n",
					boundaryString(boundary1Type),
					boundaryString(boundary2Type));
			return false;
		}
	case SPHERE: // Only one face on a sphere; no need to check faceID
		switch (boundary2Type) {
		case SPHERE:
			for (short i = 0; i < 3; i++) {
				if (boundary1[i] == boundary2[i])
					faceShared[i] = boundary1[i];
				else
					return false;
			}
			return true;
		default:
			fprintf(stderr,
					"ERROR: Boundary types %s and %s are not allowed to share a surface.\n",
					boundaryString(boundary1Type),
					boundaryString(boundary2Type));
			return false;
		}
	default:
		fprintf(stderr, "ERROR: Boundary type %s invalid to share a surface.\n",
				boundaryString(boundary1Type));
		return false;
	}
}

// Record specified face of boundary TODO add cylinders
void recordFace(const int boundary1Type, const double boundary1[],
		const short faceID, double boundaryFace[]) {

	switch (boundary1Type) {
	case RECTANGULAR_BOX:
	case RECTANGLE:
		switch (faceID) {
		case 0: // lower yz plane
			boundaryFace[0] = boundary1[0];
			boundaryFace[1] = boundary1[0];
			boundaryFace[2] = boundary1[2];
			boundaryFace[3] = boundary1[3];
			boundaryFace[4] = boundary1[4];
			boundaryFace[5] = boundary1[5];
			return;
		case 1: // upper yz plane
			boundaryFace[0] = boundary1[1];
			boundaryFace[1] = boundary1[1];
			boundaryFace[2] = boundary1[2];
			boundaryFace[3] = boundary1[3];
			boundaryFace[4] = boundary1[4];
			boundaryFace[5] = boundary1[5];
			return;
		case 2: // lower xz plane
			boundaryFace[0] = boundary1[0];
			boundaryFace[1] = boundary1[1];
			boundaryFace[2] = boundary1[2];
			boundaryFace[3] = boundary1[2];
			boundaryFace[4] = boundary1[4];
			boundaryFace[5] = boundary1[5];
			return;
		case 3: // upper xz plane
			boundaryFace[0] = boundary1[0];
			boundaryFace[1] = boundary1[1];
			boundaryFace[2] = boundary1[3];
			boundaryFace[3] = boundary1[3];
			boundaryFace[4] = boundary1[4];
			boundaryFace[5] = boundary1[5];
			return;
		case 4: // lower xy plane
			boundaryFace[0] = boundary1[0];
			boundaryFace[1] = boundary1[1];
			boundaryFace[2] = boundary1[2];
			boundaryFace[3] = boundary1[3];
			boundaryFace[4] = boundary1[4];
			boundaryFace[5] = boundary1[4];
			return;
		case 5: // lower xy plane
			boundaryFace[0] = boundary1[0];
			boundaryFace[1] = boundary1[1];
			boundaryFace[2] = boundary1[2];
			boundaryFace[3] = boundary1[3];
			boundaryFace[4] = boundary1[5];
			boundaryFace[5] = boundary1[5];
			return;
		default:
			fprintf(stderr, "ERROR: Face ID %d invalid for a %s.\n", faceID,
					boundaryString(boundary1Type));
			return;
		}
		return;
	case SPHERE:
		boundaryFace[0] = boundary1[0];
		boundaryFace[1] = boundary1[1];
		boundaryFace[2] = boundary1[2];
		boundaryFace[3] = boundary1[3];
		return;
	default:
		fprintf(stderr, "ERROR: Cannot record the face boundary of shape %s.\n",
				boundaryString(boundary1Type));
		return;
	}
}

// What is the value of the plane equation for a given point?
double planeEquation(const double point[3], const double plane[4]) {
	return point[0] * plane[0] + point[1] * plane[1] + point[2] * plane[2]
			+ plane[3];
}

// Reflect point against a boundary.
// oldPoint is used to determine direction of reflection if needed
// bReflectInside indicates whether point should be reflected into boundary
// Returns false if point did not intersect boundary
bool reflectPoint(const double oldPoint[3], const double L[3],
		const double length, const double curPoint[3], double newPoint[3],
		double intersectPoint[3], short * planeID, const int boundary1Type,
		const double boundary1[], bool bReflectInside) {
	double dist; // Distance from oldPoint to boundary along line to curPoint
	double LDotCenterToOld; // Dot product of L and centerToOld
	double pIntMinusC[3];
	double dDistance; // Distance from intersectPoint to newPoint

	newPoint[0] = curPoint[0];
	newPoint[1] = curPoint[1];
	newPoint[2] = curPoint[2];

	if (!bLineHitBoundary(oldPoint, L, length, boundary1Type, boundary1,
			planeID, *planeID, bReflectInside, &dist, intersectPoint)) // Line did not hit the boundary that it needs to reflect off of
			{
		// We should just lock boundary closest to endPoint
		if (!bLineHitBoundary(oldPoint, L, INFINITY, boundary1Type, boundary1,
				planeID, *planeID, bReflectInside, &dist, intersectPoint)) // Assume that point is already at boundary we want to reflect off of
				{
			// Just keep point at start
			intersectPoint[0] = oldPoint[0];
			intersectPoint[1] = oldPoint[1];
			intersectPoint[2] = oldPoint[2];
		} // Else line does eventually hit boundary. Just place point at that intersection
		newPoint[0] = intersectPoint[0];
		newPoint[1] = intersectPoint[1];
		newPoint[2] = intersectPoint[2];
		return false;
	}

	switch (boundary1Type) {
	case RECTANGULAR_BOX:
		switch (*planeID) {
		case 0:
			// Reflect off of lower x
			newPoint[0] = boundary1[0] + boundary1[0] - curPoint[0];
			return true;
		case 1:
			// Reflect off of upper x
			newPoint[0] = boundary1[1] + boundary1[1] - curPoint[0];
			return true;
		case 2:
			// Reflect off of lower y
			newPoint[1] = boundary1[2] + boundary1[2] - curPoint[1];
			return true;
		case 3:
			// Reflect off of upper y
			newPoint[1] = boundary1[3] + boundary1[3] - curPoint[1];
			return true;
		case 4:
			// Reflect off of lower z
			newPoint[2] = boundary1[4] + boundary1[4] - curPoint[2];
			return true;
		case 5:
			// Reflect off of upper z
			newPoint[2] = boundary1[5] + boundary1[5] - curPoint[2];
			return true;
		default:
			fprintf(stderr,
					"WARNING: Plane intersection ID %d invalid for a %s.\n",
					*planeID, boundaryString(boundary1Type));
			return false;
		}
	case SPHERE:
		*planeID = 0; // There's only one surface on a sphere

		pIntMinusC[0] = intersectPoint[0] - boundary1[0];
		pIntMinusC[1] = intersectPoint[1] - boundary1[1];
		pIntMinusC[2] = intersectPoint[2] - boundary1[2];

		dDistance = 2
				* ((curPoint[0] - intersectPoint[0]) * pIntMinusC[0]
						+ (curPoint[1] - intersectPoint[1]) * pIntMinusC[1]
						+ (curPoint[2] - intersectPoint[2]) * pIntMinusC[2])
				/ (squareDBL(pIntMinusC[0]) + squareDBL(pIntMinusC[1])
						+ squareDBL(pIntMinusC[2]));

		newPoint[0] -= dDistance * pIntMinusC[0];
		newPoint[1] -= dDistance * pIntMinusC[1];
		newPoint[2] -= dDistance * pIntMinusC[2];

		return true;
	case CYLINDER: //TODO reflections on the mantle taken from sphere and reduced to 2D, check!
		//TODO a transformation to do all the checks just once would be fine!
		if (boundary1[4] == PLANE_XY) {
			switch (*planeID) {
			case 4:
				// Reflect off of lower z
				newPoint[2] = boundary1[2] + boundary1[2] - curPoint[2];
				return true;
			case 5:
				// Reflect off of upper z
				newPoint[2] = boundary1[2] + boundary1[5] + boundary1[2]
						+ boundary1[5] - curPoint[2];
				return true;
				//all other surfaces are a reflection on the mantle
			case 0:
			case 1:
			case 2:
			case 3:
				pIntMinusC[0] = intersectPoint[0] - boundary1[0];
				pIntMinusC[1] = intersectPoint[1] - boundary1[1];
				dDistance = 2
						* ((curPoint[0] - intersectPoint[0]) * pIntMinusC[0]
								+ (curPoint[1] - intersectPoint[1])
										* pIntMinusC[1])
						/ (squareDBL(pIntMinusC[0]) + squareDBL(pIntMinusC[1]));
				newPoint[0] -= dDistance * pIntMinusC[0];
				newPoint[1] -= dDistance * pIntMinusC[1];
				return true;
			default:
				fprintf(stderr,
						"WARNING: Plane intersection ID %d invalid for a %s.\n",
						*planeID, boundaryString(boundary1Type));
				return false;
			}
		} else if (boundary1[4] == PLANE_XZ) {
			switch (*planeID) {
			case 2:
				// Reflect off of lower y
				newPoint[1] = boundary1[1] + boundary1[0] - curPoint[1];
				return true;
			case 3:
				// Reflect off of upper y
				newPoint[1] = boundary1[1] + boundary1[5] + boundary1[1]
						+ boundary1[5] - curPoint[1];
				return true;
				//all other surfaces are a reflection on the mantle
			case 0:
			case 1:
			case 4:
			case 5:
				pIntMinusC[0] = intersectPoint[0] - boundary1[0];
				pIntMinusC[2] = intersectPoint[2] - boundary1[2];
				dDistance = 2
						* ((curPoint[0] - intersectPoint[0]) * pIntMinusC[0]
								+ (curPoint[2] - intersectPoint[2])
										* pIntMinusC[2])
						/ (squareDBL(pIntMinusC[0]) + squareDBL(pIntMinusC[2]));
				newPoint[0] -= dDistance * pIntMinusC[0];
				newPoint[2] -= dDistance * pIntMinusC[2];
				return true;
			default:
				fprintf(stderr,
						"WARNING: Plane intersection ID %d invalid for a %s.\n",
						*planeID, boundaryString(boundary1Type));
				return false;
			}
		} else if (boundary1[4] == PLANE_YZ) {
			switch (*planeID) {
			case 0:
				// Reflect off of lower x
				newPoint[0] = boundary1[0] + boundary1[0] - curPoint[0];
				return true;
			case 1:
				// Reflect off of upper x
				newPoint[0] = boundary1[0] + boundary1[5] + boundary1[0]
						+ boundary1[5] - curPoint[0];
				return true;
				//all other surfaces are a reflection on the mantle
			case 2:
			case 3:
			case 4:
			case 5:
				pIntMinusC[1] = intersectPoint[1] - boundary1[1];
				pIntMinusC[2] = intersectPoint[2] - boundary1[2];
				dDistance = 2
						* ((curPoint[1] - intersectPoint[1]) * pIntMinusC[1]
								+ (curPoint[2] - intersectPoint[2])
										* pIntMinusC[2])
						/ (squareDBL(pIntMinusC[1]) + squareDBL(pIntMinusC[2]));
				newPoint[1] -= dDistance * pIntMinusC[1];
				newPoint[2] -= dDistance * pIntMinusC[2];
				return true;
			default:
				fprintf(stderr,
						"WARNING: Plane intersection ID %d invalid for a %s.\n",
						*planeID, boundaryString(boundary1Type));
				return false;
			}
		} else {
			fprintf(stderr, "ERROR: Cannot reflect a point off of a %s.\n",
					boundaryString(boundary1Type));
			return false;
		}
	default:
		fprintf(stderr,
				"ERROR: Cannot reflect a point off of an unknown region.\n");
		return false;
	}
}

// "Push" a point along a line
void pushPoint(double p1[3], double p2[3], const double dist, const double L[3]) {
	p2[0] = p1[0] + dist * L[0];
	p2[1] = p1[1] + dist * L[1];
	p2[2] = p1[2] + dist * L[2];
}

// Determine distance from point to a boundary
double distanceToBoundary(const double point[3], const int boundary1Type,
		const double boundary1[]) {
	double dist = 0.;
	double dist2;

	switch (boundary1Type) {
	case RECTANGULAR_BOX:
		if (bPointInBoundary(point, boundary1Type, boundary1)) {
			// Point is inside box find closest face
			dist = point[0] - boundary1[0];
			dist2 = boundary1[1] - point[0];
			if (dist2 < dist)
				dist = dist2;
			dist2 = point[1] - boundary1[2];
			if (dist2 < dist)
				dist = dist2;
			dist2 = boundary1[3] - point[1];
			if (dist2 < dist)
				dist = dist2;
			dist2 = point[2] - boundary1[4];
			if (dist2 < dist)
				dist = dist2;
			dist2 = boundary1[5] - point[2];
			if (dist2 < dist)
				dist = dist2;
			return dist;
		} else // Point is outside box
		{
			if (point[0] < boundary1[0])
				dist += squareDBL(boundary1[0] - point[0]);
			else if (point[0] > boundary1[1])
				dist += squareDBL(boundary1[1] - point[0]);
			if (point[1] < boundary1[2])
				dist += squareDBL(boundary1[2] - point[1]);
			else if (point[1] > boundary1[3])
				dist += squareDBL(boundary1[3] - point[1]);
			if (point[2] < boundary1[4])
				dist += squareDBL(boundary1[4] - point[2]);
			else if (point[2] > boundary1[5])
				dist += squareDBL(boundary1[5] - point[2]);
			return sqrt(dist);
		}
		return 0.;
	case SPHERE:
		dist = pointDistance(point, boundary1) - boundary1[3];
		if (dist < 0)
			dist = -dist;
		return dist;
	case CYLINDER: //TODO changed, check and finish
	default:
		fprintf(stderr,
				"ERROR: Cannot determine the distance from a point to a %s.\n",
				boundaryString(boundary1Type));
		return 0.;
	}
}

// Determine boundary of intersection of two boundaries
//  Only valid for rectangular boundaries (rectangles or boxes) or spherical intersections.
int intersectBoundary(const int boundary1Type, const double boundary1[],
		const int boundary2Type, const double boundary2[],
		double intersection[6]) {

	if ((boundary1Type == RECTANGULAR_BOX || boundary1Type == RECTANGLE)
			&& (boundary2Type == RECTANGULAR_BOX || boundary2Type == RECTANGLE)) {
		intersection[0] =
				(boundary1[0] > boundary2[0]) ? boundary1[0] : boundary2[0];
		intersection[1] =
				(boundary1[1] < boundary2[1]) ? boundary1[1] : boundary2[1];
		intersection[2] =
				(boundary1[2] > boundary2[2]) ? boundary1[2] : boundary2[2];
		intersection[3] =
				(boundary1[3] < boundary2[3]) ? boundary1[3] : boundary2[3];
		intersection[4] =
				(boundary1[4] > boundary2[4]) ? boundary1[4] : boundary2[4];
		intersection[5] =
				(boundary1[5] < boundary2[5]) ? boundary1[5] : boundary2[5];
		if (boundary1Type == RECTANGLE && boundary2Type == RECTANGLE)
			return RECTANGLE;
		else
			return RECTANGULAR_BOX;
	} else if (boundary1Type == SPHERE || boundary2Type == SPHERE) {
		// At least one of the boundaries is a sphere. One boundary must be
		// contained fully within the other boundary
		if (bBoundarySurround(boundary1Type, boundary1, boundary2Type,
				boundary2, 0.)) {
			// boundary 1 is within boundary 2
			intersection[0] = boundary1[0];
			intersection[1] = boundary1[1];
			intersection[2] = boundary1[2];
			intersection[3] = boundary1[3];
			intersection[4] = boundary1[4];
			intersection[5] = boundary1[5];
			return boundary1Type;
		} else if (bBoundarySurround(boundary2Type, boundary2, boundary1Type,
				boundary1, 0.)) {
			// boundary 2 is within boundary 1
			intersection[0] = boundary2[0];
			intersection[1] = boundary2[1];
			intersection[2] = boundary2[2];
			intersection[3] = boundary2[3];
			intersection[4] = boundary2[4];
			intersection[5] = boundary2[5];
			return boundary2Type;
		} else if (!bBoundaryIntersect(boundary2Type, boundary2, boundary1Type,
				boundary1, 0.)) {
			// Boundaries do not intersect at all
			intersection[0] = 0.;
			intersection[1] = 0.;
			intersection[2] = 0.;
			intersection[3] = 0.;
			intersection[4] = 0.;
			intersection[5] = 0.;
			return RECTANGULAR_BOX;
		} else {
			// Intersection is invalid
			fprintf(stderr,
					"ERROR: Intersection of two boundaries is invalid. At least one boundary is spherical and hits the other boundary.\n");
			return UNDEFINED_SHAPE;
		}
	} else if (boundary1Type == CYLINDER && boundary2Type == CYLINDER) {
		//intersection boundary can only be calculated if the cylinders have the same orientation
		//and the cross section of one is in the other
		if (boundary1[4] == boundary2[4]) {

			int along = 0;
			int across1 = 0;
			int across2 = 0;

			if (boundary2[4] == PLANE_XY) {
				across1 = 0;
				across2 = 1;
				along = 2;
			} else if (boundary2[4] == PLANE_XZ) {
				across1 = 0;
				along = 1;
				across2 = 2;
			} else if (boundary2[4] == PLANE_YZ) {
				along = 0;
				across1 = 1;
				across2 = 2;
			} else {
				fprintf(stderr,
						"ERROR: Cannot determine the orientation of a %s.\n",
						boundaryString(boundary1Type));
				return false;
			}

			double centerDistance = sqrt(
					squareDBL(boundary1[across1] - boundary2[across1])
							+ squareDBL(
									boundary1[across2] - boundary2[across2]));
			if (centerDistance >= boundary1[3] + boundary2[3]) {
				// no radial overlap
				intersection[0] = 0.;
				intersection[1] = 0.;
				intersection[2] = 0.;
				intersection[3] = 0.;
				intersection[4] = 0.;
				intersection[5] = 0.;
				return RECTANGULAR_BOX;
			} else if (centerDistance <= boundary1[3] - boundary2[3]) {
				//circle area of cylinder 2 inside that of cylinder 1 (or both equal)
				intersection[across1] = boundary2[across1];
				intersection[across2] = boundary2[across2];
				intersection[3] = boundary2[3];
				intersection[along] =
						(boundary1[along] > boundary2[along]) ?
								boundary1[along] : boundary2[along];
				intersection[4] = boundary2[4];
				intersection[5] = fmin(boundary1[along] + boundary1[5],
						boundary2[along] + boundary2[5]) - intersection[along];
				return CYLINDER;
			} else if (centerDistance <= boundary2[3] - boundary1[3]) {
				//circle area of cylinder 1 inside that of cylinder 2 (or both equal)
				intersection[across1] = boundary1[across1];
				intersection[across2] = boundary1[across2];
				intersection[3] = boundary1[3];
				intersection[along] =
						(boundary1[along] > boundary2[along]) ?
								boundary1[along] : boundary2[along];
				intersection[4] = boundary1[4];
				intersection[5] = fmin(boundary1[along] + boundary1[5],
						boundary2[along] + boundary2[5]) - intersection[along];
				return CYLINDER;
			} else {
				fprintf(stderr,
						"ERROR: Cannot determine the intersection of a %s and a %s if the intersection is no cylinder.\n",
						boundaryString(boundary1Type),
						boundaryString(boundary2Type));
				return false;
			}
		} else {
			fprintf(stderr,
					"ERROR: Cannot determine the intersection of a %s and a %s of different orientations.\n",
					boundaryString(boundary1Type),
					boundaryString(boundary2Type));
			return false;
		}

	} else if ((boundary1Type == CYLINDER && boundary2Type == RECTANGULAR_BOX)
			|| (boundary1Type == RECTANGULAR_BOX && boundary2Type == CYLINDER)) {
		double boundaryCylinder[6];
		double boundaryBox[6];
		//sort boundaries to unify calculations
		if (boundary1Type == CYLINDER) {
			for (int i = 0; i < 6; i++) {
				boundaryCylinder[i] = boundary1[i];
				boundaryBox[i] = boundary2[i];
			}
		} else {
			for (int i = 0; i < 6; i++) {
				boundaryCylinder[i] = boundary2[i];
				boundaryBox[i] = boundary1[i];
			}
		}
		//transform coordinates
		int along = 0;
		int across1 = 0;
		int across2 = 0;
		if (boundaryCylinder[4] == PLANE_XY) {
			across1 = 0;
			across2 = 1;
			along = 2;
		} else if (boundaryCylinder[4] == PLANE_XZ) {
			across1 = 0;
			along = 1;
			across2 = 2;
		} else if (boundaryCylinder[4] == PLANE_YZ) {
			along = 0;
			across1 = 1;
			across2 = 2;
		} else {
			fprintf(stderr,
					"ERROR: Cannot determine the orientation of a %s.\n",
					boundaryString(boundary1Type));
			return false;
		}
		// box inside circle area
		if (sqrt(
				squareDBL(boundaryBox[across1 * 2] - boundaryCylinder[across1])
						+ squareDBL(
								boundaryBox[across2 * 2]
										- boundaryCylinder[across2]))
				<= boundaryCylinder[3]
				&& sqrt(
						squareDBL(
								boundaryBox[across1 * 2 + 1]
										- boundaryCylinder[across1])
								+ squareDBL(
										boundaryBox[across2 * 2]
												- boundaryCylinder[across2]))
						<= boundaryCylinder[3]
				&& sqrt(
						squareDBL(
								boundaryBox[across1 * 2]
										- boundaryCylinder[across1])
								+ squareDBL(
										boundaryBox[across2 * 2 + 1]
												- boundaryCylinder[across2]))
						<= boundaryCylinder[3]
				&& sqrt(
						squareDBL(
								boundaryBox[across1 * 2 + 1]
										- boundaryCylinder[across1])
								+ squareDBL(
										boundaryBox[across2 * 2 + 1]
												- boundaryCylinder[across2]))
						<= boundaryCylinder[3]) {
			//cross section is a rectangle
			intersection[across1 * 2] = boundaryBox[across1 * 2];
			intersection[across1 * 2 + 1] = boundaryBox[across1 * 2 + 1];
			intersection[across2 * 2] = boundaryBox[across2 * 2];
			intersection[across2 * 2 + 1] = boundaryBox[across2 * 2 + 1];

			//length is the intersection of both lengths
			intersection[along * 2] = fmax(boundaryBox[along * 2],
					boundaryCylinder[along]);
			intersection[along * 2 + 1] = fmin(boundaryBox[along * 2 + 1],
					boundaryCylinder[along] + boundaryCylinder[5]);

			return RECTANGULAR_BOX;
		}
		//or the circle is completely in the cross section of the box
		else if (boundaryBox[across1 * 2]
				<= boundaryCylinder[across1] - boundaryCylinder[3]
				&& boundaryBox[across1 * 2 + 1]
						>= boundaryCylinder[across1] + boundaryCylinder[3]
				&& boundaryBox[across2 * 2]
						<= boundaryCylinder[across2] - boundaryCylinder[3]
				&& boundaryBox[across2 * 2 + 1]
						>= boundaryCylinder[across2] + boundaryCylinder[3]) {
			intersection[across1] = boundaryCylinder[across1];
			intersection[across2] = boundaryCylinder[across2];
			intersection[along] = fmax(boundaryCylinder[along],
					boundaryBox[along * 2]);
			intersection[3] = boundaryCylinder[3];
			intersection[4] = boundaryCylinder[4];
			intersection[5] = fmin(
					(boundaryCylinder[along] + boundaryCylinder[5]),
					boundaryBox[along * 2 + 1]) - intersection[along];
			return CYLINDER;
		} else {
			fprintf(stderr,
					"ERROR: Cannot determine the intersection of a %s and a %s if the cross section of one is not completely inside the other.\n",
					boundaryString(boundary2Type),
					boundaryString(boundary1Type));
			return UNDEFINED_SHAPE;
		}
	} else // Intersection for combination of boundary types is unknown
	{
		fprintf(stderr,
				"ERROR: Cannot determine the intersection of a %s and a %s.\n",
				boundaryString(boundary2Type), boundaryString(boundary1Type));
		return UNDEFINED_SHAPE;
	}

}

// Define unit vector pointing from one point to another
void defineLine(const double p1[3], const double p2[3], double L[3],
		double * length) {
	*length = sqrt(
			squareDBL(p2[0] - p1[0]) + squareDBL(p2[1] - p1[1])
					+ squareDBL(p2[2] - p1[2]));

	if (*length > 0.) {
		L[0] = (p2[0] - p1[0]) / (*length);
		L[1] = (p2[1] - p1[1]) / (*length);
		L[2] = (p2[2] - p1[2]) / (*length);
	} else {
		L[0] = 0.;
		L[1] = 0.;
		L[2] = 0.;
		*length = 0.;
	}
}

// Determine volume of boundary
double boundaryVolume(const int boundary1Type, const double boundary1[]) {
	switch (boundary1Type) {
	case RECTANGLE:
		if (boundary1[1] < boundary1[0] || boundary1[3] < boundary1[2]
				|| boundary1[5] < boundary1[4])
			return 0.;
		else if (boundary1[0] == boundary1[1])
			return (boundary1[5] - boundary1[4]) * (boundary1[3] - boundary1[2]);
		if (boundary1[2] == boundary1[3])
			return (boundary1[1] - boundary1[0]) * (boundary1[5] - boundary1[4]);
		if (boundary1[4] == boundary1[5])
			return (boundary1[1] - boundary1[0]) * (boundary1[3] - boundary1[2]);
	case RECTANGULAR_BOX:
		if (boundary1[1] < boundary1[0] || boundary1[3] < boundary1[2]
				|| boundary1[5] < boundary1[4])
			return 0.;
		else
			return (boundary1[1] - boundary1[0]) * (boundary1[3] - boundary1[2])
					* (boundary1[5] - boundary1[4]);
	case CIRCLE:
		return PI * squareDBL(boundary1[3]);
	case SPHERE:
		return 4 / 3 * PI * boundary1[3] * boundary1[3] * boundary1[3];
	case CYLINDER:
		return 2 * PI * boundary1[3] * boundary1[3] * boundary1[5];
	case LINE:
		return sqrt(
				squareDBL(boundary1[1] - boundary1[0])
						+ squareDBL(boundary1[3] - boundary1[2])
						+ squareDBL(boundary1[5] - boundary1[4]));
	default:
		fprintf(stderr, "ERROR: Cannot determine the volume of a %s.\n",
				boundaryString(boundary1Type));
		return 0;
	}
}

// Determine boundary surface Area
double boundarySurfaceArea(const int boundary1Type, const double boundary1[]) {
	double area = 0.;

	switch (boundary1Type) {
	case RECTANGLE:
		if (boundary1[1] < boundary1[0] || boundary1[3] < boundary1[2]
				|| boundary1[5] < boundary1[4])
			return 0.;

		area += 2 * (boundary1[1] - boundary1[0]);
		area += 2 * (boundary1[3] - boundary1[2]);
		area += 2 * (boundary1[5] - boundary1[4]);
		return area;

	case RECTANGULAR_BOX:
		if (boundary1[1] < boundary1[0] || boundary1[3] < boundary1[2]
				|| boundary1[5] < boundary1[4])
			return 0.;

		area += 2 * (boundary1[1] - boundary1[0])
				* (boundary1[3] - boundary1[2]);
		area += 2 * (boundary1[1] - boundary1[0])
				* (boundary1[5] - boundary1[4]);
		area += 2 * (boundary1[5] - boundary1[4])
				* (boundary1[3] - boundary1[2]);
		return area;

	case CIRCLE:
		return 2 * PI * boundary1[3];
	case SPHERE:
		return 4 * PI * boundary1[3] * boundary1[3];
	case CYLINDER:
		return 2 * PI * boundary1[3] * boundary1[3]
				+ 2 * PI * boundary1[3] * boundary1[5];
	default:
		fprintf(stderr, "ERROR: Boundary type %s invalid.\n",
				boundaryString(boundary1Type));
		return 0;
	}
}

// Find a random coordinate within the specified range
double uniformPoint(double rangeMin, double rangeMax) {
	return (rangeMin + (rangeMax - rangeMin) * mt_drand());
}

// Find a random coordinate within the specified boundary
void uniformPointVolume(double point[3], const int boundaryType,
		const double boundary1[], bool bSurface, const short planeID) {
	bool bNeedPoint;
	short curFace;
	double r, rSq;

	switch (boundaryType) {
	case RECTANGLE:
		if (bSurface) {
			curFace = (short) floor(4 * mt_drand());
			switch (planeID) {
			case PLANE_XY:
				point[2] = boundary1[4];
				switch (curFace) {
				case 0:
				case 1:
					point[0] = boundary1[curFace];
					point[1] = uniformPoint(boundary1[2], boundary1[3]);
					break;
				case 2:
				case 3:
					point[0] = uniformPoint(boundary1[0], boundary1[1]);
					point[1] = boundary1[curFace];
					break;
				}
				break;
			case PLANE_XZ:
				point[1] = boundary1[2];
				switch (curFace) {
				case 0:
				case 1:
					point[0] = boundary1[curFace];
					point[2] = uniformPoint(boundary1[4], boundary1[5]);
					break;
				case 2:
				case 3:
					point[0] = uniformPoint(boundary1[0], boundary1[1]);
					point[1] = boundary1[curFace + 2];
					break;
				}
				break;
			case PLANE_YZ:
				point[0] = boundary1[0];
				switch (curFace) {
				case 0:
				case 1:
					point[2] = boundary1[curFace + 4];
					point[1] = uniformPoint(boundary1[2], boundary1[3]);
					break;
				case 2:
				case 3:
					point[2] = uniformPoint(boundary1[4], boundary1[5]);
					point[1] = boundary1[curFace];
					break;
				}
				break;
			default:
				// Something went wrong
				fprintf(stderr,
						"ERROR: Cannot generate a uniform random point on plane %u of a rectangle.\n",
						planeID);
				return;
			}
			return;
		}
		switch (planeID) {
		case PLANE_XY:
			point[0] = uniformPoint(boundary1[0], boundary1[1]);
			point[1] = uniformPoint(boundary1[2], boundary1[3]);
			point[2] = boundary1[4];
			break;
		case PLANE_XZ:
			point[0] = uniformPoint(boundary1[0], boundary1[1]);
			point[1] = boundary1[2];
			point[2] = uniformPoint(boundary1[4], boundary1[5]);
			break;
		case PLANE_YZ:
			point[0] = boundary1[0];
			point[1] = uniformPoint(boundary1[2], boundary1[3]);
			point[2] = uniformPoint(boundary1[4], boundary1[5]);
			break;
		default:
			// Something went wrong
			fprintf(stderr,
					"ERROR: Cannot generate a uniform random point on plane %u of a rectangle.\n",
					planeID);
			return;
		}
	case RECTANGULAR_BOX:
		if (bSurface) {
			curFace = (short) floor(6 * mt_drand());
			switch (curFace) {
			case 0:
			case 1:
				point[0] = boundary1[curFace];
				point[1] = uniformPoint(boundary1[2], boundary1[3]);
				point[2] = uniformPoint(boundary1[4], boundary1[5]);
				break;
			case 2:
			case 3:
				point[0] = uniformPoint(boundary1[0], boundary1[1]);
				point[1] = boundary1[curFace];
				point[2] = uniformPoint(boundary1[4], boundary1[5]);
				break;
			case 4:
			case 5:
				point[0] = uniformPoint(boundary1[0], boundary1[1]);
				point[1] = uniformPoint(boundary1[2], boundary1[3]);
				point[2] = boundary1[curFace];
				break;
			}
			return;
		}
		point[0] = uniformPoint(boundary1[0], boundary1[1]);
		point[1] = uniformPoint(boundary1[2], boundary1[3]);
		point[2] = uniformPoint(boundary1[4], boundary1[5]);
		return;
	case CIRCLE:
		return;
	case SPHERE:
		// Use rejection method to create point in sphere
		bNeedPoint = true;
		while (bNeedPoint) {
			point[0] = mt_drand();
			point[1] = mt_drand();
			point[2] = mt_drand();

			rSq = squareDBL(point[0]) + squareDBL(point[1])
					+ squareDBL(point[2]);

			if (rSq < 1.) {
				// Found valid point. Scale as needed and randomize sign
				if (mt_drand() > 0.5)
					point[0] = -point[0];
				if (mt_drand() > 0.5)
					point[1] = -point[1];
				if (mt_drand() > 0.5)
					point[2] = -point[2];

				if (bSurface) {
					r = sqrt(rSq);
					point[0] /= r;
					point[1] /= r;
					point[2] /= r;
				}
				point[0] = boundary1[0] + point[0] * boundary1[3];
				point[1] = boundary1[1] + point[1] * boundary1[3];
				point[2] = boundary1[2] + point[2] * boundary1[3];

				bNeedPoint = false;
			}
		}
		return;
	default:
		fprintf(stderr,
				"ERROR: Cannot generate a uniform random point in a %s.\n",
				boundaryString(boundaryType));
		return;
	}
}

// Find distance between 2 3D points
double pointDistance(const double point1[3], const double point2[3]) {
	return sqrt(
			squareDBL(point2[0] - point1[0]) + squareDBL(point2[1] - point1[1])
					+ squareDBL(point2[2] - point1[2]));
}

// Square a double value
double squareDBL(double v) {
	return v * v;
}

// Return string with name of boundary
// Uses static memory strings in case output is not assigned
// to allocated memory
const char * boundaryString(const int boundaryType) {
	static char rectString[] = "Rectangle";
	static char boxString[] = "Rectangular Box";
	static char circleString[] = "Circle";
	static char sphereString[] = "Sphere";
	static char cylinderString[] = "Cylinder";
	static char emptyString[] = "";

	switch (boundaryType) {
	case RECTANGLE:
		return rectString;
	case RECTANGULAR_BOX:
		return boxString;
	case CIRCLE:
		return circleString;
	case SPHERE:
		return sphereString;
	case CYLINDER:
		return cylinderString;
	default:
		fprintf(stderr,
				"ERROR: Shape type %d does not have an associated name.\n",
				boundaryType);
		return emptyString;
	}
}
