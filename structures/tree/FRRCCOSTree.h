/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * FixedRadiusRootCCOTreeLegacy.h
 *
 *  Created on: Dec 01, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef FIXEDRADIUSROOTCCOTREELEGACY_H_
#define FIXEDRADIUSROOTCCOTREELEGACY_H_

#include <vector>
#include <iostream>

#include "../CCOCommonStructures.h"
#include "../domain/AbstractDomain.h"

#include <vtkPolyData.h>
#include <vtkCellLocator.h>
#include <vtkLine.h>
#include <vtkSmartPointer.h>
#include "AbstractStructuredCCOTree.h"

using namespace std;

/**
 * Models a single CCO tree which contains the vessel's data and connectivity. The implementation
 * corresponds to the algorithm 3 of Queiroz Ph.D. Thesis. Current implementation is endowed with OpenMP
 * parallelism.
 * Specifications:
 * 	- Fixed domain;
 * 	- Given flow Q;
 * 	- Given pressure gradient;
 * 	- Given N terminals;
 * 	- Constant/variable symmetry per tree level;
 * 	- Constant/variable Murray coefficient per tree level;
 * 	- Constant/variable viscosity per tree level.
 */
class FRRCCOSTree : public AbstractStructuredCCOTree{

	/**	Root radius. */
	double rootRadius;

public:

	/**
	 * Common tree creator.
	 * @param xi	Root point.
	 * @param rootRadius	Root radius.
	 * @param qi	Flow at the root.
	 * @param gam	Murray law function.
	 * @param epsLim	Sibling vessels ratio function.
	 * @param nu	Viscosity function with respect to the tree level.
	 * @param minAngle	Minimum angle allowed.
	 */
	FRRCCOSTree(point xi, double rootRadius, double qi, AbstractConstraintFunction<double,int> *gam, AbstractConstraintFunction<double,int> *epsLim, AbstractConstraintFunction<double,int> *nu, double minAngle, double refPressure, GeneratorData *instanceData);
	/**
	 * Common destructor.
	 */
	~FRRCCOSTree();

	/**
	 *	Creates a copy from the tree without the vtk structures (neither from the tree nor from the vessels).
	 * @return	Copy from the tree
	 */
	FRRCCOSTree *clone();

	/**
	 * Returns the closest point in the CCOTree with respect to @p xNew point.
	 * @param xNew	Point from which the minimum distance is computed.
	 * @param xBif	Closest point in the tree to @p xNew.
	 * @param dist	Minimum distance between @p xNew and the tree.
	 */
	void getClosestTreePoint(point xNew, point *xBif, double *dist);

	/**
	 * Return the segments in a close neighborhood of @p xNew. The neighborhood is computed based on the perfusion
	 * volume indicated by @p domain.
	 * @param xNew	Center point of the neighborhood of interest.
	 * @param domain	Domain of the segments.
	 * @param nFound	Amount of segments in the neighborhood.
	 * @return	Array of segments in the neighborhood of @p xNew.
	 */
	vector<vessel *> getCloseSegments(point xNew, AbstractDomain *domain, int *nFound);
	/**
	 * Adds a new vessel to the CCO tree. @param xProx and @param xDist are the proximal and distal nodes of the new
	 * vessel and @param parent is the attachment parent vessel.
	 * @param xProx	Proximal point of the new vessel.
	 * @param xDist Distal point of the new vessel.
	 * @param parent	Parent to the new vessel.
	 */
	void addVessel(point xProx, point xDist, vessel *parent);
	/**
	 * For a given spatial point @p xNew test its connection with @p parent vessel. It must evaluate if the restrictions
	 * of geometry and symmetry are satisfied and also if it do not intersects with other vessel of this tree. It returns
	 * in @p cost the functional value with the inclusion of such vessel.
	 * @param xNew	Distal point for the new vessel to test.
	 * @param parent	Parent vessel to test the @p xNew point connection.
	 * @param domain	Tree domain.
	 * @param neighbors	Close neighbors used for intersection test.
	 * @param dlim	Not used in the current implementation.
	 * @param xBif	Proximal point for the new vessel that present the lower function cost.
	 * @param cost	Functional cost variation for the best bifurcation position.
	 * @return	If the connection of the tree with xNew is possible. If not @p cost is INFINITY.
	 */
	int testVessel(point xNew, vessel *parent, AbstractDomain *domain, vector<vessel *> neighbors, double dlim, point *xBif, double *cost);
	/**
	 * Prints the current tree node by node.
	 */
	void print();

	/**
	 * Returns the class name identifier.
	 * @return Class name.
	 */
	string getTreeName();

	/**
	 * Returns the root radius.
	 * @return Root radius.
	 */
	double getRootRadius();

protected:
	/**
	 * Returns a string with the tree atributes to create the .cco file.
	 * @param outFile	File writer for the .cco file.
	 */
	void saveTree(ofstream *outFile);

private:
	/**
	 * Clones the current tree.
	 * @param root	Root of the tree to clone.
	 * @param segments	Segments of the tree.
	 * @return Cloned subtree.
	 */
	vessel *cloneTree(vessel *root, vector<vessel *> *segments);
	/**
	 * Returns the updated functional cost due to the new segment inclusion.
	 * @param xNew	Proximal point of the new vessel.
	 * @param xTest Distal point of the new vessel.
	 * @param parent Parent to the new vessel.
	 */
	double evaluate(point xNew, point xTest, vessel *parent);
	/**
	 * Updates the tree values for the current topology in only one tree "in order" swept (O(N)).
	 * As the recursion deepens, the level number is computed for each element. As the
	 * recursion is returning, it computes the flow and resistance for the current node and the
	 * radius ratio for its childs.
	 * @param root Root vessel for the tree to update.
	 * @param tree Tree to update.
	 */
	void updateTree(vessel *root, FRRCCOSTree *tree);
	/**
	 * For a giving pair of beta between sibling of a parent vessel, it analyze the symmetry constrain given by
	 * epsLim function.
	 *
	 * @param beta1	Beta for the 1st sibling.
	 * @param beta2	Beta for the 2nd sibling.
	 * @param nLevel Tree level of the bifurcation.
	 * @return	If the simmetry constrain is satisfied.
	 */
	int isSymmetricallyValid(double beta1, double beta2, int nLevel);
	/**
	 * Checks if the angles of the parent vessel and the new vessel do not violate the minimum angle constraint.
	 * @param xBif	Bifurcation point between the new vessel and the distal part of the parent vessel (iCon).
	 * @param xNew	Distal point of the new vessel.
	 * @param parent	Parent vessel.
	 * @return	If the angles are higher than the minimum allowed.
	 */
	int areValidAngles(point xBif, point xNew, vessel *parent);
	/**
	 * Determines if the segment xBif-xNew is closer than dLim to any other segment of the tree (without being its parent
	 * vessel). Not used.
	 * @param xBif	Proximal point for the new segment.
	 * @param xNew	Distal point for the new segment.
	 * @param parent	Parent vessel of the new segment.
	 * @param dLim	Perfusion volume for each terminal at the current state of the tree.
	 * @return	If the middle point of the segment is sufficiently distant from the tree.
	 */
	int isOverlapped(point xBif, point xNew, vessel *parent, double dLim);
	/**
	 * It returns if the segment @param p1 - @param p2 intersects any vessel of the tree beside parent.
	 * @param p1	Extreme point 1 of the line.
	 * @param p2	Extreme point 2 of the line.
	 * @param parent	Vessel excluded from the checking.
	 * @param boundaryTol	Factor of line contraction at the extremes to avoid false intersections due to contact with anastomose.
	 * @return	If the segment intersects any segment of the tree.
	 */
	int isIntersectingVessels(point p1, point p2, vessel *parent, vector<vessel *> neighbors);
};

#endif /* FIXEDRADIUSROOTCCOTREELEGACY_H_ */
