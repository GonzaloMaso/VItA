/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SingleVesselCCOOTree.h
 *
 *	Jan 28, 2021: Important! This class needs cleaning and refactoring.
 *
 *  Created on: Mar 29, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef TREE_SINGLEVESSELCCOOTREE_H_
#define TREE_SINGLEVESSELCCOOTREE_H_

#include <iostream>
#include <string>
#include <vector>

#include "../../constrains/AbstractConstraintFunction.h"
#include "../vascularElements/SingleVessel.h"
#include "AbstractObjectCCOTree.h"

using namespace std;

/**
 * N-furcation tree with only SingleVessel elements as vascular elements.
 */
class SingleVesselCCOOTree: public AbstractObjectCCOTree {
	string filenameCCO;
	/**	Root radius. */
	double rootRadius;
	/** Convergence tolerance. */
	double variationTolerance;
	/**	Amount of non-common terminals. */
	long long int nCommonTerminals;
	//	FIXME These classes should not have this kind of permissions, must rework the architecture to a POO strategy.
	friend class PruningCCOOTree;
	friend class BreadthFirstPruning;
	friend class TreeMerger;
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
	 * @param resistanceVariationTolerance	Convergence tolerance for the iterative viscosity scheme.
	 */
	SingleVesselCCOOTree(point xi, double rootRadius, double qi, AbstractConstraintFunction<double,int> *gam, AbstractConstraintFunction<double,int> *epsLim,
			AbstractConstraintFunction<double,int> *nu, double refPressure, double resistanceVariationTolerance, GeneratorData *instanceData);

//	/**
//	 * Creates a new tree from the .cco file @p filename. The obtained tree has not vtkLine objects of the vessels.
//	 * @param filenameCCO Path to the .cco file.
//	 * @param filenameVTK Path to the .vtk file.
//	 */
//	SingleVesselCCOOTree(string filenameCCO, string filenameVTK, GeneratorData *instanceData);

	/**
	 * Creates a new tree from the .cco file @p filename in VItA format.
	 * @param filenameCCO Path to the .cco file.
	 * @param gam	Murray law function.
	 * @param epsLim	Sibling vessels ratio function.
	 * @param nu	Viscosity function with respect to the tree level.
	 */
	SingleVesselCCOOTree(string filenameCCO, GeneratorData *instanceData, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim,
			AbstractConstraintFunction<double, int> *nu);
	/**
	 * Creates a new tree from the .cco file @p filename in HeMoLab format.
	 * @param filenameCCO Path to the .cco file.
	 * @param qi	Flow at the root.
	 * @param gam	Murray law function.
	 * @param epsLim	Sibling vessels ratio function.
	 * @param nu	Viscosity function with respect to the tree level.
	 * @param minAngle	Minimum angle allowed.
	 * @param refPressure	Convergence tolerance for the iterative viscosity scheme.
	 * @param viscosityTolerance	Convergence tolerance for the iterative viscosity scheme.
	 */
	SingleVesselCCOOTree(string filenameCCO, GeneratorData* instanceData, double qi, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim,
			AbstractConstraintFunction<double, int> *nu, double refPressure, double viscosityTolerance);
	/**
	 * Creates a copy of the tree only with its parameters. Does not create vessel data.
	 * @param baseTree Base tree.
	 */
	SingleVesselCCOOTree(SingleVesselCCOOTree *baseTree);
	/**
	 * Common destructor.
	 */
	~SingleVesselCCOOTree();
	/**
	 *	Creates a copy from the tree without the vtk structures (neither from the tree nor from the vessels).
	 * @return	Copy from the tree
	 */
	SingleVesselCCOOTree *clone();
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
	vector<AbstractVascularElement *> getCloseSegments(point xNew, AbstractDomain *domain, int *nFound);

	/**
	 * Adds a new vessel to the CCO tree. @param xProx and @param xDist are the proximal and distal nodes of the new
	 * vessel and @param parent is the attachment parent vessel.
	 * @param xProx	Proximal point of the new vessel.
	 * @param xDist Distal point of the new vessel.
	 * @param parent	Parent to the new vessel.
	 * @param vesselFunction Vessel function of the added vessel.
	 */
	void addVessel(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction);

	//	FIXME This function probably should be part of other class
	void addVesselMergeFast(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, unordered_map<string, SingleVessel *>* stringToPointer);

	//	FIXME This function probably should be part of other class
	void addVesselMerge(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction, unordered_map<string, SingleVessel *>* stringToPointer);
	/** 
	 * Adds a vessel that has already been validated. This function is used by BreadthFirstPrunning.
	 * @param newVessel Vessel that will be added.
	 * @param originalVessel Vessel that has been validaded in a previous tree.
	 * @param copiedTo Mapping such that copiedTo[originalVessel] = copiedVessel;
	*/
	//	FIXME This function probably should be part of other class
	void addValitatedVessel(SingleVessel *newVessel, SingleVessel *originalVessel, unordered_map<SingleVessel *, SingleVessel *>& copiedTo);

	/** 
	 * Adds a vessel that has already been validated, but do not update the tree. This function is used by BreadthFirstPrunning.
	 * @param newVessel Vessel that will be added.
	 * @param originalVessel Vessel that has been validaded in a previous tree.
	 * @param copiedTo Mapping such that copiedTo[originalVessel] = copiedVessel;
	*/
	//	FIXME This function probably should be part of other class
	void addValitatedVesselFast(SingleVessel *newVessel, SingleVessel *originalVessel, unordered_map<SingleVessel *, SingleVessel *>& copiedTo);

//	/**
//	 * Adds a new vessel to the CCO tree as continuation of the pre-existent vessel @p parent. @param xDist is the distal nodes of the new
//	 * vessel and @param parent is the proximal attachment parent vessel.
//	 * @param xDist Distal point of the new vessel.
//	 * @param parent Parent to the new vessel.
//	 * @param mode Branching mode of the added vessel.
//	 * @param vesselFunction Vessel function of the added vessel.
//	 */
//	void addVessel(point xDist, AbstractVascularElement *parent, AbstractVascularElement::BRANCHING_MODE mode, AbstractVascularElement::VESSEL_FUNCTION vesselFunction);

	/**
	 * For a given spatial point @p xNew test its connection with @p parent vessel. It must evaluate if the restrictions
	 * of geometry and symmetry are satisfied and also if it do not intersects with other vessel of this tree. It returns
	 * in @p cost the functional variation by the inclusion of such vessel.
	 * @param xNew	Distal point for the new vessel to test.
	 * @param parent	Parent vessel to test the @p xNew point connection.
	 * @param domain	Tree domain.
	 * @param neighbors	Close neighbors used for intersection test.
	 * @param dlim	Not used in the current implementation.
	 * @param xBif	Proximal point for the new vessel that present the lower function cost.
	 * @param cost	Functional cost variation for the best bifurcation position.
	 * @return	If the connection of the tree with xNew is possible. If not @p cost is INFINITY.
	 */
	int testVessel(point xNew, AbstractVascularElement *parent, AbstractDomain *domain, vector<AbstractVascularElement *> neighbors, double dlim, point *xBif, double *cost);

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

	/**
	 * Removes all branches at stage @p stage such that do not branch in the next stages.
	 * @param stage
	 */
	void removeWitheredBranches(int stage);

	/**
	 * Removes a vessel and its subtree from the current tree.
	 * @param vessel Vessel to be removed.
	 */
	void remove(SingleVessel *vessel);

	/**
	 * Returns if the @p vessel at stage @p stage has no connexions with vessels at stage @p stage + 1.
	 * @param vessel Vessel to evaluate if it is withered or not.
	 * @param stage Stage of the vessel @p vessl
	 * @return If it is withered or not.
	 */
	bool isWithered(SingleVessel *vessel);

	string getFilenameCCO();

protected:
	/**
	 * Returns a string with the tree atributes to create the .cco file.
	 * @param outFile	File writer for the .cco file.
	 */
	void saveTree(ofstream *outFile);

private:
	/**
	 * Clones the subtree with parent vessel @p levels .
	 * @param root	Root of the tree to clone.
	 * @param segments	Segments of the tree.
	 * @return Cloned subtree.
	 */
	SingleVesselCCOOTree *cloneUpTo(int levels, SingleVessel *parent);
	/**
	 * Clones the subtree with parent vessel @p root recursively.
	 * @param root	Root of the tree to clone.
	 * @param segments	Segments of the tree.
	 * @return Cloned subtree.
	 */
	SingleVessel *cloneTree(SingleVessel *root, unordered_map<long long, AbstractVascularElement *> *segments);
	/**
	 * Returns a partial variation of the cost functional due to the new segment inclusion.
	 * @param xNew	Proximal point of the new vessel.
	 * @param xTest Distal point of the new vessel.
	 * @param parent Parent to the new vessel.
	 * @param dLim Minimum distance from the new vessel to the tree.
	 */
	double evaluate(point xNew, point xTest, SingleVessel *parent, double dLim);
	/**
	 * Returns a partial variation of the cost functional due to the new segment inclusion. This method is only used for DISTAL_BRANCHING vessels.
	 * @param xNew	Proximal point of the new vessel.
	 * @param parent Parent to the new vessel.
	 * @param dLim Minimum distance from the new vessel to the tree.
	 */
	double evaluate(point xNew, SingleVessel *parent, double dLim);
	/**
	 * Updates the tree values for the current topology in only one tree "in order" swept (O(N)).
	 * As the recursion deepens, the level number is computed for each element. As the
	 * recursion is returning, it computes the flow and resistance for the current node and the
	 * radius ratio for its childs.
	 * @param root Root vessel for the tree to update.
	 * @param tree Tree to update.
	 */
	void updateTree(SingleVessel *root, SingleVesselCCOOTree *tree);
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
	 * @param minAngle 	Minimum angle constraint.
	 * @return	If the angles are higher than the minimum allowed.
	 */
	int areValidAngles(point xBif, point xNew, SingleVessel *parent, double minAngle);
	/**
	 * Checks if the angle between the plane of the parent and sibling vessel and the new vessel satisfies the opening angle constraint.
	 * @param xBif	Bifurcation point between the new vessel and the distal part of the parent vessel (iCon).
	 * @param xNew	Distal point of the new vessel.
	 * @param parent	Parent vessel.
	 * @param minPlaneAngle 	Minimum opening angle constraint.
	 * @return	If the angles satisfies the minimum plane angle.
	 */
	int isValidOpeningAngle(point xBif, point xNew, SingleVessel *parent, double minPlaneAngle);
	/**
	 * It returns if the line @param p1 - @param p2 intersects any vessel of the tree beside parent.
	 * @param p1	Extreme point 1 of the line.
	 * @param p2	Extreme point 2 of the line.
	 * @param parent	Vessel excluded from the checking.
	 * @param boundaryTol	Factor of line contraction at the extremes to avoid false intersections due to contact with anastomose.
	 * @return	If the segment intersects any segment of the tree.
	 */
	int isIntersectingVessels(point p1, point p2, SingleVessel *parent, vector<AbstractVascularElement *> neighbors);
	/**
	 * Updates the tree beta values for the current vessel diameters.
	 * @param root	Tree root.
	 * @param parentRadius	1.0 if root is actually the tree root, or parent radius of root otherwise.
	 * @param maxBetaVariation	The maximum variation of the beta value due to the update.
	 */
	void updateTreeViscositiesBeta(SingleVessel *root, double *maxBetaVariation);
	/**
	 * Returns the viscosity estimated with the Fahraeus-Lindquist model.
	 * @param radius
	 * @return Viscosity value.
	 */
	double getNuFL(double radius);

	/**
	 * Creates the VTK lines and points associated to a HeMoLab file loaded.
	 */
	void createSegmentVtkLines(AbstractVascularElement *rootVessel);

	double getVariationTolerance();

};

#endif /* TREE_SINGLEVESSELCCOOTREE_H_ */
