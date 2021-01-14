/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * GenericVesselsTree.h
 *
 *  Created on: Mar 23, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef TREE_ABSTRACTOBJECTCCOTREE_H_
#define TREE_ABSTRACTOBJECTCCOTREE_H_

#include "../CCOCommonStructures.h"
#include "../vascularElements/AbstractVascularElement.h"
#include "../vascularElements/SingleVessel.h"
#include "../domain/AbstractDomain.h"
#include "../../constrains/AbstractConstraintFunction.h"
#include "../../core/GeneratorData.h"

#include <vtkPolyData.h>
#include <vtkCellLocator.h>
#include <vtkLine.h>
#include <vtkSmartPointer.h>

#include <vector>
#include <iomanip>
#include <fstream>
#include <unordered_map>

using namespace std;

/**
 * Abstract tree with vessel structures represented by objects. Its slower than AbstractStructuredCCOTree
 * although has more functionality and is more customizable.
 */
class AbstractObjectCCOTree {
	
protected:
	/** Wrapper with parameters associated to a tree generation process. */
	GeneratorData *instanceData;

	/**	Position of the tree inlet. */
	point xPerf;
	/** Proximal flow of the tree. */
	double qProx;
	/** Percentage of the flow reserved for RESERVED terminals. */
	double qReservedFactor;
	/** Constraint for the ratio between parent and child radius based on the Murray's Law. */
	AbstractConstraintFunction<double, int> *gam;
	/** Constraint for the symmetry between sibling vessels. */
	AbstractConstraintFunction<double, int> *epsLim;
	/** Constraint for the vessel viscosity. */
	AbstractConstraintFunction<double, int> *nu;
	/** Distal reference pressure of the tree. */
	double refPressure;
	/**	Tree root. */
	AbstractVascularElement *root;
	/**	Unused, factor used to estimate pressure in the tree. */
	double psiFactor;
	/**	Unused, pressure drop along the tree. */
	double dp;
	/**	Amount of terminals within the tree. */
	long long nTerms;

	/**	VTK representation of the tree. */
	vtkSmartPointer<vtkPolyData> vtkTree;
	/**	VTK locator for the tree. */
	vtkSmartPointer<vtkCellLocator> vtkTreeLocator;
	/**	Vascular elements (vessels or set of vessels) of the tree. */
	unordered_map<long long, AbstractVascularElement *> elements;

	/**	Amount of points consumed since the beginning of the CCO generation. */
	long long int pointCounter;
	/**	Stage actually being build. */
	int currentStage;
	/**	If the tree is in cm, otherwise is assumed that is in mm (default) */
	int isInCm;

	friend class PruningCCOOTree;

public:
	/**
	 * Constructor for tree with no parameters, convenient to create a tree and load its parameters later.
	 * @param instanceData Global parameters for instance behavior (amount of tries per bifurcation, point tries before dLim reduction, etc.).
	 */
	AbstractObjectCCOTree(GeneratorData *instanceData);
	/**
	 * Common constructor with all parameters necessary for tree generation.
	 * @param xi	Root entry point.
	 * @param qi	Flow at @p xi.
	 * @param gam	Function for the exponential coefficient of the Murray Law depending on the bifurcation level.
	 * @param epsLim	Function of the symmetry ratio (ratio between the radius of the smallest over the biggest sibling) depending on the bifurcation level.
	 * @param nu	Blood viscosity depending on the bifurcation level.
	 * @param minAngle	Lowest angle allowed at the new bifurcation.
	 * @param refPressure	Distal reference pressure for the tree.
	 * @param instanceData	Parameters for this stage generation.
	 */
	AbstractObjectCCOTree(point xi, double qi, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim, AbstractConstraintFunction<double, int> *nu, double refPressure, GeneratorData *instanceData);
	/**
	 * Common destructor.
	 */
	virtual ~AbstractObjectCCOTree();
	/**
	 * Returns the closest point in the CCOTree with respect to @p xNew point.
	 * @param xNew	Point from which the minimum distance is computed.
	 * @param xBif	Closest point in the tree to @p xNew.
	 * @param dist	Minimum distance between @p xNew and the tree.
	 */
	virtual void getClosestTreePoint(point xNew, point *xBif, double *dist) = 0;
	/**
	 * Return the segments in a close neighborhood of @p xNew. The neighborhood is computed based on the perfusion
	 * volume indicated by @p domain.
	 * @param xNew	Center point of the neighborhood of interest.
	 * @param domain	Domain of the segments.
	 * @param nFound	Amount of segments in the neighborhood.
	 * @return	Array of segments in the neighborhood of @p xNew.
	 */
	virtual vector<AbstractVascularElement *> getCloseSegments(point xNew, AbstractDomain *domain, int *nFound) = 0;
	/**
	 * Adds a new vessel to the CCO tree. @param xProx and @param xDist are the proximal and distal nodes of the new
	 * vessel and @param parent is the attachment parent vessel.
	 * @param xProx	Proximal point of the new vessel.
	 * @param xDist Distal point of the new vessel.
	 * @param parent	Parent to the new vessel.
	 */
	virtual void addVessel(point xProx, point xDist, AbstractVascularElement *parent, AbstractVascularElement::VESSEL_FUNCTION vesselFunction) = 0;
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
	virtual int testVessel(point xNew, AbstractVascularElement *parent, AbstractDomain *domain, vector<AbstractVascularElement *> neighbors, double dlim, point *xBif, double *cost) = 0;
	/**
	 * Computes the pressure for the whole tree for a given reference pressure P_r (default P_r=0 Pa).
	 * @param root	Tree root.
	 */
	void computePressure(AbstractVascularElement *root);
	/**
	 * Computes the total tree volume, if the geometric constraint is not satisfied, it returns INFINITY.
	 * @param root
	 * @return Total tree volume.
	 */
	double computeTreeCost(AbstractVascularElement *root);
	/**
	 * Prints the current tree vessel by vessel.
	 */
	virtual void print() = 0;
	/**
	 * Saves the tree data with exception of the vtk objects.
	 * @param filename	File to save the tree data.
	 */
	virtual void save(string filename);
	/**
	 * Returns the name of the tree class that implements the AbstractCCOTree.
	 * @return Tree class name.
	 */
	virtual string getTreeName() = 0;
	/**
	 * Returns a vector with each point of the tree as entry. Each point is implemented as a vector of three doubles.
	 * @return
	 */
	vector< vector<double> > getVertices();

	/**
	 * Returns a vector with one entry per line element of the tree. Each entry contains the identifiers of the points
	 * that constitutes the element.
	 * @return
	 */
	vector< vector<int> > getConnectivity();
	/**
	 * Prints points and lines contained in the VTK tree structure.
	 */
	void printVtkTree();

	/**
	 * Returns the amount of terminals in the current tree.
	 * @return Amount of tree terminals.
	 */
	long long int getNTerminals();
	/**
	 * Returns the amount of terminals of a given terminal type.
	 * @param type Terminal type.
	 * @return Amount of terminals of the given terminal type.
	 */
	long long int getNTerminals(AbstractVascularElement::TERMINAL_TYPE type);

	/**
	 * Getter of @p nTerms.
	 * @return @p nTerms
	 */
	long long int getNTerms();
	/**
	 * Getter of @p xProx.
	 * @return @p xProx
	 */
	point getXProx();
	/**
	 * Getter of @p qProx.
	 * @return @p qProx
	 */
	double getQProx();
	/**
	 * Getter of @p epsLim.
	 * @return @p epsLim.
	 */
	AbstractConstraintFunction<double, int> *getEpsLim();
	/**
	 * Getter of @p gam.
	 * @return @p gam.
	 */
	AbstractConstraintFunction<double, int> *getGam();
	/**
	 * Getter of @p nu.
	 * @return @p nu.
	 */
	AbstractConstraintFunction<double, int> *getNu();
	/**
	 * Getter of @p root.
	 * @return @p root.
	 */
	AbstractVascularElement *getRoot();
	/**
	 * Getter of @p vtkTree.
	 * @return @p vtkTree.
	 */
	vtkSmartPointer<vtkPolyData> getVtkTree();
	/**
	 * Getter of @p elements.
	 * @return @p elements.
	 */
	unordered_map<long long, AbstractVascularElement*>& getSegments();
	/**
	 * Getter of vessels within @p elements.
	 * @return @p elements.
	 */
	vector<SingleVessel *> getVessels();
	/**
	 * Getter of @p dp.
	 * @return @p dp.
	 */
	double getDp() const;
	/**
	 * Getter of @p instanceData.
	 * @return @p instanceData.
	 */
	GeneratorData* getInstanceData();
	/**
	 * Getter of @p pointCounter.
	 * @return @p pointCounter.
	 */
	long long int getPointCounter() const;
	/**
	 * Getter of @p currentStage.
	 * @return @p currentStage.
	 */
	int getCurrentStage() const;
	/**
	 * Getter of @p qReservedFactor.
	 * @return @p qReservedFactor.
	 */
	double getReservedFactor() const;
	/**
	 * Getter of @p isInCm.
	 * @return @p isInCm
	 */
	int getIsInCm() const;

	/**
	 * Setter of @p epsLim.
	 * @param epsLim.
	 */
	void setEpsLim(AbstractConstraintFunction<double, int> *epsLim);
	/**
	 * Setter of @p gam.
	 * @param gam.
	 */
	void setGam(AbstractConstraintFunction<double, int> *gam);
	/**
	 * Setter of @p nu.
	 * @param nu.
	 */
	void setNu(AbstractConstraintFunction<double, int> *nu);
	/**
	 * Setter of @p pointCounter.
	 * @param pointCounter.
	 */
	void setPointCounter(long long int pointCounter);
	/**
	 * Setter of @p instanceData.
	 * @param instanceData.
	 */
	void setInstanceData(GeneratorData* instanceData);
	/**
	 * Setter of @p currentStage.
	 * @param currentStage.
	 */
	void setCurrentStage(int currentStage);
	/**
	 * Setter of @p qReservedFactor.
	 * @param qReservedFactor.
	 */
	void setReservedFactor(double reservedFactor);
	/**
	 * Setter of @p isInCm.
	 * @param isInCm.
	 */
	void setIsInCm(int isInCm);

protected:
	/**
	 * Writes a string with the vessel attributes in a .cco file.
	 * @param currentVessel	Vessel to be saved.
	 * @param outFile	File writer for the .cco file.
	 */
	void saveVessel(AbstractVascularElement *currentVessel, ofstream *outFile);
	/**
	 * Writes a string with the tree attributes in a .cco file.
	 * @param outFile	File writer for the .cco file.
	 */
	virtual void saveTree(ofstream *outFile);
	/**
	 * Updates the @p vtkSegment field of each vessel in segments using their @p vtkSegmentId from
	 * the current @p vtkTree.
	 */
	void updateSegmentVtkLines();
	/**
	 * Counts recursively all terminals of the subtree with @p root as root.
	 * @param root Root of the subtree.
	 * @return	Amount of terminals in the subtree.
	 */
	long long int countTerminals(AbstractVascularElement *root);

	/**
	 * Counts recursively the terminals with type @p type of the subtree with @p root as root.
	 * @param root Root of the subtree.
	 * @return	Amount of terminals in the subtree.
	 */
	long long int countTerminals(AbstractVascularElement* root, AbstractVascularElement::TERMINAL_TYPE type);

private:
	void saveVessels(AbstractVascularElement *root, ofstream *treeFile);
	void saveConnectivity(AbstractVascularElement *root, ofstream *treeFile);
};

#endif /* TREE_ABSTRACTOBJECTCCOTREE_H_ */
