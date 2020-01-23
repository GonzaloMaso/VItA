/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractCCOTree.h
 *
 *  Created on: Feb 5, 2018
 *      Author: gonzalo
 */

#ifndef ABSTRACTCCOTREE_H_
#define ABSTRACTCCOTREE_H_

#include <vtkPolyData.h>
#include <vtkCellLocator.h>
#include <vtkLine.h>
#include <vtkSmartPointer.h>


#include "../CCOCommonStructures.h"
#include "../domain/AbstractDomain.h"
#include "../../constrains/AbstractConstraintFunction.h"
#include "../../core/GeneratorData.h"

#include <vector>
#include <iomanip>
#include <fstream>

/**
 * Generic CCO Tree structure and required methods.
 */
class AbstractStructuredCCOTree {
protected:
	/** Wrapper with parameters associated to a tree generation process. */
	GeneratorData *instanceData;

	point xPerf;
	double qProx;
	AbstractConstraintFunction<double, int> *gam;
	AbstractConstraintFunction<double, int> *epsLim;
	AbstractConstraintFunction<double, int> *nu;
	//	Minimum aperture angle in the bifurcations (specified between 0 and pi)
	double minAngle;
	double refPressure;

	vessel *root;
	double psiFactor;
	double dp;
	long long nTerms;
	vtkSmartPointer<vtkPolyData> vtkTree;
	vtkSmartPointer<vtkCellLocator> vtkTreeLocator;
	vector<vessel *> segments;

	long long int pointCounter;

public:
	/**
	 * Constructor for tree with no parameters, convenient to create a tree and load its parameters later.
	 * @param instanceData Global parameters for instance behavior (amount of tries per bifurcation, point tries before dLim reduction, etc.).
	 */
	AbstractStructuredCCOTree(GeneratorData *instanceData);
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
	AbstractStructuredCCOTree(point xi, double qi, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim, AbstractConstraintFunction<double, int> *nu, double minAngle, double refPressure, GeneratorData *instanceData);
	/**
	 * Common destructor.
	 */
	virtual ~AbstractStructuredCCOTree();
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
	virtual vector<vessel *> getCloseSegments(point xNew, AbstractDomain *domain, int *nFound) = 0;
	/**
	 * Adds a new vessel to the CCO tree. @param xProx and @param xDist are the proximal and distal nodes of the new
	 * vessel and @param parent is the attachment parent vessel.
	 * @param xProx	Proximal point of the new vessel.
	 * @param xDist Distal point of the new vessel.
	 * @param parent	Parent to the new vessel.
	 */
	virtual void addVessel(point xProx, point xDist, vessel *parent) = 0;
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
	virtual int testVessel(point xNew, vessel *parent, AbstractDomain *domain, vector<vessel *> neighbors, double dlim, point *xBif, double *cost) = 0;
	/**
	 * Computes the pressure for the whole tree for a given reference pressure P_r (default P_r=0 Pa).
	 * @param root	Tree root.
	 */
	void computePressure(vessel *root);
	/**
	 * Prints the current tree vessel by vessel.
	 */
	virtual void print() = 0;
	/**
	 * Save in @p filename the current CCO tree in VTK format appropriate to resume execution.
	 * @param filename File name for VTP format.
	 */
	void storeVTK(string filename);
	/**
	 * Save in @p filename the current CCO tree in VTK format appropriate for visualization purpose (can not resume execution with this VTK).
	 * @param filename File name for VTP format.
	 * @param mode Establish some characteristics of the output file (e.g. the radius at the nodes instead at the cells).
	 */
	void storeVTK(string filename, unsigned int mode);
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
	vessel *getRoot();
	/**
	 * Getter of @p vtkTree.
	 * @return @p vtkTree.
	 */
	vtkSmartPointer<vtkPolyData> getVtkTree();
	/**
	 * Getter of @p segments.
	 * @return @p segments.
	 */
	const vector<vessel*>& getSegments() const;
	/**
	 * Getter of @p dp.
	 * @return @p dp.
	 */
	double getDp() const;
	/**
	 * Getter of @p minAngle.
	 * @return @p minAngle.
	 */
	double getMinAngle() const;
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
	 * Setter of @p minAngle.
	 * @param minAngle.
	 */
	void setMinAngle(double minAngle);
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

protected:
	/**
	 * Writes a string with the vessel attributes in a .cco file.
	 * @param currentVessel	Vessel to be saved.
	 * @param outFile	File writer for the .cco file.
	 */
	void saveVessel(vessel *currentVessel, ofstream *outFile);
	/**
	 * Writes a string with the tree attributes in a .cco file.
	 * @param outFile	File writer for the .cco file.
	 */
	virtual void saveTree(ofstream *outFile);
	/**
	 * Creates an ordered vector with all VTK identifiers of the vessels in the current tree.
	 * @return Vector with all VTK identifiers of the vessels in the current tree.
	 */
	vector<vessel*> createVTKIndex();
	/**
	 * Generates all VTK objects of the vessels loaded in the current tree.
	 * @param root Root of the tree.
	 */
	void generateVTKstructures(vessel *root);
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
	long long int countTerminals(vessel *root);
	/**
	 * Computes the total tree volume, if the geometric constraint is not satisfied, it returns INFINITY.
	 * @param root
	 * @param parentRadius
	 * @return Total tree volume.
	 */
	double computeTreeCost(vessel *root, double parentRadius);
};

#endif /* ABSTRACTCCOTREE_H_ */
