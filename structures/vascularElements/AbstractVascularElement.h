/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractVascularElement.h
 *
 *  Created on: Mar 23, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef VASCULARELEMENTS_ABSTRACTVASCULARELEMENT_H_
#define VASCULARELEMENTS_ABSTRACTVASCULARELEMENT_H_

#include <fstream>
#include <vector>

#include "../CCOCommonStructures.h"

using namespace std;

class SingleVessel;

/**
 * Abstract vascular element that can be composed by a set of vessels.
 */
class AbstractVascularElement {
protected:
	/**	Vessels that constitute the current vascular element. */
	vector<SingleVessel *> vessels;
public:

	/**	Indicates the behavior of the terminal for flow distribution. */
	enum TERMINAL_TYPE {COMMON, RESERVED};

	/**	Indicates the behavior of the terminal for flow distribution. */
	TERMINAL_TYPE terminalType;

	/** Fraction of the qReserved assigned to this vessel. */
	double qReservedFraction;

	/** Parent vascular structure. */
	int stage;

	/** Parent vascular structure. */
	AbstractVascularElement *parent;

	/** The children vascular structures. */
	vector<AbstractVascularElement *> children;

	/** Indicates how and where this vessel can branch when it is parent. */
	enum BRANCHING_MODE {NO_BRANCHING, RIGID_PARENT, DEFORMABLE_PARENT, DISTAL_BRANCHING, ONLY_AT_PARENT_HOTSPOTS};

	/**	Indicates the hotspots where the vessel can bifurcate. */
	BRANCHING_MODE branchingMode;

	/** Indicates the vessel function that constraints where it can bifurcates. Perforators may bifurcates while
	 * they are partially outside the domain, but its partition point must be inside the domain. Transport vessels
	 * bifurcates outside the domain and distribution vessels bifurcates inside the domain.*/
	enum VESSEL_FUNCTION {DISTRIBUTION, PERFORATOR, TRANSPORT};

	/** Type of vessel function for the current vessel. */
	VESSEL_FUNCTION vesselFunction;

	/**
	 * Common constructor.
	 */
	AbstractVascularElement();
	/**
	 * Destructor constructor.
	 */
	virtual ~AbstractVascularElement();

	//	TOPOLOGY METHODS
	/**
	 * Sets the branching mode for the current vessel.
	 * @param mode Branching mode.
	 */
	void setBranchingMode(BRANCHING_MODE mode);

	/**
	 * Fill @p branchingPoints with the potential points for branching in the current vessel, accorading to the @p branchingMode.
	 * @param branchingPoints Vector filled with the possible branching points of this vessel.
	 */
	virtual void getBranchingPoints(vector<point>* branchingPoints, point xNew) = 0;

	/**
	 * Returns all vessels in the vascular element.
	 * @return Set of all vessels in the vascular element.
	 */
	vector<SingleVessel*>& getVessels();

	/**
	 * Returns its parental AbstractVascularElement.
	 * @return AbstractVascularElement proximally attached.
	 */
	virtual AbstractVascularElement *getParent() = 0;

	/**
	 * Returns the SingleVessel that is parent for this vascular element.
	 * @param xp Point of parent vessel attachment.
	 * @return SingleVessel that is parent for this vascular element.
	 */
	virtual SingleVessel *getParentVesselTo(point xp) = 0;

	/**
	 * Returns all the children vessels to this vascular element.
	 * @return Children vessels to this vascular element.
	 */
	virtual vector<AbstractVascularElement *>& getChildren() = 0;

	/**
	 * Adds a new child to the children vector.
	 * @param newChild New child to be inserted in the structure.
	 */
	void addChild(AbstractVascularElement *newChild);

	/**
	 * Remove all children.
	 */
	void removeChildren();

	/**
	 * Returns all SingleVessel that are children of this vascular element.
	 * @param xd Point of children vessels attachment.
	 * @return SingleVessel that is parent for this vascular element.
	 */
	virtual vector<SingleVessel*> *getChildrenVesselTo(point xd) = 0;

	/**
	 * Returns all vessels in this vascular structure connected to point @p p.
	 * @param p Point at which all target vessels are connected to.
	 * @return Set of vessels connected to @p p.
	 */
	virtual vector<SingleVessel*> *getVesselsConnectedTo(point p) = 0;

	/**
	 * Returns all the inner terminals in the current vascular element.
	 * @return Terminals in the current vascular element.
	 */
	virtual long long int getTerminals() = 0;

	virtual long long int getTerminals(TERMINAL_TYPE type) = 0;

	//	FUNCTIONALITY METHODS
	/**
	 * Returns the functional cost contribution of this vascular element. If any constraint
	 * --e.g. geometric constrains-- are violated the function must return INFINITY.
	 * @return Functional cost contribution.
	 */
	virtual double getVolume() = 0;
	/**
	 * Updates vascular element pressure from children to parent direction.
	 */
	virtual void updatePressure() = 0;

	//	COMMUNICATION METHODS
	/**
	 * Returns the radius at the distal point of the vessel.
	 * @return Radius at the distal point of the vessel.
	 */
	virtual double getDistalRadius() = 0;
	/**
	 * Returns the pressure at the proximal point of the vessel.
	 * @return Pressure at the proximal point of the vessel.
	 */
	virtual double getProximalPressure() = 0;

	//	I/O METHODS
	/**
	 * Writes the vessel data at the file associated to @p treeFile stream.
	 * @param treeFile Output file stream.
	 */
	virtual void saveVesselData(ofstream *treeFile) = 0;
	/**
	 * Writes connectivity of this vessel at the file associated to @p treeFile stream.
	 * @param treeFile Output file stream.
	 */
	virtual void saveVesselConnectivity(ofstream *treeFile) = 0;

};

#endif /* VASCULARELEMENTS_ABSTRACTVASCULARELEMENT_H_ */
