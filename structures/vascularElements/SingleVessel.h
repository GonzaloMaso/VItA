/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SingleVessel.h
 *
 *  Created on: Mar 23, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef VASCULARELEMENTS_SINGLEVESSEL_H_
#define VASCULARELEMENTS_SINGLEVESSEL_H_

#include "AbstractVascularElement.h"
#include "../CCOCommonStructures.h"

#include <vtkSmartPointer.h>
#include <vtkLine.h>

#include <vector>

using namespace std;

/**
 * Composite vessel structure that for a single rectilinear vessel.
 */
class SingleVessel: public AbstractVascularElement {

public:

	/**	Points to be tested for branching in DEFORMABLE_PARENT mode. */
	static int bifurcationTests;
	/** VTK geomtric representation. */
	vtkSmartPointer<vtkLine> vtkSegment;
	/** Unique identifier for the vessel. */
	vtkIdType vtkSegmentId;

	/** Proximal position of the vessel. */
	point xProx;
	/** Distal position of the vessel. */
	point xDist;

	/** Bifurcation level from proximal to distal (Root is at level 0).*/
	int nLevel;
	/** Vessel radius. */
	double radius;
	/** Radius relative to the parent vessel (i.e. radius/parent_radius). For root, it is the radius value. */
	double beta;
	/** Distance between xProx and xDist. */
	double length;
	/** Reduced fluid-dynamic resistance. */
	double resistance;
	/** Local fluid-dynamic resistance. */
	double localResistance;
	/** Blood viscosity. */
	double viscosity;
	/** Flow in this vessel. */
	double flow;
	/** Pressure in the vessel. */
	double pressure;
	/** Generation step */
	long long int ID;
	/** Volume of this down tree branch.*/
	double treeVolume;

	SingleVessel();
	~SingleVessel();

	/**
	 * Fill @p branchingPoints with the potential points for branching in the current vessel, accorading to the @p branchingMode.
	 * @param branchingPoints Vector filled with the possible branching points of this vessel.
	 * @param xNew Terminal position of the new vessel to branch with.
	 */
	void getBranchingPoints(vector<point>* branchingPoints, point xNew);
	//	TOPOLOGY METHODS
	/**
	 * Returns its parental AbstractVascularElement.
	 * @return AbstractVascularElement proximally attached.
	 */
	AbstractVascularElement *getParent();
	/**
	 * Returns the SingleVessel that is parent for this vascular element.
	 * @param xp Point of parent vessel attachment.
	 * @return SingleVessel that is parent for this vascular element.
	 */
	SingleVessel *getParentVesselTo(point xp);
	/**
	 * Returns all the children vessels to this vascular element.
	 * @return Children vessels to this vascular element.
	 */
	vector<AbstractVascularElement *>& getChildren();

	/**
	 * Returns all SingleVessel that are children of this vascular element.
	 * @param xd Point of children vessels attachment.
	 * @return SingleVessel that is parent for this vascular element.
	 */
	vector<SingleVessel*> *getChildrenVesselTo(point xd);

	/**
	 * Returns all vessels in this vascular structure connected to point @p p.
	 * @param p Point at which all target vessels are connected to.
	 * @return Set of vessels connected to @p p.
	 */
	vector<SingleVessel*> *getVesselsConnectedTo(point p);

	/**
	 * Returns all the inner terminals in the current vascular element.
	 * @return Terminals in the current vascular element.
	 */
	long long int getTerminals();

	long long int getTerminals(TERMINAL_TYPE type);

	//	FUNCTIONALITY METHODS
	/**
	 * Returns the functional cost contribution of this vascular element. If any constraint
	 * --e.g. geometric constrains-- are violated the function must return INFINITY.
	 * @return Functional cost contribution.
	 */
	double getVolume();
	/**
	 * Updates vascular element pressure from children to parent direction.
	 */
	void updatePressure();
	/**
	 * Determine the flow at this terminal based on tree information (e.g. amount of terminals and inflow).
	 * @param tree Tree that contains this vessel.
	 * @return Amount of flow outgoing this terminal.
	 */
	double getTerminalFlow(double qProx, double qReserved, int commonTerminals);

	//	COMMUNICATION METHODS
	/**
	 * Returns the radius at the distal point of the vessel.
	 * @return Radius at the distal point of the vessel.
	 */
	double getDistalRadius();
	/**
	 * Returns the pressure at the proximal point of the vessel.
	 * @return Pressure at the proximal point of the vessel.
	 */
	double getProximalPressure();

	void saveVesselData(ofstream *treeFile);
	void saveVesselConnectivity(ofstream *treeFile);

	string coordToString();

};

#endif /* VASCULARELEMENTS_SINGLEVESSEL_H_ */
