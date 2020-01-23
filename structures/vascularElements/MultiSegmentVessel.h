/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MultiSegmentVessel.h
 *
 *  Created on: Mar 23, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef VASCULARELEMENTS_MULTISEGMENTVESSEL_H_
#define VASCULARELEMENTS_MULTISEGMENTVESSEL_H_

#include "AbstractVascularElement.h"
#include "SingleVessel.h"

/**
 * Vessel composed by a set of SingleVessel structures in serie with only
 * bifurcations at its ends.
 */
class MultiSegmentVessel: public AbstractVascularElement {
public:
	/**	Indicates the hotspots where the vessel can bifurcate. */
	BRANCHING_MODE branchingMode;

	/** Distance between xProx and xDist. */
	double length;
	/** Reduced fluid-dynamic resistance. */
	double localResistance;
	/** Volume of this down tree branch.*/
	double treeVolume;

	/**
	 * Creates a MultiSegmentVessel from a vector of contiguous vessels ordered from proximal to distal position.
	 * @param innerSegments
	 */
	MultiSegmentVessel(vector<SingleVessel *> innerSegments);
	~MultiSegmentVessel();

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

};

#endif /* VASCULARELEMENTS_MULTISEGMENTVESSEL_H_ */
