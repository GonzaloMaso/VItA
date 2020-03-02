/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_ABSTRACT1DSTATMANIPULATOR_H_
#define STATS_ABSTRACT1DSTATMANIPULATOR_H_

#include <vector>

#include "../../structures/vascularElements/SingleVessel.h"
#include "../VesselObjectHandler.h"

/**
 * Abstract class for the computation of statistical quantities of a specific attribute for an array of vessels.
 */
class Abstract1DStatManipulator {
protected:
	/**	Handler used to extract the vessel attributes. */
	VesselObjectHandler *handler;
public:
	/**
	 * Constructor.
	 */
	Abstract1DStatManipulator();
	/**
	 * Destructor.
	 */
	~Abstract1DStatManipulator();
	/**
	 * Computes the statistical quantity of interest.
	 * @param vessels	Array of vessels processed to obtain the quantity of interest.
	 * @param att		Attribute over which the quantity is computed.
	 * @return	Statistical quantity of interest.
	 */
	virtual vector<double> compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att) = 0;
};

#endif /* STATS_ABSTRACT1DSTATMANIPULATOR_H_ */
