/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_ABSTRACTSTRUCTUREDSTATMANIPULATOR_H_
#define STATS_ABSTRACTSTRUCTUREDSTATMANIPULATOR_H_

#include <vector>

#include "../../structures/CCOCommonStructures.h"
#include "../VesselStructHandler.h"

/**
 * Abstract class for the computation of statistical quantities of a specific attribute for an array of vessels.
 */
class AbstractStructuredStatManipulator {
protected:
	/**	Handler used to extract the vessel attributes. */
	VesselStructHandler *handler;
public:
	/**
	 * Constructor.
	 */
	AbstractStructuredStatManipulator();
	/**
	 * Destructor.
	 */
	~AbstractStructuredStatManipulator();
	/**
	 * Computes the statistical quantity of interest.
	 * @param vessels	Array of vessels processed to obtain the quantity of interest.
	 * @param att		Attribute over which the quantity is computed.
	 * @return	Statistical quantity of interest.
	 */
	virtual double compute(vector<vessel *> vessels, VesselStructHandler::ATTRIBUTE att) = 0;
};

#endif /* STATS_ABSTRACTSTRUCTUREDSTATMANIPULATOR_H_ */
