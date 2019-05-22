/*
 * AbstractStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_ABSTRACTSTATMANIPULATOR_H_
#define STATS_ABSTRACTSTATMANIPULATOR_H_

#include <vector>

#include "../structures/CCOCommonStructures.h"
#include "VesselHandler.h"

/**
 * Abstract class for the computation of statistical quantities of a specific attribute for an array of vessels.
 */
class AbstractStatManipulator {
protected:
	/**	Handler used to extract the vessel attributes. */
	VesselHandler *handler;
public:
	/**
	 * Constructor.
	 */
	AbstractStatManipulator();
	/**
	 * Destructor.
	 */
	~AbstractStatManipulator();
	/**
	 * Computes the statistical quantity of interest.
	 * @param vessels	Array of vessels processed to obtain the quantity of interest.
	 * @param att		Attribute over which the quantity is computed.
	 * @return	Statistical quantity of interest.
	 */
	virtual double compute(vector<vessel *> vessels, VesselHandler::ATTRIBUTE att) = 0;
};

#endif /* STATS_ABSTRACTSTATMANIPULATOR_H_ */
