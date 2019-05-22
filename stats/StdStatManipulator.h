/*
 * StdDiameterStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_STDSTATMANIPULATOR_H_
#define STATS_STDSTATMANIPULATOR_H_

#include "AbstractStatManipulator.h"
#include "MeanStatManipulator.h"

/**
 * Computes the standard deviation value of a specific attribute for an array of vessels.
 */
class StdStatManipulator: public AbstractStatManipulator {
	/**	Mean statistical manipulator used to compute the mean internally. */
	MeanStatManipulator *meanManipulator;
public:
	/**
	 * Constructor.
	 */
	StdStatManipulator();
	/**
	 * Computes the standard deviation of the attribute of interest among all @p vessels.
	 * @param vessels Array of vessels from which the statistical function is computed.
	 * @return Standard deviation value.
	 */
	double compute(vector<vessel *> vessels, VesselHandler::ATTRIBUTE att);
};

#endif /* STATS_STDSTATMANIPULATOR_H_ */
