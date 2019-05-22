/*
 * MeanStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_PERCENTILESTATMANIPULATOR_H_
#define STATS_PERCENTILESTATMANIPULATOR_H_

#include "AbstractStatManipulator.h"

/**
 * Computes the N-th percentile value of a specific attribute for an array of vessels.
 */
class PercentileStatManipulator: public AbstractStatManipulator {
	/**	Number of the percentile which the object computes. */
	int numPercentile;
public:
	/**
	 * Constructor for the statistical manipulator
	 * @param numPercentile	Number of the percentile which the object will compute.
	 */
	PercentileStatManipulator(int numPercentile);
	/**
	 * Computes the statistical value of interest.
	 * @param vessels	Array of vessels from which the statistical function is computed.
	 * @return	Value of the percentile.
	 */
	double compute(vector<vessel *> vessels, VesselHandler::ATTRIBUTE att);
};

#endif /* STATS_PERCENTILESTATMANIPULATOR_H_ */
