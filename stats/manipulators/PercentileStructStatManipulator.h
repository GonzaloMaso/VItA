/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MeanStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_PERCENTILESTRUCTSTATMANIPULATOR_H_
#define STATS_PERCENTILESTRUCTSTATMANIPULATOR_H_

#include "AbstractStructuredStatManipulator.h"

/**
 * Computes the N-th percentile value of a specific attribute for an array of vessels.
 */
class PercentileStructStatManipulator : public AbstractStructuredStatManipulator{
	/**	Number of the percentile which the object computes. */
	int numPercentile;
public:
	/**
	 * Constructor for the statistical manipulator
	 * @param numPercentile	Number of the percentile which the object will compute.
	 */
	PercentileStructStatManipulator(int numPercentile);
	/**
	 * Computes the statistical value of interest.
	 * @param vessels	Array of vessels from which the statistical function is computed.
	 * @return	Value of the percentile.
	 */
	double compute(vector<vessel *> vessels, VesselStructHandler::ATTRIBUTE att);
};

#endif /* STATS_PERCENTILESTRUCTSTATMANIPULATOR_H_ */
