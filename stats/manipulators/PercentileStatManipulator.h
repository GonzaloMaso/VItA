/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MeanStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STATS_PERCENTILESTATMANIPULATOR_H_
#define STATS_PERCENTILESTATMANIPULATOR_H_

#include "Abstract0DStatManipulator.h"

/**
 * Computes the N-th percentile value of a specific attribute for an array of vessels.
 */
class PercentileStatManipulator : public Abstract0DStatManipulator{
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
	double compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att);
};

#endif /* STATS_PERCENTILESTATMANIPULATOR_H_ */
