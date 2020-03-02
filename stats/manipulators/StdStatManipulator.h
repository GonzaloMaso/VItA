/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * StdDiameterStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STATS_STDSTATMANIPULATOR_H_
#define STATS_STDSTATMANIPULATOR_H_

#include "Abstract0DStatManipulator.h"
#include "MeanStatManipulator.h"

/**
 * Computes the standard deviation value of a specific attribute for an array of vessels.
 */
class StdStatManipulator: public Abstract0DStatManipulator {
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
	double compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att);
};

#endif /* STATS_STDSTATMANIPULATOR_H_ */
