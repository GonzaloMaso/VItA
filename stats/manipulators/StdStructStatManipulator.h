/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * StdDiameterStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_STDSTRUCTSTATMANIPULATOR_H_
#define STATS_STDSTRUCTSTATMANIPULATOR_H_

#include "AbstractStructuredStatManipulator.h"
#include "MeanStructStatManipulator.h"

/**
 * Computes the standard deviation value of a specific attribute for an array of vessels.
 */
class StdStructStatManipulator: public AbstractStructuredStatManipulator {
	/**	Mean statistical manipulator used to compute the mean internally. */
	MeanStructStatManipulator *meanManipulator;
public:
	/**
	 * Constructor.
	 */
	StdStructStatManipulator();
	/**
	 * Computes the standard deviation of the attribute of interest among all @p vessels.
	 * @param vessels Array of vessels from which the statistical function is computed.
	 * @return Standard deviation value.
	 */
	double compute(vector<vessel *> vessels, VesselStructHandler::ATTRIBUTE att);
};

#endif /* STATS_STDSTRUCTSTATMANIPULATOR_H_ */
