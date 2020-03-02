/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MeanStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_MEANSTRUCTSTATMANIPULATOR_H_
#define STATS_MEANSTRUCTSTATMANIPULATOR_H_

#include "AbstractStructuredStatManipulator.h"

/**
 * Computes the mean of a specific attribute for an array of vessels.
 */
class MeanStructStatManipulator: public AbstractStructuredStatManipulator {
public:
	/**
	 * Dummy constructor.
	 */
	MeanStructStatManipulator();
	/**
	 * Computes the mean of the attribute of interest among all @p vessels.
	 * @param vessels Array of vessels from which the statistical function is computed.
	 * @return Mean value.
	 */
	double compute(vector<vessel *> vessels, VesselStructHandler::ATTRIBUTE att);
};

#endif /* STATS_MEANDIAMETERSTATMANIPULATOR_H_ */
