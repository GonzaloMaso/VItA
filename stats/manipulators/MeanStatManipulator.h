/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MeanStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STATS_MEANSTATMANIPULATOR_H_
#define STATS_MEANSTATMANIPULATOR_H_

#include "../../structures/vascularElements/SingleVessel.h"
#include "Abstract0DStatManipulator.h"

/**
 * Computes the mean of a specific attribute for an array of vessels.
 */
class MeanStatManipulator: public Abstract0DStatManipulator {
public:
	/**
	 * Dummy constructor.
	 */
	MeanStatManipulator();
	/**
	 * Computes the mean of the attribute of interest among all @p vessels.
	 * @param vessels Array of vessels from which the statistical function is computed.
	 * @return Mean value.
	 */
	double compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att);
};

#endif /* STATS_MEANSTATMANIPULATOR_H_ */
