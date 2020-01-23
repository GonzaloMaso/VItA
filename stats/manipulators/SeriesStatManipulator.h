/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * MeanStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STATS_SERIESSTATMANIPULATOR_H_
#define STATS_SERIESSTATMANIPULATOR_H_

#include "../../structures/vascularElements/SingleVessel.h"
#include "Abstract1DStatManipulator.h"

/**
 * Extracts a specific attribute for an array of vessels.
 */
class SeriesStatManipulator: public Abstract1DStatManipulator {
public:
	/**
	 * Dummy constructor.
	 */
	SeriesStatManipulator();
	/**
	 * Extracts an attribute of interest among for each vessel in @p.
	 * @param vessels Array of vessels.
	 * @return Attribute series.
	 */
	vector<double> compute(vector<SingleVessel *> vessels, VesselObjectHandler::ATTRIBUTE att);
};

#endif /* STATS_SERIESSTATMANIPULATOR_H_ */
