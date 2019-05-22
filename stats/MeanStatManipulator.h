/*
 * MeanStatManipulator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_MEANSTATMANIPULATOR_H_
#define STATS_MEANSTATMANIPULATOR_H_

#include "AbstractStatManipulator.h"

/**
 * Computes the mean of a specific attribute for an array of vessels.
 */
class MeanStatManipulator: public AbstractStatManipulator {
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
	double compute(vector<vessel *> vessels, VesselHandler::ATTRIBUTE att);
};

#endif /* STATS_MEANDIAMETERSTATMANIPULATOR_H_ */
