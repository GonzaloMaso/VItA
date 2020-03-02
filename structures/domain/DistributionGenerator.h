/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * DistributionGenerator.h
 *
 *  Created on: 22/02/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_DISTRIBUTIONGENERATOR_H_
#define DOMAIN_DISTRIBUTIONGENERATOR_H_

#include <vector>
#include <random>
#include "../CCOCommonStructures.h"

/**
 * Abstract generator of points.
 */
class DistributionGenerator {
protected:
	/**	Bounding box of the domain where the points will be generated. */
	double *boundingBox;
	/**	Random generator. */
	mt19937 generator;
public:
	/**
	 * Constructor
	 */
	DistributionGenerator();
	/**
	 * Destructor
	 */
	virtual ~DistributionGenerator();
	/**
	 * Initializer of the generator. Execute it before any getter call is performed.
	 */
	virtual void initialize(int seed, double *boundingBox);
	/**
	 * Return a vector of @p n points of the distribution.
	 * @param n Amount of output points.
	 * @return Vector of distribution points.
	 */
	virtual vector<point> getNPoints(int n) = 0;
};

#endif /* DOMAIN_DISTRIBUTIONGENERATOR_H_ */
