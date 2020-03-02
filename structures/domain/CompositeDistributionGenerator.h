/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CompositeDistributionGenerator.h
 *
 *  Created on: 22/02/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_COMPOSITEDISTRIBUTIONGENERATOR_H_
#define DOMAIN_COMPOSITEDISTRIBUTIONGENERATOR_H_

#include "DistributionGenerator.h"

/**
 * Generate points as the union of several generators.
 */
class CompositeDistributionGenerator: public DistributionGenerator {
	/**	Generators */
	vector<DistributionGenerator *> distributions;
public:
	/**
	 * Constructor.
	 * @param distributions
	 */
	CompositeDistributionGenerator(vector<DistributionGenerator *> distributions);
	/**
	 * Destructor.
	 */
	virtual ~CompositeDistributionGenerator();
	/**
	 * Initializer of the generator. Execute it before any getter call is performed.
	 */
	virtual void initialize(int seed, double *boundingBox);
	/**
	 * Return a vector of @p n points of the distribution.
	 * @param n Amount of output points.
	 * @return Vector of distribution points.
	 */
	vector<point> getNPoints(int n);
};

#endif /* DOMAIN_COMPOSITEDISTRIBUTIONGENERATOR_H_ */
