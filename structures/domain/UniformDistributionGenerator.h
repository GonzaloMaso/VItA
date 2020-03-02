/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * UniformDistributionGenerator.h
 *
 *  Created on: 22/02/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_UNIFORMDISTRIBUTIONGENERATOR_H_
#define DOMAIN_UNIFORMDISTRIBUTIONGENERATOR_H_

#include "DistributionGenerator.h"

class UniformDistributionGenerator: public DistributionGenerator {
	/**	Generator for X component */
	uniform_real_distribution<double> distX;
	/**	Generator for Y component */
	uniform_real_distribution<double> distY;
	/**	Generator for Z component */
	uniform_real_distribution<double> distZ;
public:
	/**
	 * Constructor.
	 */
	UniformDistributionGenerator();
	/**
	 * Destructor.
	 */
	virtual ~UniformDistributionGenerator();
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

#endif /* DOMAIN_UNIFORMDISTRIBUTIONGENERATOR_H_ */
