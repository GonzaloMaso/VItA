/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * NormalDistributionGenerator.h
 *
 *  Created on: 22/02/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef DOMAIN_NORMALDISTRIBUTIONGENERATOR_H_
#define DOMAIN_NORMALDISTRIBUTIONGENERATOR_H_

#include "DistributionGenerator.h"

class NormalDistributionGenerator: public DistributionGenerator {
	/**	Generator for X component */
	normal_distribution<double> distX;
	/**	Generator for Y component */
	normal_distribution<double> distY;
	/**	Generator for Z component */
	normal_distribution<double> distZ;

	/**	Distribution mean */
	vector<double> mean;
	/**	Distribution std */
	vector<double> std;

public:
	/**
	 * Constructor.
	 * @param mean Three component vector containing the means for each point component.
	 * @param std Three component vector containing the standard deviation for each point component.
	 */
	NormalDistributionGenerator(vector<double> mean, vector<double> std);
	/**
	 * Destructor.
	 */
	virtual ~NormalDistributionGenerator();
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

#endif /* DOMAIN_NORMALDISTRIBUTIONGENERATOR_H_ */
