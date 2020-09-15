/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * NormalDistributionGenerator.cpp
 *
 *  Created on: 22/02/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "NormalDistributionGenerator.h"

NormalDistributionGenerator::NormalDistributionGenerator(vector<double> mean, vector<double> std) : DistributionGenerator(){
	this->mean = mean;
	this->std = std;
}

NormalDistributionGenerator::~NormalDistributionGenerator(){
	this->mean.clear();
	this->std.clear();
}

void NormalDistributionGenerator::initialize(int seed, double *boundingBox) {
	DistributionGenerator::initialize(seed, boundingBox);
	distX = normal_distribution<double>(mean[0], std[0]);
	distY = normal_distribution<double>(mean[1], std[1]);
	distZ = normal_distribution<double>(mean[2], std[2]);
}

vector<point> NormalDistributionGenerator::getNPoints(int n){
	vector<point> randomInnerPoints;
	for (int i = 0; i < n; ++i) {
		point p = { distX(generator), distY(generator), distZ(generator) };
		randomInnerPoints.push_back(p);
	}

	return randomInnerPoints;
}
