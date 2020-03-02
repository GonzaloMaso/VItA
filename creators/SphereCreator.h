/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SphereCreator.h
 *
 *  Created on: 5 de fev de 2018
 *      Author: gonzalo
 */

#ifndef CREATORS_SPHERECREATOR_H_
#define CREATORS_SPHERECREATOR_H_

#include "AbstractCreator.h"

#include <vector>

using namespace std;

/**
 * Creates a VTP file containing a sphere.
 */
class SphereCreator: public AbstractCreator {
	/** 3D point denoting the center of the sphere. */
	vector<double> center;
	/** Radius of the sphere. */
	double radius;
	/** Azimutal resolution of the sphere. */
	double phiResolution;
	/**	Zenital resolution of the sphere. */
	double thetaResolution;
public:
	/**
	 * Initialize the inner variables of the object.
	 * @param center 3D point denoting the center of the sphere.
	 * @param radius Radius of the sphere.
	 * @param phiResolution Azimutal resolution of the sphere.
	 * @param thetaResolution Zenital resolution of the sphere.
	 */
	SphereCreator(vector<double> center, double radius, double phiResolution, double thetaResolution);
	/**
	 * Standard destructor.
	 */
	virtual ~SphereCreator();
	/**
	 * Persist the object in a file with VTP format.
	 * @param filename	Save file for the current object.
	 */
	void create(string filename);
};

#endif /* CREATORS_SPHERECREATOR_H_ */
