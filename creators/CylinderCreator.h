/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CylinderCreator.h
 *
 *  Created on: 5 de fev de 2018
 *      Author: gonzalo
 */

#ifndef CREATORS_CYLINDERCREATOR_H_
#define CREATORS_CYLINDERCREATOR_H_

#include "AbstractCreator.h"

#include <vector>

using namespace std;

/**
 * Creates a VTP file containing a cylinder.
 */
class CylinderCreator: public AbstractCreator {
	/** 3D point denoting the center of the sphere. */
	vector<double> center;
	/** Radius of the cylinder base. */
	double radius;
	/**	Height of the cylinder. */
	double height;
	/** Circumferential resolution. */
	int resolution;
public:
	/**
	 * Initialize the inner variables of the object.
	 * @param center 3D point denoting the center of the sphere.
	 * @param radius Radius of the cylinder base.
	 * @param height Height of the cylinder.
	 * @param resolution Circumferential resolution.
	 */
	CylinderCreator(vector<double> center, double radius, double height, int resolution);
	/**
	 * Standard destructor.
	 */
	virtual ~CylinderCreator();
	/**
	 * Persist the object in a file with VTP format.
	 * @param filename	Save file for the current object.
	 */
	void create(string filename);
};

#endif /* CREATORS_CYLINDERCREATOR_H_ */
