/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * ParallelepipedCreator.h
 *
 *  Created on: 5 de fev de 2018
 *      Author: gonzalo
 */

#ifndef CREATORS_PARALLELEPIPEDCREATOR_H_
#define CREATORS_PARALLELEPIPEDCREATOR_H_

#include "AbstractCreator.h"

/**
 * Creates a VTP file containing a parallelepiped.
 */
class ParallelepipedCreator : public AbstractCreator{
	/** Vertex 1 */
	double *lb;
	/** Vertex 2 */
	double *ub;
public:
	/**
	 * Initialize the inner variables of the object. Vertexes @p lb
	 * and @p ub must be opposite vertexes of the parallelepiped.
	 * @param lb Vertex v1.
	 * @param ub Vertex v2.
	 */
	ParallelepipedCreator(double *lb, double *ub);
	/**
	 * Standard destructor.
	 */
	virtual ~ParallelepipedCreator();
	/**
	 * Persist the object in a file with VTP format.
	 * @param filename	Save file for the current object.
	 */
	void create(string filename);
};

#endif /* CREATORS_PARALLELEPIPEDCREATOR_H_ */
