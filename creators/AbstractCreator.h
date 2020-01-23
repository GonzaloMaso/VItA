/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * AbstractCreator.h
 *
 *  Created on: 5 de fev de 2018
 *      Author: gonzalo
 */

#ifndef CREATORS_ABSTRACTCREATOR_H_
#define CREATORS_ABSTRACTCREATOR_H_

#include <iostream>
using namespace std;

/**
 * Abstract class for the creation and persistance of geometric
 * objects in VTP format files.
 */
class AbstractCreator {
public:
	/**
	 * Dummy constructor.
	 */
	AbstractCreator();
	/**
	 * Dummy destructor.
	 */
	virtual ~AbstractCreator();
	/**
	 * Abstract persistent method.
	 * @param filename	Save file for the current object.
	 */
	virtual void create(string filename) = 0;
};

#endif /* CREATORS_ABSTRACTCREATOR_H_ */
