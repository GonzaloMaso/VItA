/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CheckpointSavingTask.h
 *
 *  Created on: 24/06/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef CHECKPOINTSAVINGTASK_H_
#define CHECKPOINTSAVINGTASK_H_

#include "AbstractSavingTask.h"

/**
 * Saves a SingleVesselCCOOTree structure in .cco format. Does not save dLim value, thus resume a generation process may differ from the original execution.
 */
class CheckpointSavingTask: public AbstractSavingTask {
	string path;
	string prefix;
public:
	/**
	 * Saves the @p tree in the given path.
	 * @param path	Output directory.
	 * @param prefix	Prefix for the generated files.
	 */
	CheckpointSavingTask(string path, string prefix);
	virtual ~CheckpointSavingTask();

	/**
	 * Saves the tree in its current state in a .cco format.
	 */
	void execute(long long int terminals, AbstractObjectCCOTree *tree);
};

#endif /* CHECKPOINTSAVINGTASK_H_ */
