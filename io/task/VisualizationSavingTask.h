/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * VisualizationSavingTask.h
 *
 *  Created on: 24/06/2019
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef VISUALIZATIONSAVINGTASK_H_
#define VISUALIZATIONSAVINGTASK_H_

#include "AbstractSavingTask.h"
#include "../VTKObjectTreeSplinesNodalWriter.h"

/**
 * Saves the tree in .vtk format using splines to smooth the vessels.
 */
class VisualizationSavingTask: public AbstractSavingTask {
	VTKObjectTreeSplinesNodalWriter *nodalWriter;
	string path;
	string prefix;
public:
	/**
	 * Saves the @p tree in the given path.
	 * @param path	Output directory.
	 * @param prefix	Prefix for the generated files.
	 */
	VisualizationSavingTask(string path, string prefix);
	virtual ~VisualizationSavingTask();

	/**
	 * Saves the tree in a visualization file using .vtk format.
	 */
	void execute(long long int terminals, AbstractObjectCCOTree *tree);

};

#endif /* VISUALIZATIONSAVINGTASK_H_ */
