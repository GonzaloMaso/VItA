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
	AbstractObjectCCOTree *tree;
	VTKObjectTreeSplinesNodalWriter *nodalWriter;
	string path;
	string prefix;
public:
	/**
	 * Saves the @p tree in the given path.
	 * @param path	Output directory.
	 * @param prefix	Prefix for the generated files.
	 * @param tree	Tree to be stored in the generated files.
	 */
	VisualizationSavingTask(string path, string prefix, AbstractObjectCCOTree *tree);
	virtual ~VisualizationSavingTask();

	/**
	 * Saves the tree in a visualization file using .vtk format.
	 */
	void execute(long long int terminals);

};

#endif /* VISUALIZATIONSAVINGTASK_H_ */
