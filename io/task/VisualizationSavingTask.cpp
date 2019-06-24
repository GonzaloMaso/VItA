/*
 * VisualizationSavingTask.cpp
 *
 *  Created on: 24/06/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "VisualizationSavingTask.h"

VisualizationSavingTask::VisualizationSavingTask(string path, string prefix, AbstractObjectCCOTree *tree){
	this->tree = tree;
	this->nodalWriter = new VTKObjectTreeSplinesNodalWriter();
	this->path = path;
	this->prefix = prefix;
}

VisualizationSavingTask::~VisualizationSavingTask()
{
	delete nodalWriter;
}

void VisualizationSavingTask::execute(long long int terminals){
	nodalWriter->write(path+ "/" + prefix + to_string(terminals) + "_view.vtp",tree);
}
