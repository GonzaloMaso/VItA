/*
 * CheckpointSavingTask.cpp
 *
 *  Created on: 24/06/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "CheckpointSavingTask.h"

CheckpointSavingTask::CheckpointSavingTask(string path, string prefix, SingleVesselCCOOTree *tree){
	this->tree = tree;
	this->path = path;
	this->prefix = prefix;
}

CheckpointSavingTask::~CheckpointSavingTask(){
}

void CheckpointSavingTask::execute(long long int terminals){
	tree->save(path + "/" + prefix + to_string(terminals) + ".cco");
}
