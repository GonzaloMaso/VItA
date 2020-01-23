/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * CheckpointSavingTask.cpp
 *
 *  Created on: 24/06/2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "CheckpointSavingTask.h"

CheckpointSavingTask::CheckpointSavingTask(string path, string prefix){
	this->path = path;
	this->prefix = prefix;
}

CheckpointSavingTask::~CheckpointSavingTask(){
}

void CheckpointSavingTask::execute(long long int terminals, AbstractObjectCCOTree *tree){
	tree->save(path + "/" + prefix + to_string(terminals) + ".cco");
}
