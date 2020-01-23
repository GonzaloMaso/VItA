/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeStatsManager.cpp
 *
 *  Created on: Feb 9, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "StructuredTreeStatsManager.h"

#include "manipulators/MeanStructStatManipulator.h"
#include "manipulators/StdStructStatManipulator.h"
#include "StructuredTreeIndexCreator.h"

StructuredTreeStatsManager::StructuredTreeStatsManager(AbstractStructuredCCOTree* tree){
	this->tree = tree;
}

void StructuredTreeStatsManager::getMeanPerLevel(vector<double> *levels, vector<double> *means, vector<double> *stds, VesselStructHandler::ATTRIBUTE att){

	StructuredTreeIndexCreator *index = new StructuredTreeIndexCreator(tree);
	vector<vector<vessel *> > levelSegments = index->getSegmentsByLevel();

	int level = 0;
	MeanStructStatManipulator *meanManipulator = new MeanStructStatManipulator();
	StdStructStatManipulator *stdManipulator = new StdStructStatManipulator();
	for (vector<vector<vessel *> >::iterator it = levelSegments.begin(); it < levelSegments.end(); ++it, ++level){
		means->push_back(meanManipulator->compute(*it,att));
		stds->push_back(stdManipulator->compute(*it,att));
		levels->push_back(level);
	}

}
