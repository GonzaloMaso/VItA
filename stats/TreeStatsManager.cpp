/*
 * TreeStatsManager.cpp
 *
 *  Created on: Feb 9, 2018
 *      Author: gonzalo
 */

#include "TreeStatsManager.h"

#include "TreeIndexCreator.h"
#include "MeanStatManipulator.h"
#include "StdStatManipulator.h"

TreeStatsManager::TreeStatsManager(AbstractStructuredCCOTree* tree){
	this->tree = tree;
}

void TreeStatsManager::getMeanPerLevel(vector<double> *levels, vector<double> *means, vector<double> *stds, VesselHandler::ATTRIBUTE att){

	TreeIndexCreator *index = new TreeIndexCreator(tree);
	vector<vector<vessel *> > levelSegments = index->getSegmentsByLevel();

	int level = 0;
	MeanStatManipulator *meanManipulator = new MeanStatManipulator();
	StdStatManipulator *stdManipulator = new StdStatManipulator();
	for (vector<vector<vessel *> >::iterator it = levelSegments.begin(); it < levelSegments.end(); ++it, ++level){
		means->push_back(meanManipulator->compute(*it,att));
		stds->push_back(stdManipulator->compute(*it,att));
		levels->push_back(level);
	}

}
