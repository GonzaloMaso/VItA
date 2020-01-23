/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeStatsManager.cpp
 *
 *  Created on: Feb 9, 2018
 *      Author: gonzalo
 */

#include "ObjectTreeStatsManager.h"

#include "../structures/vascularElements/SingleVessel.h"

#include "manipulators/MeanStatManipulator.h"
#include "manipulators/StdStatManipulator.h"
#include "manipulators/SeriesStatManipulator.h"

#include "ObjectTreeIndexCreator.h"
#include "VesselObjectHandler.h"

ObjectTreeStatsManager::ObjectTreeStatsManager(AbstractObjectCCOTree* tree){
	this->tree = tree;
}

void ObjectTreeStatsManager::getMeanPerLevel(vector<double> *levels, vector<double> *means, vector<double> *stds, VesselObjectHandler::ATTRIBUTE att){

	ObjectTreeIndexCreator *index = new ObjectTreeIndexCreator(tree);
	vector<vector<SingleVessel *> > levelSegments = index->getSegmentsByLevel();

	int level = 0;
	MeanStatManipulator *meanManipulator = new MeanStatManipulator();
	StdStatManipulator *stdManipulator = new StdStatManipulator();

	for (vector<vector<SingleVessel *> >::iterator it = levelSegments.begin(); it < levelSegments.end(); ++it, ++level){
		means->push_back(meanManipulator->compute(*it,att));
		stds->push_back(stdManipulator->compute(*it,att));
		levels->push_back(level);
	}

}

void ObjectTreeStatsManager::getAttributesPerLevel(vector<double> *levels, vector<vector<double>> *values, vector<VesselObjectHandler::ATTRIBUTE> att){
	ObjectTreeIndexCreator *index = new ObjectTreeIndexCreator(tree);
	vector<vector<SingleVessel *> > levelSegments = index->getSegmentsByLevel();

	int numAttributes = att.size();
	for (int idxAttribute = 0; idxAttribute < numAttributes; ++idxAttribute) {
		vector<double> currentAttribute;
		values->push_back(currentAttribute);
	}

	int level = 0;
	SeriesStatManipulator *seriesManipulator = new SeriesStatManipulator();

	for (vector<vector<SingleVessel *> >::iterator it = levelSegments.begin(); it < levelSegments.end(); ++it, ++level){
		vector<vector<double>> seriesByAttCurrentLevel;
		for (int idxAttribute = 0; idxAttribute < numAttributes; ++idxAttribute) {
			seriesByAttCurrentLevel.push_back(seriesManipulator->compute(*it,att[idxAttribute]));
		}

		for (unsigned int idxSamples = 0; idxSamples < seriesByAttCurrentLevel[0].size(); ++idxSamples) {
			levels->push_back(level);
			for (int idxAttribute = 0; idxAttribute < numAttributes; ++idxAttribute) {
				values->at(idxAttribute).push_back(seriesByAttCurrentLevel[idxAttribute][idxSamples]);
			}
		}
	}
}

void ObjectTreeStatsManager::getAttributesPerLevel(SingleVessel *root, vector<double> *levels, vector<vector<double>> *values, vector<VesselObjectHandler::ATTRIBUTE> att){
	ObjectTreeIndexCreator *index = new ObjectTreeIndexCreator(tree);
	vector<vector<SingleVessel *> > levelSegments = index->getSegmentsByLevel(root);

	int numAttributes = att.size();
	for (int idxAttribute = 0; idxAttribute < numAttributes; ++idxAttribute) {
		vector<double> currentAttribute;
		values->push_back(currentAttribute);
	}

	int level = 0;
	SeriesStatManipulator *seriesManipulator = new SeriesStatManipulator();

	for (vector<vector<SingleVessel *> >::iterator it = levelSegments.begin(); it < levelSegments.end(); ++it, ++level){
		vector<vector<double>> seriesByAttCurrentLevel;
		for (int idxAttribute = 0; idxAttribute < numAttributes; ++idxAttribute) {
			seriesByAttCurrentLevel.push_back(seriesManipulator->compute(*it,att[idxAttribute]));
		}

		for (unsigned int idxSamples = 0; idxSamples < seriesByAttCurrentLevel[0].size(); ++idxSamples) {
			levels->push_back(level);
			for (int idxAttribute = 0; idxAttribute < numAttributes; ++idxAttribute) {
				values->at(idxAttribute).push_back(seriesByAttCurrentLevel[idxAttribute][idxSamples]);
			}
		}
	}
}

void ObjectTreeStatsManager::getBranchesAttributesPerLevel(vector<double>* branches, vector<double>* levels, vector<vector<double> >* values,
		vector<VesselObjectHandler::ATTRIBUTE> att, int stage){

	ObjectTreeIndexCreator *index = new ObjectTreeIndexCreator(tree);
	vector<SingleVessel *> rootSegments = index->getRootAtStage((SingleVessel *)tree->getRoot(),stage);

	for (unsigned int idxAtt = 0; idxAtt < att.size(); ++idxAtt) {
		vector<double> currentAttribute;
		values->push_back(currentAttribute);
	}

	int idxBranch = 0;
	for (std::vector<SingleVessel *>::iterator itRoots = rootSegments.begin(); itRoots != rootSegments.end(); ++itRoots, ++idxBranch) {
		vector<double> currentLevels;
		vector<vector<double>> currentValues;
		getAttributesPerLevel(*itRoots, &currentLevels, &currentValues, att);

		int numSamples = currentLevels.size();
		for (int idxSample = 0; idxSample < numSamples; ++idxSample) {
			branches->push_back(idxBranch);
		}
		levels->insert(levels->end(),currentLevels.begin(),currentLevels.end());
		for (unsigned int idxAtt = 0; idxAtt < att.size(); ++idxAtt) {
			values->at(idxAtt).insert(values->at(idxAtt).end(),currentValues[idxAtt].begin(),currentValues[idxAtt].end());
		}
	}

}
