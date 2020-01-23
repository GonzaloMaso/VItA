/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeIndexCreator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "ObjectTreeIndexCreator.h"

/**
 * Creates an instance associated with @p tree .
 * @param tree Instance associated with the index creator.
 */
ObjectTreeIndexCreator::ObjectTreeIndexCreator(AbstractObjectCCOTree *tree){

	this->tree = tree;

}

/**
 * Returns a vector with all the tree terminals.
 * @return vector with all the tree terminals.
 */
vector<SingleVessel*> ObjectTreeIndexCreator::getTerminals(){
	vector<SingleVessel *> extractedTerminals;
	extractTerminals((SingleVessel *) tree->getRoot(), &extractedTerminals);
	return extractedTerminals;
}

/**
 * Returns a vector with all the tree segments.
 * @return vector with all the tree segments.
 */
vector<SingleVessel*> ObjectTreeIndexCreator::getAllSegments(){
	vector<SingleVessel *> extractedSegments;
	extractAllSegments((SingleVessel *) tree->getRoot(), &extractedSegments);
	return extractedSegments;
}

/**
 * Returns a vector for which the i-th entry contains a vector with all segments at the i-th tree level.
 * @return vector for which the i-th entry contains a vector with all segments at the i-th tree level.
 */
vector<vector<SingleVessel*> > ObjectTreeIndexCreator::getSegmentsByLevel(){
	vector<vector<SingleVessel*> > extractedSegments;
	extractSegmentsByLevel((SingleVessel *) tree->getRoot(), &extractedSegments, 0);
	return extractedSegments;
}

/**
 * Returns a vector for which the i-th entry contains a vector with all segments at the i-th tree level.
 * @return vector for which the i-th entry contains a vector with all segments at the i-th tree level.
 */
vector<vector<SingleVessel*> > ObjectTreeIndexCreator::getSegmentsByLevel(SingleVessel *root){
	vector<vector<SingleVessel*> > extractedSegments;
	extractSegmentsByLevel(root, &extractedSegments, 0);
	return extractedSegments;
}

vector<SingleVessel*> ObjectTreeIndexCreator::getRootAtStage(SingleVessel *root, int stage){
	vector<SingleVessel*> extractedRoots;
	extractRootsAtStage(root, &extractedRoots, stage);
	return extractedRoots;
}

/**
 * Recursively adds the tree terminals to the given vector @p extractedTerminals.
 * @param root Tree root.
 * @param extractedTerminals Extracted terminals.
 */
void ObjectTreeIndexCreator::extractTerminals(SingleVessel* root, vector<SingleVessel*> *extractedTerminals){

	if(root->getChildren().empty()){
		extractedTerminals->push_back(root);
	}else{
		for (std::vector<AbstractVascularElement *>::iterator it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
			extractTerminals((SingleVessel *) *it,extractedTerminals);
		}
	}

}

/**
 * Recursively in pre-order, adds all tree segments to the vector @p extractedSegments.
 * @param root Tree root.
 * @param extractedSegments Extracted segments.
 */
void ObjectTreeIndexCreator::extractAllSegments(SingleVessel* root, vector<SingleVessel*>* extractedSegments){
	extractedSegments->push_back(root);
	for (std::vector<AbstractVascularElement *>::iterator it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
		extractAllSegments((SingleVessel *) *it,extractedSegments);
	}
}

/**
 * Recursively in pre-order, adds all tree segments to the vector @p extractedSegments in the entry corresponding to their bifurcation level.
 * @param root Tree root.
 * @param extractedSegments Extracted segments.
 * @param level Current deepness level in the tree.
 */
void ObjectTreeIndexCreator::extractSegmentsByLevel(SingleVessel* root, vector<vector<SingleVessel*> >* extractedSegments, unsigned level){

	if(extractedSegments->size() == level){
		vector<SingleVessel *> levelVessels;
		extractedSegments->push_back(levelVessels);
	}
	int isIntermediateVessel = root->getChildren().size() == 1;

	if(!isIntermediateVessel)
		(*extractedSegments)[level].push_back(root);

	for (std::vector<AbstractVascularElement *>::iterator it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
		if(isIntermediateVessel)
			extractSegmentsByLevel((SingleVessel *) *it,extractedSegments,level);
		else
			extractSegmentsByLevel((SingleVessel *) *it,extractedSegments,level+1);
	}

}

void ObjectTreeIndexCreator::extractRootsAtStage(SingleVessel* root, vector<SingleVessel*> *extractedRoots, int stage){
	if(root->stage == stage){
		extractedRoots->push_back(root);
	}else{
		for (std::vector<AbstractVascularElement *>::iterator it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
			extractRootsAtStage((SingleVessel *) *it, extractedRoots, stage);
		}
	}


}
