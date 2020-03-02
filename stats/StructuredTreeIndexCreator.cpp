/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeIndexCreator.cpp
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#include "StructuredTreeIndexCreator.h"

/**
 * Creates an instance associated with @p tree .
 * @param tree Instance associated with the index creator.
 */
StructuredTreeIndexCreator::StructuredTreeIndexCreator(AbstractStructuredCCOTree *tree){

	this->tree = tree;

}

/**
 * Returns a vector with all the tree terminals.
 * @return vector with all the tree terminals.
 */
vector<vessel*> StructuredTreeIndexCreator::getTerminals(){
	vector<vessel *> extractedTerminals;
	extractTerminals(tree->getRoot(), &extractedTerminals);
	return extractedTerminals;
}

/**
 * Returns a vector with all the tree segments.
 * @return vector with all the tree segments.
 */
vector<vessel*> StructuredTreeIndexCreator::getAllSegments(){
	vector<vessel *> extractedSegments;
	extractAllSegments(tree->getRoot(), &extractedSegments);
	return extractedSegments;
}

/**
 * Returns a vector for which the i-th entry contains a vector with all segments at the i-th tree level.
 * @return vector for which the i-th entry contains a vector with all segments at the i-th tree level.
 */
vector<vector<vessel*> > StructuredTreeIndexCreator::getSegmentsByLevel(){
	vector<vector<vessel*> > extractedSegments;
	extractSegmentsByLevel(tree->getRoot(), &extractedSegments, 0);
	return extractedSegments;
}

/**
 * Recursively adds the tree terminals to the given vector @p extractedTerminals.
 * @param root Tree root.
 * @param extractedTerminals Extracted terminals.
 */
void StructuredTreeIndexCreator::extractTerminals(vessel* root, vector<vessel*> *extractedTerminals){

	if(root->anastomose.size()>1){
		extractTerminals(root->anastomose[1],extractedTerminals);
		extractTerminals(root->anastomose[2],extractedTerminals);
	}else{
		extractedTerminals->push_back(root);
	}

}

/**
 * Recursively in pre-order, adds all tree segments to the vector @p extractedSegments.
 * @param root Tree root.
 * @param extractedSegments Extracted segments.
 */
void StructuredTreeIndexCreator::extractAllSegments(vessel* root, vector<vessel*>* extractedSegments){
	if(root->anastomose.size()>1){
		extractedSegments->push_back(root);
		extractAllSegments(root->anastomose[1],extractedSegments);
		extractAllSegments(root->anastomose[2],extractedSegments);
	}else{
		extractedSegments->push_back(root);
	}
}

/**
 * Recursively in pre-order, adds all tree segments to the vector @p extractedSegments in the entry corresponding to their bifurcation level.
 * @param root Tree root.
 * @param extractedSegments Extracted segments.
 * @param level Current deepness level in the tree.
 */
void StructuredTreeIndexCreator::extractSegmentsByLevel(vessel* root, vector<vector<vessel*> >* extractedSegments, unsigned level){

	if(extractedSegments->size() == level){
		vector<vessel *> levelVessels;
		extractedSegments->push_back(levelVessels);
	}

	if(root->anastomose.size()>1){
		(*extractedSegments)[level].push_back(root);
		extractSegmentsByLevel(root->anastomose[1],extractedSegments,level+1);
		extractSegmentsByLevel(root->anastomose[2],extractedSegments,level+1);
	}else{
		(*extractedSegments)[level].push_back(root);
	}

}
