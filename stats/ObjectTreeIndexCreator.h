/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * TreeIndexCreator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_OBJECTTREEINDEXCREATOR_H_
#define STATS_OBJECTTREEINDEXCREATOR_H_

#include <vector>

#include "../structures/tree/AbstractObjectCCOTree.h"
#include "../structures/vascularElements/SingleVessel.h"

/**
 * Creates indexes for a given tree to ease the manipulation of its segments.
 */
class ObjectTreeIndexCreator {
	/** Tree for which the indexes are created.	 */
	AbstractObjectCCOTree *tree;

public:
	/**
	 * Constructor
	 * @param tree Tree to index.
	 */
	ObjectTreeIndexCreator(AbstractObjectCCOTree *tree);

	/**
	 * Returns all terminal segments of @p tree.
	 * @return Vector of terminal segments.
	 */
	vector<SingleVessel *> getTerminals();
	/**
	 * Returns all segments of @p tree.
	 * @return Vector of segments.
	 */
	vector<SingleVessel *> getAllSegments();
	/**
	 * Returns a vector for each tree bifurcation level, containing the segments of @p tree for such level.
	 * @return Vector of vectors of segments.
	 */
	vector<vector<SingleVessel *> > getSegmentsByLevel();

	/**
	 * Returns a vector for each tree bifurcation level, containing the segments of @p tree for such level.
	 * @param root Tree of the sub/tree to analyse.
	 * @return Vector of vectors of segments.
	 */
	vector<vector<SingleVessel *> > getSegmentsByLevel(SingleVessel *root);

	/**
	 * Returns a vector of all root for the given stage.
	 * @param root Root of the main tree.
	 * @param stage Stage to which roots and all children vessels belong.
	 * @return Vector with roots for the given stage.
	 */
	vector<SingleVessel *> getRootAtStage(SingleVessel *root, int stage);

private:
	void extractTerminals(SingleVessel *root, vector<SingleVessel *> *extractedTerminals);
	void extractAllSegments(SingleVessel *root, vector<SingleVessel *> *extractedTerminals);
	void extractSegmentsByLevel(SingleVessel *root, vector<vector<SingleVessel *> > *extractedTerminals, unsigned level);
	void extractRootsAtStage(SingleVessel *root, vector<SingleVessel *> *extractedRoots, int stage);
};

#endif /* STATS_OBJECTTREEINDEXCREATOR_H_ */
