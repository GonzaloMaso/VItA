/*
 * TreeIndexCreator.h
 *
 *  Created on: Feb 7, 2018
 *      Author: gonzalo
 */

#ifndef STATS_TREEINDEXCREATOR_H_
#define STATS_TREEINDEXCREATOR_H_

#include <vector>

#include "../structures/CCOCommonStructures.h"
#include "../structures/tree/AbstractStructuredCCOTree.h"

/**
 * Creates indexes for a given tree to ease the manipulation of its segments.
 */
class TreeIndexCreator {
	/** Tree for which the indexes are created.	 */
	AbstractStructuredCCOTree *tree;

public:
	/**
	 * Constructor
	 * @param tree Tree to index.
	 */
	TreeIndexCreator(AbstractStructuredCCOTree *tree);

	/**
	 * Returns all terminal segments of @p tree.
	 * @return Vector of terminal segments.
	 */
	vector<vessel *> getTerminals();
	/**
	 * Returns all segments of @p tree.
	 * @return Vector of segments.
	 */
	vector<vessel *> getAllSegments();
	/**
	 * Returns a vector for each tree bifurcation level, containing the segments of @p tree for such level.
	 * @return Vector of vectors of segments.
	 */
	vector<vector<vessel *> > getSegmentsByLevel();

private:
	void extractTerminals(vessel *root, vector<vessel *> *extractedTerminals);
	void extractAllSegments(vessel *root, vector<vessel *> *extractedTerminals);
	void extractSegmentsByLevel(vessel *root, vector<vector<vessel *> > *extractedTerminals, unsigned level);
};

#endif /* STATS_TREEINDEXCREATOR_H_ */
