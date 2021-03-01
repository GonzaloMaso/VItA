#ifndef BREADTH_FIRST_PRUNING_
#define BREADH_FIRT_PRUNING_H_

#include"../SingleVesselCCOOTree.h"
#include"AbstractPruningRule.h"

//	TODO The class is not following the abstract prunnig class.
//	TODO Two implementations in one class, Split this class into two.
class BreadthFirstPruning {
    public:
        SingleVesselCCOOTree *pruneTree(SingleVesselCCOOTree *tree, vector<AbstractPruningRule *>& rules);
        SingleVesselCCOOTree *pruneTreeFast(SingleVesselCCOOTree *tree, vector<AbstractPruningRule *>& rules);
};

#endif
