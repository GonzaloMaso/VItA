#ifndef BREADTH_FIRST_PRUNING_
#define BREADH_FIRT_PRUNING_H_

#include"../SingleVesselCCOOTree.h"
#include"AbstractPruningRule.h"

class BreadthFirstPruning {
    public:
        SingleVesselCCOOTree *pruneTree(SingleVesselCCOOTree *tree, vector<AbstractPruningRule *>& rules);
        SingleVesselCCOOTree *pruneTreeFast(SingleVesselCCOOTree *tree, vector<AbstractPruningRule *>& rules);
};

#endif