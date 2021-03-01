#ifndef TREE_MERGER_H_
#define TREE_MERGER_H_

#include<vector>
#include<unordered_map>
#include"../structures/tree/SingleVesselCCOOTree.h"
#include"../constrains/AbstractConstraintFunction.h"

using namespace std;

struct ReadData {
    point xBif, xNew, xPProx, xPDist;
    AbstractVascularElement::VESSEL_FUNCTION function;
};

class TreeMerger {
    SingleVesselCCOOTree *tree;
    vector<vector<ReadData>*> *vesselToMerge;
    unordered_map<string, SingleVessel *> *stringToPointer;

    public:
    TreeMerger(string baseTree, vector<string>& derivedTreePoints, GeneratorData *instanceData, AbstractConstraintFunction<double, int> *gam, AbstractConstraintFunction<double, int> *epsLim, AbstractConstraintFunction<double, int> *nu, bool isInCm);
    ~TreeMerger();
    SingleVesselCCOOTree *mergeFast();
    SingleVesselCCOOTree *merge();

    private:
    void createMapping(SingleVessel *vessel);

};

#endif