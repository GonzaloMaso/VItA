#include"BreadthFirstPruning.h"

#include<queue>

#include"../SingleVesselCCOOTree.h"
#include"../../vascularElements/SingleVessel.h"
#include"AbstractPruningRule.h"

using namespace std;

bool isVesselToBePruned(SingleVessel *vessel, vector<AbstractPruningRule *>& rules) {
    for (auto it = rules.begin(); it != rules.end(); ++it) {
        bool flag = (*it)->needsPruning(vessel);
        if (flag) {
            return true;
        }
    }
    return false;
}

queue<SingleVessel *>* vesselsToPreserve(SingleVesselCCOOTree *tree, vector<AbstractPruningRule *>& rules) {
    queue<SingleVessel *> *toPreserve = new queue<SingleVessel *>;
    queue<SingleVessel *> *searchable = new queue<SingleVessel *>;
    SingleVessel *root = (SingleVessel *) tree->getRoot();
    searchable->push(root);
    while (!searchable->empty()) {
        SingleVessel *curVessel = searchable->front();
        searchable->pop();
        if(!isVesselToBePruned(curVessel, rules)) {
            toPreserve->push(curVessel);
        }
        vector<AbstractVascularElement *> children = curVessel->getChildren();
        for (auto it = children.begin(); it != children.end(); ++it) {
            searchable->push((SingleVessel *) (*it));
        }
    }
    delete searchable;
    return toPreserve;
}

SingleVesselCCOOTree* BreadthFirstPruning::pruneTree(SingleVesselCCOOTree *tree, vector<AbstractPruningRule *>& rules) {
    unordered_map<SingleVessel*, SingleVessel*> copiedTo;
    // This represents the root parent
    copiedTo[nullptr] = nullptr;
    queue<SingleVessel *>* toPreserve = vesselsToPreserve(tree, rules);
    SingleVesselCCOOTree *newTree = new SingleVesselCCOOTree(tree);
    while(!toPreserve->empty()) {
        SingleVessel *vesselToCopy = toPreserve->front();
        toPreserve->pop();
        SingleVessel *vesselToAdd = new SingleVessel();
        copiedTo[vesselToCopy] = vesselToAdd;
        newTree->addValitatedVessel(vesselToAdd, vesselToCopy, copiedTo);
    }
    delete toPreserve;
    return newTree;
}

SingleVesselCCOOTree* BreadthFirstPruning::pruneTreeFast(SingleVesselCCOOTree *tree, vector<AbstractPruningRule *>& rules) {
    unordered_map<SingleVessel*, SingleVessel*> copiedTo;
    // This represents the root parent
    copiedTo[nullptr] = nullptr;
    queue<SingleVessel *>* toPreserve = vesselsToPreserve(tree, rules);
    SingleVesselCCOOTree *newTree = new SingleVesselCCOOTree(tree);
    while(!toPreserve->empty()) {
        SingleVessel *vesselToCopy = toPreserve->front();
        toPreserve->pop();
        SingleVessel *vesselToAdd = new SingleVessel();
        copiedTo[vesselToCopy] = vesselToAdd;
        newTree->addValitatedVesselFast(vesselToAdd, vesselToCopy, copiedTo);
    }

   	//	Update post-order nLevel, flux, pressure and determine initial resistance and beta values.
	newTree->updateTree(((SingleVessel *) newTree->root), newTree);

    //	Update resistance, pressure and betas
	double maxVariation = INFINITY;
	while (maxVariation > newTree->variationTolerance) {
	    newTree->updateTreeViscositiesBeta(((SingleVessel *) newTree->root), &maxVariation);
	}    

    delete toPreserve;
    return newTree;
}