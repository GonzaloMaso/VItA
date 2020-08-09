#include<vector>
#include<unordered_map>
#include<cstdio>

#include"TreeMerger.h"

#include"../structures/CCOCommonStructures.h"
#include"../structures/vascularElements/AbstractVascularElement.h"
#include"GeneratorData.h"

string coordToString(point &xProx, point &xDist) {
    double coordArray[6] = {xProx.p[0], xProx.p[1], xProx.p[2],
		xDist.p[0], xDist.p[1], xDist.p[2]};
    char *coordCString = (char *) malloc(6 * sizeof(double));
    memcpy(coordCString, &coordArray, 6 * sizeof(double));
    string coordString(coordCString, (6 * sizeof(double) / sizeof(char)));
    free(coordCString);
    return coordString;
}


void TreeMerger::createMapping(SingleVessel *vessel) {
    if(!vessel) {
        return;
    }
    printf("Added vessel at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
        vessel->xProx.p[0], vessel->xProx.p[1], vessel->xProx.p[2],
        vessel->xDist.p[0], vessel->xDist.p[1], vessel->xDist.p[2]);
    string coordString = vessel->coordToString();
    pair<unordered_map<string, SingleVessel *>::iterator, bool> didInsert = this->stringToPointer->insert(pair<string, SingleVessel *>(coordString, vessel));
    if(!didInsert.second) {
        printf("Failed to add element!\n");
    }
    printf("stringToPoint.size = %lu\n", this->stringToPointer->size());
    for (auto it = vessel->children.begin(); it != vessel->children.end(); ++it) {
        SingleVessel *child = (SingleVessel *) (*it);
        createMapping(child);
    }
}

TreeMerger::TreeMerger(string baseTree, vector<string>& derivedTreePoints, GeneratorData *instanceData, AbstractConstraintFunction<double, int> *gam,
    AbstractConstraintFunction<double, int> *epsLim, AbstractConstraintFunction<double, int> *nu) {
        
        this->tree = new SingleVesselCCOOTree(baseTree, instanceData, gam, epsLim, nu);
        
        this->vesselToMerge = new vector<vector<readData> *>;
        
        for (auto it = derivedTreePoints.begin(); it != derivedTreePoints.end(); ++it) {
            
            FILE *fp = fopen((*it).c_str(), "rb");
            if (!fp) {
                fprintf(stderr, "Failed to open derived tree file!\n");
                exit(EXIT_FAILURE);
            }
           
            vector<readData> *toMerge = new vector<readData>;
            this->vesselToMerge->push_back(toMerge);
            
            readData readLine;
            double tempPoints[12];
            int tempFunction;
            fread(&tempPoints, sizeof(double), 12, fp);
            fread(&tempFunction, sizeof(int), 1, fp);
            readLine.xBif.p[0] = tempPoints[0];
            readLine.xBif.p[1] = tempPoints[1];
            readLine.xBif.p[2] = tempPoints[2];
            readLine.xNew.p[0] = tempPoints[3];
            readLine.xNew.p[1] = tempPoints[4];
            readLine.xNew.p[2] = tempPoints[4];
            readLine.xPProx.p[0] = tempPoints[6];
            readLine.xPProx.p[1] = tempPoints[7];
            readLine.xPProx.p[2] = tempPoints[8];
            readLine.xPDist.p[0] = tempPoints[9];
            readLine.xPDist.p[1] = tempPoints[10];
            readLine.xPDist.p[2] = tempPoints[11];            
            while (!feof(fp)) {
                toMerge->push_back(readLine);
                fread(&tempPoints, sizeof(double), 12, fp);
                fread(&tempFunction, sizeof(int), 1, fp);
                readLine.xBif.p[0] = tempPoints[0];
                readLine.xBif.p[1] = tempPoints[1];
                readLine.xBif.p[2] = tempPoints[2];
                readLine.xNew.p[0] = tempPoints[3];
                readLine.xNew.p[1] = tempPoints[4];
                readLine.xNew.p[2] = tempPoints[4];
                readLine.xPProx.p[0] = tempPoints[6];
                readLine.xPProx.p[1] = tempPoints[7];
                readLine.xPProx.p[2] = tempPoints[8];
                readLine.xPDist.p[0] = tempPoints[9];
                readLine.xPDist.p[1] = tempPoints[10];
                readLine.xPDist.p[2] = tempPoints[11];
            }
        }

        // for (auto it1 = vesselToMerge->begin(); it1 != vesselToMerge->end(); ++it1) {
        //     for (auto it2 = (*it1)->begin(); it2 != (*it1)->end(); ++it2) {
        //         printf("xBif = (%.16e %.16e %.16e)\txNew = (%.16e %.16e %.16e)\txProx = (%.16e %.16e %.16e)\t xDist = (%.16e %.16e %.16e)\n",
        //             (*it2).xBif.p[0], (*it2).xBif.p[1], (*it2).xBif.p[2],
        //             (*it2).xNew.p[0], (*it2).xNew.p[1], (*it2).xNew.p[2],
        //             (*it2).xPProx.p[0], (*it2).xPProx.p[1], (*it2).xPProx.p[2],
        //             (*it2).xPDist.p[0], (*it2).xPDist.p[1],(*it2).xPDist.p[2]);
        //     }
        // }

        this->stringToPointer = new unordered_map<string, SingleVessel *>();
        this->createMapping((SingleVessel *) this->tree->getRoot());

        // printf("The mapping:\n");
        // for (auto it = this->stringToPointer->begin(); it != this->stringToPointer->end(); ++it) {
        //     printf("key = %s\n, pointer = %p\n", (*it).first.c_str(), (*it).second);
        // }
}

TreeMerger::~TreeMerger() {
    delete this->stringToPointer;
    for (auto it = this->vesselToMerge->begin(); it != this->vesselToMerge->end(); ++it) {
        delete (*it);
    }
    delete this->vesselToMerge;
    delete tree;
}

SingleVesselCCOOTree* TreeMerger::merge() {
    for (auto itTree = vesselToMerge->begin(); itTree != vesselToMerge->end(); ++itTree) {
        for (auto itVessels = (*itTree)->begin(); itVessels != (*itTree)->end(); ++itVessels) {
            printf("Trying to access parent at xProx = (%.16e, %.16e, %.16e) xDist = (%.16e, %.16e, %.16e)\n",
                (*itVessels).xPProx.p[0], (*itVessels).xPProx.p[1], (*itVessels).xPProx.p[2],
                (*itVessels).xPDist.p[0], (*itVessels).xPDist.p[1], (*itVessels).xPDist.p[2]);
            // printf("key_accessed = %s\n", coordToString((*itVessels).xPProx, (*itVessels).xPDist).c_str());
            SingleVessel *parentPointer = this->stringToPointer->at(coordToString((*itVessels).xPProx, (*itVessels).xPDist));
            tree->addVesselMergeFast((*itVessels).xBif, (*itVessels).xNew, parentPointer, (*itVessels).function, this->stringToPointer);
        }
    }

    // Update tree
    this->tree->updateTree(((SingleVessel *) this->tree->getRoot()), this->tree);
	
	double maxVariation = INFINITY;
	while (maxVariation > this->tree->variationTolerance) {
			this->tree->updateTreeViscositiesBeta(((SingleVessel *) this->tree->getRoot()), &maxVariation);
	}

    return this->tree;
}