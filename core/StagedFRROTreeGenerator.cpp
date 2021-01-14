/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * StagedFRROTreeGenerator.cpp
 *
 *  Created on: Mar 31, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#include "StagedFRROTreeGenerator.h"

#include <cmath>
#include <fstream>
#include<unordered_set>

#include "../io/VTKObjectTreeElementalWriter.h"
#include "../io/VTKObjectTreeSplinesNodalWriter.h"

#include "../structures/CCOCommonStructures.h"
#include "../structures/tree/SingleVesselCCOOTree.h"
#include "../structures/vascularElements/AbstractVascularElement.h"
#include "../utils/MemoryMonitor.h"

#include <omp.h>

StagedFRROTreeGenerator::StagedFRROTreeGenerator(
		StagedDomain* domain, point xi, double rootRadii, double qi,
		long long nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu,
		double refPressure, double viscosityTolerance) {

	domain->registerObserver(this);
	this->domain = domain;
	this->instanceData = domain->getInstanceData();
	SingleVessel::bifurcationTests = instanceData->nBifurcationTest;
	this->nTerminals = nTerm;
	this->stage = domain->getCurrentStage();

	this->gams = gam;
	this->epsLims = epsLim;
	this->nus = nu;

	this->tree = new SingleVesselCCOOTree(xi, rootRadii, qi, gam[0], epsLim[0],
			nu[0], refPressure, viscosityTolerance, instanceData);
	this->tree->setCurrentStage(stage);
	this->didAllocateTree = true;
	

	this->dLim = domain->getDLim(1, instanceData->perfusionAreaFactor);

	this->isGeneratingConfFile = 0;
	this->confFilename = "";

	this->dataMonitor = new GeneratorDataMonitor(domain);
	this->monitor = new MemoryMonitor(MemoryMonitor::MEGABYTE);
}

StagedFRROTreeGenerator::StagedFRROTreeGenerator(
		StagedDomain* domain, AbstractObjectCCOTree* tree, long long nTerm,
		vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu) {

	domain->registerObserver(this);
	this->domain = domain;
	this->instanceData = domain->getInstanceData();
	SingleVessel::bifurcationTests = instanceData->nBifurcationTest;
	this->nTerminals = nTerm;
	this->tree = tree;
	//	Stage can be loaded from file
	this->stage = domain->getCurrentStage();
	this->tree->setCurrentStage(domain->getCurrentStage());

	this->gams = gam;
	this->epsLims = epsLim;
	this->nus = nu;

	this->tree->setGam(gam[0]);
	this->tree->setEpsLim(epsLim[0]);
	this->tree->setNu(nu[0]);

	this->dLim = domain->getDLim(1, instanceData->perfusionAreaFactor);

	this->isGeneratingConfFile = 0;
	this->confFilename = "";

	this->dataMonitor = new GeneratorDataMonitor(domain);
	this->monitor = new MemoryMonitor(MemoryMonitor::MEGABYTE);

	this->didAllocateTree = false;
}

StagedFRROTreeGenerator::~StagedFRROTreeGenerator() {
	delete this->dataMonitor;
	delete this->monitor;
	if (didAllocateTree) {
		delete this->tree;
	}
}

AbstractObjectCCOTree *StagedFRROTreeGenerator::generate(long long int saveInterval, string tempDirectory) {

	this->beginTime = time(nullptr);
	this->dLimInitial = this->dLim;

//	VTKObjectTreeSplinesNodalWriter *nodalWriter = new VTKObjectTreeSplinesNodalWriter();
	generatesConfigurationFile(ios::out);

	point xNew;
	int i = 0;

	do {
		xNew = domain->getRandomPoint();
	} while (!isValidRootSegment(xNew, ++i));

	tree->addVessel(xNew, xNew, NULL, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction);

	for (long long i = 1; i < nTerminals; i = tree->getNTerms()) {

		dataMonitor->update();

		if (i % saveInterval == 0 ) {
			saveStatus(i);
		}

		int invalidTerminal = true;
		int iTry = 0;
		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
		while (invalidTerminal) {

			do {
				xNew = domain->getRandomPoint();
			} while (!isValidSegment(xNew, ++iTry));

			int nNeighbors;
			vector<AbstractVascularElement *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

			double minCost = INFINITY;
			point minBif;
			AbstractVascularElement *minParent = NULL;
#pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
				point xBif;
				double cost;
				tree->testVessel(xNew, neighborVessels[j], domain,
						neighborVessels, dLim, &xBif, &cost); //	Inf cost stands for invalid solution
#pragma omp critical
				{
					if (cost < minCost) {
						minCost = cost;
						minBif = xBif;
						minParent = neighborVessels[j];
					}
				}
			}
			//	end for trees

			if (minCost < INFINITY) {
				cout << "Added with a cost of " << minCost << endl;
				tree->addVessel(minBif, xNew, minParent, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction);
				invalidTerminal = false;
			}
		}
		dataMonitor->addDLimValue(dLim,i);
		domain->update();
		//	TODO Reset data monitor!!
		//	Terminal added.
	}
	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());

	this->endTime = time(nullptr);
	this->dLimLast = this->dLim;

	saveStatus(nTerminals-1);
	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	return tree;
}

AbstractObjectCCOTree *StagedFRROTreeGenerator::generateExperimental(long long int saveInterval, string tempDirectory, int maxNumOfTrials) {

	this->beginTime = time(nullptr);
	this->dLimInitial = this->dLim;

//	VTKObjectTreeSplinesNodalWriter *nodalWriter = new VTKObjectTreeSplinesNodalWriter();
	generatesConfigurationFile(ios::out);

	point xNew;
	int i = 0;

	do {
		xNew = domain->getRandomPoint();
	} while (!isValidRootSegment(xNew, ++i));

	tree->addVessel(xNew, xNew, NULL, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction);

	for (long long i = 1; i < nTerminals; i = tree->getNTerms()) {

		dataMonitor->update();

		if (i % saveInterval == 0 ) {
			saveStatus(i);
		}

		int invalidTerminal = true;
		int iTry = 0;
		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
		while (invalidTerminal) {

			do {
				xNew = domain->getRandomPoint();
			} while (!isValidSegment(xNew, ++iTry));
			if (iTry > maxNumOfTrials) {
				return NULL;
			}

			int nNeighbors;
			vector<AbstractVascularElement *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

			double minCost = INFINITY;
			point minBif;
			AbstractVascularElement *minParent = NULL;
#pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
				point xBif;
				double cost;
				tree->testVessel(xNew, neighborVessels[j], domain,
						neighborVessels, dLim, &xBif, &cost); //	Inf cost stands for invalid solution
#pragma omp critical
				{
					if (cost < minCost) {
						minCost = cost;
						minBif = xBif;
						minParent = neighborVessels[j];
					}
				}
			}
			//	end for trees

			if (minCost < INFINITY) {
				cout << "Added with a cost of " << minCost << endl;
				tree->addVessel(minBif, xNew, minParent, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction);
				invalidTerminal = false;
			}
		}
		dataMonitor->addDLimValue(dLim,i);
		domain->update();
		//	TODO Reset data monitor!!
		//	Terminal added.
	}
	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());

	this->endTime = time(nullptr);
	this->dLimLast = this->dLim;

	saveStatus(nTerminals-1);
	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	return tree;
}

int StagedFRROTreeGenerator::isValidRootSegment(point xNew,
		int iTry) {

	if (iTry % instanceData->nTerminalTrial == 0) {
		dLim *= instanceData->dLimReductionFactor;
//		cout << "DLim reduced." << endl;
	}
	if( !(domain->isSegmentInside(tree->getXProx(),xNew)) )
		return false;

	point dVect = xNew - tree->getXProx();

	return sqrt(dVect ^ dVect) > dLim;
}

int StagedFRROTreeGenerator::isValidSegment(point xNew, int iTry) {

	if (iTry % instanceData->nTerminalTrial == 0) {
		dLim *= instanceData->dLimReductionFactor;
		cout << "DLim reduced." << endl;
	}
	point xBif;
	double dist;
	tree->getClosestTreePoint(xNew, &xBif, &dist);
//	cout << iTry << ": Closest point to xNew=" << xNew << " is " << xBif << " at " << dist << " (distance limit " << dLim << ")" << endl;

	return dist > dLim;
}

StagedDomain * StagedFRROTreeGenerator::getDomain() {
	return domain;
}

vector<AbstractObjectCCOTree *> StagedFRROTreeGenerator::getTrees() {
	vector<AbstractObjectCCOTree *> trees;
	trees.push_back(tree);
	return trees;
}

void StagedFRROTreeGenerator::enableConfigurationFile(
		string filename) {
	this->isGeneratingConfFile = 1;
	this->confFilename = filename;
}

void StagedFRROTreeGenerator::generatesConfigurationFile(ios::openmode mode) {

	if (isGeneratingConfFile) {
		confFile.open(confFilename.c_str(), mode);
		confFile.setf(ios::scientific, ios::floatfield);
		confFile.precision(16);

		confFile << "*GlobalParameters" << endl;
		confFile << "N_LEVELS_SCALING_TEST " << instanceData->nLevelTest << endl;
		confFile << "N_TRIAL " << instanceData->nTerminalTrial << endl;
		confFile << "DLIM_REDUCTION_FACTOR " << instanceData->dLimReductionFactor << endl;
		confFile << "PERFUSION_AREA_FACTOR " << instanceData->perfusionAreaFactor << endl;
		confFile << "CLOSE_NEIGHBORHOOD_FACTOR " << instanceData->closeNeighborhoodFactor << endl;
		confFile << "MIDPOINT_DLIM_FACTOR " << instanceData->midPointDlimFactor << endl;
		confFile << "N_BIF_TRIES " << instanceData->nBifurcationTest << endl;
		confFile << endl;

		confFile << "*TreeParameters" << endl;
		confFile << "TREE_CLASS " << this->tree->getTreeName() << endl;
		confFile << "N_TERMINALS " << this->nTerminals << endl;
		confFile << "INPUT_POSITION " << this->tree->getXProx() << endl;
		confFile << "INPUT_FLOW " << this->tree->getQProx() << endl;
		confFile << "PRESSURE_DROP " << this->tree->getDp() << endl;
		confFile << endl;

		confFile << "*Timestamps" << endl;
		beginningTime = chrono::steady_clock::now();
	}
}

void StagedFRROTreeGenerator::markTimestampOnConfigurationFile(
		string label) {
	if (isGeneratingConfFile) {
		confFile << (chrono::duration_cast<chrono::microseconds>(
				chrono::steady_clock::now() - beginningTime).count())
				/ 1000000.0 << ": " << label << endl;
	}
}

void StagedFRROTreeGenerator::closeConfigurationFile() {

	confFile << endl << "DOMAIN_POINTS_GENERATED " << this->domain->getPointCounter() << endl;
	if (isGeneratingConfFile) {
		confFile.flush();
		confFile.close();
	}
}

AbstractObjectCCOTree *StagedFRROTreeGenerator::resume(long long int saveInterval, string tempDirectory) {

	this->beginTime = time(nullptr);
	this->dLimInitial = this->dLim;

//	VTKObjectTreeSplinesNodalWriter *nodalWriter = new VTKObjectTreeSplinesNodalWriter();
	generatesConfigurationFile(ios::out);

	for (long long int j = 0; j < tree->getPointCounter(); ++j) {
		domain->getRandomPoint();
	}

	//	Compute current nTerm
	long long currentTerminals = tree->getNTerms();
	point xNew;

	cout << "Generating from " << currentTerminals << " to " << nTerminals << "..." << endl;
	//	Be careful nTerminals may differ from the current amount of terminals since vessel-tip conexions are allowed in some cases.
	for (long long i = currentTerminals; i < nTerminals; i = tree->getNTerms()) {

		dataMonitor->update();

		if (i % saveInterval == 0 || i == currentTerminals) {
			saveStatus(i);
		}

		int invalidTerminal = true;
		int iTry = 0;
		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
		while (invalidTerminal) {

			do {
				xNew = domain->getRandomPoint();
			} while (!isValidSegment(xNew, ++iTry));

			int nNeighbors;
			vector<AbstractVascularElement *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

			double minCost = INFINITY;
			point minBif;
			AbstractVascularElement *minParent = NULL;
#pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
				point xBif;
				double cost;
				tree->testVessel(xNew, neighborVessels[j], domain,
						neighborVessels, dLim, &xBif, &cost); //	Inf cost stands for invalid solution
#pragma omp critical
				{
					if (cost < minCost) {
						minCost = cost;
						minBif = xBif;
						minParent = neighborVessels[j];
					}
				}
			}
			//	end for trees

			if (minCost < INFINITY) {
				cout << "Added with a cost of " << minCost << " with a total cost of " << ((SingleVessel *) tree->getRoot())->treeVolume << endl;
				tree->addVessel(minBif, xNew, minParent, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction);
				invalidTerminal = false;
			}
		}
		dataMonitor->addDLimValue(dLim,i);
		domain->update();
		//	TODO Reset data monitor!!
		//	Terminal added.
	}
	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());

	this->endTime = time(nullptr);
	this->dLimLast = this->dLim;

	saveStatus(nTerminals-1);
	markTimestampOnConfigurationFile("Final tree volume " + to_string(((SingleVessel *) tree->getRoot())->treeVolume));
	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	return tree;

}

AbstractObjectCCOTree *StagedFRROTreeGenerator::resumeExperimental(long long int saveInterval, string tempDirectory, int maxNumOfTrials) {

	this->beginTime = time(nullptr);
	this->dLimInitial = this->dLim;

//	VTKObjectTreeSplinesNodalWriter *nodalWriter = new VTKObjectTreeSplinesNodalWriter();
	generatesConfigurationFile(ios::out);

	for (long long int j = 0; j < tree->getPointCounter(); ++j) {
		domain->getRandomPoint();
	}

	//	Compute current nTerm
	long long currentTerminals = tree->getNTerms();
	point xNew;

	cout << "Generating from " << currentTerminals << " to " << nTerminals << "..." << endl;
	//	Be careful nTerminals may differ from the current amount of terminals since vessel-tip conexions are allowed in some cases.
	for (long long i = currentTerminals; i < nTerminals; i = tree->getNTerms()) {

		dataMonitor->update();

		if (i % saveInterval == 0 || i == currentTerminals) {
			saveStatus(i);
		}

		int invalidTerminal = true;
		int iTry = 0;
		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
		while (invalidTerminal) {

			do {
				xNew = domain->getRandomPoint();
			} while (!isValidSegment(xNew, ++iTry));
			if (iTry > maxNumOfTrials) {
				return NULL;
			}

			int nNeighbors;
			vector<AbstractVascularElement *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

			double minCost = INFINITY;
			point minBif;
			AbstractVascularElement *minParent = NULL;
#pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
				point xBif;
				double cost;
				tree->testVessel(xNew, neighborVessels[j], domain,
						neighborVessels, dLim, &xBif, &cost); //	Inf cost stands for invalid solution
#pragma omp critical
				{
					if (cost < minCost) {
						minCost = cost;
						minBif = xBif;
						minParent = neighborVessels[j];
					}
				}
			}
			//	end for trees

			if (minCost < INFINITY) {
				cout << "Added with a cost of " << minCost << " with a total cost of " << ((SingleVessel *) tree->getRoot())->treeVolume << endl;
				tree->addVessel(minBif, xNew, minParent, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction);
				invalidTerminal = false;
			}
		}
		dataMonitor->addDLimValue(dLim,i);
		domain->update();
		//	TODO Reset data monitor!!
		//	Terminal added.
	}
	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());

	this->endTime = time(nullptr);
	this->dLimLast = this->dLim;

	saveStatus(nTerminals-1);
	markTimestampOnConfigurationFile("Final tree volume " + to_string(((SingleVessel *) tree->getRoot())->treeVolume));
	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	return tree;

}

AbstractObjectCCOTree *StagedFRROTreeGenerator::resumeSavePoints(long long int saveInterval, string tempDirectory, FILE *fp) {
	this->beginTime = time(nullptr);
	this->dLimInitial = this->dLim;

//	VTKObjectTreeSplinesNodalWriter *nodalWriter = new VTKObjectTreeSplinesNodalWriter();
	generatesConfigurationFile(ios::out);

	for (long long int j = 0; j < tree->getPointCounter(); ++j) {
		domain->getRandomPoint();
	}

	//	Compute current nTerm
	long long currentTerminals = tree->getNTerms();
	point xNew;

	cout << "Generating from " << currentTerminals << " to " << nTerminals << "..." << endl;
	//	Be careful nTerminals may differ from the current amount of terminals since vessel-tip conexions are allowed in some cases.
	for (long long i = currentTerminals; i < nTerminals; i = tree->getNTerms()) {

		dataMonitor->update();

		if (i % saveInterval == 0 || i == currentTerminals) {
			saveStatus(i);
		}

		int invalidTerminal = true;
		int iTry = 0;
		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
		while (invalidTerminal) {

			do {
				xNew = domain->getRandomPoint();
			} while (!isValidSegment(xNew, ++iTry));

			int nNeighbors;
			vector<AbstractVascularElement *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

			double minCost = INFINITY;
			point minBif;
			AbstractVascularElement *minParent = NULL;
#pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
				point xBif;
				double cost;
				tree->testVessel(xNew, neighborVessels[j], domain,
						neighborVessels, dLim, &xBif, &cost); //	Inf cost stands for invalid solution
#pragma omp critical
				{
					if (cost < minCost) {
						minCost = cost;
						minBif = xBif;
						minParent = neighborVessels[j];
					}
				}
			}
			//	end for trees

			if (minCost < INFINITY) {
				cout << "Added with a cost of " << minCost << " with a total cost of " << ((SingleVessel *) tree->getRoot())->treeVolume << endl;
				SingleVessel *minParentSV = (SingleVessel *) minParent;
				fwrite(&(minBif.p[0]), sizeof(double), 1, fp);
				fwrite(&(minBif.p[1]), sizeof(double), 1, fp);
				fwrite(&(minBif.p[2]), sizeof(double), 1, fp);
				fwrite(&(xNew.p[0]), sizeof(double), 1, fp);
				fwrite(&(xNew.p[1]), sizeof(double), 1, fp);
				fwrite(&(xNew.p[2]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xProx.p[0]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xProx.p[1]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xProx.p[2]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xDist.p[0]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xDist.p[1]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xDist.p[2]), sizeof(double), 1, fp);
				fwrite(&(instanceData->vesselFunction), sizeof(int), 1, fp);
				tree->addVessel(minBif, xNew, minParent, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction);				
				invalidTerminal = false;
			}
		}
		dataMonitor->addDLimValue(dLim,i);
		domain->update();
		//	TODO Reset data monitor!!
		//	Terminal added.
	}
	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());

	this->endTime = time(nullptr);
	this->dLimLast = this->dLim;

	saveStatus(nTerminals-1);
	markTimestampOnConfigurationFile("Final tree volume " + to_string(((SingleVessel *) tree->getRoot())->treeVolume));
	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	return tree;
}

void originalVesselsRecursive(SingleVessel *root, unordered_map<string, pair<SingleVessel *, bool>>* originals, vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints) {
	if(!root) {
		return;
	}
	point midPoint = (root->xProx + root->xDist) / 2.;
	string key = root->coordToString();
	pair<SingleVessel *, bool> value = pair<SingleVessel *, bool>(root, static_cast<bool>(enclosedPoints->IsInsideSurface(midPoint.p)));
	originals->insert(pair<string, pair<SingleVessel *, bool>>(key, value));
	for(auto it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
		originalVesselsRecursive(static_cast<SingleVessel *>((*it)), originals, enclosedPoints);
	}
}


unordered_map<string, pair<SingleVessel*, bool>>* originalVessels(SingleVesselCCOOTree *tree, vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints) {
	unordered_map<string, pair<SingleVessel *, bool>>* originalSV = new unordered_map<string, pair<SingleVessel *, bool>>;
	originalVesselsRecursive(static_cast<SingleVessel *>(tree->getRoot()), originalSV, enclosedPoints);
	return originalSV;
}

AbstractObjectCCOTree *StagedFRROTreeGenerator::resumeSavePointsMidPoint(long long int saveInterval, string tempDirectory, FILE *fp) {
	this->beginTime = time(nullptr);
	this->dLimInitial = this->dLim;

	vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = domain->getEnclosedPoints();
	unordered_map<string, pair<SingleVessel *, bool>> *ogVessels = originalVessels(static_cast<SingleVesselCCOOTree *>(this->tree), enclosedPoints);

//	VTKObjectTreeSplinesNodalWriter *nodalWriter = new VTKObjectTreeSplinesNodalWriter();
	generatesConfigurationFile(ios::out);

	for (long long int j = 0; j < tree->getPointCounter(); ++j) {
		domain->getRandomPoint();
	}

	//	Compute current nTerm
	long long currentTerminals = tree->getNTerms();
	point xNew;

	cout << "Generating from " << currentTerminals << " to " << nTerminals << "..." << endl;
	//	Be careful nTerminals may differ from the current amount of terminals since vessel-tip conexions are allowed in some cases.
	for (long long i = currentTerminals; i < nTerminals; i = tree->getNTerms()) {

		dataMonitor->update();

		if (i % saveInterval == 0 || i == currentTerminals) {
			saveStatus(i);
		}

		int invalidTerminal = true;
		int iTry = 0;
		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
		while (invalidTerminal) {

			do {
				xNew = domain->getRandomPoint();
			} while (!isValidSegment(xNew, ++iTry));

			int nNeighbors;
			vector<AbstractVascularElement *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

			double minCost = INFINITY;
			point minBif;
			AbstractVascularElement *minParent = NULL;
#pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
				point xBif;
				double cost;
				auto ogIt = ogVessels->find(static_cast<SingleVessel *>(neighborVessels[j])->coordToString());
				if (ogIt != ogVessels->end() && (*ogIt).second.second == false) {
					continue;
				} 
				tree->testVessel(xNew, neighborVessels[j], domain,
						neighborVessels, dLim, &xBif, &cost); //	Inf cost stands for invalid solution
#pragma omp critical
				{
					if (cost < minCost) {
						minCost = cost;
						minBif = xBif;
						minParent = neighborVessels[j];
					}
				}
			}
			//	end for trees

			if (minCost < INFINITY) {
				cout << "Added with a cost of " << minCost << " with a total cost of " << ((SingleVessel *) tree->getRoot())->treeVolume << endl;
				SingleVessel *minParentSV = (SingleVessel *) minParent;
				fwrite(&(minBif.p[0]), sizeof(double), 1, fp);
				fwrite(&(minBif.p[1]), sizeof(double), 1, fp);
				fwrite(&(minBif.p[2]), sizeof(double), 1, fp);
				fwrite(&(xNew.p[0]), sizeof(double), 1, fp);
				fwrite(&(xNew.p[1]), sizeof(double), 1, fp);
				fwrite(&(xNew.p[2]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xProx.p[0]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xProx.p[1]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xProx.p[2]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xDist.p[0]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xDist.p[1]), sizeof(double), 1, fp);
				fwrite(&(minParentSV->xDist.p[2]), sizeof(double), 1, fp);
				fwrite(&(instanceData->vesselFunction), sizeof(int), 1, fp);
				tree->addVessel(minBif, xNew, minParent, (AbstractVascularElement::VESSEL_FUNCTION) instanceData->vesselFunction);				
				invalidTerminal = false;
			}
		}
		dataMonitor->addDLimValue(dLim,i);
		domain->update();
		//	TODO Reset data monitor!!
		//	Terminal added.
	}
	tree->computePressure(tree->getRoot());
	tree->setPointCounter(domain->getPointCounter());

	this->endTime = time(nullptr);
	this->dLimLast = this->dLim;

	ogVessels->clear();
	delete ogVessels;

	saveStatus(nTerminals-1);
	markTimestampOnConfigurationFile("Final tree volume " + to_string(((SingleVessel *) tree->getRoot())->treeVolume));
	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	return tree;
}

void StagedFRROTreeGenerator::observableModified(IDomainObservable* observableInstance) {
	cout << "Changing instance parameters from " << endl << instanceData;
	instanceData = ((AbstractDomain *) observableInstance)->getInstanceData();
	cout << "To " << endl << instanceData << endl;
	tree->setInstanceData(instanceData);

	SingleVessel::bifurcationTests = instanceData->nBifurcationTest;
	this->stage = ((StagedDomain *) observableInstance)->getCurrentStage();
	cout << "Changing to stage " << stage << endl;
	this->tree->setCurrentStage(this->stage);
	this->tree->setGam(gams[stage]);
	this->tree->setEpsLim(epsLims[stage]);
	this->tree->setNu(nus[stage]);

	this->dataMonitor->reset();
}

AbstractObjectCCOTree*& StagedFRROTreeGenerator::getTree() {
	return tree;
}

void StagedFRROTreeGenerator::setSavingTasks(const vector<AbstractSavingTask*>& savingTasks){
	this->savingTasks = savingTasks;
}

void StagedFRROTreeGenerator::saveStatus(long long int terminals){
	tree->setPointCounter(domain->getPointCounter());
	for (std::vector<AbstractSavingTask *>::iterator it = savingTasks.begin(); it != savingTasks.end(); ++it) {
		(*it)->execute(terminals,tree);
	}
//			nodalWriter->write(tempDirectory+ "/step" + to_string(i) + "_view.vtp",tree);
	markTimestampOnConfigurationFile("Generating vessel #" + to_string(terminals));
	markTimestampOnConfigurationFile("Total RAM consumption: " + to_string(monitor->getProcessMemoryConsumption()) + " MB.");
}

vector<AbstractConstraintFunction<double, int> *>* StagedFRROTreeGenerator::getGams()
{
	return &(this->gams);
}

vector<AbstractConstraintFunction<double, int> *>* StagedFRROTreeGenerator::getEpsLims()
{
	return &(this->epsLims);
}
vector<AbstractConstraintFunction<double, int> *>* StagedFRROTreeGenerator::getNus()
{
	return &(this->nus);
}

double StagedFRROTreeGenerator::getDLim() {
	return this->dLim;
}

void StagedFRROTreeGenerator::setDLim(double newDLim) {
	this->dLim = newDLim;
}

double StagedFRROTreeGenerator::getDLimInitial() {
	return this->dLimInitial;
}

double StagedFRROTreeGenerator::getDLimLast() {
	return this->dLimLast;
}

time_t StagedFRROTreeGenerator::getBeginTime() {
	return this->beginTime;
}

time_t StagedFRROTreeGenerator::getEndTime() {
	return this->endTime;
}