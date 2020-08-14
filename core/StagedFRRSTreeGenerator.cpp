/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * StagedFixedPerfusionRadiusTreeGenerator.cpp
 *
 *  Created on: Jun 2, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#include "StagedFRRSTreeGenerator.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "../constrains/AbstractConstraintFunction.h"
#include "../structures/tree/FRRCCOSTree.h"
#include "../structures/tree/FRRVariableViscosityCCOSTree.h"
#include "../structures/tree/FRRVaViOptCCOSTree.h"
#include "../utils/MemoryMonitor.h"

#include <omp.h>

StagedFRRSTreeGenerator::StagedFRRSTreeGenerator(
		StagedDomain* domain, point xi, double rootRadii, double qi,
		long long nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu,
		double minAngle, double refPressure) {

	domain->registerObserver(this);
	this->domain = domain;
	this->instanceData = domain->getInstanceData();
	this->nTerminals = nTerm;
	this->stage = 0;

	this->gams = gam;
	this->epsLims = epsLim;
	this->nus = nu;

	this->tree = new FRRCCOSTree(xi, rootRadii, qi, gam[stage],
			epsLim[stage], nu[stage], minAngle, refPressure, instanceData);
	this->dLim = domain->getDLim(1, instanceData->perfusionAreaFactor);

	this->isGeneratingConfFile = 0;
	this->confFilename = "";

	this->dataMonitor = new GeneratorDataMonitor(domain);

}

StagedFRRSTreeGenerator::StagedFRRSTreeGenerator(
		StagedDomain* domain, point xi, double rootRadii, double qi,
		long long nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu,
		double minAngle, double refPressure, double viscosityTolerance, int optimized) {

	domain->registerObserver(this);
	this->domain = domain;
	this->instanceData = domain->getInstanceData();
	this->nTerminals = nTerm;
	this->stage = 0;

	this->gams = gam;
	this->epsLims = epsLim;
	this->nus = nu;

	if (optimized) {
		this->tree = new FRRVaViOptCCOSTree(xi, rootRadii, qi, gam[stage], epsLim[stage],
				nu[stage], minAngle, refPressure, viscosityTolerance, instanceData);
	} else {
		this->tree = new FRRVariableViscosityCCOSTree(xi, rootRadii,
				qi, gam[stage], epsLim[stage], nu[stage], minAngle, refPressure, viscosityTolerance, instanceData);
	}
	this->dLim = domain->getDLim(1, instanceData->perfusionAreaFactor);

	this->isGeneratingConfFile = 0;
	this->confFilename = "";

	this->dataMonitor = new GeneratorDataMonitor(domain);
}

StagedFRRSTreeGenerator::StagedFRRSTreeGenerator(
		StagedDomain* domain, AbstractStructuredCCOTree* tree, long long nTerm,
		vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu) {

	domain->registerObserver(this);
	this->domain = domain;
	this->instanceData = domain->getInstanceData();
	this->nTerminals = nTerm;
	this->tree = tree;
	//	Stage can be loaded from file
	this->stage = 0;

	this->gams = gam;
	this->epsLims = epsLim;
	this->nus = nu;

	this->tree->setGam(gam[stage]);
	this->tree->setEpsLim(epsLim[stage]);
	this->tree->setNu(nu[stage]);

	this->dataMonitor = new GeneratorDataMonitor(domain);
}

AbstractStructuredCCOTree *StagedFRRSTreeGenerator::generate() {

	MemoryMonitor *monitor = new MemoryMonitor(MemoryMonitor::MEGABYTE);
	generatesConfigurationFile(ios::out);

	point xNew;
	int i = 0;

	do {
		xNew = domain->getRandomPoint();
	} while (!isValidRootSegment(xNew, ++i));

	tree->addVessel(xNew, xNew, NULL);

	for (long long i = 1; i < nTerminals; ++i) {

		dataMonitor->update();

		if (i % 1000 == 0) {
			tree->storeVTK("./outputs/stepStruct" + to_string(i) + ".vtp");
			markTimestampOnConfigurationFile("Generating vessel #" + to_string(i));
			markTimestampOnConfigurationFile("Total RAM consumption: " + to_string(monitor->getProcessMemoryConsumption()) + " MB.");
		}

		int invalidTerminal = true;
		int iTry = 0;
		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
		while (invalidTerminal) {

			do {
				xNew = domain->getRandomPoint();
			} while (!isValidSegment(xNew, ++iTry));

			int nNeighbors;
			vector<vessel *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

			double minCost = INFINITY;
			point minBif;
			vessel *minParent = NULL;
//#pragma omp parallel for shared(minCost, minBif, minParent), schedule(dynamic,1), num_threads(omp_get_max_threads())
			for (unsigned j = 0; j < neighborVessels.size(); ++j) {
				point xBif;
				double cost;
				tree->testVessel(xNew, neighborVessels[j], domain,
						neighborVessels, dLim, &xBif, &cost); //	Inf cost stands for invalid solution
//#pragma omp critical
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
				tree->addVessel(minBif, xNew, minParent);
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

	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	delete monitor;

	return tree;
}

int StagedFRRSTreeGenerator::isValidRootSegment(point xNew,
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

int StagedFRRSTreeGenerator::isValidSegment(point xNew, int iTry) {

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

StagedDomain * StagedFRRSTreeGenerator::getDomain() {
	return domain;
}

vector<AbstractStructuredCCOTree*> StagedFRRSTreeGenerator::getTrees() {
	vector<AbstractStructuredCCOTree*> trees;
	trees.push_back(tree);
	return trees;
}

void StagedFRRSTreeGenerator::enableConfigurationFile(
		string filename) {
	this->isGeneratingConfFile = 1;
	this->confFilename = filename;
}

void StagedFRRSTreeGenerator::generatesConfigurationFile(ios::openmode mode) {

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
		confFile << "MINIMUM_ANGLE " << this->tree->getMinAngle() << endl;
		confFile << endl;

		confFile << "*Timestamps" << endl;
		beginningTime = chrono::steady_clock::now();
	}
}

void StagedFRRSTreeGenerator::markTimestampOnConfigurationFile(
		string label) {
	if (isGeneratingConfFile) {
		confFile << (chrono::duration_cast<chrono::microseconds>(
				chrono::steady_clock::now() - beginningTime).count())
				/ 1000000.0 << ": " << label << endl;
	}
}

void StagedFRRSTreeGenerator::closeConfigurationFile() {

	confFile << endl << "DOMAIN_POINTS_GENERATED " << this->domain->getPointCounter() << endl;
	if (isGeneratingConfFile) {
		confFile.flush();
		confFile.close();
	}
}

AbstractStructuredCCOTree* StagedFRRSTreeGenerator::resume() {

	MemoryMonitor *monitor = new MemoryMonitor(MemoryMonitor::MEGABYTE);
	generatesConfigurationFile(ios::out);

	for (long long int j = 0; j < tree->getPointCounter(); ++j) {
		domain->getRandomPoint();
	}

	//	Compute current nTerm
	long long currentTerminals = tree->getNTerminals();
	point xNew;

	cout << "Generating from " << currentTerminals << " to " << nTerminals << "..." << endl;
	for (long long i = currentTerminals; i < nTerminals; ++i) {

		dataMonitor->update();

		if (i % 250 == 0) {
			tree->storeVTK("/home/gonzalo/workspace/exCCOExecuteTreeGenerator/outputs/step" + to_string(i) + ".vtp");
			markTimestampOnConfigurationFile("Generating vessel #" + to_string(i));
			markTimestampOnConfigurationFile("Total RAM consumption: " + to_string(monitor->getProcessMemoryConsumption()) + " MB.");
		}

		int invalidTerminal = true;
		int iTry = 0;
		dLim = instanceData->dLimCorrectionFactor * domain->getDLim(i, instanceData->perfusionAreaFactor);
		while (invalidTerminal) {

			do {
				xNew = domain->getRandomPoint();
			} while (!isValidSegment(xNew, ++iTry));

			int nNeighbors;
			vector<vessel *> neighborVessels = tree->getCloseSegments(xNew, domain, &nNeighbors);
			cout << "Trying segment #" << i << " at terminal point " << xNew << " with the " << nNeighbors << " closest neighbors (dLim = " << dLim << ")." << endl;

			double minCost = INFINITY;
			point minBif;
			vessel *minParent = NULL;
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
				tree->addVessel(minBif, xNew, minParent);
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

	markTimestampOnConfigurationFile("Tree successfully generated.");
	closeConfigurationFile();

	delete monitor;

	return tree;

}

void StagedFRRSTreeGenerator::observableModified(IDomainObservable* observableInstance) {
	cout << "Changing instance parameters from " << endl << instanceData;
	instanceData = ((AbstractDomain *) observableInstance)->getInstanceData();
	cout << "To " << endl << instanceData << endl;
	tree->setInstanceData(instanceData);

	this->stage = ((StagedDomain *) observableInstance)->getCurrentStage();
	this->tree->setGam(gams[stage]);
	this->tree->setEpsLim(epsLims[stage]);
	this->tree->setNu(nus[stage]);

}
