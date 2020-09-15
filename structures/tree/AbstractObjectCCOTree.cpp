/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * GenericVesselsTree.cpp
 *
 *  Created on: Mar 23, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "AbstractObjectCCOTree.h"

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkParametricSpline.h>
#include <vtkCardinalSpline.h>
#include <vtkSpline.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>

AbstractObjectCCOTree::AbstractObjectCCOTree(GeneratorData *instanceData) {

	this->instanceData = instanceData;

	this->pointCounter = 0;
	this->refPressure = 0.0;
	this->qReservedFactor = 0.0;

	this->xPerf = {0.0,0.0,0.0};
	this->qProx = 0;

	this->gam = nullptr;
	this->epsLim = nullptr;
	this->nu = nullptr;

	this->root = nullptr;

	this->psiFactor = 0.0;
	this->dp = 0.0;
	this->nTerms = 0;

	this->vtkTree = vtkSmartPointer<vtkPolyData>::New();
	this->vtkTreeLocator = vtkSmartPointer<vtkCellLocator>::New();

	this->currentStage = 0;
	this->isInCm = 0;

}

AbstractObjectCCOTree::AbstractObjectCCOTree(point xi, double qi, AbstractConstraintFunction<double, int> *gam,
		AbstractConstraintFunction<double, int> *epsLim,
		AbstractConstraintFunction<double, int> *nu, double refPressure, GeneratorData *instanceData) {

	this->instanceData = instanceData;

	this->pointCounter = 0;
	this->qReservedFactor = 0.0;

	this->xPerf = xi;
	this->qProx = qi;

	this->gam = gam;
	this->epsLim = epsLim;
	this->nu = nu;
	this->refPressure = refPressure;

	this->psiFactor = 0.0;
	this->dp = 0.0;
	this->nTerms = 0;

	this->root = NULL;

	this->vtkTree = vtkSmartPointer<vtkPolyData>::New();
	this->vtkTreeLocator = vtkSmartPointer<vtkCellLocator>::New();

	this->currentStage = 0;
	this->isInCm = 0;
}

AbstractObjectCCOTree::~AbstractObjectCCOTree() {
	for (auto it = elements.begin(); it != elements.end(); it++){
		delete it->second;
	}
	elements.clear();
}

AbstractConstraintFunction<double, int> *AbstractObjectCCOTree::getEpsLim() {
	return epsLim;
}

AbstractConstraintFunction<double, int> *AbstractObjectCCOTree::getGam() {
	return gam;
}

AbstractConstraintFunction<double, int> *AbstractObjectCCOTree::getNu() {
	return nu;
}

double AbstractObjectCCOTree::getQProx() {
	return qProx;
}

AbstractVascularElement *AbstractObjectCCOTree::getRoot() {
	return root;
}

unordered_map<long long, AbstractVascularElement*>& AbstractObjectCCOTree::getSegments() {
	return elements;
}

vtkSmartPointer<vtkPolyData> AbstractObjectCCOTree::getVtkTree() {
	return vtkTree;
}

point AbstractObjectCCOTree::getXProx() {
	return xPerf;
}

void AbstractObjectCCOTree::save(string filename) {
	ofstream treeFile;

	treeFile.open(filename.c_str(), ios::out);
	treeFile.setf(ios::scientific, ios::floatfield);
	treeFile.precision(16);

	treeFile << "*Tree" << endl;
	saveTree(&treeFile);
	treeFile << endl << endl;

	treeFile << "*Vessels" << endl << elements.size() << endl;

	saveVessels(this->root, &treeFile);
	treeFile << endl;

	treeFile << "*Connectivity" << endl;
	saveConnectivity(this->root, &treeFile);

	treeFile.flush();
	treeFile.close();
}

void AbstractObjectCCOTree::saveVessels(AbstractVascularElement * root, ofstream *treeFile){
	if(!root){
		return;
	}
	root->saveVesselData(treeFile);
	*treeFile << endl;
	for (std::vector<AbstractVascularElement *>::iterator it = root->children.begin(); it != root->children.end(); ++it) {
		saveVessels(*it,treeFile);
	}
}

void AbstractObjectCCOTree::saveConnectivity(AbstractVascularElement * root, ofstream *treeFile){
	if(!root){
		return;
	}
	root->saveVesselConnectivity(treeFile);
	*treeFile << endl;
	for (std::vector<AbstractVascularElement *>::iterator it = root->children.begin(); it != root->children.end(); ++it) {
		saveConnectivity(*it,treeFile);
	}
}

void AbstractObjectCCOTree::saveTree(ofstream* outFile) {

	*outFile << xPerf.p[0] << " " << xPerf.p[1] << " " << xPerf.p[2] << " " << qProx << " " << psiFactor << " " << dp
			<< " " << nTerms << " " << refPressure << " "
			<< (long long) pointCounter << " ";

}

double AbstractObjectCCOTree::getDp() const {
	return dp;
}

void AbstractObjectCCOTree::setEpsLim(AbstractConstraintFunction<double, int> *epsLim) {
	this->epsLim = epsLim;
}

void AbstractObjectCCOTree::setGam(AbstractConstraintFunction<double, int> *gam) {
	this->gam = gam;
}

void AbstractObjectCCOTree::setNu(AbstractConstraintFunction<double, int> *nu) {
	this->nu = nu;
}

long long int AbstractObjectCCOTree::getPointCounter() const {
	return pointCounter;
}

void AbstractObjectCCOTree::setPointCounter(long long int pointCounter) {
	this->pointCounter = pointCounter;
}

long long int AbstractObjectCCOTree::getNTerms() {
	return this->nTerms;
}

long long int AbstractObjectCCOTree::getNTerminals() {
	return countTerminals(this->root);
}

long long int AbstractObjectCCOTree::countTerminals(AbstractVascularElement* root) {
	long long int childrenTerminals = 0l;
	for (std::vector<AbstractVascularElement *>::iterator it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
		childrenTerminals += countTerminals(*it);
	}
	return root->getTerminals() + childrenTerminals;
}

long long int AbstractObjectCCOTree::countTerminals(AbstractVascularElement* root, AbstractVascularElement::TERMINAL_TYPE type) {
	long long int childrenTerminals = 0l;
	for (std::vector<AbstractVascularElement *>::iterator it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
		childrenTerminals += countTerminals(*it, type);
	}
	return root->getTerminals(type) + childrenTerminals;
}

double AbstractObjectCCOTree::computeTreeCost(AbstractVascularElement* root) {

	double currentCost = ((SingleVessel *) root)->getVolume();
	vector<AbstractVascularElement *> children = root->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		currentCost += computeTreeCost(*it);
	}
	return currentCost;
}

void AbstractObjectCCOTree::computePressure(AbstractVascularElement* root) {

	vector<AbstractVascularElement *> children = root->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = root->getChildren().begin(); it != root->getChildren().end(); ++it) {
		computePressure(*it);
	}
	root->updatePressure();
//
//	if (root->anastomose.size() > 1) {
//		computePressure(root->anastomose[1]);
//		computePressure(root->anastomose[2]);
//		root->localResistance = 8 * root->viscosity / M_PI * root->length;
//		root->pressure = root->flux * root->localResistance + root->anastomose[1]->pressure + root->anastomose[2]->pressure;
//	} else {
//		root->localResistance = 8 * root->viscosity / M_PI * root->length;
//		root->pressure = root->flux * root->localResistance;
//	}
}

void AbstractObjectCCOTree::updateSegmentVtkLines() {
	vtkTree->GetLines()->InitTraversal();
	for (auto it = elements.begin(); it != elements.end(); ++it) {
		for (auto it2 = (it->second)->getVessels().begin(); it2 != (it->second)->getVessels().end(); ++it2) {
			(*it2)->vtkSegment = vtkSmartPointer<vtkLine>::New();
			vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
			vtkTree->GetLines()->GetNextCell(pts);
			(*it2)->vtkSegment->GetPointIds()->SetId(0, pts->GetId(0));
			(*it2)->vtkSegment->GetPointIds()->SetId(1, pts->GetId(1));
		}
	}
}

void AbstractObjectCCOTree::printVtkTree() {

	vtkTree->BuildCells();
	int numLines = vtkTree->GetNumberOfLines();
	int numPoints = vtkTree->GetNumberOfPoints();

	cout << "Tree with " << numLines << " segments and " << numPoints << " points." << endl;

	vtkSmartPointer<vtkCellArray> lines = vtkTree->GetLines();
	lines->InitTraversal();
	vtkSmartPointer<vtkIdList> pts = vtkSmartPointer<vtkIdList>::New();
	int i = 0;
	while (lines->GetNextCell(pts)) {
		cout << i++ << ": ";
		for (int j = 0; j < pts->GetNumberOfIds(); ++j) {
			cout << pts->GetId(j) << " ";
		}
		cout << endl;
	}

	vtkSmartPointer<vtkPoints> points = vtkTree->GetPoints();
	for (int i = 0; i < numPoints; ++i) {
		cout << i << ": (" << points->GetPoint(i)[0] << ", " << points->GetPoint(i)[1] << ", " << points->GetPoint(i)[2] << ")" << endl;
	}

}

GeneratorData* AbstractObjectCCOTree::getInstanceData() {
	return instanceData;
}

void AbstractObjectCCOTree::setInstanceData(GeneratorData* instanceData) {
	this->instanceData = instanceData;
}

int AbstractObjectCCOTree::getCurrentStage() const
{
	return currentStage;
}

void AbstractObjectCCOTree::setCurrentStage(int currentStage) {
	this->currentStage = currentStage;
}

double AbstractObjectCCOTree::getReservedFactor() const
{
	return qReservedFactor;
}

void AbstractObjectCCOTree::setReservedFactor(double reservedFactor)
		{
	qReservedFactor = reservedFactor;
}

long long int AbstractObjectCCOTree::getNTerminals(AbstractVascularElement::TERMINAL_TYPE type) {
	return countTerminals(this->root, type);
}

vector<vector<double> > AbstractObjectCCOTree::getVertices() {
	vector<vector<double> > points;
	vtkSmartPointer<vtkPoints> vtkPointsData = vtkTree->GetPoints(); // Points
	for (unsigned int i = 0; i < vtkPointsData->GetNumberOfPoints(); ++i) {
		vector<double> point;
		double *coordinates = vtkPointsData->GetPoint(i);
		point.push_back(coordinates[0]);
		point.push_back(coordinates[1]);
		point.push_back(coordinates[2]);
		points.push_back(point);
	}

	return points;

}

vector<vector<int> > AbstractObjectCCOTree::getConnectivity() {
	vector<vector<int> > lines;
	for (auto it = elements.begin(); it != elements.end(); ++it) {
		for (vector<SingleVessel *>::iterator it2 = (it->second)->getVessels().begin(); it2 != (it->second)->getVessels().end(); ++it2) {
			SingleVessel *currentSegment = *it2;
			vector<int> line;
			vector<double> radius;

			vtkSmartPointer<vtkIdList> vtkLineData = currentSegment->vtkSegment->GetPointIds(); //	Lines
			line.push_back(vtkLineData->GetId(0) + 1);
			line.push_back(vtkLineData->GetId(1) + 1);
			lines.push_back(line);
		}
	}
	return lines;
}

vector<SingleVessel*> AbstractObjectCCOTree::getVessels(){
	vector<SingleVessel *> allVessels;
	for (auto it = elements.begin(); it != elements.end(); ++it) {
		vector<SingleVessel *> currentVessels = (it->second)->getVessels();
		allVessels.insert(allVessels.end(),currentVessels.begin(),currentVessels.end());
	}

	return allVessels;
}

int AbstractObjectCCOTree::getIsInCm() const
{
	return isInCm;
}

void AbstractObjectCCOTree::setIsInCm(int isInCm)
		{
	this->isInCm = isInCm;
}
