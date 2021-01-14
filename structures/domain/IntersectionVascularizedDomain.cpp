/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * IntersectionVascularizedDomain.cpp
 *
 *  Created on: Mar 20, 2019
 *      Author: Gonzalo D. Maso Talou
 */

#include "IntersectionVascularizedDomain.h"

#include <chrono>
#include <omp.h>

//	Model
#include <vtkPolyDataReader.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPointData.h>
#include <vtkMassProperties.h>

IntersectionVascularizedDomain::IntersectionVascularizedDomain(vector<string> filenameVascularRegions, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filenameVR = filenameVascularRegions;
	boundingBox = new double[6];

	//	Read all the data from the vascularized files
	for (unsigned int i = 0; i < filenameVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameVascularRegions[i].c_str());
		reader2->Update();
		vtkVascularizedRegions.push_back(reader2->GetOutput());

		//	Create locator for non-vascularized regions.
		vtkSmartPointer<vtkOBBTree> locatorRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorRegion->SetDataSet(vtkVascularizedRegions[i]);
		locatorRegion->BuildLocator();
		sublocators.push_back(locatorRegion);

		//	Update boundary box
		if(i==0){
			double *refBB = vtkVascularizedRegions[0]->GetBounds();
			memcpy(boundingBox,refBB,sizeof(double)*6);
		}

		double *bbCurrent = vtkVascularizedRegions[i]->GetBounds();
		boundingBox[0] = max(boundingBox[0],bbCurrent[0]);
		boundingBox[2] = max(boundingBox[2],bbCurrent[2]);
		boundingBox[4] = max(boundingBox[4],bbCurrent[4]);
		boundingBox[1] = min(boundingBox[1],bbCurrent[1]);
		boundingBox[3] = min(boundingBox[3],bbCurrent[3]);
		boundingBox[5] = min(boundingBox[5],bbCurrent[5]);
	}

	nDraw = 10000;
	this->seed = chrono::system_clock::now().time_since_epoch().count();
	generator = mt19937(this->seed);

	characteristicLength = max(max((boundingBox[1] - boundingBox[0]) / 2,
			(boundingBox[3] - boundingBox[2]) / 2),
			(boundingBox[5] - boundingBox[4]) / 2);
}

IntersectionVascularizedDomain::IntersectionVascularizedDomain(vector<string> filenameVascularRegions,
		int N, GeneratorData *instanceData) : AbstractDomain(instanceData) {
	this->filenameVR = filenameVascularRegions;
	boundingBox = new double[6];

	//	Read all the data from the vascularized files
	for (unsigned int i = 0; i < filenameVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameVascularRegions[i].c_str());
		reader2->Update();
		vtkVascularizedRegions.push_back(reader2->GetOutput());

		//	Create locator for non-vascularized regions.
		vtkSmartPointer<vtkOBBTree> locatorRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorRegion->SetDataSet(vtkVascularizedRegions[i]);
		locatorRegion->BuildLocator();
		sublocators.push_back(locatorRegion);

		//	Update boundary box
		if(i==0){
			double *refBB = vtkVascularizedRegions[0]->GetBounds();
			memcpy(boundingBox,refBB,sizeof(double)*6);
		}

		double *bbCurrent = vtkVascularizedRegions[i]->GetBounds();
		boundingBox[0] = max(boundingBox[0],bbCurrent[0]);
		boundingBox[2] = max(boundingBox[2],bbCurrent[2]);
		boundingBox[4] = max(boundingBox[4],bbCurrent[4]);
		boundingBox[1] = min(boundingBox[1],bbCurrent[1]);
		boundingBox[3] = min(boundingBox[3],bbCurrent[3]);
		boundingBox[5] = min(boundingBox[5],bbCurrent[5]);
	}

	nDraw = N;
	this->seed = chrono::system_clock::now().time_since_epoch().count();
	generator = mt19937(this->seed);

	characteristicLength = max(max((boundingBox[1] - boundingBox[0]) / 2,
			(boundingBox[3] - boundingBox[2]) / 2),
			(boundingBox[5] - boundingBox[4]) / 2);
}

IntersectionVascularizedDomain::IntersectionVascularizedDomain(vector<string> filenameVascularRegions,
		int N, int seed, GeneratorData *instanceData) : AbstractDomain(instanceData) {
	this->filenameVR = filenameVascularRegions;
	boundingBox = new double[6];

	//	Read all the data from the vascularized files
	for (unsigned int i = 0; i < filenameVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameVascularRegions[i].c_str());
		reader2->Update();
		vtkVascularizedRegions.push_back(reader2->GetOutput());

		//	Create locator for non-vascularized regions.
		vtkSmartPointer<vtkOBBTree> locatorRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorRegion->SetDataSet(vtkVascularizedRegions[i]);
		locatorRegion->BuildLocator();
		sublocators.push_back(locatorRegion);

		//	Update boundary box
		if(i==0){
			double *refBB = vtkVascularizedRegions[0]->GetBounds();
			memcpy(boundingBox,refBB,sizeof(double)*6);
		}

		double *bbCurrent = vtkVascularizedRegions[i]->GetBounds();
		boundingBox[0] = max(boundingBox[0],bbCurrent[0]);
		boundingBox[2] = max(boundingBox[2],bbCurrent[2]);
		boundingBox[4] = max(boundingBox[4],bbCurrent[4]);
		boundingBox[1] = min(boundingBox[1],bbCurrent[1]);
		boundingBox[3] = min(boundingBox[3],bbCurrent[3]);
		boundingBox[5] = min(boundingBox[5],bbCurrent[5]);
	}

	nDraw = N;
	this->seed = seed;
	generator = mt19937(seed);

	characteristicLength = max(max((boundingBox[1] - boundingBox[0]) / 2,
			(boundingBox[3] - boundingBox[2]) / 2),
			(boundingBox[5] - boundingBox[4]) / 2);

}

IntersectionVascularizedDomain::~IntersectionVascularizedDomain() {
	this->randomInnerPoints.clear();
	this->filenameVR.clear();
	this->vtkVascularizedRegions.clear();
	this->sublocators.clear();
	delete[] this->boundingBox;
}

void IntersectionVascularizedDomain::generateRandomPoints() {
	uniform_real_distribution<double> distX(boundingBox[0], boundingBox[1]);
	uniform_real_distribution<double> distY(boundingBox[2], boundingBox[3]);
	uniform_real_distribution<double> distZ(boundingBox[4], boundingBox[5]);

	for (int i = 0; i < nDraw; ++i) {
		point p = { distX(generator), distY(generator), distZ(generator) };
		randomInnerPoints.push_back(p);
	}

	removeRandomOuterPoints();
}

void IntersectionVascularizedDomain::removeRandomOuterPoints() {

	vector<vtkSmartPointer<vtkPoints>> points;
//	cout << "Testing inside domain condition for " << randomInnerPoints.size() << " points" << endl;

	vtkSmartPointer<vtkPoints> currentPoints;
	for (unsigned i = 0; i < randomInnerPoints.size(); ++i) {
		if (i % 1000 == 0) {
			currentPoints = vtkSmartPointer<vtkPoints>::New();
			points.push_back(currentPoints);
		}
		currentPoints->InsertNextPoint(randomInnerPoints[i].p);
	}

	vector<vector<vtkSmartPointer<vtkSelectEnclosedPoints> > > allEnclosedPoints;
	vector<vector<vtkSmartPointer<vtkSelectEnclosedPoints> > > allEnclosedInNVR;
#pragma omp parallel for shared(allEnclosedPoints,allEnclosedInNVR), ordered, schedule(static,1), num_threads(omp_get_max_threads())
	for (unsigned j = 0; j < points.size(); ++j) {
		vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
		pointsPolydata->SetPoints(points[j]);

		vector<vtkSmartPointer<vtkSelectEnclosedPoints> > testVR;
		for (unsigned int i = 0; i < vtkVascularizedRegions.size(); ++i) {

			//	Points inside test
			vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
			selectEnclosedPoints->SetInputData(pointsPolydata);
			selectEnclosedPoints->SetSurfaceData(vtkVascularizedRegions[i]);
			selectEnclosedPoints->SetTolerance(1E-12);
			selectEnclosedPoints->Update();

			testVR.push_back(selectEnclosedPoints);
		}

#pragma omp ordered
#pragma omp critical
		{
			allEnclosedPoints.push_back(testVR);
		}
#pragma omp flush
	}

	randomInnerPoints.clear();

	for (unsigned j = 0; j < allEnclosedPoints.size(); ++j) {
		for (vtkIdType i = 0; i < points[j]->GetNumberOfPoints(); i++) {
			int contained = true;
			for (unsigned int l = 0; l < vtkVascularizedRegions.size() && contained; ++l) {
				if (!allEnclosedPoints[j][l]->IsInside(i)) {
					contained = false;
				}
			}
			if (contained) {
				double *pd = points[j]->GetPoint(i);
				point p = { pd[0], pd[1], pd[2] };
				randomInnerPoints.push_back(p);
			}
		}
	}
}

point IntersectionVascularizedDomain::getRandomPoint() {

	if (randomInnerPoints.empty())
		generateRandomPoints();
	point p = randomInnerPoints.front();
	randomInnerPoints.pop_front();

	++pointCounter;

	return p;
}

int IntersectionVascularizedDomain::isSegmentInside(point xs, point xf) {

	vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int hasIntersected;
	double tol = 1.e-8;

	for (unsigned int j = 0; j < sublocators.size(); ++j) {
		hasIntersected = abs(sublocators[j]->IntersectWithLine(xs.p, xf.p, points, ids));
		if (hasIntersected)
			return false;
	}

	return true;
}

double IntersectionVascularizedDomain::getSize() {

	if(volume == 0.0){
		vtkMassProperties *massProperty = vtkMassProperties::New();
		for (std::vector<vtkSmartPointer<vtkPolyData> >::iterator it = vtkVascularizedRegions.begin(); it != vtkVascularizedRegions.end(); ++it) {
			massProperty->SetInputData(*it);
			massProperty->Update();
			volume += massProperty->GetVolume();
		}
	}
	return volume;
}

vtkSmartPointer<vtkOBBTree>& IntersectionVascularizedDomain::getLocator()
{
	vtkSmartPointer<vtkOBBTree> nullInstance;
	return nullInstance;
}

int IntersectionVascularizedDomain::getDraw()
{
	return nDraw;
}

deque<point>& IntersectionVascularizedDomain::getRandomInnerPoints()
{
	return randomInnerPoints;
}

vtkSmartPointer<vtkPolyData>& IntersectionVascularizedDomain::getVtkGeometry()
{
	vtkSmartPointer<vtkPolyData> nullInstance;
	return nullInstance;
}

double IntersectionVascularizedDomain::getCharacteristicLength() {
	return characteristicLength;
}

double IntersectionVascularizedDomain::getDLim(long long int nVessels, double factor) {
	return characteristicLength * cbrt(factor / nVessels);
}

double* IntersectionVascularizedDomain::getLocalNeighborhood(point p, long long int nVessels) {
	double *localBox = new double[6];
	double size = instanceData->closeNeighborhoodFactor * getDLim(nVessels, instanceData->perfusionAreaFactor);

	localBox[0] = p.p[0] - size;
	localBox[1] = p.p[0] + size;
	localBox[2] = p.p[1] - size;
	localBox[3] = p.p[1] + size;
	localBox[4] = p.p[2] - size;
	localBox[5] = p.p[2] + size;

	return localBox;
}

int IntersectionVascularizedDomain::getSeed()
{
	return this->seed;
}

vector<string> IntersectionVascularizedDomain::getFilenameVR() {
	return this->filenameVR;
}

void IntersectionVascularizedDomain::logDomainFiles(FILE *fp) {
	fprintf(fp, "IntersectionVascularizedDomain\n");
    vector<string> filenameVR = this->getFilenameVR();
    int size = filenameVR.size();
    for (int i = 0; i < size; ++i) {
        fprintf(fp, "filenameVR[%d] = %s\n", i, filenameVR[i].c_str());
    }
}

vtkSmartPointer<vtkSelectEnclosedPoints> IntersectionVascularizedDomain::getEnclosedPoints() {
	vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	enclosedPoints->Initialize(this->vtkVascularizedRegions[0]);
	return enclosedPoints;
}