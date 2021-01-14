/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * DomainNVR.cpp
 *
 *  Created on: Dec 1, 2017
 *      Author: gonzalo
 */

#include "DomainNVR.h"

#include <chrono>
#include <random>
#include <omp.h>

//	Model
#include <vtkPolyDataReader.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPointData.h>
#include <vtkMassProperties.h>

DomainNVR::DomainNVR(string filename, vector<string> filenameNonVascularRegions, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filenameHull = filename;
	this->filenameNVR = filenameNonVascularRegions;
	//	Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->BuildLocator();

	//	Read all the data from the non-vascularized files
	for (unsigned int i = 0; i < filenameNonVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameNonVascularRegions[i].c_str());
		reader2->Update();
		vtkHollowRegions.push_back(reader2->GetOutput());

		//	Create the tree
		vtkSmartPointer<vtkOBBTree> locatorHollowRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorHollowRegion->SetDataSet(vtkHollowRegions[i]);
		locatorHollowRegion->BuildLocator();
		hollowLocators.push_back(locatorHollowRegion);
	}

	nDraw = 10000;

	this->seed = chrono::system_clock::now().time_since_epoch().count();
	generator = mt19937(this->seed);

	double *bb = vtkGeometry->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);
}

DomainNVR::DomainNVR(string filename, vector<string> filenameNonVascularRegions, int N, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filenameHull = filename;
	this->filenameNVR = filenameNonVascularRegions;
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->BuildLocator();

	//	Read all the data from the non-vascularized files
	for (unsigned int i = 0; i < filenameNonVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameNonVascularRegions[i].c_str());
		reader2->Update();
		vtkHollowRegions.push_back(reader2->GetOutput());

		//	Create the tree
		vtkSmartPointer<vtkOBBTree> locatorHollowRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorHollowRegion->SetDataSet(vtkHollowRegions[i]);
		locatorHollowRegion->BuildLocator();
		hollowLocators.push_back(locatorHollowRegion);
	}

	this->seed = chrono::system_clock::now().time_since_epoch().count();
	generator = mt19937(this->seed);

	nDraw = N;
	double *bb = vtkGeometry->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);

}

DomainNVR::DomainNVR(string filename, vector<string> filenameNonVascularRegions, int N, int seed, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filenameHull = filename;
	this->filenameNVR = filenameNonVascularRegions;
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->BuildLocator();

	//	Read all the data from the non-vascularized files
	for (unsigned int i = 0; i < filenameNonVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameNonVascularRegions[i].c_str());
		reader2->Update();
		vtkHollowRegions.push_back(reader2->GetOutput());

		//	Create the tree
		vtkSmartPointer<vtkOBBTree> locatorHollowRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorHollowRegion->SetDataSet(vtkHollowRegions[i]);
		locatorHollowRegion->BuildLocator();
		hollowLocators.push_back(locatorHollowRegion);
	}

	this->seed = seed;
	generator = mt19937(seed);

	nDraw = N;
	double *bb = vtkGeometry->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);

}

DomainNVR::~DomainNVR() {
	this->randomInnerPoints.clear();
	this->filenameNVR.clear();
	this->vtkHollowRegions.clear();
	this->hollowLocators.clear();
}

void DomainNVR::generateRandomPoints() {
	double *boundingBox = vtkGeometry->GetBounds();
	uniform_real_distribution<double> distX(boundingBox[0], boundingBox[1]);
	uniform_real_distribution<double> distY(boundingBox[2], boundingBox[3]);
	uniform_real_distribution<double> distZ(boundingBox[4], boundingBox[5]);

	for (int i = 0; i < nDraw; ++i) {
		point p = { distX(generator), distY(generator), distZ(generator) };
		randomInnerPoints.push_back(p);
	}

	removeRandomOuterPoints();
}

//void DomainNVR::removeRandomOuterPointsSerial() {
//
//	vtkSmartPointer<vtkPoints> points =
//			vtkSmartPointer<vtkPoints>::New();
//	cout << "Testing inside domain condition for " << randomInnerPoints.size() << " points" << endl;
//	for (unsigned i = 0; i < randomInnerPoints.size(); ++i) {
//		points->InsertNextPoint(randomInnerPoints[i].p);
//	}
//
//	vtkSmartPointer<vtkPolyData> pointsPolydata =
//			vtkSmartPointer<vtkPolyData>::New();
//	pointsPolydata->SetPoints(points);
//
//	//	Points inside hull test
//	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints =
//			vtkSmartPointer<vtkSelectEnclosedPoints>::New();
//	selectEnclosedPoints->SetInputData(pointsPolydata);
//	selectEnclosedPoints->SetSurfaceData(vtkGeometry);
//	selectEnclosedPoints->SetTolerance(1E-6);
//	selectEnclosedPoints->Update();
//
//	vector<vtkSmartPointer<vtkSelectEnclosedPoints> > testNVR;
//	for (unsigned int i = 0; i < vtkHollowRegions.size(); ++i) {
//		cout << "Testing inside NVR #" << i << " for " << randomInnerPoints.size() << " points" << endl;
//
//		//	Points inside non-vascularized test
//		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPointsByHoles =
//				vtkSmartPointer<vtkSelectEnclosedPoints>::New();
//		selectEnclosedPointsByHoles->SetInputData(pointsPolydata);
//		selectEnclosedPointsByHoles->SetSurfaceData(vtkHollowRegions[i]);
//		selectEnclosedPointsByHoles->SetTolerance(1E-6);
//		selectEnclosedPointsByHoles->Update();
//
//		testNVR.push_back(selectEnclosedPointsByHoles);
//	}
//
//	randomInnerPoints.clear();
//
//	cout << "Evaluating conditions for " << points->GetNumberOfPoints() << " points" << endl;
//	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
//		if (selectEnclosedPoints->IsInside(i)) {
//			int outsideHoles = true;
//			for (unsigned int j = 0; j < vtkHollowRegions.size() && outsideHoles; ++j) {
//				if (testNVR[j]->IsInside(i)) {
//					outsideHoles = false;
//				}
//			}
//			if (outsideHoles) {
//				double *pd = points->GetPoint(i);
//				point p = { pd[0], pd[1], pd[2] };
//				randomInnerPoints.push_back(p);
//			}
//		}
//	}
//
//}

void DomainNVR::removeRandomOuterPoints() {

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

	vector<vtkSmartPointer<vtkSelectEnclosedPoints>> allEnclosedPoints;
	vector<vector<vtkSmartPointer<vtkSelectEnclosedPoints> > > allEnclosedInNVR;
#pragma omp parallel for shared(allEnclosedPoints,allEnclosedInNVR), ordered, schedule(static,1), num_threads(omp_get_max_threads())
	for (unsigned j = 0; j < points.size(); ++j) {
		vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
		pointsPolydata->SetPoints(points[j]);

		//	Points inside test
		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->SetInputData(pointsPolydata);
		selectEnclosedPoints->SetSurfaceData(vtkGeometry);
		selectEnclosedPoints->SetTolerance(1E-12);
		selectEnclosedPoints->Update();

		vector<vtkSmartPointer<vtkSelectEnclosedPoints> > testNVR;
		for (unsigned int i = 0; i < vtkHollowRegions.size(); ++i) {
//			cout << "Testing inside NVR #" << i + 1 << " for point batch " << j << endl;

			//	Points inside non-vascularized test
			vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPointsByHoles =
					vtkSmartPointer<vtkSelectEnclosedPoints>::New();
			selectEnclosedPointsByHoles->SetInputData(pointsPolydata);
			selectEnclosedPointsByHoles->SetSurfaceData(vtkHollowRegions[i]);
			selectEnclosedPointsByHoles->SetTolerance(1E-12);
			selectEnclosedPointsByHoles->Update();

			testNVR.push_back(selectEnclosedPointsByHoles);
		}

#pragma omp ordered
#pragma omp critical
		allEnclosedPoints.push_back(selectEnclosedPoints);
		allEnclosedInNVR.push_back(testNVR);
#pragma omp flush
	}

	randomInnerPoints.clear();

	for (unsigned j = 0; j < allEnclosedPoints.size(); ++j) {
		for (vtkIdType i = 0; i < points[j]->GetNumberOfPoints(); i++) {
			if (allEnclosedPoints[j]->IsInside(i)) {
				int outsideHoles = true;
				for (unsigned int k = 0; k < vtkHollowRegions.size() && outsideHoles; ++k) {
					if (allEnclosedInNVR[j][k]->IsInside(i)) {
						outsideHoles = false;
					}
				}
				if (outsideHoles) {
					double *pd = points[j]->GetPoint(i);
					point p = { pd[0], pd[1], pd[2] };
					randomInnerPoints.push_back(p);
				}
			}
		}
	}

}

point DomainNVR::getRandomPoint() {

	if (randomInnerPoints.empty())
		generateRandomPoints();
	point p = randomInnerPoints.front();
	randomInnerPoints.pop_front();

	++pointCounter;

	return p;
}

int DomainNVR::isSegmentInside(point xs, point xf) {

	vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	bool hasIntersected;
	double tol = 1.e-8;

	if (!isConvexDomain) {
        locator->SetTolerance(tol);
		hasIntersected = abs(locator->IntersectWithLine(xs.p, xf.p, points, ids));
		if ( hasIntersected ){
//			cout << "Out of hull - (" << xs.p[0] << ", " << xs.p[1] << ", " << xs.p[2] << ") - (" << xf.p[0] << ", " << xf.p[1] << ", " << xf.p[2] << ")" << endl;
			return false;
		}
	}

	for (unsigned int j = 0; j < hollowLocators.size(); ++j) {
		hollowLocators[j]->SetTolerance(tol);
		hasIntersected = abs(hollowLocators[j]->IntersectWithLine(xs.p, xf.p, points, ids));
		if ( hasIntersected ){
//			cout << "Within not vascularized region - (" << xs.p[0] << ", " << xs.p[1] << ", " << xs.p[2] << ") - (" << xf.p[0] << ", " << xf.p[1] << ", " << xf.p[2] << ")" << endl;
			return false;
		}
	}

	return true;
}

double DomainNVR::getSize() {
	if (volume == 0.0) {
		vtkMassProperties *massProperty = vtkMassProperties::New();
		massProperty->SetInputData(vtkGeometry);
		massProperty->Update();
		volume = massProperty->GetVolume();
	}
	return volume;
}

vtkSmartPointer<vtkOBBTree>& DomainNVR::getLocator()
{
	return locator;
}

int DomainNVR::getDraw()
{
	return nDraw;
}

deque<point>& DomainNVR::getRandomInnerPoints()
{
	return randomInnerPoints;
}

vtkSmartPointer<vtkPolyData>& DomainNVR::getVtkGeometry()
{
	return vtkGeometry;
}

double DomainNVR::getCharacteristicLength() {
	return characteristicLength;
}

double DomainNVR::getDLim(long long int nVessels, double factor) {
	return characteristicLength * cbrt(factor / nVessels);
}

double* DomainNVR::getLocalNeighborhood(point p, long long int nVessels) {
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

int DomainNVR::getSeed()
{
	return this->seed;
}

string DomainNVR::getFilenameHull() {
	return this->filenameHull;
}

vector<string> DomainNVR::getFilenameNVR() {
	return this->filenameNVR;
}

void DomainNVR::logDomainFiles(FILE *fp) {
	fprintf(fp, "DomainNVR\n");
    fprintf(fp, "filenameHull = %s\n", this->getFilenameHull().c_str());
    vector<string> filenameNVR = this->getFilenameNVR();
    int size = filenameNVR.size();
    for (int i = 0; i < size; ++i) {
        fprintf(fp, "filenameNVR[%d] = %s\n", i, filenameNVR[i].c_str());
    }
}

vtkSmartPointer<vtkSelectEnclosedPoints> DomainNVR::getEnclosedPoints() {
	vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	enclosedPoints->Initialize(this->vtkGeometry);
	return enclosedPoints;
}