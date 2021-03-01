/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * PartiallyVascularizedDomain.cpp
 *
 *  Created on: Mar 21, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "PartiallyVascularizedDomain.h"

#include <chrono>
#include <omp.h>

//	Model
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPointData.h>
#include <vtkMassProperties.h>

PartiallyVascularizedDomain::PartiallyVascularizedDomain(string filename, vector<string> filenameVascularRegions,
		vector<string> filenameNonVascularRegions, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filenameHull = filename;
	this->filenameVR = filenameVascularRegions;
	this->filenameNVR = filenameNonVascularRegions;
	//	Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkTransportRegion = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkTransportRegion);
	locator->BuildLocator();

	//	Read all the data from the vascularized files
	for (unsigned int i = 0; i < filenameVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameVascularRegions[i].c_str());
		reader2->Update();
		vtkVascularizedRegions.push_back(reader2->GetOutput());
	}

	//	Read all the data from the non-vascularized files
	for (unsigned int i = 0; i < filenameNonVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameNonVascularRegions[i].c_str());
		reader2->Update();
		vtkHollowRegions.push_back(reader2->GetOutput());

		//	Create locator for non-vascularized regions.
		vtkSmartPointer<vtkOBBTree> locatorHollowRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorHollowRegion->SetDataSet(vtkHollowRegions[i]);
		locatorHollowRegion->BuildLocator();
		hollowLocators.push_back(locatorHollowRegion);
	}

	nDraw = 10000;
	this->seed = chrono::system_clock::now().time_since_epoch().count();
	generator = mt19937(this->seed);

	double *bb = vtkTransportRegion->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);
}

PartiallyVascularizedDomain::PartiallyVascularizedDomain(string filename, vector<string> filenameVascularRegions,
		vector<string> filenameNonVascularRegions, int N, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filenameHull = filename;
	this->filenameVR = filenameVascularRegions;
	this->filenameNVR = filenameNonVascularRegions;
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkTransportRegion = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkTransportRegion);
	locator->BuildLocator();

	//	Read all the data from the vascularized files
	for (unsigned int i = 0; i < filenameVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameVascularRegions[i].c_str());
		reader2->Update();
		vtkVascularizedRegions.push_back(reader2->GetOutput());
	}

	//	Read all the data from the non-vascularized files
	for (unsigned int i = 0; i < filenameNonVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameNonVascularRegions[i].c_str());
		reader2->Update();
		vtkHollowRegions.push_back(reader2->GetOutput());

		//	Create locator for non-vascularized regions.
		vtkSmartPointer<vtkOBBTree> locatorHollowRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorHollowRegion->SetDataSet(vtkHollowRegions[i]);
		locatorHollowRegion->BuildLocator();
		hollowLocators.push_back(locatorHollowRegion);
	}

	nDraw = N;
	this->seed = chrono::system_clock::now().time_since_epoch().count();
	generator = mt19937(this->seed);

	double *bb = vtkTransportRegion->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);

}

PartiallyVascularizedDomain::PartiallyVascularizedDomain(string filename, vector<string> filenameVascularRegions,
		vector<string> filenameNonVascularRegions, int N, int seed, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filenameHull = filename;
	this->filenameVR = filenameVascularRegions;
	this->filenameNVR = filenameNonVascularRegions;
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkTransportRegion = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkTransportRegion);
	locator->BuildLocator();

	//	Read all the data from the vascularized files
	for (unsigned int i = 0; i < filenameVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameVascularRegions[i].c_str());
		reader2->Update();
		vtkVascularizedRegions.push_back(reader2->GetOutput());
	}

	//	Read all the data from the non-vascularized files
	for (unsigned int i = 0; i < filenameNonVascularRegions.size(); ++i) {
		vtkSmartPointer<vtkPolyDataReader> reader2 = vtkSmartPointer<vtkPolyDataReader>::New();
		reader2->SetFileName(filenameNonVascularRegions[i].c_str());
		reader2->Update();
		vtkHollowRegions.push_back(reader2->GetOutput());

		//	Create locator for non-vascularized regions.
		vtkSmartPointer<vtkOBBTree> locatorHollowRegion = vtkSmartPointer<vtkOBBTree>::New();
		locatorHollowRegion->SetDataSet(vtkHollowRegions[i]);
		locatorHollowRegion->BuildLocator();
		hollowLocators.push_back(locatorHollowRegion);
	}

	nDraw = N;
	this->seed = seed;
	generator = mt19937(seed);

	double *bb = vtkTransportRegion->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);

}

PartiallyVascularizedDomain::~PartiallyVascularizedDomain() {
	this->randomInnerPoints.clear();
	this->filenameVR.clear();
	this->filenameNVR.clear();
	this->vtkVascularizedRegions.clear();
	this->vtkHollowRegions.clear();
	this->hollowLocators.clear();
}

void PartiallyVascularizedDomain::generateRandomPoints() {
	double *boundingBox = vtkTransportRegion->GetBounds();
	uniform_real_distribution<double> distX(boundingBox[0], boundingBox[1]);
	uniform_real_distribution<double> distY(boundingBox[2], boundingBox[3]);
	uniform_real_distribution<double> distZ(boundingBox[4], boundingBox[5]);

	for (int i = 0; i < nDraw; ++i) {
		point p = { distX(generator), distY(generator), distZ(generator) };
		randomInnerPoints.push_back(p);
	}

	removeRandomOuterPoints();
}

void PartiallyVascularizedDomain::removeRandomOuterPoints() {

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
		{
			allEnclosedPoints.push_back(testVR);
			allEnclosedInNVR.push_back(testNVR);
		}
#pragma omp flush
	}

	randomInnerPoints.clear();

	for (unsigned j = 0; j < allEnclosedPoints.size(); ++j) {
		for (vtkIdType i = 0; i < points[j]->GetNumberOfPoints(); i++) {
			int contained = false;
			for (unsigned int l = 0; l < vtkVascularizedRegions.size() && !contained; ++l) {
				if (allEnclosedPoints[j][l]->IsInside(i)) {
					contained = true;
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
}

point PartiallyVascularizedDomain::getRandomPoint() {

	if (randomInnerPoints.empty())
		generateRandomPoints();
	point p = randomInnerPoints.front();
	randomInnerPoints.pop_front();

	++pointCounter;

	return p;
}

int PartiallyVascularizedDomain::isSegmentInside(point xs, point xf) {

	vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int hasIntersected;
	double tol = 1.e-8;

	if (!isConvexDomain) {
		hasIntersected = abs(locator->IntersectWithLine(xs.p, xf.p, points, ids));
		if (hasIntersected)
			return false;
	}

	for (unsigned int j = 0; j < hollowLocators.size(); ++j) {
		hasIntersected = abs(hollowLocators[j]->IntersectWithLine(xs.p, xf.p, points, ids));
		if (hasIntersected)
			return false;
	}

	return true;
}

double PartiallyVascularizedDomain::getSize() {

	if (volume == 0.0) {
		vtkMassProperties *massProperty = vtkMassProperties::New();
		for (std::vector<vtkSmartPointer<vtkPolyData> >::iterator it = vtkVascularizedRegions.begin(); it != vtkVascularizedRegions.end();
				++it) {
			massProperty->SetInputData(*it);
			massProperty->Update();
			volume += massProperty->GetVolume();
		}
	}
	return volume;
}

vtkSmartPointer<vtkOBBTree>& PartiallyVascularizedDomain::getLocator()
{
	return locator;
}

int PartiallyVascularizedDomain::getDraw()
{
	return nDraw;
}

deque<point>& PartiallyVascularizedDomain::getRandomInnerPoints()
{
	return randomInnerPoints;
}

vtkSmartPointer<vtkPolyData>& PartiallyVascularizedDomain::getVtkGeometry()
{
	return vtkTransportRegion;
}

double PartiallyVascularizedDomain::getCharacteristicLength() {
	return characteristicLength;
}

double PartiallyVascularizedDomain::getDLim(long long int nVessels, double factor) {
	return characteristicLength * cbrt(factor / nVessels);
}

double* PartiallyVascularizedDomain::getLocalNeighborhood(point p, long long int nVessels) {
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

void PartiallyVascularizedDomain::savePoints(string filename) {
	// Create 10 points.
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (std::deque<point>::iterator it = randomInnerPoints.begin(); it != randomInnerPoints.end(); ++it) {
		points->InsertNextPoint(it->p[0], it->p[1], it->p[2]);
	}

	// Create a polydata object and add the points to it.
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(polydata);

	writer->SetDataModeToBinary();
	writer->Write();
}

int PartiallyVascularizedDomain::getSeed()
{
	return this->seed;
}

string PartiallyVascularizedDomain::getFilenameHull() {
	return this->filenameHull;
}

vector<string> PartiallyVascularizedDomain::getFilenameVR() {
	return this->filenameVR;
}

vector<string> PartiallyVascularizedDomain::getFilenameNVR() {
	return this->filenameNVR;
}

void PartiallyVascularizedDomain::logDomainFiles(FILE *fp) {
	string filenameHull = this->getFilenameHull();
    vector<string> filenameVR = this->getFilenameVR();
    int size_vr = filenameVR.size();
    vector<string> filenameNVR = this->getFilenameNVR();
    int size_nvr = filenameNVR.size();
    fprintf(fp, "PartiallyVascularizedDomain\n");
    fprintf(fp, "filenameHull = %s\n", filenameHull.c_str());
    for (int i = 0; i < size_vr; ++i) {
        fprintf(fp, "filenameVR[%d] = %s\n", i, filenameVR[i].c_str());
    }
    for (int i = 0; i < size_nvr; ++i) {
        fprintf(fp, "filenameNVR[%d] = %s\n", i, filenameNVR[i].c_str());
    }
}

vtkSmartPointer<vtkSelectEnclosedPoints> PartiallyVascularizedDomain::getEnclosedPoints() {
	vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	enclosedPoints->Initialize(this->vtkTransportRegion);
	return enclosedPoints;
}