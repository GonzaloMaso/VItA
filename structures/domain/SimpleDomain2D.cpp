/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * Mesh.cpp
 *
 *  Created on: May 24, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#include "SimpleDomain2D.h"

#include <chrono>
#include <random>
#include <omp.h>

//	Model
#include <vtkPolyDataReader.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPointData.h>
#include <vtkMassProperties.h>

SimpleDomain2D::SimpleDomain2D(string filename, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filename = filename;
	//	Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<
			vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->BuildLocator();

	nDraw = 10000;
	this->seed = chrono::system_clock::now().time_since_epoch().count();
	generator = mt19937(this->seed);

	double *bb = vtkGeometry->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);

}

SimpleDomain2D::SimpleDomain2D(string filename, int N, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filename = filename;
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<
			vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->BuildLocator();

	nDraw = N;
	this->seed = chrono::system_clock::now().time_since_epoch().count();
	generator = mt19937(this->seed);

	double *bb = vtkGeometry->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);

}

SimpleDomain2D::SimpleDomain2D(string filename, int N, int seed, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filename = filename;			
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<
			vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->BuildLocator();

	nDraw = N;

	this->seed = seed;
	generator = mt19937(seed);

	double *bb = vtkGeometry->GetBounds();
	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);

}

SimpleDomain2D::~SimpleDomain2D() {
	this->randomInnerPoints.clear();
}

void SimpleDomain2D::generateRandomPoints() {
	if (seed == -1)
		seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 generator(seed);
	double *boundingBox = vtkGeometry->GetBounds();
	uniform_real_distribution<double> distX(boundingBox[0], boundingBox[1]);
	uniform_real_distribution<double> distY(boundingBox[2], boundingBox[3]);

	for (int i = 0; i < nDraw; ++i) {
		point p = { distX(generator), distY(generator), 0.0 };
		randomInnerPoints.push_back(p);
	}

	removeRandomOuterPoints();
}

void SimpleDomain2D::removeRandomOuterPointsSerial() {

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	cout << "Testing inside domain condition for " << randomInnerPoints.size()
			<< " points" << endl;

	for (unsigned i = 0; i < randomInnerPoints.size(); ++i) {
		points->InsertNextPoint(randomInnerPoints[i].p);
	}

	vtkSmartPointer<vtkPolyData> pointsPolydata =
			vtkSmartPointer<vtkPolyData>::New();
	pointsPolydata->SetPoints(points);

	//	Points inside test
	vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints =
			vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	selectEnclosedPoints->SetInputData(pointsPolydata);
	selectEnclosedPoints->SetSurfaceData(vtkGeometry);
	selectEnclosedPoints->SetTolerance(1E-6);
	selectEnclosedPoints->Update();

	randomInnerPoints.clear();

	for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++) {
		if (selectEnclosedPoints->IsInside(i)) {
			double *pd = points->GetPoint(i);
			point p = { pd[0], pd[1], pd[2] };
			randomInnerPoints.push_back(p);
		}
	}

}

void SimpleDomain2D::removeRandomOuterPoints() {

	vector<vtkSmartPointer<vtkPoints>> points;
	cout << "Testing inside domain condition for " << randomInnerPoints.size()
			<< " points" << endl;

	vtkSmartPointer<vtkPoints> currentPoints;
	for (unsigned i = 0; i < randomInnerPoints.size(); ++i) {
		if (i % 1000 == 0) {
			currentPoints = vtkSmartPointer<vtkPoints>::New();
			points.push_back(currentPoints);
		}
		currentPoints->InsertNextPoint(randomInnerPoints[i].p);
	}

	vector<vtkSmartPointer<vtkSelectEnclosedPoints>> allEnclosedPoints;
#pragma omp parallel for shared(allEnclosedPoints), ordered, schedule(static,1), num_threads(omp_get_max_threads())
	for (unsigned j = 0; j < points.size(); ++j) {
		vtkSmartPointer<vtkPolyData> pointsPolydata =
				vtkSmartPointer<vtkPolyData>::New();
		pointsPolydata->SetPoints(points[j]);

		//	Points inside test
		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints =
				vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->SetInputData(pointsPolydata);
		selectEnclosedPoints->SetSurfaceData(vtkGeometry);
		selectEnclosedPoints->SetTolerance(1E-6);
		selectEnclosedPoints->Update();

#pragma omp ordered
#pragma omp critical
		allEnclosedPoints.push_back(selectEnclosedPoints);
#pragma omp flush
	}

	randomInnerPoints.clear();

	for (unsigned j = 0; j < allEnclosedPoints.size(); ++j) {
		for (vtkIdType i = 0; i < points[j]->GetNumberOfPoints(); i++) {
			if (allEnclosedPoints[j]->IsInside(i)) {
				double *pd = points[j]->GetPoint(i);
				point p = { pd[0], pd[1], pd[2] };
				randomInnerPoints.push_back(p);
			}
		}
	}

}

point SimpleDomain2D::getRandomPoint() {

	if (randomInnerPoints.empty())
		generateRandomPoints();
	point p = randomInnerPoints.front();
	randomInnerPoints.pop_front();

	++pointCounter;

	return p;
}

/*
 * Return if a segment of two inside-domain points is completely inside the domain.
 */
int SimpleDomain2D::isSegmentInside(point xs, point xf) {
	double tol = 1.e-8;

	if(isConvexDomain)
		return true;

	vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int hasIntersected = abs(locator->IntersectWithLine(xs.p, xf.p, points, ids));

	return !hasIntersected;
}

double SimpleDomain2D::getCharacteristicLength() {
	return characteristicLength;
}

double SimpleDomain2D::getSize() {
	if (volume == 0.0) {
		vtkMassProperties *massProperty = vtkMassProperties::New();
		massProperty->SetInputData(vtkGeometry);
		massProperty->Update();
		volume = massProperty->GetVolume();
	}
	return volume;
}

vtkSmartPointer<vtkOBBTree>& SimpleDomain2D::getLocator() {
	return locator;
}

int SimpleDomain2D::getDraw() {
	return nDraw;
}

deque<point>& SimpleDomain2D::getRandomInnerPoints() {
	return randomInnerPoints;
}

vtkSmartPointer<vtkPolyData>& SimpleDomain2D::getVtkGeometry() {
	return vtkGeometry;
}

double SimpleDomain2D::getDLim(long long int nVessels, double factor) {
	return characteristicLength * sqrt(factor / nVessels);
}

double* SimpleDomain2D::getLocalNeighborhood(point p, long long int nVessels) {
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

int SimpleDomain2D::getSeed()
{
	return this->seed;
}

string SimpleDomain2D::getFilename()
{
	return this->filename;
}

void SimpleDomain2D::logDomainFiles(FILE *fp) {
	fprintf(fp, "SimpleDomain2D\n");
    fprintf(fp, "filename = %s\n", this->getFilename().c_str());
}

vtkSmartPointer<vtkSelectEnclosedPoints> SimpleDomain2D::getEnclosedPoints() {
	vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	enclosedPoints->Initialize(this->vtkGeometry);
	return enclosedPoints;
}