/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * Mesh.cpp
 *
 *  Created on: May 24, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#include "SimpleDomain.h"

#include <chrono>
#include <random>
#include <omp.h>

//	Model
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkMassProperties.h>
#include <vtkSmartPointer.h>
#include <vtkSelectEnclosedPoints.h>

#include "UniformDistributionGenerator.h"

SimpleDomain::SimpleDomain(string filename, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filename = filename;
	//	Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->SetTolerance(1E-12);
	locator->BuildLocator();

	nDraw = 10000;

	this->seed = chrono::system_clock::now().time_since_epoch().count();
	double *bb = vtkGeometry->GetBounds();
	distribution = new UniformDistributionGenerator();
	distribution->initialize(this->seed,bb);
	this->didAllocateDistribution = true;

	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);
}

SimpleDomain::SimpleDomain(string filename, int N, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filename = filename;
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->SetTolerance(1E-12);
	locator->BuildLocator();

	nDraw = N;

	this->seed = chrono::system_clock::now().time_since_epoch().count();
	double *bb = vtkGeometry->GetBounds();
	distribution = new UniformDistributionGenerator();
	distribution->initialize(this->seed,bb);
	this->didAllocateDistribution = true;

	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);
}

SimpleDomain::SimpleDomain(string filename, int N, int seed, GeneratorData *instanceData) :
		AbstractDomain(instanceData) {
	this->filename = filename;
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->SetTolerance(1E-12);
	locator->BuildLocator();

	nDraw = N;
	this->seed = seed;

	double *bb = vtkGeometry->GetBounds();
	distribution = new UniformDistributionGenerator();
	distribution->initialize(seed,bb);
	this->didAllocateDistribution = true;

	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);
	cout << "Characteristic length " << characteristicLength << endl;
}

SimpleDomain::SimpleDomain(string filename, int N, int seed, GeneratorData* instanceData, DistributionGenerator* distribution) :
				AbstractDomain(instanceData) {
	this->filename = filename;
	// Read all the data from the file
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
	reader->SetFileName(filename.c_str());
	reader->Update();
	vtkGeometry = reader->GetOutput();

	//	Create the tree
	locator = vtkSmartPointer<vtkOBBTree>::New();
	locator->SetDataSet(vtkGeometry);
	locator->SetTolerance(1E-12);
	locator->BuildLocator();

	nDraw = N;
	this->seed = seed;

	double *bb = vtkGeometry->GetBounds();
	this->distribution = distribution;
	this->distribution->initialize(seed,bb);
	this->didAllocateDistribution = false;

	characteristicLength = max(max((bb[1] - bb[0]) / 2, (bb[3] - bb[2]) / 2), (bb[5] - bb[4]) / 2);
	cout << "Characteristic length " << characteristicLength << endl;
}

SimpleDomain::~SimpleDomain() {
	this->randomInnerPoints.clear();
	if (this->didAllocateDistribution) {
		delete this->distribution;
	}
}

void SimpleDomain::generateRandomPoints() {
	vector<point> newPoints = distribution->getNPoints(nDraw);
	randomInnerPoints.insert(randomInnerPoints.end(), newPoints.begin(), newPoints.end());

	removeRandomOuterPoints();
}

void SimpleDomain::removeRandomOuterPoints() {

	vector<vtkSmartPointer<vtkPoints>> points;
	cout << "Testing inside domain condition for " << randomInnerPoints.size() << " points" << endl;

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
		vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
		pointsPolydata->SetPoints(points[j]);

		//	Points inside test
		vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
		selectEnclosedPoints->SetInputData(pointsPolydata);
		selectEnclosedPoints->SetSurfaceData(vtkGeometry);
		selectEnclosedPoints->SetTolerance(1E-12);
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

point SimpleDomain::getRandomPoint() {

	if (randomInnerPoints.empty())
		generateRandomPoints();
	point p = randomInnerPoints.front();
	randomInnerPoints.pop_front();

	++pointCounter;

	return p;
}

int SimpleDomain::isSegmentInside(point xs, point xf) {

	double tol = 1.e-8;

	if(isConvexDomain)
		return true;

	vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int hasIntersected = abs(locator->IntersectWithLine(xs.p, xf.p, points, ids));

	return !hasIntersected;
}

double SimpleDomain::getCharacteristicLength() {
	return characteristicLength;
}

double SimpleDomain::getSize() {
	if (volume == 0.0) {
		vtkMassProperties *massProperty = vtkMassProperties::New();
		massProperty->SetInputData(vtkGeometry);
		massProperty->Update();
		volume = massProperty->GetVolume();
	}
	return volume;
}

vtkSmartPointer<vtkOBBTree>& SimpleDomain::getLocator() {
	return locator;
}

int SimpleDomain::getDraw() {
	return nDraw;
}

deque<point>& SimpleDomain::getRandomInnerPoints() {
	return randomInnerPoints;
}

vtkSmartPointer<vtkPolyData>& SimpleDomain::getVtkGeometry() {
	return vtkGeometry;
}

double SimpleDomain::getDLim(long long int nVessels, double factor) {
	return characteristicLength * cbrt(factor / nVessels);
}

double* SimpleDomain::getLocalNeighborhood(point p, long long int nVessels) {
	double *localBox = new double[6];
	double neighborhoodRadius = instanceData->closeNeighborhoodFactor * getDLim(nVessels, instanceData->perfusionAreaFactor);

	localBox[0] = p.p[0] - neighborhoodRadius;
	localBox[1] = p.p[0] + neighborhoodRadius;
	localBox[2] = p.p[1] - neighborhoodRadius;
	localBox[3] = p.p[1] + neighborhoodRadius;
	localBox[4] = p.p[2] - neighborhoodRadius;
	localBox[5] = p.p[2] + neighborhoodRadius;

	return localBox;
}

void SimpleDomain::savePoints(string filename) {
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

int SimpleDomain::getSeed()
{
	return this->seed;
}

string SimpleDomain::getFilename()
{
	return this->filename;
}

void SimpleDomain::logDomainFiles(FILE *fp) {
	fprintf(fp, "SimpleDomain\n");
    fprintf(fp, "filename = %s\n", this->getFilename().c_str());
}

vtkSmartPointer<vtkSelectEnclosedPoints> SimpleDomain::getEnclosedPoints() {
	vtkSmartPointer<vtkSelectEnclosedPoints> enclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
	enclosedPoints->Initialize(this->vtkGeometry);
	return enclosedPoints;
}