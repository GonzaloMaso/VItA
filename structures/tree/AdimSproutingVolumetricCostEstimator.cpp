/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SproutingVolumetricCostEstimator.cpp
 *
 *  Created on: Apr 11, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "AdimSproutingVolumetricCostEstimator.h"
#include "../vascularElements/SingleVessel.h"
#include <math.h>

AdimSproutingVolumetricCostEstimator::AdimSproutingVolumetricCostEstimator(double volumeFactor, double proteolyticFactor, double diffusionFactor, double volumeRef, double radiusRef): AbstractCostEstimator(){
	this->previousVolume = 0.0;
	this->proteolyticFactor = proteolyticFactor;
	this->diffusionFactor = diffusionFactor;
	this->volumeFactor = volumeFactor;
	this->distToParent = 1.0;
	this->parentRadius = 1.0;
	this->volumeRef = volumeRef;
	this->lengthRef = cbrt(3*volumeRef/(4*M_PI));
	this->radiusRef = radiusRef;
}

AdimSproutingVolumetricCostEstimator::~AdimSproutingVolumetricCostEstimator(){
}

void AdimSproutingVolumetricCostEstimator::previousState(AbstractObjectCCOTree* tree, AbstractVascularElement* parent, point iNew, point iTest, double dLim){
	previousVolume = ((SingleVessel *) tree->getRoot())->treeVolume;

	point a = ((SingleVessel *)parent)->xProx;
	point b = ((SingleVessel *)parent)->xDist;
	//	Parent-to-iNew distance
	//	Parent vessel slope
	point m = b - a;
	//	Parameter for closer projection
	double t = (m ^ (iNew - a)) / (m^m);
	//	Confine t into [0,1] interval
	if (t < 0){
		t = 0;
	}
	else if( t > 1.0){
		t = 1.0;
	}
	//	Closest segment between iNew and parent vessel
	point proj = (iNew - a) - m * t;
	distToParent = sqrt(proj ^ proj);

	parentRadius = ((SingleVessel *)parent)->radius;

	this->dLim = dLim;
	this->bifLevel = ((SingleVessel *)parent)->nLevel;
}

double AdimSproutingVolumetricCostEstimator::computeCost(AbstractObjectCCOTree* tree){
	double volCost = volumeFactor * (computeTreeCost(tree->getRoot()) - previousVolume) / volumeRef;
	double proteolysisCost = proteolyticFactor * parentRadius / radiusRef; // 500.0
	double parentLengthRatio = distToParent / lengthRef;
	double stimulusCost = diffusionFactor * parentLengthRatio * parentLengthRatio;
//	cout << "Volume ref = " << volumeRef << " - Radius ref = " << radiusRef << " - Length ref = " << lengthRef << endl;
//	cout << "Volumetric cost = " << volCost << ", Protease degradation cost = " << proteolysisCost << ", VEGF/FGF difussion cost = " << stimulusCost << endl;
	return volCost + proteolysisCost + stimulusCost ;
}

double AdimSproutingVolumetricCostEstimator::computeTreeCost(AbstractVascularElement* root) {
	double currentCost = ((SingleVessel *)root)->getVolume();
	vector<AbstractVascularElement *> children = root->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		currentCost += computeTreeCost(*it);
	}
	return currentCost;
}

AbstractCostEstimator* AdimSproutingVolumetricCostEstimator::clone(){
	return (new AdimSproutingVolumetricCostEstimator(volumeFactor, proteolyticFactor, diffusionFactor, volumeRef, radiusRef));
}

double AdimSproutingVolumetricCostEstimator::getVolumeFactor()
{
	return this->volumeFactor;
}

double AdimSproutingVolumetricCostEstimator::getProteolyticFactor()
{
	return this->proteolyticFactor;
}

double AdimSproutingVolumetricCostEstimator::getDiffusionFactor()
{
	return this->diffusionFactor;
}

double AdimSproutingVolumetricCostEstimator::getVolumeRef()
{
	return this->volumeRef;
}

double AdimSproutingVolumetricCostEstimator::getRadiusRef()
{
	return this->radiusRef;
}

double AdimSproutingVolumetricCostEstimator::getLenghtRef()
{
	return this->lengthRef;
}

void AdimSproutingVolumetricCostEstimator::logCostEstimator(FILE *fp) {
	double v_fac, p_fac, d_fac, v_ref, r_ref, l_ref;
	v_fac = this->getVolumeFactor();
    p_fac = this->getProteolyticFactor();
    d_fac = this->getDiffusionFactor();
    v_ref = this->getVolumeRef();
    r_ref = this->getRadiusRef();
    l_ref = this->getLenghtRef();
    fprintf(fp, "This domain uses AdimSproutingVolumetricCost.\n");
    fprintf(fp, "Volume factor = %f.\n", v_fac);
    fprintf(fp, "Proteolytic factor = %f.\n", p_fac);
    fprintf(fp, "Diffusion factor = %f.\n", d_fac);
    fprintf(fp, "Reference volume = %f.\n", v_ref);
    fprintf(fp, "Reference radius = %f.\n", r_ref);
    fprintf(fp, "Reference lenght = %f.\n", l_ref);
}