/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * GeneratorDataMonitor.cpp
 *
 *  Created on: Mar 12, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#include "GeneratorDataMonitor.h"

GeneratorDataMonitor::GeneratorDataMonitor(AbstractDomain *domain)
{
	this->domain = domain;

	dLimObservations = 5;
	dLimOcurrencies = deque<double>(dLimObservations,-1.0);
}

void GeneratorDataMonitor::addDLimValue(double value, int nVessels){
	dLimOcurrencies.push_back(value / domain->getDLim(nVessels,domain->getInstanceData()->perfusionAreaFactor) );
	dLimOcurrencies.pop_front();
}

void GeneratorDataMonitor::setDLimObservations(int nObs){
	dLimObservations = nObs;
	dLimOcurrencies = deque<double>(dLimObservations,-1.0);
}

void GeneratorDataMonitor::update(){

	if(dLimOcurrencies.front() > 0){
		double maxDLim = 0.0;
		for(std::deque<double>::iterator it = dLimOcurrencies.begin(); it != dLimOcurrencies.end(); ++it) {
			if(maxDLim<*it)
				maxDLim = *it;
		}

		domain->getInstanceData()->dLimCorrectionFactor = maxDLim;
	}
}

void GeneratorDataMonitor::reset(){
	for(std::deque<double>::iterator it = dLimOcurrencies.begin(); it != dLimOcurrencies.end(); ++it) {
		*it = -1.0;
	}

}
