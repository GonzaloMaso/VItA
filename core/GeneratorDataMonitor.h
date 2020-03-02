/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * GeneratorDataMonitor.h
 *
 *  Created on: Mar 12, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef GENERATORDATAMONITOR_H_
#define GENERATORDATAMONITOR_H_

#include <deque>

#include "../structures/domain/AbstractDomain.h"
#include "GeneratorData.h"

using namespace std;

/**
 * Monitors the GeneratorData in order to adjust parameters or take other actions.
 */
class GeneratorDataMonitor {
	/**	Deque of the last @p observations ocurrences of Dlim values. */
	deque<double> dLimOcurrencies;
	/** Amount of observations stored for the dLim value. */
	int dLimObservations;
	/**	Perfusion domain. */
	AbstractDomain *domain;

public:
	/**
	 * Constructor.
	 */
	GeneratorDataMonitor(AbstractDomain *domain);

	/**
	 * Adds a new ocurrence of dLim.
	 * @param value New ocurrence of dLim.
	 * @param nVessels Amount of terminals in the domain.
	 */
	void addDLimValue(double value, int nVessel);
	/**
	 * Sets the amount of observations of dLim values that are stored in @p dLimOcurrencies.
	 * @param nObs Amount of observations of dLim values that are stored in @p dLimOcurrencies.
	 */
	void setDLimObservations(int nObs);

	/**
	 * Perform adjustments over GeneratorData based on the recorded history of the parameters.
	 */
	void update();

	/**
	 * Resets the recorded history of the parameters.
	 */
	void reset();
};

#endif /* GENERATORDATAMONITOR_H_ */
