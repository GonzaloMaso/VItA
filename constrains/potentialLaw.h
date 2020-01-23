/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * exponentialLaw.h
 *
 *  Created on: Jun 1, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef POTENTIALLAW_H_
#define POTENTIALLAW_H_

/**
 * Returns the same value for the potential law (Murray law).
 * @param treeLevel Bifurcation level.
 * @return Value for the exponential law.
 */
double constantPotentialLaw(int treeLevel);
/**
 * Returns two different values for the exponential law (Murray law), depending on the bifurcation level.
 * @param treeLevel Bifurcation level.
 * @return Value for the exponential law.
 */
double twoLevelPotentialLaw(int treeLevel);

inline double constantPotentialLaw(int treeLevel){
	return 3;
}

inline double twoLevelPotentialLaw(int treeLevel){
	if(treeLevel<6)
		return 2.55;
	else
		return 3.0;
}

#endif /* POTENTIALLAW_H_ */
