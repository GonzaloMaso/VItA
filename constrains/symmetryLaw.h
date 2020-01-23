/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * symmetryLaw.h
 *
 *  Created on: Jun 1, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef SYMMETRYLAW_H_
#define SYMMETRYLAW_H_

/**
 * Returns zero for all bifurcation levels.
 * @param treeLevel Bifurcation level.
 * @return Symmetry coefficient.
 */
double constantNullBifurcationSymmetry(int treeLevel);
/**
 * Returns the same bifurcation symmetry coeficient for all bifurcation levels.
 * @param treeLevel Bifurcation level.
 * @return	Symmetry coefficient.
 */
double constantBifurcationSymmetry(int treeLevel);
/**
 * Returns a different bifurcation symmetry coeficient depending to the current bifurcation levels.
 * @param treeLevel Current bifurcation level.
 * @return	Symmetry coefficient.
 */
double variableBifurcationSymmetry(int treeLevel);

inline double constantNullBifurcationSymmetry(int treeLevel){
	return 0.;
}

inline double constantBifurcationSymmetry(int treeLevel){
	return 0.3;
}

inline double blockBifurcationSymmetry(int treeLevel){
	if(treeLevel<6)
		return 2.0;	//
	else
		return 0.0;
}

inline double variableBifurcationSymmetry(int treeLevel){
	if(treeLevel<6)
		return 0.4;	//
	else
		return 0.0;
}

#endif /* SYMMETRYLAW_H_ */
