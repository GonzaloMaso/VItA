/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * SproutingVolumetricCostEstimator.h
 *
 *  Created on: Apr 11, 2018
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef TREE_ADIMSPROUTINGVOLUMETRICCOSTESTIMATOR_H_
#define TREE_ADIMSPROUTINGVOLUMETRICCOSTESTIMATOR_H_

#include "AbstractCostEstimator.h"

/**
 * Cost estimator that considers diffusion and vessel-wall degradation associated to
 * sprouting angiogenesis and also weighs the variation of the tree volume. Wall degradation is assumed linear w.r.t. vessel radius
 * (under the assumption that vessel thickness is linear with vessel radius) and media diffusion of biosignals is assumed to follow
 * Fick's 2nd law.
 */
class AdimSproutingVolumetricCostEstimator: public AbstractCostEstimator {
	/**	Volume at the previous step. */
	double previousVolume;
	/**	Vessel wall proteolytic degradation factor. */
	double proteolyticFactor;
	/**	Factor of media diffusion of angiogenic growth factors. */
	double diffusionFactor;
	/**	Weigh for tree volume. */
	double volumeFactor;
	/**	Distance between sprouting source and candidate parent vessel. */
	double distToParent;
	/**	Radius of the candidate parent vessel. */
	double parentRadius;
	/**	Minimum radius distance from the new vessel to the tree. */
	double dLim;
	/** Level of bifurcation. */
	double bifLevel;
	/**	Volume of the vascularized domain */
	double volumeRef;
	/**	Root radius of the vascular tree */
	double radiusRef;
	/**	Equivalent radius associated to @p vRef */
	double lengthRef;

public:
	/**
	 * Common constructor.
	 */
	AdimSproutingVolumetricCostEstimator(double volumeFactor, double proteolyticFactor, double diffusionFactor, double volumeRef, double radiusRef);
	/**
	 * Common destructor.
	 */
	virtual ~AdimSproutingVolumetricCostEstimator();
	/**
	 * Clones the current estimator instance.
	 * @return Cloned instance.
	 */
	AbstractCostEstimator *clone();

	/**
	 * Extracts information of the tree at the previous step.
	 * @param root	Root of the tree at the previous step.
	 * @param parent	Vascular element where the new vessel will be connected.
	 * @param iNew	Distal position of the new vessel.
	 * @param iTest	Proximal position of the new vessel.
	 * @param dLim	Minimum radius distance from the new vessel to the tree.
	 */
	void previousState(AbstractObjectCCOTree *tree, AbstractVascularElement *parent, point iNew, point iTest, double dLim);
	/**
	 * Computes the functional cost of the given tree.
	 * @param root	Root of the tree at the current step.
	 * @return Cost of the given tree.
	 */
	double computeCost(AbstractObjectCCOTree* tree);
	/**
	 * @return @p volumeFactor
	 */
	double getVolumeFactor();
	/**
	 * @return @p proteolyticFactor
	 */
	double getProteolyticFactor();
	/**
	 * @return @p diffusionFactor
	 */
	double getDiffusionFactor();
	/**
	 * @return @p volumeRef
	 */
	double getVolumeRef();
	/**
	 * @return @p radiusRef
	 */
	double getRadiusRef();
	/**
	 * @return @p lenghtRef
	 */
	double getLenghtRef();

	void logCostEstimator(FILE *fp);
	
private:
	/**
	 * Computes the volume for the tree with root @p root.
	 * @param root	Root of the tree.
	 * @return	Volume of the tree.
	 */
	double computeTreeCost(AbstractVascularElement* root);

};

#endif /* TREE_SPROUTINGVOLUMETRICCOSTESTIMATOR_H_ */
