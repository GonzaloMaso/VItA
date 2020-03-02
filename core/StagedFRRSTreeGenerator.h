/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * StagedFRRSTreeGenerator.h
 *
 *  Created on: Jun 2, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STAGEDFRRSTREEGENERATOR_H_
#define STAGEDFRRSTREEGENERATOR_H_

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "../constrains/AbstractConstraintFunction.h"
#include "../structures/domain/IDomainObserver.h"
#include "../structures/domain/StagedDomain.h"
#include "../structures/tree/AbstractStructuredCCOTree.h"

#include "GeneratorDataMonitor.h"

using namespace std;
/**
 * Generator for a single fixed perfusion tree with many stages.
 */
class StagedFRRSTreeGenerator : public IDomainObserver {

	/** Wrapper with parameters associated to a tree generation process. */
	GeneratorData *instanceData;
	/**	Monitor of the @p instanceData . */
	GeneratorDataMonitor *dataMonitor;

	/** Amount of terminals in the trees.*/
	long long int nTerminals;
	/** Perfusion volume.*/
	double dLim;

	/** Perfusion domain. */
	StagedDomain *domain;

	/** Generated trees. */
	AbstractStructuredCCOTree *tree;

	vector<AbstractConstraintFunction<double,int> *> gams;
	vector<AbstractConstraintFunction<double,int> *> epsLims;
	vector<AbstractConstraintFunction<double,int> *> nus;

	int stage;

	/**	If the current generation saves the configuration used.*/
	int isGeneratingConfFile;
	/**	File name for the configuration file.*/
	string confFilename;

public:
	/**
	 * Constructor for non-optimized and constant viscosity tree.
	 * @param domain	Perfusion domain.
	 * @param xi		Root proximal position.
	 * @param rootRadii	Root radius.
	 * @param qi		Inflow for the perfusion domain.
	 * @param nTerm		Total number of tree terminals.
	 * @param gam		Murray's power law function.
	 * @param epsLim	Symmetry constraint function.
	 * @param nu		Viscosity function.
	 * @param minAngle	Minimum allowed bifurcation angle.
	 */
	StagedFRRSTreeGenerator(StagedDomain *domain, point xi, double rootRadii, double qi, long long int nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu, double minAngle, double refPressure);
	/**
	 * Constructor for tree with Fahraeus-Lindquist viscosity model.
	 * @param domain	Perfusion domain.
	 * @param xi		Root proximal position.
	 * @param rootRadii	Root radius.
	 * @param qi		Inflow for the perfusion domain.
	 * @param nTerm		Total number of tree terminals.
	 * @param gam		Murray's power law function.
	 * @param epsLim	Symmetry constraint function.
	 * @param nu		Viscosity function.
	 * @param minAngle	Minimum allowed bifurcation angle.
	 * @param viscosityTolerance	Convergence tolerance for the absolute beta variation in the viscosity iterative scheme.
	 * @param optimized	If it is optimized with partial tree scaling.
	 */
	StagedFRRSTreeGenerator(StagedDomain *domain, point xi, double rootRadii, double qi, long long int nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu, double minAngle, double refPressure, double viscosityTolerance, int optimized);
	/**
	 * Constructor to resume from a pre-existent tree.
	 * @param domain	Perfusion domain.
	 * @param tree		Pre-existent tree.
	 * @param nTerm		Total number of tree terminals.
	 * @param gam		Murray's power law function.
	 * @param epsLim	Symmetry constraint function.
	 * @param nu		Viscosity function.
	 */
	StagedFRRSTreeGenerator(StagedDomain *domain, AbstractStructuredCCOTree *tree, long long int nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu);
	/**
	 * Generates the specified tree.
	 * @return	Perfusion tree.
	 */
	AbstractStructuredCCOTree * generate();
	/**
	 * Resumes the tree generation.
	 * @return	Perfusion tree.
	 */
	AbstractStructuredCCOTree *resume();
	/**
	 * Returns the perfusion domain.
	 * @return Perfusion domain.
	 */
	StagedDomain* getDomain();
	/**
	 * Returns the generated tree.
	 * @return Generated tree.
	 */
	vector<AbstractStructuredCCOTree*> getTrees();
	/**
	 * Enables the configuration file generation capabilities.
	 * @param filename	File name where the generator stores the configuration data.
	 */
	void enableConfigurationFile(string filename);
	/**
	 * Updates internal data from this instance after the domain has been modified.
	 * @param observableInstance Domain instance modified.
	 */
	void observableModified(IDomainObservable * observableInstance);

protected:
	/**	Configuration file stream. */
	ofstream confFile;
	/** Initial timestamp of the generation process. */
	chrono::steady_clock::time_point beginningTime;
	/**
	 * Returns if the length of the segment defined by xProx and xNew is higher than dLim (distance criterion).
	 * @param xNew	Proposed distal point.
	 * @param iTry	Number of trial.
	 * @return If the root segment is valid.
	 */
	int isValidRootSegment(point xNew, int iTry);
	/**
	 * Returns if the closest point to the whole try is greater than dLim (distance criterion for a new segment).
	 * @param xNew Proposed distal point.
	 * @param iTry Number of trial.
	 * @return If the segment is valid.
	 */
	int isValidSegment(point xNew, int iTry);
	/**
	 * Generates the configuration file for the current tree generation.
	 * @param mode Is the openmode used for the generated file (ios::out for generation, ios::app for resume).
	 */
	void generatesConfigurationFile(ios::openmode mode);
	/**
	 * Saves the timestamp for event @p label.
	 */
	void markTimestampOnConfigurationFile(string label);
	/**
	 * Closes the configuration file for the current tree generation.
	 */
	void closeConfigurationFile();

};

#endif /* STAGEDFRRSTREEGENERATOR_H_ */
