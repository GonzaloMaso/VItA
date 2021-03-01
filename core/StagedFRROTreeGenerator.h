/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * FixedPerfusionRadiusTreeGenerator.h
 *
 *  Created on: Jun 2, 2017
 *      Author: Gonzalo D. Maso Talou
 */

#ifndef STAGEDFRROTREEGENERATOR_H_
#define STAGEDFRROTREEGENERATOR_H_

#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include<ctime>

#include "../constrains/AbstractConstraintFunction.h"
#include "../io/task/AbstractSavingTask.h"
#include "../structures/domain/IDomainObserver.h"
#include "../structures/domain/StagedDomain.h"
#include "../structures/tree/AbstractObjectCCOTree.h"
#include "../utils/MemoryMonitor.h"

#include "GeneratorDataMonitor.h"

using namespace std;
/**
 * Generator for a single fixed perfusion tree with many stages.
 */
class StagedFRROTreeGenerator : public IDomainObserver {

	/** Time of the beggining of the tree generation process*/
	time_t beginTime;
	/** Time of the end of the tree generation process*/
	time_t endTime;
	/** Wrapper with parameters associated to a tree generation process. */
	GeneratorData *instanceData;
	/**	Monitor of the @p instanceData . */
	GeneratorDataMonitor *dataMonitor;
	/**	Monitors the memory usage. */
	MemoryMonitor *monitor;
	/**	Action executed at each save interval during generate and resume methods */
	vector<AbstractSavingTask *> savingTasks;

	/** Amount of terminals in the trees.*/
	long long int nTerminals;
	/** Perfusion volume.*/
	double dLim;
	/** Initial dLim value*/
	double dLimInitial;
	/** Final dLim value*/
	double dLimLast;
	/** Perfusion domain. */
	StagedDomain *domain;

	/** Generated trees. */
	AbstractObjectCCOTree *tree;

	vector<AbstractConstraintFunction<double,int> *> gams;
	vector<AbstractConstraintFunction<double,int> *> epsLims;
	vector<AbstractConstraintFunction<double,int> *> nus;

	/**	Current stage of generation.*/
	int stage;
	/**	If the current generation saves the configuration used.*/
	int isGeneratingConfFile;
	/**	File name for the configuration file.*/
	string confFilename;

	bool didAllocateTree;

public:
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
	 */
	StagedFRROTreeGenerator(StagedDomain *domain, point xi, double rootRadii, double qi, long long int nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu, double refPressure, double viscosityTolerance);
	/**
	 * Constructor to resume from a pre-existent tree.
	 * @param domain	Perfusion domain.
	 * @param tree		Pre-existent tree.
	 * @param nTerm		Total number of tree terminals.
	 * @param gam		Murray's power law function.
	 * @param epsLim	Symmetry constraint function.
	 * @param nu		Viscosity function.
	 */
	StagedFRROTreeGenerator(StagedDomain *domain, AbstractObjectCCOTree *tree, long long int nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu);
	/**
	 * Common destructor.
	 */
	~StagedFRROTreeGenerator();
	/**
	 * Generates the specified tree.
	 * @param saveInterval Number of iterations performed between saved steps.
	 * @param tempDirectory Directory where intermediate solutions are saved.
	 * @return	Perfusion tree.
	 */
	AbstractObjectCCOTree *generate(long long int saveInterval, string tempDirectory);
	/**
	 * Resumes the tree generation.
	 * @param saveInterval Number of iterations performed between saved steps.
	 * @param tempDirectory Directory where intermediate solutions are saved.
	 * @return	Perfusion tree.
	 */
	AbstractObjectCCOTree *resume(long long int saveInterval, string tempDirectory);
	/**
	 * Generates the specified tree with a failsafe for bad configuration.
	 * @param saveInterval Number of iterations performed between saved steps.
	 * @param tempDirectory Directory where intermediate solutions are saved.
	 * @param maxNumOfTrials Maximum number of trials before generator fails and exits.
	 * @return	Perfusion tree.
	 */
	AbstractObjectCCOTree *generateExperimental(long long int saveInterval, string tempDirectory, int maxNumOfTrials);
	/**
	 * Resumes the tree generation.
	 * @param saveInterval Number of iterations performed between saved steps.
	 * @param tempDirectory Directory where intermediate solutions are saved.
	 * @param maxNumOfTrials Maximum number of trials before generator fails and exits.
	 * @return	Perfusion tree.
	 */
	AbstractObjectCCOTree *resumeExperimental(long long int saveInterval, string tempDirectory, int maxNumOfTrials);
	/**
	 * Returns the perfusion domain.
	 * @return Perfusion domain.
	 */

	/**
	 * Resumes the tree generation process and saves the optimal xNew and xBif in @param fp.
	 */
	AbstractObjectCCOTree *resumeSavePoints(long long int saveInterval, string tempDirectory, FILE *fp);

	AbstractObjectCCOTree *resumeSavePointsMidPoint(long long int saveInterval, string tempDirectory, FILE *fp);
	
	StagedDomain* getDomain();
	/**
	 * Returns the generated tree.
	 * @return Generated tree.
	 */
	vector<AbstractObjectCCOTree *> getTrees();
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
	/**
	 * Return the currently generated @p tree
	 * @return Generated tree.
	 */
	AbstractObjectCCOTree*& getTree();
	/**
	 * Saves the current generation status using the @p savingTasks and the memory monitor.
	 * @param terminals	Current iteration number.
	 */
	void saveStatus(long long int terminals);
	/**
	 * Sets a set of saving tasks to be produced each time the generation process saves.
	 * @param savingTasks	Set of tasks to be performed each time that the tree is saved.
	 */
	void setSavingTasks(const vector<AbstractSavingTask*>& savingTasks);
	/**
	 * Returns a pointer to the vector of gams.
	 */
	vector<AbstractConstraintFunction<double, int> *>* getGams();
	/**
	 * Returns a pointer to the vector of eps_lims.
	 */
	vector<AbstractConstraintFunction<double, int> *>* getEpsLims();
	/**
	 * Returns a pointer to the vector of nus.
	 */
	vector<AbstractConstraintFunction<double, int> *>* getNus();
	/**
	 * Returns the current DLim value.
	 */
	double getDLim();
	/**
	 * Set the DLim value. Use this only to resume process from a previous generation.
	 */
	void setDLim(double newDLim);
	
	time_t getBeginTime();

	time_t getEndTime();

	double getDLimInitial();

	double getDLimLast();
	
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
#endif