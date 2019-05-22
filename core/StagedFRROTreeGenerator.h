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

#include "../constrains/AbstractConstraintFunction.h"
#include "../structures/domain/IDomainObserver.h"
#include "../structures/domain/StagedDomain.h"
#include "../structures/tree/AbstractObjectCCOTree.h"

#include "GeneratorDataMonitor.h"

using namespace std;
/**
 * Generator for a single fixed perfusion tree with many stages.
 */
class StagedFRROTreeGenerator : public IDomainObserver {

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
	AbstractObjectCCOTree *tree;

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
	StagedFRROTreeGenerator(StagedDomain *domain, point xi, double rootRadii, double qi, long long int nTerm, vector<AbstractConstraintFunction<double,int> *>gam, vector<AbstractConstraintFunction<double,int> *>epsLim, vector<AbstractConstraintFunction<double,int> *>nu, double minAngle, double refPressure, double viscosityTolerance);
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
	AbstractObjectCCOTree * generate(long long int saveInterval, string tempDirectory);
	/**
	 * Resumes the tree generation.
	 * @param saveInterval Number of iterations performed between saved steps.
	 * @param tempDirectory Directory where intermediate solutions are saved.
	 * @return	Perfusion tree.
	 */
	AbstractObjectCCOTree *resume(long long int saveInterval, string tempDirectory);
	/**
	 * Returns the perfusion domain.
	 * @return Perfusion domain.
	 */
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

#endif /* STAGEDFRROTREEGENERATOR_H_ */