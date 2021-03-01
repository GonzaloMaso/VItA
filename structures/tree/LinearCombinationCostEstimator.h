#ifndef TREE_LINEARCOMBINATIONCOSTESTIMATOR_H_
#define TREE_LINEARCOMBINATIONCOSTESTIMATOR_H_

#include<vector>
#include "AbstractCostEstimator.h"

#include "../vascularElements/AbstractVascularElement.h"

/**
 * Cost estimator that computes the tree volume.
 */
class LinearCombinationCostEstimator: public AbstractCostEstimator {
	vector<double> alpha;
	vector<unsigned int> lExp, rExp;
	/**	Cost at the previous step. */
	double previousCost;
public:
	/**
	 * Common constructor.
	 */
	LinearCombinationCostEstimator(vector<double> alpha, vector<unsigned int> lExp, vector<unsigned int> rExp);
	/**
	 * Common destructor.
	 */
	virtual ~LinearCombinationCostEstimator();
	/**
	 * Clones the current estimator instance.
	 * @return Cloned instance.
	 */
	AbstractCostEstimator *clone() override;

	/**
	 * Extracts information of the tree at the previous step.
	 * @param root	Root of the tree at the previous step.
	 * @param parent	Vascular element where the new vessel will be connected.
	 * @param iNew	Distal position of the new vessel.
	 * @param iTest	Proximal position of the new vessel.
	 * @param dLim	Minimum radius distance from the new vessel to the tree.
	 */
	void previousState(AbstractObjectCCOTree *tree, AbstractVascularElement *parent, point iNew, point iTest, double dLim) override;

	/**
	 * Computes the functional cost of the given tree.
	 * @param tree	Tree at the current step.
	 * @return Cost of the given tree.
	 */
	double computeCost(AbstractObjectCCOTree *tree) override;

	void logCostEstimator(FILE *fp) override;

private:
	/**
	 * Computes the volume for the tree with root @p root.
	 * @param root	Root of the tree.
	 * @return	Volume of the tree.
	 */
	double computeTreeCost(AbstractVascularElement* root);
	/**
     * Computes the cost of a single vessel.
     */
    double computeVesselCost(SingleVessel *vessel);
};

#endif /* TREE_LINEARCOMBINATIONCCOSTESTIMATOR_H_ */
