#include "LinearCombinationCostEstimator.h"
#include "../vascularElements/SingleVessel.h"
#include<vector>

static double fast_pow(double base,unsigned int exponent) {
    double result = 1;
    while (exponent) {
        if (exponent & 1) {
            result *= base;
            exponent -= 1;
        }
        else {
            result *= base * base;
            exponent -=2;
        }
    }
    return result;
}

LinearCombinationCostEstimator::LinearCombinationCostEstimator(vector<double> alpha, vector<unsigned int> lExp, vector<unsigned int> rExp) : AbstractCostEstimator(){
	this->previousCost = 0.0;
	this->alpha = alpha;
	this->lExp = lExp;
	this->rExp = rExp;
}

LinearCombinationCostEstimator::~LinearCombinationCostEstimator(){
}

double LinearCombinationCostEstimator::computeCost(AbstractObjectCCOTree *tree){
	return computeTreeCost(tree->getRoot()) - this->previousCost;
}

void LinearCombinationCostEstimator::previousState(AbstractObjectCCOTree *tree, AbstractVascularElement* parent, point iNew, point iTest, double dLim){
	this->previousCost = this->computeTreeCost(tree->getRoot());
}

double LinearCombinationCostEstimator::computeTreeCost(AbstractVascularElement* root) {

	double currentCost = this->computeVesselCost(static_cast<SingleVessel *>(root));
	vector<AbstractVascularElement *> children = root->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		currentCost += computeTreeCost(*it);
	}
	return currentCost;
}

AbstractCostEstimator* LinearCombinationCostEstimator::clone(){
	return (new LinearCombinationCostEstimator(this->alpha, this->lExp, this->rExp));
}

void LinearCombinationCostEstimator::logCostEstimator(FILE *fp) {	
	fprintf(fp, "This domain uses LinearCombinationCostEstimator.\n");
	fprintf(fp, "Alpha: ");
	for (size_t i = 0; i < this->alpha.size(); ++i) {
		fprintf(fp, "%lf ", this->alpha[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "lExp: ");
	for (size_t i = 0; i < this->lExp.size(); ++i) {
		fprintf(fp, "%u ", this->lExp[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "rExp: ");
	for (size_t i = 0; i < this->rExp.size(); ++i) {
		fprintf(fp, "%u ", this->rExp[i]);
	}
	fprintf(fp, "\n");
}

double LinearCombinationCostEstimator::computeVesselCost(SingleVessel *vessel) {
	double cost = 0;
	for (size_t i = 0; i < this->alpha.size(); ++i) {
		cost += this->alpha[i] * fast_pow(vessel->length, this->lExp[i]) * fast_pow(vessel->radius, this->rExp[i]);
	}
	return cost;
}