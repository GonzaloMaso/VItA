#include "PowerCostEstimator.h"
#include "../vascularElements/SingleVessel.h"
#include <math.h>

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

PowerCostEstimator::PowerCostEstimator(unsigned int lExp, unsigned int rExp): AbstractCostEstimator(){
	this->previousCost = 0.0;
    this->lExp = lExp;
    this->rExp = rExp;
}

PowerCostEstimator::~PowerCostEstimator(){
}

void PowerCostEstimator::previousState(AbstractObjectCCOTree* tree, AbstractVascularElement* parent, point iNew, point iTest, double dLim){
    this->previousCost = computeTreeCost(tree->getRoot());
}

double PowerCostEstimator::computeCost(AbstractObjectCCOTree* tree){
	return this->computeTreeCost(tree->getRoot()) - this->previousCost;
}

double PowerCostEstimator::computeTreeCost(AbstractVascularElement* root) {
	double currentCost = this->computeVesselCost(static_cast<SingleVessel *>(root));
	vector<AbstractVascularElement *> children = root->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		currentCost += computeTreeCost(*it);
	}
	return currentCost;
}

AbstractCostEstimator* PowerCostEstimator::clone() {
	return (new PowerCostEstimator(this->lExp, this->rExp));
}

double PowerCostEstimator::computeVesselCost(SingleVessel *vessel) {
    return fast_pow(vessel->length, this->lExp) * fast_pow(vessel->radius, this->rExp);
}

void PowerCostEstimator::logCostEstimator(FILE *fp) {
    fprintf(fp, "This domain uses PowerCostEstimator.\n");
    fprintf(fp, "Lenght exponent = %u.\n", this->lExp);
    fprintf(fp, "Radius exponent = %u.\n", this->rExp);
}