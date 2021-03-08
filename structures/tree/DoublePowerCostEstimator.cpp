#include "DoublePowerCostEstimator.h"
#include "../vascularElements/SingleVessel.h"
#include <cmath>

DoublePowerCostEstimator::DoublePowerCostEstimator(double lExp, double rExp): AbstractCostEstimator(){
	this->previousCost = 0.0;
    this->lExp = lExp;
    this->rExp = rExp;
}

DoublePowerCostEstimator::~DoublePowerCostEstimator(){
}

void DoublePowerCostEstimator::previousState(AbstractObjectCCOTree* tree, AbstractVascularElement* parent, point iNew, point iTest, double dLim){
    this->previousCost = computeTreeCost(tree->getRoot());
}

double DoublePowerCostEstimator::computeCost(AbstractObjectCCOTree* tree){
	return this->computeTreeCost(tree->getRoot()) - this->previousCost;
}

double DoublePowerCostEstimator::computeTreeCost(AbstractVascularElement* root) {
	double currentCost = this->computeVesselCost(static_cast<SingleVessel *>(root));
	vector<AbstractVascularElement *> children = root->getChildren();
	for (std::vector<AbstractVascularElement *>::iterator it = children.begin(); it != children.end(); ++it) {
		currentCost += computeTreeCost(*it);
	}
	return currentCost;
}

AbstractCostEstimator* DoublePowerCostEstimator::clone() {
	return (new DoublePowerCostEstimator(this->lExp, this->rExp));
}

double DoublePowerCostEstimator::computeVesselCost(SingleVessel *vessel) {
    return pow(vessel->length, this->lExp) * pow(vessel->radius, this->rExp);
}

void DoublePowerCostEstimator::logCostEstimator(FILE *fp) {
    fprintf(fp, "This domain uses PowerCostEstimator.\n");
    fprintf(fp, "Lenght exponent = %u.\n", this->lExp);
    fprintf(fp, "Radius exponent = %u.\n", this->rExp);
}