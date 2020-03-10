// Standard library dependencies
#include<cstdio>

// VItA dependencies
#include"../core/StagedFRROTreeGenerator.h"
#include"../core/GeneratorData.h"
#include"../structures/tree/AbstractCostEstimator.h"
#include"../structures/tree/VolumetricCostEstimator.h"
#include"../structures/tree/SproutingVolumetricCostEstimator.h"
#include"../structures/tree/AdimSproutingVolumetricCostEstimator.h"
#include"../structures/domain/AbstractDomain.h"
#include"../constrains/AbstractConstraintFunction.h"
#include"../structures/tree/AbstractObjectCCOTree.h"
#include"../structures/tree/SingleVesselCCOOTree.h"
#include"../constrains/AbstractConstraintFunction.h"

// This class interface
#include"StagedFRROTreeGeneratorLogger.h"

void log_cost_estimator(FILE *fp, AbstractCostEstimator *cost_estimator)
{
    if ((*cost_estimator).getWhichEstimator() == 0) {
        fprintf(fp, "This domain uses VolumetricCostEstimator.\n");
    }
    else if ((*cost_estimator).getWhichEstimator() == 1) {
        double v_fac, p_fac, d_fac;        
        SproutingVolumetricCostEstimator* sv_cost_estimator = (SproutingVolumetricCostEstimator *) cost_estimator;
        v_fac = (*sv_cost_estimator).getVolumeFactor();
        p_fac = (*sv_cost_estimator).getProteolyticFactor();
        d_fac = (*sv_cost_estimator).getDiffusionFactor();
        fprintf(fp, "This domain uses SproutingVolumetricCostEstimator.\n");
        fprintf(fp, "Volume factor = %f.\n", v_fac);
        fprintf(fp, "Proteolytic factor = %f.\n", p_fac);
        fprintf(fp, "Diffusion factor = %f.\n", d_fac);
    }
    else {
        double v_fac, p_fac, d_fac, v_ref, r_ref, l_ref;
        AdimSproutingVolumetricCostEstimator *asv_cost_estimator = (AdimSproutingVolumetricCostEstimator *) cost_estimator;
        v_fac = (*asv_cost_estimator).getVolumeFactor();
        p_fac = (*asv_cost_estimator).getProteolyticFactor();
        d_fac = (*asv_cost_estimator).getDiffusionFactor();
        v_ref = (*asv_cost_estimator).getVolumeRef();
        r_ref = (*asv_cost_estimator).getRadiusRef();
        l_ref = (*asv_cost_estimator).getLenghtRef();
        fprintf(fp, "This domain uses AdimSproutingVolumetricCost.\n");
        fprintf(fp, "Volume factor = %f.\n", v_fac);
        fprintf(fp, "Proteolytic factor = %f.\n", p_fac);
        fprintf(fp, "Diffusion factor = %f.\n", d_fac);
        fprintf(fp, "Reference volume = %f.\n", v_ref);
        fprintf(fp, "Reference radius = %f.\n", r_ref);
        fprintf(fp, "Reference lenght = %f.\n", l_ref);
    }
}

void log_gen_data(FILE *fp, GeneratorData *data)
{   
    log_cost_estimator(fp, (*data).costEstimator);
    fprintf(fp, "n_level_test = %d.\n", (*data).nLevelTest);
    fprintf(fp, "n_terminal_trial = %d.\n", (*data).nTerminalTrial);
    fprintf(fp, "d_lim_reduction_factor = %f.\n", (*data).dLimReductionFactor);
    fprintf(fp, "perfusion_area_factor = %f.\n", (*data).perfusionAreaFactor);
    fprintf(fp, "close_neighborhood factor = %f.\n", (*data).closeNeighborhoodFactor);
    fprintf(fp, "mid_point_d_lim_factor = %f.\n", (*data).midPointDlimFactor);
    fprintf(fp, "n_bifurcation_test = %d.\n", (*data).nBifurcationTest);
    fprintf(fp, "vessel_function = %d.\n", (*data).vesselFunction);
    fprintf(fp, "reset_d_lim = %d.\n", (int) (*data).resetsDLim);    
}

void log_domain(FILE *fp, AbstractDomain *domain, long long int n_term, AbstractConstraintFunction<double, int> *gam,
    AbstractConstraintFunction<double, int> *eps_lim, AbstractConstraintFunction<double, int> *nu)
{
    fprintf(fp, "n_term = %lld.\n", n_term);
    fprintf(fp, "gamma = %f.\n", (*gam).getValue(0));
    fprintf(fp, "eps_lim = %0f.\n", (*eps_lim).getValue(0));
    fprintf(fp, "nu = %f.\n", (*nu).getValue(0));
    fprintf(fp, "n_draw = %d.\n", (*domain).getDraw());
    fprintf(fp, "random seed = %d.\n", (*domain).getSeed());
    fprintf(fp, "characteristic_lenght = %f.\n", (*domain).getCharacteristicLength());
    fprintf(fp, "is_convex_domain = %d.\n", (*domain).isIsConvexDomain());
    fprintf(fp, "min_bif_angle = %f.\n", (*domain).getMinBifurcationAngle());
    fprintf(fp, "is_bif_plane_constrained = %d.\n", (int) (*domain).isIsBifPlaneContrained());
    fprintf(fp, "min_plane_angle = %f.\n", (*domain).getMinPlaneAngle());
    log_gen_data(fp, (*domain).getInstanceData());
}





StagedFRROTreeGeneratorLogger::StagedFRROTreeGeneratorLogger(FILE *fp, StagedFRROTreeGenerator *tree_gen)
{
    (*this).file_ = fp;
    (*this).tree_generator_ = tree_gen;
};

StagedFRROTreeGeneratorLogger::~StagedFRROTreeGeneratorLogger()
{

};

void StagedFRROTreeGeneratorLogger::write()
{   
    FILE *fp = (*this).file_;
    StagedFRROTreeGenerator *generator = (*this).tree_generator_;
    SingleVesselCCOOTree *tree = (SingleVesselCCOOTree *) (*generator).getTree();
    point x0 = (*tree).getXProx();
    double r0 = (*tree).getRootRadius();
    double q0 = (*tree).getQProx();
    StagedDomain* staged_domain = (*(*this).tree_generator_).getDomain();
    vector<AbstractDomain *> *domains = (*staged_domain).getDomains();
    vector<long long int> *n_terms = (*staged_domain).getNTerminals();
    vector<AbstractConstraintFunction<double, int> *> *gams = (*generator).getGams();
    vector<AbstractConstraintFunction<double, int> *> *epsLims = (*generator).getEpsLims();
    vector<AbstractConstraintFunction<double, int> *> *nus = (*generator).getNus();
    int size = (*domains).size();
    fprintf(fp, "Root position = (%f, %f, %f).\n", x0.p[0], x0.p[1], x0.p[2]);
    fprintf(fp, "Root radius = %f.\n", r0);
    fprintf(fp, "Root influx = %f.\n", q0);
    for (int i = 0; i < size; ++i) {
        fprintf(fp, "\n");
        fprintf(fp, "---Stage[%d]---\n", i);
        log_domain(fp, (*domains)[i], (*n_terms)[i], (*gams)[0], (*epsLims)[0], (*nus)[0]);
    }
}