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
#include"../structures/domain/SimpleDomain2D.h"
#include"../structures/domain/SimpleDomain.h"
#include"../structures/domain/IntersectionVascularizedDomain.h"
#include"../structures/domain/DomainNVR.h"
#include"../structures/domain/PartiallyVascularizedDomain.h"
#include"../structures/domain/DummyDomain.h"
#include"../structures/domain/StagedDomain.h"
#include"../constrains/AbstractConstraintFunction.h"
#include"../structures/tree/AbstractObjectCCOTree.h"
#include"../structures/tree/SingleVesselCCOOTree.h"
#include"../constrains/AbstractConstraintFunction.h"

// This class interface
#include"StagedFRROTreeGeneratorLogger.h"

void log_domain_files(FILE *fp, AbstractDomain *domain) 
{   
    int whichDomain = domain->getWhichDomain();
    if(whichDomain == 0) {
        SimpleDomain2D* typed_domain = (SimpleDomain2D *) domain;
        fprintf(fp, "SimpleDomain2D\n");
        fprintf(fp, "filename = %s\n", typed_domain->getFilename().c_str());
    }
    else if (whichDomain == 1) {
        SimpleDomain* typed_domain = (SimpleDomain *) domain;
        fprintf(fp, "SimpleDomain\n");
        fprintf(fp, "filename = %s\n", typed_domain->getFilename().c_str());
    }
    else if (whichDomain == 2) {
        DomainNVR* typed_domain = (DomainNVR *) domain;
        fprintf(fp, "DomainNVR\n");
        fprintf(fp, "filenameHull = %s\n", typed_domain->getFilenameHull().c_str());
        vector<string> filenameNVR = typed_domain->getFilenameNVR();
        int size = filenameNVR.size();
        for (int i = 0; i < size; ++i) {
            fprintf(fp, "filenameNVR[%d] = %s\n", i, filenameNVR[i].c_str());
        }
    }
    else if (whichDomain == 3) {
        IntersectionVascularizedDomain* typed_domain = (IntersectionVascularizedDomain *) domain;
        fprintf(fp, "IntersectionVascularizedDomain\n");
        vector<string> filenameVR = typed_domain->getFilenameVR();
        int size = filenameVR.size();
        for (int i = 0; i < size; ++i) {
            fprintf(fp, "filenameVR[%d] = %s\n", i, filenameVR[i].c_str());
        }
    }
    else if (whichDomain == 4) {
        PartiallyVascularizedDomain* typed_domain = (PartiallyVascularizedDomain *) domain;
        string filenameHull = typed_domain->getFilenameHull();
        vector<string> filenameVR = typed_domain->getFilenameVR();
        int size_vr = filenameVR.size();
        vector<string> filenameNVR = typed_domain->getFilenameNVR();
        int size_nvr = filenameNVR.size();
        fprintf(fp, "PartiallyVascularizedDomain\n");
        fprintf(fp, "filenameHull = %s\n", filenameHull.c_str());
        for (int i = 0; i < size_vr; ++i) {
            fprintf(fp, "filenameVR[%d] = %s\n", i, filenameVR[i].c_str());
        }
        for (int i = 0; i < size_nvr; ++i) {
            fprintf(fp, "filenameNVR[%d] = %s\n", i, filenameNVR[i].c_str());
        }
    }
    else if (whichDomain == 5) {
        DummyDomain* typed_domain = (DummyDomain *) domain;
        fprintf(fp, "DummyDomain\n");
        fprintf(fp, "filename = %s\n", typed_domain->getFilename().c_str());
    }
    else if (whichDomain == 6) {
        // StagedDomain* typed_domain = (StagedDomain *) domain;
        fprintf(fp, "StagedDomain\n");
    }
    else {
        fprintf(stderr, "This is not a valid domain!\n");
        exit(EXIT_FAILURE);
    }
}

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
    double q0 = (*tree).getQProx();
    StagedDomain* staged_domain = (*(*this).tree_generator_).getDomain();
    vector<AbstractDomain *> *domains = (*staged_domain).getDomains();
    vector<long long int> *n_terms = (*staged_domain).getNTerminals();
    vector<AbstractConstraintFunction<double, int> *> *gams = (*generator).getGams();
    vector<AbstractConstraintFunction<double, int> *> *epsLims = (*generator).getEpsLims();
    vector<AbstractConstraintFunction<double, int> *> *nus = (*generator).getNus();
    int size = (*domains).size();
    string filenameCCO = tree->getFilenameCCO();
    if(filenameCCO.empty()) {
        point x0 = (*tree).getXProx();
        double r0 = (*tree).getRootRadius();
        fprintf(fp, "Root position = (%f, %f, %f).\n", x0.p[0], x0.p[1], x0.p[2]);
        fprintf(fp, "Root radius = %f.\n", r0);
        fprintf(fp, "Root influx = %f.\n", q0);
    }
    else {
        fprintf(fp, "Input CCO filename = %s\n", filenameCCO.c_str());
        fprintf(fp, "Root influx = %f.\n", q0);
    }
    
    for (int i = 0; i < size; ++i) {
        fprintf(fp, "\n");
        fprintf(fp, "Stage[%d]\n", i);
        log_domain_files(fp, (*domains)[i]);
        log_domain(fp, (*domains)[i], (*n_terms)[i], (*gams)[0], (*epsLims)[0], (*nus)[0]);
    }
    
    fprintf(fp, "\n");
    fprintf(fp, "Initial dLim = %f.\n", generator->getDLimInitial());
    fprintf(fp, "Last dLim = %f.\n", generator->getDLimLast());
    time_t begin_time = generator->getBeginTime();
    time_t end_time = generator->getEndTime();
    struct tm *initial_tm = localtime(&begin_time);
    struct tm *last_tm = localtime(&end_time);
    char time_initial_c_string[21];
    char time_last_c_string[21];
    strftime(time_initial_c_string, 20, "%d_%m_%Y_%H_%M_%S", initial_tm);
    strftime(time_last_c_string, 20, "%d_%m_%Y_%H_%M_%S", last_tm);
    fprintf(fp, "\n");
    fprintf(fp, "Beginning of generation time = %s\n", time_initial_c_string);
    fprintf(fp, "End of generation time = %s\n", time_last_c_string);
}