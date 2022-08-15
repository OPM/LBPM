#include "models/MultiPhysController.h"

ScaLBL_Multiphys_Controller::ScaLBL_Multiphys_Controller(
    int RANK, int NP, const Utilities::MPI &COMM)
    : rank(RANK), nprocs(NP), Restart(0), timestepMax(0), num_iter_Stokes(0),
      num_iter_Ion(0), analysis_interval(0), visualization_interval(0),
      tolerance(0), time_conv_max(0), time_conv_MainLoop(0), comm(COMM) {}
ScaLBL_Multiphys_Controller::~ScaLBL_Multiphys_Controller() {}

void ScaLBL_Multiphys_Controller::ReadParams(string filename) {

    // read the input database
    db = std::make_shared<Database>(filename);
    study_db = db->getDatabase("MultiphysController");

    // Default parameters
    timestepMax = 10000;
    Restart = false;
    restart_interval = 100000;
    num_iter_Stokes = 1;
    num_iter_Ion.push_back(1);
    analysis_interval = 500;
    visualization_interval = 10000;
    tolerance = 1.0e-6;
    time_conv_max = 0.0;
    time_conv_MainLoop = 0.0;

    // load input parameters
    if (study_db->keyExists("timestepMax")) {
        timestepMax = study_db->getScalar<int>("timestepMax");
    }
    if (study_db->keyExists("analysis_interval")) {
        analysis_interval = study_db->getScalar<int>("analysis_interval");
    }
    if (study_db->keyExists("visualization_interval")) {
        visualization_interval =
            study_db->getScalar<int>("visualization_interval");
    }
    if (study_db->keyExists("restart_interval")) {
    	restart_interval = study_db->getScalar<int>("restart_interval");
    }
    if (study_db->keyExists("tolerance")) {
        tolerance = study_db->getScalar<double>("tolerance");
    }
    //if (study_db->keyExists( "time_conv" )){
    //	time_conv = study_db->getScalar<double>( "time_conv" );
    //}
    //if (study_db->keyExists( "Schmidt_Number" )){
    //	SchmidtNum = study_db->getScalar<double>( "Schmidt_Number" );
    //}

    // recalculate relevant parameters
    //if (SchmidtNum>1){
    //    num_iter_Stokes = int(round(SchmidtNum/2)*2);
    //    num_iter_Ion = 1;
    //}
    //else if (SchmidtNum>0 && SchmidtNum<1){
    //    num_iter_Ion = int(round((1.0/SchmidtNum)/2)*2);
    //    num_iter_Stokes = 1;
    //}
    //else if (SchmidtNum==1){
    //    num_iter_Stokes = 1;
    //    num_iter_Ion = 1;
    //}
    //else{
    //	ERROR("Error: SchmidtNum (Schmidt number) must be a positive number! \n");
    //}

    // load input parameters
    // in case user wants to have an absolute control over the iternal iteration
    if (study_db->keyExists("num_iter_Ion_List")) {
        num_iter_Ion.clear();
        num_iter_Ion = study_db->getVector<int>("num_iter_Ion_List");
    }
    if (study_db->keyExists("num_iter_Stokes")) {
        num_iter_Stokes = study_db->getScalar<int>("num_iter_Stokes");
    }
}

int ScaLBL_Multiphys_Controller::getStokesNumIter_PNP_coupling(
    double StokesTimeConv, const vector<double> &IonTimeConv) {
    //Return number of internal iterations for the Stokes solver
    int num_iter_stokes;
    vector<double> TimeConv;

    TimeConv.assign(IonTimeConv.begin(), IonTimeConv.end());
    TimeConv.insert(TimeConv.begin(), StokesTimeConv);
    vector<double>::iterator it_max =
        max_element(TimeConv.begin(), TimeConv.end());
    int idx_max = distance(TimeConv.begin(), it_max);
    if (idx_max == 0) {
        num_iter_stokes = 2;
    } else {
        double temp =
            2 * TimeConv[idx_max] /
            StokesTimeConv; //the factor 2 is the number of iterations for the element has max time_conv
        num_iter_stokes = int(round(temp / 2) * 2);
    }
    return num_iter_stokes;
}

vector<int> ScaLBL_Multiphys_Controller::getIonNumIter_PNP_coupling(
    double StokesTimeConv, const vector<double> &IonTimeConv) {
    //Return number of internal iterations for the Ion transport solver
    vector<int> num_iter_ion;
    vector<double> TimeConv;
    TimeConv.assign(IonTimeConv.begin(), IonTimeConv.end());
    TimeConv.insert(TimeConv.begin(), StokesTimeConv);
    vector<double>::iterator it_max =
        max_element(TimeConv.begin(), TimeConv.end());
    unsigned int idx_max = distance(TimeConv.begin(), it_max);
    if (idx_max == 0) {
        for (unsigned int idx = 1; idx < TimeConv.size(); idx++) {
            double temp =
                2 * StokesTimeConv /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
    } else if (idx_max == 1) {
        num_iter_ion.push_back(2);
        for (unsigned int idx = 2; idx < TimeConv.size(); idx++) {
            double temp =
                2 * TimeConv[idx_max] /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
    } else if (idx_max == TimeConv.size() - 1) {
        for (unsigned int idx = 1; idx < TimeConv.size() - 1; idx++) {
            double temp =
                2 * TimeConv[idx_max] /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
        num_iter_ion.push_back(2);
    } else {
        for (unsigned int idx = 1; idx < idx_max; idx++) {
            double temp =
                2 * TimeConv[idx_max] /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
        num_iter_ion.push_back(2);
        for (unsigned int idx = idx_max + 1; idx < TimeConv.size(); idx++) {
            double temp =
                2 * TimeConv[idx_max] /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
    }
    return num_iter_ion;
}

vector<int> ScaLBL_Multiphys_Controller::getIonNumIter_NernstPlanck_coupling(
    const vector<double> &IonTimeConv) {
    //Return number of internal iterations for the Ion transport solver
    vector<double> TimeConv;
    TimeConv.assign(IonTimeConv.begin(), IonTimeConv.end());
    vector<int> num_iter_ion;
    vector<double>::iterator it_max = max_element(TimeConv.begin(), TimeConv.end());
    unsigned int idx_max = distance(TimeConv.begin(), it_max);
    if (idx_max == 0) {
        num_iter_ion.push_back(2);
        for (unsigned int idx = 1; idx < TimeConv.size(); idx++) {
            double temp =
                2 * TimeConv[idx_max] /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
    } else if (idx_max == TimeConv.size() - 1) {
        for (unsigned int idx = 0; idx < TimeConv.size() - 1; idx++) {
            double temp =
                2 * TimeConv[idx_max] /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
        num_iter_ion.push_back(2);
    } else {
        for (unsigned int idx = 0; idx < idx_max; idx++) {
            double temp =
                2 * TimeConv[idx_max] /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
        num_iter_ion.push_back(2);
        for (unsigned int idx = idx_max + 1; idx < TimeConv.size(); idx++) {
            double temp =
                2 * TimeConv[idx_max] /
                TimeConv
                    [idx]; //the factor 2 is the number of iterations for the element has max time_conv
            num_iter_ion.push_back(int(round(temp / 2) * 2));
        }
    }
    return num_iter_ion;
}


void ScaLBL_Multiphys_Controller::getTimeConvMax_PNP_coupling(
    double StokesTimeConv, const vector<double> &IonTimeConv) {
    //Return maximum of the time converting factor from Stokes and ion solvers
    vector<double> TimeConv;

    TimeConv.assign(IonTimeConv.begin(), IonTimeConv.end());
    TimeConv.insert(TimeConv.begin(), StokesTimeConv);
    time_conv_max = *max_element(TimeConv.begin(), TimeConv.end());
}
