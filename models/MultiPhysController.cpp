#include "models/MultiPhysController.h"

ScaLBL_Multiphys_Controller::ScaLBL_Multiphys_Controller(int RANK, int NP, MPI_Comm COMM):
rank(RANK),Restart(0),timestepMax(0),num_iter_Stokes(0),num_iter_Ion(0),SchmidtNum(0),comm(COMM)
{

}
ScaLBL_Multiphys_Controller::~ScaLBL_Multiphys_Controller(){

}

void ScaLBL_Multiphys_Controller::ReadParams(string filename){
    
    // read the input database 
	db = std::make_shared<Database>( filename );
	study_db = db->getDatabase( "MultiphysController" );
    

    // Default parameters
    timestepMax = 10000;
    Restart = false;
    SchmidtNum = 1.0;
    num_iter_Stokes=1;
    num_iter_Ion=1;
	
    // load input parameters
	if (study_db->keyExists( "timestepMax" )){
		timestepMax = study_db->getScalar<int>( "timestepMax" );
	}
	if (study_db->keyExists( "Schmidt_Number" )){
		SchmidtNum = study_db->getScalar<double>( "Schmidt_Number" );
	}
    // recalculate relevant parameters
    if (SchmidtNum>1){
        num_iter_Stokes = int(round(SchmidtNum/2)*2);
        num_iter_Ion = 1;
    }
    else if (SchmidtNum>0 && SchmidtNum<1){
        num_iter_Ion = int(round((1.0/SchmidtNum)/2)*2);
        num_iter_Stokes = 1;
    }
    else{
		ERROR("Error: SchmidtNum (Schmidt number) must be a positive number! \n");
    }

    // load input parameters
    // in case user wants to have an absolute control over the iternal iteration
    if (study_db->keyExists( "num_iter_Ion" )){
        num_iter_Ion = study_db->getScalar<int>( "num_iter_Ion" );
    }
    if (study_db->keyExists( "num_iter_Stokes" )){
        num_iter_Stokes = study_db->getScalar<int>( "num_iter_Stokes" );
    }

}


