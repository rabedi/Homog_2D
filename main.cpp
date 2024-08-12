#include "SpaceMesh.h"
#include "FileIOMacros.h"
#include <experimental/filesystem>
#include <sys/stat.h>
#if !VCPP
#include <unistd.h>
#endif
#include <stdio.h>
using namespace std;
string g_rootFolder_4Output = "";
string main_dir="";
string inc_dir = "";
    

int main(int argv, char* argc[])
{
	string option_fem_input, option_bc, option_static, work_dir_base;
    bool generateFEMInputs;
	bool testdata = true;
	if(argv==1)
	{
		main_dir="/media/ghuynh/Data1/mechanics/code/project/SPDE/";
		work_dir_base = main_dir + "homo_data_disp/Res8x8/";
		g_rootFolder_4Output = work_dir_base + "Img12/SubImg34/";
		if (testdata)
		{
			main_dir = "Test/";
			work_dir_base = main_dir;
			g_rootFolder_4Output = work_dir_base;
		}
		main_dir = "Test/";
		option_fem_input = "true";
		option_bc = "mbc_DispNoShear"; //"mbc_Disp";
		option_static = "true";
	}
	else if (argv< 6)
	{
		cout<< "USAGE: " << endl;
        cout<< "option1: ./main output_directory   generateFEMInputs  bcConfig.homogType  is_static  project_directory\n";
		cout<< "generateFEMInputs=true -> generate input mesh; false -> homogenize material properties\n";
		cout<< "homogType= mbc_Disp, mbc_Trac, mbc_MixedMineDispNormal\n";
        return 0;
	}
	else
	{
		// get output
		g_rootFolder_4Output= argc[1];
		main_dir = argc[5];
		option_fem_input = argc[2];
		option_bc=argc[3];
		option_static = argc[4];
	}

	// get main directory
	//g_rootFolder_4Output= argc[1];

	//main_dir = argc[5];

	if (g_rootFolder_4Output != "")
	{

		struct stat sb;
		if (stat (g_rootFolder_4Output.c_str(), &sb) == 0)
			std::cout<< "the path exists\n";
		else
			MakeDir(g_rootFolder_4Output);;
		//g_rootFolder_4Output += "/";
	}

	// For Giang: change this boolean to true to generate LS-DYNA input, to false to read these results and homogenize stuff
	// if false -> generates input files for FEA (Abaqus/LS-DYNA), true -> reads their results to compute average stresses / srains -> compute Cs (stiffnesses0
	//option_fem_input = argc[2];
	
	if(!option_fem_input.compare("true"))
		generateFEMInputs = true; 
	else
		generateFEMInputs = false;

	string fileNameWOExtIn = "Micro";
	sMesh mesh;
	FEMSolverT solverOptionIn = fems_abaqus;
	FEMSolverT solverOptionOut = fems_abaqus; // fems_LSDYNA; // fems_LSDYNA;


	mesh.Read(fileNameWOExtIn, solverOptionOut, solverOptionIn);


	if (!generateFEMInputs)
	{
		bool onlyIncludeLastStep = true; // true for static/implicit/linear problems | false for dynamic/explicit/fracture problems
		AllLoadCases_MeanStnStrs homogResults;
		mesh.Homogenize_CFromAllLoadCases(fileNameWOExtIn, solverOptionOut, onlyIncludeLastStep, homogResults);
		cout << "successful\n";
		//getchar();
		return 0;
	}
	BCConfig bcConfig;
	bcConfig.stiffnessScale = 10.0; //1000.0;
	bcConfig.strainScale = 1.0;

	// mbc_Disp, mbc_Trac, mbc_MixedMineDispNormal, mbc_MixedMineTracNormal, mbc_MixedHazanov95DispNormal, mbc_MixedHazanov95TracNormal
	//option_bc=argc[3];
	
	if(!option_bc.compare("mbc_Disp"))
	{
		bcConfig.homogType = mbc_Disp;
		inc_dir = main_dir+"homo_data_disp/";
	}
	else if (!option_bc.compare("mbc_DispNoShear"))
	{
		bcConfig.homogType = mbc_DispNoShear;
		inc_dir = main_dir + "homo_data_dispNoShear/";
	}
	else if(!option_bc.compare("mbc_Trac"))
	{
		bcConfig.homogType =  mbc_Trac; // mbc_Disp;
		inc_dir = main_dir+"homo_data_trac/";

	}
	else
	{
		bcConfig.homogType =  mbc_MixedMineDispNormal; 
		inc_dir = main_dir+"homo_data_mixed/";

	}
	//bcConfig.img_dir += "global_params.k";

	std::cout<< "boundary type: " << bcConfig.homogType << std::endl;


	// where certain header, ... excerpts for this run exist
	bool is_static;
	//option_static = argc[4];
	if(!option_static.compare("true") )
		is_static = true;
	else
		is_static=false;
	string inputFolder4RunIn =main_dir+"generate_lsdyna_input/LS_DYNA_IMPLICIT/";


	
	bcConfig.Initialize(inputFolder4RunIn);
	bcConfig.set_problem_type(is_static);

	vector<boundary_AllPoints_uFAllDirs_1Case> bcs;
	mesh.ComputePrintBC_AllCases(bcConfig, bcs);
	cout << "run complete\n";
	//getchar();
	return 0;
}

