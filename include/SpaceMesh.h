#ifndef SPACE_MESH__H
#define SPACE_MESH__H

#include "global.h"

typedef enum {fems_none = -1, fems_abaqus, fems_LSDYNA, FEMSolverT_SIZE} FEMSolverT;
typedef enum {ELEMENT_SHELL=1, NODE_ID, NODE_COORD,NODE_VEL, SET_SHELL_LIST_TITLE, DATABASE_HISTORY_NODE_ID
, LOAD_NODE_POINT, BOUNDARY_PRESCRIBED_MOTION_NODE, SET_NODE_LIST, VEL_ID} InputPart;

// Disp: Dirichlet: full strain tensor is used -> u = eps x
// Trac: Neumann  : full stress tensor is used -> trac = sigma n
// Mixed: on different facets and components of facets uses either Dirichlet or Neumann
//			Mine:	epsX: epsilon00, sigma01, sigma11 are given	->
//						on bottom and top, only tractions consistent with sigma01, sigma11 are applied
//						on left, right shear stress is sigma01, normal displacement is epsilon00 * x0

//					epsY: epsilon11, sigma01, sigma00

//					mbc_MixedMineDispNormal:
//												case0 -> epsX with epsilon00 != 0 
//												case1 -> epsY with epsilon11 != 0
//												case2 -> epsX with epsilon00 = sigma11 = 0, sigma01 != 0
//					mbc_MixedMineTracNormal:
//												case0 -> epsY with sigma00 != 0 
//												case1 -> epsX with sigma11 != 0
//												case2 -> epsX with epsilon00 = sigma11 = 0, sigma01 != 0

//		Hazanov95: On overall properties of elastic heterogeneous bodies smaller than the representative volume.pdf
//				   section 3.
//					epsX: epsilon00, sigma01, sigma11 are given	->
//						on bottom and top, sigma11 is applied 
//						BUT u1 = eps00 x0
//						on left, right shear stress is sigma01, normal displacement is epsilon00 * x0
//					Breakdown of displacement and traction cases is similar

//		mbc_DispNoShear == Disp but does not apply "shear" displacements. That is on left and right edges ux is specified and on top and bottom uy
//		Load case 0: strainxx = 1, case 1: strainyy = 1, case 0: strainxx = strainyy = 1
typedef enum {mbc_Disp, mbc_Trac, mbc_MixedMineDispNormal, mbc_MixedMineTracNormal, mbc_MixedHazanov95DispNormal, mbc_MixedHazanov95TracNormal, mbc_DispNoShear, HomogBCType_SIZE} HomogBCType;

//typedef enum {sm_bottom, sm_right, sm_top, sm_left, sm_numSide} sideT;
typedef enum {sm_bl, sm_br, sm_tr, sm_tl, cornerT_SIZE} cornerT;

typedef enum {bt_none = -1, bt_neumann, bt_dirichlet, bt_dirichlet0, bt_SIZE} boundaryT;
void UpdateBCType(boundaryT& bt2Update, boundaryT newbt);

class sMesh;

class BCInfo
{
public:
	BCInfo();

	unsigned int caseNo;
	// indexed as [xysideIndex][minmaxIndex][Uomponent index]
	// for example side_BCType[1][mx_t][0] is for y (1) max value (mx_t) -> top face, component 0
	boundaryT side_BCType[DIM][mnmxT_SIZE][DIM_U];
	// lists corners and corresponding directions that need to be fixed, second is it's corresponding direction
	vector< pair<cornerT, int> > corner2Dir_4RigidControl;

	vector<double> voigtStrain;
	vector<double> voigtStress;

	void set_BCType(boundaryT bt = bt_none);
};

class BCConfig
{
public:
	BCConfig();
	void Initialize(string inputFolder4RunIn);
	void set_problem_type(bool& type);

	string inputFolder4Run;
	string img_dir;
	bool is_static;
	// input parameters
	HomogBCType	homogType;
	// voigt strain and stress are computed from the following two
	double strainScale;
	double stiffnessScale; // like elastic modulus

	// computed
	double stressScale;
	BCInfo	bcInfoCases[VOIGHT_STN_SZ];

	void FormCases4AllDispBC(unsigned int caseNo, BCInfo& bcInfoCase);
	void FormCases4AllTracBC(unsigned int caseNo, BCInfo& bcInfoCase);
	void FormCasesMixedMineHanzov95BC(HomogBCType homogT, unsigned int caseNo, BCInfo& bcInfoCase);
};


class smNode
{
public:
	long id;
	Vec crd;
};

class elementCP4
{
public:
	elementCP4();
	long e4_id;
	long vs4[4];
	long matID;
	double area;
};

class elset
{
public:
	string name;
	vector<long> el_ids;
};

class sSide
{
public:
	void setSize(unsigned long szIn);
	void setSegmentLengthsAndWeights();
	unsigned posCrdChanging;
	unsigned long sz;
//	vector<int> vertexIndices;
	vector<int> vertexIDs;
	vector<double> crdRel2Centroid[DIM];
	vector<double> weight;
	vector<double> segmentLengths;
	Vec normal;
};

class sPoint
{
public:
	sPoint();
//	int vertexIndex;
	int vertexID;
	Vec crdRel2Centroid;
	Vec normalTimesWeights;
	mnmxT pos[DIM];
};

// point displacement force class for 1 direction
class boundary_1Point_uF1Dir
{
public:
	boundary_1Point_uF1Dir();
	double u; // displacement
	double F; // force
	boundaryT bt;
};

class boundary_1Point_uFAllDirs
{
public:
	sPoint point;
	boundary_1Point_uF1Dir	op_ufs[DIM_U];
};

class boundary_AllPoints_uFAllDirs_1Case
{
public:
	void PrintBC(BCConfig& bcConfig, FEMSolverT solverOption, ostream& out, string inputFolder4Run, sMesh* mesh);
	vector<boundary_1Point_uFAllDirs> allDomain_uFs;
	map<int, int> vertexID2Pos;
	int sz_ptBCs;
	double maxF;
};

class OneLoadCase_MeanStnStrs
{
public:
	vector<int> timeStepNos;
	vector<double> timeValues;
	vector< vector<double> > voigtStrains;
	vector< vector<double> > voigtStresses;
};

class AllLoadCases_MeanStnStrs
{
public:
	void AllLoadCases_MeanStnStrs_Write(ostream& out);
	OneLoadCase_MeanStnStrs lc_stn_sts[VOIGHT_STN_SZ];
	vector<MATRIX> stiffnesses;
	int num_timeSteps;
};

class sMesh
{
public:
//	string fileNameWOExt;
//	vector<string> fileNamesOut;

	void Read(const string& fileNameWOExtIn, FEMSolverT solverOptionOut_In, FEMSolverT solverOptionIn_In);
	void ComputePrintBC_AllCases(BCConfig& bcConfig, vector<boundary_AllPoints_uFAllDirs_1Case>& bcs);

//	vector<smNode> nodes;
	map<long, smNode> nodes;
	int n_nodes;
	vector<elementCP4> e4;
	map<long, long> id2_e4Pos;
	int n_e4;
	double area_e4;
	vector<elset> elss;
	int num_elsets;
	Vec minCrd;
	Vec maxCrd;
	Vec meanCrd;
	Vec spanCrd;
	double volume;

	sSide sides[DIM][mnmxT_SIZE];
	// first index is on x, second is on y
	// [mn][mn] lower left
	// [mn][mx] xmin, ymax -> upper right
	sPoint corners[mnmxT_SIZE][mnmxT_SIZE];
	int boundaryPoints_cornerPoss[cornerT_SIZE];
	vector<sPoint> boundaryPoints;
	unsigned long sz_boundaryPoints;

	void Compute_AllStrain_BCs(const vector<double>& voigtStrain, boundary_AllPoints_uFAllDirs_1Case& bcs);
	void Compute_AllTraction_BCs(const vector<double>& voigtStress, boundary_AllPoints_uFAllDirs_1Case& bcs, unsigned int caseNo);

	void FormEmpty_BoundaryStorage(boundary_AllPoints_uFAllDirs_1Case& bcs);
	void Compute_GeneralMixed_BCs(BCInfo& bcInfo, boundary_AllPoints_uFAllDirs_1Case& bcs);
	void ComputePrintBC_AllCases_Aux(BCConfig& bcConfig, const vector<string>& fileNames, FEMSolverT solverOptionIn, vector<boundary_AllPoints_uFAllDirs_1Case>& bcs);
	void ComputePrintBC(BCConfig& bcConfig, unsigned int caseNo, FEMSolverT solverOptionIn, ostream& out, boundary_AllPoints_uFAllDirs_1Case& bcs, string inputFolder4Run);

	// bcs: input
	// voigt notation for 2D/3D elasticity, for anti-plane shear it's normal strain / stress vectors
	void Compute_Average_Strain_Stress_From_BoundaryNode_Disps_and_Forces(const boundary_AllPoints_uFAllDirs_1Case& bcs, vector<double>& voigtStrain, vector<double>& voigtStress);
	// onlyIncludeLastStep == true for static problems and dynamics ones [when we only want the last solution]
	int Compute_Average_Strain_Stress_From_BoundaryNode_Disps_and_Forces(const vector<boundary_AllPoints_uFAllDirs_1Case>& bcss, vector< vector<double> >& voigtStrains, vector< vector<double> >& voigtStresses, bool onlyIncludeLastStep);

	/// second stage
	bool Homogenize_CFromAllLoadCases(const string& fileNameWOExtIn, FEMSolverT solverOptionOut, bool onlyIncludeLastStep, AllLoadCases_MeanStnStrs& homogResults);
private:
	string fileNameWOExt;
	string extIn, extOut;
	FEMSolverT solverOptionIn;
	FEMSolverT solverOptionOut;
	vector<string> fileNamesOut;
	void Read(istream& in);
	void Initialize();
	void PrintNode_Element_ElementSet_LSDYNA(ostream& out);
	// return is the size of time steps

private:
	int Read_LS_DYNA_Boundary_Solutions_1LoadCase(const string& fileName_u, const string& fileName_F, 
		OneLoadCase_MeanStnStrs& lc_stn_sts, bool onlyIncludeLastStep, bool& successfulread);
	double CalculateQuadArea(const vector<smNode*>& nodePtrs);
};

// factor multiplies ten2ndOrder to get vec (e.g. 1/Volume)
double CalculateQuadAreaAux(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);
void FullElasticity_StnStnTensor2VoigtVector(const MATRIX& ten2ndOrder, vector<double>& vec, bool isStrain, double factor = 1.0);
bool DoublesAreEqual(double d1, double d2, double tol);

unsigned int BreakString(const string& inString, vector<string>& parts);
void CopyContentFileA_ostream(const string& fileA, ostream& out);

void getNumberSpaces(string& sp, const int& number, InputPart& ip);
// returns true if read is successful
bool Read_LSDYNA_NodalDisplacements(istream& in, unsigned int& timeStepNo, double& timeValue);
bool Read_LSDYNA_NodalForces(istream& in, double& timeValue);

bool checkZeroString( std::string& rstr);
string getScientificNumber( const std::string& rstr, const double& val);
#endif
