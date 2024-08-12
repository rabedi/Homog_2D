#include "SpaceMesh.h"
#include "FileIOMacros.h"
#include "commonMacros.h"
#include "bits/stdc++.h"

#define TST_HOMOG 0

using namespace std;
BCConfig::BCConfig()
{
	homogType = mbc_Disp;
	strainScale = 1.0;
	stiffnessScale = 1.0;
	stressScale = 1.0;
}

void BCConfig::Initialize(string inputFolder4RunIn)
{
	inputFolder4Run = inputFolder4RunIn;
	stressScale = stiffnessScale * strainScale;
	for (unsigned caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
		bcInfoCases[caseNo].caseNo = caseNo;

	if ((homogType == mbc_Disp) || (homogType == mbc_DispNoShear))
	{
		for (unsigned caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
			FormCases4AllDispBC(caseNo, bcInfoCases[caseNo]);
	}
	else if (homogType == mbc_Trac)
	{
		for (unsigned caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
			FormCases4AllTracBC(caseNo, bcInfoCases[caseNo]);
	}
	else if ((homogType == mbc_MixedMineDispNormal) || (homogType == mbc_MixedHazanov95DispNormal)
		|| (homogType == mbc_MixedMineTracNormal) || (homogType == mbc_MixedHazanov95TracNormal))
	{
		for (unsigned caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
			FormCasesMixedMineHanzov95BC(homogType, caseNo, bcInfoCases[caseNo]);
	}
	else
	{
		cout << "homogType\t" << homogType << '\n';
		THROW("Invalid homogType\n");
	}
}

void BCConfig::set_problem_type(bool& type)
{
	is_static = type;
}

void BCConfig::FormCases4AllDispBC(unsigned int caseNo, BCInfo& bcInfoCase)
{
	bcInfoCase.set_BCType(bt_dirichlet);
	bcInfoCase.voigtStrain[caseNo] = strainScale;
	if (homogType == mbc_DispNoShear)
	{
		if (caseNo == 2)
		{
			bcInfoCase.voigtStrain[0] = strainScale;
			bcInfoCase.voigtStrain[1] = strainScale;
			bcInfoCase.voigtStrain[2] = 0.0;
		}
		for (unsigned int epsDir = 0; epsDir < DIM; ++epsDir)
		{
			unsigned int tangentDir = 1 - epsDir;
			bcInfoCase.side_BCType[epsDir][mn_t][tangentDir] = bt_neumann;
			bcInfoCase.side_BCType[epsDir][mx_t][tangentDir] = bt_neumann;
		}
		for (unsigned int crnr = 0; crnr < cornerT_SIZE; ++crnr)
		{
			bcInfoCase.corner2Dir_4RigidControl.push_back(pair<cornerT, int>((cornerT)crnr, 0));
			bcInfoCase.corner2Dir_4RigidControl.push_back(pair<cornerT, int>((cornerT)crnr, 1));
		}
	}
}

void BCConfig::FormCases4AllTracBC(unsigned int caseNo, BCInfo& bcInfoCase)
{
	bcInfoCase.set_BCType(bt_neumann);
	bcInfoCase.voigtStress[caseNo] = stressScale;
	bcInfoCase.corner2Dir_4RigidControl.push_back(pair<cornerT, int>(sm_bl, 0));
#if !ANTI_PLANE_SHEAR
	bcInfoCase.corner2Dir_4RigidControl.push_back(pair<cornerT, int>(sm_bl, 1));
#endif
	if (caseNo != 1)
		bcInfoCase.corner2Dir_4RigidControl.push_back(pair<cornerT, int>(sm_br, 1));
	else
		bcInfoCase.corner2Dir_4RigidControl.push_back(pair<cornerT, int>(sm_tl, 0));
}

elementCP4::elementCP4()
{
	matID = -1;
	e4_id = -1;
}	

void BCConfig::FormCasesMixedMineHanzov95BC(HomogBCType homogT, unsigned int caseNo, BCInfo& bcInfoCase)
{
	// the direction at which eps is specified
	int epsDir = -1;
	if (caseNo == 2)
	{
		bcInfoCase.voigtStress[2] = stressScale;
		epsDir = 0;
	}
	else // normal modes
	{
		if ((homogT == mbc_MixedMineDispNormal) || (homogT == mbc_MixedHazanov95DispNormal))
		{
			bcInfoCase.voigtStrain[caseNo] = strainScale;
			epsDir = caseNo;
 		}
		else if ((homogT == mbc_MixedMineTracNormal) || (homogT == mbc_MixedHazanov95TracNormal))
		{
			bcInfoCase.voigtStress[caseNo] = stressScale;
			epsDir = 1 - caseNo;
		}
		else
		{
			cout << "homogT\t" << homogT << '\n';
			THROW("Invalid homogT\n");
		}
	}
	bool mine = ((homogT == mbc_MixedMineDispNormal) || (homogT == mbc_MixedMineTracNormal));
	int other_epsDir = 1 - epsDir;
	int uDir0 = epsDir;
	int uDir1 = 1 - uDir0;
#if ANTI_PLANE_SHEAR
	uDir0 = 0;
	uDir1 = 0;
#endif
	// normal components
	bcInfoCase.side_BCType[epsDir][mn_t][uDir0] = bt_dirichlet;
	bcInfoCase.side_BCType[epsDir][mx_t][uDir0] = bt_dirichlet;
	bcInfoCase.side_BCType[other_epsDir][mn_t][uDir1] = bt_neumann;
	bcInfoCase.side_BCType[other_epsDir][mx_t][uDir1] = bt_neumann;
#if !ANTI_PLANE_SHEAR
	bcInfoCase.side_BCType[epsDir][mn_t][uDir1] = bt_neumann;
	bcInfoCase.side_BCType[epsDir][mx_t][uDir1] = bt_neumann;
	boundaryT bt = bt_neumann;
	if (!mine)
		bt = bt_dirichlet;
	bcInfoCase.side_BCType[other_epsDir][mn_t][uDir0] = bt;
	bcInfoCase.side_BCType[other_epsDir][mx_t][uDir0] = bt;

	// controling rigid
	if (mine)
		bcInfoCase.corner2Dir_4RigidControl.push_back(pair<cornerT, int>(sm_bl, other_epsDir));
#endif
}

void UpdateBCType(boundaryT& bt2Update, boundaryT newbt)
{
	if (bt2Update == bt_dirichlet)
		return;
	if (newbt != bt_none)
	{
		bt2Update = newbt;
		return;
	}
}

void getNumberSpaces(string& sp, const int& number, InputPart& ip )
{
	if(ip==BOUNDARY_PRESCRIBED_MOTION_NODE||ip==DATABASE_HISTORY_NODE_ID||ip==SET_SHELL_LIST_TITLE||ip==LOAD_NODE_POINT||
	    ip==SET_NODE_LIST||ip==VEL_ID||ip==ELEMENT_SHELL||ip==NODE_ID)
	{
		if (number< 1.0e1)
			sp = "                   ";
		else if( number<1.0e2 && number>=1.0e1)
			sp = "                  ";
		else if( number <1.0e3 && number >=1.0e2)
			sp = "                 ";
		else if( number <1.0e4 && number >=1.0e3)
			sp = "                ";
		else if( number <1.0e5 && number >=1.0e4)
			sp = "               ";
		else if( number <1.0e6 && number >=1.0e5)
			sp = "              ";
		else if( number <1.0e7 && number >=1.0e6)
			sp = "             ";
		else if( number <1.0e8 && number >=1.0e7)
			sp = "            ";
		else 
		{
			std::cout<< "can't return the number of spaces as the input number is larger than 1.0e6\n";
			throw ;
		}
	}
}

string getScientificNumber( const std::string& rstr, const double& val)
{
	int iter=0;
	int posE;
	string numstr;
	while(rstr[iter])
	{
		if(iter==0)
		  numstr = rstr[0];
		// else if(iter==0&&val<0.0)
		// {
		// 	string vstr ; toString(val, vstr);
		// 	numstr = vstr[0];
		// 	if(rstr[0]!='-')
		// 		numstr = numstr+rstr[0];

		// }


		if(rstr[iter]=='e')
		{
			posE = iter;
			break;
		}
		if(iter<5&&iter>0&&rstr[iter]!='e'&&rstr[iter]!='-')
			numstr = numstr+ rstr[iter];
		iter+=1;
	}

	if(posE==1)
		numstr= numstr+".000";
	if(posE==3)
		numstr= numstr+"00";
	if(posE==4)
		numstr= numstr+"0";
	iter=0;
	while(rstr[iter])
	{
		if(iter>=posE)
			numstr = numstr+rstr[iter];
		iter+=1;
	}

	if (val<0.0)
		numstr = '-'+numstr;
	return numstr;
}

bool checkZeroString( std::string& rstr)
{
	bool IsZero;
	int iter = 0;
	int ZeroCount=0;
	if(rstr[0]!='0')
	{
		IsZero=false;
		return IsZero;
	}
		

	while (rstr[iter])
	{
		if(rstr[iter]=='0')
			ZeroCount+=1;
		iter++;
	}

	if(ZeroCount==9 || (iter==1&ZeroCount==1)) // considered to be zero if 0.0000000 or 0
		IsZero=true;
	else
		IsZero=false; // else different from 0

	return IsZero;
}

BCInfo::BCInfo()
{
	caseNo = 0;
	voigtStrain.resize(VOIGHT_STN_SZ);
	voigtStress.resize(VOIGHT_STN_SZ);
	fill(voigtStrain.begin(), voigtStrain.end(), 0.0);
	fill(voigtStress.begin(), voigtStress.end(), 0.0);
	set_BCType();
}

void BCInfo::set_BCType(boundaryT bt)
{
	for (int sidei = 0; sidei < DIM; ++sidei)
		for (int mnx = 0; mnx < mnmxT_SIZE; ++mnx)
			for (int ucompi = 0; ucompi < DIM_U; ++ucompi)
				side_BCType[sidei][mnx][ucompi] = bt;
}

void sSide::setSize(unsigned long szIn)
{
	sz = szIn;
//	vertexIndices.resize(sz);
	vertexIDs.resize(sz);
	for (int i = 0; i < DIM; ++i)
		crdRel2Centroid[i].resize(sz);
	weight.resize(sz);
	if (sz > 1)
		segmentLengths.resize(sz - 1);
}

void sSide::setSegmentLengthsAndWeights()
{
	if (sz < 2)
		return;
	for (int i = 0; i < sz - 1; ++i)
		segmentLengths[i] = crdRel2Centroid[posCrdChanging][i + 1] - crdRel2Centroid[posCrdChanging][i];
	weight[0] = 0.5 * segmentLengths[0];
	for (int i = 1; i < sz - 1; ++i)
		weight[i] = 0.5 * (segmentLengths[i - 1] + segmentLengths[i]);
	weight[sz - 1] = 0.5 * segmentLengths[sz - 2];

//	double szs = 0.0;
//	for (int i = 0; i < weight.size(); ++i)
//		szs += weight[i];
//	cout << szs << '\n';
}

sPoint::sPoint()
{
	pos[0] = mnx_none_t;
	pos[1] = mnx_none_t;
}

boundary_1Point_uF1Dir::boundary_1Point_uF1Dir()
{
	u = 0.0;
	F = 0.0;
	bt = bt_none;
}

void boundary_AllPoints_uFAllDirs_1Case::PrintBC(BCConfig& bcConfig, FEMSolverT solverOption, ostream & out, string inputFolder4Run, sMesh* mesh)
{
	maxF = 0;
	for (unsigned int i = 0; i < sz_ptBCs; ++i)
	{
		int id = allDomain_uFs[i].point.vertexID;
		if (solverOption == fems_abaqus)
			out << "*Nset, nset=NODE_" << id << '\n' << id << '\n';
		for (int j = 0; j < DIM_U; ++j)
		{
			if (allDomain_uFs[i].op_ufs[j].bt == bt_neumann)
				maxF = MAX(maxF, fabs(allDomain_uFs[i].op_ufs[j].F));
		}
	}
	double tolF = 1e-15 * maxF;


	if (solverOption == fems_abaqus)
	{
		out << "*Nset, nset=NODE_BOUNDARY\n";
		unsigned long cntr = 0;
		unsigned int sz_ptBCsm1 = sz_ptBCs - 1;
		for (unsigned int i = 0; i < sz_ptBCs; ++i)
		{
			int id = allDomain_uFs[i].point.vertexID;
			out << id;
			if ((++cntr % 16 == 0) || (i == sz_ptBCsm1))
				out << '\n';
			else
				out << ",\t";
		}

		// plotting displacements
//		if (solverOption == fems_abaqus)
			out << "*BOUNDARY\n";
		for (unsigned int i = 0; i < sz_ptBCs; ++i)
		{
			for (int j = 0; j < DIM_U; ++j)
			{
				unsigned jp1 = j + 1;
				if (allDomain_uFs[i].op_ufs[j].bt == bt_dirichlet)
				{
					double u = allDomain_uFs[i].op_ufs[j].u;
					int id = allDomain_uFs[i].point.vertexID;
//					if (solverOption == fems_abaqus)
						out << "NODE_" << id << ", " << jp1 << ", " << jp1 << ", " << u << '\n';
				}
			}
		}
		if (maxF == 0)
			return;
		if (solverOption == fems_abaqus)
			out << "*CLOAD\n";
		for (unsigned int i = 0; i < sz_ptBCs; ++i)
		{
			for (int j = 0; j < DIM_U; ++j)
			{
				unsigned jp1 = j + 1;
				if (allDomain_uFs[i].op_ufs[j].bt == bt_neumann)
				{
					double F = allDomain_uFs[i].op_ufs[j].F;
					if (fabs(F) < tolF)
						continue;
					int id = allDomain_uFs[i].point.vertexID;
//					if (solverOption == fems_abaqus)
						out << "NODE_" << id << ", " << jp1 << ", " << F << '\n';
				}
			}
		}
	}
	// abaqus - finished
	else if (solverOption == fems_LSDYNA)
	{
		double maxU = 0;
		for (unsigned int i = 0; i < sz_ptBCs; ++i)
		{
			int id = allDomain_uFs[i].point.vertexID;
			for (int j = 0; j < DIM_U; ++j)
			{
				if (allDomain_uFs[i].op_ufs[j].bt == bt_dirichlet)
					maxU = MAX(maxU, fabs(allDomain_uFs[i].op_ufs[j].u));
			}
		}
		double tolU = 1e-15 * maxU;
		if (tolU < 10.0 * DBL_MIN)
			tolU = tolF / bcConfig.stiffnessScale;

		map<int, vector<int> > vert_pos_2_indices_fixed;
		int id_lb = mesh->corners[mn_t][mn_t].vertexID;
		int id_lt = mesh->corners[mn_t][mx_t].vertexID;
		int id_rb = mesh->corners[mx_t][mn_t].vertexID;

		for (unsigned int i = 0; i < sz_ptBCs; ++i)
		{
			bool add2List = false;
			int id = allDomain_uFs[i].point.vertexID;
			vector<int> poss(3);
			fill(poss.begin(), poss.end(), 0);
			for (int j = 0; j < DIM_U; ++j)
			{
				if ((allDomain_uFs[i].op_ufs[j].bt == bt_dirichlet) && (fabs(allDomain_uFs[i].op_ufs[j].u) < tolU))
				{
					allDomain_uFs[i].op_ufs[j].bt = bt_dirichlet0;
					poss[j] = 1;
					add2List = true;
				}
			}
			if (true) //((id == id_lb) || (id == id_rb) || (id == id_lt))
			{
				poss[2] = 1;
				add2List = true;
			}
			if (add2List)
				vert_pos_2_indices_fixed[i] = poss;
		}
		//4Giang
		out << "*BOUNDARY_SPC_NODE\n";
		out <<"$#               nid                 cid                dofx                dofy                dofz               dofrx               dofry               dofrz\n";
		string delim_id;
		map<int, vector<int> >::iterator it, itb = vert_pos_2_indices_fixed.begin(), ite = vert_pos_2_indices_fixed.end();
		InputPart ip = BOUNDARY_PRESCRIBED_MOTION_NODE;
		int i, id;
		for (it = itb; it != ite; ++it)
		{
			i = it->first;
			id = allDomain_uFs[i].point.vertexID;
			getNumberSpaces(delim_id, id, ip);
			out << delim_id << id;
			getNumberSpaces(delim_id, 0, ip);
			out<< delim_id << "0" << delim_id<< it->second[0] << delim_id << it->second[1] << delim_id << it->second[2] << delim_id<< "0"<< delim_id<<"0"<< delim_id<<"0\n";
//			out << delim_id << id << "         0         " << it->second[0] << "         " << it->second[1] << "         " << it->second[2] << "         1         1         1\n";
		}
//		int id_left_bottom = mesh->corners[mn_t][mn_t].vertexID;
//		getNumberSpaces(delim_id, id_left_bottom, ip);
//		out << delim_id << id_left_bottom << "         0         1         1         1         1         1         1\n";
//		int id_right_bottom = mesh->corners[mx_t][mn_t].vertexID;
//		getNumberSpaces(delim_id, id_right_bottom, ip);
//		out << delim_id << id_right_bottom << "         0         0         1         1         1         1         1\n";
		
		string delim = " ";
		string delim3 = "   ";
		string delim4 = "    ";
		string delim19="                   ";
		string curve_name = bcConfig.is_static==true ? "_LS_DYNA_disp_load_curves.txt":"_LS_DYNA_disp_load_curves_dynamics.txt";
		string file_disp_load_curves = inputFolder4Run + curve_name;
		CopyContentFileA_ostream(file_disp_load_curves, out);
//		this->allDomain_uFs

		ip = BOUNDARY_PRESCRIBED_MOTION_NODE;
		out << "*BOUNDARY_PRESCRIBED_MOTION_NODE\n";
		for (unsigned int i = 0; i < sz_ptBCs; ++i)
		{
			for (int j = 0; j < DIM_U; ++j)
			{
				unsigned jp1 = j + 1;
				if (allDomain_uFs[i].op_ufs[j].bt == bt_dirichlet)
				{
					int id = allDomain_uFs[i].point.vertexID;
					// if(id==147)
					// 	std::cout<<"debug\n";
					double u = allDomain_uFs[i].op_ufs[j].u;

					string number_str;
					double abF = abs(u);
					toString(abF, number_str);
					bool IsZero = checkZeroString(number_str);
					if (IsZero)
						continue;
					if(abF<1.0e-9)
						continue;

	
					getNumberSpaces(delim_id, id,ip);

	
						out << "$#               nid                 dof                 vad                lcid                  sf                 vid               death               birth\n";
						// out << std::fixed;
						// out << std::setprecision(7); // show 16 digits of precision
						//if (abs(u) < 1.0e-10)
						//	out << delim_id << id << delim9 << jp1 << delim9 << "2" << delim9 << "1" << delim7 << "0.0" << delim9 << "01.00000E28" << delim7 << "0.0\n"; 0    9.9999994421e+27
						if (u>=0.0)
							out << delim_id << id << delim19 << jp1 << delim19 << "2" << delim19 << "1" << delim4 << std::fixed<< std::setprecision(14) <<u<< delim19<< "0"<< delim4<< "9.9999994421e+27"<< delim4 << "0.0000000000e+00\n";
						else 
							out << delim_id << id << delim19 << jp1 << delim19 << "2" << delim19 << "1" << delim3 <<std::fixed<< std::setprecision(14) <<u << delim19<< "0"<< delim4<<"9.9999994421e+27"<< delim4 << "0.0000000000e+00\n";

				}
			}
		}
		if (maxF != 0)
		{
			ip=LOAD_NODE_POINT;
			for (unsigned int i = 0; i < sz_ptBCs; ++i)
			{
				for (int j = 0; j < DIM_U; ++j)
				{
					unsigned jp1 = j + 1;
					if (allDomain_uFs[i].op_ufs[j].bt == bt_neumann)
					{
						int id = allDomain_uFs[i].point.vertexID;
						// if(id==11)
						// 	cout<<"debug\n";
						string number_str;
						double F = allDomain_uFs[i].op_ufs[j].F;
						double abF = abs(F);
						toString(abF, number_str);
						bool IsZero = checkZeroString(number_str);
						if (IsZero)
							continue;
						if(abF<1.0e-9)
							continue;

						out << "*LOAD_NODE_POINT\n";
						out << "$#     nid       dof      lcid        sf       cid        m1        m2        m3\n";
						// out << std::fixed;
						// out << std::setprecision(7); // show 16 digits of precision
						getNumberSpaces(delim_id, id, ip);
						out << delim_id << id;
						if(F>=0.0)
							out << delim19<< jp1 << delim19 << "2" << delim4 <<std::fixed<< std::setprecision(14)<< F << delim19 << "0" << delim19 << "0" << delim19 << "0" << delim19 << "0\n";
						else
							out << delim19<< jp1 << delim19 << "2" << delim3 <<std::fixed<< std::setprecision(14)<<  F << delim19 << "0" << delim19 << "0" << delim19 << "0" << delim19 << "0\n";
	
					}
				}
			}
		}

		// printing nodes for which the output is needed
		ip = DATABASE_HISTORY_NODE_ID;
		out << "*DATABASE_HISTORY_NODE_ID\n";
		out << "$#               id1\n";
//		out << "$# id1  heading\n";
		for (unsigned int i = 0; i < sz_ptBCs; ++i)
		{
			int id = allDomain_uFs[i].point.vertexID;
			getNumberSpaces(delim_id, id, ip);
			out << delim_id << id;
			out << '\n';
		}

		out << "*SET_NODE_LIST\n";
		out << "$#               sid                 da1                 da2                 da3                 da4              solver                 its                   -\n";
		out << delim19 << 1 << delim4 << "0.0000000000e+00" << delim4<< "0.0000000000e+00" << delim4 << "0.0000000000e+00" << delim4<< "0.0000000000e+00"
		              << "MECH" << "                "<< 1 << "\n";
		out << "$#              nid1                nid2                nid3                nid4                nid5                nid6                nid7                nid8\n";

//		out << "$# id1  heading\n";
        int quotient, remainder, szBC_tmp;
		quotient = sz_ptBCs/8;
		szBC_tmp = quotient*8;
		remainder = sz_ptBCs%8;
		ip = SET_NODE_LIST;
		for (unsigned int i = 0; i < szBC_tmp; i=i+8)
		{
			int id;
			for(int ii=0; ii< 8; ++ii)
			{
				id = allDomain_uFs[i+ii].point.vertexID;
				getNumberSpaces(delim_id, id, ip);
				out << delim_id << id;
			}
			out << '\n';
		}
		if (szBC_tmp< sz_ptBCs)
		{
			for(int i=szBC_tmp; i< sz_ptBCs; ++i)
			{
				id = allDomain_uFs[i].point.vertexID;
				getNumberSpaces(delim_id, id, ip);
				out << delim_id << id;			
			}

			for(int i=sz_ptBCs; i< (quotient+1)*8;++i)
				out<< delim19 << 0;
			out<< '\n';
		}

	}
}

void AllLoadCases_MeanStnStrs::AllLoadCases_MeanStnStrs_Write(ostream& out)
{
	out << "num_timeSteps\t" << num_timeSteps << '\n';
	for (unsigned int ti = 0; ti < num_timeSteps; ++ti)
	{
		out << "timeIndex\t" << ti << '\n';
		for (unsigned int caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
		{
			vector<double>*strn = &lc_stn_sts[caseNo].voigtStrains[ti];
			vector<double>*strs = &lc_stn_sts[caseNo].voigtStresses[ti];

			out << "\tcaseNo\t" << caseNo;
			out << "\tstrains";
			for (unsigned int i = 0; i < VOIGHT_STN_SZ; ++i)
				out << "\t" << (*strn)[i];
			out << "\tstresses";
			for (unsigned int i = 0; i < VOIGHT_STN_SZ; ++i)
				out << "\t" << (*strs)[i];
			out << "\n";
		}
		out << "stiffness\n";
		out << stiffnesses[ti] << '\n';
	}
}

void sMesh::Read(const string & fileNameWOExtIn, FEMSolverT solverOptionOut_In, FEMSolverT solverOptionIn_In)
{
	fileNameWOExt = fileNameWOExtIn;
	solverOptionIn = solverOptionIn_In;
	solverOptionOut = solverOptionOut_In;
	fileNamesOut.resize(VOIGHT_STN_SZ);
	string name;
	if (solverOptionIn == fems_abaqus)
		extIn = "inp";
	else if (solverOptionIn == fems_LSDYNA)
		extIn = "k";
	if (solverOptionOut == fems_abaqus)
		extOut = "inp";
	else if (solverOptionOut == fems_LSDYNA)
		extOut = "k";

	name = g_rootFolder_4Output +fileNameWOExt + "." + extIn;
	for (unsigned int caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
	{
		string str;
		toString(caseNo, str);
//		string dir;
		string subfolderName = g_rootFolder_4Output + fileNameWOExt + str;

		MakeDir(subfolderName);
		fileNamesOut[caseNo] = subfolderName + "/input." + extOut;
//		fileNamesOut[caseNo] = fileNameWOExt + "_out_" + str + "." + extOut;
	}
	fstream in(name.c_str(), ios::in);
	if (!in.is_open())
	{
		cout << name << '\n';
		THROW("Cannot open the file\n");
	}
	Read(in);
	Initialize();
}

void sMesh::Read(istream& in)
{
	string buf = "";
	while (buf != "*NODE")
		READ_NSTRING(in, buf, buf);
	in >> buf;
	long id;
	nodes.clear();
	while (fromString(buf, id))
	{
		smNode node;
		node.id = id;
		in >> node.crd[0];
		READ_NSTRING(in, buf, buf);
		in >> node.crd[1];
//		for (int j = 0; j < DIM; ++j)
//		{
//			in >> buf;
//			in >> node.crd[j];
//		}
		READ_NSTRING(in, buf, buf);
		nodes[node.id] = node;
//		nodes.push_back(node);
	}
	n_nodes = nodes.size();
//	READ_NSTRING(in, buf, buf);

	elss.clear();
	num_elsets = 0;

	while (!in.eof())
	{
		unsigned int sz = buf.size();
		while (sz == 0)
			READ_NSTRING(in, buf, buf);
		bool hasStar = buf[0] = '*';
		bool has2Star = (hasStar && (sz >= 2) && (buf[1] == '*'));
		if (has2Star)
		{
			THROW("Error READ_NSTRING should pass comment lines\n");
			getline(in, buf);
			READ_NSTRING(in, buf, buf);
		}
		if (hasStar)
		{
			if (buf.compare(0, 8, "*ELEMENT") == 0)
			{
				READ_NSTRING(in, buf, buf);
				if (buf.find("CPS4"))
				{
					in >> buf;
					long id;
					while (fromString(buf, id))
					{
						elementCP4 ele;
						ele.e4_id = id;
						in >> ele.vs4[0];
						in >> buf;
						in >> ele.vs4[1];
						in >> buf;
						in >> ele.vs4[2];
						in >> buf;
						in >> ele.vs4[3];
						READ_NSTRING(in, buf, buf);
						e4.push_back(ele);
						id2_e4Pos[id] = e4.size() - 1;
					}
					n_e4 = e4.size();
				}
			}
			else if (buf.compare(0, 5, "*NSET") == 0)
			{
				getline(in, buf);
				READ_NSTRING(in, buf, buf);
				while (fromString(buf, id))
					READ_NSTRING(in, buf, buf);
			}
			else if (buf.compare(0, 6, "*ELSET") == 0)
			{
				READ_NSTRING(in, buf, buf);
				std::string delimiter = "=";
				size_t pos = 0;
				std::string token;
				pos = buf.find(delimiter);
				if (pos < 0)
				{
					THROW("Cannot find =\n");
				}
				token = buf.substr(pos + 1, buf.size());// buf.substr(0, pos);
				elset els;
				transform(token.begin(), token.end(), token.begin(), ::toupper);
				els.name = token;
				READ_NSTRING(in, buf, buf);
				while (fromString(buf, id))
				{
					els.el_ids.push_back(id);
					READ_NSTRING(in, buf, buf);
				}
				elss.push_back(els);
				num_elsets = elss.size();
				unsigned int elSet_index = num_elsets - 1;
				unsigned int sz_set = els.el_ids.size();
				for (unsigned ii = 0; ii < sz_set; ++ii)
				{
					long id = els.el_ids[ii];
					map<long, long>::iterator it = id2_e4Pos.find(id);
					if (it == id2_e4Pos.end())
					{
						THROW("This is not necessarily an error as we may have triangular elements too\n");
					}
					int pos = it->second;
					e4[pos].matID = elSet_index;
				}
			}
			else if (buf.compare(0, 6, "*SOLID") == 0)
			{
				getline(in, buf);
				READ_NSTRING(in, buf, buf);
			}
			else if (buf.compare(0, 9, "*MATERIAL") == 0)
			{
				getline(in, buf);
				getline(in, buf);
				READ_NSTRING(in, buf, buf);
			}
		}
		else
		{
			cout << "buf\n" << buf << '\n';
			THROW("Invalid buf\n");
		}
	}
}

void sMesh::Initialize()
{
	long sz = e4.size();
	area_e4 = 0.0;
	for (unsigned int i = 0; i < sz; ++i)
	{
		elementCP4* ePtr = &e4[i];
		vector<smNode*> nodePtrs(4);
		for (int vi = 0; vi < 4; ++vi)
			nodePtrs[vi] = &nodes[ePtr->vs4[vi]];
		ePtr->area = CalculateQuadArea(nodePtrs);
		area_e4 += ePtr->area;
	}
	string area_out = g_rootFolder_4Output + "/areas.txt";
	fstream aout(area_out.c_str(), ios::out);
	aout << "total_area_e4\t" << area_e4 << '\n';
	aout << "num_e4\t" << sz << '\n';
	for (unsigned int i = 0; i < sz; ++i)
	{
		elementCP4* ePtr = &e4[i];
		aout << ePtr->e4_id << '\t' << ePtr->area << '\n';
	}
	for (int i = 0; i < DIM; ++i)
	{
		minCrd[i] = P_INFINITY;
		maxCrd[i] = M_INFINITY;
	}

	map<long, smNode>::iterator niter, niter0 = nodes.begin(), nitere = nodes.end();
	for (niter = niter0; niter != nitere; ++niter)
	{
		smNode* nPtr = &niter->second;
		for (int j = 0; j < DIM; ++j)
		{
			minCrd[j] = MIN(minCrd[j], nPtr->crd[j]);
			maxCrd[j] = MAX(maxCrd[j], nPtr->crd[j]);
		}
	}
	volume = 1.0;
	double minEdge = P_INFINITY;
	for (int j = 0; j < DIM; ++j)
	{
		meanCrd[j] = 0.5 * (maxCrd[j] + minCrd[j]);
		spanCrd[j] = maxCrd[j] - minCrd[j];
		minEdge = MIN(minEdge, spanCrd[j]);
		volume *= spanCrd[j];
	}
	double tol = 1e-6 * minEdge;
	double mnmx[DIM][mnmxT_SIZE];
	for (int i = 0; i < DIM; ++i)
	{
		mnmx[i][mn_t] = minCrd[i];
		mnmx[i][mx_t] = maxCrd[i];
	}

	map<double, int> crd2VertInd[DIM][mnmxT_SIZE];
	//map<long, smNode>::iterator 
	niter0 = nodes.begin(), nitere = nodes.end();
	for (niter = niter0; niter != nitere; ++niter)
	{
		smNode* nPtr = &niter->second;
		for (int j = 0; j < DIM; ++j)
		{
			double crdv = nPtr->crd[j];
			// only for 2D
			double crdOther = nPtr->crd[1 - j];
			for (int k = 0; k < mnmxT_SIZE; ++k)
			{
				if (fabs(crdv - mnmx[j][k]) < tol)
				{
					crd2VertInd[j][k][crdOther] = niter->first;
				}
			}
		}
	}
	for (int j = 0; j < DIM; ++j)
	{
		unsigned int posCrdChanging = 1 - j;
		for (int k = 0; k < mnmxT_SIZE; ++k)
		{
			map<double, int>* mp = &crd2VertInd[j][k];
			sSide *side = &sides[j][k];
			side->normal[j] = 1.0;
			if (k == mn_t)
				side->normal[j] = -1.0;
			side->normal[posCrdChanging] = 0.0;
			side->posCrdChanging = posCrdChanging;
			side->setSize(mp->size());
			map<double, int>::iterator it, itb = mp->begin(), ite = mp->end();
			unsigned int cntr = 0;
			for (it = itb; it != ite; ++it)
			{
//				int index = it->second;
				int nodeID = it->second;
//				side->vertexIndices[cntr] = nodeID;
				side->vertexIDs[cntr] = nodeID; // nodes[index].id;
				for (int m = 0; m < DIM; ++m)
					side->crdRel2Centroid[m][cntr] = nodes[nodeID].crd[m] - meanCrd[m];
				++cntr;
			}
			side->setSegmentLengthsAndWeights();
		}
	}
	for (unsigned int xmM = 0; xmM < mnmxT_SIZE; ++xmM)
	{
		sSide *side_x = &sides[0][xmM];
		for (unsigned int ymM = 0; ymM < mnmxT_SIZE; ++ymM)
		{
			sSide *side_y = &sides[1][ymM];
			sPoint *corner = &corners[xmM][ymM];
			corner->pos[0] = mnmxT(xmM);
			corner->pos[1] = mnmxT(ymM);
			unsigned pos_x = 0;
			if (ymM == mx_t)
				pos_x = side_x->sz - 1;
			unsigned pos_y = 0;
			if (xmM == mx_t)
				pos_y = side_y->sz - 1;
//			corner->vertexIndex = side_x->vertexIndices[pos_x];
			corner->vertexID = side_x->vertexIDs[pos_x];
			corner->crdRel2Centroid[0] = side_x->crdRel2Centroid[0][pos_x];
			corner->crdRel2Centroid[1] = side_x->crdRel2Centroid[1][pos_x];
			corner->normalTimesWeights[0] = side_x->weight[pos_x];
			if (pos_y == 0)
				corner->normalTimesWeights[0] *= -1;

			corner->normalTimesWeights[1] = side_y->weight[pos_y];
			if (pos_x == 0)
				corner->normalTimesWeights[1] *= -1;
		}
	}

	// forming boundary points
	sz_boundaryPoints = sides[0][mn_t].sz + sides[0][mx_t].sz + sides[1][mn_t].sz + sides[1][mx_t].sz - 4;
	boundaryPoints.resize(sz_boundaryPoints);
	unsigned int cntr = 0;
	boundaryPoints_cornerPoss[sm_bl] = cntr;
	boundaryPoints[cntr++] = corners[mn_t][mn_t];
	sSide* side = &sides[1][mn_t];
	for (int i = 1; i < side->sz - 1; ++i)
	{
		boundaryPoints[cntr].crdRel2Centroid[0] = side->crdRel2Centroid[0][i];
		boundaryPoints[cntr].crdRel2Centroid[1] = side->crdRel2Centroid[1][i];
		boundaryPoints[cntr].vertexID = side->vertexIDs[i];
//		boundaryPoints[cntr].vertexIndex = side->vertexIndices[i];
		boundaryPoints[cntr].normalTimesWeights[0] = side->normal[0] * side->weight[i];
		boundaryPoints[cntr].normalTimesWeights[1] = side->normal[1] * side->weight[i];
		boundaryPoints[cntr].pos[0] = mnx_none_t;
		boundaryPoints[cntr].pos[1] = mn_t;
		cntr++;
	}
	boundaryPoints_cornerPoss[sm_br] = cntr;
	boundaryPoints[cntr++] = corners[mx_t][mn_t];
	side = &sides[0][mx_t];
	for (int i = 1; i < side->sz - 1; ++i)
	{
		boundaryPoints[cntr].crdRel2Centroid[0] = side->crdRel2Centroid[0][i];
		boundaryPoints[cntr].crdRel2Centroid[1] = side->crdRel2Centroid[1][i];
		boundaryPoints[cntr].vertexID = side->vertexIDs[i];
//		boundaryPoints[cntr].vertexIndex = side->vertexIndices[i];
		boundaryPoints[cntr].normalTimesWeights[0] = side->normal[0] * side->weight[i];
		boundaryPoints[cntr].normalTimesWeights[1] = side->normal[1] * side->weight[i];
		boundaryPoints[cntr].pos[0] = mx_t;
		boundaryPoints[cntr].pos[1] = mnx_none_t;
		cntr++;
	}
	boundaryPoints_cornerPoss[sm_tr] = cntr;
	boundaryPoints[cntr++] = corners[mx_t][mx_t];
	side = &sides[1][mx_t];
	for (int i = side->sz - 2; i > 0; --i)
	{
		boundaryPoints[cntr].crdRel2Centroid[0] = side->crdRel2Centroid[0][i];
		boundaryPoints[cntr].crdRel2Centroid[1] = side->crdRel2Centroid[1][i];
		boundaryPoints[cntr].vertexID = side->vertexIDs[i];
//		boundaryPoints[cntr].vertexIndex = side->vertexIndices[i];
		boundaryPoints[cntr].normalTimesWeights[0] = side->normal[0] * side->weight[i];
		boundaryPoints[cntr].normalTimesWeights[1] = side->normal[1] * side->weight[i];
		boundaryPoints[cntr].pos[0] = mnx_none_t;
		boundaryPoints[cntr].pos[1] = mx_t;
		cntr++;
	}
	boundaryPoints_cornerPoss[sm_tl] = cntr;
	boundaryPoints[cntr++] = corners[mn_t][mx_t];
	side = &sides[0][mn_t];
	for (int i = side->sz - 2; i > 0; --i)
	{
		boundaryPoints[cntr].crdRel2Centroid[0] = side->crdRel2Centroid[0][i];
		boundaryPoints[cntr].crdRel2Centroid[1] = side->crdRel2Centroid[1][i];
		boundaryPoints[cntr].vertexID = side->vertexIDs[i];
//		boundaryPoints[cntr].vertexIndex = side->vertexIndices[i];
		boundaryPoints[cntr].normalTimesWeights[0] = side->normal[0] * side->weight[i];
		boundaryPoints[cntr].normalTimesWeights[1] = side->normal[1] * side->weight[i];
		boundaryPoints[cntr].pos[0] = mn_t;
		boundaryPoints[cntr].pos[1] = mnx_none_t;
		cntr++;
	}
	sz_boundaryPoints = boundaryPoints.size();
}

void sMesh::Compute_AllStrain_BCs(const vector<double>& voigtStrain, boundary_AllPoints_uFAllDirs_1Case& bcs)
{
	bcs.allDomain_uFs.resize(sz_boundaryPoints);
	bcs.sz_ptBCs = sz_boundaryPoints;
	double eps00 = voigtStrain[0];
#if !ANTI_PLANE_SHEAR
	double eps01 = 0.5 * voigtStrain[2];
	double eps11 = voigtStrain[1];
#else
	double eps01 = voigtStrain[1];
#endif
	for (int i = 0; i < sz_boundaryPoints; ++i)
	{
		bcs.allDomain_uFs[i].point = boundaryPoints[i];
		double x0 = bcs.allDomain_uFs[i].point.crdRel2Centroid[0];
		double x1 = bcs.allDomain_uFs[i].point.crdRel2Centroid[1];
		bcs.allDomain_uFs[i].op_ufs[0].u = x0 * eps00 + x1 * eps01;
		bcs.allDomain_uFs[i].op_ufs[0].bt = bt_dirichlet;
#if !ANTI_PLANE_SHEAR
		bcs.allDomain_uFs[i].op_ufs[1].u = x0 * eps01 + x1 * eps11;
		bcs.allDomain_uFs[i].op_ufs[1].bt = bt_dirichlet;
#endif
	}
}

void sMesh::Compute_AllTraction_BCs(const vector<double>& voigtStress, boundary_AllPoints_uFAllDirs_1Case& bcs, unsigned int caseNo)
{
	bcs.allDomain_uFs.resize(sz_boundaryPoints);
	bcs.sz_ptBCs = sz_boundaryPoints;
	double sig00 = voigtStress[0];
#if !ANTI_PLANE_SHEAR
	double sig01 = voigtStress[2];
	double sig11 = voigtStress[1];
#else
	double sig01 = voigtStress[1];
#endif
	for (int i = 0; i < sz_boundaryPoints; ++i)
	{
		bcs.allDomain_uFs[i].point = boundaryPoints[i];
		double nTilde0 = bcs.allDomain_uFs[i].point.normalTimesWeights[0];
		double nTilde1 = bcs.allDomain_uFs[i].point.normalTimesWeights[1];
		bcs.allDomain_uFs[i].op_ufs[0].F = nTilde0 * sig00 + nTilde1 * sig01;
		bcs.allDomain_uFs[i].op_ufs[0].bt = bt_neumann;
#if !ANTI_PLANE_SHEAR
		bcs.allDomain_uFs[i].op_ufs[1].F = nTilde0 * sig01 + nTilde1 * sig11;
		bcs.allDomain_uFs[i].op_ufs[1].bt = bt_neumann;
#endif
	}
	// handing rigid motions
	if (caseNo < 0)
		return;
#if !ANTI_PLANE_SHEAR
	if ((caseNo == 0) || (caseNo == 1)) // loading in x or y direction
	{
		// first position: bottom edge first after left corner
		int posFixed = boundaryPoints_cornerPoss[sm_bl];
		//		int posFixed = 1;
		//		if (caseNo == 1)
		//		{
					// first position: left edge first above bottom corner
		//			posFixed = sz_boundaryPoints - 1;
		//		}
		bcs.allDomain_uFs[posFixed].op_ufs[0].u = 0.0;
		bcs.allDomain_uFs[posFixed].op_ufs[1].u = 0.0;
		bcs.allDomain_uFs[posFixed].op_ufs[0].bt = bt_dirichlet;
		bcs.allDomain_uFs[posFixed].op_ufs[1].bt = bt_dirichlet;
		bcs.allDomain_uFs[posFixed].op_ufs[0].F = 0.0;
		bcs.allDomain_uFs[posFixed].op_ufs[1].F = 0.0;

		// now adding the roller
		if (caseNo == 0)
		{
			int posRoller = boundaryPoints_cornerPoss[sm_br]; //sides[1][mn_t].sz - 1;
			bcs.allDomain_uFs[posRoller].op_ufs[1].u = 0.0;
			bcs.allDomain_uFs[posRoller].op_ufs[1].bt = bt_dirichlet;
			bcs.allDomain_uFs[posRoller].op_ufs[1].F = 0.0;
		}
		else
		{
			int posRoller = boundaryPoints_cornerPoss[sm_tl];//bcs.allDomain_uFs.size() - sides[0][mn_t].sz + 1;
			bcs.allDomain_uFs[posRoller].op_ufs[0].u = 0.0;
			bcs.allDomain_uFs[posRoller].op_ufs[0].bt = bt_dirichlet;
			bcs.allDomain_uFs[posRoller].op_ufs[0].F = 0.0;
		}
	}
	else
	{
		// shear loading
		int posFixed = boundaryPoints_cornerPoss[sm_bl];

		bcs.allDomain_uFs[posFixed].op_ufs[0].u = 0.0;
		bcs.allDomain_uFs[posFixed].op_ufs[1].u = 0.0;
		bcs.allDomain_uFs[posFixed].op_ufs[0].bt = bt_dirichlet;
		bcs.allDomain_uFs[posFixed].op_ufs[1].bt = bt_dirichlet;
		// now adding the corresponding force to neighbor node
//		bcs.allDomain_uFs[1].op_ufs[0].F += bcs.allDomain_uFs[posFixed].op_ufs[0].F;
		bcs.allDomain_uFs[posFixed].op_ufs[0].F = 0.0;
		//		bcs.allDomain_uFs[sz_boundaryPoints - 1].op_ufs[1].F += bcs.allDomain_uFs[posFixed].op_ufs[1].F;
		bcs.allDomain_uFs[posFixed].op_ufs[1].F = 0.0;

		int posRoller = boundaryPoints_cornerPoss[sm_br]; //sides[1][mn_t].sz - 1;
		bcs.allDomain_uFs[posRoller].op_ufs[1].u = 0.0;
		bcs.allDomain_uFs[posRoller].op_ufs[1].bt = bt_dirichlet;
		bcs.allDomain_uFs[posRoller].op_ufs[1].F = 0.0;
	}
#else
	{
		int posFixed = boundaryPoints_cornerPoss[sm_bl];
		//		int posFixed = 1;
		//		if (caseNo == 1)
		//		{
					// first position: left edge first above bottom corner
		//			posFixed = sz_boundaryPoints - 1;
		//		}
		bcs.allDomain_uFs[posFixed].op_ufs[0].u = 0.0;
		bcs.allDomain_uFs[posFixed].op_ufs[0].bt = bt_dirichlet;
		bcs.allDomain_uFs[posFixed].op_ufs[0].F = 0.0;
	}
#endif
}

void sMesh::FormEmpty_BoundaryStorage(boundary_AllPoints_uFAllDirs_1Case& bcs)
{
	bcs.sz_ptBCs = sz_boundaryPoints;
	bcs.allDomain_uFs.resize(sz_boundaryPoints);

	for (int i = 0; i < sz_boundaryPoints; ++i)
	{
		bcs.allDomain_uFs[i].point = boundaryPoints[i];
		bcs.vertexID2Pos[boundaryPoints[i].vertexID] = i;
	}
}

void sMesh::Compute_GeneralMixed_BCs(BCInfo& bcInfo, boundary_AllPoints_uFAllDirs_1Case& bcs)
{
	bcs.sz_ptBCs = sz_boundaryPoints;
	bcs.allDomain_uFs.resize(sz_boundaryPoints);
	double eps00 = bcInfo.voigtStrain[0];
	double sig00 = bcInfo.voigtStress[0];
#if !ANTI_PLANE_SHEAR
	double eps01 = 0.5 * bcInfo.voigtStrain[2];
	double eps11 = bcInfo.voigtStrain[1];
	double sig01 = bcInfo.voigtStress[2];
	double sig11 = bcInfo.voigtStress[1];
#else
	double eps01 = bcInfo.voigtStrain[1];
	double sig01 = bcInfo.voigtStress[1];
#endif
	double u[DIM_U], F[DIM_U];
	boundaryT bc[DIM_U];

	for (int i = 0; i < sz_boundaryPoints; ++i)
	{
		bcs.allDomain_uFs[i].point = boundaryPoints[i];
		boundary_1Point_uFAllDirs* b1p = &bcs.allDomain_uFs[i];
		sPoint* pointPtr = &b1p->point;
		double x0 = pointPtr->crdRel2Centroid[0];
		double x1 = pointPtr->crdRel2Centroid[1];
		double nTilde0 = pointPtr->normalTimesWeights[0];
		double nTilde1 = pointPtr->normalTimesWeights[1];
		u[0] = x0 * eps00 + x1 * eps01;
		F[0] = nTilde0 * sig00 + nTilde1 * sig01;
#if !ANTI_PLANE_SHEAR
		u[1] = x0 * eps01 + x1 * eps11;
		F[1] = nTilde0 * sig01 + nTilde1 * sig11;
#endif

		for (int ucompi = 0; ucompi < DIM_U; ++ucompi)
		{
			bc[ucompi] = bt_none;
			for (int sidei = 0; sidei < DIM; ++sidei)
			{
				mnmxT mnx = pointPtr->pos[sidei];
				if (mnx != mnx_none_t)
					UpdateBCType(bc[ucompi], bcInfo.side_BCType[sidei][mnx][ucompi]);
			}
			b1p->op_ufs[ucompi].bt = bc[ucompi];
//			if (bc[ucompi] == bt_dirichlet)
				b1p->op_ufs[ucompi].u = u[ucompi];
//			else
				b1p->op_ufs[ucompi].F = F[ucompi];
		}
	}
	int sz = bcInfo.corner2Dir_4RigidControl.size();
	for (int i = 0; i < sz; ++i)
	{
		cornerT cor = bcInfo.corner2Dir_4RigidControl[i].first;
#if !ANTI_PLANE_SHEAR
		int dir = bcInfo.corner2Dir_4RigidControl[i].second;
#else
		int dir = 0;
#endif
		int pos = boundaryPoints_cornerPoss[cor];
		bcs.allDomain_uFs[pos].op_ufs[dir].bt = bt_dirichlet;
	}
}

void sMesh::ComputePrintBC_AllCases_Aux(BCConfig& bcConfig, const vector<string>& fileNames, FEMSolverT solverOptionIn, vector<boundary_AllPoints_uFAllDirs_1Case>& bcs)
{
	bcs.resize(VOIGHT_STN_SZ);

	if (solverOptionOut == fems_abaqus)
	{
		op command = copyF;
		string source = fileNameWOExt + "." + extIn;
		string extinpt = "inp";
		string file_before = "before_BC." + extinpt;
		string file_after = "after_BC." + extinpt;
		for (unsigned int caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
		{
			fileOperation(command, source, fileNames[caseNo]);
			fstream out(fileNames[caseNo].c_str(), ios::app);
			CopyContentFileA_ostream(file_before, out);
			ComputePrintBC(bcConfig, caseNo, solverOptionIn, out, bcs[caseNo], bcConfig.inputFolder4Run);
			CopyContentFileA_ostream(file_after, out);
		}
	}


	else if (solverOptionOut == fems_LSDYNA)
	{
		string header_name = bcConfig.is_static==true ? "_LS_DYNA_header.txt":"_LS_DYNA_header_dynamics.txt";
		string delim19="                   ";
		string delim18="                  ";
		string delim17="                 ";
		string delim16="                ";
		string delim15="               ";
		string delim13="             ";
		string delim4="    ";

		for (unsigned int caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
		{
			string glob_name = inc_dir+"global_params.k" ;
			fstream out_glob(glob_name.c_str(), ios::out);
			out_glob<< "*PARAMETER\n";
			out_glob<< "$#"<<delim18<<"0"<<delim19<<"0"<< delim19<<"0"<<delim19<<"0"<<delim19<<"0"<< delim19<<"0"<<delim19<<"0"<< delim19<<"\n"; 
			out_glob<< "$#             prmr1                val1               prmr2                val2               prmr3                val3               prmr4                val4\n";
			out_glob<< "RAX"<< delim17<<"10.0"<< delim16<<"RAY"<< delim17<<"10.0"<<delim16<< "RST"<< delim17<< "1.0"<< delim17<<"RTF"<< delim17<< "0.10"<< "\n";
			out_glob<< "RSO"<< delim17<<"0.050"<<delim15<<"RDT"<<delim17<<  "0.001"<< delim15<< "RDTPLOT"<< delim13<<"0.001"<<"\n";
			out_glob.close();


			string out_name = inc_dir+"output_params.k" ;
			fstream out_file(out_name.c_str(), ios::out);
			out_file<< "*DATABASE_BNDOUT\n";
			out_file<<"$#                dt              binary                lcur               ioopt             option1             option2             option3             option4\n";
			out_file<<delim17<<"&DT"<<delim19<<"3"<<delim19<<"0"<<delim19<<"1"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<"\n";

			out_file<< "*DATABASE_ELOUT\n";
			out_file<<"$#                dt              binary                lcur               ioopt             option1             option2             option3             option4\n";
			out_file<<delim17<<"&DT"<<delim19<<"3"<<delim19<<"0"<<delim19<<"1"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<"\n";

			out_file<< "*DATABASE_BINARY_D3PLOT\n";
			out_file<<"$#                dt                lcdt                beam               npltc              psetid\n";
			out_file<<delim13<<"&DTPLOT"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0\n";
			out_file<<"$#             ioopt                rate              cutoff              window                type                pset\n";
	        out_file<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0\n";                                                                              

			out_file<<"*DATABASE_GLSTAT\n";
			out_file<<"$#                dt              binary                lcur               ioopt\n";
			out_file<< delim17<<"&DT"<<delim19<<"3"<<delim19<<"0"<<delim19<<"1\n";

			out_file<<"*DATABASE_NODFOR\n";
			out_file<<"$#                dt              binary                lcur               ioopt\n";
			out_file<< delim17<<"&DT"<<delim19<<"3"<<delim19<<"0"<<delim19<<"1\n";

			out_file<<"*DATABASE_NODOUT\n";
			out_file<<"$#                dt              binary                lcur               ioopt             option1             option2\n";
			out_file<< delim17<<"&DT"<<delim19<<"3"<<delim19<<"0"<<delim19<<"1"<<delim4<<"0.0000000000e+00"<<delim19<< "0" <<"\n";

			out_file<<"*DATABASE_SPCFORC\n";
			out_file<<"$#                dt              binary                lcur               ioopt\n";
			out_file<< delim17<<"&DT"<<delim19<<"3"<<delim19<<"0"<<delim19<<"1\n";

			out_file<<"*DATABASE_HISTORY_SHELL_SET\n";
			out_file<<"$#               id1                 id2                 id3                 id4                 id5                 id6                 id7                 id8\n";
			out_file<<delim19<<"1"<<delim19<<"2\n";  


			out_file<<"*DATABASE_EXTENT_BINARY\n";
			out_file<<"$#             neiph               neips              maxint              strflg              sigflg              epsflg              rltflg              engflg\n";
			out_file<<delim19<<"0"<<delim19<<"4"<<delim19<<"3"<<delim19<<"0"<<delim19<<"1"<<delim19<<"1"<<delim19<<"1"<<delim19<<"1"<<"\n";
			out_file<<"$#            cmpflg              ieverp              beamip               dcomp                shge               stssz              n3thdt             ialemat\n";
			out_file<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"1"<<delim19<<"1"<<delim19<<"1"<<delim19<<"2"<<delim19<<"1\n";
			out_file<<"$#           nintsld             pkp_sen                sclp               hydro               msscl               therm              intout              nodout\n";
			out_file<<delim19<<"0"<<delim19<<"0"<<delim17<<"1.0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0\n";                                        
			out_file<<"$#              dtdt              resplt               neipb             quadsld              cubsld             deleres\n";
			out_file<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0\n";

			out_file.close();


			fstream out(fileNames[caseNo].c_str(), ios::out);		
			out<<"*KEYWORD    LONG=Y    memory=30000000\n";
			out<<"*INCLUDE\n";   
			out<< "../../../../global_params.k" << "\n";
			out<< "../../../../output_params.k" << "\n";
			//out<< inc_dir+"initial_veloc.k\n";
			out<<"*TITLE\n";
			string fileHeader = bcConfig.inputFolder4Run + header_name;
			CopyContentFileA_ostream(fileHeader, out);


			ComputePrintBC(bcConfig, caseNo, solverOptionIn, out, bcs[caseNo], bcConfig.inputFolder4Run);
			PrintNode_Element_ElementSet_LSDYNA(out);
			out << "*END\n";
		}
	}
}

void sMesh::ComputePrintBC(BCConfig& bcConfig, unsigned int caseNo, FEMSolverT solverOptionIn, ostream& out, boundary_AllPoints_uFAllDirs_1Case& bcs, string inputFolder4Run)
{
	Compute_GeneralMixed_BCs(bcConfig.bcInfoCases[caseNo], bcs);
	if (solverOptionIn == fems_none)
		return;
	bcs.PrintBC(bcConfig, solverOptionIn, out, bcConfig.inputFolder4Run, this);
}

void sMesh::Compute_Average_Strain_Stress_From_BoundaryNode_Disps_and_Forces(const boundary_AllPoints_uFAllDirs_1Case& bcs, vector<double>& voigtStrain, vector<double>& voigtStress)
{
	unsigned int sz = bcs.sz_ptBCs;
	double inverseVol = 1.0 / volume;
	voigtStrain.resize(VOIGHT_STN_SZ);
	voigtStress.resize(VOIGHT_STN_SZ);
#if ANTI_PLANE_SHEAR
		THROW("implement this part\n");
#else
	MATRIX stn(DIM_U, DIM_U), sts(DIM_U, DIM_U);
	stn = 0.0;
	sts = 0.0;
	for (unsigned int ni = 0; ni < sz; ++ni)
	{
		const boundary_1Point_uFAllDirs* pt = &bcs.allDomain_uFs[ni];
		for (int i = 0; i < DIM_U; ++i)
		{
			for (int j = 0; j < DIM_U; ++j)
			{
				stn[i][j] += pt->op_ufs[i].u * pt->point.normalTimesWeights[j];
				sts[i][j] += pt->op_ufs[i].F * pt->point.crdRel2Centroid[j];
			}
		}
//		cout << "stn\n" << stn << '\n';
//		cout << "sts\n" << sts << '\n';
		FullElasticity_StnStnTensor2VoigtVector(stn, voigtStrain, true, inverseVol);
		FullElasticity_StnStnTensor2VoigtVector(sts, voigtStress, false, inverseVol);
	}
#endif
}

int sMesh::Compute_Average_Strain_Stress_From_BoundaryNode_Disps_and_Forces(const vector<boundary_AllPoints_uFAllDirs_1Case>& bcss, vector< vector<double> >& voigtStrains, vector< vector<double> >& voigtStresses, bool onlyIncludeLastStep)
{
	unsigned int sz = bcss.size();
	if (sz == 0)
		return 0;
	if (onlyIncludeLastStep || (sz == 1))
	{
		voigtStrains.resize(1);
		voigtStresses.resize(1);
		Compute_Average_Strain_Stress_From_BoundaryNode_Disps_and_Forces(bcss[sz - 1], voigtStrains[0], voigtStresses[0]);
		return 1;
	}
	voigtStrains.resize(sz);
	voigtStresses.resize(sz);
	for (unsigned int i = 0; i < sz; ++i)
		Compute_Average_Strain_Stress_From_BoundaryNode_Disps_and_Forces(bcss[i], voigtStrains[i], voigtStresses[i]);
	return sz;
}

void sMesh::PrintNode_Element_ElementSet_LSDYNA(ostream& out)
{
	string delim4 = "    ";
	string delim6 = "    ";
	string delim8 = "      ";
	string delim2="  ";

	string delim19="                   ";
	string delim18="                  ";
	string delim17="                 ";
	string delim_tmp;

    out<<"*PART\n";
	out<<"$#\n";                  
	out<<"BLACK\n";
	out<<"$#               pid               secid                 mid               eosid                hgid                grav              adpopt                tmid\n";
    out<<delim19<<"1"<<delim19<<"1"<<delim19<<"1"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0\n";
	out<<"*SECTION_SHELL_TITLE\n";
	out<<"BLACK\n";
	out<<"$#             secid              elform                shrf                 nip               propt             qr/irid               icomp               setyp\n";
    out<<delim19<<"1"<<delim18<<"13"<<delim4<<"0.0000000000e+00"<<delim19<<"4"<<delim17<<"1.0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"1\n";
	out<<"$#                t1                  t2                  t3                  t4                nloc               marea                idof              edgset\n";
    out<<delim17<<"1.0"<<delim17<<"1.0"<<delim17<<"1.0"<<delim17<<"1.0"<<delim4<<"0.0000000000e+00"<<delim4<<"0.0000000000e+00"<<delim4<<"0.0000000000e+00"<<delim19<<"0\n";



    out<<"*PART\n";
	out<<"$#\n";                  
	out<<"GREY\n";
	out<<"$#               pid               secid                 mid               eosid                hgid                grav              adpopt                tmid\n";
    out<<delim19<<"2"<<delim19<<"2"<<delim19<<"2"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"0\n";
	out<<"*SECTION_SHELL_TITLE\n";
	out<<"GREY\n";
	out<<"$#             secid              elform                shrf                 nip               propt             qr/irid               icomp               setyp\n";
    out<<delim19<<"2"<<delim18<<"13"<<delim4<<"0.0000000000e+00"<<delim19<<"4"<<delim17<<"1.0"<<delim19<<"0"<<delim19<<"0"<<delim19<<"1\n";
	out<<"$#                t1                  t2                  t3                  t4                nloc               marea                idof              edgset\n";
    out<<delim17<<"1.0"<<delim17<<"1.0"<<delim17<<"1.0"<<delim17<<"1.0"<<delim4<<"0.0000000000e+00"<<delim4<<"0.0000000000e+00"<<delim4<<"0.0000000000e+00"<<delim19<<"0\n";


	out << "*ELEMENT_SHELL\n";
//	out << "$# eid pid n1 n2 n3 n4 n5 n6 n7 n8\n";
	out << "$#               eid                 pid                  n1                  n2                  n3                  n4                  n5                  n6                  n7                  n8\n";
	InputPart ip = ELEMENT_SHELL;
	long sz = e4.size();
	for (unsigned int i = 0; i < sz; ++i)
	{
		elementCP4* ePtr = &e4[i];
		int matID = ePtr->matID;
		if (matID == -1)
		{
			cout << "i\t" << i << "\tid\t" << ePtr->e4_id << '\n';
			THROW("Material ID is not set for this element\n");
		}
		getNumberSpaces(delim_tmp, ePtr->e4_id , ip);
		out << delim_tmp << ePtr->e4_id ;

		getNumberSpaces(delim_tmp, matID  , ip);
		out << delim_tmp << matID + 1;

		for (unsigned int j = 0; j < 4; ++j)
		{
			getNumberSpaces(delim_tmp,ePtr->vs4[j] , ip);
			out << delim_tmp << ePtr->vs4[j];
		}
			
		out << delim19 << "0" << delim19 << "0" << delim19 << "0" << delim19 << "0\n";
	}
	out << "*NODE\n";
//	out << "$# nid x y z tc rc\n";
	out << "$#               nid                   x                   y                   z                  tc                  rc\n";

	
	// print velocity
	//string vel_file =inc_dir+"initial_veloc.k";
	map<long, smNode>::iterator niter, niter0 = nodes.begin(), nitere = nodes.end();
	for (niter = niter0; niter != nitere; ++niter)
	{
		// out << std::fixed;
		// out << std::setprecision(8); // show 16 digits of precision

		smNode* nPtr = &niter->second;
		ip=NODE_ID;
		getNumberSpaces(delim_tmp, nPtr->id,  ip);
		out << delim_tmp << nPtr->id ;
		out<< delim8 <<std::fixed<< std::setprecision(12)<<  nPtr->crd[0];
		out<< delim8 <<std::fixed<< std::setprecision(12)<<  nPtr->crd[1];
		out << delim4 << "0.0000000000e+00" << delim19 << "0" << delim19 << "0\n";
	}


	//out_vel.close();
// fstream out_vel(vel_file.c_str(), ios::out);
// 	{
		out<< "*INITIAL_VELOCITY_NODE\n";
		out<< "$#               nid                  vx                  vy                  vz                 vxr                 vyr                 vzr                icid\n";
	//}

	niter0 = nodes.begin(), nitere = nodes.end();
	for (niter = niter0; niter != nitere; ++niter)
	{
		// out << std::fixed;
		// out << std::setprecision(8); // show 16 digits of precision

		smNode* nPtr = &niter->second;
		ip=VEL_ID;
		getNumberSpaces(delim_tmp, nPtr->id,  ip);
		out<< delim_tmp<< nPtr->id;
		ip=NODE_VEL;
		out<<delim2<< "&AX*"<< nPtr->crd[0];

		getNumberSpaces(delim_tmp, nPtr->crd[1],  ip);
		out <<delim2<<"&AY*"<< nPtr->crd[1];
		out<< delim4<< "0.0000000000e+00"<< delim4<< "0.0000000000e+00"<< delim4<< "0.0000000000e+00"<< delim4<< "0.0000000000e+00"<< delim19<< "0"<<"\n";

		
	}


	ip=SET_SHELL_LIST_TITLE;
	for (unsigned int esi = 0; esi < num_elsets; ++esi)
	{
		elset* eSet = &elss[esi];
		out << "*SET_SHELL_LIST_TITLE\n";
		out << eSet->name << '\n';
		out << "$#               sid                 da1                 da2                 da3                 da4\n";
		out << delim19 << esi + 1 << delim4 << "0.0000000000e+00" << delim4 << "0.0000000000e+00" << delim4 << "0.0000000000e+00" << delim4 << "0.0000000000e+00\n";
		out << "$#              eid1                eid2                eid3                eid4                eid5                eid6                eid7                eid8\n";
		sz = eSet->el_ids.size();
		for (unsigned int i = 0; i < sz; ++i)
		{
			getNumberSpaces(delim_tmp, eSet->el_ids[i],  ip);
			out << delim_tmp << eSet->el_ids[i];
			if ((i + 1) % 8 == 0)
				out << '\n';
		}
		unsigned int tmpi = sz  % 8;
		if (tmpi > 0)
		{
			for (unsigned int i = tmpi; i < 8; ++i)
				out << delim19 << 0;
			out << '\n';
		}
	}
}

bool sMesh::Homogenize_CFromAllLoadCases(const string& fileNameWOExtIn, FEMSolverT solverOptionOut, bool onlyIncludeLastStep, AllLoadCases_MeanStnStrs& homogResults)
{
	int numActiveSteps = 0;
	for (unsigned caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
	{
		string str;
		toString(caseNo, str);
		string subfolderName = g_rootFolder_4Output + fileNameWOExt + str + "/";
		bool successfulread;
		if (solverOptionOut == fems_LSDYNA)
		{
			string fileName_u = subfolderName + "nodout";
//			string fileName_F = subfolderName + "bndout";
			string fileName_F = subfolderName + "nodfor";

			numActiveSteps = Read_LS_DYNA_Boundary_Solutions_1LoadCase(fileName_u, fileName_F,
				homogResults.lc_stn_sts[caseNo], onlyIncludeLastStep, successfulread);
			if (!successfulread)
				return false;
		}
		else if (solverOptionOut == fems_abaqus)
		{
			THROW("Add a similar read of nodal u and F function\n");
		}
	}
	homogResults.num_timeSteps = numActiveSteps;
	homogResults.stiffnesses.resize(homogResults.num_timeSteps);


	for (unsigned int ti = 0; ti < numActiveSteps; ++ti)
	{
		MATRIX strns(VOIGHT_STN_SZ, VOIGHT_STN_SZ), strss(VOIGHT_STN_SZ, VOIGHT_STN_SZ), *stiff = &homogResults.stiffnesses[ti];
		stiff->resize(VOIGHT_STN_SZ, VOIGHT_STN_SZ);
		(*stiff) = 0.0;
		for (unsigned i = 0; i < VOIGHT_STN_SZ; ++i)
		{
			for (unsigned caseNo = 0; caseNo < VOIGHT_STN_SZ; ++caseNo)
			{
				strns[i][caseNo] = homogResults.lc_stn_sts[caseNo].voigtStrains[ti][i];
				strss[i][caseNo] = homogResults.lc_stn_sts[caseNo].voigtStresses[ti][i];
			}
		}

		MATRIX strnsInv;
		strns.ComputeInverse(strnsInv);
		Product(*stiff, strss, strnsInv);
	}
	string fileName = g_rootFolder_4Output+fileNameWOExtIn + ".homog";
	fstream out(fileName.c_str(), ios::out);
	homogResults.AllLoadCases_MeanStnStrs_Write(out);
	return true;
}

int sMesh::Read_LS_DYNA_Boundary_Solutions_1LoadCase(const string& fileName_u, const string& fileName_F,
	OneLoadCase_MeanStnStrs& lc_stn_sts, bool onlyIncludeLastStep, bool& successfulread)
{
	successfulread = true;	vector<boundary_AllPoints_uFAllDirs_1Case> bcss;
#if TST_HOMOG
	double eps00 = 1.24, eps11 = 2.34, eps01 = -1.75;
	double sig00 = -4.3, sig11 = 2.22, sig01 = -12.1;
#endif
	fstream inu(fileName_u.c_str(), ios::in);
	if (!inu.is_open())
	{
		cout << fileName_u << '\n';
		successfulread = false;
		return 0;
	//		THROW("Cannot open file\n");
	}
	unsigned int timeStepNoRead;
	double timeValueRead;
	double ux, uy;
	int vertID, vertIndex;
	string buf;
	while (!inu.eof())
	{
		if (Read_LSDYNA_NodalDisplacements(inu, timeStepNoRead, timeValueRead))
		{
			boundary_AllPoints_uFAllDirs_1Case bcs;
			FormEmpty_BoundaryStorage(bcs);
			map<int, int>::iterator it, ite = bcs.vertexID2Pos.end();

			for (unsigned int i = 0; i < sz_boundaryPoints; ++i)
			{
				inu >> vertID >> ux >> uy;
				getline(inu, buf);
				it = bcs.vertexID2Pos.find(vertID);
				if (it == ite)
				{
					cout << "vertID\t" << vertID << " this vertex does not belong to boundary set\n";
//					successfulread = false;
//					return 0;
					THROW("Invalid vertexID\n");
				}
				vertIndex = it->second;
				boundary_1Point_uFAllDirs* pt = &bcs.allDomain_uFs[vertIndex];
				pt->op_ufs[0].bt = bt_dirichlet;
#if TST_HOMOG
				ux = eps00 * pt->point.crdRel2Centroid[0] + eps01 * pt->point.crdRel2Centroid[1];
				uy = eps01 * pt->point.crdRel2Centroid[0] + eps11 * pt->point.crdRel2Centroid[1];
#endif
				pt->op_ufs[0].u = ux;
				pt->op_ufs[1].u = uy;
			}
			for (unsigned int i = 0; i < sz_boundaryPoints; ++i)
			{
				if (bcs.allDomain_uFs[i].op_ufs[0].bt != bt_dirichlet)
				{
					cout << "bcs.allDomain_uFs[i]\n" << bcs.allDomain_uFs[i].point.vertexID << '\n';
					THROW("ux, uy for this vertex not provided\n");
				}
			}
			bcss.push_back(bcs);
			lc_stn_sts.timeStepNos.push_back(timeStepNoRead);
			lc_stn_sts.timeValues.push_back(timeValueRead);
		}
		else
			break;
	}
	unsigned int numTimeSteps = bcss.size();

	// nodal forces
	fstream inF(fileName_F.c_str(), ios::in);
	if (!inF.is_open())
	{
		cout << fileName_F << '\n';
		successfulread = false;
		return 0;
		THROW("Cannot open file\n");
	}
	double Fx, Fy;
	unsigned int cntr = 0;
	while (!inF.eof())
	{
		if (Read_LSDYNA_NodalForces(inF, timeValueRead))
		{
			boundary_AllPoints_uFAllDirs_1Case* bcs = &bcss[cntr];
			map<int, int>::iterator it, ite = bcs->vertexID2Pos.end();

			for (unsigned int i = 0; i < sz_boundaryPoints; ++i)
			{
				inF >> buf >> vertID >> buf >> Fx >> buf >> Fy;
				getline(inF, buf);
				it = bcs->vertexID2Pos.find(vertID);
				if (it == ite)
				{
					cout << "vertID\t" << vertID << " this vertex does not belong to boundary set\n";
					THROW("Invalid vertexID\n");
				}
				vertIndex = it->second;
				boundary_1Point_uFAllDirs* pt = &bcs->allDomain_uFs[vertIndex];
				pt->op_ufs[0].bt = bt_neumann;

#if TST_HOMOG
				Fx = sig00 * pt->point.normalTimesWeights[0] + sig01 * pt->point.normalTimesWeights[1];
				Fy = sig01 * pt->point.normalTimesWeights[0] + sig11 * pt->point.normalTimesWeights[1];
#endif
				pt->op_ufs[0].F = Fx;
				pt->op_ufs[1].F = Fy;
			}
			for (unsigned int i = 0; i < sz_boundaryPoints; ++i)
			{
				if (bcs->allDomain_uFs[i].op_ufs[0].bt != bt_neumann)
				{
					cout << "bcs->allDomain_uFs[i]\n" << bcs->allDomain_uFs[i].point.vertexID << '\n';
					THROW("ux, uy for this vertex not provided\n");
				}
			}
			if (!DoublesAreEqual(lc_stn_sts.timeValues[cntr], timeValueRead, 1e-4))
			{
				cout << "timeValues[cntr]\t" << lc_stn_sts.timeValues[cntr] << '\n';
				cout << "timeValueRead\t" << timeValueRead << '\n';
				THROW("time values don't match\n");
			}
			++cntr;

			lc_stn_sts.timeStepNos.push_back(timeStepNoRead);
			lc_stn_sts.timeValues.push_back(timeValueRead);
		}
		else
			break;
	}
	int numActiveSteps = Compute_Average_Strain_Stress_From_BoundaryNode_Disps_and_Forces(bcss, lc_stn_sts.voigtStrains, lc_stn_sts.voigtStresses, onlyIncludeLastStep);
	if (onlyIncludeLastStep)
	{
		double tmp = lc_stn_sts.timeValues[numTimeSteps - 1];
		lc_stn_sts.timeValues.resize(1); lc_stn_sts.timeValues[0] = tmp;
		int tmpi = lc_stn_sts.timeStepNos[numTimeSteps - 1];
		lc_stn_sts.timeStepNos.resize(1); lc_stn_sts.timeStepNos[0] = tmpi;
	}
	return numActiveSteps;
}


void sMesh::ComputePrintBC_AllCases(BCConfig& bcConfig, vector<boundary_AllPoints_uFAllDirs_1Case>& bcs)
{
	ComputePrintBC_AllCases_Aux(bcConfig, fileNamesOut, solverOptionOut, bcs);
}

double sMesh::CalculateQuadArea(const vector<smNode*>& nodePtrs)
{
	double x0 = nodePtrs[0]->crd[0], y0 = nodePtrs[0]->crd[1];
	double x1 = nodePtrs[1]->crd[0], y1 = nodePtrs[1]->crd[1];
	double x2 = nodePtrs[2]->crd[0], y2 = nodePtrs[2]->crd[1];
	double x3 = nodePtrs[3]->crd[0], y3 = nodePtrs[3]->crd[1];
	return CalculateQuadAreaAux(x0, y0, x1, y1, x2, y2, x3, y3);
}

double CalculateQuadAreaAux(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3)
{
	double v01_x = x1 - x0, v01_y = y1 - y0;
	double v02_x = x2 - x0, v02_y = y2 - y0;
	double v03_x = x3 - x0, v03_y = y3 - y0;
	double area_a = fabs(v02_y * v01_x - v02_x * v01_y);
	double area_b = fabs(v02_y * v03_x - v02_x * v03_y);
	return 0.5 * (area_a + area_b);
}

bool Read_LSDYNA_NodalDisplacements(istream& in, unsigned int& timeStepNo, double& timeValue)
{
	string buf;
	unsigned int sz;
	string s2 = "time";
	size_t loc;
	while (!in.eof())
	{
		getline(in, buf);
		loc = buf.find(s2);
		if (loc != std::string::npos)
		{
			vector<string> parts;
			int numF = BreakString(buf, parts);
			for (unsigned int i = 0; i < numF; ++i)
			{
				if (fromString(parts[i], timeStepNo) == true)
				{
					for (unsigned int j = i + 1; j < numF; ++j)
					{
						if (fromString(parts[j], timeValue))
						{
							buf = "";
							while (!in.eof())
							{
								getline(in, buf);
								loc = buf.find("nodal");
								if (loc != std::string::npos)
								{
									loc = buf.find("x-disp");
									if (loc != std::string::npos)
										return true;
									else
									{
										loc = buf.find("x-rot");
										if (loc == std::string::npos)
										{
											THROW("error reading the file\n");
										}
										return Read_LSDYNA_NodalDisplacements(in, timeStepNo, timeValue);
									}
								}
							}
							return false;
						}
					}
					return false;
				}
			}
			return false;
		}
	}
	return false;
}

bool Read_LSDYNA_NodalForces(istream& in, double& timeValue)
{
	string buf;
	unsigned int sz;
	string s2 = "t=";
	size_t loc;
	while (!in.eof())
	{
		getline(in, buf);
		loc = buf.find(s2);
		if (loc != std::string::npos)
		{
			vector<string> parts;
			int numF = BreakString(buf, parts);
			for (unsigned int i = 0; i < numF; ++i)
			{
				if (fromString(parts[i], timeValue) == true)
				{
					if (fromString(parts[i], timeValue))
					{
						buf = "";
						while (!in.eof())
						{
							in >> buf;
//							if (buf == "velocity")
							if (buf == "nodal")
							{
								getline(in, buf);
								return true;
							}
						}
						return false;
					}
					return false;
				}
			}
			return false;
		}
	}
	return false;
}

void FullElasticity_StnStnTensor2VoigtVector(const MATRIX& ten2ndOrder, vector<double>& vec, bool isStrain, double factor)
{
	double factorShear = factor;
	if (!isStrain)
		factorShear *= 0.5;
	vec.resize(VOIGHT_STN_SZ);
	vec[0] = factor * ten2ndOrder[0][0];
#if !DIM1
	vec[1] = factor * ten2ndOrder[1][1];
#if DIM2
	vec[2] = factorShear * (ten2ndOrder[0][1] + ten2ndOrder[1][0]);
#else
	vec[2] = factor * ten2ndOrder[2][2];
	vec[3] = factorShear * (ten2ndOrder[1][2] + ten2ndOrder[2][1]);
	vec[4] = factorShear * (ten2ndOrder[0][2] + ten2ndOrder[2][0]);
	vec[5] = factorShear * (ten2ndOrder[0][1] + ten2ndOrder[1][0]);
#endif
#endif
}

bool DoublesAreEqual(double d1, double d2, double tol)
{
	double den = fabs(d1);
	if (fabs(d2) > den)	den = fabs(d2);
	if (den < 1.0) den = 1.0;
	return (fabs((d1 - d2) / tol) < den);
}

unsigned int BreakString(const string& inString, vector<string>& parts)
{
	std::stringstream ss(inString);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> vstrings(begin, end);
	unsigned int numF = vstrings.size();
	parts.resize(numF);
	for (unsigned int i = 0; i < numF; ++i)
		parts[i] = vstrings[i];
	return numF;
}

void CopyContentFileA_ostream(const string& fileA, ostream& out) 
{
	fstream in(fileA.c_str(), ios::in);
	if (!in.is_open())
	{
		return;
	}
		
	string buf;
	getline(in, buf);
	while (!in.eof())
	{
		out << buf << '\n';
		getline(in, buf);
	}
}
