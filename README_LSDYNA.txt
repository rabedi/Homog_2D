To run the code:

1. Open Target.sln (or open your own version of the project with h and cpp files)

2. In main.cpp line 16 change the homogType to mbc_Disp, mbc_Trac, mbc_MixedMineDispNormal for displacement, traction, and mixed BCs, respectively. Right now, the state is for Displacement BC which is fine.

3. Run the code

4. Micro_out_0.k, Micro_out_1.k, Micro_out_2.k are the generated files for load cases 0 to 2 (pull in x, pull in y, shear in xy)

5. Try running Micro_out_0.k

---------------------------------------------------

In order to fix any remaining format errors:

A. Header, ... stuff
In the code folder, you'll find the two files:
LS_DYNA_IMPLICIT/_LS_DYNA_header.txt
LS_DYNA_IMPLICIT/_LS_DYNA_disp_load_curves.txt

You can directly open these two files, make correction, and redo steps 1 to 4 above. Make sure you leave one empty line at the end of the files (as is now).

B. Code generated parts:
If needed, edit functions
SpaceMesh.cpp

void sMesh::PrintNode_Element_ElementSet_LSDYNA(ostream& out)
void boundary_AllPoints_uFAllDirs_1Case::PrintBC(FEMSolverT solverOption, ostream & out, string inputFolder4Run) - line 266 to the end of functions

these functions are called from 
void sMesh::ComputePrintBC_AllCases_Aux(BCConfig& bcConfig, const vector<string>& fileNames, FEMSolverT solverOptionIn, vector<boundary_AllPoints_uFAllDirs_1Case>& bcs)
Again, if you make any changes here, clearly steps 1 to 4 should be repeated.
------------
