#ifndef FSI_GLOBAL__H
#define FSI_GLOBAL__H

#include "Dims.h"
#include "commonMacros.h"

// See the private class in SL_interfacePPtData_Time_Seq, the size - timeSeqPtrs
#define SIZE_PT_TIME_SEQUENCE	202

// Ring problem is 1D ring fragmentation problem. It has a source terms from v_r
#define RING_PROBLEM	0

#define USE_ISO_ASSUMPTION_PRE	0	// if iso assumption is used, normal and shear modes (2D, 3D) are decoupled and C, Z, Y, ..., in matrix form are not used
#define USE_PLANE_STRAIN	1	// for 2D isotropic model uses plane strain
// turn the flag below one, if the problem has source term that is 0th order in q = [v, sigma]
// For example, when damping is nonzero (D_vv = damping / rho) the flag below is 1
// see D matrix components in SL_Bulk_Properties
#define HAVE_SOURCE_ORDER0_q	0 //1
#define HAVE_SOURCE_TERM_OTHER	0	// source term other than 0th order q and ring problem source (e.g. constant body force ...)


#if HAVE_SOURCE_ORDER0_q
#define HAVE_SOURCE	1
#else
	#if RING_PROBLEM
	#define HAVE_SOURCE	1
	#else
		#if HAVE_SOURCE_TERM_OTHER
		#define HAVE_SOURCE	1
		#else
		#define HAVE_SOURCE	0
		#endif
	#endif
#endif

#if DiM1
#define USE_ISO_ASSUMPTION	1
#else
#define USE_ISO_ASSUMPTION	USE_ISO_ASSUMPTION_PRE
#endif

// C, Z, Y, etc. are only necessary for anisotropic C, but still can be computed for isotropic C
// for debugging, we can still compute them (second line), but first line should be the default
//#define COMPUTE_ANISO_BULK	!USE_ISO_ASSUMPTION
#define COMPUTE_ANISO_BULK	!USE_ISO_ASSUMPTION

#if USE_DEALII
#include <deal.II/base/tensor.h>
using namespace dealii ;
// uses deal.ii I/O format for vectors and matrices
// include the correct deal.ii files
typedef Tensor<1, DIM> VEC;
typedef Tensor<2, DIM> MAT;

#if DiM2a3_F
typedef Tensor<1, DIMm1> VECm1;
typedef Tensor<2, DIMm1> MATm1;
#endif

#else
//#include "LAfuncs.h"
//typedef VECTOR VEC;
//typedef MATRIX MAT;

#include "LAfuncsFixed.h"
typedef Vc_dd VEC;
typedef Mtrx_dd MAT;

#if 1 //DiM2a3_F
typedef Vc_dm1 VECm1;
typedef Mtrx_dm1 MATm1;
#endif

typedef long GID;
#endif
#endif

extern string g_prefileName;

// indices for printing for directions
extern vector<string> pINDs;
// this one has a comma too (if needed)
extern vector<string> pINDcommas;
void setGlobalMembers();
