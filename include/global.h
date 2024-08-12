#ifndef FSI_GLOBAL__H
#define FSI_GLOBAL__H

#define USE_DEALII	0
#define USE_DEALII_VECMATIO 1

#include <iostream> /* for input and output using >> and << operator */
#include <fstream>  /* file streams for input and output */
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>
#include <assert.h>     /* assert */
#include <math.h>     /* math functions */
#include <array>
#include <map>

#if USE_DEALII
#include <deal.II/base/tensor.h>
using namespace dealii ;
// uses deal.ii I/O format for vectors and matrices
// include the correct deal.ii files
typedef Tensor<1, DIM> VEC;
typedef Tensor<2, DIM> MAT;
#else
#include "LAfuncs.h"
typedef VECTOR VEC;
typedef MATRIX MAT;

typedef double Vec[2];
#endif

#define DIM1	0
#define DIM2	1
#define DIM3	0

#define ANTI_PLANE_SHEAR 0

#if ANTI_PLANE_SHEAR
#define DIM_U 1
#else
#define DIM_U DIM
#endif


#if DIM2
#define DIM	2
	#if ANTI_PLANE_SHEAR
	#define VOIGHT_STN_SZ 2
	#else
	#define VOIGHT_STN_SZ 3
	#endif
#else
#if DIM3
#define DIM 3
#else
#define DIM 1
#endif
#endif


#include <iostream> /* for input and output using >> and << operator */
#include <fstream>  /* file streams for input and output */
#include <assert.h>     /* assert */
#include <math.h>     /* math functions */
#include <array>

using namespace std;

#define MAX(a,b) ((a) > (b) ? a : b )
#define MIN(a,b) ((a) > (b) ? b : a )

#ifndef EXIT
#define EXIT exit(1);getchar();getchar();
#endif

// to exit the code with useful information when something is wrong. It provides file and line number with a message if a check is incorrect.
// example
//					if (denominator == 0) THROW("cannot divide be zero!\n");
#ifndef THROW
#define THROW(msg){{char tmpstr[255];sprintf(tmpstr, "In %s, line %d : %s \n", __FILE__, __LINE__, msg);	cerr<<tmpstr;getchar(); getchar(); throw(tmpstr); }}
#endif

#define P_INFINITY 0.1 * DBL_MAX
#define M_INFINITY -0.1 * DBL_MAX

typedef enum {mnx_none_t = -1, mn_t, mx_t, mnmxT_SIZE} mnmxT;
extern string g_rootFolder_4Output;
extern string main_dir;
extern string inc_dir;

#endif
