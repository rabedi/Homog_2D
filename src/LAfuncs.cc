#include "LAfuncs.h"
#include <string>
#include <assert.h>     /* assert */
#include <math.h>       /* fabs */
#include "global.h"

using namespace std;

#if !USE_BLITZ


istream& operator>>(istream& in, VECTOR& dat)
{
	string buf;
	unsigned int sz;
#if USE_DEALII_VECMATIO
	sz = DIM;
#else
	in >> buf >> sz;
#endif
	dat.resize(sz);
	for (unsigned int i = 0; i < sz; ++i)
		in >> dat[i];
	return in;
}

ostream& operator<<(ostream& out, const VECTOR& dat)
{
	unsigned int sz = dat.size();
#if !USE_DEALII_VECMATIO
	out << "size\t" << sz << '\t';
#endif
	for (unsigned int i = 0; i < sz; ++i)
		out << dat[i] << '\t';
	return out;
}


VECTOR::VECTOR(unsigned int sizeIn)
{
	resize(sizeIn);
}

void VECTOR::resize(unsigned int sizeIn)
{
	vec.resize(sizeIn);
}

unsigned int VECTOR::size() const
{
	return vec.size();
}

double& VECTOR::operator[](unsigned int i)
{
	return vec[i];
}

double VECTOR::operator[](unsigned int i) const
{
	return vec[i];
}

VECTOR& VECTOR::operator=(double val)
{
	for (unsigned int i = 0; i < size(); ++i)
		vec[i] = val;
	return *this;
}


//// MATRIX CLASS

MATRIX::MATRIX(unsigned int rowsIn, unsigned int colsIn)
{
	resize(rowsIn, colsIn);
}
void MATRIX::resize(unsigned int rowsIn, unsigned int colsIn)
{
	nrows = rowsIn;
	ncols = colsIn;
	if (nrows < 0)
	{
		cout << "nrows\t" << nrows << '\n';
		THROW("nrows < 0");
	}
	matx.resize(rowsIn);
	for (unsigned int i = 0; i < nrows; ++i)
	{
		matx[i].resize(ncols);
	}
}

unsigned int MATRIX::rows() const
{
	return nrows;
}

unsigned int MATRIX::columns() const
{
	return ncols;
}

vector<double>& MATRIX::operator[](unsigned int i)
{
	return matx[i];
}

const vector<double>& MATRIX::operator[](unsigned int i) const
{
	return matx[i];
}

void MATRIX::ComputeInverse(MATRIX& invMatrix)
{
	unsigned int dim = nrows;
	invMatrix.resize(nrows, nrows);
	if (dim == 1)
	{
		invMatrix[0][0] = 1.0 / matx[0][0];
		return;
	}
	if (dim == 2)
	{
		double det = matx[0][0] * matx[1][1] - matx[0][1] * matx[1][0];
		double detInv = 1.0 / det;
		invMatrix[0][0] =  detInv * matx[1][1];
		invMatrix[1][1] =  detInv * matx[0][0];
		invMatrix[0][1] = -detInv * matx[0][1];
		invMatrix[1][0] = -detInv * matx[1][0];
		return;
	}
	if (dim == 3)
	{
		double Ap00 =   matx[1][1] * matx[2][2] - matx[2][1] * matx[1][2];
		double Ap01 = -(matx[1][0] * matx[2][2] - matx[2][0] * matx[1][2]);
		double Ap02 =   matx[1][0] * matx[2][1] - matx[2][0] * matx[1][1];

		double Ap10 = -(matx[0][1] * matx[2][2] - matx[2][1] * matx[0][2]);
		double Ap11 =   matx[0][0] * matx[2][2] - matx[2][0] * matx[0][2];
		double Ap12 = -(matx[0][0] * matx[2][1] - matx[2][0] * matx[0][1]);

		double Ap20 =   matx[0][1] * matx[1][2] - matx[1][1] * matx[0][2];
		double Ap21 = -(matx[0][0] * matx[1][2] - matx[1][0] * matx[0][2]);
		double Ap22 =   matx[0][0] * matx[1][1] - matx[1][0] * matx[0][1];
		double det = matx[0][0] * Ap00 + matx[0][1] * Ap01 + matx[0][2] * Ap02;
		double tmp, tol = 0.0;
		for (int i = 0; i < nrows; ++i)
			for (int j = 0; j < nrows; ++j)
			{
				tmp = fabs(matx[i][j]);
				if (tol < tmp)
					tol = tmp;
			}
		tol *= 1e-16;
		if (fabs(det) < tol)
		{
			cout << "det\t" << det;
			THROW("zero determinant matrix\n");
		}
		double detInv = 1.0 / det;
		invMatrix[0][0] = detInv * Ap00;
		invMatrix[0][1] = detInv * Ap10;
		invMatrix[0][2] = detInv * Ap20;

		invMatrix[1][0] = detInv * Ap01;
		invMatrix[1][1] = detInv * Ap11;
		invMatrix[1][2] = detInv * Ap21;

		invMatrix[2][0] = detInv * Ap02;
		invMatrix[2][1] = detInv * Ap12;
		invMatrix[2][2] = detInv * Ap22;
	}
}

#if 0
double& MATRIX::operator()(unsigned int i, unsigned int j)
{
	// example of using macro DB defined before for debugging
#if DEBUG_MODE
 if ((i >= rows()) || (j >> columns()))
   THROW("index output of bound\n");
#endif
	return matx[i][j];
}

double MATRIX::operator()(unsigned int i, unsigned int j) const
{
	return matx[i][j];
}
#endif

MATRIX& MATRIX::operator=(double val)
{
	for (unsigned int i = 0; i < nrows; ++i)
	{
		for (unsigned int j = 0; j < ncols; ++j)		
			matx[i][j] = val;
	}
	return *this;
}

istream& operator>>(istream& in, MATRIX& dat)
{
	string buf;
	unsigned int nrow, ncol;
#if USE_DEALII_VECMATIO
	nrow = DIM, ncol = DIM;
#else
	in >> buf >> nrow >> buf >> ncol;
#endif
	dat.resize(nrow, ncol);
	for (unsigned int i = 0; i < nrow; ++i)
	{
		for (unsigned int j = 0; j < ncol; ++j)
			in >> dat[i][j];
	}
	return in;
}

ostream& operator<<(ostream& out, const MATRIX& dat)
{
	unsigned int nrow = dat.rows(), ncol = dat.columns();
#if !USE_DEALII_VECMATIO
	out << "rows\t" << nrow << "\tcols\t" << ncol << '\n';
#endif
	for (unsigned int i = 0; i < nrow; ++i)
	{
		for (unsigned int j = 0; j < ncol; ++j)
			out << dat[i][j] << '\t';
#if USE_DEALII_VECMATIO
		out << '\t';
#else
		out << '\n';
#endif
	}
	return out;
}



// 	This is the LU decomposition and solve implementation that we used 
// 	with TNT.  I's put here so people don't need to use CLAPACK to run 
// 	the code.
int LUsolve(MATRIX &K, VECTOR &F)
{
    int m, n, pivsign;
    m = K.rows();
    n = K.columns();
    MATRIX LU_(m, n);
	LU_ = K;	// CHECK: can we just get along with referencing? It seems to be the case, then we save some mempory copy here
//    MATRIX LU_(K); // refernce option

    int i=0;
    int j=0;
    int k=0;


    pivsign = 1;
    VECTOR LUcolj(m);
    vector<int> piv(m);

    for (i = 0; i < m; i++)
        piv[i] = i;

    for (j = 0; j < n; j++) 
    {
        // Make a copy of the j-th column to localize references.
        for (i = 0; i < m; i++)
            LUcolj[i] = LU_[i][j];

        // Apply previous transformations.
        for (i = 0; i < m; i++) 
        {
            int kmax = i<j?i:j; //min(i,j);
            double s = 0.0;
            for (k = 0; k < kmax; k++)
                s += LU_[i][k]*LUcolj[k];

            LU_[i][j] = LUcolj[i] -= s;
        }

        // Find pivot and exchange if necessary.
        int p = j;
        for (int i = j+1; i < m; i++)
            if (fabs(LUcolj[i] ) > fabs(LUcolj[p] ) )
                p = i;

        if (p != j) 
        {
            for (k = 0; k < n; k++) 
            {
                double t = LU_[p][k]; 
                LU_[p][k] = LU_[j][k]; 
                LU_[j][k] = t;
            }
            k = piv[p]; 
			piv[p] = piv[j]; 
            piv[j] = k;
            pivsign = -pivsign;
        }

        // Compute multipliers.
        if ((j < m) && (LU_[j][j] != 0.0))
            for (int i = j+1; i < m; i++) //... [thite, 6/30/2003]
                LU_[i][j] /= LU_[j][j];
    }

    int issingular=0;
    for(int j=0; j<n; j++)
    { 
        if(LU_[j][j] == 0)
            issingular = 1;
    }
#if 0
    if( issingular != 0 )
    {
        std::cout<<"Solve procedure failed:  FILE:  "<<__FILE__

			<<", LINE:  "<<__LINE__<<std::endl;
	    return issingular;
//        exit(1);
    }
#endif
    int piv_length = piv.size();
    VECTOR x(piv_length);
    for (int I = 0; I < piv_length; I++) 
        x[I] = F[piv[I]];

    // Solve L*Y = B(piv)
    for (k = 0; k < n; k++) 
        for (i = k+1; i < n; i++) 
            x[i] -= x[k]*LU_[i][k];

    // Solve U*X = Y;
    for (k = n-1; k >= 0; k--) 
    {
        x[k] /= LU_[k][k];
        for (i = 0; i < k; i++) 
            x[i] -= x[k]*LU_[i][k];
    }

    F = x; 
    K = LU_;
    return issingular;
} // end TNT version of LUsolve with 1D F

#endif
