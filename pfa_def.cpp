#include "pfa_def.h"

#include <memory>






cmat cnew(int r, int c)
{
	cmat mat;
	
	mat.row = r;
	mat.col = c;

	mat.nmat.real = new double[r * c];
	mat.nmat.imag = new double[r * c];

	memset(mat.nmat.real, 0, r*c * sizeof(double));
	memset(mat.nmat.imag, 0, r*c * sizeof(double));

	return mat;
}



cmat cnew(double * real, double * imag, int r, int c)
{
	cmat mat;

	mat.row = r;
	mat.col = c;

	mat.nmat.real = new double[mat.row * mat.col];
	mat.nmat.imag = new double[mat.row * mat.col];

	memcpy(mat.nmat.real, real, mat.row * mat.col * sizeof(double));
	memcpy(mat.nmat.imag, imag, mat.row * mat.col * sizeof(double));

	return mat;
}






void cdel(cmat a)
{
	delete[] a.nmat.real;
	delete[] a.nmat.imag;
}







void cshow(cmat a)
{
	for (int i = 0; i < a.row; i++)
	{
		for (int j = 0; j < a.col; j++)
		{
			double r = a.nmat.real[i*a.col + j];
			double c = a.nmat.imag[i*a.col + j];
			
			
			
			if(c >= 0)
			    printf("%.4lf+j%.4lf ", r, c);
			else
				printf("%.4lf-j%.4lf ", r, -c);
		}
		printf("\n");
	}

}







cplx cget(cmat a, int r, int c)
{
	cplx v;

	v.real = a.nmat.real[r*a.col + c];
	v.imag = a.nmat.imag[r*a.col + c];

	return v;
}





void cset(cmat a, cplx v, int r, int c)
{
	
	a.nmat.real[r*a.col + c] = v.real;
	a.nmat.imag[r*a.col + c] = v.imag;
}







cmat ccpy(cmat a)
{
	cmat mat;

	mat.row = a.row;
	mat.col = a.col;

	mat.nmat.real = a.nmat.real;
	mat.nmat.imag = a.nmat.imag;

	return mat;
}










cmat cdup(cmat a)
{
	cmat mat;

	mat.row = a.row;
	mat.col = a.col;

	mat.nmat.real = new double[mat.row * mat.col];
	mat.nmat.imag = new double[mat.row * mat.col];

	memcpy(mat.nmat.real, a.nmat.real, mat.row * mat.col * sizeof(double));
	memcpy(mat.nmat.imag, a.nmat.imag, mat.row * mat.col * sizeof(double));

	return mat;
}










cmat creal(cmat a)
{
	cmat mat;

	mat.row = a.row;
	mat.col = a.col;

	mat.nmat.real = new double[mat.row * mat.col];
	mat.nmat.imag = new double[mat.row * mat.col];

	memcpy(mat.nmat.real, a.nmat.real, mat.row * mat.col * sizeof(double));
	memset(mat.nmat.imag, 0, mat.row * mat.col * sizeof(double));

	return mat;
}

cmat cimag(cmat a)
{
	cmat mat;

	mat.row = a.row;
	mat.col = a.col;

	mat.nmat.real = new double[mat.row * mat.col];
	mat.nmat.imag = new double[mat.row * mat.col];

	memset(mat.nmat.real, 0, mat.row * mat.col * sizeof(double));
	memcpy(mat.nmat.imag, a.nmat.imag, mat.row * mat.col * sizeof(double));

	return mat;
}

cmat cconj(cmat a)
{
	cmat mat;

	mat = cdup(a);

	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < mat.col; j++)
		{
			mat.nmat.imag[i*mat.col + j] = -1 * mat.nmat.imag[i*mat.col + j];
		}
	}

	return mat;
}


cmat ctrans(cmat a)
{
	cmat mat;

	mat = cdup(a);

	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < i; j++)
		{
			double real = mat.nmat.real[i*mat.col + j];
			mat.nmat.real[i*mat.col + j] = mat.nmat.real[j*mat.col + i];
			mat.nmat.real[j*mat.col + i] = real;

			double imag = mat.nmat.imag[i*mat.col + j];
			mat.nmat.imag[i*mat.col + j] = mat.nmat.imag[j*mat.col + i];
			mat.nmat.imag[j*mat.col + i] = imag;
		}
	}

	return mat;
}

cmat cadd(cmat a, cmat b)
{
	cmat mat;
	int row;
	int col;

	row = a.row;
	col = a.col;

	mat = cnew(row, col);

	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < mat.col; j++)
		{
			double real = a.nmat.real[i*mat.col + j] + b.nmat.real[i*mat.col + j];
			mat.nmat.real[i*mat.col + j] = real;

			double imag = a.nmat.imag[i*mat.col + j] + b.nmat.imag[i*mat.col + j];
			mat.nmat.imag[i*mat.col + j] = imag;
		}
	}

	return mat;
}

cmat csub(cmat a, cmat b)
{
	cmat mat;
	int row;
	int col;

	row = a.row;
	col = a.col;

	mat = cnew(row, col);

	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < mat.col; j++)
		{
			double real = a.nmat.real[i*mat.col + j] - b.nmat.real[i*mat.col + j];
			mat.nmat.real[i*mat.col + j] = real;

			double imag = a.nmat.imag[i*mat.col + j] - b.nmat.imag[i*mat.col + j];
			mat.nmat.imag[i*mat.col + j] = imag;
		}
	}

	return mat;
}

cmat cmul(cmat a, cmat b)
{
	cmat mat;
	int row;
	int col;

	row = a.row;
	col = b.col;

	mat = cnew(row, col);

	for (int i = 0; i < mat.row; i++)
	{
		for (int j = 0; j < mat.col; j++)
		{
			for (int k = 0; k < a.col; k++)
			{
				double ar = a.nmat.real[i*a.col + k];
				double ai = a.nmat.imag[i*a.col + k];
				double br = b.nmat.real[k*b.col + j];
				double bi = b.nmat.imag[k*b.col + j];

				mat.nmat.real[i*mat.col + j] += ar * br - ai * bi;
				mat.nmat.imag[i*mat.col + j] += ar * bi + ai * br;
			}
		}
	}

	return mat;
}




cmat cdiv(cmat a, cmat b)
{
	cmat mat;

	if (a.row != b.row || b.row != b.col)
		return mat;

	cmat inv = cinv(b);

	mat = cmul(inv, a);

	cdel(inv);

	return mat;
}







void cadds(cmat a, cplx v, int r, int c)
{
	cplx ve = cget(a, r, c);
	ve.real = ve.real + v.real;
	ve.imag = ve.imag + v.imag;
	cset(a, ve, r, c);
}

void csubs(cmat a, cplx v, int r, int c)
{
	cplx ve = cget(a, r, c);
	ve.real = ve.real - v.real;
	ve.imag = ve.imag - v.imag;
	cset(a, ve, r, c);
}

void cmuls(cmat a, cplx v, int r, int c)
{
	cplx vr;
	cplx ve = cget(a, r, c);
	vr.real = ve.real * v.real - ve.imag * v.imag;
	vr.imag = ve.real * v.imag + ve.imag * v.real;
	cset(a, vr, r, c);
}

void cdivs(cmat a, cplx v, int r, int c)
{
	cplx vr;
	cplx ve = cget(a, r, c);
    double den = v.real * v.real + v.imag * v.imag;
	vr.real = (ve.real * v.real + ve.imag * v.imag) / den;
	vr.imag = (0 - ve.real * v.imag + ve.imag * v.real) / den;
	cset(a, vr, r, c);
}









cmat cmnr(cmat a, int r, int c)
{
	cmat mat;
	cplx v;
	int idx = 0;
	int n = a.row;

	if (a.row != a.col)
		return mat;

	mat = cnew(n - 1, n - 1);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i != r && j != c)
			{
				v = cget(a, i, j);
				mat.nmat.real[idx] = v.real;
				mat.nmat.imag[idx] = v.imag;
				idx++;
			}
		}
	}

	return mat;
}


double cdet(cmat a)
{
	double det = 0.0;
	double si = 1.0;
	cmat mnr;
	int n = a.row;


	if (a.row != a.col)
		return det;

	if (n == 1)
		return a.nmat.real[0];

	for (int k = 0; k < n; k++)
	{
		mnr = cmnr(a, 0, k);

		si = (k % 2) ? -1 : 1;
		det += si * cget(a, 0, k).real * cdet(mnr);
		cdel(mnr);
	}

	return det;
}


cmat cadjg(cmat a)
{
	double det = 0.0;
	double si = 1.0;
	int n = a.row;
	cmat adjg;
	cmat mnr;

	if (a.row != a.col)
		return adjg;

	if (n == 1)
		return a;

	adjg = cnew(n, n);

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			mnr = cmnr(a, i, j);

			si = ((i + j) % 2) ? -1 : 1;

			cplx v;
			v.real = si * cdet(mnr);
			v.imag = 0.0;
			cset(adjg, v, j, i);

			cdel(mnr);
		}
	}

	return adjg;
}

cmat cinv(cmat a)
{
	cmat mat;
	double det = 0.0;

	int n = a.row;

	if (a.row != a.col)
		return mat;

	mat = cadjg(a);
	det = cdet(a);

    cplx v;
	v.real = det;
	v.imag = 0.0;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cdivs(mat, v, i, j);
		}
	}

	return mat;
}











