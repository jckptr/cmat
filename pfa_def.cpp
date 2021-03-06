#include "pfa_def.h"

#include <memory>
#include <cmath>
#include <random>
#include <ctime>





static struct _csys_mem_ _m;




cmat cnew(int r, int c)
{
    cmat mat;
    int i;

    mat.row = r;
    mat.col = c;

    if (r*c == 0)
        return mat;

    if (_m.nnmat >= NMAX_NORMALMAT)
    {
        printf("cnew: Memory allocate failed\n");
        return mat;
    }

    for (i = 0; i < NMAX_NORMALMAT; i++)
    {
        if (_m._nmat_pos_[i] == 0)
        {
            _m._nmat_pos_[i] = 1;
            _m.nnmat++;
            mat.nmat = &(_m._nmat_[i]);
            break;
        }
    }

    if (i >= NMAX_NORMALMAT)
    {
        printf("cnew: Memory allocate failed\n");
        return mat;
    }


    mat.nmat -> real = new double[r * c];
    mat.nmat -> imag = new double[r * c];

    memset(mat.nmat -> real, 0, r * c * sizeof(double));
    memset(mat.nmat -> imag, 0, r * c * sizeof(double));

    //printf("cnew: new cmat allocated successfully\n");

    return mat;
}






void cdel(cmat a)
{
    if (_m.nnmat <= 0)
        return;

    for (int i = 0; i < NMAX_NORMALMAT; i++)
    {
        if (_m._nmat_pos_[i] == 1)
        {
            if (&(_m._nmat_[i]) == a.nmat)
            {
                _m._nmat_pos_[i] = 0;
                _m.nnmat--;

                delete[] _m._nmat_[i].real;
                delete[] _m._nmat_[i].imag;

                //printf("cdel: a cmat deleted successfully\n");

                break;
            }
        }
    } 
}




cmat cdup(cmat a)
{
    cmat mat;

    mat = cnew(a.row, a.col);

    memcpy(mat.nmat -> real, a.nmat -> real, mat.row * mat.col * sizeof(double));
    memcpy(mat.nmat -> imag, a.nmat -> imag, mat.row * mat.col * sizeof(double));

    return mat;
}



void cclr()
{
    if (_m.nnmat <= 0)
        return;

    //printf("remain: %d\n", _m.nnmat);

    for (int i = 0; i < NMAX_NORMALMAT; i++)
    {
        if (_m._nmat_pos_[i] == 1)
        {
            _m._nmat_pos_[i] = 0;
            _m.nnmat--;

            delete[] _m._nmat_[i].real;
            delete[] _m._nmat_[i].imag;

            //printf("cclr: a cmat deleted successfully\n");

            if (_m.nnmat <= 0)
                break;
        }
    }
}




int cisempt(cmat a)
{
    return ((a.col * a.row) ? 1 : 0);
}


int ciszero(cplx v)
{
    return (abs(v.real) < CDBL_EPSILON && abs(v.imag) < CDBL_EPSILON);
}

int cisreal(cplx v)
{
    return (abs(v.real) > CDBL_EPSILON && abs(v.imag) < CDBL_EPSILON);
}

int cisimag(cplx v)
{
    return (abs(v.real) < CDBL_EPSILON && abs(v.imag) > CDBL_EPSILON);
}

int cisint(cplx v)
{
    cplx vr;
    double null;

    vr = cround(v);
    vr = csub(v, vr);

    return ciszero(vr);
}

int ciszero(cmat a)
{
    cplx v;

    int row = a.row;
    int col = a.col;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v = cget(a, i, j);

            if (!ciszero(v))
            {
                return 0;
            }
        }
    }

    return 1;
}

int cisreal(cmat a)
{
    cplx v;
    int allzero = 1;

    int row = a.row;
    int col = a.col;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v = cget(a, i, j);

            if (cisimag(v))
            {
                return 0;
            }
            else if (cisreal(v))
            {
                allzero = 0;
            }
        }
    }

    return (!allzero);
}

int cisimag(cmat a)
{
    cplx v;
    int allzero = 1;

    int row = a.row;
    int col = a.col;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v = cget(a, i, j);

            if (cisreal(v))
            {
                return 0;
            }
            else if (cisimag(v))
            {
                allzero = 0;
            }
        }
    }

    return (!allzero);
}

int cisint(cmat a)
{
    cplx v;

    int row = a.row;
    int col = a.col;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v = cget(a, i, j);

            if (!cisint(v))
            {
                return 0;
            }
        }
    }

    return 1;
}






void cshow(cplx v)
{
    if (ciszero(v))
    {
        printf("0\n");
    }
    else if (cisreal(v))
    {
        if (cisint(v))
            printf("%.0lf\n", v.real);
        else
            printf("%.4lf\n", v.real);
    }
    else if (cisimag(v))
    {
        if (cisint(v))
            printf("%.0lf\n", v.imag);
        else
        {
            if (v.imag >= 0)
                printf("j%.4lf\n", v.imag);
            else
                printf("-j%.4lf\n", -v.imag);
        }
    }
    else
    {
        if (v.imag >= 0)
            printf("%.4lf+j%.4lf\n", v.real, v.imag);
        else
            printf("%.4lf-j%.4lf\n", v.real, -v.imag);
    }
}


void cshow(cmat a)
{
    cplx v;

    printf("[");
    for (int i = 0; i < a.row; i++)
    {
        if (i != 0)
            printf(" ");

        printf("[");
        for (int j = 0; j < a.col; j++)
        {
            v.real = a.nmat -> real[i*a.col + j];
            v.imag = a.nmat -> imag[i*a.col + j];

            if (ciszero(v))
            {
                printf("0,");
            }
            else if (cisreal(v))
            {
                if (cisint(v))
                    printf("%.0lf,", v.real);
                else
                    printf("%.4lf,", v.real);
            }
            else if (cisimag(v))
            {
                if (cisint(v))
                    printf("j%.0lf,", v.imag);
                else
                {
                    if (v.imag >= 0)
                        printf("j%.4lf,", v.imag);
                    else
                        printf("-j%.4lf,", -v.imag);
                }
            }
            else
            {
                if (v.imag >= 0)
                    printf("%.4lf+j%.4lf,", v.real, v.imag);
                else
                    printf("%.4lf-j%.4lf,", v.real, -v.imag);
            }
        }
        if (i >= a.row - 1)
            printf("\b]]\n");
        else
            printf("\b],\n");
    }
}







cplx cget(cmat a, int r, int c)
{
    cplx v;

    v.real = a.nmat -> real[r*a.col + c];
    v.imag = a.nmat -> imag[r*a.col + c];

    return v;
}





void cset(cmat a, cplx v, int r, int c)
{
    
    a.nmat -> real[r*a.col + c] = v.real;
    a.nmat -> imag[r*a.col + c] = v.imag;
}




int ceql(cplx a, cplx b)
{
    return (abs(a.real - b.real) < DBL_EPSILON && abs(a.imag - b.imag) < DBL_EPSILON);
}

int ceql(cmat a, cmat b)
{
    if (a.row != b.row || a.col != b.col)
        return 0;

    int row = a.row;
    int col = a.col;

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (!ceql(cget(a, i, j), cget(b, i, j)))
                return 0;
        }
    }

    return 1;
}







cplx cceil(cplx v)
{
    cplx vr;

    vr.real = ceil(v.real);
    vr.imag = ceil(v.imag);

    return vr;
}

cplx cfloor(cplx v)
{
    cplx vr;

    vr.real = floor(v.real);
    vr.imag = floor(v.imag);

    return vr;
}

cplx cround(cplx v)
{
    cplx vr;

    vr.real = round(v.real);
    vr.imag = round(v.imag);

    return vr;
}


cmat cceil(cmat a)
{
    cmat mat;
    cplx v;

    int row = a.row;
    int col = a.col;

    mat = cnew(row, col);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v = cget(a, i, j);
            v = cceil(v);
            cset(mat, v, i, j);
        }
    }

    return mat;
}

cmat cfloor(cmat a)
{
    cmat mat;
    cplx v;

    int row = a.row;
    int col = a.col;

    mat = cnew(row, col);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v = cget(a, i, j);
            v = cfloor(v);
            cset(mat, v, i, j);
        }
    }

    return mat;
}

cmat cround(cmat a)
{
    cmat mat;
    cplx v;

    int row = a.row;
    int col = a.col;

    mat = cnew(row, col);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v = cget(a, i, j);
            v = cround(v);
            cset(mat, v, i, j);
        }
    }

    return mat;
}




cmat creal(cmat a)
{
    cmat mat;

    mat = cnew(a.row, a.col);

    memcpy(mat.nmat -> real, a.nmat -> real, mat.row * mat.col * sizeof(double));

    return mat;
}



cmat cimag(cmat a)
{
    cmat mat;

    mat = cnew(a.row, a.col);

    memcpy(mat.nmat -> imag, a.nmat -> imag, mat.row * mat.col * sizeof(double));

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
            mat.nmat -> imag[i*mat.col + j] = -1 * mat.nmat -> imag[i*mat.col + j];
        }
    }

    return mat;
}






cmat ctrans(cmat a)
{
    cmat mat;
    cplx v;

    int row = a.row;
    int col = a.col;

    mat = cnew(col, row);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            v = cget(a, i, j);
            cset(mat, v, j, i);
        }
    }

    return mat;
}




cplx ctrace(cmat a)
{
    cplx v;
    cplx temp;
    int n = a.row;

    if (a.row != a.col)
        return v;

    for (int i = 0; i < n; i++)
    {
        temp = cget(a, i, i);
        v = cadd(v, temp);
    }

    return v;
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
            double real = a.nmat -> real[i*mat.col + j] + b.nmat -> real[i*mat.col + j];
            mat.nmat -> real[i*mat.col + j] = real;

            double imag = a.nmat -> imag[i*mat.col + j] + b.nmat -> imag[i*mat.col + j];
            mat.nmat -> imag[i*mat.col + j] = imag;
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
            double real = a.nmat -> real[i*mat.col + j] - b.nmat -> real[i*mat.col + j];
            mat.nmat -> real[i*mat.col + j] = real;

            double imag = a.nmat -> imag[i*mat.col + j] - b.nmat -> imag[i*mat.col + j];
            mat.nmat -> imag[i*mat.col + j] = imag;
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
                double ar = a.nmat -> real[i*a.col + k];
                double ai = a.nmat -> imag[i*a.col + k];
                double br = b.nmat -> real[k*b.col + j];
                double bi = b.nmat -> imag[k*b.col + j];

                mat.nmat -> real[i*mat.col + j] += ar * br - ai * bi;
                mat.nmat -> imag[i*mat.col + j] += ar * bi + ai * br;
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



cmat cadds(cmat a, cplx v)
{
    cmat mat;
    cplx ve;
    int row;
    int col;

    row = a.row;
    col = a.col;

    mat = cnew(row, col);

    for (int i = 0; i < mat.row; i++)
    {
        for (int j = 0; j < mat.col; j++)
        {
            ve = cget(a, i, j);
            cset(mat, cadd(ve, v), i, j);
        }
    }

    return mat;
}

cmat csubs(cmat a, cplx v)
{
    cmat mat;
    cplx ve;
    int row;
    int col;

    row = a.row;
    col = a.col;

    mat = cnew(row, col);

    for (int i = 0; i < mat.row; i++)
    {
        for (int j = 0; j < mat.col; j++)
        {
            ve = cget(a, i, j);
            cset(mat, csub(ve, v), i, j);
        }
    }

    return mat;
}

cmat cmuls(cmat a, cplx v)
{
    cmat mat;
    cplx ve;
    int row;
    int col;

    row = a.row;
    col = a.col;

    mat = cnew(row, col);

    for (int i = 0; i < mat.row; i++)
    {
        for (int j = 0; j < mat.col; j++)
        {
            ve = cget(a, i, j);
            cset(mat, cmul(ve, v), i, j);
        }
    }

    return mat;
}


cmat cdivs(cmat a, cplx v)
{
    cmat mat;
    cplx ve;
    int row;
    int col;

    row = a.row;
    col = a.col;

    mat = cnew(row, col);

    for (int i = 0; i < mat.row; i++)
    {
        for (int j = 0; j < mat.col; j++)
        {
            ve = cget(a, i, j);
            cset(mat, cdiv(ve, v), i, j);
        }
    }

    return mat;
}



cmat cmuls(cmat a, cmat b)
{
    cmat mat;
    int row;
    int col;

    row = a.row;
    col = a.col;

    if (a.row != b.row || a.col != b.col)
        return mat;

    mat = cnew(row, col);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            cset(mat, cmul(cget(a, i, j), cget(b, i, j)), i, j);
        }
    }

    return mat;
}


cmat cdivs(cmat a, cmat b)
{
    cmat mat;
    int row;
    int col;

    row = a.row;
    col = a.col;

    if (a.row != b.row || a.col != b.col)
        return mat;

    mat = cnew(row, col);

    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            cset(mat, cdiv(cget(a, i, j), cget(b, i, j)), i, j);
        }
    }

    return mat;
}



cmat cone(int n)
{
    cmat mat;
    cplx v;

    mat = cnew(n, n);
    v.real = 1.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cset(mat, v, i, j);
        }
    }

    return mat;
}

cmat cone(int r, int c)
{
    cmat mat;
    cplx v;

    mat = cnew(r, c);
    v.real = 1.0;

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            cset(mat, v, i, j);
        }
    }

    return mat;
}

cmat ceye(int n)
{
    cmat mat;
    cplx v;

    mat = cnew(n, n);
    v.real = 1.0;

    for (int i = 0; i < n; i++)
    {
        cset(mat, v, i, i);
    }

    return mat;
}

cplx crand()
{
    return cplx();
}

cmat crand(int n)
{
    return cmat();
}

cmat crand(int r, int c)
{
    return cmat();
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
                mat.nmat -> real[idx] = v.real;
                mat.nmat -> imag[idx] = v.imag;
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
        return a.nmat -> real[0];

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








cplx cadd(cplx a, cplx b)
{
    cplx v;
    
    v.real = a.real + b.real;
    v.imag = a.imag + b.imag;

    return v;
}

cplx csub(cplx a, cplx b)
{
    cplx v;

    v.real = a.real - b.real;
    v.imag = a.imag - b.imag;

    return v;
}

cplx cmul(cplx a, cplx b)
{
    cplx v;

    v.real = a.real * b.real - a.imag * b.imag;
    v.imag = a.real * b.imag + a.imag * b.real;

    return v;
}

cplx cdiv(cplx a, cplx b)
{
    cplx v;

    double den = b.real * b.real + b.imag * b.imag;
    v.real = (a.real * b.real + a.imag * b.imag) / den;
    v.imag = (0 - a.real * b.imag + a.imag * b.real) / den;

    return v;
}











