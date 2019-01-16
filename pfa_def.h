

#ifndef CMAT_PFA_DEF_H
#define CMAT_PFA_DEF_H






typedef struct cplx      cplx;
typedef struct cmat      cmat;
typedef struct normalmat normalmat;
typedef struct sparsemat sparsemat;



struct cplx
{
    double real;
    double imag;
};

struct normalmat
{
    double * real;
    double * imag;
};


struct sparsemat
{
    
};

struct cmat
{
    int row = 0;
    int col = 0;
    
    int sparse = 0;

    normalmat nmat;
    sparsemat smat;
};





cmat cnew(int r, int c);
cmat cnew(double * real, double * imag, int r, int c);
void cdel(cmat a);

void cshow(cmat a);

cplx cget(cmat a, int r, int c);
void cset(cmat a, cplx v, int r, int c);


cmat ccpy(cmat a);
cmat cdup(cmat a);

cmat creal(cmat a);
cmat cimag(cmat a);
cmat cconj(cmat a);

cmat ctrans(cmat a);

cmat cadd(cmat a, cmat b);
cmat csub(cmat a, cmat b);
cmat cmul(cmat a, cmat b);
cmat cdiv(cmat a, cmat b);

void cadds(cmat a, cplx v, int r, int c);
void csubs(cmat a, cplx v, int r, int c);
void cmuls(cmat a, cplx v, int r, int c);
void cdivs(cmat a, cplx v, int r, int c);

cmat   cmnr(cmat a, int r, int c);
double cdet(cmat a);
cmat   cadjg(cmat a);
cmat   cinv(cmat a);



#endif