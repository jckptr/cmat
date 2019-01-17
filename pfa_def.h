

#ifndef CMAT_PFA_DEF_H
#define CMAT_PFA_DEF_H






typedef struct cplx      cplx;
typedef struct cmat      cmat;
typedef struct normalmat normalmat;
typedef struct sparsemat sparsemat;



struct cplx
{
    double real = 0.0;
    double imag = 0.0;
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
cplx ctrace(cmat a);

cmat   cmnr(cmat a, int r, int c);
double cdet(cmat a);
cmat   cadjg(cmat a);
cmat   cinv(cmat a);


cplx cadd(cplx a, cplx b);
cplx csub(cplx a, cplx b);
cplx cmul(cplx a, cplx b);
cplx cdiv(cplx a, cplx b);

cmat cadd(cmat a, cmat b);
cmat csub(cmat a, cmat b);
cmat cmul(cmat a, cmat b);
cmat cdiv(cmat a, cmat b);

void cadds(cmat a, cplx v, int r, int c);
void csubs(cmat a, cplx v, int r, int c);
void cmuls(cmat a, cplx v, int r, int c);
void cdivs(cmat a, cplx v, int r, int c);

cmat cadds(cmat a, cplx v);
cmat csubs(cmat a, cplx v);
cmat cmuls(cmat a, cplx v);
cmat cdivs(cmat a, cplx v);

cmat cadds(cmat a, cmat b);
cmat csubs(cmat a, cmat b);
cmat cmuls(cmat a, cmat b);
cmat cdivs(cmat a, cmat b);


#endif