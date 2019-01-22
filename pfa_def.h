

#ifndef CMAT_PFA_DEF_H
#define CMAT_PFA_DEF_H




#define CDBL_EPSILON 1.0e-09



#define NMAX_NORMALMAT 50
#define NMAX_SPARSEMAT 50



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

    normalmat * nmat;
    sparsemat * smat;
};





struct _csys_mem_
{
    int nnmat = 0;
    int _nmat_pos_[NMAX_NORMALMAT] = { 0 };
    normalmat _nmat_[NMAX_NORMALMAT];
};






cmat cnew(int r, int c);
cmat cdup(cmat a);
void cdel(cmat a);
void cclr();

int cisempt(cmat a);
int ciszero(cplx v);
int cisreal(cplx v);
int cisimag(cplx v);
int cisint(cplx v);
int ciszero(cmat a);
int cisreal(cmat a);
int cisimag(cmat a);
int cisint(cmat a);
int ceql(cplx a, cplx b);
int ceql(cmat a, cmat b);

void cshow(cplx v);
void cshow(cmat a);

cplx cget(cmat a, int r, int c);
void cset(cmat a, cplx v, int r, int c);

// d == 'r' concat as row order
// d == 'c' concat as col order
cmat ccat(cmat a, cmat b, char d);


cplx cceil(cplx v);
cplx cfloor(cplx v);
cplx cround(cplx v);
cmat cceil(cmat a);
cmat cfloor(cmat a);
cmat cround(cmat a);

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

cmat cmuls(cmat a, cmat b);
cmat cdivs(cmat a, cmat b);


cmat cone(int n);
cmat cone(int r, int c);
cmat ceye(int n);
cplx crand();
cmat crand(int n);
cmat crand(int r, int c);

cmat clns(double from, double to, int n);





#endif