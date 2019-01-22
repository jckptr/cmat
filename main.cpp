



#include <iostream>


#include "pfa_def.h"
#include <cmath>





int main()
{
    using namespace std;

    cmat a = cnew(2, 2);

    cmat b = cnew(2, 1);

    cplx t;
    t.real = 1;
    t.imag = 0;
    cset(a, t, 0, 0);
    t.real = 5;
    t.imag = 0;
    cset(a, t, 0, 1);
    t.real = 7;
    t.imag = 0;
    cset(a, t, 1, 0);
    t.real = 1;
    t.imag = 0;
    cset(a, t, 1, 1);

    t.real = 13;
    cset(b, t, 0, 0);
    t.real = 23;
    cset(b, t, 1, 0);

    cshow(a);
    cshow(b);

    cout << a.row << ' ' << a.col << endl;
    cout << b.row << ' ' << b.col << endl;

    cshow(ctrans(a));
    cshow(ctrace(a));


    cmat x = cdiv(b, a);
    cout << x.row << ' ' << x.col << endl;
    cshow(x);
    
    cshow(cmuls(b, x));
    cshow(cmuls(a, t));

    cmat o = cone(3);
    cmat p = ceye(3);

    cplx v = { 0.0, 1.0 };

    cshow(o);
    cshow(p);
    cshow(cmuls(o, v));
    cshow(cmuls(o, p));

    //cdel(a);
    //cdel(b);
    //cdel(x);

    cclr();

    system("pause");
    return 0;
}